"""
03_integration.py — Multi-sample integration via Harmony.

Loads eligible QC files from an explicit manifest, computes HVGs, runs Harmony,
generates UMAP and PaCMAP embeddings, and clusters at two Leiden resolutions.
Doublets are annotated in obs["scDblFinder.class"] but not removed here.
"""

import argparse
import datetime
import json
import logging
import os

import harmonypy as hm
import numpy as np
import pandas as pd
import pacmap
from processing_utils import select_hvgs, run_pca_neighbors
import scanpy as sc

os.environ["PYTHONHASHSEED"] = "0"
import random

np.random.seed(0)
random.seed(0)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def load_samples(manifest):
    included = manifest.loc[manifest["status"] == "included"]
    if included.empty:
        raise ValueError("No samples are eligible for integration; see sample_inclusion.csv")
    adatas = []
    for row in included.itertuples(index=False):
        p = row.h5ad_path
        logger.info("  Loading %s", p)
        a = sc.read_h5ad(p)
        status = a.uns.get("sample_status", {})
        if not status.get("integration_eligible", False) or status.get("sample") != row.sample:
            raise ValueError(f"QC status does not match manifest: {p}")
        if a.n_obs != row.n_cells_retained:
            raise ValueError(f"QC cell count does not match manifest: {p}")
        if set(a.obs["Sample"].astype(str)) != {row.sample}:
            raise ValueError(f"Sample labels do not match manifest: {p}")
        sample_name = row.sample
        if "after_soupx" in a.layers:
            a.X = a.layers["after_soupx"].copy()
            a.obs["SoupX_run"] = True
            logger.info(f"    [{sample_name}] Using 'after_soupx' layer")
        elif "b4_soupx" in a.layers:
            a.X = a.layers["b4_soupx"].copy()
            a.obs["SoupX_run"] = False
            logger.info(f"    [{sample_name}] Using 'b4_soupx' layer (SoupX not run)")
        elif "counts" in a.layers:
            a.X = a.layers["counts"].copy()
            a.obs["SoupX_run"] = False
        else:
            raise ValueError(f"[{sample_name}] QC output is missing its count layer")

        adatas.append(a)

    return adatas


def main(args):
    leiden_resolutions = [float(r) for r in args.leiden_resolutions.split()]

    manifest = pd.read_csv(args.manifest, keep_default_na=False, dtype={"sample": str})
    if manifest["sample"].duplicated().any():
        raise ValueError("Duplicate samples in manifest")
    if not manifest["status"].isin(["included", "excluded"]).all():
        raise ValueError("Invalid sample status in manifest")
    adatas = load_samples(manifest)
    n_samples = len(adatas)
    plot_dir = os.path.dirname(args.output_file) or "."
    os.makedirs(plot_dir, exist_ok=True)
    total_before = sum(a.n_obs for a in adatas)
    logger.info(f"Loaded {len(adatas)} samples, {total_before} total cells")

    # Harvest per-step software versions from the first sample before concat
    # drops them (sc.concat merge="same" discards uns keys that differ across samples).
    _r_step_logs = {}
    for _step in ("soupx", "scdblfinder"):
        _entry = adatas[0].uns.get("pipeline_log", {}).get(_step)
        if _entry and "software" in _entry:
            _r_step_logs[_step] = {"software": _entry["software"]}

    adata = sc.concat(adatas, join="outer", merge="same", fill_value=0)
    adata.uns["sample_inclusion"] = manifest.drop(columns=["h5ad_path"]).copy()
    del adatas
    logger.info(f"Merged: {adata.n_obs} cells")

    # Re-attach R software version stubs so they survive into the final h5ad.
    for _step, _info in _r_step_logs.items():
        adata.uns.setdefault("pipeline_log", {}).setdefault(_step, {}).update(_info)

    if "scDblFinder.class" in adata.obs.columns:
        n_dbl = (adata.obs["scDblFinder.class"] == "doublet").sum()
        logger.info(
            f"  {n_dbl} doublets annotated in obs['scDblFinder.class'] "
            f"({100 * n_dbl / adata.n_obs:.1f}%) — kept, not removed"
        )

    logger.info(f"Computing {args.n_top_genes} highly variable genes …")
    select_hvgs(adata, n_top_genes=args.n_top_genes, batch_key="Sample")
    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"  {n_hvg} HVGs selected")

    logger.info("Normalising, log-transforming, scaling …")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False)
    run_pca_neighbors(adata)

    if n_samples > 1:
        logger.info(f"Running Harmony (max_iter={args.max_iter_harmony}) …")
        ho = hm.run_harmony(
            adata.obsm["X_pca"], adata.obs, "Sample",
            max_iter_harmony=args.max_iter_harmony, random_state=0,
            nclust=min(100, max(2, adata.n_obs // 30)),
        )
        corrected = ho.Z_corr
        # harmonypy versions differ in orientation.
        if corrected.shape != adata.obsm["X_pca"].shape:
            corrected = corrected.T
        adata.obsm["X_pca_harmony"] = corrected
    else:
        logger.info("One eligible sample: skipping Harmony batch correction")
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"].copy()
    sc.pp.neighbors(adata, use_rep="X_pca_harmony",
                    n_neighbors=min(15, adata.n_obs - 1))

    logger.info("Leiden clustering …")
    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        logger.info(f"  resolution={res}: {adata.obs[key].nunique()} clusters")

    logger.info("Cell cycle scoring …")
    cell_cycle_scored = False
    if args.cell_cycle_genes and os.path.exists(args.cell_cycle_genes):
        with open(args.cell_cycle_genes) as fh:
            cc = json.load(fh)
        s_genes = [g for g in cc["s_genes"] if g in adata.var_names]
        g2m_genes = [g for g in cc["g2m_genes"] if g in adata.var_names]
        if s_genes and g2m_genes and adata.n_vars >= 50:
            sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)
            cell_cycle_scored = True
        else:
            logger.warning("Skipping cell cycle scoring: insufficient matching genes")
            adata.obs["phase"] = "unassigned"
        logger.info(
            f"  Phase distribution: {adata.obs['phase'].value_counts().to_dict()}"
        )
    else:
        logger.warning(
            f"Cell cycle genes file not found: {args.cell_cycle_genes!r}. Skipping."
        )

    logger.info("Computing UMAP …")
    sc.tl.umap(adata, init_pos="random" if adata.n_obs < 5 else "spectral")

    logger.info("Computing PaCMAP …")
    if adata.n_obs >= 20:
        embedding = pacmap.PaCMAP(n_neighbors=min(10, adata.n_obs - 1))
        adata.obsm["X_pacmap"] = embedding.fit_transform(adata.obsm["X_pca_harmony"])
        adata.uns["pacmap_status"] = "completed"
    else:
        adata.uns["pacmap_status"] = "skipped: fewer than 20 cells"
        logger.info(adata.uns["pacmap_status"])

    logger.info("Saving UMAP plots …")
    sc.settings.figdir = plot_dir

    ri = np.random.permutation(adata.n_obs)
    sc.pl.umap(adata[ri, :], color="Sample", save="_by_sample.png", show=False)

    cluster_cols = [f"leiden_{str(r).replace('.', '_')}" for r in leiden_resolutions]
    sc.pl.umap(
        adata,
        color=cluster_cols,
        legend_loc="on data",
        save="_clusters.png",
        show=False,
    )

    qc_cols = [c for c in ["scDblFinder.class", "cell_quality"] if c in adata.obs.columns]
    if qc_cols:
        sc.pl.umap(adata[ri, :], color=qc_cols, save="_qc.png", show=False)

    # Reproducibility log
    import anndata as ad

    adata.uns.setdefault("pipeline_log", {})["integration"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "n_samples": n_samples,
        "n_samples_excluded": int((manifest["status"] == "excluded").sum()),
        "harmony_run": n_samples > 1,
        "n_cells_output": int(adata.n_obs),
        "n_hvg": int(n_hvg),
        "leiden_resolutions": leiden_resolutions,
        "max_iter_harmony": args.max_iter_harmony,
        "cell_cycle_scoring": cell_cycle_scored,
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
            "pacmap": pacmap.__version__,
        },
    }

    os.makedirs(plot_dir, exist_ok=True)
    adata.write_h5ad(args.output_file)
    logger.info(f"Saved integrated data → {args.output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Integrate scRNA-seq samples via Harmony."
    )
    parser.add_argument(
        "--manifest", required=True, help="Explicit sample inclusion CSV from QC"
    )
    parser.add_argument("--output_file", required=True, help="Output .h5ad path")
    parser.add_argument("--n_top_genes", type=int, default=4000)
    parser.add_argument(
        "--leiden_resolutions",
        default="1.5 3.0",
        help="Space-separated list of Leiden resolutions",
    )
    parser.add_argument("--max_iter_harmony", type=int, default=100)
    parser.add_argument(
        "--cell_cycle_genes",
        default="config/cell_cycle_genes.json",
        help="Path to cell_cycle_genes.json",
    )
    args = parser.parse_args()
    main(args)
