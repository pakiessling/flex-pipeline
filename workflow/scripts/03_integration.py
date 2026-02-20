"""
03_integration.py — Multi-sample integration via Harmony.

Loads all per-sample h5ad files, computes HVGs, runs Harmony batch correction,
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
import pacmap
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


def load_samples(input_dir: str):
    paths = [
        os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".h5ad")
    ]
    if not paths:
        raise FileNotFoundError(
            f"No .h5ad files found in {input_dir!r}. "
            "Check that the QC step completed successfully."
        )

    adatas = []
    for p in paths:
        logger.info(f"  Loading {p}")
        a = sc.read_h5ad(p)

        sample_name = a.obs.get("Sample", ["?"])[0]
        if "after_soupx" in a.layers:
            a.X = a.layers["after_soupx"].copy()
            a.obs["SoupX_run"] = True
            logger.info(f"    [{sample_name}] Using 'after_soupx' layer")
        elif "b4_soupx" in a.layers:
            a.X = a.layers["b4_soupx"].copy()
            a.obs["SoupX_run"] = False
            logger.info(f"    [{sample_name}] Using 'b4_soupx' layer (SoupX not run)")
        else:
            a.obs["SoupX_run"] = False
            logger.warning(f"    [{sample_name}] No SoupX layer found; using .X as-is")

        adatas.append(a)

    return adatas


def main(args):
    leiden_resolutions = [float(r) for r in args.leiden_resolutions.split()]

    logger.info(f"Loading samples from {args.input_dir} …")
    adatas = load_samples(args.input_dir)
    total_before = sum(a.n_obs for a in adatas)
    logger.info(f"Loaded {len(adatas)} samples, {total_before} total cells")

    adata = sc.concat(adatas, join="outer", merge="same")
    del adatas
    logger.info(f"Merged: {adata.n_obs} cells")

    if "scDblFinder.class" in adata.obs.columns:
        n_dbl = (adata.obs["scDblFinder.class"] == "doublet").sum()
        logger.info(
            f"  {n_dbl} doublets annotated in obs['scDblFinder.class'] "
            f"({100 * n_dbl / adata.n_obs:.1f}%) — kept, not removed"
        )

    logger.info(f"Computing {args.n_top_genes} highly variable genes …")
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3_paper",
        n_top_genes=args.n_top_genes,
        batch_key="Sample",
    )
    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"  {n_hvg} HVGs selected")

    logger.info("Normalising, log-transforming, scaling …")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False)
    sc.pp.pca(adata)

    logger.info(f"Running Harmony (max_iter={args.max_iter_harmony}) …")
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        "Sample",
        max_iter_harmony=args.max_iter_harmony,
        random_state=0,
    )
    adata.obsm["X_pca_harmony"] = ho.Z_corr
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")

    logger.info("Leiden clustering …")
    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        logger.info(f"  resolution={res}: {adata.obs[key].nunique()} clusters")

    logger.info("Cell cycle scoring …")
    if args.cell_cycle_genes and os.path.exists(args.cell_cycle_genes):
        with open(args.cell_cycle_genes) as fh:
            cc = json.load(fh)
        sc.tl.score_genes_cell_cycle(adata, cc["s_genes"], cc["g2m_genes"])
        logger.info(
            f"  Phase distribution: {adata.obs['phase'].value_counts().to_dict()}"
        )
    else:
        logger.warning(
            f"Cell cycle genes file not found: {args.cell_cycle_genes!r}. Skipping."
        )

    logger.info("Computing UMAP …")
    sc.tl.umap(adata)

    logger.info("Computing PaCMAP …")
    embedding = pacmap.PaCMAP()
    adata.obsm["X_pacmap"] = embedding.fit_transform(adata.obsm["X_pca_harmony"])

    logger.info("Saving UMAP plots …")
    plot_dir = os.path.dirname(args.output_file)
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
        "n_samples": len(adata.obs["Sample"].unique()),
        "n_cells_output": int(adata.n_obs),
        "n_hvg": int(n_hvg),
        "leiden_resolutions": leiden_resolutions,
        "max_iter_harmony": args.max_iter_harmony,
        "cell_cycle_scoring": os.path.exists(args.cell_cycle_genes),
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
        "--input_dir", required=True, help="Directory with per-sample .h5ad files"
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
