"""
02_qc.py — Per-sample quality control.

Cells are *marked* (not removed). Doublets are annotated by scDblFinder in the
preceding SoupX step (.obs["scDblFinder.class"] / .obs["scDblFinder.score"]).
Low-quality cells are flagged here in .obs["cell_quality"]. Only zero-total
count cells are removed. Optional diagnostics run separately with --diagnostics.
"""

import argparse
import datetime
import json
import logging
import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from processing_utils import (InsufficientData, prepare_counts, select_hvgs, run_pca_neighbors)
from scipy.stats import median_abs_deviation as mad
from statsmodels.stats.multitest import multipletests

os.environ["PYTHONHASHSEED"] = "0"
import random

np.random.seed(0)
random.seed(0)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def illico_to_rank_genes_groups(de_df: pd.DataFrame, group_col: str, n_top: int) -> dict:
    """Convert flat illico output to scanpy's rank_genes_groups structured-recarray format."""
    groups = sorted(de_df[group_col].astype(str).unique())
    gene_dtype    = [(g, "U200")     for g in groups]
    float32_dtype = [(g, np.float32) for g in groups]
    float64_dtype = [(g, np.float64) for g in groups]

    names_arr     = np.empty(n_top, dtype=gene_dtype)
    scores_arr    = np.zeros(n_top, dtype=float32_dtype)
    pvals_arr     = np.ones(n_top,  dtype=float64_dtype)
    pvals_adj_arr = np.ones(n_top,  dtype=float64_dtype)
    logfcs_arr    = np.zeros(n_top, dtype=float32_dtype)

    for g in groups:
        grp = de_df[de_df[group_col].astype(str) == g].copy()
        grp = grp.sort_values("statistic", ascending=False).reset_index(drop=True)
        n = min(len(grp), n_top)
        _, padj, _, _ = multipletests(grp["p_value"].values, method="fdr_bh")
        names_arr[g][:n]     = grp["gene"].values[:n]
        scores_arr[g][:n]    = grp["statistic"].values[:n].astype(np.float32)
        pvals_arr[g][:n]     = grp["p_value"].values[:n]
        pvals_adj_arr[g][:n] = padj[:n]
        logfcs_arr[g][:n]    = np.log2(
            np.maximum(grp["fold_change"].values[:n], 1e-10)
        ).astype(np.float32)

    return {
        "params": {
            "groupby": group_col, "reference": "rest",
            "method": "wilcoxon", "use_raw": False,
            "layer": None, "corr_method": "benjamini-hochberg",
        },
        "names": names_arr, "scores": scores_arr,
        "pvals": pvals_arr, "pvals_adj": pvals_adj_arr,
        "logfoldchanges": logfcs_arr,
    }


def calculate_qc(adata):
    adata.var["MT"] = adata.var_names.str.startswith("MT-")
    if not adata.n_obs or not adata.n_vars:
        for col in ("total_counts", "n_genes_by_counts", "log1p_total_counts",
                    "log1p_n_genes_by_counts", "pct_counts_MT", "pct_counts_in_top_20_genes"):
            adata.obs[col] = np.zeros(adata.n_obs)
        return adata
    top = min(20, adata.n_vars)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["MT"], inplace=True, percent_top=[top], log1p=True
    )
    if top != 20:
        adata.obs["pct_counts_in_top_20_genes"] = adata.obs[f"pct_counts_in_top_{top}_genes"]
    return adata


def flag_outliers(adata, nmads: int, min_genes: int = 2):
    """Mark cells as low-quality based on MAD thresholds. Does NOT remove cells."""

    def is_outlier(metric, upper_only=False):
        M = adata.obs[metric]
        if M.empty:
            return pd.Series(False, index=M.index)
        if upper_only:
            return M > np.median(M) + nmads * mad(M)
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    outlier_mask = (
        (adata.obs["n_genes_by_counts"] < min_genes)
        | is_outlier("log1p_total_counts")
        | is_outlier("log1p_n_genes_by_counts")
        | is_outlier("pct_counts_in_top_20_genes")
    )

    adata.obs["cell_quality"] = "high-quality"
    adata.obs.loc[outlier_mask, "cell_quality"] = "low-quality"

    n_low = outlier_mask.sum()
    logger.info(
        f"  Low-quality cells flagged: {n_low} / {adata.n_obs} "
        f"({100 * n_low / max(adata.n_obs, 1):.1f}%)"
    )
    return adata


def cluster_and_embed(adata, leiden_resolutions):
    select_hvgs(adata)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    run_pca_neighbors(adata)

    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        logger.info(f"  Leiden resolution={res}: {adata.obs[key].nunique()} clusters")

    # Rank genes before scaling (scaling can introduce negative values)
    try:
        from illico import asymptotic_wilcoxon
    except ImportError:
        raise ImportError(
            "illico is not installed. Install with: pip install illico\n"
            "See https://github.com/remydubois/illico"
        )
    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        rg_key = f"rank_genes_groups_{str(res).replace('.', '_')}"

        if adata.obs[key].nunique() < 2:
            logger.info("Skipping markers for %s: only one cluster", key)
            continue
        try:
            de_df = asymptotic_wilcoxon(adata, group_keys=key, reference=None, is_log1p=True)
        except np.linalg.LinAlgError as error:
            logger.warning("Skipping markers for %s: %s", key, error)
            adata.uns.setdefault("diagnostic_failures", {})[rg_key] = str(error)
            continue
        de_df = de_df.reset_index()
        de_df = de_df.rename(columns={"pert": key, "feature": "gene"})
        adata.uns[rg_key] = illico_to_rank_genes_groups(de_df, key, n_top=100)

    sc.pp.scale(adata, zero_center=False)
    run_pca_neighbors(adata)
    try:
        sc.tl.umap(adata, init_pos="random" if adata.n_obs < 5 else "spectral")
    except np.linalg.LinAlgError as error:
        logger.warning("Skipping UMAP: %s", error)
        adata.uns.setdefault("diagnostic_failures", {})["umap"] = str(error)
    return adata


def save_qc_plots(adata, sample: str, qc_folder: str, leiden_resolutions):
    sample_dir = os.path.join(qc_folder, sample)
    os.makedirs(sample_dir, exist_ok=True)

    color_cols = [
        "scDblFinder.class",
        "scDblFinder.score",
        "cell_quality",
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_MT",
    ] + [f"leiden_{str(r).replace('.', '_')}" for r in leiden_resolutions]

    # Filter to columns that exist
    color_cols = [c for c in color_cols if c in adata.obs.columns]

    if "X_umap" in adata.obsm:
        fig = sc.pl.umap(
            adata,
            color=color_cols,
            legend_loc="on data",
            return_fig=True,
            vmax="p99",
            sort_order=False,
        )
        umap_path = os.path.join(sample_dir, f"{sample}_umap.png")
        fig.savefig(umap_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
        logger.info(f"  UMAP saved → {umap_path}")

    else:
        logger.info("No UMAP available for %s", sample)

    # Top marker gene CSVs
    for res in leiden_resolutions:
        rg_key = f"rank_genes_groups_{str(res).replace('.', '_')}"
        if rg_key in adata.uns:
            df = pd.DataFrame(adata.uns[rg_key]["names"]).head(100)
            csv_path = os.path.join(
                sample_dir, f"{sample}_markers_{str(res).replace('.', '_')}.csv"
            )
            df.to_csv(csv_path)
            logger.info(f"  Marker genes saved → {csv_path}")


def main(args):
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input not found: {args.input}")

    leiden_resolutions = [float(r) for r in args.leiden_resolutions.split()]

    logger.info(f"[{args.sample}] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    if args.diagnostics:
        status = adata.uns["sample_status"]
        diagnostic_status = {"status": "skipped", "reason": str(status["reason"])}
        if status["integration_eligible"]:
            try:
                adata = cluster_and_embed(adata, leiden_resolutions)
                save_qc_plots(adata, args.sample, args.qc_folder, leiden_resolutions)
                failures = adata.uns.get("diagnostic_failures", {})
                diagnostic_status = {
                    "status": "partial" if failures else "completed",
                    "reason": "; ".join(f"{k}: {v}" for k, v in failures.items()),
                }
            except (InsufficientData, np.linalg.LinAlgError) as error:
                logger.warning("Skipping diagnostics: %s", error)
                diagnostic_status = {"status": "skipped", "reason": str(error)}
        adata.uns["diagnostics"] = diagnostic_status
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        adata.write_h5ad(args.output)
        return
    n_start = adata.n_obs
    adata, status = prepare_counts(adata, args.sample)
    logger.info("[%s] Eligibility: %s; removed %d zero-count cells; %s",
                args.sample, status["integration_eligible"],
                status["n_zero_count_cells_removed"], status["reason"])
    # JSON keeps a complete audit of removed barcode IDs.
    adata.uns["sample_status"] = {k: v for k, v in status.items() if k != "removed_cell_ids"}
    logger.info(f"[{args.sample}] Loaded {n_start} cells × {adata.n_vars} genes")

    # Append sample suffix to cell barcodes to ensure uniqueness after merging
    adata.obs["Sample"] = args.sample
    adata.obs.index = adata.obs.index + "_" + args.sample

    logger.info(f"[{args.sample}] Calculating QC metrics …")
    adata = calculate_qc(adata)

    logger.info(
        f"[{args.sample}] Flagging outlier cells (MAD threshold={args.mad_threshold}) …"
    )
    adata = flag_outliers(adata, nmads=args.mad_threshold, min_genes=args.min_genes)

    # Reproducibility log
    import anndata as ad

    n_dbl = (
        int((adata.obs["scDblFinder.class"] == "doublet").sum())
        if "scDblFinder.class" in adata.obs
        else 0
    )
    adata.uns.setdefault("pipeline_log", {})["qc"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "n_cells_input": int(n_start),
        "n_cells_output": int(adata.n_obs),
        "n_low_quality": int((adata.obs["cell_quality"] == "low-quality").sum()),
        "n_scdblfinder_doublets": n_dbl,
        "mad_threshold": args.mad_threshold,
        "min_genes_annotation_only": args.min_genes,
        "cell_removal_policy": "zero_total_counts_only",
        "integration_min_cells": 3,
        "integration_min_variable_genes": 3,
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output)
    status_path = args.status_output or args.output + ".status.json"
    os.makedirs(os.path.dirname(status_path) or ".", exist_ok=True)
    with open(status_path, "w") as fh:
        json.dump(status, fh, indent=2)
    logger.info(f"[{args.sample}] Saved → {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Per-sample QC for single-cell RNA-seq."
    )
    parser.add_argument("--diagnostics", action="store_true", help="Run optional exploration on an existing QC result")
    parser.add_argument("--status_output", help="Sample eligibility JSON")
    parser.add_argument("--input", required=True, help="Input .h5ad path")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, help="Output .h5ad path")
    parser.add_argument(
        "--qc_folder", required=True, help="Folder for QC plots and CSVs"
    )
    parser.add_argument(
        "--min_genes", type=int, default=2, help="Minimum genes for quality annotation only; cells are retained"
    )
    parser.add_argument(
        "--mad_threshold",
        type=int,
        default=5,
        help="MAD multiplier for outlier detection",
    )
    parser.add_argument(
        "--leiden_resolutions",
        default="1.5 3.0",
        help="Space-separated Leiden clustering resolutions",
    )
    args = parser.parse_args()
    main(args)
