"""
02_qc.py — Per-sample quality control.

Cells are *marked* (not removed). Doublets are annotated by scDblFinder in the
preceding SoupX step (.obs["scDblFinder.class"] / .obs["scDblFinder.score"]).
Low-quality cells are flagged here in .obs["cell_quality"].
"""

import argparse
import datetime
import logging
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation as mad

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


def calculate_qc(adata):
    adata.var["MT"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["MT"], inplace=True, percent_top=[20], log1p=True
    )
    return adata


def flag_outliers(adata, nmads: int):
    """Mark cells as low-quality based on MAD thresholds. Does NOT remove cells."""
    def is_outlier(metric, upper_only=False):
        M = adata.obs[metric]
        if upper_only:
            return M > np.median(M) + nmads * mad(M)
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    sc.pp.filter_cells(adata, min_genes=2)

    outlier_mask = (
        is_outlier("log1p_total_counts")
        | is_outlier("log1p_n_genes_by_counts")
        | is_outlier("pct_counts_in_top_20_genes")
    )

    adata.obs["cell_quality"] = "high-quality"
    adata.obs.loc[outlier_mask, "cell_quality"] = "low-quality"

    n_low = outlier_mask.sum()
    logger.info(
        f"  Low-quality cells flagged: {n_low} / {adata.n_obs} "
        f"({100 * n_low / adata.n_obs:.1f}%)"
    )
    return adata



def cluster_and_embed(adata, leiden_resolutions):
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata)

    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        logger.info(
            f"  Leiden resolution={res}: "
            f"{adata.obs[key].nunique()} clusters"
        )

    # Rank genes before scaling (scaling can introduce negative values)
    for res in leiden_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        rg_key = f"rank_genes_groups_{str(res).replace('.', '_')}"
        sc.tl.rank_genes_groups(adata, key, method="wilcoxon", tie_correct=True, key_added=rg_key)

    sc.pp.scale(adata, zero_center=False)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
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

    # Top marker gene CSVs
    for res in leiden_resolutions:
        rg_key = f"rank_genes_groups_{str(res).replace('.', '_')}"
        if rg_key in adata.uns:
            df = pd.DataFrame(adata.uns[rg_key]["names"]).head(100)
            csv_path = os.path.join(sample_dir, f"{sample}_markers_{str(res).replace('.', '_')}.csv")
            df.to_csv(csv_path)
            logger.info(f"  Marker genes saved → {csv_path}")


def main(args):
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input not found: {args.input}")

    leiden_resolutions = [float(r) for r in args.leiden_resolutions.split()]

    logger.info(f"[{args.sample}] Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    n_start = adata.n_obs
    logger.info(f"[{args.sample}] Loaded {n_start} cells × {adata.n_vars} genes")

    # Append sample suffix to cell barcodes to ensure uniqueness after merging
    adata.obs["Sample"] = args.sample
    adata.obs.index = adata.obs.index + "_" + args.sample

    logger.info(f"[{args.sample}] Calculating QC metrics …")
    adata = calculate_qc(adata)

    logger.info(f"[{args.sample}] Flagging outlier cells (MAD threshold={args.mad_threshold}) …")
    adata = flag_outliers(adata, nmads=args.mad_threshold)

    logger.info(f"[{args.sample}] Clustering and embedding …")
    adata = cluster_and_embed(adata, leiden_resolutions)

    logger.info(f"[{args.sample}] Saving QC plots …")
    save_qc_plots(adata, args.sample, args.qc_folder, leiden_resolutions)

    # Reproducibility log
    import anndata as ad
    n_dbl = int((adata.obs["scDblFinder.class"] == "doublet").sum()) \
        if "scDblFinder.class" in adata.obs else 0
    adata.uns.setdefault("pipeline_log", {})["qc"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "n_cells_input": int(n_start),
        "n_cells_output": int(adata.n_obs),
        "n_low_quality": int((adata.obs["cell_quality"] == "low-quality").sum()),
        "n_scdblfinder_doublets": n_dbl,
        "mad_threshold": args.mad_threshold,
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    adata.write_h5ad(args.output)
    logger.info(f"[{args.sample}] Saved → {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Per-sample QC for single-cell RNA-seq.")
    parser.add_argument("--input",   required=True, help="Input .h5ad path")
    parser.add_argument("--sample",  required=True, help="Sample identifier")
    parser.add_argument("--output",  required=True, help="Output .h5ad path")
    parser.add_argument("--qc_folder", required=True, help="Folder for QC plots and CSVs")
    parser.add_argument("--min_genes", type=int, default=2, help="Minimum genes per cell")
    parser.add_argument("--mad_threshold", type=int, default=5, help="MAD multiplier for outlier detection")
    parser.add_argument("--leiden_resolutions", default="1.5 3.0",
                        help="Space-separated Leiden clustering resolutions")
    args = parser.parse_args()
    main(args)
