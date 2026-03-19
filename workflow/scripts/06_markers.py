"""
06_markers.py — Differential expression / marker genes via illico.

Uses illico.asymptotic_wilcoxon for fast Wilcoxon rank-sum tests and stores
the result in both:
  - adata.uns["rank_genes_groups"]  (scanpy-compatible format for CyteType)
  - CSV files in --markers_dir
"""

import argparse
import datetime
import logging
import os
import random

import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multitest import multipletests

os.environ["PYTHONHASHSEED"] = "0"
np.random.seed(0)
random.seed(0)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def add_adjusted_pvalues(de_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """FDR correction per group, adds a pval_adj column."""
    padj_parts = []
    for _, grp in de_df.groupby(group_col, sort=False):
        _, padj, _, _ = multipletests(grp["p_value"].values, method="fdr_bh")
        padj_parts.append(pd.Series(padj, index=grp.index))
    de_df["pval_adj"] = pd.concat(padj_parts)
    return de_df


def illico_to_rank_genes_groups(
    de_df: pd.DataFrame, group_col: str, n_top: int, pval_cutoff: float
) -> dict:
    """
    Convert the flat illico output DataFrame to scanpy's rank_genes_groups
    structured-recarray format.

    Genes are filtered to significant (pval_adj <= cutoff) and upregulated
    (fold_change > 1), then ranked by descending log2 fold change.
    """
    groups = sorted(de_df[group_col].astype(str).unique())
    gene_dtype = [(g, "U200") for g in groups]
    float32_dtype = [(g, np.float32) for g in groups]
    float64_dtype = [(g, np.float64) for g in groups]

    names_arr = np.full(n_top, "", dtype=gene_dtype)
    scores_arr = np.zeros(n_top, dtype=float32_dtype)
    pvals_arr = np.ones(n_top, dtype=float64_dtype)
    pvals_adj_arr = np.ones(n_top, dtype=float64_dtype)
    logfcs_arr = np.zeros(n_top, dtype=float32_dtype)

    for g in groups:
        grp = de_df[de_df[group_col].astype(str) == g].copy()
        grp = grp[(grp["pval_adj"] <= pval_cutoff) & (grp["fold_change"] > 1.0)]
        grp["log2fc"] = np.log2(grp["fold_change"].values)
        grp = grp.sort_values("log2fc", ascending=False).reset_index(drop=True)

        n = min(len(grp), n_top)
        if n == 0:
            logger.warning(f"Group '{g}': no significant upregulated genes at FDR <= {pval_cutoff}")
            continue

        names_arr[g][:n] = grp["gene"].values[:n]
        scores_arr[g][:n] = grp["statistic"].values[:n].astype(np.float32)
        pvals_arr[g][:n] = grp["p_value"].values[:n]
        pvals_adj_arr[g][:n] = grp["pval_adj"].values[:n]
        logfcs_arr[g][:n] = grp["log2fc"].values[:n].astype(np.float32)

    return {
        "params": {
            "groupby": group_col,
            "reference": "rest",
            "method": "wilcoxon",
            "use_raw": False,
            "layer": None,
            "corr_method": "benjamini-hochberg",
        },
        "names": names_arr,
        "scores": scores_arr,
        "pvals": pvals_arr,
        "pvals_adj": pvals_adj_arr,
        "logfoldchanges": logfcs_arr,
    }


def main(args):
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input not found: {args.input_file}")

    logger.info(f"Loading {args.input_file}")
    adata = sc.read_h5ad(args.input_file)

    if args.group_key not in adata.obs.columns:
        raise KeyError(
            f"Group key '{args.group_key}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    logger.info(
        f"Running illico.asymptotic_wilcoxon on '{args.group_key}' "
        f"({adata.obs[args.group_key].nunique()} groups, {adata.n_obs} cells)"
    )

    try:
        from illico import asymptotic_wilcoxon
    except ImportError:
        raise ImportError(
            "illico is not installed. Install with: pip install illico\n"
            "See https://github.com/remydubois/illico"
        )

    de_df = asymptotic_wilcoxon(
        adata,
        group_keys=args.group_key,
        reference=None,
        is_log1p=args.is_log1p,
    )

    de_df = de_df.reset_index()
    de_df = de_df.rename(columns={"pert": args.group_key, "feature": "gene"})
    group_col = args.group_key

    logger.info(f"illico returned {len(de_df)} rows for {de_df[group_col].nunique()} groups")

    de_df = add_adjusted_pvalues(de_df, group_col)
    de_df["log2fc"] = np.where(
        de_df["fold_change"] > 0,
        np.log2(de_df["fold_change"]),
        np.nan,
    )

    os.makedirs(args.markers_dir, exist_ok=True)

    all_genes_path = os.path.join(args.markers_dir, f"all_genes_{args.group_key}.csv")
    de_df.to_csv(all_genes_path, index=False)
    logger.info(f"All DE genes saved -> {all_genes_path}")

    sig_up = de_df[(de_df["pval_adj"] <= args.pval_cutoff) & (de_df["fold_change"] > 1.0)].copy()
    n_sig = len(sig_up)
    logger.info(
        f"{n_sig} significant upregulated genes at FDR <= {args.pval_cutoff} "
        f"(out of {len(de_df)} total)"
    )

    top_markers = (
        sig_up.sort_values([group_col, "fold_change"], ascending=[True, False])
        .groupby(group_col, sort=False)
        .head(args.n_top)
        .reset_index(drop=True)
    )
    markers_path = os.path.join(args.markers_dir, f"marker_genes_{args.group_key}.csv")
    top_markers.to_csv(markers_path, index=False)
    logger.info(f"Top marker genes saved -> {markers_path}")

    logger.info("Converting to scanpy rank_genes_groups format")
    adata.uns["rank_genes_groups"] = illico_to_rank_genes_groups(
        de_df, group_col, n_top=args.n_top, pval_cutoff=args.pval_cutoff
    )

    import anndata as ad

    adata.uns.setdefault("pipeline_log", {})["markers"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "tool": "illico",
        "group_key": args.group_key,
        "n_top": args.n_top,
        "is_log1p": args.is_log1p,
        "pval_cutoff": args.pval_cutoff,
        "n_groups": int(de_df[group_col].nunique()),
        "n_significant_upregulated": n_sig,
        "filter": f"pval_adj <= {args.pval_cutoff} & fold_change > 1, sorted by log2fc desc",
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }

    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    adata.write_h5ad(args.output_file)
    logger.info(f"Saved h5ad with rank_genes_groups -> {args.output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Marker gene analysis via illico.")
    parser.add_argument("--input_file", required=True, help="Input .h5ad path")
    parser.add_argument("--output_file", required=True, help="Output .h5ad path")
    parser.add_argument(
        "--markers_dir", required=True, help="Directory for marker CSV files"
    )
    parser.add_argument(
        "--group_key", default="leiden_3_0", help="Obs column with cluster labels"
    )
    parser.add_argument(
        "--n_top", type=int, default=100, help="Top N markers to report per cluster"
    )
    parser.add_argument(
        "--pval_cutoff", type=float, default=0.05, help="FDR-adjusted p-value cutoff"
    )
    parser.add_argument(
        "--is_log1p",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether expression matrix is already log1p-transformed",
    )
    args = parser.parse_args()
    main(args)