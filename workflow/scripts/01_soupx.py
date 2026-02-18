"""
01_soupx.py — Ambient RNA correction using SoupX.

If --skip is passed (or SoupX errors), the script falls back to the
unmodified filtered matrix so downstream steps always receive a valid h5ad.
"""

import argparse
import datetime
import logging
import os

import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)
rcb.logger.setLevel(logging.WARNING)

ro.pandas2ri.activate()
anndata2ri.activate()


def _load_filtered(filtered_matrix: str, sample: str):
    logger.info(f"[{sample}] Loading filtered matrix: {filtered_matrix}")
    adata = sc.read_10x_h5(filtered_matrix)
    adata.var_names_make_unique()
    return adata


def run_soupx(adata, raw_matrix: str, sample: str):
    """Run SoupX ambient RNA correction via rpy2. Returns (adata, success)."""
    logger.info(f"[{sample}] Preprocessing for SoupX clustering …")
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    soupx_groups = adata_pp.obs["soupx_groups"]
    del adata_pp

    logger.info(f"[{sample}] Loading raw matrix: {raw_matrix}")
    adata_raw = sc.read_10x_h5(raw_matrix)
    adata_raw.var_names_make_unique()
    common_genes = adata_raw.var_names.intersection(adata.var_names).sort_values()
    adata = adata[:, common_genes].copy()
    adata_raw = adata_raw[:, common_genes].copy()
    data_tod = adata_raw.X.T
    del adata_raw

    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["genes"] = genes
    ro.globalenv["cells"] = cells
    ro.globalenv["soupx_groups"] = soupx_groups

    logger.info(f"[{sample}] Running SoupX in R …")
    try:
        ro.r("""
            library(SoupX)
            library(Matrix)
            rownames(data) <- genes
            colnames(data) <- cells
            data     <- as(data,     "sparseMatrix")
            data_tod <- as(data_tod, "sparseMatrix")
            sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)
            soupProf <- data.frame(
                row.names = rownames(data),
                est    = rowSums(data) / sum(data),
                counts = rowSums(data)
            )
            sc <- setSoupProfile(sc, soupProf)
            sc <- setClusters(sc, soupx_groups)
            sc <- autoEstCont(sc, doPlot = FALSE)
            out <- adjustCounts(sc, roundToInt = TRUE)
        """)

        out = ro.globalenv["out"]
        before_sum = float(adata.X.sum())
        adata.layers["b4_soupx"] = adata.X.copy()
        adata.layers["after_soupx"] = out.T
        adata.X = adata.layers["after_soupx"].copy()
        after_sum = float(adata.X.sum())
        pct = (after_sum - before_sum) / before_sum * 100

        logger.info(
            f"[{sample}] SoupX complete. "
            f"Total counts before={before_sum:.0f}, after={after_sum:.0f} ({pct:+.2f}%)"
        )
        return adata, True

    except Exception as exc:
        logger.warning(f"[{sample}] SoupX failed ({exc}); falling back to raw filtered matrix.")
        adata.layers["b4_soupx"] = adata.X.copy()
        return adata, False


def main(args):
    adata = _load_filtered(args.filtered_matrix, args.sample)
    n_cells_raw = adata.n_obs

    soupx_success = False
    if args.skip:
        logger.info(f"[{args.sample}] SoupX skipped (--skip flag).")
        adata.layers["b4_soupx"] = adata.X.copy()
    else:
        adata, soupx_success = run_soupx(adata, args.raw_matrix, args.sample)

    # Add sample metadata
    adata.obs["Sample"] = args.sample

    # Reproducibility log
    import anndata as ad
    adata.uns.setdefault("pipeline_log", {})["soupx"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "soupx_run": soupx_success,
        "skipped": args.skip,
        "n_cells": int(n_cells_raw),
        "n_genes": int(adata.n_vars),
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    adata.write_h5ad(args.output)
    logger.info(
        f"[{args.sample}] Saved → {args.output} "
        f"(SoupX: {'yes' if soupx_success else 'no/skipped'})"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SoupX ambient RNA correction.")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, help="Output .h5ad path")
    parser.add_argument("--filtered_matrix", required=True, help="Filtered 10x H5 matrix path")
    parser.add_argument("--raw_matrix", default="", help="Raw 10x H5 matrix path (required unless --skip)")
    parser.add_argument("--skip", action="store_true", help="Skip SoupX, pass through filtered matrix only")
    args = parser.parse_args()

    if not args.skip and not args.raw_matrix:
        parser.error("--raw_matrix is required unless --skip is set")

    main(args)
