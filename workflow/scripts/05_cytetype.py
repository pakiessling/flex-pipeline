"""
05_cytetype.py — LLM-based cluster annotation via CyteType.

Requires:
  - adata.uns["rank_genes_groups"] in scanpy format (produced by 06_markers.py)
  - adata.var["gene_symbol"] column
  - CyteType installed: pip install cytetype

The CyteType annotation JSON is written to --output_json so it can be used
independently by the report step.
"""

import argparse
import datetime
import json
import logging
import os

import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


def _plot_cytetype_umaps(adata, results: dict, group_key: str, plot_dir: str):
    """Save UMAP PNGs coloured by CyteType annotations."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm:
        logger.warning("No UMAP coordinates found; skipping CyteType UMAP plots.")
        return

    os.makedirs(plot_dir, exist_ok=True)

    # Detect any obs columns the library added (e.g. cytetype_label)
    cytetype_cols = [c for c in adata.obs.columns if c.startswith("cytetype_")]

    # If no obs columns, build a mapping from the JSON and add it
    if not cytetype_cols and group_key in adata.obs.columns:
        annotations = results.get("annotations", [])
        if annotations:
            cluster_to_label = {str(a["clusterId"]): a["annotation"] for a in annotations}
            adata.obs["cytetype_label"] = (
                adata.obs[group_key].astype(str).map(cluster_to_label)
            )
            cytetype_cols = ["cytetype_label"]
            logger.info("Mapped CyteType cluster annotations to adata.obs['cytetype_label']")

    if not cytetype_cols:
        logger.warning("No CyteType annotation columns found; skipping UMAP plots.")
        return

    for col in cytetype_cols:
        fig, ax = plt.subplots(figsize=(11, 7))
        sc.pl.umap(
            adata, color=col, ax=ax, show=False,
            title=col.replace("_", " ").title(),
        )
        out_png = os.path.join(plot_dir, f"{col}_umap.png")
        fig.savefig(out_png, bbox_inches="tight", dpi=120)
        plt.close(fig)
        logger.info(f"UMAP plot saved → {out_png}")


def main(args):
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input not found: {args.input_file}")

    logger.info(f"Loading {args.input_file} …")
    adata = sc.read_h5ad(args.input_file)

    if "rank_genes_groups" not in adata.uns:
        raise KeyError(
            "adata.uns['rank_genes_groups'] is missing. "
            "Run the markers step (06_markers.py) before CyteType."
        )

    if args.group_key not in adata.obs.columns:
        raise KeyError(
            f"Group key '{args.group_key}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    # CyteType requires a gene_symbol column in .var
    if "gene_symbol" not in adata.var.columns:
        logger.info("Adding 'gene_symbol' column from var_names …")
        adata.var["gene_symbol"] = adata.var_names

    logger.info(
        f"Running CyteType on group_key='{args.group_key}', "
        f"n_top_genes={args.n_top_genes} …"
    )

    try:
        from cytetype import CyteType
    except ImportError:
        raise ImportError(
            "cytetype is not installed. Install with: pip install cytetype\n"
            "See https://github.com/NygenAnalytics/CyteType"
        )

    annotator = CyteType(
        adata,
        group_key=args.group_key,
        rank_key="rank_genes_groups",
        n_top_genes=args.n_top_genes,
        gene_symbols_column="gene_symbol",
    )

    adata = annotator.run(study_context=args.study_context)
    logger.info("CyteType annotation complete.")

    # Extract and save the annotation JSON
    results = json.loads(adata.uns["cytetype_results"]["result"])

    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
    with open(args.output_json, "w") as fh:
        json.dump(results, fh, indent=4)
    logger.info(f"Annotation JSON saved → {args.output_json}")

    # Reproducibility log
    import anndata as ad
    adata.uns.setdefault("pipeline_log", {})["cytetype"] = {
        "completed_at": datetime.datetime.now().isoformat(),
        "group_key": args.group_key,
        "n_top_genes": args.n_top_genes,
        "study_context": args.study_context,
        "software": {
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }

    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    adata.write_h5ad(args.output_file)
    logger.info(f"Annotated h5ad saved → {args.output_file}")

    _plot_cytetype_umaps(adata, results, args.group_key, args.plot_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LLM-based cluster annotation via CyteType.")
    parser.add_argument("--input_file",    required=True, help="Input .h5ad (must have rank_genes_groups in .uns)")
    parser.add_argument("--output_file",   required=True, help="Output .h5ad path")
    parser.add_argument("--output_json",   required=True, help="Output JSON path for CyteType annotation")
    parser.add_argument("--group_key",     default="leiden_3", help="Obs column with cluster labels")
    parser.add_argument("--n_top_genes",   type=int, default=50, help="Top marker genes sent to LLM per cluster")
    parser.add_argument("--study_context", default="Single-cell RNA-seq data from heart tissue.",
                        help="Biological context description for the LLM")
    parser.add_argument("--plot_dir",      default="",
                        help="Directory to save UMAP plots (default: same dir as output_file)")
    args = parser.parse_args()
    if not args.plot_dir:
        args.plot_dir = os.path.dirname(args.output_file)
    main(args)
