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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LLM-based cluster annotation via CyteType.")
    parser.add_argument("--input_file",    required=True, help="Input .h5ad (must have rank_genes_groups in .uns)")
    parser.add_argument("--output_file",   required=True, help="Output .h5ad path")
    parser.add_argument("--output_json",   required=True, help="Output JSON path for CyteType annotation")
    parser.add_argument("--group_key",     default="leiden_3", help="Obs column with cluster labels")
    parser.add_argument("--n_top_genes",   type=int, default=50, help="Top marker genes sent to LLM per cluster")
    parser.add_argument("--study_context", default="Single-cell RNA-seq data from heart tissue.",
                        help="Biological context description for the LLM")
    args = parser.parse_args()
    main(args)
