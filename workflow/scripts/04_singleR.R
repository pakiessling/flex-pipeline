#!/usr/bin/env Rscript
# 04_singleR.R — Label transfer from a reference dataset using SingleR.
#
# Reads both query and reference as h5ad files via anndataR, runs SingleR,
# and writes the query h5ad back with SingleR labels added to .obs.

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleR)
  library(BiocParallel)
  library(argparser)
})

parser <- arg_parser("SingleR label transfer")
parser <- add_argument(parser, "--query",        help = "Path to query .h5ad (integrated)")
parser <- add_argument(parser, "--reference",    help = "Path to reference .h5ad")
parser <- add_argument(parser, "--label_column", help = "Column in reference .obs to use", default = "cell_type")
parser <- add_argument(parser, "--output",       help = "Output .h5ad path")
args <- parse_args(parser)

cat("Loading query:", args$query, "\n")
query_ad <- read_h5ad(args$query, to = "SingleCellExperiment")

cat("Loading reference:", args$reference, "\n")
ref_ad <- read_h5ad(args$reference, to = "SingleCellExperiment")

# Validate that the label column exists in the reference
if (!(args$label_column %in% colnames(colData(ref_ad)))) {
  stop(paste0(
    "Label column '", args$label_column, "' not found in reference .obs. ",
    "Available columns: ", paste(colnames(colData(ref_ad)), collapse = ", ")
  ))
}

ref_labels <- colData(ref_ad)[[args$label_column]]
cat("Reference cell types:", paste(unique(ref_labels), collapse = ", "), "\n")

cat("Running SingleR …\n")
pred <- SingleR(
  test      = query_ad,
  ref       = ref_ad,
  labels    = ref_labels,
  BPPARAM   = MulticoreParam(workers = max(1L, parallel::detectCores() - 1L))
)

# --- Write labels back to the h5ad ----------------------------------------
# Re-load as AnnData to modify .obs and write back
query_py <- read_h5ad(args$query)

query_py$obs[["singler_label"]]       <- pred$labels
query_py$obs[["singler_pruned_label"]] <- pred$pruned.labels
query_py$obs[["singler_score_delta"]] <- pred$delta.next

cat("SingleR label counts:\n")
print(table(pred$labels))

# Reproducibility info
query_py$uns[["pipeline_log"]][["label_transfer"]] <- list(
  tool           = "SingleR",
  reference_path = args$reference,
  label_column   = args$label_column,
  n_reference    = nrow(ref_ad),
  completed_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)

dir.create(dirname(args$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(query_py, args$output)
cat("Saved →", args$output, "\n")
