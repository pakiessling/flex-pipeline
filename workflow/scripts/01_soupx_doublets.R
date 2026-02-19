#!/usr/bin/env Rscript
# 01_soupx.R — Ambient RNA correction (SoupX) + doublet detection (scDblFinder).
#
# Pure R implementation. Loads filtered and raw 10X H5 files directly via
# Seurat::Read10X_h5, builds a SoupChannel, and runs autoEstCont/adjustCounts.
# Falls back to the unmodified filtered matrix on failure or when --skip is set.
# Runs scDblFinder on the (corrected) counts and stores scDblFinder.score and
# scDblFinder.class in obs. Writes output as .h5ad via anndataR.

if (!requireNamespace("anndataR", quietly = TRUE)) {
  pak::pak("scverse/anndataR@v0.1.0")
}

suppressPackageStartupMessages({
  library(SoupX)
  library(Matrix)
  library(anndataR)
  library(argparser)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
})

parser <- arg_parser("SoupX ambient RNA correction + scDblFinder doublet detection")
parser <- add_argument(parser, "--sample",                help = "Sample identifier")
parser <- add_argument(parser, "--output",                help = "Output .h5ad path")
parser <- add_argument(parser, "--filtered_matrix",       help = "Filtered 10x H5 matrix path")
parser <- add_argument(parser, "--raw_matrix",            help = "Raw 10x H5 matrix path", default = NA)
parser <- add_argument(parser, "--skip",                  help = "Skip SoupX correction", flag = TRUE)
parser <- add_argument(parser, "--expected_doublet_rate", help = "Expected doublet rate (dbr for scDblFinder)", default = NA)
args <- parse_args(parser)

if (!args$skip && (is.na(args$raw_matrix) || args$raw_matrix == "")) {
  stop("--raw_matrix is required unless --skip is set")
}

# ---------------------------------------------------------------------------
# Load filtered matrix
# ---------------------------------------------------------------------------
cat(sprintf("[%s] Loading filtered matrix: %s\n", args$sample, args$filtered_matrix))
filt_mat <- Seurat::Read10X_h5(args$filtered_matrix)
# FLEX outputs can return a named list of modalities; keep Gene Expression only
if (is.list(filt_mat)) {
  filt_mat <- filt_mat[["Gene Expression"]]
}

# ---------------------------------------------------------------------------
# Run SoupX (or skip)
# ---------------------------------------------------------------------------
soupx_run <- FALSE

if (args$skip) {
  cat(sprintf("[%s] SoupX skipped (--skip flag).\n", args$sample))
  corrected_mat <- filt_mat

} else {
  cat(sprintf("[%s] Loading raw matrix: %s\n", args$sample, args$raw_matrix))
  raw_mat <- Seurat::Read10X_h5(args$raw_matrix)
  if (is.list(raw_mat)) {
    raw_mat <- raw_mat[["Gene Expression"]]
  }

  # Align genes (raw may contain more features than filtered)
  common_genes <- intersect(rownames(filt_mat), rownames(raw_mat))
  filt_mat <- filt_mat[common_genes, ]
  raw_mat  <- raw_mat[common_genes, ]

  # autoEstCont requires cluster labels. Since we load H5 files directly
  # (bypassing load10X), we must supply them ourselves via Seurat.
  cat(sprintf("[%s] Clustering filtered matrix for SoupX …\n", args$sample))
  seu <- Seurat::CreateSeuratObject(counts = filt_mat)
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu), verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 15, verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:15, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, verbose = FALSE)
  soupx_clusters <- setNames(as.character(Seurat::Idents(seu)), colnames(filt_mat))
  rm(seu)

  result <- tryCatch({
    cat(sprintf("[%s] Building SoupChannel and running autoEstCont …\n", args$sample))
    sc <- SoupChannel(raw_mat, filt_mat)
    sc <- setClusters(sc, soupx_clusters)
    sc <- autoEstCont(sc, doPlot = FALSE)
    out <- adjustCounts(sc, roundToInt = TRUE)
    list(corrected = out, success = TRUE)
  }, error = function(e) {
    cat(sprintf("[%s] SoupX failed (%s); falling back to filtered matrix.\n",
                args$sample, conditionMessage(e)))
    list(corrected = filt_mat, success = FALSE)
  })

  corrected_mat <- result$corrected
  soupx_run     <- result$success

  if (soupx_run) {
    before_sum <- sum(filt_mat)
    after_sum  <- sum(corrected_mat)
    pct        <- (after_sum - before_sum) / before_sum * 100
    cat(sprintf("[%s] SoupX complete. Counts before=%.0f, after=%.0f (%+.2f%%)\n",
                args$sample, before_sum, after_sum, pct))
  }
}

# ---------------------------------------------------------------------------
# scDblFinder doublet detection
# ---------------------------------------------------------------------------
n_cells <- ncol(filt_mat)
n_genes <- nrow(filt_mat)
cat(sprintf("[%s] Cells: %d, Genes: %d\n", args$sample, n_cells, n_genes))

cat(sprintf("[%s] Running scDblFinder …\n", args$sample))
set.seed(0)

dbr_arg <- if (!is.na(args$expected_doublet_rate)) as.numeric(args$expected_doublet_rate) else NULL

dbl_df <- tryCatch({
  sce <- SingleCellExperiment(assays = list(counts = corrected_mat))
  sce <- scDblFinder(sce, dbr = dbr_arg, BPPARAM = SerialParam())
  data.frame(
    scDblFinder.score = sce$scDblFinder.score,
    scDblFinder.class = as.character(sce$scDblFinder.class),
    row.names         = colnames(sce),
    stringsAsFactors  = FALSE
  )
}, error = function(e) {
  cat(sprintf("[%s] scDblFinder failed (%s); filling with NA.\n",
              args$sample, conditionMessage(e)))
  data.frame(
    scDblFinder.score = rep(NA_real_,      n_cells),
    scDblFinder.class = rep(NA_character_, n_cells),
    row.names         = colnames(filt_mat),
    stringsAsFactors  = FALSE
  )
})

n_doublets <- sum(dbl_df$scDblFinder.class == "doublet", na.rm = TRUE)
cat(sprintf("[%s] scDblFinder: %d / %d doublets (%.1f%%)\n",
            args$sample, n_doublets, n_cells, 100 * n_doublets / n_cells))

# ---------------------------------------------------------------------------
# Build AnnData and write h5ad
# ---------------------------------------------------------------------------
obs_df <- data.frame(
  Sample                = rep(args$sample, n_cells),
  scDblFinder.score     = dbl_df[colnames(filt_mat), "scDblFinder.score"],
  scDblFinder.class     = dbl_df[colnames(filt_mat), "scDblFinder.class"],
  row.names             = colnames(filt_mat),
  stringsAsFactors      = FALSE
)
var_df <- data.frame(row.names = rownames(filt_mat),
                     stringsAsFactors = FALSE)

# AnnData expects cells × genes; Read10X_h5 returns genes × cells
layers <- list(b4_soupx = t(filt_mat))
if (soupx_run) {
  layers[["after_soupx"]] <- t(corrected_mat)
}

adata <- AnnData(
  X      = t(corrected_mat),
  obs    = obs_df,
  var    = var_df,
  layers = layers,
  uns    = list(
    pipeline_log = list(
      soupx = list(
        completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        soupx_run    = soupx_run,
        skipped      = args$skip,
        n_cells      = n_cells,
        n_genes      = n_genes,
        software     = list(
          SoupX    = as.character(packageVersion("SoupX")),
          anndataR = as.character(packageVersion("anndataR"))
        )
      ),
      scdblfinder = list(
        completed_at         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        n_doublets           = n_doublets,
        n_cells              = n_cells,
        doublet_rate         = n_doublets / n_cells,
        expected_doublet_rate = if (!is.na(args$expected_doublet_rate)) as.numeric(args$expected_doublet_rate) else NULL,
        software             = list(
          scDblFinder = as.character(packageVersion("scDblFinder"))
        )
      )
    )
  )
)

dir.create(dirname(args$output), showWarnings = FALSE, recursive = TRUE)
write_h5ad(adata, args$output)
cat(sprintf("[%s] Saved → %s (SoupX: %s)\n",
            args$sample, args$output,
            ifelse(soupx_run, "yes", "no/skipped")))
