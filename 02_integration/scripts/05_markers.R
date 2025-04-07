if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# Check for required packages and install if needed
required_packages <- c("anndataR", "presto", "Matrix", "SingleCellExperiment", "argparser")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "anndataR") {
      pak::pak("scverse/anndataR")
    } else if (pkg == "presto") {
      pak::pak("immunogenomics/presto")
    } else {
      install.packages(pkg)
    }
  }
}





suppressPackageStartupMessages({
  library(anndataR)
  library(Matrix)
  library(presto)
  library(SingleCellExperiment)
  library(argparser, quietly = TRUE)
})

parser <- arg_parser("Find markers for each cluster with presto")
parser <- add_argument(parser, "--input", help = "Path to h5ad")
parser <- add_argument(parser, "--layer", help = "Which matrix to use", default = "X")
parser <- add_argument(parser, "--ct_column",
  help = "Which cell type column to use",
  default = "leiden_3"
)
parser <- add_argument(parser, "--output", help = "Path to output folder for files")
args <- parse_args(parser)


sce <- read_h5ad(args$input, to = "SingleCellExperiment")

print(sce)

print(assay(sce))
# for whatever reason Presto needs an assay named counts or logcounts
assay(sce, "logcounts") <- as(assay(sce, args$layer), "CsparseMatrix")

cat("Running Wilcoxon rank-sum test (wilcoxauc) for cluster markers.\n")
res <- wilcoxauc(sce, group_by = args$ct_column, assay = "logcounts")

top <- top_markers(res, n = 100, auc_min = .5, pct_in_min = 1, padj_max = 0.05)

cat("Writing ranked genes")
all_genes_file <- file.path(args$output, paste0("all_genes_ranked_", args$ct_column, ".csv"))
marker_genes_file <- file.path(args$output, paste0("marker_genes_ranked_", args$ct_column, ".csv"))
write.csv(res, file = all_genes_file, row.names = FALSE)
write.csv(top, file = marker_genes_file, row.names = FALSE)
