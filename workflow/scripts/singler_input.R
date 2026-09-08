# Explicit expression selection for SingleR. Inputs must be unscaled logcounts.
singler_expression <- function(adata, layer, role) {
  mat <- if (identical(layer, "X")) adata$X else adata$layers[[layer]]
  if (is.null(mat)) {
    stop(sprintf("%s is missing expression layer '%s'", role, layer))
  }
  values <- if (inherits(mat, "sparseMatrix")) mat@x else as.vector(mat)
  if (any(!is.finite(values)) || any(values < 0)) {
    stop(sprintf("%s layer '%s' must contain finite, nonnegative logcounts", role, layer))
  }
  genes <- rownames(adata$var)
  cells <- rownames(adata$obs)
  if (is.null(genes) || is.null(cells) || anyDuplicated(genes) || anyDuplicated(cells)) {
    stop(sprintf("%s must have unique gene and cell names", role))
  }
  if (!identical(dim(mat), c(length(cells), length(genes)))) {
    stop(sprintf("%s expression dimensions do not match obs/var", role))
  }
  mat <- t(mat)
  dimnames(mat) <- list(genes, cells)
  mat
}
