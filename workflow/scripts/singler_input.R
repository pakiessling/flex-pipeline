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

singler_workers <- function(requested, allocated = Sys.getenv("SLURM_CPUS_PER_TASK", "")) {
  positive_integer <- function(x) {
    n <- suppressWarnings(as.numeric(x))
    if (length(n) != 1L || !is.finite(n) || n < 1 || n != floor(n)) {
      stop("SingleR CPU limits must be positive integers")
    }
    as.integer(n)
  }
  workers <- positive_integer(requested)
  if (nzchar(allocated)) workers <- min(workers, positive_integer(allocated))
  workers
}
