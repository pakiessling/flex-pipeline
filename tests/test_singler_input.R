library(Matrix)
source("workflow/scripts/singler_input.R")
x <- matrix(c(10, 20, 30, 40), nrow = 2)
a <- list(X = x, layers = list(logcounts = log1p(x)),
          obs = data.frame(row.names = c("cell2", "cell1")),
          var = data.frame(row.names = c("gene2", "gene1")))
out <- singler_expression(a, "logcounts", "Query")
stopifnot(all(out == t(log1p(x))), identical(colnames(out), c("cell2", "cell1")),
          identical(rownames(out), c("gene2", "gene1")))
a$layers$logcounts <- as(a$layers$logcounts, "sparseMatrix")
stopifnot(all(singler_expression(a, "logcounts", "Query") == out))
stopifnot(inherits(try(singler_expression(a, "missing", "Query"), silent = TRUE), "try-error"))
a$X[1, 1] <- -1
stopifnot(inherits(try(singler_expression(a, "X", "Reference"), silent = TRUE), "try-error"))
cat("SingleR expression selection tests passed\n")
