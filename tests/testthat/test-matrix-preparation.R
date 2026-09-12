test_that("matrix preparation exactly preserves historical moments and normalization", {
  # Independent test-only reference for the previous three-pass preparation.
  reference <- function(mat) {
    if (nrow(mat) == 0L || ncol(mat) == 0L) {
      return(list(mean = rep(0, ncol(mat)), sd = rep(1, ncol(mat)),
                   normalized = matrix(numeric(0), nrow(mat), ncol(mat))))
    }
    center <- apply(mat, 2, function(column) mean(column, na.rm = TRUE))
    scale <- apply(mat, 2, function(column) {
      value <- sd(column, na.rm = TRUE)
      if (is.na(value) || value == 0) 1 else value
    })
    normalized <- apply(mat, 2, function(column) {
      center <- mean(column, na.rm = TRUE)
      scale <- sd(column, na.rm = TRUE)
      if (is.na(scale) || scale == 0) rep(0, length(column)) else (column - center) / scale
    })
    if (is.null(dim(normalized))) normalized <- matrix(normalized, nrow(mat), ncol(mat))
    list(mean = center, sd = scale, normalized = normalized)
  }
  set.seed(909)
  for (index in seq_len(300L)) {
    n_rows <- sample(0:12, 1L)
    n_columns <- sample(0:8, 1L)
    mat <- matrix(rnorm(n_rows * n_columns), n_rows, n_columns)
    if (n_rows > 0L && n_columns > 0L && index %% 3L == 0L) mat[, 1L] <- 2
    if (n_rows > 0L && n_columns > 0L && index %% 5L == 0L) {
      mat[sample(length(mat), 1L)] <- sample(c(NA_real_, NaN, Inf, -Inf, 1e308, 1e-308), 1L)
    }
    if (index %% 2L == 0L) dimnames(mat) <- list(
      if (n_rows) paste0("r", seq_len(n_rows)) else character(),
      if (n_columns) paste0("c", seq_len(n_columns)) else character())
    original <- mat
    rng <- .Random.seed
    expect_identical(IntegMultiReg:::.imr_prepare_matrix(mat), reference(mat))
    expect_identical(mat, original)
    expect_identical(.Random.seed, rng)
  }
  for (mat in list(data.frame(a = 1:3, b = c(2, 1, 4)),
                   data.frame(a = 1:3, block = I(matrix(1:6, 3))))) {
    expect_identical(IntegMultiReg:::.imr_prepare_matrix(mat), reference(mat))
  }
})
