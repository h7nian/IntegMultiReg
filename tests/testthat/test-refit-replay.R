test_that("saved refit folds replay predictions and fitting seeds", {
  for (type in c("continuous", "binary", "right.censored")) {
    fit <- fit_demo(type, total = 15, burn = 5, seed = 53)
    original <- cv_imr(fit, cv_method = "refit", k = 3, rounds = 2)
    set.seed(817)
    state <- .Random.seed
    for (workers in c(1L, 2L)) {
      folds <- original$control$folds
      replay <- cv_imr(fit, cv_method = "refit",
        folds = folds[rev(seq_len(nrow(folds))), ], workers = workers)
      expect_identical(replay[1:5], original[1:5], info = type)
      expect_identical(replay$control$refit_seeds, original$control$refit_seeds)
      expect_identical(.Random.seed, state)
    }
    expect_identical(dim(original$control$refit_seeds), c(2L, 3L))
  }
})
