test_that("cv_imr returns labelled accuracy matrices", {
  cv <- cv_imr(cv_method = "refit", fit_bin, k = 5, rounds = 2)
  expect_named(cv, c("pooled", "fold_mean", "predictions", "metric", "validation"))
  expect_equal(ncol(cv$pooled), length(fit_bin$model$subgroup_names) + 1L)
  expect_equal(nrow(cv$pooled), 2L)
  expect_equal(tail(colnames(cv$pooled), 1), "all")
  expect_identical(cv$metric, "AUC")
  expect_true(all(is.finite(cv$pooled)))
})

test_that("the cross-validated metric reflects the outcome type", {
  fc <- fit_demo("continuous", total = 200, burn = 100, seed = 4)
  expect_identical(cv_imr(cv_method = "refit", fc, k = 5, rounds = 1)$metric, "MSE")
  fs <- fit_demo("right.censored", total = 200, burn = 100, seed = 4)
  expect_identical(cv_imr(cv_method = "refit", fs, k = 5, rounds = 1)$metric, "C-index")
})

test_that("survival cross-validated C-index beats chance", {
  fs <- fit_demo("right.censored", total = 800, burn = 300, seed = 8)
  cv <- cv_imr(cv_method = "refit", fs, k = 5, rounds = 3)
  expect_gt(mean(cv$pooled[, "all"]), 0.6)
})

test_that("cv_imr validates cross-validation controls", {
  expect_error(cv_imr(cv_method = "refit", fit_bin, k = 1), "`k`")
  expect_error(cv_imr(cv_method = "refit", fit_bin, k = max(fit_bin$model$sample_sizes) + 1), "sample size")
  expect_error(cv_imr(cv_method = "refit", fit_bin, rounds = 0), "`rounds`")
  expect_error(cv_imr(cv_method = "refit", fit_bin, max_models = 0), "`max_models`")
  expect_error(cv_imr(cv_method = "refit", fit_bin, method = "lasso"), "should be one of")
  expect_error(cv_imr(cv_method = "refit", fit_bin, method = "bms"), "must match the fitted object")
})

test_that("cv_imr accepts an explicitly matching fitted method", {
  fit_bms <- imr(
    simIMR$platforms, simIMR$outcome.binary, covariates = simIMR$covariates,
    outcome_type = "binary", method = "bms", nu = c(-4, -3, -4),
    draws = 80, burnin = 40, min_subgroup_size = 30, seed = 14
  )
  cv <- cv_imr(cv_method = "refit", fit_bms, k = 5, rounds = 1, method = "bms")
  expect_identical(cv$metric, "AUC")
  expect_true(all(vapply(fit_bms$posterior$interaction_draws, is.null, logical(1))))
})
