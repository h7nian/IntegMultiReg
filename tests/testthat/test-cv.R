test_that("cv_imr returns labelled accuracy matrices", {
  cv <- cv_imr(cv_method = "refit", fit_bin, k = 5, rounds = 2)
  expect_named(cv, c("pooled", "fold_mean", "predictions", "metric", "validation", "control"))
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
  expect_error(cv_imr(cv_method = "refit", fit_bin, method = "lasso"), "unused argument")
  expect_error(cv_imr(cv_method = "refit", fit_bin, model_variant = "bms"), "unused argument")
})

test_that("CV retains the model variant of its fitted object", {
  fit_bms <- imr(
    simIMR$platforms, simIMR$outcome.binary, covariates = simIMR$covariates,
    outcome_type = "binary", model_variant = "bms", nu = c(-4, -3, -4),
    draws = 80, burnin = 40, min_subgroup_size = 30, seed = 14
  )
  cv <- cv_imr(cv_method = "refit", fit_bms, k = 5, rounds = 1)
  expect_identical(cv, cv_imr(fit_bms, k = 5, rounds = 1, cv_method = "refit"))
  expect_error(cv_imr(fit_bms, model_variant = "imr"), "unused argument")
  expect_identical(cv$metric, "AUC")
  expect_true(all(vapply(fit_bms$posterior$interaction_draws, is.null, logical(1))))
})

test_that("cross-validation results carry a class and print a summary", {
  f <- fit_bin
  cv <- cv_imr(f, k = 2, rounds = 1)
  expect_s3_class(cv, "imr_cv")
  # the object is still an ordinary list for existing code
  expect_true(is.list(cv))
  expect_false(is.null(cv$pooled))
  expect_false(is.null(cv$control$folds))

  out <- capture.output(print(cv))
  expect_true(any(grepl("IMR cross-validation", out)))
  expect_true(any(grepl(cv$metric, out, fixed = TRUE)))
  expect_true(any(grepl("1 round of 2-fold", out, fixed = TRUE)))
  expect_true(any(grepl("refit preprocessing", out)))
  expect_identical(cv$control$model_variant, "imr")
  capture.output(expect_invisible(print(cv)))

  refit <- cv_imr(f, k = 2, rounds = 1, cv_method = "refit")
  expect_s3_class(refit, "imr_cv")
  expect_false(any(grepl("States:", capture.output(print(refit)))))
  weighted <- cv_imr(f, k = 2, rounds = 1, cv_method = "reweight")
  expect_true(any(grepl("all retained draws", capture.output(print(weighted)))))
})

test_that("the printed scores keep significant digits, not decimal places", {
  cv <- cv_imr(fit_bin, k = 2, rounds = 1)
  cv$pooled[] <- 0.000123456
  out <- capture.output(print(cv))
  # round() would have shown 0 here; signif() keeps the leading digits
  expect_true(any(grepl("0.000123", out, fixed = TRUE)))
  expect_false(any(grepl("^\\s*\\[1,\\]\\s+0\\s*$", out)))
})
