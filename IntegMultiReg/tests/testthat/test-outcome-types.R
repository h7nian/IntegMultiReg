test_that("all three outcome types fit without error and return valid output", {
  for (type in c("binary", "continuous", "right.censored")) {
    fit <- fit_demo(type, total = 200, burn = 100, seed = 3)
    expect_s3_class(fit, "imr")
    expect_identical(fit$control$outcome_type, type)
    expect_true(all(is.finite(fit$posterior$log_posterior)))
    expect_true(all(vapply(coef(fit), function(m) all(m >= 0 & m <= 1), logical(1))))
  }
})

test_that("the BMS (non-integrative) method also runs", {
  fit <- imr(
    simIMR$platforms, simIMR$outcome.binary, covariates = simIMR$covariates,
    outcome_type = "binary", method = "bms",
    nu = c(-4, -3, -4), draws = 150, burnin = 75, min_subgroup_size = 30, seed = 1)
  expect_s3_class(fit, "imr")
  expect_identical(fit$control$method, "bms")
})

test_that("the model runs without clinical covariates", {
  fit <- imr(
    simIMR$platforms, simIMR$outcome.binary, covariates = NULL,
    outcome_type = "binary",
    nu = c(-4, -3, -4), draws = 150, burnin = 75, min_subgroup_size = 30, seed = 1)
  expect_s3_class(fit, "imr")
  expect_length(fit$model$covariate_names, 0L)
  expect_true(all(vapply(fit$preprocessing$covariates, ncol, integer(1L)) == 0L))
})

test_that("integer-valued outcomes are handled (coerced to double)", {
  oc <- simIMR$outcome.binary
  oc$y <- as.integer(oc$y)
  expect_silent(
    suppressWarnings(imr(
      simIMR$platforms, oc, covariates = simIMR$covariates, outcome_type = "binary",
      nu = c(-4, -3, -4), draws = 120, burnin = 60, min_subgroup_size = 30, seed = 1))
  )
})
