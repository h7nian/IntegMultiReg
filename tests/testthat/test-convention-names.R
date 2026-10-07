test_that("saved convention labels upgrade without changing numerical results", {
  args <- list(x = simIMR$platforms, outcome = simIMR$outcome.continuous,
               covariates = simIMR$covariates, outcome_type = "continuous",
               nu = c(-4, -3, -4), draws = 40, burnin = 20, seed = 91)
  for (sampler in c("corrected", "original")) {
    current <- do.call(imr, c(args, list(sampler_method = sampler,
                                       prior_indexing = "original")))
    saved <- current
    saved$control$sampler_method <- if (sampler == "corrected") "paper" else "legacy"
    saved$control$numerical$prior_indexing <- "code2017"
    expect_error(validate_imr(saved))
    upgraded <- upgrade_imr_fit(saved)
    expect_identical(upgraded, current)
    expect_identical(upgraded$posterior, saved$posterior)
    expect_identical(predict(upgraded, simIMR$platforms, covariates = simIMR$covariates),
                     predict(current, simIMR$platforms, covariates = simIMR$covariates))
    expect_identical(upgrade_imr_fit(upgraded), upgraded)
  }
})

test_that("retired argument values are rejected rather than treated as aliases", {
  args <- list(x = simIMR$platforms, outcome = simIMR$outcome.continuous,
               outcome_type = "continuous", draws = 10, burnin = 5)
  expect_error(do.call(imr, c(args, list(sampler_method = "paper"))), "arg")
  expect_error(do.call(imr, c(args, list(sampler_method = "legacy"))), "arg")
  expect_error(do.call(imr, c(args, list(prior_indexing = "code2017"))), "arg")
  expect_error(cv_imr(fit_bin, cv_method = "legacy"), "arg")
  expect_error(cv_imr(fit_bin, cv_method = "importance", df_method = "legacy_integer"),
               "df_method")
  expect_error(cv_imr(fit_bin, cv_method = "importance", score_method = "legacy"),
               "score_method")
})
