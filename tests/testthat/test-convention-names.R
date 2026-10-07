old_schema_fit <- function(fit, sampler = "paper", indexing = "standard") {
  fit$schema_version <- 2L
  names(fit$control)[names(fit$control) == "model_variant"] <- "method"
  names(fit$control)[names(fit$control) == "selection_update"] <- "sampler_method"
  fit$control$sampler_method <- sampler
  fit$control$numerical$prior_indexing <- indexing
  fit
}

test_that("explicit metadata conversion preserves stored draws", {
  current <- fit_demo("continuous", total = 20, burn = 10, seed = 91)
  saved <- old_schema_fit(current)
  expect_error(validate_imr_object(saved), "upgrade_imr_object")
  upgraded <- upgrade_imr_object(saved)
  expect_identical(upgraded$posterior, current$posterior)
  expect_identical(upgraded$control$model_variant, "imr")
  expect_identical(upgraded$control$selection_update, "symmetric_mrf_hastings")
  expect_identical(predict(upgraded, simIMR$platforms, covariates = simIMR$covariates),
                   predict(current, simIMR$platforms, covariates = simIMR$covariates))
  expect_identical(upgrade_imr_object(upgraded), upgraded)
  historical <- upgrade_imr_object(old_schema_fit(current, "legacy", "code2017"))
  expect_true(validate_imr_object(historical))
  expect_identical(historical$posterior, current$posterior)
  expect_identical(inclusion_probabilities(historical), inclusion_probabilities(current))
  expect_error(predict(historical, simIMR$platforms), "archived source")
  expect_error(cv_imr(historical), "archived source")
  expect_error(sample_regression_posterior(historical), "archived source")
})

test_that("retired values do not become silent aliases", {
  expect_error(cv_imr(fit_bin, cv_method = "importance"), "arg")
  expect_error(cv_imr(fit_bin, cv_method = "legacy"), "arg")
  expect_error(cv_imr(fit_bin, cv_method = "postfit_original"), "arg")
  expect_error(cv_imr(fit_bin, cv_method = "reweight", model_set = "draws"), "model_set")
  expect_error(cv_imr(fit_bin, cv_method = "reweight", df_method = "legacy_integer"), "df_method")
  expect_error(cv_imr(fit_bin, score_method = "legacy"), "unused argument")
  expect_error(imr(simIMR$platforms, simIMR$outcome.continuous,
                   outcome_type = "continuous", method = "imr"), "Unused argument")
})

test_that("coefficient extraction never returns inclusion probabilities or starts MCMC", {
  set.seed(213)
  rng <- .Random.seed
  expect_error(coef(fit_bin), "does not store regression coefficients")
  expect_identical(.Random.seed, rng)
  expect_true(all(vapply(inclusion_probabilities(fit_bin), function(m) all(m >= 0 & m <= 1), TRUE)))
  expect_error(confint(fit_bin), "Specify.*parm")
})

test_that("prediction labels identify subgroups and their platforms", {
  p <- predict(fit_bin, simIMR$platforms, covariates = simIMR$covariates)
  expect_s3_class(p, "imr_predictions")
  expect_named(p, paste0("subgroup:", fit_bin$model$subgroup_names))
  text <- capture.output(visible <- withVisible(print(p)))
  expect_false(visible$visible)
  expect_identical(visible$value, p)
  expect_true(any(grepl("genomic", text, fixed = TRUE)))
  expect_error(predict(fit_bin, simIMR$platforms, type = "response"), "Unused argument")
})
