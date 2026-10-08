test_that("outcome dimension is validated per outcome type", {
  bad <- data.frame(id = simIMR$outcome.binary$id, y = simIMR$outcome.binary$y, extra = 1)
  expect_error(imr(simIMR$platforms, bad, outcome_type = "binary", min_subgroup_size = 30, mcmc = imr_mcmc(
    draws = 50,
    burnin = 25, chains = 2L, diagnostics = FALSE
  )), "id")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "right.censored", min_subgroup_size = 30,
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "right-censored")
})

test_that("retained-draw and burn-in counts are validated", {
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", min_subgroup_size = 30,
    mcmc = imr_mcmc(draws = 0, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "`draws`")
  expect_s3_class(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", min_subgroup_size = 30,
    mcmc = imr_mcmc(draws = 60, burnin = 80, seed = 1, chains = 2L, diagnostics = FALSE)
  ), "imr")
})

test_that("an min_subgroup_size that excludes every subgroup is an error", {
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", min_subgroup_size = 1e+06,
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "No availability subgroup")
})

test_that("outcome_type and model_variant are matched against their choices", {
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary, outcome_type = "poisson", mcmc = imr_mcmc(
    draws = 50,
    burnin = 25, chains = 2L, diagnostics = FALSE
  )), "should be one of")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", model_variant = "lasso",
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "should be one of")
})

test_that("downstream functions reject non-imr input", {
  expect_error(cv_imr(cv_method = "refit", list(1)), "imr")
  expect_error(IntegMultiReg:::predict.imr(list(1), newdata = list()), "imr")
})

test_that("input data frames have standard id and numeric-column validation", {
  no_id <- simIMR$platforms
  names(no_id[[1]])[1] <- "sample_id"
  expect_error(
    imr(no_id, simIMR$outcome.binary,
      covariates = simIMR$covariates, outcome_type = "binary",
      min_subgroup_size = 30, mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
    ),
    "id.*first column"
  )
  dup_id <- simIMR$platforms
  dup_id[[1]]$id[2] <- dup_id[[1]]$id[1]
  expect_error(
    imr(dup_id, simIMR$outcome.binary,
      covariates = simIMR$covariates, outcome_type = "binary",
      min_subgroup_size = 30, mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
    ),
    "unique subject identifiers"
  )
  non_numeric <- simIMR$platforms
  non_numeric[[1]]$G01 <- as.character(non_numeric[[1]]$G01)
  expect_error(
    imr(non_numeric, simIMR$outcome.binary,
      covariates = simIMR$covariates, outcome_type = "binary",
      min_subgroup_size = 30, mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
    ),
    "must be numeric"
  )
  bad_binary <- simIMR$outcome.binary
  bad_binary$y[1] <- 2
  expect_error(
    imr(simIMR$platforms, bad_binary,
      covariates = simIMR$covariates, outcome_type = "binary",
      min_subgroup_size = 30, mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
    ),
    "0/1"
  )
  bad_survival <- simIMR$outcome.survival
  bad_survival$time[1] <- 0
  expect_error(
    imr(simIMR$platforms, bad_survival,
      covariates = simIMR$covariates, outcome_type = "right.censored",
      min_subgroup_size = 30, mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
    ),
    "positive"
  )
})

test_that("scalar and hyper-parameter arguments are validated before sampling", {
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary, outcome_type = "binary", priors = imr_priors(nu = c(
    -3,
    -3
  )), mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)), "one per platform")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", priors = imr_priors(molecular_scale = -1),
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "`molecular_scale`")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary, outcome_type = "binary", mcmc = imr_mcmc(
    draws = 50.5,
    burnin = 25, chains = 2L, diagnostics = FALSE
  )), "`draws`")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", min_subgroup_size = -1,
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "`min_subgroup_size`")
  expect_error(imr(simIMR$platforms, simIMR$outcome.binary,
    outcome_type = "binary", verbose = NA,
    mcmc = imr_mcmc(draws = 50, burnin = 25, chains = 2L, diagnostics = FALSE)
  ), "`verbose`")
})
