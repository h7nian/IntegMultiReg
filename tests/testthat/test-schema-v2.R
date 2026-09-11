test_that("schema v2 has a fixed named structure", {
  expect_named(fit_bin, c(
    "schema_version", "control", "model", "preprocessing", "posterior"
  ))
  expect_identical(fit_bin$schema_version, 2L)
  expect_true(validate_imr(fit_bin))

  bad <- fit_bin
  bad$unexpected <- TRUE
  expect_error(validate_imr(bad), "unsupported top-level fields")

  bad <- fit_bin
  bad$model$subgroup_platforms[[1L]] <- 99L
  expect_error(validate_imr(bad), "platform indices")

  bad <- fit_bin
  storage.mode(bad$preprocessing$features[[1L]][[1L]]) <- "integer"
  expect_error(validate_imr(bad), "platform 1 is invalid")
})

legacy_fit_fixture <- function(f) {
  priors <- f$control$priors
  legacy <- list(
    gam_mean = f$posterior$inclusion_probabilities,
    theta_mean = f$posterior$interaction_means,
    estimate_latent_y = f$posterior$latent_response_mean,
    log_posterior = f$posterior$log_posterior,
    gam_sample = f$posterior$selection_draws,
    theta_sample = f$posterior$interaction_draws,
    list_hyperpara = c(priors$forced_scale, priors$molecular_scale,
      priors$residual, priors$interaction, f$control$seed, priors$nu),
    data1 = list(f$model$n_platforms,
      lapply(f$model$platform_subgroups, function(i) i - 1L),
      lapply(f$model$subgroup_platforms, function(i) i - 1L),
      length(f$model$subgroup_names), f$model$sample_sizes,
      lengths(f$model$feature_names), length(f$model$covariate_names),
      match(f$control$outcome_type,
            c("right.censored", "binary", "continuous")),
      f$control$mcmc$draws),
    data2 = list(f$preprocessing$features, f$preprocessing$response,
      f$preprocessing$covariates, f$preprocessing$feature_center,
      f$preprocessing$feature_scale, f$preprocessing$covariate_center,
      f$preprocessing$covariate_scale),
    call = f$control$call, type_outcome = f$control$outcome_type,
    response_scale = f$control$response_scale,
    method = toupper(f$control$method), n_platform = f$model$n_platforms,
    platform_names = f$model$platform_names,
    feature_names = f$model$feature_names,
    covariate_names = f$model$covariate_names,
    model_bitstrings = f$model$subgroup_names,
    sample_size = f$model$sample_sizes,
    model_platforms = f$model$subgroup_platforms,
    platform_models = f$model$platform_subgroups, nu = priors$nu,
    ssize = f$control$min_subgroup_size,
    sample_mcmc = c(total = f$control$mcmc$draws,
                    burnin = f$control$mcmc$burnin),
    input_data = f$preprocessing$input_data,
    formula = f$preprocessing$formula,
    formula_data = f$preprocessing$formula_data,
    terms = f$preprocessing$terms, contrasts = f$preprocessing$contrasts,
    xlevels = f$preprocessing$xlevels, formula_id = f$preprocessing$id
  )
  class(legacy) <- "imr"
  legacy
}

test_that("complete 0.1.x fits upgrade explicitly and damaged fits do not", {
  legacy <- legacy_fit_fixture(fit_bin)
  upgraded <- upgrade_imr_fit(legacy)
  expect_true(validate_imr(upgraded))
  expect_equal(upgraded$posterior, fit_bin$posterior)
  expect_equal(upgrade_imr_fit(upgraded), upgraded)

  legacy$data2 <- legacy$data2[-1L]
  expect_error(upgrade_imr_fit(legacy), "incomplete; refit")
})

test_that("MRF subgroup capacity is enforced at 16 before native code", {
  expect_true(IntegMultiReg:::.imr_check_mrf_capacity(list(seq_len(16L))))
  expect_error(
    IntegMultiReg:::.imr_check_mrf_capacity(list(seq_len(17L))),
    "at most 16"
  )
  bad <- fit_bin
  bad$model$platform_subgroups[[1L]] <- seq_len(17L)
  expect_error(validate_imr(bad), "at most 16")
})

test_that("extreme sparsity priors keep the MRF log normalizer finite", {
  ids <- seq_len(24L)
  platform <- data.frame(id = ids, marker = sin(ids))
  fit <- imr(
    list(assay = platform), data.frame(id = ids, y = cos(ids)),
    outcome_type = "continuous", nu = 1000,
    min_subgroup_size = 0, draws = 2, burnin = 0, seed = 91
  )
  expect_true(all(is.finite(fit$posterior$log_posterior)))
})
