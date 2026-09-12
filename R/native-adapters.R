## The only R functions that know the native argument layout. Public methods
## validate and prepare named data, then hand it to these narrow adapters.

.imr_call_cv_postfit_native <- function(object, k, rounds, max_models,
                                        verbose, importance,
                                        cache_bytes = 128 * 1024^2,
                                        cache_hash_bits = 64L,
                                        stage = c("full", "plan", "predict", "score"),
                                        tasks = NULL, predictions = NULL) {
  # Internal controls exercise bounded-cache and collision paths in tests.
  # They never change which draws contribute to the statistical calculation.
  cache_bytes <- .imr_check_integer_scalar(cache_bytes, "cache_bytes", min = 0)
  cache_hash_bits <- .imr_check_integer_scalar(cache_hash_bits, "cache_hash_bits", min = 0)
  if (cache_bytes > 128 * 1024^2 || cache_hash_bits > 64L)
    .imr_abort("Invalid internal post-fit cache controls.")
  control <- object$control
  model <- object$model
  preprocessing <- object$preprocessing
  posterior <- object$posterior
  priors <- control$priors
  stage <- match.arg(stage)
  if (is.null(tasks)) tasks <- matrix(TRUE, k, rounds)
  if (!is.logical(tasks) || anyNA(tasks) ||
      !identical(dim(tasks), as.integer(c(k, rounds))) ||
      (stage != "predict" && !all(tasks)))
    .imr_abort("Invalid internal post-fit task mask.")
  if (stage == "score") {
    if (!is.double(predictions) || any(!is.finite(predictions)) ||
        !identical(dim(predictions), as.integer(c(sum(model$sample_sizes), rounds))))
      .imr_abort("Invalid internal post-fit scoring predictions.")
  } else if (!is.null(predictions)) {
    .imr_abort("Predictions may only be supplied to the internal scoring stage.")
  }
  .quietly(verbose, .Call(
    "imr_cv_postfit",
    as.double(priors$forced_scale), as.double(priors$molecular_scale),
    as.double(priors$residual[["shape"]]), as.double(priors$residual[["rate"]]),
    as.double(priors$interaction[["shape"]]), as.double(priors$interaction[["rate"]]),
    as.double(control$seed), as.double(priors$nu),
    posterior$latent_response_mean, posterior$selection_draws,
    posterior$interaction_means,
    as.integer(model$n_platforms),
    lapply(model$platform_subgroups, function(i) as.integer(i - 1L)),
    lapply(model$subgroup_platforms, function(i) as.integer(i - 1L)),
    as.integer(length(model$subgroup_names)), as.integer(model$sample_sizes),
    as.integer(lengths(model$feature_names)), as.integer(length(model$covariate_names)),
    preprocessing$features, preprocessing$response, preprocessing$covariates,
    as.integer(match(control$outcome_type,
                     c("right.censored", "binary", "continuous"))),
    as.integer(control$mcmc$draws), as.integer(k), as.integer(rounds),
    as.integer(max_models), importance, as.double(c(cache_bytes, cache_hash_bits)),
    list(match(stage, c("full", "plan", "predict", "score")) - 1L,
         tasks, predictions),
    PACKAGE = "IntegMultiReg"
  ))
}

.imr_call_fit_native <- function(priors, seed, nu, method, n_platforms,
                                 platform_subgroups, subgroup_platforms,
                                 sample_sizes, n_features, n_covariates,
                                 features, response, outcome_type, covariates,
                                 draws, burnin, verbose) {
  .quietly(verbose, .Call(
    "imr_fit",
    as.double(priors$forced_scale), as.double(priors$molecular_scale),
    as.double(priors$residual[["shape"]]),
    as.double(priors$residual[["rate"]]),
    as.double(priors$interaction[["shape"]]),
    as.double(priors$interaction[["rate"]]), as.double(seed), as.double(nu),
    toupper(method), as.integer(n_platforms), platform_subgroups,
    subgroup_platforms, as.integer(length(sample_sizes)),
    as.integer(sample_sizes), as.integer(n_features), as.integer(n_covariates),
    features, response,
    as.integer(match(outcome_type,
                     c("right.censored", "binary", "continuous"))),
    covariates, as.integer(draws), as.integer(burnin),
    PACKAGE = "IntegMultiReg"
  ))
}

.imr_call_predict_native <- function(control, model, posterior,
                                     features, covariates, test_features,
                                     test_covariates, test_sample_sizes,
                                     max_models, verbose) {
  priors <- control$priors
  .quietly(verbose, .Call(
    "imr_predict",
    as.double(priors$forced_scale), as.double(priors$molecular_scale),
    as.double(priors$residual[["shape"]]),
    as.double(priors$residual[["rate"]]),
    as.double(priors$interaction[["shape"]]),
    as.double(priors$interaction[["rate"]]), as.double(control$seed),
    as.double(priors$nu), posterior$latent_response_mean,
    posterior$selection_draws, posterior$interaction_means,
    toupper(control$method), as.integer(model$n_platforms),
    lapply(model$platform_subgroups, function(index) as.integer(index - 1L)),
    lapply(model$subgroup_platforms, function(index) as.integer(index - 1L)),
    as.integer(length(model$subgroup_names)), as.integer(model$sample_sizes),
    as.integer(lengths(model$feature_names)),
    as.integer(length(model$covariate_names)), features, covariates,
    as.integer(control$mcmc$draws), test_features, test_covariates,
    as.integer(test_sample_sizes), as.integer(max_models),
    as.integer(match(control$outcome_type,
                     c("right.censored", "binary", "continuous"))),
    PACKAGE = "IntegMultiReg"
  ))
}
