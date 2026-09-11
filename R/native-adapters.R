## The only R functions that know the native argument layout. Public methods
## validate and prepare named data, then hand it to these narrow adapters.

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
