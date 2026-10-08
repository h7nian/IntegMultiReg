# Adapt the shared, validated preprocessing to the retained Laplace engine.
# Fixed choices here use the symmetric MRF, the boundary Hastings correction,
# and standard clinical/molecular prior indexing. Historical unadjusted updates
# are not another posterior sampler offered by the new interface.
.imr_call_collapsed <- function(model, prep, control, seed, initial = NULL,
                                verbose = FALSE,
                                numerical = NULL) {
  priors <- control$priors
  if (is.null(numerical)) {
    settings <- control$numerical %||% imr_control()
    numerical <- c(0, settings$laplace_max_iter, settings$laplace_tolerance)
  }
  .imr_quietly(verbose, .Call(
    "imr_collapsed_sample",
    as.double(priors$forced_scale), as.double(priors$molecular_scale),
    as.double(priors$residual[["shape"]]), as.double(priors$residual[["rate"]]),
    as.double(priors$interaction[["shape"]]), as.double(priors$interaction[["rate"]]),
    as.double(seed), as.double(priors$nu), toupper(control$model_variant),
    as.integer(model$n_platforms),
    lapply(model$platform_subgroups, function(i) as.integer(i - 1L)),
    lapply(model$subgroup_platforms, function(i) as.integer(i - 1L)),
    as.integer(length(model$subgroup_names)), as.integer(model$sample_sizes),
    as.integer(lengths(model$feature_names)), as.integer(length(model$covariate_names)),
    prep$features, prep$response,
    as.integer(match(control$outcome_type, c("right.censored", "binary", "continuous"))),
    prep$covariates, as.integer(control$mcmc$draws), as.integer(control$mcmc$burnin),
    1L, as.double(numerical), initial,
    as.integer(c(control$mcmc$thin, control$mcmc$keep_latent)),
    PACKAGE = "IntegMultiReg"
  ))
}
