.imr_selection_history <- function(object, indices = NULL) {
  draws <- object$control$mcmc$draws
  if (is.null(indices)) indices <- seq_len(draws * object$control$mcmc$chains)
  lapply(indices, function(index) {
    row <- (index - 1L) %% draws + 1L
    chain <- (index - 1L) %/% draws + 1L
    lapply(seq_len(object$model$n_platforms), function(l) {
      matrix(as.integer(object$posterior$selection[[l]][row, chain, ]),
        nrow = length(object$model$platform_subgroups[[l]]),
        ncol = length(object$model$feature_names[[l]]), byrow = TRUE
      )
    })
  })
}

.imr_predict_selection <- function(object, inputs, max_models, verbose, type) {
  if (!is.null(inputs$empty)) {
    return(.imr_prediction_result(inputs$empty, object))
  }
  control <- object$control
  model <- object$model
  post <- object$posterior
  priors <- control$priors
  numerical <- control$numerical %||% imr_control()
  history <- .imr_selection_history(object)
  value <- .imr_quietly(verbose, .Call("imr_predict_selection",
    as.double(priors$forced_scale), as.double(priors$molecular_scale),
    as.double(priors$residual[["shape"]]), as.double(priors$residual[["rate"]]),
    as.double(priors$interaction[["shape"]]), as.double(priors$interaction[["rate"]]),
    as.double(control$seed), as.double(priors$nu),
    post$latent_mean, history, post$interaction_mean, toupper(control$model_variant),
    as.integer(model$n_platforms),
    lapply(model$platform_subgroups, function(i) as.integer(i - 1L)),
    lapply(model$subgroup_platforms, function(i) as.integer(i - 1L)),
    as.integer(length(model$subgroup_names)), as.integer(model$sample_sizes),
    as.integer(lengths(model$feature_names)), as.integer(length(model$covariate_names)),
    inputs$x_train, object$preprocessing$covariates, as.integer(length(history)),
    inputs$x_test, inputs$cova_test, as.integer(inputs$samplesize_test),
    as.integer(max_models),
    as.integer(if (type == "link") 3L else match(control$outcome_type, c("right.censored", "binary", "continuous"))),
    as.double(c(numerical$laplace_max_iter, numerical$laplace_tolerance)),
    PACKAGE = "IntegMultiReg"
  ))
  result <- Map(function(id, prediction) {
    data.frame(id = id, prediction = prediction, row.names = NULL, stringsAsFactors = FALSE)
  }, inputs$sample_ids, value)
  .imr_prediction_result(result, object)
}
