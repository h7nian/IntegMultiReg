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

#' Original Point Predictions from a Laplace Selection Fit
#'
#' Preserves the original ranked-model calculation for a fit that marginalizes
#' both regression coefficients and residual variance. It uses approximate
#' coefficient modes and rescored distinct selection states, rather than MCMC
#' visit-frequency weights. This method returns point predictions without
#' additional posterior sampling.
#'
#' @section Prediction rule:
#' Candidate states are collected from the end of the stored history until
#' all draws have been inspected or `100 * max_models` distinct states have
#' been collected. The highest-scoring `max_models` states receive normalized
#' exponential weights. For transformed predictor row x, prediction averages
#' \deqn{\sum_m w_m x^T\widehat\beta_m.}{sum_m weight_m * x-transpose * approximate_beta_m.}
#' Binary predictions instead average
#' \deqn{\sum_m w_m\Phi(x^T\widehat\beta_m).}{sum_m weight_m * Phi(x-transpose * approximate_beta_m).}
#' This is the original unit-variance probit approximation. For a log-time
#' survival fit, the returned prediction is on the log-time scale. Exponentiating
#' it gives a transformed point prediction, not a posterior mean survival time.
#'
#' @param object An `imr_selection` fit from
#'   `imr(..., marginalize = "coefficients_and_variance")`.
#' @param newdata,platform_names,covariates New-subject inputs as in [predict.imr()].
#' @param max_models Maximum retained distinct models in the point approximation.
#' @param verbose Print the native model-scoring diagnostics.
#' @param ... Unused arguments are rejected.
#' @return An `imr_predictions` list. Binary predictions are probabilities;
#'   continuous and survival predictions use the working response scale.
#' @examples
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   marginalize = "coefficients_and_variance", min_subgroup_size = 0,
#'   priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 40, burnin = 20, chains = 1,
#'     initial = NULL, seed = 1, diagnostics = FALSE
#'   )
#' )
#' predict(fit, list(assay = x))
#' @export
predict.imr_selection <- function(object, newdata, platform_names = NULL,
                                  covariates = NULL, max_models = 100L,
                                  verbose = FALSE, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(object)
  .imr_check_flag(verbose, "verbose")
  max_models <- .imr_check_integer_scalar(max_models, "max_models", min = 1)
  inputs <- .imr_prediction_inputs(object, newdata, platform_names, covariates)
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
    as.integer(match(control$outcome_type, c("right.censored", "binary", "continuous"))),
    as.double(c(0, numerical$laplace_max_iter, numerical$laplace_tolerance)),
    PACKAGE = "IntegMultiReg"
  ))
  result <- Map(function(id, prediction) {
    data.frame(id = id, prediction = prediction, row.names = NULL, stringsAsFactors = FALSE)
  }, inputs$sample_ids, value)
  .imr_prediction_result(result, object)
}
