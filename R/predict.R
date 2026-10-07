#' Predict Outcomes for New Subjects
#'
#' @description
#' `predict()` method for objects of class `"imr"` produced by [imr()].  New
#' subjects are matched to the availability subgroup models learned during
#' training, missing platforms are handled automatically, and the
#' test features are standardized with the training set's centring and scaling
#' factors before Bayesian model averaging.
#'
#' @section Ranked-model point prediction:
#' The retained joint selection states are scanned in reverse iteration order
#' until the chain is exhausted or `100 * max_models` distinct states have
#' been collected. These states are ranked by recomputed approximate
#' log-posterior scores. The calculation
#' uses the full-fit latent-response means and mean MRF interactions. Let
#' \eqn{\mathcal M}{M} be the retained set of at most `max_models` states and
#' \eqn{\ell_m}{ell_m} their scores. The prediction weights are
#' \deqn{w_m=\frac{\exp(\ell_m-\ell_{\max})}
#'                  {\sum_{h\in\mathcal M}\exp(\ell_h-\ell_{\max})}.}{w_m = exp(ell_m - ell_max) / sum_h exp(ell_h - ell_max), over retained models.}
#' They are normalized scores of the ranked distinct states, rather than their
#' empirical MCMC visit frequencies. Let \eqn{\widehat b_{s,m}}{approximate_beta_s,m} be the
#' approximate coefficient mode from the iteration-limited coordinate-ascent
#' calculation, and let
#' \eqn{z} be the new subject's transformed design row. Then
#' \deqn{\widehat y=\sum_{m\in\mathcal M}w_m z^T\widehat b_{s,m}}{predicted working response = sum_m w_m * z-transpose * approximate_beta_s,m.}
#' for continuous and working-scale survival responses. Binary predictions
#' apply the probit link within each model:
#' \deqn{\widehat p=\sum_{m\in\mathcal M}w_m
#'       \Phi(z^T\widehat b_{s,m}).}{predicted event probability = sum_m w_m * Phi(z-transpose * approximate_beta_s,m).}
#' The binary point-prediction path uses the unit-variance probit approximation;
#' it does not average the conditional variance draws from [sample_regression_posterior()].
#' For a default survival fit, \eqn{\widehat y}{predicted working response} is on the log-time scale.
#'
#' These are plug-in point predictions. For coefficient uncertainty and
#' predictive intervals, use [sample_regression_posterior()] followed by
#' [predict.imr_posterior()], which also uses a different model-averaging
#' construction. Increasing `max_models` changes this finite ranked-model
#' approximation; it does not lengthen the fitted MCMC chain.
#'
#' @param object A fitted object of class `"imr"` returned by [imr()].
#' @param newdata A list of data frames with the new platform measurements, or
#'   an [imr_data()] object created for prediction. When an `imr_data` object is
#'   supplied, its covariates are used automatically.
#'   Each data frame must include `id` as the first column, followed by finite
#'   numeric feature columns matching the corresponding training platform.
#' @param platform_names A character vector giving, for each element of
#'   `newdata`, the index (`"1"`, `"2"`, ...) of the corresponding training
#'   platform.  Defaults to `NULL`, meaning the elements of `newdata` are taken
#'   to be in the same order as the platforms supplied to [imr()].
#' @param covariates An optional data frame of clinical covariates for the test
#'   subjects, including `id` as the first column.  Required when the model was
#'   fitted with covariates and ignored with a warning when it was not.
#'   For a formula fit, supply the original predictor columns (including factors
#'   and variables used in transformations). The training terms, factor levels
#'   and contrasts are reused. The training identifier name is also accepted.
#'   Alternatively, an already encoded numeric model matrix may be supplied as
#'   a data frame with `id` and exactly the fitted covariate column names.
#' @param max_models Integer; the maximum number of distinct selection models
#'   (gamma configurations) used for Bayesian model averaging.  Default `100`.
#' @param verbose Logical; if `TRUE`, print the C routine's diagnostics.
#'   Defaults to `FALSE`.
#' @param ... Unused arguments are rejected.
#'
#' @details
#' Predictions are only produced for subjects observed on at least one platform
#' (and, when `covariates` is supplied, with covariate data).  For `"binary"`
#' outcomes the returned `prediction` column is a probability obtained through the
#' probit link (`pnorm`) within each model before averaging; for `"continuous"` and `"right.censored"` outcomes it
#' is the predicted working response. For default log-time survival fits it is
#' on the log-time scale; exponentiating gives a transformed point prediction,
#' not a posterior mean survival time. Old and identity-scale fits retain their
#' historical response scale.
#'
#' @return A named list with one data frame per active availability subgroup
#'   model. Each data frame has columns `id` (subject identifier) and `prediction`
#'   (predicted value or probability at full numeric precision).
#'
#' @seealso [imr()], [cv_imr()]
#'
#' @examples
#' \donttest{
#' data("simIMR", package = "IntegMultiReg")
#' fit <- imr(
#'   x = simIMR$platforms, outcome = simIMR$outcome,
#'   covariates = simIMR$covariates, outcome_type = "binary",
#'   nu = c(-4, -3, -4), draws = 200, burnin = 100,
#'   min_subgroup_size = 5, seed = 1
#' )
#' new_x <- simIMR$platforms[[1]][1:10, ]
#' new_z <- simIMR$platforms[[2]][1:7, ]
#' predict(fit, newdata = list(new_x, new_z), covariates = simIMR$covariates)
#' }
#' @export
predict.imr <- function(object, newdata, platform_names = NULL,
                        covariates = NULL, max_models = 100,
                        verbose = FALSE, ...) {
  .imr_reject_dots(...)
  validate_imr_object(object)
  .imr_require_current_updates(object)
  .imr_check_flag(verbose, "verbose")
  max_models <- .imr_check_integer_scalar(max_models, "max_models", min = 1)
  inputs <- .imr_prediction_inputs(object, newdata, platform_names, covariates)
  if (!is.null(inputs$empty)) return(.imr_prediction_result(inputs$empty, object))
  x_train <- inputs$x_train
  x_test <- inputs$x_test
  cova_test <- inputs$cova_test
  sample_ids <- inputs$sample_ids
  samplesize_test <- inputs$samplesize_test
  subgroup_names <- inputs$subgroup_names
  n_platforms <- inputs$n_platforms
  control <- object$control
  model <- object$model
  prep <- object$preprocessing
  posterior <- object$posterior
  results <- .imr_call_predict_native(
    control = control, model = model, posterior = posterior,
    features = x_train, covariates = prep$covariates,
    test_features = x_test, test_covariates = cova_test,
    test_sample_sizes = samplesize_test, max_models = max_models,
    verbose = verbose
  )
  names(results) <- subgroup_names
  res <- mapply(function(x, y) {
    data.frame(id = x, prediction = y, row.names = NULL, stringsAsFactors = FALSE)
  }, sample_ids, results, SIMPLIFY = FALSE)

  .imr_prediction_result(res, object)
}

.imr_prediction_result <- function(result, fit) {
  names(result) <- paste0("subgroup:", fit$model$subgroup_names)
  attr(result, "platforms") <- lapply(fit$model$subgroup_platforms, function(index)
    fit$model$platform_names[index])
  class(result) <- c("imr_predictions", "list")
  result
}

#' Print Predictions by Availability Subgroup
#'
#' Labels each prediction table by its availability subgroup and the measured
#' platform names. List keys such as `subgroup:011` identify availability
#' patterns, not variable-selection models.
#'
#' @param x An `imr_predictions` object returned by [predict.imr()] or
#'   [predict.imr_posterior()].
#' @param row.names Whether to print row names in the subgroup tables.
#' @param ... Arguments passed to the data-frame printing method.
#' @return `x`, invisibly.
#' @export
print.imr_predictions <- function(x, row.names = FALSE, ...) {
  platforms <- attr(x, "platforms")
  for (g in seq_along(x)) {
    cat(sprintf("%s (%s)\n", names(x)[g], paste(platforms[[g]], collapse = " + ")))
    print(x[[g]], row.names = row.names, ...)
    if (g < length(x)) cat("\n")
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.imr_prediction_covariates <- function(object, covariates) {
  if (!is.data.frame(covariates)) {
    .imr_abort("`covariates` must be a data frame.")
  }
  id <- if ("id" %in% names(covariates)) "id" else object$preprocessing$id
  if (is.null(id) || !id %in% names(covariates)) {
    .imr_abort("`covariates` must contain the subject identifier column.")
  }
  ids <- covariates[[id]]
  .imr_check_id_frame(data.frame(id = ids), "covariates",
                      require_rows = FALSE, require_features = FALSE)
  tt <- stats::delete.response(object$preprocessing$terms)
  variables <- all.vars(tt)
  if (!all(variables %in% names(covariates))) {
    # Preserve the component-wise interface for explicitly encoded matrices.
    if (identical(setdiff(names(covariates), id), object$model$covariate_names)) {
      names(covariates)[names(covariates) == id] <- "id"
      return(covariates[, c("id", object$model$covariate_names), drop = FALSE])
    }
    .imr_abort(sprintf("`covariates` is missing formula variable(s): %s.",
                       paste(setdiff(variables, names(covariates)), collapse = ", ")))
  }
  mf <- stats::model.frame(tt, data = covariates, na.action = stats::na.fail,
                            xlev = object$preprocessing$xlevels)
  mm <- stats::model.matrix(tt, mf, contrasts.arg = object$preprocessing$contrasts)
  mm <- mm[, attr(mm, "assign") != 0L, drop = FALSE]
  if (nrow(mm) != length(ids) ||
      !identical(colnames(mm), object$model$covariate_names)) {
    .imr_abort("Formula covariates do not match the training model matrix.")
  }
  data.frame(id = ids, mm, check.names = FALSE, row.names = NULL)
}

# Shared routing/standardization for point and posterior predictions.
.imr_prediction_inputs <- function(object, newdata, platform_names, covariates) {
  if (!inherits(object, "imr")) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  model <- object$model
  prep <- object$preprocessing
  n_platforms <- model$n_platforms
  n_covariates <- length(model$covariate_names)

  if (inherits(newdata, "imr_data")) {
    validate_imr_data(newdata)
    if (!is.null(covariates)) {
      .imr_abort(
        "Supply covariates either inside `newdata` or through `covariates`, not both."
      )
    }
    covariates <- newdata$covariates
    if (is.null(platform_names)) {
      matched <- match(names(newdata$platforms), model$platform_names)
      if (anyNA(matched) || anyDuplicated(matched)) {
        .imr_abort(paste0(
          "Platform names in `newdata` must uniquely match the fitted model: ",
          paste(model$platform_names, collapse = ", "), "."
        ))
      }
      platform_names <- as.character(matched)
    }
    newdata <- newdata$platforms
  }

  if (!is.list(newdata) || length(newdata) == 0L) {
    .imr_abort("`newdata` must be a non-empty list of data frames.")
  }
  if (length(newdata) > n_platforms) {
    .imr_abort("`newdata` cannot contain more platforms than the fitted model.")
  }

  ## Default: the new platforms are supplied in the same order as at training.
  if (is.null(platform_names)) {
    platform_names <- as.character(seq_along(newdata))
  } else {
    if (length(platform_names) != length(newdata)) {
      .imr_abort("`platform_names` must have the same length as `newdata`.")
    }
    platform_names <- as.character(platform_names)
    if (anyNA(platform_names) || anyDuplicated(platform_names)) {
      .imr_abort("`platform_names` must not contain missing or duplicated values.")
    }
    valid_platforms <- as.character(seq_len(n_platforms))
    if (any(!platform_names %in% valid_platforms)) {
      .imr_abort(sprintf(
        "`platform_names` must use training platform indices in {%s}.",
        paste(valid_platforms, collapse = ", ")
      ))
    }
  }
  names(newdata) <- platform_names
  subgroup_names <- model$subgroup_names

  for (i in seq_along(newdata)) {
    arg <- sprintf("newdata[[%d]]", i)
    platform_index <- as.integer(platform_names[i])
    .imr_check_id_frame(newdata[[i]], arg, require_rows = FALSE)
    .imr_check_numeric_columns(newdata[[i]], arg)
    expected_n <- length(model$feature_names[[platform_index]])
    if (ncol(newdata[[i]]) - 1L != expected_n) {
      .imr_abort(sprintf(
        "`%s` must contain %d feature column(s) for training platform %s.",
        arg, expected_n, platform_names[i]
      ))
    }
    expected_names <- model$feature_names[[platform_index]]
    if (!is.null(expected_names) &&
        !identical(colnames(newdata[[i]])[-1], expected_names)) {
      .imr_abort(sprintf(
        "Feature columns in `%s` must match the training feature names.",
        arg
      ))
    }
  }

  if (n_covariates > 0L && is.null(covariates)) {
    .imr_abort("`covariates` is required because the model was fitted with covariates.")
  }
  if (n_covariates == 0L && !is.null(covariates)) {
    .imr_warn(
      "`covariates` was supplied, but the model was fitted without covariates; ignoring it."
    )
    covariates <- NULL
  }
  if (!is.null(covariates)) {
    if (!is.null(prep$terms)) {
      covariates <- .imr_prediction_covariates(object, covariates)
    }
    .imr_check_id_frame(covariates, "covariates", require_rows = FALSE)
    .imr_check_numeric_columns(covariates, "covariates")
    if (ncol(covariates) - 1L != n_covariates) {
      .imr_abort(sprintf(
        "`covariates` must contain %d covariate column(s).", n_covariates
      ))
    }
    expected_covariates <- model$covariate_names
    if (!is.null(expected_covariates) && length(expected_covariates) > 0L &&
        !identical(colnames(covariates)[-1], expected_covariates)) {
      .imr_abort("Columns in `covariates` must match the training covariate names.")
    }
  }

  all_ids <- unique(unlist(lapply(newdata, function(df) df$id)))

  ## We only predict subjects with covariate data that are observed in at least one omics platform
  if (!is.null(covariates)) {
    all_ids <- intersect(all_ids, covariates$id)
  }

  if (length(all_ids) == 0) {
    .imr_warn(
      "No subjects have both the required covariates and at least one platform."
    )
    return(list(empty = .imr_empty_predictions(subgroup_names)))
  }
  # Rows correspond to subjects and columns correspond to platforms.
  presence <- data.frame(do.call(cbind, lapply(newdata, function(df) {
    all_ids %in% df$id
  })))
  names(presence) <- platform_names
  not_active_platform <- setdiff(as.character(seq_len(n_platforms)), names(presence))
  for (l in not_active_platform) {
    presence[[l]] <- rep(FALSE, length(all_ids))
  }
  presence <- presence[, as.character(seq_len(n_platforms)), drop = FALSE]

  bitstrings <- apply(
    presence, 1,
    function(x) paste(as.integer(rev(x)), collapse = "")
  )

  x_train <- prep$features
  unique_patterns <- subgroup_names
  platforms <- newdata
  for (l in not_active_platform) {
    platforms[[l]] <- data.frame(matrix(nrow = 0, ncol = ncol(x_train[[1]][[as.numeric(l)]])))
  }

  platforms <- platforms[as.character(seq_len(n_platforms))]
  x_test <- list()
  sample_ids <- list()

  for (pat in unique_patterns) {
    subgroup_ids <- all_ids[bitstrings == pat]
    subgroup_list <- vector("list", n_platforms)
    names(subgroup_list) <- paste0("platform", seq_len(n_platforms))
    for (i in seq_len(n_platforms)) {
      bit <- substr(pat, n_platforms - i + 1, n_platforms - i + 1)
      if (bit == "1") {
        subgroup_list[[i]] <- .imr_match_rows(
          platforms[[i]], subgroup_ids, sprintf("newdata platform %d", i)
        )
      } else {
        subgroup_list[[i]] <- platforms[[i]][FALSE, , drop = FALSE]
      }
    }
    sample_ids[[pat]] <- subgroup_ids
    x_test[[pat]] <- subgroup_list
  }

  if (length(unlist(sample_ids)) == 0) {
    .imr_warn(
      "No new subjects belong to availability subgroup models retained during training."
    )
    return(list(empty = .imr_empty_predictions(subgroup_names)))
  }
  routed_ids <- unique(unlist(sample_ids, use.names = FALSE))
  dropped_ids <- setdiff(as.character(all_ids), as.character(routed_ids))
  if (length(dropped_ids) > 0L) {
    .imr_warn(sprintf(
      "%d subject(s) do not belong to availability subgroup models retained during training and were not predicted.",
      length(dropped_ids)
    ))
  }

  ### Create a vector that shows the presence of models in the test data
  samplesize_test <- unlist(lapply(sample_ids, length))
  if (!is.null(covariates)) {
    cova_test <- lapply(sample_ids, function(x) {
      .imr_match_rows(covariates, x, "covariates")
    })
  } else {
    cova_test <- lapply(sample_ids, function(x) matrix(numeric(0), nrow = length(x), ncol = 0))
  }
  x_test <- .imr_platform_matrices(x_test)

  if (!is.null(covariates)) {
    cova_test <- lapply(
      seq_along(cova_test),
      function(i) as.matrix(cova_test[[i]][, -1, drop = FALSE])
    )
  }

  norm_mean <- prep$feature_center
  norm_sd <- prep$feature_scale

  x_test <- .imr_standardize_platforms(x_test, norm_mean, norm_sd)

  mean_cov <- prep$covariate_center
  sd_cov <- prep$covariate_scale

  if (!is.null(covariates)) {
    cova_test <- mapply(.imr_standardize_matrix,
      cova_test, mean_cov, sd_cov,
      SIMPLIFY = FALSE
    )
  }

  list(x_train = x_train, x_test = x_test, cova_test = cova_test,
       sample_ids = sample_ids, samplesize_test = samplesize_test,
       subgroup_names = subgroup_names, n_platforms = n_platforms)
}
