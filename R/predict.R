#' Predict from a Joint IMR Posterior
#'
#' Integrates over the fitted joint coefficient and variance samples, routing
#' subjects by their measured platforms and reusing training transformations.
#' `quantity = "conditional_mean"` describes the conditional expected response;
#' `"new_observation"` additionally simulates outcome noise.
#'
#' @section Predictive quantities:
#' For a transformed design row x, each draw gives
#' \deqn{\eta^{(b,c)}=x^T\beta^{(b,c)}.}{eta = transpose(x) * beta for each joint draw.}
#' Continuous conditional means are eta, and binary probabilities are
#' \deqn{p^{(b,c)}=\Phi(\eta^{(b,c)}/\sqrt{v^{(b,c)}}).}{Binary probability = NormalCDF(eta / sqrt(variance)).}
#' For log-time survival, the draw-wise conditional mean time is
#' \deqn{\exp(\eta^{(b,c)}+v^{(b,c)}/2).}{Conditional mean survival time = exp(eta + variance/2).}
#' `new_observation` samples Gaussian working responses, then applies the binary
#' threshold or survival exponential where needed. Future censoring is not simulated.
#' Response-scale log-time survival summaries use posterior medians because an
#' inverse-gamma variance mixture need not have a finite time-scale mean.
#' Other point summaries are posterior means. Intervals use empirical inverse
#' CDF quantiles (type 1), preserving binary values and point masses.
#'
#' @param object A joint `imr` fit.
#' @param newdata An [imr_data()] object or list of platform data frames.
#' @param platform_names Optional training platform names or indices for that list.
#'   Named lists are matched to fitted names when possible.
#' @param covariates Clinical predictors; formula fits require the original
#'   formula columns so the training transformation can be reused.
#' @param quantity `"conditional_mean"` or `"new_observation"`.
#' @param type `"response"` returns the outcome scale; `"link"` returns the
#'   Gaussian working scale (latent utility for binary data, log time for AFT).
#' @param interval Include equal-tail intervals in the prediction tables.
#' @param level Equal-tail probability level.
#' @param seed Prediction-noise seed. It does not refit or resample parameters.
#' @param ... Unused arguments are rejected.
#' @return An `imr_predictions` list of data frames by availability subgroup,
#'   each containing `id`, `prediction` and optional `lower`/`upper` columns.
#' @export
predict.imr <- function(object, newdata, platform_names = NULL, covariates = NULL,
                        quantity = c("conditional_mean", "new_observation"),
                        type = c("response", "link"), interval = FALSE,
                        level = .95, seed = 1L, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(object)
  quantity <- match.arg(quantity)
  type <- match.arg(type)
  .imr_check_flag(interval, "interval")
  level <- .imr_check_level(level)
  seed <- .imr_check_integer_scalar(seed, "seed", min = 0)
  inputs <- .imr_prediction_inputs(object, newdata, platform_names, covariates)
  saved <- .imr_save_rng()
  on.exit(.imr_restore_rng(saved), add = TRUE)
  if (quantity == "new_observation") set.seed(seed)
  .imr_predict_joint(object, inputs, quantity, type, interval, level)
}

.imr_flat_draws <- function(x) {
  matrix(x,
    nrow = dim(x)[1L] * dim(x)[2L], ncol = dim(x)[3L],
    dimnames = list(NULL, dimnames(x)[[3L]])
  )
}


.imr_predict_joint <- function(object, inputs, quantity = "conditional_mean", type = "response",
                               interval = FALSE, level = .95, weights = NULL) {
  if (!is.null(inputs$empty)) {
    empty <- inputs$empty
    if (interval) {
      empty <- lapply(empty, function(x) {
        x$lower <- numeric()
        x$upper <- numeric()
        x
      })
    }
    return(.imr_prediction_result(empty, object))
  }
  total <- as.double(object$control$mcmc$draws) * object$control$mcmc$chains
  if (!is.null(weights) && (length(weights) != total || any(!is.finite(weights)) ||
    any(weights < 0) || sum(weights) <= 0)) {
    .imr_abort("Invalid posterior prediction weights.")
  }
  if (!is.null(weights)) weights <- weights / sum(weights)
  outcome <- object$control$outcome_type
  log_time <- outcome == "right.censored" && object$control$response_scale == "log" && type == "response"
  result <- lapply(seq_along(object$model$subgroup_names), function(g) {
    ids <- inputs$sample_ids[[g]]
    out <- data.frame(id = ids, prediction = rep(NA_real_, length(ids)))
    if (interval) {
      out$lower <- rep(NA_real_, length(ids))
      out$upper <- rep(NA_real_, length(ids))
    }
    if (!length(ids)) {
      return(out)
    }
    columns <- c(
      list(rep(1, length(ids)), inputs$cova_test[[g]]),
      inputs$x_test[[g]][object$model$subgroup_platforms[[g]]]
    )
    design <- do.call(cbind, columns)
    beta <- .imr_flat_draws(object$posterior$coefficients[[g]])
    variance <- as.vector(object$posterior$variance[, , g])
    # A Gaussian expected value is linear: avoid a subject-by-draw matrix when
    # only its mean is requested. Nonlinear links and intervals use batches.
    if (quantity == "conditional_mean" && !interval && (outcome == "continuous" || type == "link" ||
      (outcome == "right.censored" && !log_time))) {
      coefficients <- if (is.null(weights)) colMeans(beta) else drop(crossprod(weights, beta))
      out$prediction <- drop(design %*% coefficients)
      return(out)
    }
    batch_size <- max(1L, floor(1e6 / total))
    for (start in seq.int(1L, length(ids), by = batch_size)) {
      rows <- seq.int(start, min(length(ids), start + batch_size - 1L))
      values <- beta %*% t(design[rows, , drop = FALSE])
      if (quantity == "new_observation") {
        values <- values +
          matrix(stats::rnorm(length(values)), nrow = total) * sqrt(variance)
      }
      if (type == "response" && outcome == "binary") {
        if (quantity == "new_observation") {
          values[] <- as.numeric(values > 0)
        } else {
          values[] <- stats::pnorm(values / sqrt(variance))
        }
      } else if (log_time && quantity == "conditional_mean") values <- values + variance / 2
      if (any(!is.finite(values))) .imr_abort("Non-finite posterior predictions; inspect predictor scales and tails.")
      out$prediction[rows] <- if (log_time) apply(values, 2L, .imr_quantile, probability = .5, weights = weights) else if (is.null(weights)) colMeans(values) else drop(crossprod(weights, values))
      if (interval) {
        q <- vapply(seq_along(rows), function(j) {
          .imr_quantile(
            values[, j],
            c((1 - level) / 2, (1 + level) / 2), weights
          )
        }, numeric(2L))
        out$lower[rows] <- q[1L, ]
        out$upper[rows] <- q[2L, ]
      }
    }
    if (log_time) {
      for (name in setdiff(names(out), "id")) out[[name]] <- exp(out[[name]])
      if (any(!is.finite(out$prediction))) .imr_abort("Time-scale prediction overflows; inspect type = 'link' and posterior tails.")
    }
    out
  })
  .imr_prediction_result(result, object)
}

.imr_prediction_result <- function(result, fit) {
  names(result) <- paste0("subgroup:", fit$model$subgroup_names)
  attr(result, "platforms") <- lapply(fit$model$subgroup_platforms, function(index) {
    fit$model$platform_names[index]
  })
  class(result) <- c("imr_predictions", "list")
  result
}

#' Print Predictions by Availability Subgroup
#'
#' Labels each prediction table by its availability subgroup and the measured
#' platform names. List keys such as `subgroup:011` identify availability
#' patterns, not variable-selection models.
#'
#' @param x An `imr_predictions` object returned by [predict.imr()].
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
    require_rows = FALSE, require_features = FALSE
  )
  tt <- stats::delete.response(object$preprocessing$terms)
  variables <- all.vars(tt)
  if (!all(variables %in% names(covariates))) {
    # Preserve the component-wise interface for explicitly encoded matrices.
    if (identical(setdiff(names(covariates), id), object$model$covariate_names)) {
      names(covariates)[names(covariates) == id] <- "id"
      return(covariates[, c("id", object$model$covariate_names), drop = FALSE])
    }
    .imr_abort(sprintf(
      "`covariates` is missing formula variable(s): %s.",
      paste(setdiff(variables, names(covariates)), collapse = ", ")
    ))
  }
  mf <- stats::model.frame(tt,
    data = covariates, na.action = stats::na.fail,
    xlev = object$preprocessing$xlevels
  )
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
    supplied_names <- names(newdata)
    if (!is.null(supplied_names) && all(nzchar(supplied_names)) && all(supplied_names %in% model$platform_names)) {
      platform_names <- as.character(match(supplied_names, model$platform_names))
    } else {
      platform_names <- as.character(seq_along(newdata))
    }
  } else {
    if (is.character(platform_names) && all(platform_names %in% model$platform_names)) {
      platform_names <- as.character(match(platform_names, model$platform_names))
    }
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

  list(
    x_train = x_train, x_test = x_test, cova_test = cova_test,
    sample_ids = sample_ids, samplesize_test = samplesize_test,
    subgroup_names = subgroup_names, n_platforms = n_platforms
  )
}
