#' Cross-Validated Predictive Accuracy of an IMR Fit
#'
#' @description
#' `cv_imr()` evaluates predictive accuracy of a model fitted with [imr()] using
#' repeated \eqn{K}-fold splits. Each training fold independently standardizes
#' predictors, runs the MCMC sampler and estimates model-averaging weights.
#' The original fit supplies model settings and raw data, not posterior draws.
#' Fit a separate object with `method = "BMS"` for a no-borrowing comparison.
#' The accuracy measure depends on the outcome type: the concordance index
#' (C-index) for right-censored outcomes, the area under the ROC curve (AUC) for
#' binary outcomes, and the mean squared error (MSE) for continuous outcomes.
#'
#' @param object A fitted object of class `"imr"` returned by [imr()].
#' @param k Integer number of cross-validation folds per round (default `5`).
#'   Must be at least `2`, and each availability subgroup must contain at least
#'   `k` subjects.
#' @param rounds Integer number of independent cross-validation rounds to
#'   average over (default `2`).  Must be positive.
#' @param method Optional compatibility argument.  If supplied, it must match
#'   the method stored in `object`; `cv_imr()` cannot turn an IMR fit into a BMS
#'   fit or vice versa.
#' @param max_models Integer maximum number of selection models used for
#'   Bayesian model averaging (default `100`).  Must be positive.
#' @param verbose Logical; if `TRUE`, print fold progress and sampler
#'   diagnostics.  Defaults to `FALSE`.
#'
#' @return A list with two numeric matrices, each of dimension
#'   `rounds` x `(n_subgroups + 1)`, whose last column corresponds to all
#'   subjects pooled and whose remaining columns are named by the availability
#'   subgroup bitstrings:
#'   \item{total_cindex}{Accuracy computed on the pooled cross-validated
#'     predictions within each availability subgroup (and overall).}
#'   \item{subset_cindex}{Accuracy computed fold-by-fold within each
#'     availability subgroup.}
#'   The accuracy metric (`"C-index"`, `"AUC"` or `"MSE"`) is recorded in the
#'   `"metric"` attribute of the returned list. The `"predictions"` attribute
#'   contains subject IDs, subgroups, rounds, folds and out-of-fold predictions.
#'   Undefined AUCs (single class) and C-indices (no comparable pairs) are `NA`.
#'
#' @details
#' Only subgroups retained in `object` are evaluated. Their training subsets
#' remain included even when they are smaller than the original `ssize`.
#' Formula transformations are rebuilt using the training rows. Raw formula
#' data must be stored in the fit; older formula fits must be refitted.
#' Each fold uses the original prior and MCMC settings, so runtime is roughly
#' `k * rounds` full fits. Folds and sampler seeds are reproducible from the
#' fitted seed, and the caller's R RNG state is restored. Hyperparameters are
#' treated as fixed; tuning them requires an additional validation layer.
#'
#' @seealso [imr()], [predict.imr()]
#'
#' @examples
#' \donttest{
#' data("simIMR", package = "IntegMultiReg")
#' fit <- imr(
#'   platform_data_list = simIMR$platforms, outcome = simIMR$outcome,
#'   cov = simIMR$covariates, type_outcome = "binary",
#'   nu = c(-4, -3, -4), sample_mcmc = c(200, 100), ssize = 5, seed = 1
#' )
#' cv <- cv_imr(fit, k = 5, rounds = 2)
#' cv$total_cindex
#' }
#' @export
cv_imr <- function(object, k = 5, rounds = 2,
                   method = NULL,
                   max_models = 100, verbose = FALSE) {
  if (!inherits(object, "imr")) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  .imr_check_flag(verbose, "verbose")
  k <- .imr_check_integer_scalar(k, "k", min = 2)
  rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  max_models <- .imr_check_integer_scalar(max_models, "max_models", min = 1)
  object_method <- object$method
  if (is.null(object_method)) object_method <- "IMR"
  if (is.null(method)) {
    method <- object_method
  } else {
    method <- match.arg(method, c("IMR", "BMS"))
    if (!identical(method, object_method)) {
      .imr_abort(
        "`method` must match the fitted object; fit a separate `imr(..., method = \"BMS\")` object for BMS validation."
      )
    }
  }
  if (any(object$sample_size < k)) {
    .imr_abort(
      "`k` must not exceed the sample size of any modelled availability subgroup."
    )
  }
  if (is.null(object$input_data)) {
    .imr_abort("Refit this model to retain the raw inputs required for cross-validation.")
  }
  validate_imr_data(object$input_data)
  if (object$type_outcome == "right.censored" && is.null(object$response_scale)) {
    .imr_abort("Refit this survival model with an explicit `survival_scale`.")
  }
  if (!is.null(object$terms) && is.null(object$formula_data)) {
    .imr_abort("Refit this formula model to retain the raw formula data required for cross-validation.")
  }
  dat <- object$input_data
  subjects <- dat$availability[
    dat$availability$subgroup %in% object$model_bitstrings, c("id", "subgroup"), drop = FALSE]
  outcome <- .imr_match_rows(dat$outcome, subjects$id, "outcome")
  groups <- lapply(object$model_bitstrings, function(g) which(subjects$subgroup == g))
  if (any(lengths(groups) < k)) .imr_abort("`k` must not exceed the sample size of any modelled availability subgroup.")
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(object$list_hyperpara[7L])
  # Generate every split and sampler seed before fitting, since imr() seeds R.
  partitions <- lapply(seq_len(rounds), function(r) {
    folds <- integer(nrow(subjects))
    for (idx in groups) {
      strata <- if (object$type_outcome == "continuous") rep(1, length(idx)) else {
        outcome[idx, if (object$type_outcome == "binary") 2L else 3L]
      }
      folds[idx] <- .imr_cv_folds(strata, k)
    }
    folds
  })
  seeds <- matrix(sample.int(.Machine$integer.max, rounds * k, replace = TRUE), rounds, k)
  labels <- c(object$model_bitstrings, "all")
  total <- subset <- matrix(NA_real_, rounds, length(labels), dimnames = list(NULL, labels))
  records <- vector("list", rounds)
  score <- function(idx, prediction) .imr_cv_accuracy(
    object$type_outcome, prediction[idx], outcome[idx, , drop = FALSE])
  for (r in seq_len(rounds)) {
    folds <- partitions[[r]]
    prediction <- rep(NA_real_, nrow(subjects))
    fold_scores <- matrix(NA_real_, k, length(labels))
    for (fold in seq_len(k)) {
      if (verbose) cat(sprintf("CV round %d/%d, fold %d/%d\n", r, rounds, fold, k))
      test <- which(folds == fold)
      train <- which(folds != fold)
      fitted <- .imr_cv_refit(object, subjects$id[train], seeds[r, fold], verbose)
      test_platforms <- lapply(dat$platforms, function(x) x[x$id %in% subjects$id[test], , drop = FALSE])
      present <- which(vapply(test_platforms, nrow, integer(1L)) > 0L)
      covariates <- if (!is.null(object$terms)) {
        object$formula_data[object$formula_data[[object$formula_id]] %in% subjects$id[test], , drop = FALSE]
      } else if (!is.null(dat$covariates)) {
        .imr_match_rows(dat$covariates, subjects$id[test], "covariates")
      } else NULL
      if (length(fitted$covariate_names) == 0L) covariates <- NULL
      predicted <- do.call(rbind, stats::predict(fitted, test_platforms[present],
        platform_names = as.character(present), covariates = covariates,
        max_models = max_models, verbose = verbose))
      prediction[test] <- predicted$predict[match(subjects$id[test], predicted$id)]
      if (any(!is.finite(prediction[test]))) {
        .imr_abort("A training-fold model did not return finite predictions for every held-out subject.")
      }
      for (g in seq_along(groups)) fold_scores[fold, g] <- score(intersect(test, groups[[g]]), prediction)
      fold_scores[fold, length(labels)] <- score(test, prediction)
    }
    for (g in seq_along(groups)) total[r, g] <- score(groups[[g]], prediction)
    total[r, length(labels)] <- score(seq_len(nrow(subjects)), prediction)
    # Undefined fold metrics remain NA; do not silently discard difficult folds.
    subset[r, ] <- colMeans(fold_scores)
    records[[r]] <- data.frame(round = r, fold = folds, subjects,
                               predict = prediction, row.names = NULL)
  }
  results <- list(total_cindex = total, subset_cindex = subset)
  attr(results, "metric") <- switch(object$type_outcome,
    right.censored = "C-index", binary = "AUC", continuous = "MSE")
  attr(results, "predictions") <- do.call(rbind, records)
  attr(results, "validation") <- "training-fold refits"
  results
}

# Keep folds nonempty and balanced even when separate strata contain < k rows.
.imr_cv_folds <- function(strata, k) {
  out <- integer(length(strata))
  offset <- 0L
  for (level in unique(strata)) {
    idx <- which(strata == level)
    idx <- idx[sample.int(length(idx))]
    out[idx] <- (offset + seq_along(idx) - 1L) %% k + 1L
    offset <- offset + length(idx)
  }
  out
}

.imr_cv_refit <- function(object, train_ids, seed, verbose = FALSE) {
  dat <- object$input_data
  # Keep rows outside the eligible cohort for otherwise unused platforms;
  # subgroup_data() intersects these with the training outcome IDs before any
  # normalization or sampling. No held-out outcome is passed to the refit.
  train_platforms <- lapply(dat$platforms, function(x) {
    keep <- x$id %in% train_ids | !x$id %in% dat$availability$id
    x[keep, , drop = FALSE]
  })
  # A platform occurring only in filtered-out subgroups still needs a valid
  # input frame, but none of its rows will enter a training subgroup.
  for (p in seq_along(train_platforms)) {
    if (nrow(train_platforms[[p]]) == 0L && length(object$platform_models[[p]]) == 0L)
      train_platforms[[p]] <- dat$platforms[[p]]
  }
  h <- object$list_hyperpara
  control <- list(type_outcome = object$type_outcome, method = object$method,
    ssize = 0L, nu = object$nu, h0 = h[1L], hh = h[2L],
    sig_alpha_psi = h[3:4], thet_alph_bet = h[5:6],
    sample_mcmc = unname(object$sample_mcmc[c("total", "burnin")]),
    survival_scale = if (object$type_outcome == "right.censored") object$response_scale else "identity",
    seed = seed, verbose = verbose)
  if (!is.null(object$terms)) {
    id <- object$formula_id
    train_platforms <- lapply(train_platforms, function(x) {names(x)[1L] <- id; x})
    return(do.call(imr, c(list(formula = object$formula,
      data = object$formula_data[object$formula_data[[id]] %in% train_ids, , drop = FALSE],
      platforms = train_platforms, id = id), control)))
  }
  do.call(imr, c(list(platform_data_list = train_platforms,
    outcome = .imr_match_rows(dat$outcome, train_ids, "outcome"),
    cov = if (!is.null(dat$covariates)) .imr_match_rows(dat$covariates, train_ids, "covariates") else NULL), control))
}

.imr_cv_accuracy <- function(type, prediction, outcome) {
  if (!length(prediction)) return(NA_real_)
  y <- outcome[[2L]]
  if (type == "continuous") return(mean((prediction - y)^2))
  if (type == "binary") {
    positive <- sum(y == 1)
    negative <- sum(y == 0)
    if (!positive || !negative) return(NA_real_)
    return((sum(rank(prediction, ties.method = "average")[y == 1]) -
              positive * (positive + 1) / 2) / (positive * negative))
  }
  .Call("imrConcordance", as.double(prediction), as.double(y),
        as.integer(outcome[[3L]]), PACKAGE = "IntegMultiReg")
}
