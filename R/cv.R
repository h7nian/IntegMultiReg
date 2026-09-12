#' Cross-Validated Predictive Accuracy of an IMR Fit
#'
#' @description
#' `cv_imr()` evaluates predictive accuracy of a model fitted with [imr()] using
#' repeated \eqn{K}-fold splits. `cv_method` selects historical post-fit
#' validation (`"legacy"`, the default), training-fold refitting (`"refit"`),
#' or a paper-derived importance average over retained draws (`"importance"`).
#' Fit a separate object with `method = "bms"` for a no-borrowing comparison.
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
#'   Bayesian model averaging in legacy and refit modes (default `100`). Must
#'   be positive. Importance mode uses every retained draw, including repeated
#'   states, and does not truncate using this argument.
#' @param cv_method Validation algorithm: `"legacy"` (default), `"refit"`,
#'   or `"importance"`. This is independent of the fitted IMR/BMS `method`.
#' @param verbose Logical; if `TRUE`, print fold progress and sampler
#'   diagnostics.  Defaults to `FALSE`.
#'
#' @return A named list. `pooled` and `fold_mean` are numeric matrices of dimension
#'   `rounds` x `(n_subgroups + 1)`, whose last column corresponds to all
#'   subjects pooled and whose remaining columns are named by the availability
#'   subgroup bitstrings:
#'   \item{pooled}{Accuracy computed on the pooled cross-validated
#'     predictions within each availability subgroup (and overall).}
#'   \item{fold_mean}{Accuracy computed fold-by-fold within each
#'     availability subgroup.}
#'   `predictions` contains subject IDs, subgroups, rounds, folds and out-of-fold
#'   predictions. `metric` names the accuracy measure and `validation` records
#'   the selected `cv_method` as a character string.
#'   Undefined AUCs (single class) and C-indices (no comparable pairs) are `NA`.
#'
#' @details
#' Only subgroups retained in `object` are evaluated. Their training subsets
#' remain included even when they are smaller than the original
#' `min_subgroup_size`.
#' In refit mode, preprocessing, formula transformations, selection MCMC and
#' prediction weights are recomputed on training rows. Runtime is roughly
#' `k * rounds` full fits. This is the algorithm used by version 0.1.4.
#'
#' Legacy mode restores the post-fit CV algorithm from version 0.1.0: ranked
#' distinct full-fit selection models, full-fit standardization and augmented
#' response means, fold-specific ridge coefficient estimates (penalty 0.001),
#' inverse predictive-density weights and historical AUC/concordance scoring.
#' It uses the supplied fit; reproducing a 0.1.0 run also requires its original
#' fit and preprocessing. Undefined metrics return NA instead of NaN.
#'
#' Importance mode follows the empirical importance average in equations 6--7
#' of Chekouo et al. (2017), retaining MCMC state multiplicities rather than
#' selecting equally weighted distinct models. It uses fractional predictive
#' degrees of freedom and the current scoring definitions. The full-fit
#' augmented response mean plug-in is explicitly part of Section 4.1 of the
#' paper. The 0.001 ridge stabilization differs from its unpenalized coefficient
#' estimate; this mode is not an exact reproduction of the original application.
#' Binary and continuous outcomes extend the survival procedure.
#'
#' Both post-fit modes use the historical GSL fold generator, stratified by
#' event status for survival. Refit uses R's generator and additionally
#' stratifies binary outcomes by class. Identical seeds therefore do not give
#' identical partitions across refit and post-fit modes. Within each mode,
#' results are reproducible and the caller's R RNG state is preserved.
#' Hyperparameters are fixed; data-driven tuning needs an outer validation layer.
#'
#' @seealso [imr()], [predict.imr()]
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
#' cv <- cv_imr(fit, k = 5, rounds = 2, cv_method = "legacy")
#' cv$pooled
#' }
#' @export
cv_imr <- function(object, k = 5, rounds = 2,
                   method = NULL,
                   max_models = 100, verbose = FALSE,
                   cv_method = c("legacy", "refit", "importance")) {
  if (!inherits(object, "imr")) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  validate_imr(object)
  cv_method <- match.arg(cv_method)
  .imr_check_flag(verbose, "verbose")
  k <- .imr_check_integer_scalar(k, "k", min = 2)
  rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  max_models <- .imr_check_integer_scalar(max_models, "max_models", min = 1)
  control <- object$control
  model <- object$model
  preprocessing <- object$preprocessing
  object_method <- control$method
  if (is.null(method)) {
    method <- object_method
  } else {
    method <- match.arg(method, c("imr", "bms"))
    if (!identical(method, object_method)) {
      .imr_abort(
        "`method` must match the fitted object; fit a separate `imr(..., method = \"bms\")` object for BMS validation."
      )
    }
  }
  if (any(model$sample_sizes < k)) {
    .imr_abort(
      "`k` must not exceed the sample size of any modelled availability subgroup."
    )
  }
  if (is.null(preprocessing$input_data)) {
    .imr_abort("Refit this model to retain the raw inputs required for cross-validation.")
  }
  validate_imr_data(preprocessing$input_data)
  if (control$outcome_type == "right.censored" && is.null(control$response_scale)) {
    .imr_abort("Refit this survival model with an explicit `survival_scale`.")
  }
  if (!is.null(preprocessing$terms) && is.null(preprocessing$formula_data)) {
    .imr_abort("Refit this formula model to retain the raw formula data required for cross-validation.")
  }
  if (cv_method != "refit") {
    return(.imr_cv_postfit(object, k, rounds, max_models, verbose, cv_method))
  }
  dat <- preprocessing$input_data
  subjects <- dat$availability[
    dat$availability$subgroup %in% model$subgroup_names, c("id", "subgroup"), drop = FALSE]
  outcome <- .imr_match_rows(dat$outcome, subjects$id, "outcome")
  groups <- lapply(model$subgroup_names, function(g) which(subjects$subgroup == g))
  if (any(lengths(groups) < k)) .imr_abort("`k` must not exceed the sample size of any modelled availability subgroup.")
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(control$seed)
  # Generate every split and sampler seed before fitting, since imr() seeds R.
  partitions <- lapply(seq_len(rounds), function(r) {
    folds <- integer(nrow(subjects))
    for (idx in groups) {
      strata <- if (control$outcome_type == "continuous") rep(1, length(idx)) else {
        outcome[idx, if (control$outcome_type == "binary") 2L else 3L]
      }
      folds[idx] <- .imr_cv_folds(strata, k)
    }
    folds
  })
  seeds <- matrix(sample.int(.Machine$integer.max, rounds * k, replace = TRUE), rounds, k)
  labels <- c(model$subgroup_names, "all")
  total <- subset <- matrix(NA_real_, rounds, length(labels), dimnames = list(NULL, labels))
  records <- vector("list", rounds)
  score <- function(idx, prediction) .imr_cv_accuracy(
    control$outcome_type, prediction[idx], outcome[idx, , drop = FALSE])
  for (r in seq_len(rounds)) {
    folds <- partitions[[r]]
    prediction <- rep(NA_real_, nrow(subjects))
    fold_scores <- matrix(NA_real_, k, length(labels))
    for (fold in seq_len(k)) {
      if (verbose) cat(sprintf("CV round %d/%d, fold %d/%d\n", r, rounds, fold, k))
      test <- which(folds == fold)
      train <- which(folds != fold)
      prediction[test] <- .imr_cv_refit_fold(
        object, subjects$id[train], subjects$id[test], seeds[r, fold],
        max_models, verbose)
      for (g in seq_along(groups)) fold_scores[fold, g] <- score(intersect(test, groups[[g]]), prediction)
      fold_scores[fold, length(labels)] <- score(test, prediction)
    }
    for (g in seq_along(groups)) total[r, g] <- score(groups[[g]], prediction)
    total[r, length(labels)] <- score(seq_len(nrow(subjects)), prediction)
    # Undefined fold metrics remain NA; do not silently discard difficult folds.
    subset[r, ] <- colMeans(fold_scores)
    records[[r]] <- data.frame(round = r, fold = folds, subjects,
                               prediction = prediction, row.names = NULL)
  }
  list(
    pooled = total,
    fold_mean = subset,
    predictions = do.call(rbind, records),
    metric = switch(control$outcome_type,
      right.censored = "C-index", binary = "AUC", continuous = "MSE"),
    validation = "refit"
  )
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
  control <- object$control
  model <- object$model
  preprocessing <- object$preprocessing
  priors <- control$priors
  dat <- preprocessing$input_data
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
    if (nrow(train_platforms[[p]]) == 0L && length(model$platform_subgroups[[p]]) == 0L)
      train_platforms[[p]] <- dat$platforms[[p]]
  }
  refit_control <- list(outcome_type = control$outcome_type, method = control$method,
    min_subgroup_size = 0L, nu = priors$nu,
    forced_prior_scale = priors$forced_scale,
    molecular_prior_scale = priors$molecular_scale,
    residual_prior = priors$residual,
    interaction_prior = priors$interaction,
    draws = control$mcmc$draws, burnin = control$mcmc$burnin,
    survival_scale = if (control$outcome_type == "right.censored") control$response_scale else "identity",
    seed = seed, verbose = verbose)
  if (!is.null(preprocessing$terms)) {
    id <- preprocessing$id
    train_platforms <- lapply(train_platforms, function(x) {names(x)[1L] <- id; x})
    return(do.call(imr, c(list(x = preprocessing$formula,
      data = preprocessing$formula_data[preprocessing$formula_data[[id]] %in% train_ids, , drop = FALSE],
      platforms = train_platforms, id = id), refit_control)))
  }
  do.call(imr, c(list(x = train_platforms,
    outcome = .imr_match_rows(dat$outcome, train_ids, "outcome"),
    covariates = if (!is.null(dat$covariates)) .imr_match_rows(dat$covariates, train_ids, "covariates") else NULL), refit_control))
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
  .Call("imr_concordance", as.double(prediction), as.double(y),
        as.integer(outcome[[3L]]), PACKAGE = "IntegMultiReg")
}
