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
#'   Bayesian model averaging with `model_set = "ranked_unique"` and for refit
#'   prediction (default `100`). Must be positive. `model_set = "draws"` uses
#'   every retained draw, including repeated states, and ignores this limit.
#' @param cv_method Validation algorithm: `"legacy"` (default), `"refit"`,
#'   or `"importance"`. This is independent of the fitted IMR/BMS `method`.
#' @param verbose Logical; if `TRUE`, print fold progress and sampler
#'   diagnostics.  Defaults to `FALSE`.
#' @param workers Positive integer number of PSOCK worker processes (default
#'   `1L`, serial). Applies to all three CV methods. At most `k * rounds`
#'   processes are used; no automatic CPU detection is performed.
#' @param ridge Nonnegative post-fit diagonal penalty, including the intercept.
#'   `NULL` uses `0.001`; `0` requests the unpenalized estimate in the paper and
#'   original code. A singular solve stops with fold and subgroup context; no
#'   penalty or generalized inverse is silently substituted. Not used by refit.
#' @param model_set Post-fit states to average: `"ranked_unique"` (legacy
#'   default) ranks distinct states and keeps at most `max_models`; `"draws"`
#'   (importance default) retains empirical multiplicities and draw order.
#' @param df_method Predictive-density degrees of freedom: `"legacy_integer"`
#'   (legacy default) truncates `2 * shape + n_train` to an integer;
#'   `"fractional"` (importance default) keeps its numeric value. Post-fit only.
#' @param score_method Metric convention, independent of predictions:
#'   `"legacy"` (legacy default) retains historical AUC/concordance rules;
#'   `"standard"` (refit and importance default) uses current tie and comparable
#'   pair rules. MSE has the same definition in both. `NULL` selects the default.
#' @param folds Optional data frame with `id`, `round`, and `fold`. Every
#'   modelled subject must occur once per round; labels are consecutive integers
#'   starting at one. Every subgroup-fold needs training and test subjects.
#'   Omitted `k` and `rounds` are inferred; explicit values must match. Optional
#'   `row_order` is a permutation within each subgroup-round that preserves
#'   post-fit summation order. Reuse `result$control$folds` for exact replay.
#' @param fold_rng Post-fit fold generator starting point: `NULL` or `"reset"`
#'   initializes GSL from the fit seed (existing default); `"continue"` resumes
#'   the saved state at the end of fitting, as the original program did.
#'   Continuation requires a new fit saved on a compatible platform. Cannot
#'   be supplied with explicit `folds` or with refit CV.
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
#'   `control` records effective options, sampler/numerical conventions, version,
#'   seed and actual fold membership/order. It is additional metadata; the first
#'   five fields retain their existing meanings.
#'   Refit results also record `control$refit_seeds`, a rounds-by-fold matrix.
#'   Supplied folds use the same fitting seed plan as generated folds, so saved
#'   folds replay predictions under the same runtime and RNG kind.
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
#' paper. Its default 0.001 ridge stabilization differs from the unpenalized
#' coefficient estimate; set `ridge = 0` to select that estimate.
#' Binary and continuous outcomes extend the survival procedure.
#'
#' Mode names provide defaults. `ridge`, `model_set`, `df_method`, and
#' `score_method` can be chosen independently. Explicit post-fit-only arguments
#' are rejected in refit mode. The original released CV code combines `ridge = 0`,
#' `model_set = "ranked_unique"`, `df_method = "legacy_integer"`,
#' `score_method = "legacy"` and `max_models = 100`. Paper equations 6--7 use
#' `ridge = 0`, `model_set = "draws"` and `df_method = "fractional"`.
#' Neither combination alone reproduces a historical experiment: sampler,
#' numerical controls, data, initial states and random stream also matter.
#'
#' Both post-fit modes use the historical GSL fold generator, stratified by
#' event status for survival. Refit uses R's generator and additionally
#' stratifies binary outcomes by class. Identical seeds therefore do not give
#' identical partitions across refit and post-fit modes. Within each mode,
#' results are reproducible and the caller's R RNG state is preserved.
#' Hyperparameters are fixed; data-driven tuning needs an outer validation layer.
#'
#' Parallel execution preserves each method's partitions, seeds and output
#' ordering. Each process holds its own fit and training-fold workspace, with
#' single-threaded mathematical libraries; memory use grows with `workers`.
#' Post-fit matrix products may use up to 8 MiB of temporary training-column
#' storage per worker, in addition to the model-index cache and fitted data.
#' Process startup can make short jobs slower. Use only one parallel layer
#' when running multiple experiments. Worker failures stop the entire call;
#' no partial result or silent serial fallback is returned. Worker warnings
#' are relayed in task order after computation. Verbose worker output is not
#' guaranteed to appear in the calling console.
#'
#' For parallel formula refits, custom functions must be available in a
#' serializable formula environment (for example a local closure), or use a
#' package-qualified function name. The caller's global workspace is not
#' exported to workers. `workers = 1L` retains ordinary formula evaluation.
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
                   cv_method = c("legacy", "refit", "importance"),
                   workers = 1L, ridge = NULL, model_set = NULL,
                   df_method = NULL, score_method = NULL, folds = NULL,
                   fold_rng = NULL) {
  if (!inherits(object, "imr")) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  validate_imr(object)
  cv_method <- match.arg(cv_method)
  settings <- .imr_cv_settings(cv_method, ridge, model_set, df_method, score_method, fold_rng)
  if (!is.null(folds) && !is.null(fold_rng))
    .imr_abort("`fold_rng` cannot be supplied with explicit `folds`.")
  if (cv_method != "refit") .imr_cv_rng_state(object, settings$fold_rng)
  if (!is.null(folds)) {
    settings$fold_rng <- NULL
    supplied <- .imr_cv_validate_folds(folds, object,
      if (missing(k)) NULL else k, if (missing(rounds)) NULL else rounds)
    k <- supplied$k
    rounds <- supplied$rounds
    folds <- supplied$folds
  }
  .imr_check_flag(verbose, "verbose")
  k <- .imr_check_integer_scalar(k, "k", min = 2)
  rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  max_models <- .imr_check_integer_scalar(max_models, "max_models", min = 1)
  workers <- .imr_check_integer_scalar(workers, "workers", min = 1)
  if (as.double(k) * rounds > .Machine$integer.max)
    .imr_abort("The requested number of CV tasks exceeds the supported index limit.")
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
    return(.imr_cv_postfit(object, k, rounds, max_models, verbose, cv_method, workers,
                            settings, folds))
  }
  .imr_cv_refit_result(object, k, rounds, max_models, verbose, workers, settings, folds)
}

.imr_cv_refit_result <- function(object, k, rounds, max_models, verbose,
                                 workers = 1L, settings = .imr_cv_settings("refit"),
                                 supplied_folds = NULL) {
  control <- object$control
  model <- object$model
  preprocessing <- object$preprocessing
  dat <- preprocessing$input_data
  subjects <- dat$availability[
    dat$availability$subgroup %in% model$subgroup_names, c("id", "subgroup"), drop = FALSE]
  outcome <- .imr_match_rows(dat$outcome, subjects$id, "outcome")
  groups <- lapply(model$subgroup_names, function(g) which(subjects$subgroup == g))
  if (any(lengths(groups) < k)) .imr_abort("`k` must not exceed the sample size of any modelled availability subgroup.")
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(control$seed)
  plan <- .imr_cv_refit_plan(groups, outcome, control$outcome_type, k, rounds)
  partitions <- if (is.null(supplied_folds)) plan$partitions else lapply(seq_len(rounds), function(round) {
    rows <- supplied_folds[supplied_folds$round == round, , drop = FALSE]
    rows$fold[match(subjects$id, rows$id)]
  })
  seeds <- plan$seeds
  tasks <- unlist(lapply(seq_len(rounds), function(r) {
    lapply(seq_len(k), function(fold) {
      list(round = r, fold = fold, seed = seeds[r, fold],
           train_ids = subjects$id[partitions[[r]] != fold],
           test_ids = subjects$id[partitions[[r]] == fold])
    })
  }), recursive = FALSE)
  fold_predictions <- .imr_cv_map(tasks, .imr_cv_refit_task, workers,
                                 object = object, max_models = max_models,
                                 verbose = verbose)
  labels <- c(model$subgroup_names, "all")
  total <- subset <- matrix(NA_real_, rounds, length(labels), dimnames = list(NULL, labels))
  records <- vector("list", rounds)
  score <- function(idx, prediction) .imr_cv_accuracy(
    control$outcome_type, prediction[idx], outcome[idx, , drop = FALSE],
    settings$score_method)
  for (r in seq_len(rounds)) {
    folds <- partitions[[r]]
    prediction <- rep(NA_real_, nrow(subjects))
    fold_scores <- matrix(NA_real_, k, length(labels))
    for (fold in seq_len(k)) {
      test <- which(folds == fold)
      prediction[test] <- fold_predictions[[(r - 1L) * k + fold]]
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
    validation = "refit",
    control = c(.imr_cv_control(settings, object, k, rounds, max_models,
      if (is.null(supplied_folds)) do.call(rbind, lapply(seq_len(rounds), function(r) {
        data.frame(id = subjects$id, round = r, fold = partitions[[r]],
          row_order = as.integer(stats::ave(seq_len(nrow(subjects)), subjects$subgroup, FUN = seq_along)))
      })) else supplied_folds,
      if (is.null(supplied_folds)) "r" else "supplied"),
      list(refit_seeds = seeds))
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
    sampler_method = control$sampler_method %||% "legacy",
    standardize = control$standardize %||% TRUE, initial = control$initial,
    min_subgroup_size = 0L, nu = priors$nu,
    forced_prior_scale = priors$forced_scale,
    molecular_prior_scale = priors$molecular_scale,
    residual_prior = priors$residual,
    interaction_prior = priors$interaction,
    draws = control$mcmc$draws, burnin = control$mcmc$burnin,
    survival_scale = if (control$outcome_type == "right.censored") control$response_scale else "identity",
    seed = seed, verbose = verbose)
  refit_control <- c(refit_control, .imr_fit_numerical_control(control))
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

.imr_cv_accuracy <- function(type, prediction, outcome, score_method = "standard") {
  if (!length(prediction)) return(NA_real_)
  y <- outcome[[2L]]
  if (type == "continuous") return(mean((prediction - y)^2))
  if (score_method == "legacy") return(.imr_call_cv_legacy_score(type, prediction, outcome))
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
