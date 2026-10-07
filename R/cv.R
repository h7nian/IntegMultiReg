#' Cross-Validated Predictive Accuracy of an IMR Fit
#'
#' @description
#' `cv_imr()` evaluates predictive accuracy of a model fitted with [imr()] using
#' repeated \eqn{K}-fold splits. `cv_method = "refit"` (the default) fits each
#' training fold independently. `"reweight"` reuses a full-data fit and applies
#' inverse-density importance weights to its selection states.
#' Fit a separate object with `model_variant = "bms"` for a no-borrowing comparison.
#' The accuracy measure depends on the outcome type: the concordance index
#' (C-index) for right-censored outcomes, the area under the ROC curve (AUC) for
#' binary outcomes, and the mean squared error (MSE) for continuous outcomes.
#'
#' @section Refit and post-fit calculations:
#' For refit validation, a subject \eqn{i} in fold \eqn{f} receives
#' \deqn{\widehat y_i^{\mathrm{CV}}=
#'       \widehat g_{-f}(x_i),}{CV prediction for subject i in fold f = prediction from the model fitted without fold f.}
#' where \eqn{\widehat g_{-f}}{g_hat_minus_f} is fitted using only the training subjects.
#' Each fold reruns [imr()] and [predict.imr()] on training data, including
#' preprocessing, selection sampling and ranked-model prediction weights;
#' the point-prediction rule is described in [predict.imr()].
#'
#' Reweighting reuses full-fit selection states, transformed predictors and
#' augmented-response means. Within each state, its fold
#' coefficient estimate is
#' \deqn{\widehat b_m=(Z_{-f,m}^{T}Z_{-f,m}+\lambda I)^{-1}
#'                      Z_{-f,m}^{T}\bar y_{-f}^{*},}{beta_m = inverse(transpose(Z_train,m) * Z_train,m + ridge * I) * transpose(Z_train,m) * mean_working_y_train.}
#' where \eqn{\lambda}{lambda} is `ridge` and the diagonal penalty includes the
#' intercept. An unpenalized singular system is an error. If \eqn{a_m} is the
#' implemented inverse predictive-density log score for the held-out working
#' responses, the normalized fold weights are
#' \deqn{w_m=\frac{\exp(a_m-a_{\max})}
#'                  {\sum_h\exp(a_h-a_{\max})}.}{w_m = exp(a_m - a_max) / sum_h exp(a_h - a_max).}
#' Predictions average the fold model predictions with these weights, applying
#' the probit link before averaging for binary outcomes. `model_set = "all_draws"`
#' retains every state occurrence; `model_set = "top_unique"` uses at most
#' `max_models` ranked distinct states. The empirical `"all_draws"` calculation
#' is motivated by equations 6--7 of Chekouo et al. (2017). Both collections
#' reuse the full-data fit and therefore differ from training-fold refitting.
#'
#' @section Accuracy measures:
#' For a scored set of \eqn{n} subjects, continuous-outcome error is
#' \deqn{\mathrm{MSE}=\frac{1}{n}\sum_i(y_i-\widehat y_i)^2.}{MSE = sum_i (y_i - prediction_i)^2 / n.}
#' Binary AUC gives half credit to tied predictions:
#' \deqn{\mathrm{AUC}=\frac{1}{n_1n_0}
#'       \sum_{i:y_i=1}\sum_{j:y_j=0}
#'       \{I(\widehat p_i>\widehat p_j)+\frac{1}{2} I(\widehat p_i=\widehat p_j)\}.}{AUC = sum over positive-negative pairs of [I(p_positive > p_negative) + 0.5 * I(equal predictions)], divided by n_positive * n_negative.}
#' Here \eqn{n_1} and \eqn{n_0} are the class counts. For survival, let
#' \eqn{\mathcal C}{C} contain each comparable pair once, with the earlier
#' subject \eqn{i} first. For unequal times, the pair is comparable when
#' \eqn{t_i<t_j} and subject \eqn{i} has an observed event. For equal times,
#' put the event subject first and the censored subject second. Two tied
#' events are excluded.
#' The standard C-index is
#' \deqn{C=\frac{1}{|\mathcal C|}\sum_{(i,j)\in\mathcal C}
#'       \{I(\widehat y_i<\widehat y_j)+\frac{1}{2} I(\widehat y_i=\widehat y_j)\}.}{C-index = sum over ordered comparable pairs of [I(prediction_i < prediction_j) + 0.5 * I(equal predictions)], divided by the number of comparable pairs.}
#' Larger predicted working survival times mean longer survival. AUC is
#' undefined for a single class, and C is undefined without comparable pairs;
#' these return `NA`. The same scoring rules apply to both validation methods.
#'
#' @section Pooled scores and fold means:
#' Each output row is one validation round. If \eqn{A} is the chosen scoring
#' function and \eqn{D_f} contains fold \eqn{f}'s outcomes and predictions,
#' \deqn{A_{\mathrm{pooled}}=A\Bigl(\bigcup_{f=1}^{K}D_f\Bigr),}{Pooled score = score of all folds combined,}
#' \deqn{A_{\mathrm{fold\ mean}}=\frac{1}{K}\sum_{f=1}^{K}A(D_f).}{Fold mean = sum of the K individual fold scores / K.}
#' The pooled score uses all subject predictions together. The fold mean gives
#' equal weight to folds; an undefined fold makes the fold mean `NA`. These
#' quantities can differ, especially for AUC and concordance. Both are reported
#' per availability subgroup and overall; rounds remain separate rows.
#'
#' @param object A fitted object of class `"imr"` returned by [imr()].
#' @param k Integer number of cross-validation folds per round (default `5`).
#'   Must be at least `2`, and each availability subgroup must contain at least
#'   `k` subjects.
#' @param rounds Integer number of independent cross-validation rounds to
#'   average over (default `2`).  Must be positive.
#' @param max_models Integer maximum number of selection models used for
#'   Bayesian model averaging with `model_set = "top_unique"` and for refit
#'   prediction (default `100`). Must be positive. `model_set = "all_draws"` uses
#'   every retained draw, including repeated states, and ignores this limit.
#' @param cv_method `"refit"` (default) rebuilds preprocessing and MCMC inside
#'   every training fold. `"reweight"` reuses the fitted selection draws with
#'   importance weights. To compare IMR and BMS, first fit separate
#'   `model_variant = "imr"` and `"bms"` objects.
#' @param verbose Logical; if `TRUE`, print fold progress and sampler
#'   diagnostics.  Defaults to `FALSE`.
#' @param workers Positive integer number of PSOCK worker processes (default
#'   `1L`, serial). Applies to both CV methods. At most `k * rounds`
#'   processes are used; no automatic CPU detection is performed.
#' @param ridge Nonnegative post-fit diagonal penalty, including the intercept.
#'   `NULL` uses `0.001`; `0` requests the unpenalized estimate in the paper and
#'   original code. A singular solve stops with fold and subgroup context; no
#'   penalty or generalized inverse is silently substituted. Not used by refit.
#' @param model_set For reweighting, `"all_draws"` (default) retains every
#'   selection draw, including repeated states. `"top_unique"` ranks distinct
#'   selection models by their observed frequency and retains at most
#'   `max_models`. These are state collections, not different CV algorithms.
#' @param df_method For reweighting, `"fractional"` (default) keeps the numeric
#'   predictive degrees of freedom `2 * shape + n_train`; `"integer"` truncates
#'   them to an integer. Not used by refit.
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
#' @return An `imr_cv` object: a named list with a compact [print.imr_cv()]
#'   method. `pooled` and `fold_mean` are numeric matrices of dimension
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
#' Reweighting with `model_set = "all_draws"` follows the empirical importance
#' average in equations 6--7 of the paper. It conditions on full-fit predictor
#' transformations and posterior-mean augmented responses. The response plug-in
#' is part of Section 4.1; the default ridge penalty of 0.001 differs from the
#' paper's unpenalized coefficient estimate. Set `ridge = 0` explicitly for that
#' estimate, provided every training design has full column rank.
#'
#' Held-out outcomes enter inverse-density weights as an importance correction;
#' their use alone does not identify an error. This approximation has a different
#' scope from training-fold refitting. Choosing `model_set = "top_unique"`
#' replaces the empirical draw collection with selected distinct states.
#' Historical bundles of settings are documented with their archived sources;
#' a mode label alone never establishes historical numerical reproduction.
#'
#' Supplied folds replay actual partitions. Refit uses R's generator and
#' stratifies binary outcomes by class. Generated reweighting folds use GSL and
#' stratify survival outcomes by event status. Matching seed numbers across the
#' two generators does not imply matching partitions. Within each method the
#' caller's R RNG state is preserved. Hyperparameters are fixed; tuning needs an
#' outer validation layer.
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
#' @references
#' Chekouo et al. (2017), Section 4.1, equations 6--7. \doi{10.1111/biom.12587}.
#' [Read paper](https://academic.oup.com/biometrics/article/73/2/615/7537638) |
#' [Publisher PDF](https://academic.oup.com/biometrics/article-pdf/73/2/615/55973435/biometrics_73_2_615.pdf).
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
#' cv <- cv_imr(fit, k = 5, rounds = 2, cv_method = "reweight")
#' cv$pooled
#' }
#' @export
cv_imr <- function(object, k = 5, rounds = 2,
                   max_models = 100, verbose = FALSE,
                   cv_method = c("refit", "reweight"),
                   workers = 1L, ridge = NULL, model_set = NULL,
                   df_method = NULL, folds = NULL,
                   fold_rng = NULL) {
  validate_imr_object(object)
  .imr_require_current_updates(object)
  cv_method <- match.arg(cv_method)
  settings <- .imr_cv_settings(cv_method, ridge, model_set, df_method, fold_rng)
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

#' Print Method for IMR Cross-Validation Results
#'
#' @section Interpreting the scores:
#' The display reads the stored `pooled` and `fold_mean` matrices and reports
#' the validation algorithm, fold design and scoring rule. For reweighting it
#' also shows the resolved penalty, state collection and degrees of freedom.
#' The model cap is shown when refitting or using top-ranked distinct states.
#' Consult `x$control` for all recorded settings. [cv_imr()] defines the
#' scoring formulas and explains why a pooled score can differ from the mean
#' fold score. `digits` changes significant digits in the display only;
#' stored values and subject-level predictions retain full precision.
#'
#' @param x An `imr_cv` object returned by [cv_imr()].
#' @param digits Number of significant digits for the reported scores.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the printed summary.
#' @seealso [cv_imr()]
#' @export
#' @examples
#' \donttest{
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            min_subgroup_size = 5, forced_prior_scale = 1,
#'            draws = 500, burnin = 250, seed = 1)
#' cv_imr(fit, k = 2, rounds = 1)
#' }
print.imr_cv <- function(x, digits = 3, ...) {
  control <- x$control
  label <- if (is.null(control$model_variant)) "IntegMultiReg" else toupper(control$model_variant)
  cat(label, "cross-validation:", x$metric, "by availability subgroup\n")
  if (identical(x$validation, "refit")) {
    cat("  Algorithm: refit preprocessing and MCMC within each training fold.\n")
  } else if (identical(x$validation, "reweight")) {
    cat("  Algorithm: reuse the full-data fit and reweight selection states.\n")
  } else {
    cat("  Archived validation configuration:", x$validation, "\n")
  }
  cat(sprintf("  %d %s of %d-fold validation; %s folds.\n",
              control$rounds, if (control$rounds == 1L) "round" else "rounds", control$k,
              switch(control$fold_source, supplied = "supplied", gsl = "GSL-generated",
                     r = "R-generated", control$fold_source)))
  if (identical(x$validation, "reweight")) {
    if (identical(control$model_set, "all_draws")) {
      cat("  States: all retained draws, preserving repeated-state counts.\n")
    } else {
      cat(sprintf("  States: at most %d top-ranked distinct selection models.\n", control$max_models))
    }
    cat(sprintf("  Ridge penalty: %s; predictive df: %s.\n", format(control$ridge),
                if (identical(control$df_method, "fractional")) "unrounded" else "integer-truncated"))
  } else if (identical(x$validation, "refit")) {
    cat(sprintf("  Prediction cap per fit: %d distinct selection models.\n", control$max_models))
  }
  if (!is.null(control$score_rule)) {
    description <- switch(control$score_rule,
      auc_half_ties = "pairwise AUC, with half credit for tied predictions",
      concordance_comparable_pairs = "C-index over comparable pairs, with half credit for tied predictions",
      mean_squared_error = "mean squared prediction error", control$score_rule)
    cat("  Scoring:", description, "\n")
  } else if (identical(control$score_method, "standard")) {
    cat("  Scoring: standard pairwise AUC/C-index or mean squared error.\n")
  } else {
    cat("  Scoring: archived rules; consult the recorded source and controls.\n")
  }
  cat("\nPooled out-of-fold score:\n")
  print(signif(x$pooled, digits))
  cat("\nMean fold score:\n")
  print(signif(x$fold_mean, digits))
  cat(sprintf("\n%d subject predictions in `predictions`;", nrow(x$predictions)))
  cat(" effective settings and a reusable fold table in `control`.\n")
  invisible(x)
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
  structure(class = "imr_cv", list(
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
  ))
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
  # .imr_subgroup_data() intersects these with the training outcome IDs before any
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
  refit_control <- list(outcome_type = control$outcome_type, model_variant = control$model_variant,
    standardize = control$standardize %||% TRUE, initial = control$initial,
    min_subgroup_size = 0L, nu = priors$nu,
    forced_prior_scale = priors$forced_scale,
    molecular_prior_scale = priors$molecular_scale,
    residual_prior = priors$residual,
    interaction_prior = priors$interaction,
    draws = control$mcmc$draws, burnin = control$mcmc$burnin,
    survival_scale = if (control$outcome_type == "right.censored") control$response_scale else "identity",
    seed = seed, verbose = verbose)
  refit_control <- c(refit_control, .imr_fit_numerical_control(control)[c("laplace_max_iter", "laplace_tolerance")])
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
  if (score_method == "original") return(.imr_call_cv_legacy_score(type, prediction, outcome))
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
