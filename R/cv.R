#' Cross-Validate a Joint IMR Model
#'
#' The default refits the complete model within each training fold, including
#' preprocessing and formula transformations. Reweighting uses Pareto-smoothed
#' importance ratios from the stored joint posterior instead of refitting.
#'
#' @section Algorithms:
#' For refitting, every held-out prediction uses a model fitted without that
#' fold. Hyperparameters stay fixed; selecting them from the data requires an
#' outer validation layer. Matching actual folds pairs candidate comparisons.
#'
#' Reweighting uses the joint observed-response likelihood of the held-out fold:
#' \deqn{\log r_{bc,f}=-\sum_{i\in f}\log p(y_i\mid\beta^{(b,c)},(\sigma^2)^{(b,c)}).}{Log importance ratio for a fold = negative sum of its observed-data log likelihood over each posterior draw.}
#' Binary likelihoods integrate out latent utilities; censored survival uses
#' survival probabilities. The loo package smooths these ratios by PSIS and
#' estimates their Pareto k and effective sample size. All retained joint draws
#' contribute; no model ranking, ridge solve or latent-response plug-in is used.
#' Reweighting conditions on full-fit preprocessing. It estimates a different
#' procedure from training-fold refitting, and large held-out folds often have
#' unstable importance weights. Inspect diagnostics; use refit when unreliable.
#'
#' @section Scores:
#' Continuous outcomes use mean squared error. Binary AUC gives half credit to
#' tied predictions. Survival concordance counts comparable pairs, with the
#' earlier event first; tied events are excluded and an event tied with censoring
#' is comparable. An undefined fold makes its round's mean fold score `NA`.
#' Pooled scores instead combine all subject predictions within a round.
#'
#' @param object An `imr` fit. Refit CV retains its marginalization choice.
#'   Reweighting requires regression draws and is unavailable for
#'   `marginalize = "coefficients_and_variance"` fits.
#' @param k Folds per round, at least two. Inferred from supplied folds if omitted.
#' @param rounds Repetitions, inferred from supplied folds if omitted.
#' @param cv_method `"refit"` (default) or posterior `"reweight"`.
#' @param folds Optional data frame with `id`, `round` and `fold`; each modelled
#'   subject must appear exactly once per round.
#' @param seed Seed for folds and refit seed plans; defaults to the fit's seed.
#' @param workers Processes across folds. Each training fit runs chains serially.
#' @param verbose Print fold and sampling progress.
#' @return An `imr_cv` object with pooled/mean-fold score matrices, subject-level
#'   predictions, actual folds and seeds, and per-fold MCMC or PSIS diagnostics.
#'   Refit diagnostics give the maximum defined R-hat and minimum defined
#'   bulk/tail ESS, with counts of unassessed constant parameters. Reweighting
#'   records warnings from relative-efficiency estimation and PSIS by fold.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' cv_imr(fit, k = 2, rounds = 1)
#' @export
cv_imr <- function(object, k = 5L, rounds = 2L,
                   cv_method = c("refit", "reweight"), folds = NULL,
                   seed = NULL, workers = 1L, verbose = FALSE) {
  .imr_check_fit(object)
  cv_method <- match.arg(cv_method)
  if (cv_method == "reweight" && is.null(object$posterior$coefficients)) {
    .imr_abort("PSIS reweighting requires stored regression parameter draws; use cv_method = 'refit' for a Laplace selection fit.")
  }
  .imr_check_flag(verbose, "verbose")
  workers <- .imr_check_integer_scalar(workers, "workers", min = 1)
  seed <- .imr_check_integer_scalar(seed %||% object$control$seed, "seed", min = 0)
  if (!is.null(folds)) {
    checked <- .imr_cv_validate_folds(
      folds, object, if (missing(k)) NULL else k,
      if (missing(rounds)) NULL else rounds
    )
    k <- checked$k
    rounds <- checked$rounds
    folds <- checked$folds
  } else {
    k <- .imr_check_integer_scalar(k, "k", min = 2)
    rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  }
  if (as.double(k) * rounds > .Machine$integer.max) .imr_abort("Too many CV tasks.")
  data <- object$preprocessing$input_data
  validate_imr_data(data)
  subjects <- data$availability[data$availability$subgroup %in% object$model$subgroup_names,
    c("id", "subgroup"),
    drop = FALSE
  ]
  outcome <- .imr_match_rows(data$outcome, subjects$id, "outcome")
  groups <- lapply(object$model$subgroup_names, function(g) which(subjects$subgroup == g))
  if (any(lengths(groups) < k)) .imr_abort("Every modelled subgroup needs at least k subjects.")
  saved <- .imr_save_rng()
  on.exit(.imr_restore_rng(saved), add = TRUE)
  set.seed(seed)
  seeds <- matrix(sample.int(.Machine$integer.max, k * rounds, replace = TRUE), rounds, k)
  partitions <- if (is.null(folds)) {
    lapply(seq_len(rounds), function(r) {
      result <- integer(nrow(subjects))
      for (g in groups) {
        strata <- if (object$control$outcome_type == "continuous") {
          rep(1, length(g))
        } else {
          outcome[g, if (object$control$outcome_type == "binary") 2L else 3L]
        }
        result[g] <- .imr_cv_folds(strata, k)
      }
      result
    })
  } else {
    lapply(seq_len(rounds), function(r) {
      rows <- folds[folds$round == r, ]
      rows$fold[match(subjects$id, rows$id)]
    })
  }
  tasks <- unlist(lapply(seq_len(rounds), function(r) {
    lapply(seq_len(k), function(f) {
      list(
        round = r, fold = f, seed = seeds[r, f],
        train_ids = subjects$id[partitions[[r]] != f], test_ids = subjects$id[partitions[[r]] == f]
      )
    })
  }), recursive = FALSE)
  # Training-fold fits need raw inputs and settings, not full-data posterior
  # arrays or transformed design copies. Avoid serializing those to every worker.
  task_object <- if (cv_method == "refit") {
    list(
      control = object$control,
      model = object$model,
      preprocessing = list(
        input_data = object$preprocessing$input_data,
        terms = object$preprocessing$terms,
        formula = object$preprocessing$formula,
        formula_data = object$preprocessing$formula_data,
        id = object$preprocessing$id
      )
    )
  } else {
    object
  }
  answers <- .imr_map_tasks(tasks, .imr_cv_task, workers, object = task_object, method = cv_method, verbose = verbose)
  labels <- c(object$model$subgroup_names, "all")
  pooled <- fold_mean <- matrix(NA_real_, rounds, length(labels), dimnames = list(NULL, labels))
  records <- vector("list", rounds)
  all_groups <- c(groups, list(seq_len(nrow(subjects))))
  for (r in seq_len(rounds)) {
    prediction <- numeric(nrow(subjects))
    per_fold <- matrix(NA_real_, k, length(labels))
    for (f in seq_len(k)) {
      idx <- which(partitions[[r]] == f)
      prediction[idx] <- answers[[(r - 1L) * k + f]]$prediction
      for (g in seq_along(all_groups)) {
        rows <- intersect(idx, all_groups[[g]])
        per_fold[f, g] <- .imr_cv_accuracy(object$control$outcome_type, prediction[rows], outcome[rows, , drop = FALSE])
      }
    }
    for (g in seq_along(all_groups)) {
      pooled[r, g] <- .imr_cv_accuracy(
        object$control$outcome_type,
        prediction[all_groups[[g]]], outcome[all_groups[[g]], , drop = FALSE]
      )
    }
    fold_mean[r, ] <- colMeans(per_fold)
    records[[r]] <- data.frame(round = r, fold = partitions[[r]], subjects, prediction = prediction, row.names = NULL)
  }
  diagnostics <- do.call(rbind, lapply(seq_along(tasks), function(i) {
    cbind(round = tasks[[i]]$round, fold = tasks[[i]]$fold, answers[[i]]$diagnostics)
  }))
  if (cv_method == "refit" && any(diagnostics$rhat_failed > 0 | diagnostics$constant_chain_parameters > 0, na.rm = TRUE)) {
    .imr_warn("Some training-fold chains have poor MCMC diagnostics; inspect result$diagnostics before interpreting CV.")
  }
  if (cv_method == "reweight" && any(!diagnostics$reliable)) {
    .imr_warn("Some PSIS fold estimates are unreliable; inspect result$diagnostics and use cv_method = 'refit'.")
  }
  effective_mcmc <- object$control$mcmc
  if (cv_method == "refit") {
    effective_mcmc$workers <- 1L
    effective_mcmc$keep_latent <- FALSE
    effective_mcmc$seed <- NULL
  }
  structure(list(
    pooled = pooled, fold_mean = fold_mean, predictions = do.call(rbind, records),
    metric = switch(object$control$outcome_type,
      continuous = "MSE",
      binary = "AUC",
      right.censored = "C-index"
    ),
    validation = cv_method, diagnostics = diagnostics,
    control = list(
      k = k, rounds = rounds, seed = seed, model_variant = object$control$model_variant,
      folds = do.call(rbind, lapply(seq_len(rounds), function(r) {
        data.frame(id = subjects$id, round = r, fold = partitions[[r]], row.names = NULL)
      })),
      refit_seeds = if (cv_method == "refit") seeds else NULL,
      refit_chain_seeds = if (cv_method == "refit") {
        array(unlist(lapply(answers, `[[`, "chain_seeds")),
          dim = c(object$control$mcmc$chains, k, rounds), dimnames = list(chain = NULL, fold = NULL, round = NULL)
        )
      } else {
        NULL
      },
      refit_initial_seeds = if (cv_method == "refit") {
        array(unlist(lapply(answers, `[[`, "initial_seeds")),
          dim = c(object$control$mcmc$chains, k, rounds)
        )
      } else {
        NULL
      },
      rng_kind = RNGkind(), mcmc = effective_mcmc, priors = object$control$priors,
      marginalize = .imr_marginalize(object$control), numerical = object$control$numerical %||% imr_control(),
      preprocessing = if (cv_method == "refit") "training_fold" else "full_fit",
      package_version = "0.3.0"
    )
  ), class = "imr_cv")
}

.imr_cv_test_inputs <- function(object, test_ids) {
  prep <- object$preprocessing
  data <- prep$input_data
  platforms <- lapply(data$platforms, function(x) x[x$id %in% test_ids, , drop = FALSE])
  present <- which(vapply(platforms, nrow, 1L) > 0L)
  covariates <- if (!is.null(prep$terms)) prep$formula_data[prep$formula_data[[prep$id]] %in% test_ids, , drop = FALSE] else if (!is.null(data$covariates)) .imr_match_rows(data$covariates, test_ids, "covariates") else NULL
  list(platforms = platforms[present], platform_names = as.character(present), covariates = covariates)
}

.imr_cv_task <- function(task, object, method, verbose) {
  tryCatch(
    {
      if (verbose) cat(sprintf("CV round %d, fold %d\n", task$round, task$fold))
      data <- .imr_cv_test_inputs(object, task$test_ids)
      if (method == "refit") {
        fitted <- withCallingHandlers(.imr_cv_refit(object, task$train_ids, task$seed, verbose),
          imr_mcmc_warning = function(w) invokeRestart("muffleWarning")
        )
        result <- stats::predict(fitted, data$platforms,
          platform_names = data$platform_names,
          covariates = data$covariates
        )
        d <- fitted$diagnostics
        diagnostics <- data.frame(
          rhat_max = if (is.null(d) || !any(!is.na(d$rhat))) NA_real_ else max(d$rhat, na.rm = TRUE),
          rhat_failed = if (is.null(d)) NA_integer_ else sum(d$rhat > 1.01, na.rm = TRUE),
          ess_bulk_min = if (is.null(d) || !any(is.finite(d$ess_bulk))) NA_real_ else min(d$ess_bulk, na.rm = TRUE),
          ess_tail_min = if (is.null(d) || !any(is.finite(d$ess_tail))) NA_real_ else min(d$ess_tail, na.rm = TRUE),
          constant_chain_parameters = if (is.null(d)) NA_integer_ else sum(d$status == "constant_in_chain"),
          constant_parameters = if (is.null(d)) NA_integer_ else sum(d$status == "constant")
        )
      } else {
        log_likelihood <- .imr_holdout_log_likelihood(object, task$test_ids)
        scaled_likelihood <- exp(log_likelihood - max(log_likelihood))
        chain_id <- rep(seq_len(object$control$mcmc$chains), each = object$control$mcmc$draws)
        messages <- character()
        record_warning <- function(w) {
          messages <<- c(messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
        relative <- withCallingHandlers(
          loo::relative_eff(scaled_likelihood, chain_id = chain_id, cores = 1),
          warning = record_warning
        )
        if (!is.finite(relative) || relative <= 0) .imr_abort("Cannot estimate relative efficiency for importance weights.")
        psis <- withCallingHandlers(loo::psis(-log_likelihood, r_eff = relative, cores = 1),
          warning = record_warning
        )
        weights <- as.vector(stats::weights(psis, normalize = TRUE, log = FALSE))
        inputs <- .imr_prediction_inputs(object, data$platforms, data$platform_names, data$covariates)
        result <- .imr_predict_joint(object, inputs, weights = weights)
        k <- as.numeric(loo::pareto_k_values(psis))
        ess <- as.numeric(loo::psis_n_eff_values(psis))
        threshold <- min(.7, 1 - 1 / log10(length(weights)))
        diagnostics <- data.frame(
          pareto_k = k, pareto_k_threshold = threshold,
          importance_ess = ess, relative_efficiency = relative,
          reliable = is.finite(k) && k < threshold && !length(messages),
          warning = paste(unique(messages), collapse = "; ")
        )
      }
      prediction <- do.call(rbind, result)
      prediction <- prediction$prediction[match(task$test_ids, prediction$id)]
      if (any(!is.finite(prediction))) .imr_abort("Missing or non-finite held-out predictions.")
      list(
        prediction = prediction, diagnostics = diagnostics,
        chain_seeds = if (method == "refit") fitted$control$chain_seeds else NULL,
        initial_seeds = if (method == "refit") fitted$control$initial_seeds else NULL
      )
    },
    error = function(e) .imr_abort(sprintf("CV round %d, fold %d failed: %s", task$round, task$fold, conditionMessage(e)))
  )
}

.imr_holdout_log_likelihood <- function(object, ids) {
  samples <- as.double(object$control$mcmc$draws) * object$control$mcmc$chains
  result <- numeric(samples)
  outcome <- object$control$outcome_type
  for (g in seq_along(object$model$subgroup_names)) {
    rows <- which(object$preprocessing$subject_ids[[g]] %in% ids)
    if (!length(rows)) next
    design <- .imr_joint_design(object$model, object$preprocessing, g)
    beta <- .imr_flat_draws(object$posterior$coefficients[[g]])
    sd <- sqrt(as.vector(object$posterior$variance[, , g]))
    response <- object$preprocessing$response[[g]]
    # One observation at a time bounds temporary memory independently of fold size.
    for (i in rows) {
      mean <- drop(beta %*% design[i, ])
      value <- if (outcome == "binary") {
        direction <- if (response[i, 1L] == 1) 1 else -1
        stats::pnorm(direction * mean / sd, log.p = TRUE)
      } else if (outcome == "right.censored" && response[i, 2L] == 0) {
        stats::pnorm((response[i, 1L] - mean) / sd, lower.tail = FALSE, log.p = TRUE)
      } else {
        density <- stats::dnorm(response[i, 1L], mean, sd, log = TRUE)
        if (outcome == "right.censored" && object$control$response_scale == "log") density <- density - response[i, 1L]
        density
      }
      result <- result + value
    }
  }
  if (any(!is.finite(result))) .imr_abort("Non-finite held-out log likelihood; use training-fold refits.")
  result
}

#' Print IMR Cross-Validation Results
#'
#' @param x An `imr_cv` result.
#' @param digits Significant digits for scores.
#' @param ... Unused arguments are rejected.
#' @return `x`, invisibly.
#' @export
print.imr_cv <- function(x, digits = 3L, ...) {
  .imr_reject_dots(...)
  cat(toupper(x$control$model_variant), "cross-validation:", x$metric, "\n")
  cat(sprintf(
    "%d round(s) of %d folds; %s.\n", x$control$rounds, x$control$k,
    if (x$validation == "refit") "independent training-fold MCMC" else "PSIS reweighting of the full joint posterior"
  ))
  cat("Fitting algorithm:", .imr_sampler_description(x$control), "\n")
  cat("Pooled out-of-fold score:\n")
  print(signif(x$pooled, digits))
  cat("Mean fold score:\n")
  print(signif(x$fold_mean, digits))
  cat("Per-fold diagnostics are in $diagnostics; actual partitions are in $control$folds.\n")
  invisible(x)
}
