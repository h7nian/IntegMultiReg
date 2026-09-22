# Post-fit modes use the full-fit transformed data and augmented responses.
# Reconstruct row IDs with the same subgroup builder used by imr().
.imr_cv_postfit <- function(object, k, rounds, max_models, verbose, cv_method,
                            workers = 1L, settings = .imr_cv_settings(cv_method),
                            supplied_folds = NULL) {
  workers <- .imr_check_integer_scalar(workers, "workers", min = 1)
  model <- object$model
  control <- object$control
  dat <- object$preprocessing$input_data
  if (as.double(sum(model$sample_sizes)) * rounds > .Machine$integer.max) {
    .imr_abort("The requested post-fit CV prediction matrix exceeds the native index limit.")
  }
  if (2 * control$priors$residual[["shape"]] + max(model$sample_sizes) >
      .Machine$integer.max) {
    .imr_abort("Residual prior shape exceeds the supported post-fit CV range.")
  }
  grouped <- subgroup_data(dat$outcome, dat$covariates, dat$platforms)
  response <- grouped[[1L]][model$subgroup_names]
  ids <- unlist(lapply(response, function(x) x[[1L]]), use.names = FALSE)
  subgroup <- rep(model$subgroup_names, model$sample_sizes)
  if (length(ids) != sum(model$sample_sizes)) {
    .imr_abort("Stored inputs do not match the fitted subgroup rows.")
  }
  max_models <- min(max_models, control$mcmc$draws)
  partitions <- .imr_cv_fold_matrices(supplied_folds, ids, rounds)
  result <- if (workers == 1L) {
    # Current scores are computed below; skip unused historical scores.
    .imr_call_cv_postfit_native(object, k, rounds, max_models,
                                verbose, settings$model_set == "draws",
                                stage = if (settings$score_method == "standard") "predict" else "full",
                                settings = settings, folds = partitions$folds,
                                row_order = partitions$row_order)
  } else {
    .imr_cv_postfit_parallel(object, k, rounds, max_models, verbose,
                              settings, workers, partitions)
  }
  labels <- c(model$subgroup_names, "all")
  colnames(result$total_cindex) <- colnames(result$subset_cindex) <- labels
  records <- do.call(rbind, lapply(seq_len(rounds), function(round) {
    data.frame(round = round, fold = result$folds[, round], id = ids,
               subgroup = subgroup, prediction = result$predictions[, round],
               row.names = NULL)
  }))
  if (any(!is.finite(records$prediction))) {
    .imr_abort("Post-fit CV produced non-finite predictions; inspect the fit and importance weights.")
  }
  # Scoring is independent of the model collection and predictive density.
  # Current rules preserve NA for undefined folds, as in refit.
  if (settings$score_method == "standard") {
    outcome <- .imr_match_rows(dat$outcome, ids, "outcome")
    groups <- c(lapply(model$subgroup_names, function(g) which(subgroup == g)),
                list(seq_along(ids)))
    for (round in seq_len(rounds)) {
      prediction <- result$predictions[, round]
      folds <- result$folds[, round]
      score <- function(index) .imr_cv_accuracy(
        control$outcome_type, prediction[index], outcome[index, , drop = FALSE])
      result$total_cindex[round, ] <- vapply(groups, score, numeric(1))
      result$subset_cindex[round, ] <- vapply(groups, function(index) {
        mean(vapply(seq_len(k), function(fold) score(index[folds[index] == fold]),
                    numeric(1)))
      }, numeric(1))
    }
  }
  list(pooled = result$total_cindex, fold_mean = result$subset_cindex,
       predictions = records,
       metric = switch(control$outcome_type, right.censored = "C-index",
                       binary = "AUC", continuous = "MSE"),
       validation = cv_method,
       control = .imr_cv_control(settings, object, k, rounds, max_models,
         do.call(rbind, lapply(seq_len(rounds), function(round) {
           data.frame(id = ids, round = round, fold = result$folds[, round],
                      row_order = result$row_order[, round])
         })), if (is.null(supplied_folds)) "gsl" else "supplied"))
}
