# Batch independent folds so each worker indexes full-fit states only once.
# The GSL stream is replayed locally to retain both fold membership and the
# historical within-fold row order; every replay is checked against the plan.
.imr_cv_postfit_parallel <- function(object, k, rounds, max_models, verbose,
                                      settings, workers, partitions) {
  plan <- .imr_call_cv_postfit_native(object, k, rounds, max_models, FALSE,
                                      settings$model_set == "draws", stage = "plan",
                                      settings = settings, folds = partitions$folds,
                                      row_order = partitions$row_order)
  n_tasks <- k * rounds
  batches <- unname(split(seq_len(n_tasks),
                          rep(seq_len(min(workers, n_tasks)), length.out = n_tasks)))
  values <- .imr_cv_map(batches, .imr_cv_postfit_task, workers,
                        object = object, k = k, rounds = rounds,
                        max_models = max_models, verbose = verbose,
                        settings = settings, partitions = partitions,
                        folds = plan$folds, row_order = plan$row_order)
  predictions <- plan$predictions
  filled <- matrix(FALSE, nrow(predictions), ncol(predictions))
  for (batch in seq_along(batches)) {
    active <- .imr_cv_postfit_active(batches[[batch]], k, rounds, plan$folds)
    if (any(filled[active]) || length(values[[batch]]) != sum(active))
      .imr_abort("Post-fit worker results contain overlapping or missing tasks.")
    predictions[active] <- values[[batch]]
    filled[active] <- TRUE
  }
  if (!all(filled)) .imr_abort("Post-fit workers did not cover every held-out prediction.")
  if (settings$score_method == "standard") {
    # R's existing public scoring path replaces both native metric matrices.
    plan$predictions <- predictions
    return(plan)
  }
  result <- .imr_call_cv_postfit_native(object, k, rounds, max_models, verbose,
                                        settings$model_set == "draws", stage = "score",
                                        predictions = predictions, settings = settings,
                                        folds = partitions$folds, row_order = partitions$row_order)
  if (!identical(result$folds, plan$folds))
    .imr_abort("Post-fit scoring did not reproduce the planned partitions.")
  result
}

.imr_cv_postfit_active <- function(batch, k, rounds, folds) {
  tasks <- matrix(FALSE, k, rounds)
  tasks[batch] <- TRUE
  active <- matrix(FALSE, nrow(folds), ncol(folds))
  for (round in seq_len(rounds)) active[, round] <- tasks[folds[, round], round]
  active
}

.imr_cv_postfit_task <- function(batch, object, k, rounds, max_models, verbose,
                                 settings, partitions, folds, row_order) {
  tryCatch({
    tasks <- matrix(FALSE, k, rounds)
    tasks[batch] <- TRUE
    result <- .imr_call_cv_postfit_native(object, k, rounds, max_models, verbose,
                                          settings$model_set == "draws", stage = "predict", tasks = tasks,
                                          settings = settings, folds = partitions$folds,
                                          row_order = partitions$row_order)
    if (!identical(result$folds, folds) || !identical(result$row_order, row_order))
      .imr_abort("Worker did not reproduce the planned partitions.")
    active <- .imr_cv_postfit_active(batch, k, rounds, folds)
    bad <- which(active & !is.finite(result$predictions), arr.ind = TRUE)
    if (nrow(bad)) {
      round <- bad[1L, 2L]
      .imr_abort(sprintf("Non-finite prediction at round %d, fold %d.",
                         round, folds[bad[1L, 1L], round]))
    }
    if (any(!is.na(result$predictions[!active])))
      .imr_abort("Worker computed a prediction outside its assigned folds.")
    # Return only the assigned cells, in the native matrix's column-major order.
    result$predictions[active]
  }, error = function(error) {
    labels <- sprintf("round %d/fold %d", (batch - 1L) %/% k + 1L,
                       (batch - 1L) %% k + 1L)
    .imr_abort(sprintf("Post-fit worker assigned %s failed: %s",
                       paste(labels, collapse = ", "), conditionMessage(error)))
  })
}
