# A deterministic unit of refit CV once its row IDs and sampler seed are fixed.
# Keep all training-dependent preparation inside .imr_cv_refit(). The caller
# owns partition generation, RNG restoration and ordered metric aggregation.
.imr_cv_refit_task <- function(task, object, max_models, verbose) {
  tryCatch({
    if (verbose) cat(sprintf("CV round %d, fold %d\n", task$round, task$fold))
    .imr_cv_refit_fold(object, task$train_ids, task$test_ids, task$seed,
                       max_models, verbose)
  }, error = function(error) {
    .imr_abort(sprintf("CV round %d, fold %d failed: %s",
                       task$round, task$fold, conditionMessage(error)))
  })
}

.imr_cv_refit_fold <- function(object, train_ids, test_ids, seed,
                               max_models, verbose = FALSE) {
  preprocessing <- object$preprocessing
  dat <- preprocessing$input_data
  fitted <- .imr_cv_refit(object, train_ids, seed, verbose)
  test_platforms <- lapply(dat$platforms, function(x) {
    x[x$id %in% test_ids, , drop = FALSE]
  })
  present <- which(vapply(test_platforms, nrow, integer(1L)) > 0L)
  covariates <- if (!is.null(preprocessing$terms)) {
    preprocessing$formula_data[
      preprocessing$formula_data[[preprocessing$id]] %in% test_ids, , drop = FALSE]
  } else if (!is.null(dat$covariates)) {
    .imr_match_rows(dat$covariates, test_ids, "covariates")
  } else NULL
  if (length(fitted$model$covariate_names) == 0L) covariates <- NULL
  predicted <- do.call(rbind, stats::predict(fitted, test_platforms[present],
    platform_names = as.character(present), covariates = covariates,
    max_models = max_models, verbose = verbose))
  prediction <- predicted$prediction[match(test_ids, predicted$id)]
  if (any(!is.finite(prediction))) {
    .imr_abort("A training-fold model did not return finite predictions for every held-out subject.")
  }
  prediction
}
