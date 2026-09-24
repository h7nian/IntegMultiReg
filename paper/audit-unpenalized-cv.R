# Diagnose a failed strict paper CV without substituting another estimator.
# Internal fold planning is used only by this version-pinned research audit.
audit_unpenalized_failure <- function(path, settings, error) {
  args <- settings$reference_arguments$cv
  checkpoint <- paste0(path, ".fit.rds")
  if (!identical(args$ridge, 0) || !identical(args$model_set, "draws") ||
      !file.exists(checkpoint) ||
      !grepl("Cholesky|singular|positive definite", conditionMessage(error), ignore.case = TRUE))
    stop(error)
  fit <- experiment_fit_checkpoint(checkpoint, stop("Missing completed fit"))
  controls <- do.call(IntegMultiReg:::.imr_cv_settings, args)
  plan <- IntegMultiReg:::.imr_call_cv_postfit_native(fit, settings$k,
    settings$rounds, 100L, FALSE, TRUE, stage = "plan", settings = controls)
  reports <- list()
  example <- NULL
  offset <- 0L
  for (g in seq_along(fit$model$subgroup_names)) {
    n <- fit$model$sample_sizes[g]
    platforms <- fit$model$subgroup_platforms[[g]]
    masks <- lapply(fit$posterior$selection_draws, function(state) {
      c(rep(TRUE, 1L + length(fit$model$covariate_names)),
        unlist(lapply(platforms, function(p) {
          state[[p]][match(g, fit$model$platform_subgroups[[p]]), ] == 1
        }), use.names = FALSE))
    })
    columns <- vapply(masks, sum, 0)
    for (round in seq_len(settings$rounds)) for (fold in seq_len(settings$k)) {
      train <- plan$folds[offset + seq_len(n), round] != fold
      bad <- which(columns > sum(train))
      reports[[length(reports) + 1L]] <- data.frame(
        subgroup = fit$model$subgroup_names[g], round, fold,
        training_rows = sum(train), states = length(columns),
        dimensionally_singular = length(bad), fraction = length(bad) / length(columns))
      if (length(bad) && is.null(example)) {
        X <- IntegMultiReg:::.imr_posterior_design(fit, g)
        X <- X[train, masks[[bad[1L]]], drop = FALSE]
        example <- list(subgroup = fit$model$subgroup_names[g], round = round,
          fold = fold, state = bad[1L], design = X, rank = qr(X)$rank)
        stopifnot(example$rank < ncol(X))
      }
    }
    offset <- offset + n
  }
  if (is.null(example)) stop(error) # unexplained numerical failures block acceptance
  report <- do.call(rbind, reports)
  write.csv(report, paste0(path, ".rank-audit.csv"), row.names = FALSE)
  saveRDS(list(error = conditionMessage(error), report = report, example = example,
    note = "Structural failure retained. Example is not asserted to be the first native failing state."),
    paste0(path, ".rank-audit.rds"))
  cat("AUDITED structural rank deficiency:", path, "\n")
  invisible(NULL)
}
