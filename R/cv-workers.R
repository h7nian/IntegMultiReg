# Task order is fixed before dispatch. Static chunks transmit the fit once per
# worker, and parLapply returns results in input order, not completion order.
.imr_cv_map <- function(tasks, evaluate, workers, ...) {
  workers <- .imr_check_integer_scalar(workers, "workers", min = 1)
  if (!length(tasks)) return(list())
  workers <- min(workers, length(tasks))
  if (workers == 1L) return(lapply(tasks, evaluate, ...))

  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  # Children inherit these limits at process startup, before loading a math
  # library. Restore the parent's environment on success and on every error.
  thread_variables <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                        "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                        "BLIS_NUM_THREADS")
  thread_settings <- Sys.getenv(thread_variables, unset = NA_character_)
  on.exit({
    missing <- is.na(thread_settings)
    Sys.unsetenv(thread_variables[missing])
    if (any(!missing)) do.call(Sys.setenv, as.list(thread_settings[!missing]))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(thread_variables)),
                                            thread_variables)))
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  package_path <- getNamespaceInfo(asNamespace("IntegMultiReg"), "path")
  library_paths <- unique(c(dirname(package_path), .libPaths()))
  initialize <- function(paths, expected_path, rng_kind, model_options) {
    .libPaths(paths)
    do.call(RNGkind, as.list(rng_kind))
    options(model_options)
    namespace <- loadNamespace("IntegMultiReg")
    if (!identical(normalizePath(getNamespaceInfo(namespace, "path")),
                   normalizePath(expected_path))) {
      stop("CV worker loaded a different IntegMultiReg installation.")
    }
    NULL
  }
  # Bootstrap must deserialize before the package's library path is known.
  # A base-only closure also avoids serializing the parent's task/cluster frame.
  environment(initialize) <- baseenv()
  parallel::clusterCall(cluster, initialize, library_paths, package_path,
                        RNGkind(), options()[c("contrasts", "na.action")])
  results <- parallel::parLapply(cluster, tasks, .imr_cv_worker_result,
                                 evaluate = evaluate, ...)
  lapply(results, function(result) {
    for (condition in result$warnings) warning(condition)
    result$value
  })
}

# Child stderr is not a reliable diagnostic channel. Relay conditions through
# the result transport; the parent emits them in the original task order.
.imr_cv_worker_result <- function(task, evaluate, ...) {
  warnings <- list()
  value <- withCallingHandlers(evaluate(task, ...), warning = function(condition) {
    warnings[[length(warnings) + 1L]] <<- condition
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
