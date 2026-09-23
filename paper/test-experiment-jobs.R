script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "original-experiment-helpers.R"))
directory <- tempfile("job-records-")
dir.create(directory)
path <- file.path(directory, "result.rds")
checkpoint <- paste0(path, ".fit.rds")
calls <- 0L
saved_fit <- experiment_fit_checkpoint(checkpoint, {calls <- calls + 1L; list(value = 42L)})
tryCatch(record_experiment_job(path, stop("CV failed after fitting")), error = identity)
restored_fit <- experiment_fit_checkpoint(checkpoint, stop("Fit must not rerun"))
stopifnot(identical(saved_fit, restored_fit), calls == 1L)
# The following job-history checks use a separate result path.
path <- file.path(directory, "record-test.rds")
saveRDS("incomplete output", path)
stopifnot(!experiment_job_complete(path))
set.seed(53)
rng <- .Random.seed
observed <- character()
value <- withCallingHandlers(record_experiment_job(path, {
  warning("test convergence warning")
  42L
}), warning = function(w) {
  observed <<- c(observed, conditionMessage(w))
  invokeRestart("muffleWarning")
})
status <- readRDS(paste0(path, ".status.rds"))
stopifnot(experiment_job_complete(path))
stopifnot(identical(value, 42L), identical(.Random.seed, rng),
  identical(observed, "test convergence warning"), status$status == "completed",
  length(status$warnings) == 1L, status$ended >= status$started,
  status$elapsed[["elapsed"]] >= 0, !is.null(status$heap_end))
error <- tryCatch(record_experiment_job(path, stop("test failure")), error = identity)
status <- readRDS(paste0(path, ".status.rds"))
stopifnot(!experiment_job_complete(path))
stopifnot(inherits(error, "error"), status$status == "failed",
  status$error$message == "test failure", status$attempt == 2L,
  readRDS(paste0(path, ".attempt-1.rds"))$status == "completed")
cat("PASS: warning propagation, timing, RNG preservation, failure and attempt history\n")
unlink(directory, recursive = TRUE)

# Resume reads the small fold sidecar, not a potentially multi-GiB result.
path <- tempfile(fileext = ".rds")
folds <- data.frame(id = 1:4, round = 1L, fold = c(1L, 2L, 1L, 2L))
saveRDS(list(cv = list(control = list(folds = folds))), path)
stopifnot(identical(experiment_result_folds(path), folds))
saveRDS(folds, paste0(path, ".folds.rds"))
writeLines("deliberately unreadable result payload", path)
stopifnot(identical(experiment_result_folds(path), folds))
unlink(c(path, paste0(path, ".folds.rds")))
cat("PASS: fold sidecar avoids large result reads; old results remain compatible\n")
