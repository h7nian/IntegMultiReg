# Serial, isolated warm-up plus five measurements; no concurrent timing jobs.
# Arguments: output root, library, fixture, outcome, method, draws, burnin.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L, !dir.exists(args[1L]), dir.exists(args[2L]))
dir.create(args[1L], recursive = TRUE)
Sys.setenv(R_LIBS = normalizePath(args[2L]), OPENBLAS_NUM_THREADS = "1",
           OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
# Freeze the worker before launch: Rscript may read a source file incrementally.
# Repository edits during a long run must not alter an in-flight measurement.
worker <- file.path(normalizePath(args[1L]), "performance-run.R")
stopifnot(file.copy(".github/scripts/performance-run.R", worker))
rscript <- file.path(R.home("bin"), "Rscript")
for (index in 0:5) {
  name <- if (index == 0L) "warmup" else paste0("repeat-", index)
  output <- file.path(normalizePath(args[1L]), name)
  arguments <- c(shQuote(worker), shQuote(output), shQuote(args[3:7]))
  log <- file.path(args[1L], paste0(name, ".log"))
  if (Sys.info()[["sysname"]] == "Darwin") {
    status <- system2("/usr/bin/time", c("-l", shQuote(rscript), arguments),
                      stdout = log, stderr = log)
  } else {
    status <- system2(rscript, arguments, stdout = log, stderr = log)
  }
  if (status != 0L) stop("Measurement failed; inspect ", log)
  cat("Completed", name, "\n")
}
timings <- do.call(rbind, lapply(1:5, function(index) {
  data <- read.csv(file.path(args[1L], paste0("repeat-", index), "timings.csv"))
  data$repeat_index <- index
  data
}))
write.csv(timings, file.path(args[1L], "timings.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(timings, timings$stage), function(data) {
  data.frame(stage = data$stage[1L], median = median(data$elapsed),
             minimum = min(data$elapsed), maximum = max(data$elapsed))
}))
write.csv(summary, file.path(args[1L], "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
