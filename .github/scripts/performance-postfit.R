# Isolate post-fit CV timing from the cost of fitting; one warm-up, five repeats.
# Arguments: frozen fit.rds, new output directory, all_draws|top_unique, [workers].
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(3L, 4L), !file.exists(args[2L]))
library(IntegMultiReg)
fit <- upgrade_imr_object(readRDS(args[1L]))
model_set <- match.arg(args[3L], c("all_draws", "top_unique"))
workers <- if (length(args) == 4L) as.integer(args[4L]) else 1L
stopifnot(!is.na(workers), workers >= 1L)
controls <- list(object = fit, k = 5L, rounds = 2L, cv_method = "reweight", model_set = model_set)
if ("workers" %in% names(formals(cv_imr))) {
  controls$workers <- workers
} else if (workers != 1L) stop("Installed baseline does not support workers")
dir.create(args[2L], recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(args[2L], "sessionInfo.txt"))
times <- numeric(5L)
reference <- NULL
for (iteration in 0:5) {
  gc()
  elapsed <- system.time(result <- do.call(cv_imr, controls))[["elapsed"]]
  if (iteration == 0L) {
    reference <- result
    saveRDS(result, file.path(args[2L], "result.rds"))
  } else {
    stopifnot(identical(reference, result))
    times[iteration] <- elapsed
  }
}
write.csv(data.frame(repeat_index = 1:5, elapsed = times),
          file.path(args[2L], "timings.csv"), row.names = FALSE)
saveRDS(list(source_sha = Sys.getenv("IMR_SOURCE_SHA", unset = NA_character_),
             input_md5 = tools::md5sum(args[1L]), workers = workers, cv_method = "reweight", model_set = model_set),
        file.path(args[2L], "provenance.rds"))
print(c(median = median(times), minimum = min(times), maximum = max(times)))
