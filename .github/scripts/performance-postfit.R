# Isolate post-fit CV timing from the cost of fitting; one warm-up, five repeats.
# Arguments: frozen fit.rds, new output directory, legacy|importance.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, !dir.exists(args[2L]))
library(IntegMultiReg)
fit <- readRDS(args[1L])
method <- match.arg(args[3L], c("legacy", "importance"))
dir.create(args[2L], recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(args[2L], "sessionInfo.txt"))
times <- numeric(5L)
reference <- NULL
for (iteration in 0:5) {
  gc()
  elapsed <- system.time(result <- cv_imr(fit, k = 5L, rounds = 2L,
                                         cv_method = method))[["elapsed"]]
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
print(c(median = median(times), minimum = min(times), maximum = max(times)))
