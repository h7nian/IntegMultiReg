# Development-only timing for the staged, internal refit task dispatcher.
# Run in a fresh process with single-threaded math libraries. Do not run other
# benchmarks concurrently. This is not a paper experiment or a public API.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: fit.rds new-output-directory folds rounds")
library(IntegMultiReg)
fit <- readRDS(args[[1L]])
validate_imr(fit)
output <- args[[2L]]
if (file.exists(output)) stop("Refusing to overwrite benchmark evidence")
dir.create(output, recursive = TRUE)
k <- as.integer(args[[3L]])
rounds <- as.integer(args[[4L]])
stopifnot(!is.na(k), k >= 2L, !is.na(rounds), rounds >= 1L)
reference <- NULL
timings <- list()
for (workers in c(1L, 2L)) {
  for (iteration in 0:5) {
    gc()
    set.seed(739)
    rng <- .Random.seed
    elapsed <- system.time({
      result <- IntegMultiReg:::.imr_cv_refit_result(
        fit, k, rounds, 100L, FALSE, workers)
    })[["elapsed"]]
    stopifnot(identical(.Random.seed, rng))
    if (is.null(reference)) reference <- result
    stopifnot(identical(result, reference))
    timings[[length(timings) + 1L]] <- data.frame(
      workers = workers, iteration = iteration, seconds = elapsed)
    cat(workers, iteration, elapsed, "EXACT\n")
  }
}
write.csv(do.call(rbind, timings), file.path(output, "timings.csv"), row.names = FALSE)
saveRDS(reference, file.path(output, "result.rds"))
saveRDS(list(input_md5 = tools::md5sum(args[[1L]]), folds = k, rounds = rounds,
             source_sha = Sys.getenv("IMR_SOURCE_SHA", unset = NA_character_),
             library = getNamespaceInfo(asNamespace("IntegMultiReg"), "path")),
        file.path(output, "provenance.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
