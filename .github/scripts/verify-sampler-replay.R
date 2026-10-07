# Compile the supplied pre-optimization kernel with the current R toolchain.
# Compare complete retained states and RNG state, not just posterior averages.
library(IntegMultiReg)
source(".github/reference/native-bridge.R")

verify_replay <- function(source) {
  directory <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "sampler-replay-evidence")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  copied <- file.path(directory, "imr_replay_baseline.c")
  stopifnot(file.copy(source, copied, overwrite = TRUE))
  library <- file.path(directory, paste0("imr_replay_baseline", .Platform$dynlib.ext))
  status <- system2(file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "-o", shQuote(library), shQuote(copied)),
    stdout = file.path(directory, "compile.log"), stderr = file.path(directory, "compile.log"))
  if (status != 0L) stop("Could not compile baseline kernel; see compile.log.")
  dll <- dyn.load(library)
  on.exit(dyn.unload(library), add = TRUE)
  baseline <- getNativeSymbolInfo("imr_joint_sample", dll)$address
  reports <- list()
  for (case in c("continuous", "correlated-swap", "binary", "right.censored", "availability-imr", "availability-bms")) {
    availability <- startsWith(case, "availability-")
    groups <- readRDS(file.path(".github/reference/fixtures",
      paste0(if (availability) "availability" else case, "-data.rds")))
    for (chain in 1:4) for (thin in c(1L, 3L)) for (latent in c(FALSE, TRUE)) {
      args <- list(groups = groups,
        feature_platform = if (availability) c(1L, 2L) else c(1L, 1L),
        nu = if (availability) c(-1, -1.5) else -1,
        draws = 1000L, burnin = 200L, seed = 8100L + chain,
        model_variant = if (case == "availability-bms") "bms" else "imr",
        initial = c("empty", "full", "alternating", "empty")[chain],
        thin = thin, keep_latent = latent)
      old <- do.call(native_chain, c(args, list(kernel = baseline)))
      rng <- .Random.seed
      current <- do.call(native_chain, args)
      result <- data.frame(case, chain, thin, latent,
        samples_identical = identical(old, current, num.eq = FALSE),
        rng_identical = identical(rng, .Random.seed))
      reports[[length(reports) + 1L]] <- result
    }
  }
  report <- do.call(rbind, reports)
  write.csv(report, file.path(directory, "equivalence.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(directory, "sessionInfo.txt"))
  writeLines(capture.output(tools::md5sum(copied)), file.path(directory, "baseline-checksum.txt"))
  stopifnot(all(report$samples_identical), all(report$rng_identical))
  cat("Verified", nrow(report), "complete sampler/RNG replays with the current toolchain.\n")
}

arguments <- commandArgs(TRUE)
if (length(arguments) != 1L || !file.exists(arguments[[1L]])) {
  stop("Supply the path to the pre-optimization joint_sampler.c file.")
}
verify_replay(arguments[[1L]])
