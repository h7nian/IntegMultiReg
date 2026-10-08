# Compare public serial/parallel calls under the same numerical runtime.
# Set these limits before launching R; changing env vars after BLAS startup
# does not necessarily change an existing thread pool.
library(IntegMultiReg)
limits <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
if (any(Sys.getenv(limits) != "1")) {
  stop("Start R with the documented one-thread numerical environment for this exact replay.")
}
root <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "marginal-worker-evidence")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
id <- 1:20
platforms <- list(assay = data.frame(id, marker = sin(id)))
reports <- list()
for (type in c("continuous", "binary", "right.censored")) {
  outcome <- switch(type,
    continuous = data.frame(id, y = .7 + sin(id) + cos(id) / 3),
    binary = data.frame(id, y = rep(0:1, 10)),
    right.censored = data.frame(id, time = exp(.7 + sin(id)),
      status = rep(c(1, 1, 0), length.out = 20)))
  for (choice in c("variance", "coefficients", "coefficients_and_variance")) {
    fit <- function(workers) imr(platforms, outcome,
      outcome_type = type, marginalize = choice, min_subgroup_size = 0,
      priors = imr_priors(forced_scale = 1),
      mcmc = imr_mcmc(draws = 8, burnin = 4, chains = 2, seed = 58,
        workers = workers, keep_latent = TRUE, diagnostics = FALSE))
    set.seed(517)
    rng <- .Random.seed
    serial <- fit(1)
    parallel <- fit(2)
    reports[[length(reports) + 1L]] <- data.frame(type, choice,
      posterior_identical = identical(serial$posterior, parallel$posterior, num.eq = FALSE),
      initial_identical = identical(serial$control$initial, parallel$control$initial),
      rng_identical = identical(rng, .Random.seed))
  }
}
report <- do.call(rbind, reports)
write.csv(report, file.path(root, "equivalence.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
print(report)
stopifnot(all(report$posterior_identical), all(report$initial_identical), all(report$rng_identical))
