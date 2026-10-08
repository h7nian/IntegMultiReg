# Preserve the frozen R transition and RNG semantics across numerical policies.
library(IntegMultiReg)
reference_root <- file.path("tests", "reference", "marginal-r")
reference <- new.env(parent = asNamespace("IntegMultiReg"))
for (file in c("helpers.R", "variance.R", "coefficients.R")) {
  sys.source(file.path(reference_root, file), reference)
}
output <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "native-marginal-evidence")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
original_rng <- RNGkind()
original_options <- options("matprod")
reports <- list()
for (policy in c("default", "internal", "blas")) {
  options(matprod = policy)
  for (normal in c("Inversion", "Box-Muller")) {
    RNGkind("Mersenne-Twister", normal, "Rejection")
    for (choice in c("variance", "coefficients")) {
      name <- paste0(".imr_", choice, "_marginal_sample")
      for (case in c("continuous", "correlated-swap", "binary", "right.censored", "availability-imr", "availability-bms")) {
        availability <- startsWith(case, "availability")
        groups <- readRDS(file.path(reference_root, paste0(if (availability) "availability" else case, "-data.rds")))
        arguments <- list(
          groups = groups, feature_platform = if (availability) c(1L, 2L) else c(1L, 1L),
          nu = if (availability) c(-1, -1.5) else -1,
          draws = 300L, burnin = 30L, thin = 2L, seed = 923L,
          model_variant = if (case == "availability-bms") "bms" else "imr",
          initial = "full", keep_latent = TRUE, swap_rate = 1.5
        )
        expected <- do.call(get(name, reference), arguments)
        expected_rng <- .Random.seed
        actual <- do.call(get(name, asNamespace("IntegMultiReg")), arguments)
        samples_identical <- identical(expected, actual, num.eq = FALSE)
        rng_identical <- identical(expected_rng, .Random.seed)
        row <- data.frame(policy, normal, choice, case, samples_identical, rng_identical)
        print(row)
        reports[[length(reports) + 1L]] <- row
        if (!samples_identical || !rng_identical) {
          saveRDS(
            list(arguments = arguments, expected = expected, actual = actual),
            file.path(output, paste(policy, normal, choice, case, "failure.rds", sep = "-"))
          )
        }
      }
    }
  }
}
options(original_options)
do.call(RNGkind, as.list(original_rng))
result <- do.call(rbind, reports)
write.csv(result, file.path(output, "replay.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
stopifnot(all(result$samples_identical), all(result$rng_identical))
