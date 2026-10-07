# Independent small-model integrals are the reference, not an older sampler.
library(IntegMultiReg)
source(".github/reference/native-bridge.R")
root <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "joint-reference-evidence")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
read <- function(name) readRDS(file.path(".github/reference/fixtures", paste0(name, ".rds")))
reports <- list()
for (case in c("continuous", "correlated-swap", "binary", "right.censored", "availability-imr", "availability-bms")) {
  missing_platform <- startsWith(case, "availability-")
  groups <- read(paste0(if (missing_platform) "availability" else case, "-data"))
  observed <- case %in% c("binary", "right.censored")
  if (observed) {
    quadrature <- read(paste0(case, "-quadrature"))
    stopifnot(quadrature$max_difference < .001)
    reference <- quadrature$fine
  } else {
    reference <- read(paste0(case, "-reference"))
  }
  variant <- if (case == "availability-bms") "bms" else "imr"
  features <- if (missing_platform) c(1L, 2L) else c(1L, 1L)
  nu <- if (missing_platform) c(-1, -1.5) else -1
  elapsed <- system.time(chains <- lapply(1:4, function(i) {
    native_chain(groups, features, nu,
      seed = 8100L + i, model_variant = variant,
      initial = c("empty", "full", "alternating", "empty")[i]
    )
  }))
  free <- if (missing_platform) reference$free else seq_len(4L)
  values <- lapply(chains, function(x) {
    cbind(
      x$selection[, free], do.call(cbind, x$beta),
      do.call(cbind, lapply(x$beta, function(b) b^2)), x$variance, x$theta,
      if (observed) do.call(cbind, x$latent)
    )
  })
  target <- c(
    reference$inclusion, unlist(reference$beta_mean), unlist(reference$beta_second),
    reference$variance_mean, if (variant == "imr") reference$theta_mean,
    if (observed) unlist(reference$latent_mean)
  )
  stopifnot(length(target) == ncol(values[[1]]))
  batch <- do.call(rbind, lapply(values, function(m) {
    do.call(rbind, lapply(
      1:60,
      function(b) colMeans(m[(b - 1) * 500 + 1:500, , drop = FALSE])
    ))
  }))
  mcse <- apply(batch, 2, sd) / sqrt(nrow(batch))
  rhat <- vapply(seq_along(target), function(j) {
    posterior::rhat(vapply(values, function(m) m[, j], numeric(nrow(values[[1]]))))
  }, 0)
  estimate <- colMeans(do.call(rbind, values))
  moments <- data.frame(
    parameter = seq_along(target), target, estimate, mcse, rhat,
    error = estimate - target
  )
  states <- do.call(rbind, lapply(chains, function(x) x$selection[, free, drop = FALSE]))
  frequency <- tabulate(1L + drop(states %*% 2^(seq_along(free) - 1)), nbins = nrow(reference$states)) / nrow(states)
  tv <- sum(abs(frequency - reference$probability)) / 2
  passed <- all(abs(moments$error) < .01 + 5 * mcse) && max(rhat, na.rm = TRUE) < 1.01 && tv < .02
  out <- file.path(root, case)
  dir.create(out, showWarnings = FALSE)
  saveRDS(chains, file.path(out, "chains.rds"))
  write.csv(moments, file.path(out, "moments.csv"), row.names = FALSE)
  write.csv(data.frame(reference$states, reference = reference$probability, sampled = frequency),
    file.path(out, "states.csv"),
    row.names = FALSE
  )
  reports[[case]] <- data.frame(case,
    total_variation = tv, max_rhat = max(rhat, na.rm = TRUE),
    passed, seconds = elapsed[["elapsed"]]
  )
  print(reports[[case]])
}
write.csv(do.call(rbind, reports), file.path(root, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
stopifnot(all(vapply(reports, function(x) x$passed, TRUE)))

source(".github/scripts/verify-joint-quantiles.R")
