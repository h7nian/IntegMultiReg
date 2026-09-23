# Research diagnostics for matched IMR chains. No automatic convergence claim.
# coda's classical Gelman-Rubin PSRF is reported explicitly, not rank-normalized Rhat.
appendix_chain_diagnostics <- function(fits, out_dir, lag_max = 50L) {
  stopifnot(requireNamespace("coda", quietly = TRUE), length(fits) >= 2L,
            length(lag_max) == 1L, is.finite(lag_max), lag_max >= 1L)
  ref <- fits[[1L]]
  for (fit in fits) {
    stopifnot(identical(fit$model, ref$model),
              identical(fit$control$mcmc, ref$control$mcmc))
    for (p in seq_along(ref$posterior$inclusion_probabilities)) {
      a <- fit$posterior$inclusion_probabilities[[p]]
      b <- ref$posterior$inclusion_probabilities[[p]]
      stopifnot(identical(dim(a), dim(b)), identical(dimnames(a), dimnames(b)),
                all(is.finite(a)))
      theta <- fit$posterior$interaction_draws[[p]]
      stopifnot(is.matrix(theta), all(is.finite(theta)),
                nrow(theta) == ref$control$mcmc$draws,
                ncol(theta) == choose(nrow(a), 2L))
    }
  }
  stopifnot(ref$control$mcmc$draws >= 4L)
  correlations <- parameters <- acfs <- list()
  for (p in seq_along(ref$posterior$inclusion_probabilities)) {
    template <- ref$posterior$inclusion_probabilities[[p]]
    platform <- ref$model$platform_names[p]
    subgroups <- ref$model$subgroup_names[ref$model$platform_subgroups[[p]]]
    stopifnot(length(subgroups) == nrow(template))
    for (s in seq_len(nrow(template))) {
      for (i in seq_len(length(fits) - 1L)) for (j in seq.int(i + 1L, length(fits))) {
        a <- fits[[i]]$posterior$inclusion_probabilities[[p]][s, ]
        b <- fits[[j]]$posterior$inclusion_probabilities[[p]][s, ]
        defined <- length(a) > 1L && sd(a) > 0 && sd(b) > 0
        correlations[[length(correlations) + 1L]] <- data.frame(
          platform, subgroup = subgroups[s], chain1 = i, chain2 = j,
          n_features = length(a), pearson = if (defined) cor(a, b) else NA_real_,
          status = if (defined) "computed" else "undefined_constant_or_single_feature")
      }
    }
    if (nrow(template) < 2L) next
    # Native export order: (1,2), (1,3), (2,3), (1,4), ... .
    pairs <- do.call(cbind, lapply(seq.int(2L, nrow(template)),
                                  function(j) rbind(seq_len(j - 1L), j)))
    for (k in seq_len(ncol(pairs))) {
      draws <- lapply(fits, function(f) f$posterior$interaction_draws[[p]][, k])
      chains <- coda::mcmc.list(lapply(draws, coda::mcmc))
      nonconstant <- all(vapply(draws, sd, numeric(1L)) > 0)
      psrf <- if (nonconstant) coda::gelman.diag(chains, autoburnin = FALSE,
        multivariate = FALSE, transform = FALSE)$psrf[1L, ] else c(NA_real_, NA_real_)
      id <- data.frame(platform, subgroup1 = subgroups[pairs[1L, k]],
                       subgroup2 = subgroups[pairs[2L, k]])
      parameters[[length(parameters) + 1L]] <- cbind(id, data.frame(
        chains = length(fits), retained_per_chain = length(draws[[1L]]),
        classical_psrf = psrf[1L], classical_psrf_upper95 = psrf[2L],
        status = if (nonconstant) "computed" else "undefined_constant_chain"))
      for (i in seq_along(draws)) {
        nlag <- min(as.integer(lag_max), length(draws[[i]]) - 1L)
        values <- if (sd(draws[[i]]) > 0) as.numeric(stats::acf(draws[[i]],
          lag.max = nlag, plot = FALSE)$acf) else rep(NA_real_, nlag + 1L)
        acfs[[length(acfs) + 1L]] <- cbind(id, data.frame(chain = i,
          lag = seq.int(0L, nlag), acf = values))
      }
    }
  }
  result <- list(mpip_correlations = do.call(rbind, correlations),
                 theta_psrf = do.call(rbind, parameters), theta_acf = do.call(rbind, acfs))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result)) write.csv(result[[name]],
    file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  writeLines(c("Classical coda Gelman-Rubin PSRF; autoburnin=FALSE; transform=FALSE.",
    "Interaction matrices contain retained draws only; no second burn-in is discarded.",
    "Pearson correlations compare feature mPIPs within one platform and subgroup.",
    "Undefined correlations/PSRF are retained as NA, never treated as passing.",
    "Short chains test this exporter only. These tables do not establish convergence.",
    "Trace plots and scientific review remain necessary; historical starts are unavailable."),
    file.path(out_dir, "INTERPRETATION.txt"))
  capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
  invisible(result)
}
