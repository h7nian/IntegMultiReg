# Compare all eight completed chains; top-list overlap is descriptive, not a convergence test.
# Usage: Rscript chain-ranking-stability.R RUN_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, requireNamespace('coda', quietly = TRUE))
root <- normalizePath(args[1])
out <- args[2]
settings <- readRDS(file.path(root, 'settings.rds'))
stopifnot(length(settings$seeds) == 8L)
paths <- file.path(root, sprintf('chain-%02d-diagnostic.rds', 1:8))
for (i in 1:8) stopifnot(identical(readRDS(file.path(root,
  sprintf('chain-%02d.rds.status.rds', i)))$status, 'completed'))
fits <- lapply(paths, readRDS)
for (chain in 1:8) {
  fit <- fits[[chain]]
  stopifnot(identical(fit$model, fits[[1]]$model),
    identical(fit$control$mcmc, fits[[1]]$control$mcmc),
    fit$control$seed == settings$seeds[chain], fit$control$mcmc$draws == settings$draws)
  for (p in seq_along(fit$posterior$inclusion_probabilities)) {
    probabilities <- fit$posterior$inclusion_probabilities[[p]]
    stopifnot(identical(dimnames(probabilities),
        dimnames(fits[[1]]$posterior$inclusion_probabilities[[p]])),
      all(is.finite(probabilities)), all(probabilities >= 0 & probabilities <= 1))
  }
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- list()
ess <- list()
for (p in seq_along(fits[[1]]$posterior$inclusion_probabilities)) {
  platform <- fits[[1]]$model$platform_names[p]
  reference_prob <- fits[[1]]$posterior$inclusion_probabilities[[p]]
  groups <- fits[[1]]$model$subgroup_names[fits[[1]]$model$platform_subgroups[[p]]]
  k <- if (platform == 'mirna') 6L else 10L
  for (g in seq_len(nrow(reference_prob))) for (chain1 in 1:7) for (chain2 in (chain1 + 1):8) {
    prob1 <- fits[[chain1]]$posterior$inclusion_probabilities[[p]][g, ]
    prob2 <- fits[[chain2]]$posterior$inclusion_probabilities[[p]][g, ]
    top1 <- order(-prob1, seq_along(prob1))[1:k]
    top2 <- order(-prob2, seq_along(prob2))[1:k]
    rows[[length(rows) + 1L]] <- data.frame(platform, subgroup = groups[g],
     chain1 = chain1, chain2 = chain2,
     mean_abs_difference = mean(abs(prob1 - prob2)),
     max_abs_difference = max(abs(prob1 - prob2)),
     top_k = k, top_overlap = length(intersect(top1, top2)),
     boundary_tie_x = sum(prob1 == min(prob1[top1])) > sum(prob1[top1] == min(prob1[top1])),
     boundary_tie_y = sum(prob2 == min(prob2[top2])) > sum(prob2[top2] == min(prob2[top2])))
  }
  pairs <- do.call(cbind, lapply(2:nrow(reference_prob),
    function(upper) rbind(1:(upper - 1), upper)))
  for (chain in 1:8) for (pair in seq_len(ncol(pairs))) {
    interaction <- fits[[chain]]$posterior$interaction_draws[[p]][, pair]
    stopifnot(length(interaction) == settings$draws, all(is.finite(interaction)))
    ess_value <- as.numeric(coda::effectiveSize(coda::mcmc(interaction)))
    ess[[length(ess) + 1L]] <- data.frame(platform, subgroup1 = groups[pairs[1, pair]],
     subgroup2 = groups[pairs[2, pair]], chain = chain, retained = length(interaction),
     mean = mean(interaction), sd = sd(interaction), classical_ess = ess_value,
     mean_mcse = if (ess_value > 0) sd(interaction) / sqrt(ess_value) else NA_real_)
  }
}
rows <- do.call(rbind, rows)
ess <- do.call(rbind, ess)
write.csv(rows, file.path(out, 'mpip-differences.csv'), row.names = FALSE)
write.csv(ess, file.path(out, 'theta-ess.csv'), row.names = FALSE)
print(aggregate(cbind(mean_abs_difference, max_abs_difference)~platform + subgroup, rows, max), row.names = FALSE)
print(aggregate(cbind(top_overlap, top_k)~platform + subgroup, rows, min), row.names = FALSE)
print(aggregate(classical_ess~platform, ess, min), row.names = FALSE)
cat('ESS uses coda spectral estimation; top-list overlap is descriptive, with boundary ties flagged. All eight requested chains included; no convergence certification.\n')
saveRDS(list(input_hashes = tools::md5sum(paths),
  script_hash = tools::md5sum(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))),
  settings = settings, chains = 1:8, session = sessionInfo(),
  ranking = 'Descending mPIP, original feature order for ties; boundary ties explicitly flagged.',
  scope = 'All eight completed chains; per-chain classical spectral ESS differs from joint rank-normalized bulk/tail ESS. Neither numerical completion nor top-list overlap certifies convergence.'),
  file.path(out, 'ranking-stability-provenance.rds'))
