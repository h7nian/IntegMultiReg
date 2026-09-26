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
for (i in 1:8) {
  f <- fits[[i]]
  stopifnot(identical(f$model, fits[[1]]$model),
    identical(f$control$mcmc, fits[[1]]$control$mcmc),
    f$control$seed == settings$seeds[i], f$control$mcmc$draws == settings$draws)
  for (p in seq_along(f$posterior$inclusion_probabilities)) {
    m <- f$posterior$inclusion_probabilities[[p]]
    stopifnot(identical(dimnames(m), dimnames(fits[[1]]$posterior$inclusion_probabilities[[p]])),
      all(is.finite(m)), all(m >= 0 & m <= 1))
  }
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- list()
 ess <- list()
for (p in seq_along(fits[[1]]$posterior$inclusion_probabilities)) {
  platform <- fits[[1]]$model$platform_names[p]
  a <- fits[[1]]$posterior$inclusion_probabilities[[p]]
  groups <- fits[[1]]$model$subgroup_names[fits[[1]]$model$platform_subgroups[[p]]]
  k <- if (platform == 'mirna') 6L else 10L
  for (g in seq_len(nrow(a))) for (i in 1:7) for (j in (i + 1):8) {
    x <- fits[[i]]$posterior$inclusion_probabilities[[p]][g, ]
    y <- fits[[j]]$posterior$inclusion_probabilities[[p]][g, ]
    ix <- order(-x, seq_along(x))[1:k]
    iy <- order(-y, seq_along(y))[1:k]
    rows[[length(rows) + 1L]] <- data.frame(platform, subgroup = groups[g], chain1 = i, chain2 = j,
     mean_abs_difference = mean(abs(x - y)), max_abs_difference = max(abs(x - y)),
     top_k = k, top_overlap = length(intersect(ix, iy)),
     boundary_tie_x = sum(x == min(x[ix])) > sum(x[ix] == min(x[ix])),
     boundary_tie_y = sum(y == min(y[iy])) > sum(y[iy] == min(y[iy])))
  }
  pairs <- do.call(cbind, lapply(2:nrow(a), function(j)rbind(1:(j - 1), j)))
  for (i in 1:8) for (j in seq_len(ncol(pairs))) {
    d <- fits[[i]]$posterior$interaction_draws[[p]][, j]
    stopifnot(length(d) == settings$draws, all(is.finite(d)))
    n <- as.numeric(coda::effectiveSize(coda::mcmc(d)))
    ess[[length(ess) + 1L]] <- data.frame(platform, subgroup1 = groups[pairs[1, j]],
     subgroup2 = groups[pairs[2, j]], chain = i, retained = length(d),
     mean = mean(d), sd = sd(d), classical_ess = n, mean_mcse = if (n > 0) sd(d) / sqrt(n) else NA_real_)
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
