# Complement the historical classical PSRF with posterior's rank-based metrics.
# Usage: Rscript rank-chain-diagnostics.R RUN_DIRECTORY OUTPUT_DIRECTORY [CHAINS]
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 2:3, requireNamespace('posterior', quietly = TRUE))
root <- normalizePath(args[1])
out <- args[2]
nchains <- if (length(args) == 3L) as.numeric(args[3]) else 8L
stopifnot(length(nchains) == 1L, is.finite(nchains), nchains == as.integer(nchains), nchains >= 2L, nchains <= 8L)
paths <- file.path(root, sprintf('chain-%02d-diagnostic.rds', seq_len(nchains)))
for (i in seq_len(nchains))stopifnot(identical(readRDS(file.path(root,
 sprintf('chain-%02d.rds.status.rds', i)))$status, 'completed'))
fits <- lapply(paths, readRDS)
ref <- fits[[1]]
for (f in fits)stopifnot(identical(f$model, ref$model), identical(f$control$mcmc, ref$control$mcmc))
rows <- list()
for (p in seq_along(ref$posterior$interaction_draws)) {
  groups <- ref$model$subgroup_names[ref$model$platform_subgroups[[p]]]
  if (length(groups) < 2L) next
  pairs <- do.call(cbind, lapply(2:length(groups), function(j)rbind(seq_len(j - 1), j)))
  for (j in seq_len(ncol(pairs))) {
    draws <- vapply(fits, function(f)f$posterior$interaction_draws[[p]][, j],
     numeric(ref$control$mcmc$draws))
    stopifnot(is.matrix(draws), ncol(draws) == nchains, all(is.finite(draws)))
    rows[[length(rows) + 1L]] <- data.frame(platform = ref$model$platform_names[p],
     subgroup1 = groups[pairs[1, j]], subgroup2 = groups[pairs[2, j]],
     chains = nchains, retained_per_chain = nrow(draws),
     rank_normalized_split_rhat = posterior::rhat(draws),
     bulk_ess = posterior::ess_bulk(draws), tail_ess = posterior::ess_tail(draws),
     mean_mcse = posterior::mcse_mean(draws))
  }
}
result <- do.call(rbind, rows)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(out, 'theta-rank-diagnostics.csv'), row.names = FALSE)
saveRDS(list(result = result, input_hashes = tools::md5sum(paths),
 script_hash = tools::md5sum(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE))),
 posterior_version = as.character(utils::packageVersion('posterior')), session = sessionInfo(),
 chains = seq_len(nchains), interim = nchains != 8L,
 interpretation = 'Rank-normalized split/folded Rhat and bulk/tail ESS supplement classical PSRF. All retained draws are used without another burn-in or thinning. Constant-chain results remain NA. These diagnostics alone cannot certify convergence.',
 sources = c('https://mc-stan.org/posterior/reference/rhat.html',
 'https://mc-stan.org/posterior/reference/ess_bulk.html',
 'https://mc-stan.org/posterior/reference/ess_tail.html')),
 file.path(out, 'theta-rank-diagnostics.rds'))
print(result, row.names = FALSE)
cat(if (nchains == 8L) 'All eight requested chains included.\n' else 'INTERIM: requested eight-chain study is incomplete.\n')
