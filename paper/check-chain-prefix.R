# Compare overlapping trajectories without reconstructing all selection matrices.
# This establishes deterministic prefix agreement, not sampler convergence.
check_chain_prefix <- function(old, new) {
  stopifnot(inherits(old, 'imr_experiment_fit_checkpoint_v1'),
            inherits(new, 'imr_experiment_fit_checkpoint_v1'))
  a <- old$fit
   b <- new$fit
  n <- a$control$mcmc$draws
   m <- b$control$mcmc$draws
  burn <- a$control$mcmc$burnin
  stopifnot(m > n, identical(burn, b$control$mcmc$burnin),
            identical(a$model, b$model))
  # Calls and final RNG states necessarily differ after a longer run.
  keys <- setdiff(names(a$control), c('call', 'mcmc', 'rng_state'))
  stopifnot(identical(keys, setdiff(names(b$control), c('call', 'mcmc', 'rng_state'))),
            identical(a$control[keys], b$control[keys]),
            length(old$layout) == n, length(new$layout) == m,
            identical(old$layout, new$layout[seq_len(n)]),
            length(a$posterior$log_posterior) == n + burn,
            length(b$posterior$log_posterior) == m + burn,
            all(is.finite(a$posterior$log_posterior)),
            all(is.finite(b$posterior$log_posterior)),
            identical(a$posterior$log_posterior,
                      b$posterior$log_posterior[seq_len(n + burn)]),
            length(old$pools) == length(new$pools),
            length(old$indices) == length(old$pools),
            length(new$indices) == length(new$pools),
            length(a$posterior$interaction_draws) == length(b$posterior$interaction_draws))
  for (p in seq_along(a$posterior$interaction_draws)) {
    x <- a$posterior$interaction_draws[[p]]
    y <- b$posterior$interaction_draws[[p]]
    stopifnot(nrow(x) == n, nrow(y) == m, all(is.finite(x)), all(is.finite(y)),
              identical(x, y[seq_len(n), , drop = FALSE]))
  }
  for (p in seq_along(old$pools)) {
    i <- old$indices[[p]]
     j <- new$indices[[p]]
    stopifnot(length(i) == n, length(j) == m, !anyNA(i), !anyNA(j),
              all(i >= 1L & i <= length(old$pools[[p]])),
              all(j >= 1L & j <= length(new$pools[[p]])))
    pairs <- unique(data.frame(i = i, j = j[seq_len(n)]))
    # Pool numbering is an implementation detail; compare actual matrices.
    stopifnot(all(vapply(seq_len(nrow(pairs)), function(k)
      identical(old$pools[[p]][[pairs$i[k]]], new$pools[[p]][[pairs$j[k]]]), logical(1))))
  }
  data.frame(seed = a$control$seed, original_retained = n, longer_retained = m,
             burnin = burn, exact_log_posterior = TRUE, exact_interaction = TRUE,
             exact_selection = TRUE, exact_initial = TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 4L) stop('Usage: Rscript check-chain-prefix.R OLD_ROOT NEW_ROOT OUT_DIR CHAIN')
  chain <- suppressWarnings(as.integer(args[4L]))
  stopifnot(!is.na(chain), chain %in% 1:8, as.character(chain) == args[4L])
  roots <- args[1:2]
   settings <- lapply(file.path(roots, 'settings.rds'), readRDS)
  a <- settings[[1]]
   b <- settings[[2]]
  stopifnot(b$draws > a$draws,
    identical(a[c('burnin', 'seeds', 'sampler_method', 'nu', 'quick')],
              b[c('burnin', 'seeds', 'sampler_method', 'nu', 'quick')]),
    identical(unname(a$data_hash), unname(b$data_hash)),
    identical(unname(a$source_hashes), unname(b$source_hashes)))
  for (s in settings) stopifnot(identical(tools::md5sum(names(s$data_hash)), s$data_hash),
    identical(tools::md5sum(names(s$source_hashes)), s$source_hashes))
  paths <- file.path(roots, sprintf('chain-%02d.rds', chain))
  statuses <- paste0(paths, '.status.rds')
  stopifnot(all(vapply(statuses, function(p) identical(readRDS(p)$status, 'completed'), logical(1))))
  result <- check_chain_prefix(readRDS(paths[1]), readRDS(paths[2]))
  stopifnot(result$seed == a$seeds[chain], result$original_retained == a$draws,
            result$longer_retained == b$draws, result$burnin == a$burnin)
  dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(args[3], sprintf('chain-%02d-prefix', chain))
  write.csv(result, paste0(stem, '.csv'), row.names = FALSE)
  script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
  saveRDS(list(result = result, input_hashes = tools::md5sum(c(paths, statuses,
    file.path(roots, 'settings.rds'))), script_hash = tools::md5sum(script),
    checked_at = Sys.time(), session = sessionInfo(),
    scope = 'Exact overlapping trajectories; not convergence acceptance.'), paste0(stem, '.rds'))
  print(result)
}
