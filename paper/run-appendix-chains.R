# Eight sequential, resumable original-scale chains for Appendix F diagnostics.
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1L]))
materials <- dirname(script)
source(file.path(materials, 'original-experiment-helpers.R'))
source(file.path(materials, 'appendix-chain-diagnostics.R'))
args <- parse_experiment_args(commandArgs(TRUE))
allowed <- c('--out-dir', '--data', '--draws', '--burnin', '--seed', '--quick', '--chain')
if (any(!names(args) %in% allowed)) stop('Unsupported chain-runner argument')
value <- function(flag, fallback) if (is.null(args[[flag]])) fallback else args[[flag]]
int <- function(flag, fallback, minimum) {
  x <- suppressWarnings(as.numeric(value(flag, fallback)))
  if (length(x) != 1L || !is.finite(x) || x != floor(x) || x < minimum || x > .Machine$integer.max - 8L)
    stop('Invalid ', flag)
  as.integer(x)
}
library(IntegMultiReg)
stopifnot(requireNamespace('coda', quietly = TRUE), requireNamespace('digest', quietly = TRUE))
quick <- isTRUE(args[['--quick']])
data_path <- normalizePath(value('--data', file.path(materials, 'data/kirc_table1_full.rda')))
out <- value('--out-dir', file.path(dirname(materials), 'output/appendix-chains', if (quick) 'quick' else 'full'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
settings <- list(draws = int('--draws', if (quick) 100L else 350000L, 4L),
  burnin = int('--burnin', if (quick) 20L else 50000L, 0L),
  seeds = int('--seed', 100L, 0L) + 0:7, sampler_method = 'paper', nu = c(-4, -3, -4),
  quick = quick, data_hash = tools::md5sum(data_path),
  source_hashes = tools::md5sum(c(script,
    file.path(materials, c('original-experiment-helpers.R', 'appendix-chain-diagnostics.R')),
    system.file('R', 'IntegMultiReg.rdb', package = 'IntegMultiReg'),
    system.file('libs', paste0('IntegMultiReg', .Platform$dynlib.ext), package = 'IntegMultiReg'))))
manifest <- file.path(out, 'settings.rds')
if (file.exists(manifest) && !identical(readRDS(manifest), settings)) stop('Settings/source/data changed; use a new output directory')
saveRDS(settings, manifest)
writeLines(capture.output(dput(settings)), file.path(out, 'settings.txt'))
capture.output(sessionInfo(), file = file.path(out, 'sessionInfo.txt'))
e <- new.env()
 load(data_path, e)
 d <- e$kirc_full
# A one-draw fit obtains public, correctly named initial-state templates only.
# Its result is never included in diagnostics.
template_fit <- imr(d$platforms, d$outcome.raw, covariates = d$covariates,
  outcome_type = 'right.censored', sampler_method = 'paper', nu = settings$nu, draws = 1, burnin = 0, seed = 99)
selection <- lapply(coef(template_fit), function(m) { storage.mode(m) <- "integer"
 m })
rm(template_fit)
 gc()
atomic_save <- function(x, path) {
  tmp <- paste0(path, '.tmp-', Sys.getpid())
   on.exit(unlink(tmp))
  saveRDS(x, tmp)
   if (!file.rename(tmp, path)) stop('Could not commit ', path)
}
chains <- if (is.null(args[['--chain']])) 1:8 else int('--chain', 1L, 1L)
stopifnot(all(chains <= 8L))
record_task_selection(out, list(kind = 'chains', chains = chains))
slim <- vector('list', 8L)
for (chain in chains) {
  path <- file.path(out, sprintf('chain-%02d.rds', chain))
  initial <- list(selection = selection, interaction = lapply(selection, function(m)
    matrix(0, nrow(m), nrow(m), dimnames = list(rownames(m), rownames(m)))))
  for (p in seq_along(selection)) {
    initial$selection[[p]][] <- 0L
    if (chain > 1L) initial$selection[[p]][, seq_len(min(chain - 1L, ncol(selection[[p]])))] <- 1L
    initial$interaction[[p]][] <- .025 * chain
    diag(initial$interaction[[p]]) <- 0
  }
  initial_path <- file.path(out, sprintf('chain-%02d-initial.rds', chain))
  if (file.exists(initial_path)) stopifnot(identical(readRDS(initial_path), initial))
  else atomic_save(initial, initial_path)
  diagnostic_path <- file.path(out, sprintf('chain-%02d-diagnostic.rds', chain))
  if (experiment_job_complete(path) && file.exists(diagnostic_path)) {
    cat(format(Sys.time()), 'chain', chain, 'resumed completed checkpoint\n')
    slim[[chain]] <- readRDS(diagnostic_path)
    next
  }
  cat(format(Sys.time()), 'chain', chain, 'starting; seed', settings$seeds[chain], '\n')
  record_experiment_job(path, {
    fit <- experiment_fit_checkpoint(path,
      imr(d$platforms, d$outcome.raw, covariates = d$covariates,
        outcome_type = 'right.censored', nu = settings$nu, sampler_method = 'paper',
        draws = settings$draws, burnin = settings$burnin, seed = settings$seeds[chain], initial = initial))
    stopifnot(identical(fit$control$initial, initial))
    # Retain full fit losslessly on disk; read only diagnostic fields across chains.
    summary <- list(model = fit$model, control = fit$control, posterior = fit$posterior[
      c('inclusion_probabilities', 'interaction_draws', 'log_posterior')])
    atomic_save(summary, diagnostic_path)
    slim[[chain]] <- summary
    rm(fit)
     gc()
    TRUE
  })
  cat(format(Sys.time()), 'chain', chain, 'completed\n')
}
if (length(chains) == 8L) appendix_chain_diagnostics(slim, file.path(out, 'diagnostics'))
writeLines(c('Eight distinct, explicitly recorded starts and distinct seeds.',
  'Starts are new documented choices, not recovered historical values.',
  'Full fits use lossless compact checkpoints; diagnostic files omit selection histories.',
  'No automatic convergence decision: inspect diagnostics and full log-posterior traces.',
  if (quick) 'SHORT INTERFACE CHECK ONLY; not scientific convergence evidence.' else
    'Original iteration budget; assess mixing and sensitivity before interpretation.'), file.path(out, 'INTERPRETATION.txt'))
cat(format(Sys.time()), 'Requested chains completed; collect all eight before interpretation\n')
