args <- commandArgs(TRUE)
 stopifnot(length(args) == 1L)
stopifnot(capabilities('png'))
options(bitmapType = if (capabilities('aqua')) 'quartz' else 'cairo')
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
materials <- file.path(dirname(dirname(script)), 'paper')
source(file.path(materials, 'appendix-chain-diagnostics.R'))
source(file.path(materials, 'plot-appendix-chains.R'))
root <- normalizePath(args[1])
fits <- lapply(file.path(root, sprintf('chain-%02d-diagnostic.rds', 1:8)), readRDS)
settings <- readRDS(file.path(root, 'settings.rds'))
for (f in fits) {
  stopifnot(length(f$posterior$log_posterior) == settings$draws+settings$burnin)
  for (m in f$posterior$interaction_draws) stopifnot(nrow(m) == settings$draws, all(is.finite(m)))
  for (m in f$posterior$inclusion_probabilities) stopifnot(all(is.finite(m)), all(m >= 0 & m <= 1))
}
appendix_chain_diagnostics(fits, file.path(root, 'diagnostics'))
withCallingHandlers(plot_appendix_chains(root), warning = function(w)stop(conditionMessage(w)))
plots <- list.files(file.path(root, 'figures'), pattern = '[.]png$', full.names = TRUE)
stopifnot(length(plots) == 12L, all(file.info(plots)$size > 1000))
warnings <- lapply(file.path(root, sprintf('chain-%02d.rds.status.rds', 1:8)),
  function(p)readRDS(p)$warnings)
saveRDS(warnings, file.path(root, 'diagnostics/chain-warnings.rds'))
writeLines(c('All eight chains collected. Diagnostics do not certify convergence.',
  'Review rank Rhat, ESS, MCSE, traces, half drift, inclusion rankings and warnings.',
  'New Linux runs are not an exact continuation of the archived macOS chains.'),
  file.path(root, 'INTERPRETATION.txt'))
