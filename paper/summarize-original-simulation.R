# Summarize independent simulation replicates, never individual CV folds.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1])
 s <- readRDS(file.path(root, 'settings.rds'))
stopifnot(s$experiment %in% c('simulation', 'correlated'))
for (name in c('simulation-acceptance.rds', 'baseline-selection-acceptance.rds')) {
  audit <- readRDS(file.path(root, name))
  stopifnot(identical(audit$settings, s))
}
stopifnot(identical(tools::md5sum(names(s$source_hashes)), s$source_hashes),
  identical(tools::md5sum(names(s$data_hash)), s$data_hash))
configs <- if (s$experiment == 'simulation') data.frame(scenario = 1:3, rho = 0, half = FALSE) else
  expand.grid(scenario = 1, rho = c(.2, .5, .8), half = c(FALSE, TRUE))
rows <- list()
 hashes <- character()
read_csv <- function(path) {
  hashes <<- c(hashes, tools::md5sum(path))
  read.csv(path, colClasses = c(subgroup = 'character'))
}
for (configuration in seq_len(nrow(configs))) for (replicate in seq_len(s$replicates)) {
  job <- file.path(root, sprintf('configuration-%02d-replicate-%03d', configuration, replicate))
  for (method in c('imr-molecular', 'bms-molecular', 'l1-cph', 'uni-cph')) {
    selection <- read_csv(file.path(job, paste0(method, '-selection.csv')))
    if (method != 'uni-cph') {
      result_path <- file.path(job, paste0(method, '.rds'))
      hashes <- c(hashes, tools::md5sum(result_path))
      result <- readRDS(result_path)
      if (inherits(result, 'imr_experiment_result_v1')) result <- result$result
      result$fit <- NULL
      stopifnot(isTRUE(all.equal(selection, result$selection,
        tolerance = 1e-12, check.attributes = FALSE)))
    }
    stopifnot(all(is.finite(selection$auc)), all(selection$auc >= 0 & selection$auc <= 1))
    rows[[length(rows) + 1L]] <- data.frame(configuration, scenario = configs$scenario[configuration],
      rho = configs$rho[configuration], half = configs$half[configuration], replicate,
      method, metric = 'selection_auc', platform = selection$platform,
      subgroup = selection$subgroup, value = selection$auc,
      folds_valid = NA_integer_, folds_expected = NA_integer_)
    if (method == 'uni-cph') next
    prediction <- read_csv(file.path(job, paste0(method, '-scores.csv')))
    stopifnot(isTRUE(all.equal(prediction, result$summary,
      tolerance = 1e-12, check.attributes = FALSE)), all(is.finite(prediction$mean)),
      all(prediction$n_valid == prediction$n_expected),
      all(prediction$n_expected == s$k * s$rounds))
    rows[[length(rows) + 1L]] <- data.frame(configuration, scenario = configs$scenario[configuration],
      rho = configs$rho[configuration], half = configs$half[configuration], replicate,
      method, metric = 'prediction_score', platform = 'combined',
      subgroup = prediction$subgroup, value = prediction$mean,
      folds_valid = prediction$n_valid, folds_expected = prediction$n_expected)
  }
}
values <- do.call(rbind, rows)
 rownames(values) <- NULL
keys <- c('configuration', 'method', 'metric', 'platform', 'subgroup')
key <- do.call(paste, c(values[keys], sep = ':'))
stopifnot(!anyDuplicated(paste(key, values$replicate, sep = ':')))
summary <- do.call(rbind, lapply(split(values, key), function(x) {
  stopifnot(nrow(x) == s$replicates, setequal(x$replicate, seq_len(s$replicates)))
  n <- nrow(x)
   deviation <- if (n > 1L) stats::sd(x$value) else NA_real_
  data.frame(x[1, c('configuration', 'scenario', 'rho', 'half', 'method', 'metric', 'platform', 'subgroup')],
    mean = mean(x$value), sd = deviation, se = deviation / sqrt(n),
    n_valid = n, n_requested = s$replicates,
    n_original_design = if (s$experiment == 'simulation') 50L else 30L)
}))
rownames(summary) <- NULL
write.csv(values, file.path(root, 'replicate-metrics.csv'), row.names = FALSE)
write.csv(summary, file.path(root, 'across-replicate-summary.csv'), row.names = FALSE)
saveRDS(list(settings = s, source_files = hashes, values = values, summary = summary,
  created_at = Sys.time(),
  uncertainty_unit = 'Independent simulation replicate; predictions first averaged over the recorded validation folds.',
  score_caution = 'Bayesian predictions use the configured score method; L1 uses standard sample-size-weighted subgroup concordance for the all row. Do not label different definitions as interchangeable.',
  scope = 'Executed manifest only; one replicate has undefined between-replicate SD and SE. Acceptance does not certify convergence or warned Cox inference.'),
  file.path(root, 'across-replicate-summary.rds'))
cat('Summarized', nrow(values), 'replicate-level values into', nrow(summary), 'groups; replicates per group:', s$replicates, '\n')
