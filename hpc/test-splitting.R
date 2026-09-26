# Exact local smoke comparison of serial vs independently selected real tasks.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
out <- args[2]
stopifnot(!dir.exists(out))
dir.create(out, recursive = TRUE)
out <- normalizePath(out)
run_study_script <- function(script, flags, success = TRUE) {
  log <- tempfile('command-', tmpdir = out, fileext = '.log')
  code <- system2(file.path(R.home('bin'), 'Rscript'),
    c('--vanilla', shQuote(file.path(root, script)), vapply(flags, shQuote, '')), stdout = log, stderr = log)
  if (success && code != 0L) stop('Command failed: ', log)
  if (!success && code == 0L) stop('Negative control incorrectly succeeded: ', log)
  invisible(log)
}
chains <- file.path(out, 'chain-tasks')
serial <- file.path(out, 'chain-serial')
run_study_script('paper/run-appendix-chains.R', c('--quick', '--out-dir', serial))
for (i in 1:8) {
  task <- file.path(chains, sprintf('%03d', i))
  run_study_script('paper/run-appendix-chains.R', c('--quick', '--chain', i, '--out-dir', task))
  for (name in sprintf(c('chain-%02d.rds', 'chain-%02d-diagnostic.rds', 'chain-%02d-initial.rds'), i))
    stopifnot(identical(readRDS(file.path(serial, name)), readRDS(file.path(task, name))))
  cat('Exact chain comparison', i, 'PASS\n')
}
run_study_script('hpc/collect-legacy-layout.R', c('chains', chains, file.path(out, 'chains-combined')))
# Wrong task selection and changed settings must reject reuse.
run_study_script('paper/run-appendix-chains.R', c('--quick', '--chain', '2', '--out-dir', file.path(chains, '001')), FALSE)
run_study_script('paper/run-appendix-chains.R', c('--quick', '--chain', '1', '--seed', '101', '--out-dir', file.path(chains, '001')), FALSE)
status_path <- file.path(chains, '001/chain-01.rds.status.rds')
status <- readRDS(status_path)
changed <- status
changed$status <- 'running'
saveRDS(changed, status_path)
run_study_script('hpc/collect-legacy-layout.R', c('chains', chains, file.path(out, 'must-reject-chains')), FALSE)
saveRDS(status, status_path)
sim_serial <- file.path(out, 'simulation-serial')
sim_tasks <- file.path(out, 'simulation-tasks')
flags <- c('--quick', '--experiment', 'simulation', '--reference', 'code2017', '--replicates', '2')
run_study_script('paper/original-experiments.R', c(flags, '--out-dir', sim_serial))
id <- 0L
for (config in 1:3) for (rep in 1:2) {
  id <- id + 1L
  task <- file.path(sim_tasks, sprintf('%03d', id))
  run_study_script('paper/original-experiments.R', c(flags, '--configuration', config, '--replicate', rep, '--out-dir', task))
  job <- sprintf('configuration-%02d-replicate-%03d', config, rep)
  for (name in c('data.rds', 'imr-molecular.rds', 'bms-molecular.rds', 'l1-cph.rds')) {
    a <- readRDS(file.path(sim_serial, job, name))
    b <- readRDS(file.path(task, job, name))
    # Runtime measurements are not scientific outputs.
    if (inherits(a, 'imr_experiment_result_v1')) {a$result$elapsed <- NULL
    b$result$elapsed <- NULL}
    stopifnot(identical(a, b))
  }
  a <- readLines(file.path(sim_serial, job, 'uni-cph-selection.csv'))
  b <- readLines(file.path(task, job, 'uni-cph-selection.csv'))
  stopifnot(identical(a, b))
  cat('Exact simulation comparison', config, rep, 'PASS\n')
}
combined <- file.path(out, 'simulation-combined')
run_study_script('hpc/collect-legacy-layout.R', c('simulation', sim_tasks, combined))
run_study_script('paper/validate-original-simulation.R', combined)
report <- read.csv(file.path(combined, 'simulation-acceptance.csv'))
stopifnot(nrow(report) == 18L, setequal(report$config, 1:3), setequal(report$rep, 1:2))
# A late configuration must be checked, not silently skipped after an earlier fit.
late_status <- file.path(combined, 'configuration-03-replicate-002/bms-molecular.rds.status.rds')
original_status <- readBin(late_status, 'raw', n = file.info(late_status)$size)
bad <- readRDS(late_status)
bad$status <- 'failed'
saveRDS(bad, late_status)
run_study_script('paper/validate-original-simulation.R', combined, FALSE)
writeBin(original_status, late_status)
run_study_script('paper/validate-simulation-baselines.R', combined)
run_study_script('paper/summarize-original-simulation.R', combined)
run_study_script('paper/validate-original-simulation.R', c(combined, '--require-complete-study'), FALSE)
run_study_script('hpc/collect-legacy-layout.R', c('simulation', sim_tasks, combined))
stopifnot(file.rename(file.path(sim_tasks, '006'), file.path(sim_tasks, 'held-006')))
run_study_script('hpc/collect-legacy-layout.R', c('simulation', sim_tasks, file.path(out, 'must-reject-missing')), FALSE)
stopifnot(file.rename(file.path(sim_tasks, 'held-006'), file.path(sim_tasks, '006')))
run_study_script('paper/original-experiments.R', c(flags, '--configuration', '1', '--replicate', '3', '--out-dir', file.path(out, 'must-reject-replicate')), FALSE)
writeLines('PASS: real serial/split chain and simulation equality; collection; scientific validators; negative controls.', file.path(out, 'PASS'))
run_study_script('hpc/check-smoke.R', c(file.path(chains, '001'), file.path(sim_tasks, '003')))
