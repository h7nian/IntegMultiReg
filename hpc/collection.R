# Shared collector for verified task directories in manifest order.
collect_tasks <- function(kind, tasks, out) {
settings <- readRDS(file.path(tasks[1], 'settings.rds'))
stopifnot(identical(tools::md5sum(names(settings$source_hashes)), settings$source_hashes),
  identical(tools::md5sum(names(settings$data_hash)), settings$data_hash))
if (kind != 'chains') stopifnot(identical(settings$experiment, kind))
n <- if (kind == 'chains') 8L else settings$replicates * if (kind == 'simulation') 3L else 6L
stopifnot(length(tasks) == n)
from <- to <- character()
add <- function(src, dst) {from <<- c(from, src)
  to <<- c(to, dst)}
for (id in seq_len(n)) {
    task <- tasks[id]
    stopifnot(identical(readRDS(file.path(task, 'settings.rds')), settings))
    selection <- readRDS(file.path(task, 'task.rds'))
    if (kind == 'chains') {
    stopifnot(identical(selection, list(kind = 'chains', chains = as.integer(id))))
    paths <- file.path(task, sprintf(c('chain-%02d.rds', 'chain-%02d-diagnostic.rds',
      'chain-%02d-initial.rds', 'chain-%02d.rds.status.rds'), id))
    stopifnot(all(file.exists(paths)), identical(readRDS(paths[4])$status, 'completed'))
    d <- readRDS(paths[2])
     s <- readRDS(paths[1])
     f <- s$fit
    stopifnot(inherits(s, 'imr_experiment_fit_checkpoint_v1'), length(s$layout) == settings$draws,
      f$control$seed == settings$seeds[id], f$control$mcmc$draws == settings$draws,
      f$control$mcmc$burnin == settings$burnin, identical(f$control$initial, readRDS(paths[3])),
      identical(f$control$sampler_method, settings$sampler_method),
      identical(d, list(model = f$model, control = f$control, posterior = f$posterior[
        c('inclusion_probabilities', 'interaction_draws', 'log_posterior')])),
      all(is.finite(d$posterior$log_posterior)))
    for (p in seq_along(s$pools)) {
      index <- s$indices[[p]]
      stopifnot(length(index) == settings$draws, !anyNA(index),
        all(index >= 1L & index <= length(s$pools[[p]])),
        all(vapply(s$pools[[p]], function(m)all(m %in% 0:1), TRUE)))
    }
    rm(s, f, d)
     gc()
    files <- list.files(task, pattern = sprintf('^chain-%02d', id), full.names = TRUE)
    add(files, basename(files))
    } else {
    config <- as.integer((id - 1L) %/% settings$replicates + 1L)
    rep <- as.integer((id - 1L) %% settings$replicates + 1L)
    stopifnot(identical(selection, list(kind = kind, configurations = config, replicates = rep)))
    job <- sprintf('configuration-%02d-replicate-%03d', config, rep)
    required <- file.path(task, job, c('imr-molecular.rds', 'bms-molecular.rds',
      'l1-cph.rds', 'uni-cph-selection.csv'))
    for (path in required) stopifnot(file.exists(path),
      identical(readRDS(paste0(path, '.status.rds'))$status, 'completed'))
    for (file in list.files(file.path(task, job), recursive = TRUE))
      add(file.path(task, job, file), file.path(job, file))
    }
    add(file.path(task, 'task.rds'), file.path('task-provenance', sprintf('%03d-task.rds', id)))
    add(file.path(task, 'sessionInfo.txt'), file.path('task-provenance', sprintf('%03d-sessionInfo.txt', id)))
}
stopifnot(!anyDuplicated(to))
for (file in c('settings.rds', 'settings.txt', 'sessionInfo.txt')) add(file.path(tasks[1], file), file)
hashes <- tools::md5sum(from)
if (dir.exists(out)) {
    saved <- readRDS(file.path(out, 'COLLECTION.rds'))
    stopifnot(identical(saved$input_hashes, hashes), identical(saved$destinations, to),
      identical(unname(tools::md5sum(file.path(out, to))), unname(hashes)))
    cat('Existing complete collection verified; continuing validation.\n')
} else {
    temporary <- paste0(out, '.partial-', Sys.getpid())
    dir.create(temporary, recursive = TRUE)
    for (i in seq_along(from)) {
    target <- file.path(temporary, to[i])
     dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(from[i], target), identical(unname(tools::md5sum(target)), unname(hashes[i])))
    }
    saveRDS(list(input_hashes = hashes, destinations = to, kind = kind, tasks = n,
      collected_at = Sys.time(), note = 'File collection only; scientific validation follows.'),
      file.path(temporary, 'COLLECTION.rds'))
    stopifnot(file.rename(temporary, out))
    cat('Collected', n, 'complete tasks without changing their scientific settings.\n')
}

}
