script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), 'original-experiment-helpers.R'))
library(IntegMultiReg)
id <- 1:30
fit <- imr(list(a = data.frame(id, x = sin(id)), b = data.frame(id, x = cos(id))),
 data.frame(id, y = sin(id/3)), outcome_type = 'continuous', draws = 30, burnin = 5,
 min_subgroup_size = 0, seed = 71)
path <- tempfile(fileext = '.rds')
set.seed(42)
 rng <- .Random.seed
saved <- experiment_fit_checkpoint(path, fit)
stopifnot(identical(.Random.seed, rng), identical(saved, fit))
container <- readRDS(path)
stopifnot(inherits(container, 'imr_experiment_fit_checkpoint_v1'))
restored <- experiment_fit_checkpoint(path, stop('must reuse checkpoint'))
stopifnot(identical(restored, fit), identical(.Random.seed, rng))
stopifnot(identical(cv_imr(fit, k = 2, rounds = 1), cv_imr(restored, k = 2, rounds = 1)))
states <- lapply(restored$posterior$selection_draws, `[[`, 1L)
j <- which(duplicated(states))[1L]
i <- which(vapply(states, identical, logical(1), states[[j]]))[1L]
if (capabilities('profmem')) {
  a <- tracemem(states[[i]])
  b <- tracemem(states[[j]])
  untracemem(states[[i]])
  untracemem(states[[j]])
  stopifnot(identical(a, b))
}
old <- states[[i]][1L]
states[[j]][1L] <- 1L-old
stopifnot(identical(states[[i]][1L], old))
saveRDS(fit, path) # legacy checkpoint remains readable
stopifnot(identical(experiment_fit_checkpoint(path, stop('must reuse')), fit))
unlink(path)
cat('PASS: lossless full fit, CV replay, RNG, shared reload, copy-on-modify, legacy checkpoint\n')
