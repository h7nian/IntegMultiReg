test_that("refit task dispatch preserves complete results for each outcome and model", {
  outcomes <- list(binary = simIMR$outcome.binary,
                   continuous = simIMR$outcome.continuous,
                   right.censored = simIMR$outcome.survival)
  evaluate <- IntegMultiReg:::.imr_cv_refit_result
  # R CMD check can restrict tests to two processes. Standalone verification
  # also exercises three workers; do not override the checker's core budget.
  limited <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  worker_counts <- if (nzchar(limited) && limited != "false") 2L else 2:3
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      fit <- imr(simIMR$platforms, outcomes[[outcome_type]],
                 covariates = simIMR$covariates, outcome_type = outcome_type,
                 method = method, draws = 12, burnin = 6,
                 min_subgroup_size = 30, seed = 53)
      reference <- cv_imr(fit, k = 3, rounds = 2, max_models = 4,
                          cv_method = "refit")
      for (workers in worker_counts) {
        set.seed(871)
        rng <- .Random.seed
        result <- evaluate(fit, 3L, 2L, 4L, FALSE, workers)
        expect_identical(result, reference, info = paste(outcome_type, method, workers))
        expect_identical(.Random.seed, rng)
      }
    }
  }
})

test_that("refit tasks preserve formula transformations and serialization", {
  previous <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(previous))
  clinical <- merge(simIMR$outcome.continuous, simIMR$covariates,
                    by = "id", sort = FALSE)
  clinical$group <- factor(ifelse(clinical$age > 0, "older", "younger"))
  fit <- imr(y ~ poly(age, 2) + group * sex, data = clinical,
             platforms = simIMR$platforms, outcome_type = "continuous",
             draws = 12, burnin = 6, min_subgroup_size = 30, seed = 73)
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(IntegMultiReg:::.imr_cv_refit_result(restored, 2L, 1L, 4L, FALSE, 2L),
                   cv_imr(fit, k = 2, rounds = 1, max_models = 4, cv_method = "refit"))
})

test_that("worker dispatch is bounded, ordered and restores process state", {
  map <- IntegMultiReg:::.imr_cv_map
  rng <- IntegMultiReg:::.imr_save_rng()
  on.exit(IntegMultiReg:::.imr_restore_rng(rng))
  IntegMultiReg:::.imr_restore_rng(NULL)
  settings <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                           "BLIS_NUM_THREADS"), unset = NA_character_)
  connections <- showConnections(all = TRUE)
  result <- map(as.list(1:3), function(x) {
    list(index = x, threads = Sys.getenv("OMP_NUM_THREADS"))
  }, 2L)
  expect_identical(vapply(result, `[[`, integer(1), "index"), 1:3)
  expect_identical(vapply(result, `[[`, character(1), "threads"), rep("1", 3))
  expect_identical(map(list(1L), identity, 20L), list(1L))
  expect_identical(map(list(), identity, 2L), list())
  expect_identical(Sys.getenv(names(settings), unset = NA_character_), settings)
  expect_identical(showConnections(all = TRUE), connections)
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  for (invalid in list(0, -1, 1.5, NA, Inf, TRUE, "2", c(1, 2))) {
    expect_error(map(list(1L), identity, invalid), "workers")
  }
})

test_that("worker errors identify the fold and close sockets without partial results", {
  map <- IntegMultiReg:::.imr_cv_map
  task <- list(round = 2L, fold = 3L, seed = 1L,
               train_ids = "missing", test_ids = "missing")
  connections <- showConnections(all = TRUE)
  settings <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  set.seed(15)
  rng <- .Random.seed
  expect_error(map(list(task, task), IntegMultiReg:::.imr_cv_refit_task,
                   2L, object = fit_bin, max_models = 4L, verbose = FALSE),
                "CV round 2, fold 3 failed")
  expect_identical(showConnections(all = TRUE), connections)
  expect_identical(Sys.getenv("OMP_NUM_THREADS", unset = NA_character_), settings)
  expect_identical(.Random.seed, rng)
  expect_identical(map(list(1L, 2L), identity, 2L), list(1L, 2L))
})

test_that("workers inherit the caller's RNG kind without consuming its state", {
  previous_kind <- RNGkind()
  previous_rng <- IntegMultiReg:::.imr_save_rng()
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    IntegMultiReg:::.imr_restore_rng(previous_rng)
  })
  RNGkind("L'Ecuyer-CMRG")
  set.seed(183)
  rng <- .Random.seed
  kinds <- IntegMultiReg:::.imr_cv_map(list(1L, 2L), function(task) RNGkind(), 2L)
  expect_identical(kinds, rep(list(RNGkind()), 2L))
  expect_identical(.Random.seed, rng)
})
