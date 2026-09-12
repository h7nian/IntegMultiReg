test_that("parallel post-fit batches preserve complete results and caller RNG", {
  outcomes <- list(binary = simIMR$outcome.binary,
                   continuous = simIMR$outcome.continuous,
                   right.censored = simIMR$outcome.survival)
  limited <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  worker_counts <- if (nzchar(limited) && limited != "false") 2L else 2:3
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      fit <- imr(simIMR$platforms, outcomes[[outcome_type]],
                 covariates = simIMR$covariates, outcome_type = outcome_type,
                 method = method, draws = 12, burnin = 6,
                 min_subgroup_size = 30, seed = 53)
      for (mode in c("legacy", "importance")) {
        reference <- cv_imr(fit, k = 3L, rounds = 2L, max_models = 4L, cv_method = mode)
        for (workers in worker_counts) {
          set.seed(931)
          rng <- .Random.seed
          result <- IntegMultiReg:::.imr_cv_postfit(fit, 3L, 2L, 4L, FALSE, mode, workers)
          expect_identical(result, reference, info = paste(outcome_type, method, mode, workers))
          expect_identical(.Random.seed, rng)
        }
      }
    }
  }
})

test_that("post-fit worker failures report assigned folds and close connections", {
  call_native <- IntegMultiReg:::.imr_call_cv_postfit_native
  plan <- call_native(fit_bin, 2L, 1L, 4L, FALSE, FALSE, stage = "plan")
  damaged <- plan$folds
  damaged[1L] <- 3L - damaged[1L]
  connections <- showConnections(all = TRUE)
  expect_error(IntegMultiReg:::.imr_cv_map(list(1L, 2L),
    IntegMultiReg:::.imr_cv_postfit_task, 2L, object = fit_bin,
    k = 2L, rounds = 1L, max_models = 4L, verbose = FALSE, importance = FALSE,
    folds = damaged), "round 1/fold [12].*planned partitions")
  expect_identical(showConnections(all = TRUE), connections)
})
