test_that("post-fit task stages preserve full native results and partition order", {
  call_native <- IntegMultiReg:::.imr_call_cv_postfit_native
  outcomes <- list(binary = simIMR$outcome.binary,
                   continuous = simIMR$outcome.continuous,
                   right.censored = simIMR$outcome.survival)
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      fit <- imr(simIMR$platforms, outcomes[[outcome_type]],
                 covariates = simIMR$covariates, outcome_type = outcome_type,
                 method = method, draws = 12, burnin = 6,
                 min_subgroup_size = 30, seed = 53)
      for (importance in c(FALSE, TRUE)) {
        context <- paste(outcome_type, method, importance)
        run <- function(...) call_native(fit, 3L, 2L, 4L, FALSE, importance, ...)
        reference <- run()
        plan <- run(stage = "plan")
        expect_identical(plan$folds, reference$folds, info = context)
        expect_true(all(is.na(plan$predictions)), info = context)
        expect_true(all(is.na(plan$total_cindex)), info = context)
        expect_true(all(is.na(plan$subset_cindex)), info = context)
        merged <- plan$predictions
        # Nonconsecutive batches deliberately include the second round first.
        # Survival's censored shuffle is carried across rounds, not reset.
        for (indices in list(c(6L, 2L, 4L), c(1L, 3L, 5L))) {
          tasks <- matrix(FALSE, 3L, 2L)
          tasks[indices] <- TRUE
          partial <- run(stage = "predict", tasks = tasks)
          expect_identical(partial$folds, reference$folds, info = context)
          expected <- reference$predictions
          for (round in 1:2) {
            expected[!tasks[reference$folds[, round], round], round] <- NA_real_
          }
          expect_identical(partial$predictions, expected, info = context)
          expect_true(all(is.na(partial$total_cindex)), info = context)
          use <- !is.na(partial$predictions)
          merged[use] <- partial$predictions[use]
        }
        expect_identical(merged, reference$predictions, info = context)
        expect_identical(run(stage = "score", predictions = merged), reference,
                         info = context)
        empty <- run(stage = "predict", tasks = matrix(FALSE, 3L, 2L))
        expect_true(all(is.na(empty$predictions)), info = context)
      }
    }
  }
})

test_that("internal post-fit stage controls reject malformed tasks and scores", {
  run <- function(...) IntegMultiReg:::.imr_call_cv_postfit_native(
    fit_bin, 2L, 1L, 4L, FALSE, FALSE, ...)
  expect_error(run(stage = "unknown"), "arg")
  for (tasks in list(c(TRUE, TRUE), matrix(1, 2, 1), matrix(NA, 2, 1),
                     matrix(TRUE, 1, 2))) {
    expect_error(run(stage = "predict", tasks = tasks), "task mask")
  }
  expect_error(run(stage = "plan", tasks = matrix(FALSE, 2, 1)), "task mask")
  expect_error(run(stage = "score"), "scoring predictions")
  expect_error(run(stage = "score", predictions = matrix(Inf, sum(fit_bin$model$sample_sizes), 1)),
                "scoring predictions")
  expect_error(run(stage = "plan", predictions = matrix(0, 1, 1)), "only.*scoring")
})
