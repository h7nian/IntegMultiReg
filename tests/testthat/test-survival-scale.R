test_that("survival time transformation agrees with the original ReadData.c", {
  x <- data.frame(id = 1:14, x = seq(-1, 1, length.out = 14))
  time <- exp(seq(-2, 2, length.out = 14))
  y <- data.frame(id = x$id, time = time, status = 1)
  survival <- imr(list(assay = x), y, outcome_type = "right.censored",
                  min_subgroup_size = 2, draws = 12, burnin = 6, seed = 23)
  continuous <- imr(list(assay = x), data.frame(id = x$id, y = log(time)),
                    outcome_type = "continuous", min_subgroup_size = 2,
                    draws = 12, burnin = 6, seed = 23)
  expect_equal(as.numeric(survival$preprocessing$response[[1]][, 1]), log(time))
  expect_equal(survival$posterior$latent_response_mean, continuous$posterior$latent_response_mean)
  expect_equal(survival$posterior$selection_draws, continuous$posterior$selection_draws)
  expect_identical(survival$control$response_scale, "log")
  expect_equal(predict(survival, list(assay = x)), predict(continuous, list(assay = x)))
  old <- imr(list(assay = x), y, outcome_type = "right.censored",
             survival_scale = "identity", min_subgroup_size = 2,
             draws = 12, burnin = 6, seed = 23)
  expect_equal(as.numeric(old$preprocessing$response[[1]][, 1]), time)
  expect_identical(old$control$response_scale, "identity")
  expect_error(imr(list(assay = x), y, survival_scale = "invalid"), "arg")
})

test_that("censored times below one use finite negative log bounds", {
  x <- data.frame(id = 1:14, x = seq(-1, 1, length.out = 14))
  y <- data.frame(id = x$id, time = exp(seq(-3, 2, length.out = 14)),
                  status = rep(c(0, 1), 7))
  f <- imr(list(assay = x), y, min_subgroup_size = 2, draws = 20, burnin = 10, seed = 23)
  expect_true(all(is.finite(f$posterior$log_posterior)))
  expect_true(all(f$posterior$latent_response_mean[[1]][y$status == 0] >= log(y$time[y$status == 0])))
})
