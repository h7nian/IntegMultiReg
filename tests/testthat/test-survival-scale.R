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
  expect_error(compare_fit_summaries(survival, old), "same outcome type, response scale")
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

test_that("censored working times respect the historical proposal support", {
  x <- data.frame(id = 1:14, x = seq(-1, 1, length.out = 14))
  y <- data.frame(id = x$id, time = seq(1, 3, length.out = 14), status = 1)
  y$status[1] <- 0
  run <- function(scale) imr(list(assay = x), y,
    outcome_type = "right.censored", survival_scale = scale,
    min_subgroup_size = 2, draws = 4, burnin = 2, seed = 23)
  for (time in c(1000.49, 1000.5, 1200)) {
    y$time[1] <- time
    expect_error(run("identity"), "0.01 initialization below 1000.5", fixed = TRUE)
  }
  # The same raw time is supported on the default log scale.
  log_fit <- run("log")
  expect_true(all(is.finite(log_fit$posterior$log_posterior)))
  expect_gt(log_fit$posterior$latent_response_mean[[1]][1], log(y$time[1]))
  # The cap applies to augmented censored values, not observed event times.
  y$status[1] <- 1
  observed_fit <- run("identity")
  expect_equal(observed_fit$posterior$latent_response_mean[[1]][1], y$time[1])
})
