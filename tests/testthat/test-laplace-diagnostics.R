test_that("fitting records Laplace limits by subgroup and stage", {
  id <- 1:40
  fit <- imr(list(a = data.frame(id, a = sin(id), b = cos(id))),
    data.frame(id, y = sin(id / 4)), outcome_type = "continuous", method = "bms",
    draws = 10, burnin = 5, laplace_max_iter = 1L,
    laplace_tolerance = 1e-12, seed = 31)
  d <- fit$control$laplace_diagnostics
  expect_identical(d$stage, c("initial", "selection", "latent"))
  expect_identical(d$subgroup, rep("1", 3))
  expect_equal(d$calls, c(1, 15, 0))
  expect_gt(sum(d$iteration_limit), 0)
  expect_true(all(d$iteration_limit <= d$calls))
  expect_true(all(d$nonfinite == 0))
  expect_true(all(d$factorization_failures == 0))
  expect_true(validate_imr(fit))
})
