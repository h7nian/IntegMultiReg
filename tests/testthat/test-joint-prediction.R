test_that("point and interval predictions agree with independent stored-draw calculations", {
  f <- small_fit()
  d <- small_data()
  standardized <- scale(d$platforms$assay$marker,
    center = f$preprocessing$feature_center[[1]][[1]], scale = f$preprocessing$feature_scale[[1]][[1]]
  )
  X <- cbind(1, standardized)
  b <- matrix(f$posterior$coefficients[[1]], nrow = 80, ncol = 2)
  eta <- b %*% t(X)
  p <- predict(f, d$platforms, interval = TRUE)[[1]]
  expect_named(p, c("id", "prediction", "lower", "upper"))
  expect_equal(p$prediction, colMeans(eta))
  expect_equal(p$lower, apply(eta, 2, quantile, .025, type = 1, names = FALSE))
  expect_equal(p$upper, apply(eta, 2, quantile, .975, type = 1, names = FALSE))
  expect_equal(predict(f, d$platforms, type = "link")[[1]]$prediction, p$prediction)
  set.seed(112)
  rng <- .Random.seed
  a <- predict(f, d$platforms, quantity = "new_observation", interval = TRUE, seed = 33)
  expect_identical(.Random.seed, rng)
  expect_identical(a, predict(f, d$platforms, quantity = "new_observation", interval = TRUE, seed = 33))
  expect_gt(mean(a[[1]]$upper - a[[1]]$lower), mean(p$upper - p$lower))
})

test_that("binary prediction averages probabilities before applying observation noise", {
  f <- small_fit("binary")
  d <- small_data("binary")
  X <- cbind(1, as.numeric(scale(d$platforms$assay$marker)))
  b <- matrix(f$posterior$coefficients[[1]], 80, 2)
  eta <- b %*% t(X)
  v <- as.vector(f$posterior$variance)
  expected <- colMeans(pnorm(eta / sqrt(v)))
  p <- predict(f, d$platforms, interval = TRUE)[[1]]
  expect_equal(p$prediction, expected)
  y <- predict(f, d$platforms, quantity = "new_observation", interval = TRUE)[[1]]
  expect_true(all(c(y$lower, y$upper) %in% 0:1))
  expect_equal(predict(f, d$platforms, type = "link")[[1]]$prediction, colMeans(eta))
})

test_that("survival medians transform draw-wise conditional means", {
  f <- small_fit("right.censored")
  d <- small_data("right.censored")
  b <- matrix(f$posterior$coefficients[[1]], 80, 2)
  X <- cbind(1, as.numeric(scale(d$platforms$assay$marker)))
  eta <- b %*% t(X)
  v <- as.vector(f$posterior$variance)
  expected <- exp(apply(eta + v / 2, 2, quantile, .5, type = 1, names = FALSE))
  p <- predict(f, d$platforms, interval = TRUE)[[1]]
  expect_equal(p$prediction, expected)
  expect_equal(predict(f, d$platforms, type = "link")[[1]]$prediction, colMeans(eta))
  expect_true(all(p$lower > 0 & p$upper >= p$lower))
})

test_that("available-platform routing uses names and stored training transforms", {
  f <- fit_bin
  one <- simIMR$platforms[c("metabolomic")]
  named <- predict(f, one, covariates = simIMR$covariates)
  indexed <- predict(f, unname(one), platform_names = 3, covariates = simIMR$covariates)
  expect_identical(named, indexed)
  expect_s3_class(named, "imr_predictions")
  expect_output(print(named), "metabolomic")
  expect_error(predict(f, one, platform_names = "unknown", covariates = simIMR$covariates), "platform")
  expect_error(predict(f, one, max_models = 3), "only to a Laplace")
})
