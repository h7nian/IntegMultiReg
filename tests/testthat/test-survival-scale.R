test_that("survival time transformation agrees with the original ReadData.c", {
  x <- data.frame(id = 1:14, x = seq(-1, 1, length.out = 14))
  time <- exp(seq(-2, 2, length.out = 14))
  y <- data.frame(id = x$id, time = time, status = 1)
  survival <- imr(list(assay = x), y, type_outcome = "right.censored",
                  ssize = 2, sample_mcmc = c(12, 6), seed = 23)
  continuous <- imr(list(assay = x), data.frame(id = x$id, y = log(time)),
                    type_outcome = "continuous", ssize = 2,
                    sample_mcmc = c(12, 6), seed = 23)
  expect_equal(as.numeric(survival$data2$yy[[1]][, 1]), log(time))
  expect_equal(survival$estimate_latent_y, continuous$estimate_latent_y)
  expect_equal(survival$gam_sample, continuous$gam_sample)
  expect_identical(survival$response_scale, "log")
  expect_equal(predict(survival, list(assay = x)), predict(continuous, list(assay = x)))
  old <- imr(list(assay = x), y, type_outcome = "right.censored",
             survival_scale = "identity", ssize = 2,
             sample_mcmc = c(12, 6), seed = 23)
  expect_equal(as.numeric(old$data2$yy[[1]][, 1]), time)
  expect_identical(old$response_scale, "identity")
  expect_error(imr(list(assay = x), y, survival_scale = "invalid"), "arg")
})

test_that("censored times below one use finite negative log bounds", {
  x <- data.frame(id = 1:14, x = seq(-1, 1, length.out = 14))
  y <- data.frame(id = x$id, time = exp(seq(-3, 2, length.out = 14)),
                  status = rep(c(0, 1), 7))
  f <- imr(list(assay = x), y, ssize = 2, sample_mcmc = c(20, 10), seed = 23)
  expect_true(all(is.finite(f$log_posterior)))
  expect_true(all(f$estimate_latent_y[[1]][y$status == 0] >= log(y$time[y$status == 0])))
})
