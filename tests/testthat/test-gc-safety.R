test_that("native return objects survive frequent garbage collection", {
  x <- list(assay = data.frame(id = seq_len(12), x = seq_len(12) / 12))
  y <- data.frame(id = seq_len(12), y = sin(seq_len(12)))
  exercise <- function() {
    previous <- gctorture2(10L)
    on.exit(gctorture2(previous), add = TRUE)
    fit <- imr(x, y, type_outcome = "continuous", ssize = 1,
               sample_mcmc = c(3, 1), seed = 11)
    predictions <- predict(fit, x, max_models = 2)
    cv <- cv_imr(fit, k = 2, rounds = 1, max_models = 2)
    list(fit = fit, predictions = predictions, cv = cv)
  }
  result <- exercise()
  expect_true(validate_imr(result$fit))
  expect_true(all(is.finite(result$predictions[[1L]]$predict)))
  expect_true(all(is.finite(result$cv$total_cindex)))
})
