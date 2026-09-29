test_that("native return objects survive frequent garbage collection", {
  x <- list(assay = data.frame(id = seq_len(12), x = seq_len(12) / 12))
  y <- data.frame(id = seq_len(12), y = sin(seq_len(12)))
  exercise <- function() {
    previous <- gctorture2(10L)
    on.exit(gctorture2(previous), add = TRUE)
    fit <- imr(x, y, outcome_type = "continuous", min_subgroup_size = 1,
               draws = 3, burnin = 1, seed = 11)
    predictions <- predict(fit, x, max_models = 2)
    cv <- cv_imr(cv_method = "refit", fit, k = 2, rounds = 1, max_models = 2)
    postfit <- lapply(c("legacy", "importance"), function(mode) {
      cv_imr(fit, k = 2, rounds = 1, max_models = 2, cv_method = mode)
    })
    list(fit = fit, predictions = predictions, cv = cv, postfit = postfit)
  }
  result <- exercise()
  expect_true(validate_imr(result$fit))
  expect_true(all(is.finite(result$predictions[[1L]]$prediction)))
  expect_true(all(is.finite(result$cv$pooled)))
  expect_true(all(vapply(result$postfit, function(x) all(is.finite(x$pooled)),
                         logical(1))))
})
