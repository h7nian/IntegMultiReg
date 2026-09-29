test_that("importance avoids historical scores while legacy retains them", {
  importance_log <- capture.output(importance <- cv_imr(
    fit_bin, k = 2L, rounds = 1L, cv_method = "importance", verbose = TRUE))
  legacy_log <- capture.output(legacy <- cv_imr(
    fit_bin, k = 2L, rounds = 1L, max_models = 2L,
    cv_method = "legacy", verbose = TRUE))
  expect_false(any(grepl("Average of", importance_log, fixed = TRUE)))
  expect_true(any(grepl("Average of", legacy_log, fixed = TRUE)))
  expect_identical(importance, cv_imr(fit_bin, k = 2L, rounds = 1L,
                                     cv_method = "importance"))
  expect_identical(legacy, cv_imr(fit_bin, k = 2L, rounds = 1L,
                                max_models = 2L, cv_method = "legacy"))
})
