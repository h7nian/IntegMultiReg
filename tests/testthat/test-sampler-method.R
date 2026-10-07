test_that("fitting and refitting use symmetric MRF Hastings updates", {
  fit <- fit_demo("continuous", total = 12, burn = 6, seed = 53)
  expect_identical(fit$control$selection_update, "symmetric_mrf_hastings")
  expect_identical(fit$control$numerical$prior_indexing, "coefficient_blocks")
  expect_true(validate_imr_object(fit))
  expect_true(all(is.finite(fit$posterior$log_posterior)))
  train_ids <- unlist(lapply(fit$model$subgroup_names, function(group) {
    ids <- fit$preprocessing$input_data$availability
    head(ids$id[ids$subgroup == group], -1L)
  }), use.names = FALSE)
  refit <- IntegMultiReg:::.imr_cv_refit(fit, train_ids, seed = 71)
  expect_identical(refit$control$selection_update, fit$control$selection_update)
  bad <- fit; bad$control$selection_update <- "unknown"
  expect_error(validate_imr_object(bad), "selection-update")
  expect_error(imr(simIMR$platforms, simIMR$outcome.continuous, outcome_type = "continuous", sampler_method = "paper"), "Unused argument")
  expect_error(imr(simIMR$platforms, simIMR$outcome.continuous, outcome_type = "continuous", prior_indexing = "code2017"), "Unused argument")
})
