test_that("the paper sampler is explicit, recorded, and preserved by refits", {
  args <- list(x = simIMR$platforms, outcome = simIMR$outcome.continuous,
    outcome_type = "continuous", covariates = simIMR$covariates,
    draws = 12, burnin = 6, min_subgroup_size = 30, seed = 53)
  legacy <- do.call(imr, args)
  explicit <- do.call(imr, c(args, list(sampler_method = "legacy")))
  expect_identical(legacy$posterior, explicit$posterior)
  paper <- do.call(imr, c(args, list(sampler_method = "paper")))
  expect_identical(paper$control$sampler_method, "paper")
  expect_true(validate_imr(paper))
  expect_true(all(is.finite(paper$posterior$log_posterior)))
  train_ids <- unlist(lapply(paper$model$subgroup_names, function(group) {
    ids <- paper$preprocessing$input_data$availability
    head(ids$id[ids$subgroup == group], -1L)
  }), use.names = FALSE)
  refit <- IntegMultiReg:::.imr_cv_refit(paper, train_ids, seed = 71)
  expect_identical(refit$control$sampler_method, "paper")
  old <- legacy; old$control$sampler_method <- NULL
  expect_true(validate_imr(old))
  expect_identical(IntegMultiReg:::.imr_cv_refit(old, train_ids, seed = 71)$control$sampler_method,
                   "legacy")
  bad <- paper; bad$control$sampler_method <- "unknown"
  expect_error(validate_imr(bad), "sampler_method")
  expect_error(do.call(imr, c(args, list(sampler_method = "unknown"))), "arg")
})
