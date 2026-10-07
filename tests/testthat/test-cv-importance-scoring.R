test_that("both reweighting collections use standard scores and avoid historical output", {
  for (collection in c("all_draws", "top_unique")) {
    log <- capture.output(result <- cv_imr(fit_bin, k = 2L, rounds = 1L,
      max_models = 2L, cv_method = "reweight", model_set = collection, verbose = TRUE))
    expect_false(any(grepl("Average of", log, fixed = TRUE)))
    expect_identical(result$control$score_rule, "auc_half_ties")
    expect_identical(result, cv_imr(fit_bin, k = 2L, rounds = 1L,
      max_models = 2L, cv_method = "reweight", model_set = collection))
  }
})
