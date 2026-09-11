test_that("validate_imr checks fitted object integrity", {
  expect_true(validate_imr(fit_bin))
  bad <- fit_bin
  bad$posterior$inclusion_probabilities[[1L]][1L, 1L] <- 2
  expect_error(validate_imr(bad), "invalid")

  bad_draws <- fit_bin
  bad_draws$posterior$selection_draws[[1L]][[1L]] <- matrix(0, 1, 1)
  expect_error(validate_imr(bad_draws), "inconsistent")

  missing_metadata <- fit_bin
  missing_metadata$control <- NULL
  expect_error(validate_imr(missing_metadata), "missing `control`")

  malformed_mcmc <- fit_bin
  malformed_mcmc$control$mcmc <- list(burnin = 2L)
  expect_error(validate_imr(malformed_mcmc), "valid `draws` and `burnin`")

  bad_theta <- fit_bin
  bad_theta$posterior$interaction_draws[[1L]] <- bad_theta$posterior$interaction_draws[[1L]][, -1L, drop = FALSE]
  expect_error(validate_imr(bad_theta), "Interaction draws")
})

test_that("posterior summaries include uncertainty for gamma and theta", {
  out <- posterior_summary(fit_bin, level = 0.9)
  expect_s3_class(out, "posterior_summary.imr")
  expect_named(out, c("level", "selection", "theta"))
  expect_true(all(c("mean", "sd", "lower", "median", "upper") %in%
                    names(out$selection$genomic)))
  expect_true(nrow(out$selection$genomic) > 0L)
  expect_true(nrow(out$theta$genomic) > 0L)
  expect_output(print(out), "posterior summary")
  ci <- confint(fit_bin, parm = "theta", level = 0.9)
  expect_named(ci, fit_bin$model$platform_names)
})

test_that("compare_imr creates a common descriptive table", {
  out <- compare_imr(binary = fit_bin, copy = fit_bin, threshold = 0.5)
  expect_equal(nrow(out), 2L)
  expect_named(out, c(
    "fit", "outcome", "method", "platforms", "subgroups",
    "retained_draws", "selected_features"
  ))

  incompatible <- fit_bin
  incompatible$control$outcome_type <- "continuous"
  expect_error(compare_imr(fit_bin, incompatible), "same outcome type")
})

test_that("additional parameter trace plots run", {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  expect_invisible(plot(fit_bin, type = "theta_trace", platform = 1,
                        parameter = 1))
  expect_invisible(plot(fit_bin, type = "selection_trace", platform = 1,
                        subgroup = 1, feature = 1))
  grDevices::dev.off()
  expect_true(file.exists(file))
})
