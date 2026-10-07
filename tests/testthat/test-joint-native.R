test_that("native numerical errors restore RNG ownership and permit another fit", {
  dat <- small_data()
  set.seed(1904)
  before <- .Random.seed
  expect_error(
    imr(dat$platforms, dat$outcome,
      outcome_type = "continuous", min_subgroup_size = 0,
      priors = imr_priors(molecular_scale = 1e-320, forced_scale = 1),
      mcmc = imr_mcmc(draws = 4, burnin = 0, chains = 1, seed = 1, diagnostics = FALSE)
    ),
    "Non-finite conditional Bayes factor"
  )
  expect_identical(.Random.seed, before)
  expect_identical(small_fit()$posterior, small_fit()$posterior)
})

test_that("short-chain progress reports every iteration without altering draws", {
  dat <- small_data()
  args <- list(
    x = dat$platforms, outcome = dat$outcome, outcome_type = "continuous",
    min_subgroup_size = 0, mcmc = imr_mcmc(draws = 4, burnin = 2, chains = 1, seed = 2, diagnostics = FALSE)
  )
  quiet <- do.call(imr, args)
  output <- capture.output(loud <- do.call(imr, c(args, list(verbose = TRUE))))
  expect_identical(grep("^Iteration", output, value = TRUE), paste("Iteration", 1:6, "of 6"))
  expect_identical(quiet$posterior, loud$posterior)
})

test_that("R allocations remain protected under frequent garbage collection", {
  dat <- small_data("binary", n = 12)
  args <- list(
    x = dat$platforms, outcome = dat$outcome, outcome_type = "binary",
    min_subgroup_size = 0, mcmc = imr_mcmc(
      draws = 4, burnin = 0, chains = 1,
      keep_latent = TRUE, seed = 73, diagnostics = FALSE
    )
  )
  baseline <- do.call(imr, args)
  previous <- gctorture2(10L)
  on.exit(gctorture2(previous), add = TRUE)
  frequent_gc <- do.call(imr, args)
  gctorture2(previous)
  expect_identical(frequent_gc$posterior, baseline$posterior)
  expect_true(validate_imr_object(frequent_gc))
})

test_that("identity-scale censoring is not bounded by the retired truncation constant", {
  dat <- small_data("right.censored")
  dat$outcome[[2L]] <- dat$outcome[[2L]] + 2000
  fit <- imr(dat$platforms, dat$outcome,
    outcome_type = "right.censored",
    survival_scale = "identity", min_subgroup_size = 0,
    mcmc = imr_mcmc(
      draws = 4, burnin = 4, chains = 1, keep_latent = TRUE,
      seed = 171, diagnostics = FALSE
    )
  )
  expect_true(validate_imr_object(fit))
  expect_true(all(posterior_draws(fit, "latent")[[1L]] > 1000.5))
})
