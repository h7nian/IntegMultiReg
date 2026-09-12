test_that("interleaved prediction aliases preserve draw-wise results", {
  for (outcome_type in c("continuous", "binary", "right.censored")) {
    fit <- fit_demo(outcome_type, total = 20, burn = 10, seed = 71)
    first <- fit$posterior$selection_draws[[1L]]
    second <- first
    second[[1L]][1L] <- 1L - second[[1L]][1L]
    # An artificial state sequence, not a new fitted posterior: exercise aliases
    # separated by owned rows while retaining each draw's averaging contribution.
    fit$posterior$selection_draws <- rep(list(first, first, second, first), 5L)
    native <- function(cache_bytes) {
      IntegMultiReg:::.imr_call_cv_postfit_native(
        fit, k = 3L, rounds = 2L, max_models = 20L,
        verbose = FALSE, importance = TRUE, cache_bytes = cache_bytes,
        stage = "predict")
    }
    set.seed(719)
    rng <- .Random.seed
    expected <- native(0)
    expect_identical(native(128 * 1024^2), expected)
    expect_identical(native(1024), expected)
    expect_identical(.Random.seed, rng)
  }
})

test_that("a single owned prediction row supports zero burn-in", {
  for (outcome_type in c("continuous", "binary", "right.censored")) {
    fit <- fit_demo(outcome_type, total = 1, burn = 0, seed = 71)
    native <- function(cache_bytes) {
      IntegMultiReg:::.imr_call_cv_postfit_native(
        fit, k = 3L, rounds = 1L, max_models = 1L,
        verbose = FALSE, importance = TRUE, cache_bytes = cache_bytes,
        stage = "predict")
    }
    expect_identical(native(128 * 1024^2), native(0))
  }
})
