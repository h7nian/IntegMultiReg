test_that("bounded post-fit caches preserve every draw contribution", {
  for (outcome_type in c("continuous", "binary", "right.censored")) {
    fit <- fit_demo(outcome_type, total = 20, burn = 10, seed = 71)
    native <- function(bytes, bits = 64L) {
      IntegMultiReg:::.imr_call_cv_postfit_native(
        fit, k = 3L, rounds = 2L, max_models = 20L,
        verbose = FALSE, importance = TRUE,
        cache_bytes = bytes, cache_hash_bits = bits)
    }
    reference <- native(0)
    expect_identical(native(128 * 1024^2), reference)
    expect_identical(native(32), reference) # insufficient index budget: no cache
    expect_identical(native(256), reference) # bounded table, then uncached misses
    expect_identical(native(1024, 0L), reference) # deliberate hash collisions
    fit$posterior$selection_draws <- rep(fit$posterior$selection_draws[1L], 20L)
    expect_identical(native(128 * 1024^2), native(0))
  }
})

test_that("internal cache controls cannot exceed the memory bound", {
  call <- function(bytes, bits = 64L) {
    IntegMultiReg:::.imr_call_cv_postfit_native(
      fit_bin, 2L, 1L, 2L, FALSE, TRUE,
      cache_bytes = bytes, cache_hash_bits = bits)
  }
  expect_error(call(-1), "cache_bytes")
  expect_error(call(128 * 1024^2 + 1), "cache controls")
  expect_error(call(32, 65L), "cache controls")
})
