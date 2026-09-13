test_that("post-fit solve failures stop cleanly without partial results", {
  id <- seq_len(120L)
  fit <- imr(list(assay = data.frame(id = id, a = sin(id), b = cos(id))),
    data.frame(id = id, y = sin(id / 3)), outcome_type = "continuous",
    method = "bms", draws = 4L, burnin = 2L,
    min_subgroup_size = 0L, seed = 31L)
  zero <- fit$posterior$selection_draws[[1L]]
  zero[[1L]][] <- 0L
  both <- zero
  both[[1L]][] <- 1L
  fit$posterior$selection_draws <- list(zero, both, zero, both)
  broken <- fit
  # Deliberately corrupt saved preprocessing while preserving schema dimensions.
  # At this scale the duplicated columns lose the ridge in floating point.
  broken$preprocessing$features[[1L]][[1L]][, 2L] <-
    broken$preprocessing$features[[1L]][[1L]][, 1L]
  broken$preprocessing$features[[1L]][[1L]] <-
    broken$preprocessing$features[[1L]][[1L]] * 1e8
  expect_true(validate_imr(broken))
  set.seed(731)
  rng <- .Random.seed
  # Force the bounded path as well as the ordinary row path on the same object.
  for (importance in c(FALSE, TRUE)) {
    native_failure <- function(cache_bytes) tryCatch(
      IntegMultiReg:::.imr_call_cv_postfit_native(
        broken, k = 2L, rounds = 2L, max_models = 4L, verbose = FALSE,
        importance = importance, cache_bytes = cache_bytes), error = identity)
    # Floating-point tools can reject different folds of this near-singular
    # fixture. Require every budget to match this platform's ordinary path.
    reference_failure <- native_failure(128 * 1024^2)
    expect_s3_class(reference_failure, "error")
    reference_message <- conditionMessage(reference_failure)
    expect_match(reference_message,
      "CV Cholesky solve failed \\(round [12], fold [12], subgroup 1\\)")
    for (cache_bytes in c(0, 256, 4096, 128 * 1024^2)) {
      expect_identical(conditionMessage(native_failure(cache_bytes)), reference_message)
      expect_identical(.Random.seed, rng)
    }
  }
  # R CMD check permits at most two simultaneous workers. Unrestricted tests
  # (including the native sanitizer suite) also exercise three workers.
  limited <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  worker_counts <- if (nzchar(limited) && limited != "false") 1:2 else 1:3
  for (mode in c("legacy", "importance")) {
    expected <- cv_imr(fit, k = 2L, rounds = 2L, cv_method = mode)
    for (workers in worker_counts) {
      connections <- rownames(showConnections(all = TRUE))
      expect_error(cv_imr(broken, k = 2L, rounds = 2L,
        cv_method = mode, workers = workers),
        "CV Cholesky solve failed \\(round [12], fold [12], subgroup 1\\)")
      expect_identical(.Random.seed, rng)
      expect_identical(rownames(showConnections(all = TRUE)), connections)
      expect_identical(cv_imr(fit, k = 2L, rounds = 2L,
        cv_method = mode, workers = workers), expected)
      expect_identical(.Random.seed, rng)
    }
  }
})
