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
