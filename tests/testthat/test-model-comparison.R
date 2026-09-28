test_that("prediction counts complete selection states and preserves its cutoff", {
  state <- fit_bin$posterior$selection_draws[[1L]]
  expect_gte(length(state[[1L]]), 8L)
  unique_states <- lapply(0:255, function(index) {
    value <- state
    value[[1L]][1:8] <- as.integer(intToBits(index)[1:8])
    value
  })
  patterns <- list(
    repeated = rep(list(state), 300L),
    unique = unique_states[rep(1:256, length.out = 300L)],
    mixed = rep(unique_states[c(1L, 128L, 256L)], 100L)
  )
  for (states in patterns) {
    fit <- fit_bin
    fit$posterior$selection_draws <- states
    keys <- vapply(states, function(value) paste(unlist(value), collapse = ""), "")
    for (limit in c(1L, 2L, 3L)) {
      set.seed(792)
      rng <- .Random.seed
      transcript <- capture.output(prediction <- predict(
        fit, simIMR$platforms, covariates = simIMR$covariates,
        max_models = limit, verbose = TRUE
      ))
      count_line <- grep("nbr of models", transcript, value = TRUE, fixed = TRUE)
      expect_length(count_line, 1L)
      count <- as.integer(sub(".*= *", "", count_line))
      expect_identical(count, as.integer(min(length(unique(keys)), 100L * limit)))
      expect_identical(.Random.seed, rng)
      expect_true(all(is.finite(unlist(lapply(prediction, `[[`, "prediction")))))
    }
  }
})
