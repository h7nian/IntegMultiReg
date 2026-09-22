test_that("repeated selection states retain multiplicity and copy-on-modify", {
  id <- 1:24
  fit <- imr(list(a = data.frame(id, x = sin(id))),
    data.frame(id, y = cos(id)), outcome_type = "continuous",
    draws = 30, burnin = 5, min_subgroup_size = 0, seed = 71)
  states <- fit$posterior$selection_draws
  expect_length(states, 30L)
  matrices <- lapply(states, `[[`, 1L)
  repeated <- which(duplicated(matrices))[1L]
  first <- which(vapply(matrices, identical, logical(1), matrices[[repeated]]))[1L]
  expect_lt(first, repeated)
  if (capabilities("profmem")) {
    a <- tracemem(matrices[[first]])
    b <- tracemem(matrices[[repeated]])
    untracemem(matrices[[first]])
    untracemem(matrices[[repeated]])
    expect_identical(a, b)
  }
  original <- states[[first]][[1L]][1L, 1L]
  states[[repeated]][[1L]][1L, 1L] <- 1L - original
  expect_identical(states[[first]][[1L]][1L, 1L], original)
  expect_identical(fit$posterior$selection_draws[[repeated]][[1L]][1L, 1L], original)
  expect_identical(states[[repeated]][[1L]][1L, 1L], 1L - original)
  # Serialization remains a complete, ordinary list of integer matrices.
  expect_identical(unserialize(serialize(fit$posterior$selection_draws, NULL)),
    fit$posterior$selection_draws)
})
