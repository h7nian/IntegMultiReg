test_that("explicit starts reach the sampler and preserve default behavior", {
  id <- 1:40
  x <- list(a = data.frame(id, v1 = sin(id), v2 = cos(id), v3 = sin(id/2),
    v4 = cos(id/2), v5 = sin(id/3), v6 = cos(id/3)))
  y <- data.frame(id, y = sin(id/4))
  run <- function(...) imr(x, y, outcome_type = "continuous", method = "bms",
    draws = 1L, burnin = 0L, min_subgroup_size = 0L, seed = 31L, ...)
  base <- run()
  expect_identical(run(initial = NULL)$posterior, base$posterior)
  start <- coef(base)
  start[[1]][] <- 0L
  zero <- run(initial = list(selection = start))
  start[[1]][] <- 1L
  full <- run(initial = list(selection = start))
  expect_lte(sum(zero$posterior$selection_draws[[1]][[1]]), 1L)
  expect_gte(sum(full$posterior$selection_draws[[1]][[1]]), 5L)
  expect_identical(full$posterior, run(initial = list(selection = start))$posterior)
  expect_true(all(full$control$initial$selection[[1]] == 1L))
  expect_error(run(initial = list(unknown = start)), "initial")
  invalid <- start; invalid[[1]][1, 1] <- .5
  expect_error(run(initial = list(selection = invalid)), "zero or one")
  invalid <- start; colnames(invalid[[1]]) <- rev(colnames(invalid[[1]]))
  expect_error(run(initial = list(selection = invalid)), "dimnames")
  expect_error(run(initial = list(interaction = start)), "BMS")
})

test_that("interaction starts are validated and refits inherit explicit starts", {
  id <- 1:40
  x <- list(a = data.frame(id, a = sin(id)), b = data.frame(id = 1:20, b = cos(1:20)))
  y <- data.frame(id, y = sin(id/4))
  run <- function(...) imr(x, y, outcome_type = "continuous", draws = 4L,
    burnin = 2L, min_subgroup_size = 0L, seed = 31L, ...)
  template <- run()
  initial <- list(selection = coef(template))
  initial$interaction <- lapply(initial$selection, function(m) {
    matrix(0, nrow(m), nrow(m), dimnames = list(rownames(m), rownames(m)))
  })
  initial$selection <- lapply(initial$selection, function(m) {m[] <- 0; m})
  initial$interaction <- lapply(initial$interaction, function(m) {m[] <- .2; diag(m) <- 0; m})
  fit <- run(initial = initial)
  expect_identical(fit$control$initial$interaction, initial$interaction)
  expect_identical(fit$posterior, run(initial = initial)$posterior)
  bad <- initial; bad$interaction[[1]][1,2] <- -.1
  expect_error(run(initial = bad), "symmetric")
  expected <- cv_imr(fit, k = 2, rounds = 1, cv_method = "refit")
  expect_identical(cv_imr(fit, k = 2, rounds = 1, cv_method = "refit", workers = 2), expected)
})
