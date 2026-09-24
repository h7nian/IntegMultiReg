test_that("tiny positive interaction starts cannot create zero theta states", {
  id <- 1:40
  x <- list(a = data.frame(id, a = sin(id)),
            b = data.frame(id = 1:20, b = cos(1:20)))
  y <- data.frame(id, y = sin(id / 4))
  run <- function(...) imr(x, y, outcome_type = "continuous", draws = 100,
    burnin = 0, min_subgroup_size = 0, seed = 31, ...)
  template <- run()
  initial <- list(interaction = lapply(coef(template), function(m) {
    z <- matrix(1e-8, nrow(m), nrow(m),
                dimnames = list(rownames(m), rownames(m)))
    diag(z) <- 0
    z
  }))
  for (method in c("legacy", "paper")) {
    fit <- run(sampler_method = method, initial = initial)
    expect_true(validate_imr(fit))
    expect_true(all(fit$posterior$interaction_draws[[1]] > 0))
    expect_true(all(is.finite(fit$posterior$log_posterior)))
    bad <- fit
    bad$posterior$interaction_draws[[1]][1, 1] <- 0
    expect_error(validate_imr(bad), "Interaction draws")
    bad <- fit
    bad$posterior$interaction_means[[1]][1, 2] <- 0
    expect_error(validate_imr(bad), "Interaction means")
  }
})
