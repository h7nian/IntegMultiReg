test_that("one fit supplies joint parameters and coherent exclusion zeros", {
  f <- fit_bin
  expect_true(validate_imr_object(f))
  expect_identical(f$schema_version, 4L)
  expect_named(f$posterior, c("coefficients", "variance", "interaction", "latent", "log_density"))
  expect_false(any(c("selection_draws", "inclusion_probabilities") %in% names(f$posterior)))
  for (g in seq_along(f$posterior$coefficients)) {
    b <- f$posterior$coefficients[[g]]
    expect_identical(dim(b)[1:2], c(40L, 2L))
    expect_true(all(b[, , 1:4] != 0))
    expect_equal(unname(coef(f)[[g]]), vapply(seq_len(dim(b)[3L]), function(j) mean(b[, , j]), 0))
  }
  selection <- posterior_draws(f, parm = "selection")
  probabilities <- inclusion_probabilities(f)
  for (l in seq_along(selection)) {
    expect_true(all(selection[[l]] %in% 0:1))
    calculated <- vapply(seq_len(dim(selection[[l]])[3L]), function(j) mean(selection[[l]][, , j]), 0)
    expect_equal(as.vector(t(probabilities[[l]])), unname(calculated))
  }
})

test_that("same seeds and recorded initial states replay the entire chain", {
  f <- small_fit(chains = 4)
  same <- small_fit(chains = 4)
  expect_identical(f$posterior, same$posterior)
  replay <- small_fit(chains = 4, initial = f$control$initial)
  expect_identical(f$posterior, replay$posterior)
  expect_identical(f$control$chain_seeds, replay$control$chain_seeds)
  expect_identical(f$control$initial_seeds, replay$control$initial_seeds)
  set.seed(178)
  before <- .Random.seed
  invisible(small_fit())
  expect_identical(.Random.seed, before)
  expect_false(identical(f$posterior, small_fit(chains = 4, seed = 18)$posterior))
})

test_that("latent storage changes storage only and respects all likelihood supports", {
  for (type in c("binary", "right.censored")) {
    a <- small_fit(type)
    b <- small_fit(type, keep_latent = TRUE)
    expect_identical(a$posterior$coefficients, b$posterior$coefficients)
    expect_identical(a$posterior$variance, b$posterior$variance)
    expect_true(validate_imr_object(b))
    z <- posterior_draws(b, "latent")[[1L]]
    y <- b$preprocessing$response[[1L]]
    for (i in seq_len(nrow(y))) {
      if (type == "binary") {
        expect_true(all(if (y[i, 1] == 1) z[, , i] >= 0 else z[, , i] <= 0))
      } else if (y[i, 2] == 1) {
        expect_true(all(z[, , i] == y[i, 1]))
      } else {
        expect_true(all(z[, , i] >= y[i, 1]))
      }
    }
  }
})

test_that("thinning preserves every retained state for all outcome families", {
  take <- seq(3L, 12L, 3L)
  retained <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.list(x)) {
      return(lapply(x, retained))
    }
    if (length(dim(x)) == 2L) x[take, , drop = FALSE] else x[take, , , drop = FALSE]
  }
  for (type in c("continuous", "binary", "right.censored")) {
    d <- small_data(type)
    common <- list(
      x = d$platforms, outcome = d$outcome, outcome_type = type,
      min_subgroup_size = 0, priors = imr_priors(forced_scale = 1)
    )
    a <- do.call(imr, c(common, list(mcmc = imr_mcmc(
      draws = 12, burnin = 4, chains = 2, seed = 22,
      keep_latent = TRUE, diagnostics = FALSE
    ))))
    b <- do.call(imr, c(common, list(mcmc = imr_mcmc(
      draws = 4, burnin = 4, chains = 2, thin = 3, seed = 22,
      keep_latent = TRUE, diagnostics = FALSE
    ))))
    expect_identical(retained(a$posterior), b$posterior)
  }
})

test_that("draw budgets, priors and old APIs fail before expensive computation", {
  expect_error(imr_mcmc(draws = 3), "draws")
  expect_error(imr_mcmc(draws = .Machine$integer.max, chains = 2), "draw count")
  expect_error(imr_mcmc(burnin = .Machine$integer.max), "iteration limit")
  expect_error(imr_priors(molecular_scale = 0), "molecular_scale")
  expect_error(imr_priors(residual = c(a = 1, b = 1)), "shape")
  d <- small_data()
  expect_error(imr(d$platforms, d$outcome, outcome_type = "continuous", draws = 10), "Unused argument")
  expect_error(imr(d$platforms, d$outcome,
    outcome_type = "continuous", min_subgroup_size = 0,
    mcmc = imr_mcmc(max_draw_memory_mb = 1e-5)
  ), "Retained posterior arrays")
  old <- fit_bin
  old$schema_version <- 3L
  expect_error(coef(old), "Earlier selection-only fits")
  bad <- fit_bin
  bad$posterior$variance[1] <- -1
  expect_error(validate_imr_object(bad), "variance draw values")
  bad <- fit_bin
  bad$posterior$coefficients[[1]][1, 1, 1] <- 0
  expect_error(validate_imr_object(bad), "forced coefficient")
  retired <- c("sample_regression_posterior", "selection_summary", "upgrade_imr_object")
  expect_false(any(retired %in% getNamespaceExports("IntegMultiReg")))
})


test_that("integer and numeric prior inputs specify the same native model", {
  integer_priors <- imr_priors(
    molecular_scale = 1L, forced_scale = 1L,
    residual = c(shape = 100000L, rate = 100000L), interaction = c(shape = 2L, rate = 2L)
  )
  numeric_priors <- imr_priors(
    molecular_scale = 1, forced_scale = 1,
    residual = c(shape = 1e5, rate = 1e5), interaction = c(shape = 2, rate = 2)
  )
  expect_identical(integer_priors, numeric_priors)
  d <- small_data("binary")
  args <- list(
    x = d$platforms, outcome = d$outcome, outcome_type = "binary",
    min_subgroup_size = 0, mcmc = imr_mcmc(
      draws = 4, burnin = 0, chains = 1,
      seed = 70, diagnostics = FALSE
    )
  )
  a <- do.call(imr, c(args, list(priors = integer_priors)))
  b <- do.call(imr, c(args, list(priors = numeric_priors)))
  expect_identical(a$posterior, b$posterior)
})


test_that("explicit starts and saved metadata reject ambiguous subgroup labels", {
  fit <- small_fit()
  initial <- fit$control$initial
  dimnames(initial[[1]]$interaction[[1]]) <- list("wrong", "wrong")
  expect_error(small_fit(initial = initial), "Initial interaction matrices")
  initial <- fit$control$initial
  initial[[1]]$coefficients[[1]][] <- "1"
  expect_error(small_fit(initial = initial), "initial coefficients")
  bad <- fit
  bad$control$initial_seeds <- NULL
  expect_error(validate_imr_object(bad), "chain seeds")
  bad <- fit
  bad$control$priors$nu <- c(-3, -4)
  expect_error(validate_imr_object(bad), "resolved")
  bad <- fit
  dimnames(bad$posterior$variance)[[3L]] <- "wrong"
  expect_error(validate_imr_object(bad), "Variance parameter labels")
  bad <- c(fit, fit["posterior"])
  class(bad) <- "imr"
  expect_error(validate_imr_object(bad), "top-level fields")
})
