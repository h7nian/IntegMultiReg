marginal_r_reference <- function() {
  reference <- new.env(parent = asNamespace("IntegMultiReg"))
  root <- testthat::test_path("..", "reference", "marginal-r")
  for (file in c("helpers.R", "variance.R", "coefficients.R")) {
    sys.source(file.path(root, file), reference)
  }
  reference
}

test_that("native marginal transitions replay the frozen R algorithms", {
  reference <- marginal_r_reference()
  for (choice in c("variance", "coefficients")) {
    name <- paste0(".imr_", choice, "_marginal_sample")
    for (case in c("continuous", "correlated-swap", "binary", "right.censored", "availability-imr", "availability-bms")) {
      availability <- startsWith(case, "availability")
      groups <- readRDS(testthat::test_path(
        "..", "reference", "marginal-r",
        paste0(if (availability) "availability" else case, "-data.rds")
      ))
      arguments <- list(
        groups = groups, feature_platform = if (availability) c(1L, 2L) else c(1L, 1L),
        nu = if (availability) c(-1, -1.5) else -1,
        draws = 80L, burnin = 20L, thin = 2L, seed = 412L,
        model_variant = if (case == "availability-bms") "bms" else "imr",
        initial = "alternating", keep_latent = TRUE
      )
      expected <- do.call(get(name, reference), arguments)
      expected_rng <- .Random.seed
      actual <- do.call(get(name, asNamespace("IntegMultiReg")), arguments)
      expect_true(identical(actual, expected, num.eq = FALSE), info = paste(choice, case))
      expect_identical(.Random.seed, expected_rng, info = paste(choice, case))
    }
  }
})

test_that("both native marginal samplers report progress without changing samples", {
  d <- small_data("binary")
  for (choice in c("variance", "coefficients")) {
    fit <- function(verbose) {
      imr(d$platforms, d$outcome,
        outcome_type = "binary", marginalize = choice, min_subgroup_size = 0,
        priors = imr_priors(forced_scale = 1), verbose = verbose,
        mcmc = imr_mcmc(
          draws = 4, burnin = 2, chains = 1, seed = 117,
          diagnostics = FALSE, keep_latent = TRUE
        )
      )
    }
    quiet <- fit(FALSE)
    progress <- capture.output(reported <- fit(TRUE))
    expect_identical(progress, sprintf("Iteration %d of 6", 1:6))
    expect_true(identical(reported$posterior, quiet$posterior, num.eq = FALSE))
  }
})

test_that("native preparation preserves named priors and numeric index counts", {
  reference <- marginal_r_reference()
  groups <- readRDS(testthat::test_path("..", "reference", "marginal-r", "continuous-data.rds"))
  groups <- lapply(groups, function(g) {
    g$n_forced <- as.double(g$n_forced)
    g$residual_prior <- g$residual_prior[c("rate", "shape")]
    g
  })
  for (choice in c("variance", "coefficients")) {
    name <- paste0(".imr_", choice, "_marginal_sample")
    arguments <- list(groups = groups, feature_platform = c(1L, 1L), nu = -1,
      draws = 20L, burnin = 5L, seed = 99L)
    expected <- do.call(get(name, reference), arguments)
    expected_rng <- .Random.seed
    actual <- do.call(get(name, asNamespace("IntegMultiReg")), arguments)
    expect_true(identical(actual, expected, num.eq = FALSE))
    expect_identical(.Random.seed, expected_rng)
  }
})
