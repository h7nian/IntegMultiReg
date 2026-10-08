test_that("the four targets expose the parameter draws they actually store", {
  d <- small_data()
  for (choice in c("none", "variance", "coefficients", "coefficients_and_variance")) {
    fit <- imr(d$platforms, d$outcome,
      outcome_type = "continuous",
      marginalize = choice, min_subgroup_size = 0,
      priors = imr_priors(forced_scale = 1),
      mcmc = imr_mcmc(
        draws = 40, burnin = 20, chains = 2, seed = 18,
        initial = NULL, diagnostics = FALSE
      )
    )
    expect_true(validate_imr_object(fit))
    expect_true(all(is.finite(unlist(inclusion_probabilities(fit)))))
    expect_true(all(is.finite(unlist(lapply(predict(fit, d$platforms), function(x) x$prediction)))))
    if (choice == "coefficients_and_variance") {
      expect_s3_class(fit, "imr_selection")
      malformed <- fit
      class(malformed) <- "imr"
      expect_error(validate_imr_object(malformed), "class and marginalization")
      malformed <- fit
      malformed$posterior$latent_mean[[1]] <- numeric()
      expect_error(predict(malformed, d$platforms), "response means")
      malformed <- fit
      malformed$posterior$selection_mean[[1]] <- matrix(0, 1, 0)
      expect_error(validate_imr_object(malformed), "mean dimensions")
      expect_null(fit$posterior$coefficients)
      expect_null(fit$posterior$variance)
      expect_error(coef(fit), "no stored regression coefficients")
      expect_identical(posterior_draws(fit), posterior_draws(fit, "selection"))
      expect_error(cv_imr(fit, cv_method = "reweight"), "requires stored regression")
    } else {
      expect_true(all(is.finite(unlist(coef(fit)))))
      expect_true(nrow(confint(fit)) > 0)
      expect_identical(
        inclusion_probabilities(fit)[[1]][1, 1],
        mean(posterior_draws(fit)[[1]][, , 2] != 0)
      )
    }
  }
})

test_that("marginal samplers keep thinning independent of auxiliary recovery", {
  d <- small_data("binary")
  retained <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.list(x)) {
      return(lapply(x, retained))
    }
    if (length(dim(x)) == 2L) x[seq(3L, 12L, 3L), , drop = FALSE] else x[seq(3L, 12L, 3L), , , drop = FALSE]
  }
  for (choice in c("variance", "coefficients")) {
    fit <- function(draws, thin) {
      imr(d$platforms, d$outcome,
        outcome_type = "binary",
        marginalize = choice, min_subgroup_size = 0,
        priors = imr_priors(forced_scale = 1),
        mcmc = imr_mcmc(
          draws = draws, thin = thin, burnin = 4, chains = 2,
          seed = 25, keep_latent = TRUE, diagnostics = FALSE
        )
      )
    }
    expect_identical(retained(fit(12, 1)$posterior), fit(4, 3)$posterior)
  }
})

test_that("unsupported integration sizes and irrelevant proposal overrides fail early", {
  d <- small_data()
  expect_error(imr(d$platforms, d$outcome,
    outcome_type = "continuous",
    marginalize = "coefficients", min_subgroup_size = 0,
    control = imr_control(max_integration_nodes = 2)
  ), "exceeding")
  expect_error(imr(d$platforms, d$outcome,
    outcome_type = "continuous",
    mcmc = imr_mcmc(variance_step = .2)
  ), "only when coefficients")
  expect_error(imr(d$platforms, d$outcome,
    outcome_type = "continuous",
    marginalize = "coefficients_and_variance", mcmc = imr_mcmc(theta_step = .8)
  ), "retains its proposal")
})

test_that("conditional regression retains source selection states and separate diagnostics", {
  d <- small_data()
  fit <- imr(d$platforms, d$outcome,
    outcome_type = "continuous",
    marginalize = "coefficients_and_variance", min_subgroup_size = 0,
    priors = imr_priors(forced_scale = 1),
    mcmc = imr_mcmc(
      draws = 100, burnin = 40, chains = 1, initial = NULL,
      seed = 17, diagnostics = FALSE
    )
  )
  set.seed(513)
  rng <- .Random.seed
  # This short run tests representation and pairing, not convergence.
  messages <- character()
  posterior <- withCallingHandlers(
    sample_regression_posterior(fit,
      output_draws = 60, burnin = 40, min_draws_per_model_chain = 40, seed = 19
    ),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(all(grepl("Conditional-chain diagnostics|ESS has been capped", messages)))
  expect_output(print(posterior), "Conditional regression posterior")
  expect_true(validate_imr_object(posterior))
  source <- fit$posterior$selection[[1]][posterior$selection_draw_index, 1, 1]
  expect_identical(as.integer(posterior$coefficients[[1]][, 2] != 0), as.integer(source))
  expect_true(all(c("selection_model", "draws_per_model_chain", "rhat", "ess_bulk") %in% names(mcmc_diagnostics(posterior))))
  expect_identical(dimnames(posterior_draws(posterior)[[1]])[[2]], "mixture")
  expect_true(all(is.finite(unlist(coef(posterior)))))
  expect_true(nrow(confint(posterior)) > 0)
  expect_true(all(is.finite(predict(posterior, d$platforms)[[1]]$prediction)))
  expect_identical(.Random.seed, rng)
})

test_that("marginal samplers replay starts and parallel chains", {
  d <- small_data("right.censored")
  for (choice in c("variance", "coefficients", "coefficients_and_variance")) {
    fit <- function(workers = 1L, initial = "dispersed") {
      imr(d$platforms, d$outcome,
        outcome_type = "right.censored", marginalize = choice,
        min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
        mcmc = imr_mcmc(
          draws = 8, burnin = 4, chains = 2,
          seed = 58, workers = workers, initial = initial,
          keep_latent = TRUE, diagnostics = FALSE
        )
      )
    }
    serial <- fit()
    expect_identical(serial$posterior, fit(workers = 2)$posterior)
    expect_identical(serial$posterior, fit(initial = serial$control$initial)$posterior)
    cv <- cv_imr(serial, k = 2, rounds = 1, seed = 71)
    expect_identical(cv$control$marginalize, choice)
    expect_identical(cv$control$numerical, serial$control$numerical)
    expect_true(all(is.finite(cv$predictions$prediction)))
    initial <- serial$control$initial
    initial[[1]]$coefficients[[1]][1] <- 0
    expect_error(fit(initial = initial), "Initial forced coefficients")
    initial <- serial$control$initial
    initial[[1]]$interaction[[1]][1, 1] <- 1
    expect_error(fit(initial = initial), "symmetric with zero diagonal")
  }
})
