consistency_fit <- function(choice, type = "continuous", control = imr_control(), ...) {
  d <- small_data(type)
  imr(d$platforms, d$outcome,
    outcome_type = type, marginalize = choice,
    min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
    mcmc = imr_mcmc(draws = 20, burnin = 10, chains = 2, seed = 317, diagnostics = FALSE),
    control = control, ...
  )
}

test_that("all samplers expose the same draw-family structure without implicit sampling", {
  for (choice in c("none", "variance", "coefficients", "coefficients_and_variance")) {
    fit <- consistency_fit(choice)
    set.seed(913)
    seed <- .Random.seed
    draws <- posterior_draws(fit)
    expect_named(draws, c("coefficients", "variance", "selection", "interaction", "latent"))
    expect_identical(.Random.seed, seed)
    expect_identical(draws$selection, posterior_draws(fit, "selection"))
    expect_null(draws$latent)
    if (choice == "coefficients_and_variance") {
      expect_null(draws$coefficients)
      expect_null(draws$variance)
      expect_error(posterior_draws(fit, "coefficients"), "require conditional")
    } else {
      expect_identical(draws$coefficients, posterior_draws(fit, "coefficients"))
    }
    # Tiny chains exercise the interface, not inferential precision.
    result <- withCallingHandlers(summary(fit), warning = function(w) {
      expect_match(conditionMessage(w), "ESS has been capped")
      invokeRestart("muffleWarning")
    })
    expect_s3_class(result, "summary.imr")
    expect_s3_class(confint(fit, parm = "all"), "data.frame")
    sampler <- mcmc_diagnostics(fit, type = "sampler")
    expect_named(sampler, c("chain", "update", "method", "proposed", "accepted", "rate"))
    expect_identical(unique(sampler$update), c("selection", "swap", "interaction", "variance", "latent"))
    expect_true(all(sampler$accepted <= sampler$proposed, na.rm = TRUE))
    expect_true(all(is.na(sampler$rate[sampler$method != "metropolis"])))
    if (choice == "coefficients_and_variance") {
      selection <- sampler[sampler$update == "selection", ]
      expect_identical(selection$proposed, rep(30, 2))
    }
  }
})

test_that("inapplicable numerical controls fail rather than being silently ignored", {
  for (choice in c("none", "variance", "coefficients")) {
    expect_error(consistency_fit(choice, control = imr_control(laplace_tolerance = 1e-8)), "laplace_tolerance.*not applicable")
  }
  for (choice in c("none", "variance", "coefficients_and_variance")) {
    expect_error(consistency_fit(choice, control = imr_control(max_integration_nodes = 100)), "max_integration_nodes.*not applicable")
  }
  fit <- consistency_fit("coefficients", type = "binary")
  expect_named(fit$control$effective$numerical, "max_integration_nodes")
  expect_equal(unname(fit$control$effective$mcmc$variance_step), 1.5 / sqrt(1e5 + 10))
  fit <- consistency_fit("coefficients_and_variance")
  expect_false(any(c("theta_step", "swap_rate", "variance_step") %in% names(fit$control$effective$mcmc)))
  expect_named(fit$control$effective$numerical, c("laplace_max_iter", "laplace_tolerance"))
})

test_that("prediction reports its scale and calculation under the shared interface", {
  d <- small_data("right.censored")
  for (choice in c("none", "variance", "coefficients", "coefficients_and_variance")) {
    fit <- consistency_fit(choice, type = "right.censored")
    response <- predict(fit, d$platforms)
    link <- predict(fit, d$platforms, type = "link")
    expect_identical(attr(response, "prediction")$type, "response")
    expect_identical(attr(link, "prediction")$type, "link")
    expect_false(attr(response, "prediction")$interval)
    expect_true(all(response[[1]]$prediction > 0))
    if (choice == "coefficients_and_variance") {
      expect_identical(response[[1]]$prediction, exp(link[[1]]$prediction))
      expect_identical(attr(response, "prediction")$summary, "transformed_mode_average")
      expect_error(predict(fit, d$platforms, interval = TRUE), "sample_regression_posterior")
      expect_error(predict(fit, d$platforms, quantity = "new_observation"), "sample_regression_posterior")
    } else {
      expect_error(predict(fit, d$platforms, quantity = "model_average"), "only to a Laplace")
    }
  }
  expect_identical(formals(getS3method("predict", "imr")), formals(getS3method("predict", "imr_selection")))
  expect_identical(formals(getS3method("predict", "imr")), formals(getS3method("predict", "imr_posterior")))
})

test_that("conditional computation uses common MCMC and result conventions", {
  fit <- consistency_fit("coefficients_and_variance", type = "binary")
  post <- sample_regression_posterior(fit,
    output_draws = 20,
    mcmc = imr_mcmc(
      draws = 10, burnin = 5, chains = 1, seed = 48,
      thin = 2, keep_latent = TRUE, diagnostics = FALSE
    )
  )
  expect_true(validate_imr_object(post))
  draws <- posterior_draws(post)
  expect_named(draws, c("coefficients", "variance", "selection", "interaction", "latent"))
  expect_null(draws$selection)
  expect_null(draws$interaction)
  expect_identical(dimnames(draws$coefficients[[1]])[[2]], "mixture")
  result <- summary(post)
  expect_s3_class(result, "summary.imr")
  expect_true(is.na(result$chains))
  expect_true(all(c("coefficients", "variance", "latent") %in% result$parameters$family))
  expect_s3_class(confint(post, parm = "all"), "data.frame")
  expect_true(all(is.na(mcmc_diagnostics(post, type = "sampler")$rate)))
  expect_error(mcmc_diagnostics(post), "not retained")
  prediction <- predict(post, small_data("binary")$platforms)
  expect_false(attr(prediction, "prediction")$interval)
  expect_identical(attr(prediction, "prediction")$method, "conditional_posterior_draws")
  expect_error(sample_regression_posterior(fit, mcmc = imr_mcmc(theta_step = .7)), "not applicable")
  expect_error(sample_regression_posterior(fit, mcmc = imr_mcmc(initial = "full")), "fixed-model")
})

test_that("conditional thinning preserves update and RNG order", {
  X <- cbind(1, sin(1:12))
  run <- function(draws, thin) {
    set.seed(121)
    IntegMultiReg:::.imr_conditional_chain(X, as.double(rep(0:1, 6)), c(1, .5), 1e5, 1e5,
      draws = draws, burnin = 3, outcome_type = "binary", keep_latent = TRUE, thin = thin
    )
  }
  full <- run(12, 1)
  seed <- .Random.seed
  thinned <- run(4, 3)
  expect_identical(unname(thinned[, ]), unname(full[c(3, 6, 9, 12), ]))
  expect_identical(attr(thinned, "latent"), attr(full, "latent")[c(3, 6, 9, 12), ])
  expect_identical(.Random.seed, seed)
})
