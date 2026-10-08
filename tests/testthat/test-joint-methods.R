test_that("summaries and intervals all refer to the stored joint draws", {
  f <- small_fit()
  set.seed(90)
  rng <- .Random.seed
  beta <- posterior_draws(f, "coefficients")[[1L]]
  ci <- confint(f)
  for (j in seq_len(dim(beta)[3L])) {
    q <- quantile(beta[, , j], c(.025, .5, .975), type = 1, names = FALSE)
    expect_equal(unname(as.numeric(ci[j, c("lower", "median", "upper")])), q)
  }
  variance <- confint(f, parm = "variance")
  expect_equal(variance$mean, mean(f$posterior$variance))
  expect_identical(confint(f, parm = "assay:marker")$parameter, "assay:marker")
  indicators <- confint(f, parm = "selection")
  expect_true(all(c(indicators$lower, indicators$median, indicators$upper) %in% 0:1))
  expect_s3_class(summary(f), "summary.imr")
  expect_identical(.Random.seed, rng)
  expect_output(print(f), "joint posterior fit")
  expect_output(print(summary(f), max_rows = 1), "additional rows")
  expect_error(coef(f, unused = 1), "Unused")
  expect_error(confint(f, parm = "unknown"), "Unknown coefficient")
  expect_error(summary(f, level = 1), "less than one")
})

test_that("rank diagnostics use real chains and leave constant states undefined", {
  f <- small_fit(draws = 101, chains = 4)
  diagnostic <- mcmc_diagnostics(f)
  beta <- matrix(f$posterior$coefficients[[1]][, , 1], 101, 4)
  row <- diagnostic[diagnostic$family == "coefficients" & diagnostic$parameter == "(Intercept)", ]
  expect_equal(row$rhat, posterior::rhat(beta))
  expect_equal(row$ess_bulk, posterior::ess_bulk(beta))
  f$posterior$coefficients[[1]][, , 2] <- 0
  diagnostic <- mcmc_diagnostics(f)
  selection <- diagnostic[diagnostic$family == "selection", ]
  expect_identical(selection$status, "constant")
  expect_true(is.na(selection$rhat))
  f$posterior$coefficients[[1]][, 1, 2] <- 1
  diagnostic <- mcmc_diagnostics(f)
  selection <- diagnostic[diagnostic$family == "selection", ]
  expect_identical(selection$status, "constant_in_chain")
  expect_warning(IntegMultiReg:::.imr_warn_diagnostics(diagnostic), class = "imr_mcmc_warning")
})

test_that("feature ranking and descriptive comparisons have explicit meanings", {
  f <- small_fit(draws = 1000, chains = 1)
  f$posterior$coefficients[[1]][, , 2] <- 0
  f$posterior$coefficients[[1]][1:29, , 2] <- 1
  expect_equal(inclusion_probabilities(f)[[1]][1, 1], .029)
  expect_output(print(f, rank = TRUE, top = 1), "0.029")
  expect_output(print(f, threshold = .5), "0/1 features")
  table <- compare_fit_summaries(first = f, second = f)
  expect_identical(table$selected_features, c(0L, 0L))
  expect_identical(table$chains, c(1L, 1L))
  bad <- f
  bad$control$response_scale <- "log"
  expect_error(compare_fit_summaries(f, bad), "response scale")
})

test_that("all plot types retain the chain dimension", {
  f <- small_fit()
  path <- tempfile(fileext = ".pdf")
  pdf(path)
  on.exit(
    {
      dev.off()
      unlink(path)
    },
    add = TRUE
  )
  for (type in c("selection", "interaction", "trace", "coefficient_trace", "variance_trace", "selection_trace")) {
    expect_invisible(plot(f, type = type))
  }
  expect_invisible(plot_top_features(f))
  expect_invisible(plot_subgroup_sizes(f))
  expect_error(plot(f, type = "interaction_trace"), "no interaction")
})


test_that("empirical quantile boundaries are stable and respect point masses", {
  q <- IntegMultiReg:::.imr_quantile
  expect_identical(q(1:80, (1 - .95) / 2), q(1:80, .025))
  expect_identical(q(1:80, .025), 2L)
  expect_identical(q(1:80, .975), 78L)
  expect_equal(q(1:80, c(.025, .5, .975), rep(1 / 80, 80)), c(2, 40, 78))
  expect_identical(q(c(rep(0, 79), 1), c(.025, .5, .975)), c(0, 0, 0))
})


test_that("infinite R-hat remains a reported diagnostic failure", {
  fit <- small_fit()
  fit$posterior$coefficients[[1]][, , 1] <- c(rep(1, 20), rep(2, 20), rep(1, 20), rep(4, 20))
  d <- mcmc_diagnostics(fit)
  target <- d$family == "coefficients" & d$parameter == "(Intercept)"
  expect_identical(d$rhat[target], Inf)
  expect_warning(IntegMultiReg:::.imr_warn_diagnostics(d), "R-hat above 1.01")
  expect_output(IntegMultiReg:::.imr_print_diagnostics(d), "Inf")
})

test_that("parallel diagnostics preserve constants, numerical values, RNG and sockets", {
  # Many excluded variables exercise large-block dispatch without requiring a
  # long inferential run. Constant indicators must still remain undefined.
  x <- data.frame(id = 1:20, matrix(0, 20, 200))
  fit <- imr(list(assay = x), data.frame(id = x$id, y = cos(x$id)),
    outcome_type = "continuous", min_subgroup_size = 0,
    priors = imr_priors(nu = -100, forced_scale = 1),
    mcmc = imr_mcmc(
      draws = 1600, burnin = 20, chains = 2,
      seed = 91, diagnostics = FALSE
    )
  )
  sockets <- function() sum(showConnections(all = TRUE)[, "class"] == "sockconn")
  before <- sockets()
  set.seed(199)
  rng <- .Random.seed
  serial <- mcmc_diagnostics(fit)
  fit$control$mcmc$workers <- 2L
  parallel <- mcmc_diagnostics(fit)
  expect_identical(parallel, serial)
  expect_true(all(parallel$status[parallel$family == "selection"] == "constant"))
  expect_identical(.Random.seed, rng)
  expect_identical(sockets(), before)
})
