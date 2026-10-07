# Short API fixtures may legitimately trigger the convergence diagnostic.
# Accept that specific warning while failing on any unrelated warning.
small_posterior <- function(...) withCallingHandlers(sample_regression_posterior(...),
  warning = function(w) {
    expect_match(conditionMessage(w), "Conditional split R-hat")
    invokeRestart("muffleWarning")
  })

posterior_fixture <- function(type = "continuous") {
  x <- data.frame(id = 1:16, x = seq(-1, 1, length.out = 16))
  y <- .5 + x$x + sin(seq_len(16)) / 3
  outcome <- if (type == "binary") data.frame(id = x$id, y = as.integer(y > .5)) else
    if (type == "right.censored") data.frame(id = x$id, time = exp(y), status = rep(c(0, 1), 8)) else
      data.frame(id = x$id, y = y)
  imr(list(assay = x), outcome, outcome_type = type, min_subgroup_size = 2,
      draws = 16, burnin = 8, forced_prior_scale = 1, seed = 12)
}

test_that("coefficient draws preserve zero selection mass and caller RNG", {
  f <- posterior_fixture()
  f$posterior$selection_draws <- lapply(seq_along(f$posterior$selection_draws), function(i) list(matrix(i %% 2L, 1, 1)))
  set.seed(926)
  before <- .Random.seed
  d <- small_posterior(f, output_draws = 100, burnin = 50,
                       min_draws_per_model_chain = 50, seed = 13)
  expect_identical(.Random.seed, before)
  expect_s3_class(d, "imr_posterior")
  expect_equal(dim(d$beta[[1]]), c(100, 2))
  expect_true(all(d$variance[[1]] > 0))
  inactive <- d$selection_draw_index %% 2L == 0L
  expect_true(all(d$beta[[1]][inactive, 2] == 0))
  expect_true(all(d$beta[[1]][!inactive, 2] != 0))
  again <- small_posterior(f, output_draws = 100, burnin = 50,
                           min_draws_per_model_chain = 50, seed = 13)
  expect_identical(d, again)
  expect_named(summary(d)[[1]], c("term", "mean", "sd", "lower", "median", "upper", "probability_nonzero"))
  expect_equal(confint(d), summary(d))
  expect_equal(coef(d)[[1]], colMeans(d$beta[[1]]))
  expect_output(print(d), "Regression posterior")
  expect_error(sample_regression_posterior(f, chains = 1), "chains")
  expect_error(sample_regression_posterior(f, output_draws = NA), "draws")
  expect_error(summary(d, level = 1), "level")
  expect_error(confint(d, parm = "bad"), "parm")
})

test_that("posterior prediction includes residual variance and preserves routing", {
  f <- posterior_fixture()
  d <- small_posterior(f, output_draws = 300, burnin = 100,
                       min_draws_per_model_chain = 150, seed = 31)
  new <- list(assay = data.frame(id = c(19, 18), x = c(.1, -.5)))
  mean <- predict(d, new, quantity = "conditional_mean", seed = 1)
  response <- predict(d, new, quantity = "new_observation", seed = 1)
  expect_identical(mean[[1]]$id, c(19, 18))
  expect_true(all(mean[[1]]$lower < mean[[1]]$upper))
  expect_true(all(response[[1]]$upper - response[[1]]$lower > mean[[1]]$upper - mean[[1]]$lower))
  reversed <- predict(d, list(assay = new$assay[2:1, ]))
  expect_equal(reversed[[1]]$prediction, rev(mean[[1]]$prediction))
  expect_warning(empty <- predict(d, list(assay = new$assay[FALSE, ])), "No subjects")
  expect_equal(nrow(empty[[1]]), 0L)
  expect_named(empty[[1]], c("id", "prediction", "lower", "upper"))
})

test_that("binary and survival posterior predictions use their response scales", {
  for (type in c("binary", "right.censored")) {
    f <- posterior_fixture(type)
    d <- small_posterior(f, output_draws = 80, burnin = 80,
                         min_draws_per_model_chain = 40, seed = 18)
    new <- f$preprocessing$input_data$platforms
    p <- predict(d, new)
    expect_true(all(is.finite(p[[1]]$prediction)))
    expect_true(all(p[[1]]$prediction > 0))
    if (type == "binary") {
      expect_true(all(p[[1]]$upper <= 1))
      r <- predict(d, new, quantity = "new_observation")
      expect_true(all(r[[1]]$lower %in% c(0, 1)))
    } else {
      f$control$response_scale <- NA_character_
      expect_error(sample_regression_posterior(f), "Refit")
    }
  }
})

test_that("conditional kernels have independent mathematical oracles", {
  set.seed(923)
  x <- replicate(5000, IntegMultiReg:::.imr_pmom_normal(0, 1))
  expect_lt(abs(mean(x^2) - 3), .2)
  expect_lt(abs(mean(x > 0) - .5), .04)
  y <- c(-.3, .4, .9, 1.4, -.2, .8)
  A <- length(y) + .5
  mu <- sum(y) / A
  shape <- 3 + length(y) / 2
  rate <- 2 + (sum(y^2) - A * mu^2) / 2
  s <- IntegMultiReg:::.imr_conditional_chain(matrix(1, length(y), 1), y, 2, 0, 3, 2,
                              draws = 4000, burnin = 200)
  expected <- mu + qt(c(.05, .95), 2 * shape) * sqrt(rate / shape / A)
  expect_lt(max(abs(quantile(s[, 1], c(.05, .95)) - expected)), .09)
  tail <- IntegMultiReg:::.imr_lower_normal(rep(-40, 200), 1, 0)
  expect_true(all(is.finite(tail) & tail > 0))
})

test_that("posterior predictions reuse formula encoding across availability groups", {
  ids <- 1:36
  clinical <- data.frame(id = ids, y = 1 + sin(ids),
                         age = ids / 10, group = factor(ids %% 2))
  platforms <- list(a = data.frame(id = 1:24, a = cos(1:24)),
                    b = data.frame(id = 13:36, b = sin(13:36)))
  f <- imr(y ~ age + group, data = clinical, platforms = platforms,
           outcome_type = "continuous", forced_prior_scale = 1, min_subgroup_size = 2,
           draws = 20, burnin = 10, seed = 2)
  original <- f
  d <- small_posterior(f, output_draws = 20, burnin = 50, min_draws_per_model_chain = 20)
  expect_identical(f, original)
  expect_length(d$beta, 3L)
  expect_true(all(vapply(d$beta, function(b) all(b[, 1:3] != 0), TRUE)))
  expect_true(all(vapply(d$beta, function(b) identical(colnames(b)[1:3],
    c("(Intercept)", "clinical:age", "clinical:group1")), TRUE)))
  encoded <- data.frame(id = ids, age = clinical$age, group1 = ids %% 2)
  p <- predict(d, platforms, covariates = clinical)
  expect_equal(p, predict(d, platforms, covariates = encoded))
  expect_equal(sum(vapply(p, nrow, 1L)), length(ids))
  # Identical joint model rows must route the same retained mask in every group.
  for (g in seq_along(d$beta)) {
    offset <- 3L
    for (platform in f$model$subgroup_platforms[[g]]) {
      active <- vapply(d$selection_draw_index, function(s) f$posterior$selection_draws[[s]][[platform]][
        match(g, f$model$platform_subgroups[[platform]]), 1], 0)
      expect_identical(d$beta[[g]][, offset + 1L] != 0, active == 1)
      offset <- offset + 1L
    }
  }
  bad <- clinical
  bad$group <- "unseen"
  expect_error(predict(d, platforms, covariates = bad), "new level")
})

test_that("survival point summaries are medians and prediction restores RNG", {
  f <- posterior_fixture("right.censored")
  d <- small_posterior(f, output_draws = 30, burnin = 40, min_draws_per_model_chain = 20)
  # Independent oracle at a training row: its standardized marker is known.
  row <- f$preprocessing$input_data$platforms[[1]][1, , drop = FALSE]
  x <- (row$x - mean(f$preprocessing$input_data$platforms[[1]]$x)) /
    sd(f$preprocessing$input_data$platforms[[1]]$x)
  value <- exp(d$beta[[1]][, 1] + d$beta[[1]][, 2] * x + d$variance[[1]] / 2)
  set.seed(393)
  rng <- .Random.seed
  result <- predict(d, list(assay = row))[[1]]
  expect_identical(.Random.seed, rng)
  expect_equal(result$prediction, median(value))
  expect_equal(c(result$lower, result$upper), unname(quantile(value, c(.025, .975))))
})

test_that("latent draws are opt-in and leave the coefficient draws untouched", {
  f <- posterior_fixture("right.censored")
  without <- small_posterior(f, output_draws = 12, burnin = 6, chains = 2, seed = 3)
  with <- small_posterior(f, output_draws = 12, burnin = 6, chains = 2, seed = 3,
                          latent = TRUE)
  expect_null(without$latent)
  expect_identical(without$beta, with$beta)
  expect_identical(without$variance, with$variance)
  expect_equal(dim(with$latent[[1L]]),
               c(12L, nrow(f$preprocessing$response[[1L]])))
  expect_false(anyNA(with$latent[[1L]]))
})

test_that("the augmented response moves only where the outcome is censored", {
  f <- posterior_fixture("right.censored")
  p <- small_posterior(f, output_draws = 12, burnin = 6, chains = 2, seed = 4,
                       latent = TRUE)
  for (g in seq_along(p$latent)) {
    response <- f$preprocessing$response[[g]]
    status <- response[, 2L]
    moves <- apply(p$latent[[g]], 2L, stats::sd) > 0
    expect_identical(unname(moves), unname(status == 0))
    # an event time is observed, so its augmented value is that time
    expect_equal(unname(p$latent[[g]][1L, status == 1]),
                 unname(response[status == 1, 1L]))
    # a censored time is only a lower bound, so every draw sits above it
    censored <- which(status == 0)
    if (length(censored)) {
      bound <- matrix(response[censored, 1L], nrow(p$latent[[g]]),
                      length(censored), byrow = TRUE)
      expect_true(all(p$latent[[g]][, censored, drop = FALSE] >= bound))
    }
  }
})

test_that("latent summaries are per subject and refuse the cases without them", {
  f <- posterior_fixture("binary")
  p <- small_posterior(f, output_draws = 12, burnin = 6, chains = 2, seed = 5,
                       latent = TRUE)
  s <- summary(p, parm = "latent")
  expect_named(s[[1L]], c("id", "mean", "sd", "lower", "median", "upper"))
  expect_identical(nrow(s[[1L]]), nrow(f$preprocessing$response[[1L]]))
  expect_identical(s, confint(p, parm = "latent"))
  expect_error(summary(small_posterior(f, output_draws = 12, burnin = 6, chains = 2,
                                       seed = 5), parm = "latent"),
               "sample_regression_posterior\\(latent = TRUE\\)")
  continuous <- small_posterior(posterior_fixture(), output_draws = 12, burnin = 6,
                                chains = 2, seed = 5, latent = TRUE)
  expect_null(continuous$latent)
  expect_error(summary(continuous, parm = "latent"), "no latent response")
})

test_that("variance intervals summarize variance draws, not coefficients", {
  d <- small_posterior(posterior_fixture(), output_draws = 20, burnin = 10,
                       min_draws_per_model_chain = 10)
  out <- summary(d, parm = "variance", level = .8)[[1]]
  expect_identical(out$parameter, "residual_variance")
  expect_equal(out$mean, mean(d$variance[[1]]))
  expect_equal(unname(unlist(out[c("lower", "median", "upper")])),
               unname(quantile(d$variance[[1]], c(.1,.5,.9), type = 7)))
  expect_identical(confint(d, parm = "variance", level = .8),
                   summary(d, parm = "variance", level = .8))
  expect_error(predict(d, list(assay = posterior_fixture()$preprocessing$input_data$platforms[[1]]),
                       type = "response"), "Unused argument")
})

test_that("stored regression posterior controls upgrade without resampling", {
  d <- small_posterior(posterior_fixture(), output_draws = 20, burnin = 10,
                       min_draws_per_model_chain = 10)
  old <- d
  old$fit$schema_version <- 2L
  names(old$fit$control)[names(old$fit$control) == "model_variant"] <- "method"
  names(old$fit$control)[names(old$fit$control) == "selection_update"] <- "sampler_method"
  old$fit$control$sampler_method <- "paper"
  old$fit$control$numerical$prior_indexing <- "standard"
  names(old$control)[names(old$control) == "output_draws"] <- "draws"
  names(old$control)[names(old$control) == "min_draws_per_model_chain"] <- "conditional_draws"
  names(old)[names(old) == "selection_draw_index"] <- "model_draw"
  names(old$diagnostics)[names(old$diagnostics) == "selection_model"] <- "model"
  names(old$diagnostics)[names(old$diagnostics) == "output_draws"] <- "returned_draws"
  names(old$diagnostics)[names(old$diagnostics) == "draws_per_model_chain"] <- "conditional_draws"
  expect_error(validate_imr_object(old), "upgrade_imr_object")
  set.seed(872)
  rng <- .Random.seed
  restored <- upgrade_imr_object(old)
  expect_identical(.Random.seed, rng)
  expect_identical(restored$beta, d$beta)
  expect_identical(restored$variance, d$variance)
  expect_identical(restored$selection_draw_index, d$selection_draw_index)
  expect_identical(restored$diagnostics, d$diagnostics)
  expect_true(validate_imr_object(restored))
  broken <- restored; broken$variance[[1]][1] <- -1
  expect_error(validate_imr_object(broken), "variance draws")
})
