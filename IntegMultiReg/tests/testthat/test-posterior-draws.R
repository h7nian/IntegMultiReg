# Short API fixtures may legitimately trigger the convergence diagnostic.
# Accept that specific warning while failing on any unrelated warning.
small_posterior <- function(...) withCallingHandlers(posterior_draws(...),
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
  d <- small_posterior(f, draws = 100, burnin = 50,
                       conditional_draws = 50, seed = 13)
  expect_identical(.Random.seed, before)
  expect_s3_class(d, "imr_posterior")
  expect_equal(dim(d$beta[[1]]), c(100, 2))
  expect_true(all(d$variance[[1]] > 0))
  inactive <- d$model_draw %% 2L == 0L
  expect_true(all(d$beta[[1]][inactive, 2] == 0))
  expect_true(all(d$beta[[1]][!inactive, 2] != 0))
  again <- small_posterior(f, draws = 100, burnin = 50,
                           conditional_draws = 50, seed = 13)
  expect_identical(d, again)
  expect_named(summary(d)[[1]], c("term", "mean", "sd", "lower", "median", "upper", "probability_nonzero"))
  expect_equal(confint(d), summary(d))
  expect_equal(coef(d)[[1]], colMeans(d$beta[[1]]))
  expect_output(print(d), "coefficient posterior")
  expect_error(posterior_draws(f, chains = 1), "chains")
  expect_error(posterior_draws(f, draws = NA), "draws")
  expect_error(summary(d, level = 1), "level")
  expect_error(confint(d, parm = "bad"), "parm")
})

test_that("posterior prediction includes residual variance and preserves routing", {
  f <- posterior_fixture()
  d <- small_posterior(f, draws = 300, burnin = 100,
                       conditional_draws = 150, seed = 31)
  new <- list(assay = data.frame(id = c(19, 18), x = c(.1, -.5)))
  mean <- predict(d, new, type = "mean", seed = 1)
  response <- predict(d, new, type = "response", seed = 1)
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
    d <- small_posterior(f, draws = 80, burnin = 80,
                         conditional_draws = 40, seed = 18)
    new <- f$preprocessing$input_data$platforms
    p <- predict(d, new)
    expect_true(all(is.finite(p[[1]]$prediction)))
    expect_true(all(p[[1]]$prediction > 0))
    if (type == "binary") {
      expect_true(all(p[[1]]$upper <= 1))
      r <- predict(d, new, type = "response")
      expect_true(all(r[[1]]$lower %in% c(0, 1)))
    } else {
      f$control$response_scale <- NA_character_
      expect_error(posterior_draws(f), "Refit")
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
  d <- small_posterior(f, draws = 20, burnin = 50, conditional_draws = 20)
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
      active <- vapply(d$model_draw, function(s) f$posterior$selection_draws[[s]][[platform]][
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
  d <- small_posterior(f, draws = 30, burnin = 40, conditional_draws = 20)
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
