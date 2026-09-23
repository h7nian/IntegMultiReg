test_that("concordance scores ties only among comparable pairs", {
  score <- function(p, t, e) IntegMultiReg:::.imr_cv_accuracy(
    "right.censored", p, data.frame(id = seq_along(p), time = t, status = e))
  expect_equal(score(c(0, 0), c(1, 2), c(1, 1)), .5)
  expect_equal(score(c(0, 0, 0), c(1, 2, 3), c(0, 1, 1)), .5)
  expect_equal(score(c(1, 2), c(1, 2), c(1, 1)), 1)
  expect_equal(score(c(2, 1), c(1, 2), c(1, 1)), 0)
  expect_equal(score(c(1, 2), c(1, 1), c(1, 0)), 1)
  expect_true(is.na(score(c(1, 2), c(1, 2), c(0, 0))))
  if (requireNamespace("survival", quietly = TRUE)) {
    set.seed(321)
    for (i in 1:20) {
      t <- sample(1:10, 30, replace = TRUE)
      e <- sample(0:1, 30, replace = TRUE)
      p <- sample(1:8, 30, replace = TRUE)
      expected <- survival::concordance(survival::Surv(t, e) ~ p)$concordance
      expect_equal(score(p, t, e), unname(expected))
    }
  }
})

test_that("CV refits without using held-out responses or full-fit posterior draws", {
  x <- data.frame(id = 1:16, marker = sin(1:16))
  y <- data.frame(id = 1:16, y = 1 + x$marker + cos(1:16) / 5)
  fit <- imr(list(assay = x), y, outcome_type = "continuous", min_subgroup_size = 10,
             forced_prior_scale = 1, draws = 20, burnin = 10, seed = 72)
  set.seed(123)
  rng <- .Random.seed
  base <- cv_imr(cv_method = "refit", fit, k = 2, rounds = 1)
  expect_identical(.Random.seed, rng)
  records <- base$predictions
  expect_setequal(records$id, x$id)
  expect_true(all(table(records$fold) == 8L))
  changed <- fit
  heldout <- records$id[records$fold == 1]
  changed$preprocessing$input_data$outcome$y[changed$preprocessing$input_data$outcome$id %in% heldout] <- 1000
  changed$posterior$latent_response_mean <- lapply(changed$posterior$latent_response_mean, function(y) y + 100)
  changed$posterior$selection_draws <- lapply(changed$posterior$selection_draws, function(s) lapply(s, function(g) 1 - g))
  altered <- cv_imr(cv_method = "refit", changed, k = 2, rounds = 1)$predictions
  expect_identical(records$fold, altered$fold)
  expect_equal(records$prediction[records$fold == 1], altered$prediction[altered$fold == 1])
  expect_identical(base$validation, "refit")
  train_ids <- records$id[records$fold != 1]
  changed$preprocessing$input_data$platforms[[1]]$marker[changed$preprocessing$input_data$platforms[[1]]$id %in% heldout] <- 1e6
  refit1 <- IntegMultiReg:::.imr_cv_refit(fit, train_ids, 42)
  refit2 <- IntegMultiReg:::.imr_cv_refit(changed, train_ids, 42)
  expect_equal(refit1$preprocessing[c("features", "response", "covariates")],
               refit2$preprocessing[c("features", "response", "covariates")])
  expect_equal(refit1$posterior$selection_draws, refit2$posterior$selection_draws)

  # The metric is calculated from the returned out-of-fold predictions.
  expect_equal(unname(base$pooled[1, "all"]),
               mean((records$prediction - y$y[match(records$id, y$id)])^2))
})

test_that("CV balances small strata and preserves formula preprocessing", {
  folds <- IntegMultiReg:::.imr_cv_folds(c(0, 0, 1, 1, 1), 5)
  expect_equal(sort(folds), 1:5)
  x <- data.frame(subject = 1:20, marker = sin(1:20))
  d <- data.frame(subject = 1:20, age = 1:20, y = cos(1:20))
  fit <- imr(y ~ poly(age, 2), data = d, platforms = list(assay = x),
    id = "subject", outcome_type = "continuous", min_subgroup_size = 1,
    forced_prior_scale = 1, draws = 10, burnin = 5, seed = 73)
  refit <- IntegMultiReg:::.imr_cv_refit(fit, 1:12, 10)
  expect_equal(refit$preprocessing$formula_data$subject, 1:12)
  expect_equal(unname(as.matrix(refit$preprocessing$input_data$covariates[, -1])),
               unname(poly(1:12, 2)), ignore_attr = TRUE)
  expect_true(all(is.finite(cv_imr(cv_method = "refit", fit, k = 2, rounds = 1)$pooled)))
  old <- fit
  old$preprocessing$formula_data <- NULL
  expect_error(cv_imr(cv_method = "refit", old, k = 2), "raw formula data|missing `formula_data`")
})

test_that("theta uncertainty labels follow the native pair ordering", {
  x <- data.frame(id = 1:32, x = sin(1:32))
  fit <- imr(list(A = x, B = x[17:32, ], C = x[c(9:16, 25:32), ]),
    data.frame(id = 1:32, y = cos(1:32)), outcome_type = "continuous",
    min_subgroup_size = 1, forced_prior_scale = 1, nu = rep(-10, 3), draws = 4, burnin = 2, seed = 1)
  for (j in 1:6) fit$posterior$interaction_draws[[1]][, j] <- j
  result <- posterior_summary(fit)$theta$A
  expect_equal(paste(result$subgroup1, result$subgroup2),
    c("001 011", "001 101", "011 101", "001 111", "011 111", "101 111"))
  expect_equal(result$mean, 1:6)
  expect_equal(result$lower, 1:6)
})

test_that("binary BMA averages probabilities within each selection model", {
  set.seed(71)
  x <- data.frame(id = 1:30, marker = rnorm(30))
  fit <- imr(list(assay = x), data.frame(id = x$id, y = rep(0:1, 15)),
    outcome_type = "binary", forced_prior_scale = 1, nu = 0,
    draws = 20, burnin = 10, min_subgroup_size = 1, seed = 71)
  fit$posterior$selection_draws <- list(list(matrix(0L, 1, 1)), list(matrix(1L, 1, 1)))
  fit$control$mcmc$draws <- 2L
  fit$control$mcmc$burnin <- 0L
  fit$posterior$log_posterior <- fit$posterior$log_posterior[1:2]
  fit$posterior$interaction_draws[[1]] <- fit$posterior$interaction_draws[[1]][1:2, , drop = FALSE]
  new <- list(data.frame(id = 31:32, marker = c(-3, 3)))
  latent <- function(f) {
    f$control$outcome_type <- "continuous"
    predict(f, new)[[1]]$prediction
  }
  modes <- lapply(1:2, function(i) {
    f <- fit
    f$posterior$selection_draws <- f$posterior$selection_draws[i]
    f$control$mcmc$draws <- 1L
    f$posterior$log_posterior <- f$posterior$log_posterior[1L]
    f$posterior$interaction_draws[[1]] <- f$posterior$interaction_draws[[1]][1L, , drop = FALSE]
    latent(f)
  })
  eta <- latent(fit)
  w <- (eta - modes[[1]]) / (modes[[2]] - modes[[1]])
  expected <- (1 - w) * pnorm(modes[[1]]) + w * pnorm(modes[[2]])
  expect_equal(predict(fit, new)[[1]]$prediction, expected, tolerance = 1e-10)
  expect_gt(max(abs(expected - pnorm(eta))), 1e-6)
})
