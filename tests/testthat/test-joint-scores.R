test_that("concordance scores ties only among comparable pairs", {
  score <- function(p, t, e) {
    IntegMultiReg:::.imr_cv_accuracy(
      "right.censored", p, data.frame(id = seq_along(p), time = t, status = e)
    )
  }
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


test_that("binary scores equal pairwise tie-aware AUC and undefined scores stay NA", {
  accuracy <- IntegMultiReg:::.imr_cv_accuracy
  y <- c(0, 0, 1, 1, 1)
  for (p in list(c(.1, .2, .3, .4, .5), c(.5, .5, .2, .7, .5), rep(.5, 5))) {
    difference <- outer(p[y == 1], p[y == 0], "-")
    expected <- mean((difference > 0) + .5 * (difference == 0))
    expect_equal(accuracy("binary", p, data.frame(id = 1:5, y = y)), expected)
  }
  expect_true(is.na(accuracy("binary", c(.2, .8), data.frame(id = 1:2, y = c(1, 1)))))
  expect_true(is.na(accuracy(
    "right.censored", 1:2,
    data.frame(id = 1:2, time = 1:2, status = c(0, 0))
  )))
  expect_true(is.na(accuracy("continuous", numeric(), data.frame(id = integer(), y = numeric()))))
})

test_that("training-fold formula bases use only training covariates", {
  x <- data.frame(subject = 1:20, marker = sin(1:20))
  d <- data.frame(subject = 1:20, age = 1:20, y = cos(1:20))
  fit <- imr(y ~ poly(age, 2),
    data = d, platforms = list(assay = x), id = "subject",
    outcome_type = "continuous", min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
    mcmc = imr_mcmc(draws = 10, burnin = 5, chains = 2, seed = 73, diagnostics = FALSE)
  )
  refit <- IntegMultiReg:::.imr_cv_refit(fit, 1:12, 10)
  expect_equal(refit$preprocessing$formula_data$subject, 1:12)
  expect_equal(unname(as.matrix(refit$preprocessing$input_data$covariates[, -1])),
    unname(poly(1:12, 2)),
    ignore_attr = TRUE
  )
  changed <- fit
  changed$preprocessing$formula_data$age[13:20] <- 1e6
  changed$preprocessing$input_data$platforms[[1]]$marker[13:20] <- -1e6
  again <- IntegMultiReg:::.imr_cv_refit(changed, 1:12, 10)
  expect_identical(refit$posterior, again$posterior)
  expect_identical(refit$preprocessing$features, again$preprocessing$features)
})

test_that("interaction summaries retain native pair labels with four subgroups", {
  x <- data.frame(id = 1:32, x = sin(1:32))
  fit <- imr(list(A = x, B = x[17:32, ], C = x[c(9:16, 25:32), ]),
    data.frame(id = 1:32, y = cos(1:32)),
    outcome_type = "continuous", min_subgroup_size = 1,
    priors = imr_priors(forced_scale = 1, nu = rep(-10, 3)),
    mcmc = imr_mcmc(draws = 4, burnin = 2, chains = 2, seed = 1, diagnostics = FALSE)
  )
  expected <- c("001:011", "001:101", "011:101", "001:111", "011:111", "101:111")
  expect_identical(dimnames(fit$posterior$interaction$A)[[3L]], expected)
  for (j in 1:6) fit$posterior$interaction$A[, , j] <- j
  table <- subset(confint(fit, parm = "interaction"), group == "A")
  expect_identical(table$parameter, expected)
  expect_equal(table$mean, 1:6)
  expect_equal(table$lower, 1:6)
})
