test_that("refit CV records real partitions, seeds and independently reconstructed scores", {
  f <- small_fit(draws = 20, burnin = 10)
  set.seed(176)
  rng <- .Random.seed
  cv <- cv_imr(f, k = 2, rounds = 2, seed = 56)
  expect_s3_class(cv, "imr_cv")
  expect_identical(.Random.seed, rng)
  expect_identical(dim(cv$control$refit_chain_seeds), c(2L, 2L, 2L))
  expect_identical(dim(cv$control$refit_initial_seeds), c(2L, 2L, 2L))
  expect_identical(cv$control$mcmc$workers, 1L)
  expect_false(cv$control$mcmc$keep_latent)
  y <- f$preprocessing$input_data$outcome
  for (round in 1:2) {
    p <- cv$predictions[cv$predictions$round == round, ]
    error <- (p$prediction - y$y[match(p$id, y$id)])^2
    expect_equal(unname(cv$pooled[round, "all"]), mean(error))
    expect_equal(unname(cv$fold_mean[round, "all"]), mean(tapply(error, p$fold, mean)))
  }
  replay <- cv_imr(f, folds = cv$control$folds, seed = 56)
  expect_identical(cv$predictions, replay$predictions)
  expect_identical(cv$control$refit_chain_seeds, replay$control$refit_chain_seeds)
  expect_output(print(cv), "training-fold joint MCMC")
})

test_that("held-out values cannot alter their training-fold fit or seed", {
  d <- small_data("binary")
  fit <- imr(d$platforms, d$outcome,
    outcome_type = "binary", min_subgroup_size = 0,
    priors = imr_priors(forced_scale = 1), mcmc = imr_mcmc(draws = 20, burnin = 10, chains = 2, seed = 10, diagnostics = FALSE)
  )
  base <- cv_imr(fit, k = 2, rounds = 1, seed = 71)
  heldout <- base$control$folds$id[base$control$folds$fold == 1]
  altered <- fit
  idx <- match(heldout, altered$preprocessing$input_data$outcome$id)
  altered$preprocessing$input_data$outcome$y[idx] <- 1 - altered$preprocessing$input_data$outcome$y[idx]
  check <- cv_imr(altered, folds = base$control$folds, seed = 71)
  expect_identical(base$control$refit_seeds, check$control$refit_seeds)
  expect_identical(
    base$predictions$prediction[base$predictions$fold == 1],
    check$predictions$prediction[check$predictions$fold == 1]
  )
})

test_that("observed likelihood ratios integrate binary utilities and censored times", {
  for (type in c("continuous", "binary", "right.censored")) {
    f <- small_fit(type)
    g <- 1L
    id <- f$preprocessing$subject_ids[[g]][1:3]
    X <- IntegMultiReg:::.imr_joint_design(f$model, f$preprocessing, g)
    b <- matrix(f$posterior$coefficients[[g]], 80, 2)
    v <- as.vector(f$posterior$variance)
    y <- f$preprocessing$response[[g]]
    expected <- numeric(80)
    for (i in 1:3) {
      mu <- drop(b %*% X[i, ])
      ll <- if (type == "binary") {
        pnorm(if (y[i, 1] == 1) mu / sqrt(v) else -mu / sqrt(v), log.p = TRUE)
      } else if (type == "right.censored" && y[i, 2] == 0) {
        pnorm((y[i, 1] - mu) / sqrt(v), lower.tail = FALSE, log.p = TRUE)
      } else {
        dnorm(y[i, 1], mu, sqrt(v), log = TRUE) - if (type == "right.censored") y[i, 1] else 0
      }
      expected <- expected + ll
    }
    expect_equal(IntegMultiReg:::.imr_holdout_log_likelihood(f, id), expected)
  }
})

test_that("reweighting exposes unstable weights instead of silently treating them as refits", {
  f <- small_fit(draws = 100, burnin = 100)
  warnings <- character()
  cv <- withCallingHandlers(cv_imr(f, k = 2, rounds = 1, cv_method = "reweight"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(cv$control$preprocessing, "full_fit")
  expect_true(all(c("pareto_k", "importance_ess", "reliable", "warning") %in% names(cv$diagnostics)))
  expect_true(all(is.finite(cv$pooled)))
  if (any(!cv$diagnostics$reliable)) expect_true(any(grepl("unreliable", warnings)))
  expect_error(cv_imr(f, ridge = .001), "unused argument")
  expect_error(cv_imr(f, cv_method = "legacy"), "arg")
})

test_that("parallel chains and folds preserve numerical results and close sockets", {
  sockets <- function() sum(showConnections(all = TRUE)[, "class"] == "sockconn")
  before <- sockets()
  serial <- small_fit(draws = 20, burnin = 10, workers = 1)
  parallel <- small_fit(draws = 20, burnin = 10, workers = 2)
  expect_identical(serial$posterior, parallel$posterior)
  expect_identical(serial$control$chain_seeds, parallel$control$chain_seeds)
  a <- cv_imr(serial, k = 2, rounds = 1, workers = 1)
  b <- cv_imr(serial, k = 2, rounds = 1, workers = 2)
  expect_identical(a, b)
  expect_identical(sockets(), before)
  expect_error(IntegMultiReg:::.imr_map_tasks(list(1, 2), function(task) stop("worker failure"), 2), "worker failure")
  expect_identical(sockets(), before)
})
