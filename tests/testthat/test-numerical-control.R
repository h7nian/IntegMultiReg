test_that("numerical conventions are validated and travel with the fit", {
  args <- list(x = simIMR$platforms, outcome = simIMR$outcome.continuous,
    covariates = simIMR$covariates, outcome_type = "continuous",
    draws = 8, burnin = 4, min_subgroup_size = 30, seed = 53)
  baseline <- do.call(imr, args)
  explicit <- do.call(imr, c(args, list(prior_indexing = "standard",
    laplace_max_iter = c(initial = 25, selection = 40, latent = 25, prediction = 40),
    laplace_tolerance = .001)))
  expect_identical(baseline$posterior, explicit$posterior)
  historical <- do.call(imr, c(args, list(prior_indexing = "code2017",
    laplace_max_iter = c(initial = 25, selection = 40, latent = 25, prediction = 25))))
  expect_false(identical(historical$posterior, baseline$posterior))
  expect_true(validate_imr(historical))
  prediction <- predict(historical, simIMR$platforms, covariates = simIMR$covariates)
  values <- unlist(lapply(prediction, `[[`, "prediction"), use.names = FALSE)
  expect_gt(length(values), 0)
  expect_true(all(is.finite(values)))
  ids <- historical$preprocessing$input_data$availability$id
  refit <- IntegMultiReg:::.imr_cv_refit(historical, ids, seed = 53)
  expect_identical(refit$control$numerical, historical$control$numerical)
  for (mode in c("legacy", "importance", "refit")) {
    cv <- cv_imr(historical, k = 2, rounds = 1, cv_method = mode)
    expect_identical(cv$control$numerical, historical$control$numerical)
  }
  old <- baseline; old$control$numerical <- NULL
  expect_identical(predict(old, simIMR$platforms, covariates = simIMR$covariates),
                   predict(baseline, simIMR$platforms, covariates = simIMR$covariates))
  scalar <- IntegMultiReg:::.imr_numerical_control(laplace_max_iter = 20)
  expect_identical(unname(scalar$laplace_max_iter), rep(20L, 4))
  for (bad in list(0, -1, Inf, NA_real_, 1.5, c(25, 40, 25, 40),
                  c(selection = 40, initial = 25, latent = 25, prediction = 40)))
    expect_error(IntegMultiReg:::.imr_numerical_control(laplace_max_iter = bad),
                 "laplace_max_iter")
  for (bad in c(0, -1, Inf, NA_real_))
    expect_error(IntegMultiReg:::.imr_numerical_control(laplace_tolerance = bad),
                 "laplace_tolerance")
})

test_that("continued fold RNG is reproducible, replayable, and parallel-safe", {
  fit <- fit_demo("continuous", total = 8, burn = 4)
  reset <- cv_imr(fit, k = 3, rounds = 2)
  continued <- cv_imr(fit, k = 3, rounds = 2, fold_rng = "continue")
  expect_false(identical(reset$control$folds, continued$control$folds))
  expect_identical(continued, cv_imr(fit, k = 3, rounds = 2, fold_rng = "continue"))
  parallel <- cv_imr(fit, k = 3, rounds = 2, fold_rng = "continue", workers = 2)
  expect_identical(parallel, continued)
  replay <- cv_imr(fit, folds = continued$control$folds)
  expect_identical(replay[1:5], continued[1:5])
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_identical(continued, cv_imr(readRDS(path), k = 3, rounds = 2, fold_rng = "continue"))
  old <- fit; old$control$rng_state <- NULL
  expect_identical(cv_imr(old, k = 3, rounds = 2), reset)
  expect_error(cv_imr(old, fold_rng = "continue"), "saved GSL state")
  expect_error(cv_imr(fit, folds = continued$control$folds, fold_rng = "continue"), "explicit.*folds")
  expect_error(cv_imr(fit, cv_method = "refit", fold_rng = "continue"), "only to post-fit")
  bad <- fit; bad$control$rng_state$state <- as.raw(1)
  expect_error(cv_imr(bad, fold_rng = "continue"), "incompatible")
})

test_that("external preprocessing retains supplied units throughout prediction and refit", {
  id <- seq_len(40)
  x <- data.frame(id, a = 2 + sin(id), b = 5 * cos(id))
  covariates <- data.frame(id, c = seq_len(40) / 20)
  y <- data.frame(id, y = sin(id / 3))
  fit <- imr(list(assay = x), y, covariates = covariates,
    outcome_type = "continuous", standardize = FALSE,
    draws = 8, burnin = 4, min_subgroup_size = 0, seed = 3)
  expect_equal(unname(fit$preprocessing$features[[1]][[1]]), unname(as.matrix(x[-1])))
  expect_equal(unname(fit$preprocessing$covariates[[1]]), unname(as.matrix(covariates[-1])))
  expect_equal(unname(fit$preprocessing$feature_center[[1]][[1]]), c(0, 0))
  expect_equal(unname(fit$preprocessing$feature_scale[[1]][[1]]), c(1, 1))
  whole <- predict(fit, list(assay = x), covariates = covariates)[[1]]
  subset <- predict(fit, list(assay = x[1:4, ]), covariates = covariates[1:4, ])[[1]]
  expect_equal(subset$prediction, whole$prediction[1:4], tolerance = 1e-12)
  refit <- IntegMultiReg:::.imr_cv_refit(fit, id[1:30], seed = 3)
  expect_false(refit$control$standardize)
  expect_equal(unname(refit$preprocessing$features[[1]][[1]]),
               unname(as.matrix(x[1:30, -1])))
  expect_false(cv_imr(fit, k = 2, rounds = 1, cv_method = "refit")$control$standardize)
  expect_error(imr(list(assay = x), y, standardize = NA), "standardize")
})
