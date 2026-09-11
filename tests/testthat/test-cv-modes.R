test_that("legacy is the explicit and implicit default", {
  implicit <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 4)
  explicit <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 4,
                     cv_method = "legacy")
  expect_identical(implicit, explicit)
  expect_identical(implicit$validation, "legacy")
  expect_error(cv_imr(fit_bin, cv_method = "unknown"), "should be one of")
})

test_that("all CV modes return complete reproducible subject records", {
  set.seed(271)
  state <- .Random.seed
  for (mode in c("legacy", "importance", "refit")) {
    first <- cv_imr(fit_bin, k = 3, rounds = 2, max_models = 3, cv_method = mode)
    second <- cv_imr(fit_bin, k = 3, rounds = 2, max_models = 3, cv_method = mode)
    expect_identical(first, second)
    expect_identical(.Random.seed, state)
    expect_identical(first$validation, mode)
    expect_named(first, c("pooled", "fold_mean", "predictions", "metric", "validation"))
    expect_true(all(is.finite(first$predictions$prediction)))
    expect_true(all(first$predictions$prediction >= 0 & first$predictions$prediction <= 1))
    for (round in 1:2) {
      records <- first$predictions[first$predictions$round == round, ]
      expect_equal(nrow(records), sum(fit_bin$model$sample_sizes))
      expect_equal(anyDuplicated(records$id), 0L)
      expect_true(all(records$fold %in% 1:3))
    }
  }
})

test_that("post-fit modes share historical partitions but importance uses all draws", {
  legacy <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 1)
  importance <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 1,
                       cv_method = "importance")
  unlimited <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 100,
                      cv_method = "importance")
  expect_identical(importance, unlimited)
  expect_identical(legacy$predictions[c("id", "fold", "subgroup")],
                   importance$predictions[c("id", "fold", "subgroup")])
})

test_that("importance prediction agrees with a direct empirical weight calculation", {
  ids <- seq_len(24)
  fit <- imr(list(assay = data.frame(id = ids, marker = sin(ids))),
             data.frame(id = ids, y = cos(ids)), outcome_type = "continuous",
             method = "bms", draws = 8, burnin = 2, seed = 13,
             min_subgroup_size = 0)
  cv <- cv_imr(fit, k = 2, rounds = 1, cv_method = "importance")
  y <- fit$posterior$latent_response_mean[[1]]
  x <- fit$preprocessing$features[[1]][[1]]
  alpha <- fit$control$priors$residual[["shape"]]
  psi <- fit$control$priors$residual[["rate"]]
  folds <- cv$predictions$fold[match(ids, cv$predictions$id)]
  for (fold in 1:2) {
    train <- which(folds != fold)
    test <- which(folds == fold)
    predictions <- matrix(0, length(test), length(fit$posterior$selection_draws))
    log_weights <- numeric(ncol(predictions))
    for (draw in seq_along(log_weights)) {
      selected <- as.logical(fit$posterior$selection_draws[[draw]][[1]][1, ])
      design <- cbind(1, x[, selected, drop = FALSE])
      beta <- solve(crossprod(design[train, , drop = FALSE]) +
                      diag(.001, ncol(design)),
                    crossprod(design[train, , drop = FALSE], y[train]))
      predictions[, draw] <- design[test, , drop = FALSE] %*% beta
      df <- 2 * alpha + length(train)
      scale <- (psi + sum((y[train] - design[train, , drop = FALSE] %*% beta)^2)) / df
      log_weights[draw] <- length(test) / 2 * log(scale) +
        (2 * alpha + length(y)) / 2 *
        log1p(sum((y[test] - predictions[, draw])^2) / (df * scale))
    }
    weights <- exp(log_weights - max(log_weights))
    expected <- drop(predictions %*% (weights / sum(weights)))
    observed <- cv$predictions$prediction[match(ids[test], cv$predictions$id)]
    expect_equal(observed, expected, tolerance = 1e-10)
  }
})

test_that("post-fit modes support zero burn-in and one retained draw", {
  ids <- seq_len(12)
  fit <- imr(list(assay = data.frame(id = ids, marker = sin(ids))),
             data.frame(id = ids, y = cos(ids)), outcome_type = "continuous",
             draws = 1, burnin = 0, seed = 11, min_subgroup_size = 0)
  for (mode in c("legacy", "importance")) {
    result <- cv_imr(fit, k = 2, rounds = 1, cv_method = mode)
    expect_true(all(is.finite(result$pooled)))
  }
})
