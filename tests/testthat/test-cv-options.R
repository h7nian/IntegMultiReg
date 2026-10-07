test_that("reweighting exposes independent state and density choices", {
  fit <- fit_demo("continuous", total = 12, burn = 6, seed = 53)
  defaults <- cv_imr(fit, k = 3, rounds = 2, cv_method = "reweight")
  explicit <- cv_imr(fit, k = 3, rounds = 2, cv_method = "reweight", ridge = .001,
                    model_set = "all_draws", df_method = "fractional")
  expect_identical(explicit, defaults)
  expect_null(defaults$control$max_models)
  expect_identical(defaults$control$model_set, "all_draws")
  ranked <- cv_imr(fit, k = 3, rounds = 2, cv_method = "reweight",
                  model_set = "top_unique", max_models = 2, df_method = "integer")
  expect_identical(ranked$control$max_models, 2L)
  expect_identical(ranked$control$df_method, "integer")
  expect_identical(ranked$control$score_rule, "mean_squared_error")
})

test_that("saved post-fit membership and row order replay every prediction", {
  for (type in c("continuous", "binary", "right.censored")) {
    fit <- fit_demo(type, total = 12, burn = 6, seed = 53)
    for (mode in "reweight") {
      original <- cv_imr(fit, k = 3, rounds = 2, cv_method = mode)
      folds <- original$control$folds
      folds <- folds[rev(seq_len(nrow(folds))), ]
      set.seed(71)
      rng <- .Random.seed
      replay <- cv_imr(fit, cv_method = mode, folds = folds)
      expect_identical(.Random.seed, rng)
      expect_identical(replay[1:5], original[1:5], info = paste(type, mode))
      expect_identical(replay$control$k, 3L)
      expect_identical(replay$control$rounds, 2L)
    }
  }
})

test_that("ridge and predictive degrees of freedom match independent linear algebra", {
  id <- seq_len(40)
  x <- data.frame(id, a = sin(id), b = cos(id / 3))
  y <- data.frame(id, y = sin(id / 2) + cos(id / 5))
  fit <- imr(list(assay = x), y, outcome_type = "continuous", model_variant = "bms",
    draws = 6, burnin = 2, seed = 17, min_subgroup_size = 0,
    residual_prior = c(shape = .37, rate = .21))
  empty <- fit$posterior$selection_draws[[1]]
  empty[[1]][] <- 0L
  first <- both <- empty
  first[[1]][1, 1] <- 1L
  both[[1]][] <- 1L
  fit$posterior$selection_draws <- list(empty, first, first, both, empty, first)
  for (ridge in c(0, .001, .3)) for (df_method in c("fractional", "integer")) {
    result <- cv_imr(fit, k = 2, rounds = 1, cv_method = "reweight",
                     ridge = ridge, df_method = df_method)
    records <- result$predictions
    observed <- fit$posterior$latent_response_mean[[1]]
    for (fold in 1:2) {
      test <- which(records$fold == fold)
      train <- which(records$fold != fold)
      eta <- matrix(0, length(test), 6)
      log_weight <- numeric(6)
      for (draw in 1:6) {
        selected <- as.logical(fit$posterior$selection_draws[[draw]][[1]][1, ])
        design <- cbind(1, fit$preprocessing$features[[1]][[1]][, selected, drop = FALSE])
        beta <- solve(crossprod(design[train, , drop = FALSE]) + diag(ridge, ncol(design)),
                      crossprod(design[train, , drop = FALSE], observed[train]))
        eta[, draw] <- drop(design[test, , drop = FALSE] %*% beta)
        train_sse <- sum((observed[train] - design[train, , drop = FALSE] %*% beta)^2)
        test_sse <- sum((observed[test] - eta[, draw])^2)
        df <- 2 * .37 + length(train)
        if (df_method == "integer") df <- trunc(df)
        scale <- (.21 + train_sse) / df
        log_weight[draw] <- length(test) / 2 * log(scale) +
          (2 * .37 + length(id)) / 2 * log1p(test_sse / (df * scale))
      }
      weights <- exp(log_weight - max(log_weight))
      expected <- drop(eta %*% (weights / sum(weights)))
      expect_equal(records$prediction[test], expected, tolerance = 1e-11,
                   info = paste(ridge, df_method, fold))
    }
  }
})

test_that("explicit options work through parallel and bounded-cache paths", {
  fit <- fit_demo("binary", total = 12, burn = 6, seed = 53)
  settings <- IntegMultiReg:::.imr_cv_settings("reweight", ridge = .2,
    model_set = "all_draws", df_method = "integer")
  result <- cv_imr(fit, k = 2, rounds = 1, cv_method = "reweight", ridge = .2, model_set = "all_draws",
                   df_method = "integer")
  replay <- cv_imr(fit, cv_method = "reweight", ridge = .2, model_set = "all_draws", df_method = "integer",
    folds = result$control$folds, workers = 2)
  expect_identical(replay[1:5], result[1:5])
  for (bytes in c(0, 256, 128 * 1024^2)) {
    native <- IntegMultiReg:::.imr_call_cv_postfit_native(fit, 2L, 1L, 12L,
      FALSE, FALSE, settings = settings, cache_bytes = bytes, cache_hash_bits = 0L)
    expect_identical(as.vector(native$predictions), result$predictions$prediction)
  }
  refit <- cv_imr(fit, cv_method = "refit", folds = result$control$folds,
    max_models = 4)
  expect_identical(refit$predictions[c("id", "round", "fold")],
    cv_imr(fit, cv_method = "refit", folds = result$control$folds,
      max_models = 4)$predictions[c("id", "round", "fold")])
})

test_that("invalid CV options and fold tables fail before computation", {
  for (value in list(-1, Inf, NA_real_, numeric(), c(0, 1), "0"))
    expect_error(cv_imr(fit_bin, cv_method = "reweight", ridge = value), "ridge")
  for (name in c("model_set", "df_method")) {
    args <- list(object = fit_bin, cv_method = "reweight")
    args[[name]] <- "unknown"
    expect_error(do.call(cv_imr, args), name)
  }
  for (name in c("ridge", "model_set", "df_method")) {
    args <- list(object = fit_bin, cv_method = "refit")
    args[[name]] <- if (name == "ridge") 0 else "all_draws"
    expect_error(do.call(cv_imr, args), "only to reweighting")
  }
  folds <- cv_imr(fit_bin, k = 2, rounds = 1)$control$folds
  expect_error(cv_imr(fit_bin, folds = folds, k = 3), "agree")
  expect_error(cv_imr(fit_bin, folds = folds, rounds = 2), "agree")
  expect_error(cv_imr(fit_bin, folds = folds[-1, ]), "each modelled")
  expect_error(cv_imr(fit_bin, folds = rbind(folds, folds[1, ])), "each modelled")
  bad <- folds; bad$id[1] <- "unknown"
  expect_error(cv_imr(fit_bin, folds = bad), "each modelled")
  bad <- folds; bad$row_order[] <- 1L
  expect_error(cv_imr(fit_bin, folds = bad), "permutation")
  bad <- folds; bad$fold[1] <- NA_integer_
  expect_error(cv_imr(fit_bin, folds = bad), "positive finite integers")
})
