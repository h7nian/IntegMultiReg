# Independent R linear algebra for the documented empirical importance average.
# This conditions on the fitted posterior and is NOT a Table 1 replication.
# Unlike the native implementation, normalization below uses a softmax directly.
test_that("importance matches a multi-platform R reference for every outcome", {
  outcomes <- list(binary = simIMR$outcome.binary,
                   continuous = simIMR$outcome.continuous,
                   right.censored = simIMR$outcome.survival)
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      fit <- imr(simIMR$platforms, outcomes[[outcome_type]],
                 covariates = simIMR$covariates, outcome_type = outcome_type,
                 method = method, nu = c(-4, -3, -4),
                 residual_prior = c(shape = 0.37, rate = 0.21),
                 draws = 8, burnin = 4, min_subgroup_size = 30, seed = 47)
      result <- cv_imr(fit, k = 2, rounds = 1, cv_method = "importance")
      for (subgroup in seq_along(fit$model$subgroup_names)) {
        name <- fit$model$subgroup_names[subgroup]
        records <- result$predictions[result$predictions$subgroup == name, ]
        y <- fit$posterior$latent_response_mean[[subgroup]]
        alpha <- fit$control$priors$residual[["shape"]]
        rate <- fit$control$priors$residual[["rate"]]
        for (fold in 1:2) {
          train <- which(records$fold != fold)
          test <- which(records$fold == fold)
          predictions <- matrix(NA_real_, length(test), fit$control$mcmc$draws)
          log_weights <- numeric(ncol(predictions))
          for (draw in seq_along(log_weights)) {
            design <- cbind(1, fit$preprocessing$covariates[[subgroup]])
            for (platform in fit$model$subgroup_platforms[[subgroup]]) {
              row <- match(subgroup, fit$model$platform_subgroups[[platform]])
              selected <- as.logical(fit$posterior$selection_draws[[draw]][[platform]][row, ])
              design <- cbind(design,
                fit$preprocessing$features[[subgroup]][[platform]][, selected, drop = FALSE])
            }
            beta <- solve(crossprod(design[train, , drop = FALSE]) +
                            diag(0.001, ncol(design)),
                          crossprod(design[train, , drop = FALSE], y[train]))
            eta <- drop(design[test, , drop = FALSE] %*% beta)
            train_sse <- sum((y[train] - design[train, , drop = FALSE] %*% beta)^2)
            test_sse <- sum((y[test] - eta)^2)
            df <- 2 * alpha + length(train)
            scale <- (rate + train_sse) / df
            log_weights[draw] <- length(test) / 2 * log(scale) +
              (2 * alpha + length(y)) / 2 * log1p(test_sse / (df * scale))
            predictions[, draw] <- if (outcome_type == "binary") pnorm(eta) else eta
          }
          weights <- exp(log_weights - max(log_weights))
          expected <- drop(predictions %*% (weights / sum(weights)))
          expect_equal(records$prediction[test], expected, tolerance = 1e-9,
                       info = paste(outcome_type, method, name, fold))
        }
      }
    }
  }
})
