# Short chains exercise all routes; these are not convergence/performance claims.
test_that("every CV mode works for all outcome and model combinations", {
  outcomes <- list(binary = simIMR$outcome.binary,
                   continuous = simIMR$outcome.continuous,
                   right.censored = simIMR$outcome.survival)
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      fit <- imr(simIMR$platforms, outcomes[[outcome_type]],
                 covariates = simIMR$covariates, outcome_type = outcome_type,
                 method = method, nu = c(-4, -3, -4),
                 draws = 12, burnin = 6, min_subgroup_size = 30, seed = 53)
      for (cv_method in c("legacy", "refit", "importance")) {
        context <- paste(outcome_type, method, cv_method)
        result <- cv_imr(fit, k = 2, rounds = 1, max_models = 4,
                         cv_method = cv_method)
        expect_identical(result$validation, cv_method, info = context)
        expect_true(all(is.finite(result$predictions$prediction)), info = context)
        expect_true(all(is.finite(result$pooled)), info = context)
        expect_true(all(is.finite(result$fold_mean)), info = context)
        expect_equal(nrow(result$predictions), sum(fit$model$sample_sizes),
                     info = context)
        expect_equal(anyDuplicated(result$predictions$id), 0L, info = context)
        if (outcome_type == "continuous") {
          rows <- match(result$predictions$id, outcomes[[outcome_type]][[1]])
          mse <- mean((result$predictions$prediction -
                        outcomes[[outcome_type]][[2]][rows])^2)
          expect_equal(unname(result$pooled[1, "all"]), mse, info = context)
        } else {
          expect_true(all(result$pooled >= 0 & result$pooled <= 1), info = context)
        }
      }
    }
  }
})
