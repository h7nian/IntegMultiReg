# Independently installed frozen candidate and repaired candidate, same runtime.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, args[3] %in% c("capture", "compare"))
.libPaths(c(args[1], .libPaths()))
library(IntegMultiReg)
data(simIMR)
results <- list()
for (type in c("binary", "continuous", "right.censored")) {
  y <- switch(type, binary = simIMR$outcome.binary,
    continuous = simIMR$outcome.continuous, right.censored = simIMR$outcome.survival)
  for (method in c("imr", "bms")) {
    fit <- imr(simIMR$platforms, y, covariates = simIMR$covariates,
      outcome_type = type, method = method, nu = c(-4, -3, -4),
      draws = 100, burnin = 50, seed = 53)
    cv <- lapply(c("legacy", "importance", "refit"), function(mode)
      cv_imr(fit, k = 3, rounds = 1, cv_method = mode)[1:5])
    results[[paste(type, method)]] <- list(posterior = fit$posterior,
      model = fit$model, preprocessing = fit$preprocessing,
      predictions = predict(fit, simIMR$platforms,
                            covariates = simIMR$covariates), cv = cv)
  }
}
if (args[3] == "capture") saveRDS(results, args[2]) else {
  stopifnot(identical(results, readRDS(args[2])))
  cat("PASS: exact default numerical regression, six fits and eighteen CV combinations.\n")
}
