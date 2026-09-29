## Shared fixtures for the test suite.  Small MCMC runs keep R CMD check fast
## while still exercising the full fit -> predict -> cv pipeline.

data("simIMR", package = "IntegMultiReg")

fit_demo <- function(type = c("binary", "continuous", "right.censored"),
                     total = 300, burn = 150, seed = 42) {
  type <- match.arg(type)
  outcome <- switch(type,
    binary = simIMR$outcome.binary,
    continuous = simIMR$outcome.continuous,
    right.censored = simIMR$outcome.survival
  )
  imr(
    x = simIMR$platforms,
    outcome = outcome,
    covariates = simIMR$covariates,
    outcome_type = type,
    nu = c(-4, -3, -4),
    draws = total, burnin = burn,
    min_subgroup_size = 30,
    seed = seed
  )
}

## A binary fit reused across several test files.
fit_bin <- fit_demo("binary", seed = 42)
