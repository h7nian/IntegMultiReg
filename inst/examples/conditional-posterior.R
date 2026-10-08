\donttest{
x <- data.frame(id = 1:20, marker = sin(1:20))
y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
fit <- imr(list(assay = x), y,
  outcome_type = "continuous", min_subgroup_size = 0,
  marginalize = "coefficients_and_variance",
  priors = imr_priors(forced_scale = 1),
  mcmc = imr_mcmc(draws = 1000, burnin = 1000, chains = 2, seed = 1)
)
posterior <- sample_regression_posterior(fit,
  output_draws = 1000,
  mcmc = imr_mcmc(draws = 1000, burnin = 1000, chains = 2, seed = 2)
)
summary(posterior)
coef(posterior)
confint(posterior)
predict(posterior, list(assay = x), interval = TRUE)
# Inspect both the conditional chains and their source selection fit.
mcmc_diagnostics(posterior)
mcmc_diagnostics(fit)
}
