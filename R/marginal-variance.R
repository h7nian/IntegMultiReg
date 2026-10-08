# Variance-marginal posterior updates and exact conditional recovery.
.imr_variance_marginal_sample <- function(groups, feature_platform, nu, draws = 30000L,
                                          burnin = 2000L, seed = 1L, model_variant = "imr",
                                          initial = "empty", keep_latent = TRUE,
                                          interaction_prior = c(shape = 2, rate = 2),
                                          thin = 1L, theta_step = .4, swap_rate = .5, initial_state = NULL, verbose = FALSE) {
  .imr_marginal_chain(groups, feature_platform, nu,
    draws = draws, burnin = burnin, seed = seed, model_variant = model_variant,
    initial = initial, keep_latent = keep_latent, interaction_prior = interaction_prior,
    thin = thin, theta_step = theta_step, swap_rate = swap_rate,
    initial_state = initial_state, verbose = verbose,
    marginalize = "variance"
  )
}
