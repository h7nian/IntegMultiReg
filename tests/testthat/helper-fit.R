data("simIMR", package = "IntegMultiReg")

fit_demo <- function(type = c("binary", "continuous", "right.censored"),
                     total = 40L, burn = 20L, seed = 42L, variant = "imr",
                     chains = 2L, latent = FALSE, workers = 1L) {
  type <- match.arg(type)
  outcome <- switch(type,
    binary = simIMR$outcome.binary,
    continuous = simIMR$outcome.continuous,
    right.censored = simIMR$outcome.survival
  )
  imr(simIMR$platforms, outcome,
    covariates = simIMR$covariates,
    outcome_type = type, model_variant = variant,
    priors = imr_priors(nu = c(-4, -3, -4)),
    mcmc = imr_mcmc(
      draws = total, burnin = burn, seed = seed,
      chains = chains, workers = workers, keep_latent = latent, diagnostics = FALSE
    )
  )
}
fit_bin <- fit_demo("binary")

small_data <- function(type = "continuous", n = 20L) {
  id <- seq_len(n)
  list(
    platforms = list(assay = data.frame(id = id, marker = sin(id))),
    outcome = switch(type,
      continuous = data.frame(id = id, y = .7 + .8 * sin(id) + cos(id) / 3),
      binary = data.frame(id = id, y = rep(0:1, length.out = n)),
      right.censored = data.frame(id = id, time = exp(.7 + sin(id)), status = rep(c(1, 1, 0), length.out = n))
    )
  )
}
small_fit <- function(type = "continuous", draws = 40L, burnin = 20L,
                      chains = 2L, seed = 17L, workers = 1L, keep_latent = FALSE, initial = "dispersed") {
  d <- small_data(type)
  imr(d$platforms, d$outcome,
    outcome_type = type, min_subgroup_size = 0,
    priors = imr_priors(forced_scale = 1),
    mcmc = imr_mcmc(
      draws = draws, burnin = burnin, chains = chains, seed = seed,
      workers = workers, keep_latent = keep_latent, initial = initial, diagnostics = FALSE
    )
  )
}
