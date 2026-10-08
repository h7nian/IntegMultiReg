.imr_marginalize <- function(control) control$marginalize %||% "none"

.imr_sampler_description <- function(control) {
  switch(.imr_marginalize(control),
    none = "joint posterior fit",
    variance = "variance-marginal posterior fit",
    coefficients = "coefficient-marginal posterior fit",
    coefficients_and_variance = "Laplace-marginal selection fit"
  )
}

.imr_sample_chain <- function(task, spec, model, control, verbose = FALSE) {
  marginalize <- .imr_marginalize(control)
  if (marginalize == "none") {
    return(.imr_joint_chain(task, spec, model, control, verbose))
  }
  saved_rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(saved_rng), add = TRUE)
  requested <- control$mcmc$initial
  if (is.list(requested)) requested <- requested[[task$chain]]
  set.seed(task$initial_seed)
  initial <- if (marginalize == "coefficients_and_variance" && is.null(requested)) {
    NULL
  } else {
    .imr_chain_start(task$chain, spec, model, control)
  }
  set.seed(task$seed)
  if (marginalize == "coefficients_and_variance") {
    native_initial <- if (is.null(initial)) {
      NULL
    } else {
      selected <- lapply(spec$platforms, function(platform) {
        value <- matrix(0L, length(platform$groups), ncol(platform$columns))
        for (i in seq_along(platform$groups)) {
          value[i, ] <- as.integer(initial$coefficients[[platform$groups[i] + 1L]][platform$columns[i, ] + 1L] != 0)
        }
        value
      })
      list(selection = selected, interaction = if (control$model_variant == "imr") initial$interaction else NULL)
    }
    raw <- .imr_call_collapsed(model, spec$preprocessing, control, task$seed,
      initial = native_initial, verbose = verbose
    )
    # Preserve the original public fit's guard against running-average drift.
    raw$gam_mean <- lapply(raw$gam_mean, function(x) {
      x[x > 1] <- 1
      x[x < 0] <- 0
      x
    })
    selection <- lapply(seq_len(model$n_platforms), function(l) {
      count <- length(model$platform_subgroups[[l]]) * length(model$feature_names[[l]])
      out <- matrix(0L, control$mcmc$draws, count)
      for (i in seq_len(nrow(out))) out[i, ] <- as.vector(t(raw$gam_sample[[i]][[l]]))
      out
    })
    interaction <- lapply(raw$theta_sample, function(x) {
      if (is.null(x)) matrix(numeric(), control$mcmc$draws, 0L) else x
    })
    return(list(
      coefficients = NULL, variance = NULL, selection = selection,
      interaction = interaction, latent = raw$latent_sample,
      log_density = utils::tail(raw$log_posterior, control$mcmc$draws),
      selection_mean = raw$gam_mean, interaction_mean = raw$theta_mean,
      latent_mean = raw$estimate_latent_y, laplace_diagnostics = raw$laplace_diagnostics,
      initial = initial, acceptance = numeric()
    ))
  }
  offsets <- c(0L, cumsum(lengths(model$feature_names)))
  groups <- lapply(seq_along(spec$groups), function(s) {
    g <- spec$groups[[s]]
    features <- unlist(lapply(model$subgroup_platforms[[s]], function(l) {
      as.integer(offsets[l] + seq_along(model$feature_names[[l]]))
    }), use.names = FALSE)
    list(
      X = g$design, y = g$response, status = g$status,
      prior_scale = g$prior_scale, n_forced = g$forced, feature_index = features,
      outcome_type = control$outcome_type,
      residual_prior = stats::setNames(g$residual_prior, c("shape", "rate"))
    )
  })
  names(groups) <- model$subgroup_names
  arguments <- list(
    groups = groups,
    feature_platform = rep(seq_len(model$n_platforms), lengths(model$feature_names)),
    nu = control$priors$nu, draws = control$mcmc$draws, burnin = control$mcmc$burnin,
    thin = control$mcmc$thin, seed = task$seed, model_variant = control$model_variant,
    keep_latent = control$mcmc$keep_latent, interaction_prior = control$priors$interaction,
    theta_step = control$mcmc$theta_step, swap_rate = control$mcmc$swap_rate,
    initial_state = initial
  )
  engine <- .imr_variance_marginal_sample
  if (marginalize == "coefficients") {
    engine <- .imr_coefficients_marginal_sample
    arguments$max_nodes <- control$numerical$max_integration_nodes
    arguments$variance_step <- control$mcmc$variance_step
  }
  value <- do.call(engine, arguments)
  value$initial <- initial
  value
}
