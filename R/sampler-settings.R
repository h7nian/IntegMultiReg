# One description of the settings each transition kernel actually uses.
.imr_sampler_settings <- function(marginalize, mcmc, numerical, model_variant) {
  laplace <- marginalize == "coefficients_and_variance"
  mcmc_fields <- c(
    "draws", "burnin", "chains", "thin", "seed", "workers",
    "keep_latent", "initial", "diagnostics", "max_draw_memory_mb"
  )
  if (!laplace) mcmc_fields <- c(mcmc_fields, "swap_rate")
  if (!laplace && model_variant == "imr") mcmc_fields <- c(mcmc_fields, "theta_step")
  if (marginalize == "coefficients") mcmc_fields <- c(mcmc_fields, "variance_step")
  numerical_fields <- if (laplace) c("laplace_max_iter", "laplace_tolerance") else if (marginalize == "coefficients") "max_integration_nodes" else character()
  .imr_check_applicable(mcmc, imr_mcmc(), mcmc_fields, paste(marginalize, model_variant, sep = "/"))
  .imr_check_applicable(numerical, imr_control(), numerical_fields, marginalize)
  list(mcmc = mcmc[mcmc_fields], numerical = numerical[numerical_fields])
}

.imr_check_applicable <- function(settings, defaults, applicable, context) {
  unused <- setdiff(names(settings), applicable)
  changed <- unused[!vapply(unused, function(name) identical(settings[[name]], defaults[[name]]), TRUE)]
  if (length(changed)) {
    .imr_abort(sprintf(
      "Setting%s %s %s not applicable to '%s'. Leave %s at the default.",
      if (length(changed) > 1L) "s" else "", paste(sprintf("`%s`", changed), collapse = ", "),
      if (length(changed) > 1L) "are" else "is", context, if (length(changed) > 1L) "them" else "it"
    ))
  }
  invisible(NULL)
}

.imr_acceptance_table <- function(chains, marginalize, outcome, sharing) {
  updates <- c("selection", "swap", "interaction", "variance", "latent")
  laplace <- marginalize == "coefficients_and_variance"
  methods <- c(
    if (laplace) "metropolis" else "gibbs",
    if (laplace) "included_in_selection" else "metropolis",
    if (sharing) "metropolis" else "fixed",
    if (marginalize %in% c("variance", "coefficients_and_variance")) "integrated" else if (marginalize == "coefficients") "metropolis" else "gibbs",
    if (outcome == "continuous") "observed" else if (marginalize %in% c("coefficients", "coefficients_and_variance")) "metropolis" else "gibbs"
  )
  result <- lapply(seq_along(chains), function(i) {
    counts <- chains[[i]]$acceptance
    if (length(counts) && is.null(names(counts))) {
      names(counts) <- c("swap_proposals", "swap_accepts", "interaction_proposals", "interaction_accepts")
    }
    get_count <- function(update, suffix, method) {
      name <- paste0(update, suffix)
      if (method != "metropolis" || !name %in% names(counts)) NA_real_ else unname(counts[[name]])
    }
    proposed <- mapply(get_count, updates, "_proposals", methods, USE.NAMES = FALSE)
    accepted <- mapply(get_count, updates, "_accepts", methods, USE.NAMES = FALSE)
    data.frame(
      chain = i, update = updates, method = methods, proposed = proposed,
      accepted = accepted, rate = ifelse(proposed > 0, accepted / proposed, NA_real_),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}
