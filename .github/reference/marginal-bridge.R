# Adapt the public sampler specification to independently integrated fixtures.
native_chain <- function(groups, feature_platform, nu, ...) {
  name <- if (Sys.getenv("IMR_REFERENCE_MARGINALIZE") == "variance") {
    ".imr_variance_marginal_sample"
  } else ".imr_coefficients_marginal_sample"
  result <- get(name, asNamespace("IntegMultiReg"))(groups, feature_platform, nu, ...)
  selection <- matrix(NA_integer_, nrow(result$variance), length(groups) * length(feature_platform))
  for (g in seq_along(groups)) {
    columns <- (g - 1L) * length(feature_platform) + groups[[g]]$feature_index
    selection[, columns] <- result$coefficients[[g]][,
      groups[[g]]$n_forced + seq_along(groups[[g]]$feature_index), drop = FALSE] != 0
  }
  list(beta = result$coefficients, variance = result$variance, selection = selection,
    theta = do.call(cbind, result$interaction), latent = result$latent,
    log_density = result$log_density, diagnostics = result$acceptance)
}
