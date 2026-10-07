#' Priors for an Integrative Regression Model
#'
#' Groups the statistical prior choices independently of MCMC settings.
#' @param nu Baseline inclusion log-odds, conditional on the other subgroup
#'   indicators being zero: one number or one per platform. MRF interactions
#'   also affect marginal selection probabilities.
#' @param molecular_scale pMOM scale for molecular regression coefficients.
#' @param forced_scale pMOM scale for the intercept and clinical coefficients.
#' @param residual Inverse-gamma shape and rate, or `NULL` for family defaults.
#'   Gaussian and survival defaults are `c(shape = .001, rate = .001)`.
#'   Binary models anchor the latent scale with `c(shape = 1e5, rate = 1e5)`;
#'   a conflicting explicit binary variance prior is rejected.
#' @param interaction Gamma shape and rate for positive MRF interactions.
#' @return A validated `imr_priors` specification for [imr()].
#' @export
#' @examples
#' imr_priors(nu = c(-4, -3, -4))
imr_priors <- function(nu = -3, molecular_scale = .087, forced_scale = 10000,
                       residual = NULL,
                       interaction = c(shape = 40, rate = 10)) {
  if (!is.numeric(nu) || !length(nu) || any(!is.finite(nu))) {
    .imr_abort("`nu` must contain finite prior log-odds.")
  }
  structure(list(
    nu = as.double(nu),
    molecular_scale = as.double(.imr_check_numeric_vector(molecular_scale, "molecular_scale", 1, positive = TRUE)),
    forced_scale = as.double(.imr_check_numeric_vector(forced_scale, "forced_scale", 1, positive = TRUE)),
    residual = if (is.null(residual)) NULL else .imr_check_named_pair(residual, "residual"),
    interaction = .imr_check_named_pair(interaction, "interaction")
  ), class = "imr_priors")
}

#' MCMC Settings for Joint Posterior Sampling
#'
#' `draws` is the number of retained iterations per chain. Each chain performs
#' `burnin + draws * thin` updates. Thinning reduces stored draws, not the cost
#' of their intervening updates, and is not a remedy for poor mixing.
#' @param draws Retained iterations per chain, at least four.
#' @param burnin Initial iterations discarded from each chain.
#' @param chains Number of separately seeded chains. Four is the default;
#'   R-hat requires at least two.
#' @param thin Updates between retained iterations.
#' @param seed Nonnegative seed, or `NULL` to draw a seed from R's current RNG.
#' @param workers Number of worker processes for chains. CV distributes folds
#'   instead and runs their chains serially, avoiding nested parallelism.
#' @param keep_latent Retain augmented responses for binary/censored models.
#' @param initial `"dispersed"`, `"empty"`, or `"full"` selection starts, or
#'   a list of explicitly named starting states, one per chain.
#' @param theta_step Standard deviation of the log-interaction proposal.
#' @param swap_rate Whole-feature pair proposal rate. Each update attempts
#'   `ceiling(swap_rate * number_of_platform_features)` exchanges per platform;
#'   zero disables exchanges.
#' @param diagnostics Compute rank R-hat, effective sample sizes and MCSE after
#'   fitting. These can also be computed later by [mcmc_diagnostics()].
#' @param max_draw_memory_mb Limit on estimated retained-array size in MiB.
#'   This is not a limit on total process memory; fitting/combining draws needs
#'   additional workspace. Increase explicitly for large stored posteriors.
#' @return A validated `imr_mcmc` specification for [imr()].
#' @export
#' @examples
#' imr_mcmc(draws = 4000, burnin = 2000, chains = 4, seed = 123)
imr_mcmc <- function(draws = 2000L, burnin = 1000L, chains = 4L, thin = 1L,
                     seed = NULL, workers = 1L, keep_latent = FALSE,
                     initial = "dispersed", theta_step = .4, swap_rate = .5,
                     diagnostics = TRUE, max_draw_memory_mb = 1024) {
  draws <- .imr_check_integer_scalar(draws, "draws", min = 4)
  burnin <- .imr_check_integer_scalar(burnin, "burnin", min = 0)
  chains <- .imr_check_integer_scalar(chains, "chains", min = 1)
  thin <- .imr_check_integer_scalar(thin, "thin", min = 1)
  workers <- .imr_check_integer_scalar(workers, "workers", min = 1)
  if (as.double(draws) * chains > .Machine$integer.max) {
    .imr_abort("The total retained draw count exceeds the supported matrix row limit.")
  }
  if (as.double(draws) * thin + burnin > .Machine$integer.max) {
    .imr_abort("`burnin + draws * thin` exceeds the supported iteration limit.")
  }
  if (!is.null(seed)) seed <- .imr_check_integer_scalar(seed, "seed", min = 0)
  .imr_check_flag(keep_latent, "keep_latent")
  .imr_check_flag(diagnostics, "diagnostics")
  if (is.character(initial)) {
    initial <- match.arg(initial, c("dispersed", "empty", "full"))
  } else if (!is.list(initial) || length(initial) != chains) {
    .imr_abort("Explicit `initial` states must be a list with one element per chain.")
  }
  theta_step <- .imr_check_numeric_vector(theta_step, "theta_step", 1, positive = TRUE)
  swap_rate <- .imr_check_numeric_vector(swap_rate, "swap_rate", 1, nonnegative = TRUE)
  if (swap_rate > 10) .imr_abort("`swap_rate` must not exceed 10.")
  max_draw_memory_mb <- .imr_check_numeric_vector(max_draw_memory_mb, "max_draw_memory_mb", 1, positive = TRUE)
  structure(list(
    draws = draws, burnin = burnin, chains = chains, thin = thin,
    seed = seed, workers = workers, keep_latent = keep_latent, initial = initial,
    theta_step = theta_step, swap_rate = swap_rate, diagnostics = diagnostics,
    max_draw_memory_mb = max_draw_memory_mb
  ), class = "imr_mcmc")
}

.imr_validate_specification <- function(x, constructor, label) {
  expected <- names(formals(constructor))
  if (!is.list(x) || is.null(names(x)) || !setequal(names(x), expected) || anyDuplicated(names(x))) {
    .imr_abort(sprintf("Use `%s()` to construct `%s`.", label, label))
  }
  do.call(constructor, unclass(x))
}
