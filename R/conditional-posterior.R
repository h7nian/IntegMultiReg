#' Sample Regression Parameters Conditional on Stored Selection States
#'
#' Retains the conditional pMOM Gibbs procedure for a Laplace selection fit.
#' One source selection-state index is sampled for every output row and shared
#' across subgroups, preserving model-selection dependence. Conditional chains
#' use the observed data within each selected subgroup model.
#'
#' @section Interpretation:
#' This is a mixture of conditional regression posteriors with empirical model
#' weights from the original selection fit. Those weights retain the Laplace
#' approximation and the selection chain's finite-sample limitations. Conditional
#' diagnostics assess the separate fixed-model chains; the returned mixture rows
#' must not be treated as additional chains of the selection sampler.
#' For a fixed model with d active coefficients, the variance update is
#' \deqn{\sigma^2\mid\cdots\sim IG\left(a+\frac n2+\frac{3d}{2},
#' b+\frac12\left[\|z-X\beta\|^2+\sum_j\beta_j^2/\tau_j\right]\right).}{sigma^2 | rest ~ inverse-Gamma(a + n/2 + 3d/2, b + (residual sum of squares + sum_j beta_j^2/tau_j)/2).}
#'
#' @param object An `imr_selection` fit. The other three marginalization choices
#'   already store regression posterior draws, available with [posterior_draws()].
#' @param output_draws Number of returned mixture rows.
#' @param burnin Discarded iterations of each conditional chain.
#' @param chains Conditional chains per subgroup selection model, at least two.
#' @param min_draws_per_model_chain Minimum retained length of each conditional
#'   chain; frequent models receive additional draws when needed for output.
#' @param seed Seed for the conditional computation; the caller's RNG is restored.
#' @param latent Retain augmented binary/censored responses paired with coefficients.
#' @return An `imr_posterior` object with coefficient, variance and optional latent
#'   draws, source selection indices, conditional diagnostics and the originating fit.
#' @example inst/examples/conditional-posterior.R
#' @export
sample_regression_posterior <- function(object, output_draws = 1000L, burnin = 1000L,
                                        chains = 2L, min_draws_per_model_chain = 200L, seed = 1L,
                                        latent = FALSE) {
  validate_imr_object(object)
  if (!inherits(object, "imr_selection")) .imr_abort("This fit already contains regression posterior draws; use posterior_draws().")
  output_draws <- .imr_check_integer_scalar(output_draws, "output_draws", min = 2L)
  burnin <- .imr_check_integer_scalar(burnin, "burnin", min = 0L)
  chains <- .imr_check_integer_scalar(chains, "chains", min = 2L)
  min_draws_per_model_chain <- .imr_check_integer_scalar(
    min_draws_per_model_chain, "min_draws_per_model_chain",
    min = 4L
  )
  seed <- .imr_check_integer_scalar(seed, "seed", min = 0L)
  .imr_check_flag(latent, "latent")
  # For a continuous outcome the response is observed, so there is nothing
  # latent to return.
  keep_latent <- latent && object$control$outcome_type != "continuous"
  if (object$control$outcome_type == "right.censored" &&
    (length(object$control$response_scale) != 1L ||
      !object$control$response_scale %in% c("log", "identity"))) {
    .imr_abort("Refit this survival model with an explicit `survival_scale`.")
  }
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(seed)
  selection_draw_index <- sample.int(object$control$mcmc$draws * object$control$mcmc$chains, output_draws, replace = TRUE)
  selected_states <- .imr_selection_history(object, selection_draw_index)
  beta <- variance <- diagnostics <- augmented <-
    vector("list", length(object$model$subgroup_names))
  priors <- object$control$priors
  for (g in seq_along(beta)) {
    design <- .imr_joint_design(object$model, object$preprocessing, g)
    masks <- lapply(selected_states, function(state) {
      c(rep(TRUE, 1L + length(object$model$covariate_names)), unlist(lapply(
        object$model$subgroup_platforms[[g]], function(p) {
          state[[p]][match(g, object$model$platform_subgroups[[p]]), ] == 1
        }
      ), use.names = FALSE))
    })
    keys <- vapply(masks, function(x) paste(as.integer(x), collapse = ""), "")
    beta[[g]] <- matrix(0, output_draws, ncol(design), dimnames = list(NULL, colnames(design)))
    variance[[g]] <- numeric(output_draws)
    if (keep_latent) {
      augmented[[g]] <- matrix(NA_real_, output_draws,
        nrow(object$preprocessing$response[[g]]),
        dimnames = list(NULL, rownames(object$preprocessing$response[[g]]))
      )
    }
    records <- list()
    for (key in unique(keys)) {
      positions <- which(keys == key)
      active <- masks[[positions[1L]]]
      X <- design[, active, drop = FALSE]
      h <- c(
        rep(priors$forced_scale, 1L + length(object$model$covariate_names)),
        rep(priors$molecular_scale, ncol(design) - 1L - length(object$model$covariate_names))
      )[active]
      n <- max(min_draws_per_model_chain, ceiling(length(positions) / chains))
      y <- object$preprocessing$response[[g]][, 1L]
      status <- if (object$control$outcome_type == "right.censored") object$preprocessing$response[[g]][, 2L] else NULL
      samples <- lapply(seq_len(chains), function(chain) {
        .imr_conditional_chain(X, y, h,
          priors$residual[["shape"]], priors$residual[["rate"]],
          draws = n, burnin = burnin,
          initial_beta = rep(if (chain %% 2L) -.5 else .5, ncol(X)),
          initial_variance = 1, outcome_type = object$control$outcome_type,
          status = status, keep_latent = keep_latent
        )
      })
      pool <- do.call(rbind, samples)
      chosen <- sample.int(nrow(pool), length(positions), replace = FALSE)
      beta[[g]][positions, active] <- pool[chosen, seq_len(ncol(X)), drop = FALSE]
      variance[[g]][positions] <- pool[chosen, ncol(X) + 1L]
      if (keep_latent) {
        # Index the pooled chains with the same rows, so that a latent draw and
        # the coefficient draw beside it come from one sweep of the sampler.
        pooled_latent <- do.call(rbind, lapply(samples, attr, "latent"))
        augmented[[g]][positions, ] <- pooled_latent[chosen, , drop = FALSE]
      }
      records[[length(records) + 1L]] <- .imr_conditional_diagnostics(
        samples, object, g, key, length(positions)
      )
    }
    diagnostics[[g]] <- do.call(rbind, records)
  }
  names(beta) <- names(variance) <- object$model$subgroup_names
  if (keep_latent) names(augmented) <- object$model$subgroup_names
  diagnostics <- do.call(rbind, diagnostics)
  if (any(diagnostics$rhat > 1.01, na.rm = TRUE) ||
    any(diagnostics$status == "constant_in_chain")) {
    .imr_warn("Conditional-chain diagnostics need inspection; increase the per-model budget and inspect the source selection chains separately.")
  }
  structure(
    list(
      coefficients = beta, variance = variance,
      latent = if (keep_latent) augmented else NULL,
      selection_draw_index = selection_draw_index,
      diagnostics = diagnostics, fit = object,
      control = list(
        output_draws = output_draws, burnin = burnin, chains = chains,
        min_draws_per_model_chain = min_draws_per_model_chain, seed = seed
      ),
      approximation = "Empirical selection weights from the Laplace-based fit; conditional pMOM Gibbs draws."
    ),
    class = "imr_posterior"
  )
}
