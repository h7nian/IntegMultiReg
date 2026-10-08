.imr_conditional_diagnostics <- function(samples, fit, subgroup, model, output_draws, workers = 1L) {
  rows <- nrow(samples[[1L]])
  parameters <- ncol(samples[[1L]]) - 1L
  array_for <- function(matrices, labels) {
    out <- aperm(array(
      unlist(matrices, use.names = FALSE),
      c(rows, length(labels), length(matrices))
    ), c(1L, 3L, 2L))
    dimnames(out) <- list(NULL, paste0("chain", seq_along(matrices)), labels)
    out
  }
  blocks <- list(
    list(
      family = "coefficients", group = fit$model$subgroup_names[subgroup],
      draws = array_for(lapply(samples, function(x) x[, seq_len(parameters), drop = FALSE]), colnames(samples[[1L]])[seq_len(parameters)]),
      observed = rep(FALSE, parameters)
    ),
    list(
      family = "variance", group = fit$model$subgroup_names[subgroup],
      draws = array_for(lapply(samples, function(x) x[, parameters + 1L, drop = FALSE]), "residual_variance"), observed = FALSE
    )
  )
  if (!is.null(attr(samples[[1L]], "latent"))) {
    latent <- array_for(lapply(samples, attr, "latent"), as.character(fit$preprocessing$subject_ids[[subgroup]]))
    observed <- if (fit$control$outcome_type == "right.censored") {
      fit$preprocessing$response[[subgroup]][, 2L] == 1
    } else {
      rep(FALSE, dim(latent)[3L])
    }
    blocks[[3L]] <- list(family = "latent", group = fit$model$subgroup_names[subgroup], draws = latent, observed = observed)
  }
  result <- .imr_diagnose_blocks(blocks, workers)
  result$selection_model <- model
  result$output_draws <- output_draws
  result$draws_per_model_chain <- rows
  rownames(result) <- NULL
  result
}

.imr_check_regression_posterior <- function(object) {
  if (!inherits(object, "imr_posterior") || !is.list(object)) .imr_abort("An imr_posterior object is required.")
  .imr_check_fit(object$fit)
  groups <- object$fit$model$subgroup_names
  n <- object$control$output_draws
  if (!is.list(object$coefficients) || !identical(names(object$coefficients), groups) ||
    !is.list(object$variance) || !identical(names(object$variance), groups) ||
    length(n) != 1L || !.imr_is_integerish(n) || n < 2L) {
    .imr_abort("Invalid conditional-posterior structure.")
  }
  for (g in seq_along(groups)) {
    coefficients <- object$coefficients[[g]]
    labels <- colnames(.imr_joint_design(object$fit$model, object$fit$preprocessing, g))
    if (!is.matrix(coefficients) || nrow(coefficients) != n || !identical(colnames(coefficients), labels) ||
      any(!is.finite(coefficients)) || length(object$variance[[g]]) != n ||
      any(!is.finite(object$variance[[g]]) | object$variance[[g]] <= 0)) {
      .imr_abort("Invalid conditional coefficient or variance draws.")
    }
    if (!is.null(object$latent)) {
      x <- object$latent[[g]]
      if (!is.matrix(x) || !identical(dim(x), c(as.integer(n), as.integer(object$fit$model$sample_sizes[g]))) ||
        any(!is.finite(x))) {
        .imr_abort("Invalid conditional latent draws.")
      }
    }
  }
  total <- as.double(object$fit$control$mcmc$draws) * object$fit$control$mcmc$chains
  if (length(object$selection_draw_index) != n || !.imr_is_integerish(object$selection_draw_index) ||
    any(object$selection_draw_index < 1 | object$selection_draw_index > total)) {
    .imr_abort("Invalid source selection indices.")
  }
  invisible(TRUE)
}

# A computation view for shared summaries/prediction, not a fitted Markov chain.
# Its second array dimension is explicitly called mixture. Diagnostics always
# come from the separate fixed-model chains, never from this view.
.imr_conditional_view <- function(object) {
  fit <- object$fit
  n <- object$control$output_draws
  as_draws <- function(x, labels = colnames(x)) {
    array(x, c(n, 1L, ncol(x)), dimnames = list(NULL, "mixture", labels))
  }
  variance <- array(unlist(object$variance, use.names = FALSE),
    c(n, 1L, length(fit$model$subgroup_names)),
    dimnames = list(NULL, "mixture", fit$model$subgroup_names)
  )
  control <- fit$control
  control$conditional_regression <- TRUE
  control$mcmc$draws <- n
  control$mcmc$chains <- 1L
  list(
    model = fit$model, preprocessing = fit$preprocessing, control = control,
    posterior = list(
      coefficients = lapply(object$coefficients, as_draws), variance = variance,
      latent = if (is.null(object$latent)) {
        NULL
      } else {
        stats::setNames(lapply(seq_along(object$latent), function(g) {
          as_draws(object$latent[[g]], as.character(fit$preprocessing$subject_ids[[g]]))
        }), fit$model$subgroup_names)
      }
    )
  )
}

#' Summarize Conditional Regression Posterior Draws
#'
#' Summarizes an `imr_posterior` object without further sampling. Intervals use
#' the same empirical inverse-CDF convention as [confint.imr()] and include
#' exclusion zeros. Conditional-chain diagnostics are available separately with
#' [mcmc_diagnostics()]; mixture output rows are not additional MCMC chains.
#' @param object,x An `imr_posterior` object from [sample_regression_posterior()].
#' @param parm `"all"` (summary default), `"coefficients"` (interval default),
#'   `"variance"`, or `"latent"`.
#' @param level Equal-tail probability level.
#' @param ... Unused arguments are rejected.
#' @return `summary()` returns a `summary.imr` object with a parameters table;
#'   `confint()` returns the interval table. `coef()` returns
#'   named coefficient means by subgroup; `print()` returns its input invisibly.
#' @example inst/examples/conditional-posterior.R
#' @name imr_posterior_methods
NULL

#' @rdname imr_posterior_methods
#' @export
summary.imr_posterior <- function(object, parm = "all", level = .95, ...) {
  .imr_reject_dots(...)
  .imr_check_regression_posterior(object)
  level <- .imr_check_level(level)
  parm <- match.arg(parm, c("all", "coefficients", "variance", "latent"))
  if (parm == "latent" && is.null(object$latent)) .imr_abort("No conditional latent draws were retained.")
  .imr_summary(.imr_conditional_view(object), parm, level, conditional = TRUE)
}

#' @rdname imr_posterior_methods
#' @export
confint.imr_posterior <- function(object, parm = "coefficients", level = .95, ...) {
  .imr_reject_dots(...)
  .imr_check_regression_posterior(object)
  .imr_intervals(.imr_conditional_view(object), parm, .imr_check_level(level))
}

#' @rdname imr_posterior_methods
#' @export
coef.imr_posterior <- function(object, ...) {
  .imr_reject_dots(...)
  .imr_check_regression_posterior(object)
  lapply(object$coefficients, colMeans)
}

#' @rdname imr_posterior_methods
#' @export
print.imr_posterior <- function(x, ...) {
  .imr_reject_dots(...)
  .imr_check_regression_posterior(x)
  cat("Conditional regression posterior:", x$control$output_draws, "mixture draws\n")
  cat(x$approximation, "\n")
  if (is.null(x$diagnostics)) cat("Conditional diagnostics were not computed.\n") else .imr_print_diagnostics(x$diagnostics)
  cat("Diagnostics refer to fixed-model conditional chains; inspect the source fit separately.\n")
  invisible(x)
}
