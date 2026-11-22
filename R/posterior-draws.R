#' Sample Regression Coefficients and Residual Variances
#'
#' Augments the retained selection models with conditional coefficient and
#' variance draws. The fitted sampler and cross-validation results are unchanged.
#'
#' @param object An `imr` fit with stored data and selection draws.
#' @param draws Number of model-averaged draws to return (default `1000`).
#' @param burnin Conditional Gibbs burn-in iterations for each distinct
#'   subgroup selection model and each chain (default `1000`).
#' @param chains Number of conditional chains, at least two (default `2`).
#' @param conditional_draws Minimum retained iterations per conditional chain
#'   used for sampling and split R-hat (default `200`). Increase with `burnin`
#'   if the reported conditional diagnostics are poor.
#' @param seed Integer seed. The caller's random-number state is restored.
#'
#' @details
#' All active coefficients, including the always-included intercept and clinical
#' effects, follow the original first-order pMOM prior. Inactive molecular
#' coefficients are exactly zero. Molecular and clinical slopes are on the
#' subgroup-standardized predictor scale; responses are not standardized.
#'
#' Whole selection draws are resampled uniformly from `object$gam_sample`,
#' preserving their empirical joint distribution across subgroups. For each
#' distinct subgroup model, a Gibbs sampler draws coefficients and variance
#' conditional on the observed data. Binary and censored responses are augmented
#' anew; posterior mean latent responses are not treated as observed data.
#' The binary variance uses the same highly concentrated inverse-gamma prior
#' as the fitted model (approximately, not exactly, one).
#'
#' This is a two-stage posterior approximation: model weights inherit the
#' original sampler's Laplace approximation and finite-chain exploration. The
#' conditional Gibbs draws do not make those weights exact. Classical split
#' R-hat is reported for each conditional model and does not diagnose the
#' original selection chain. Inspect both stages and increase simulation effort
#' before interpreting intervals. Runtime grows with the number of distinct
#' subgroup selection models, not just `draws`.
#'
#' Old survival fits without explicit response-scale metadata must be refitted.
#' Log-time fits return coefficients for log time; identity-scale fits are
#' retained only for historical compatibility.
#'
#' @return An `imr_posterior` object containing `beta` (one draws-by-coefficients
#'   matrix per subgroup), `variance`, source `model_draw` indices, conditional
#'   `diagnostics`, and the originating `fit`. Coefficient column names distinguish
#'   clinical variables from platform features. Use `summary()`, `confint()`,
#'   `coef()` and `predict()` on this object.
#' @references
#' Chekouo et al. (2017). \doi{10.1111/biom.12587}, Section 3.1 and Web Appendix C.
#' @export
#' @examples
#' \donttest{
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, type_outcome = "continuous",
#'            h0 = 1, sample_mcmc = c(500, 250), seed = 1)
#' draws <- posterior_draws(fit, draws = 1000, burnin = 1000,
#'                          conditional_draws = 1000, seed = 2)
#' summary(draws)
#' }
posterior_draws <- function(object, draws = 1000L, burnin = 1000L,
                            chains = 2L, conditional_draws = 200L, seed = 1L) {
  validate_imr(object)
  draws <- .imr_check_integer_scalar(draws, "draws", min = 2L)
  burnin <- .imr_check_integer_scalar(burnin, "burnin", min = 0L)
  chains <- .imr_check_integer_scalar(chains, "chains", min = 2L)
  conditional_draws <- .imr_check_integer_scalar(
    conditional_draws, "conditional_draws", min = 4L)
  seed <- .imr_check_integer_scalar(seed, "seed", min = 0L)
  if (object$type_outcome == "right.censored" &&
      (length(object$response_scale) != 1L ||
       !object$response_scale %in% c("log", "identity"))) {
    .imr_abort("Refit this survival model with an explicit `survival_scale`.")
  }
  .imr_check_posterior_data(object)
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(seed)
  model_draw <- sample.int(length(object$gam_sample), draws, replace = TRUE)
  beta <- variance <- diagnostics <- vector("list", length(object$model_bitstrings))
  hyper <- object$list_hyperpara
  for (g in seq_along(beta)) {
    design <- .imr_posterior_design(object, g)
    masks <- lapply(model_draw, function(s) {
      c(rep(TRUE, 1L + length(object$covariate_names)), unlist(lapply(
        object$model_platforms[[g]], function(p) {
          object$gam_sample[[s]][[p]][match(g, object$platform_models[[p]]), ] == 1
        }), use.names = FALSE))
    })
    keys <- vapply(masks, function(x) paste(as.integer(x), collapse = ""), "")
    beta[[g]] <- matrix(0, draws, ncol(design), dimnames = list(NULL, colnames(design)))
    variance[[g]] <- numeric(draws)
    records <- list()
    for (key in unique(keys)) {
      positions <- which(keys == key)
      active <- masks[[positions[1L]]]
      X <- design[, active, drop = FALSE]
      h <- c(rep(hyper[1L], 1L + length(object$covariate_names)),
             rep(hyper[2L], ncol(design) - 1L - length(object$covariate_names)))[active]
      n <- max(conditional_draws, ceiling(length(positions) / chains))
      y <- object$data2$yy[[g]][, 1L]
      status <- if (object$type_outcome == "right.censored") object$data2$yy[[g]][, 2L] else NULL
      samples <- lapply(seq_len(chains), function(chain) {
        .imr_conditional_chain(X, y, h, rep(1, ncol(X)), hyper[3L], hyper[4L],
          draws = n, burnin = burnin,
          initial_beta = rep(if (chain %% 2L) -.5 else .5, ncol(X)),
          initial_variance = 1, outcome_type = object$type_outcome, status = status)
      })
      rhat <- .imr_split_rhat(samples)
      pool <- do.call(rbind, samples)
      chosen <- sample.int(nrow(pool), length(positions), replace = FALSE)
      beta[[g]][positions, active] <- pool[chosen, seq_len(ncol(X)), drop = FALSE]
      variance[[g]][positions] <- pool[chosen, ncol(X) + 1L]
      records[[length(records) + 1L]] <- data.frame(
        subgroup = object$model_bitstrings[g], model = key,
        returned_draws = length(positions), conditional_draws = n,
        max_split_rhat = max(rhat), row.names = NULL)
    }
    diagnostics[[g]] <- do.call(rbind, records)
  }
  names(beta) <- names(variance) <- object$model_bitstrings
  diagnostics <- do.call(rbind, diagnostics)
  if (any(!is.finite(diagnostics$max_split_rhat) | diagnostics$max_split_rhat > 1.05)) {
    .imr_warn("Conditional split R-hat exceeds 1.05 or is undefined; inspect `diagnostics` and increase burn-in/draws before using intervals.")
  }
  structure(list(beta = beta, variance = variance, model_draw = model_draw,
    diagnostics = diagnostics, fit = object,
    control = list(draws = draws, burnin = burnin, chains = chains,
                   conditional_draws = conditional_draws, seed = seed),
    approximation = "Empirical selection weights from the Laplace-based fit; conditional pMOM Gibbs draws."),
    class = "imr_posterior")
}

.imr_posterior_design <- function(fit, g, xx = fit$data2$xx, cc = fit$data2$cc) {
  n <- nrow(cc[[g]])
  X <- cbind(rep(1, n), cc[[g]])
  nm <- c("(Intercept)", paste0("clinical:", fit$covariate_names))
  if (!length(fit$covariate_names)) nm <- "(Intercept)"
  for (p in fit$model_platforms[[g]]) {
    X <- cbind(X, xx[[g]][[p]])
    nm <- c(nm, paste0(fit$platform_names[p], ":", fit$feature_names[[p]]))
  }
  colnames(X) <- make.unique(nm)
  X
}

.imr_check_posterior_data <- function(fit) {
  if (length(fit$list_hyperpara) < 4L || any(!is.finite(fit$list_hyperpara[1:4])) ||
      any(fit$list_hyperpara[1:4] <= 0) ||
      !all(c("xx", "yy", "cc") %in% names(fit$data2))) {
    .imr_abort("The fit is missing valid posterior data or hyperparameters.")
  }
  for (g in seq_along(fit$model_bitstrings)) {
    X <- .imr_posterior_design(fit, g)
    y <- fit$data2$yy[[g]]
    if (!is.matrix(y) || nrow(X) != nrow(y) || nrow(y) < 1L ||
        any(!is.finite(X)) || any(!is.finite(y))) {
      .imr_abort("The fit contains inconsistent or non-finite posterior data.")
    }
  }
  invisible(TRUE)
}

.imr_save_rng <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
}
.imr_restore_rng <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", state, envir = .GlobalEnv)
}

#' Summarize Regression Posterior Draws
#' @param object,x An `imr_posterior` object returned by [posterior_draws()].
#' @param level Equal-tail credible level, between zero and one.
#' @param parm Currently only `"coefficients"`.
#' @param ... Unused.
#' @return `summary()` and `confint()` return coefficient tables by subgroup.
#'   `coef()` returns posterior mean coefficient vectors. `print()` returns
#'   its argument invisibly. Effects are on subgroup-standardized predictor
#'   scales. Intervals include zero-mass from inactive selection models.
#' @name imr_posterior_methods
NULL

#' @rdname imr_posterior_methods
#' @export
summary.imr_posterior <- function(object, level = .95, ...) {
  .imr_interval_level(level)
  lapply(object$beta, function(x) {
    q <- t(apply(x, 2L, stats::quantile, probs = c((1-level)/2, .5, (1+level)/2), names = FALSE))
    data.frame(term = colnames(x), mean = colMeans(x), sd = apply(x, 2L, stats::sd),
      lower = q[, 1L], median = q[, 2L], upper = q[, 3L],
      probability_nonzero = colMeans(x != 0), row.names = NULL)
  })
}

#' @rdname imr_posterior_methods
#' @export
confint.imr_posterior <- function(object, parm = "coefficients", level = .95, ...) {
  if (!identical(parm, "coefficients")) .imr_abort("`parm` must be 'coefficients'.")
  summary(object, level = level)
}

#' @rdname imr_posterior_methods
#' @export
coef.imr_posterior <- function(object, ...) lapply(object$beta, colMeans)

#' @rdname imr_posterior_methods
#' @export
print.imr_posterior <- function(x, ...) {
  cat("IMR coefficient posterior:", length(x$model_draw), "draws;",
      length(x$beta), "availability subgroups\n")
  cat(x$approximation, "\n")
  cat("Maximum conditional split R-hat:", max(x$diagnostics$max_split_rhat), "\n")
  invisible(x)
}

.imr_interval_level <- function(level) {
  .imr_check_numeric_vector(level, "level", length = 1L, positive = TRUE)
  if (level >= 1) .imr_abort("`level` must be less than one.")
  invisible(level)
}

#' Predict Using Coefficient Posterior Draws
#' @param object An `imr_posterior` object.
#' @param newdata,platform_names,covariates As in [predict.imr()].
#' @param type `"mean"` returns uncertainty in the conditional response mean
#'   (event probability for binary data). `"response"` additionally generates
#'   new outcomes, including residual variability. Binary response intervals
#'   are discrete 0/1; use `"mean"` for event-probability intervals.
#' @param level Equal-tail interval level.
#' @param seed Integer simulation seed; the caller's RNG state is restored.
#' @param ... Unused.
#' @details Log-time survival fits return time-scale results: for `"mean"`,
#'   each draw is exp(eta + variance/2); for `"response"`, each draw is
#'   exp(eta + error). Censoring times for future observations are not generated.
#'   Identity-scale fits return their historical working scale. These results
#'   integrate conditional parameter uncertainty and may differ from the
#'   existing plug-in `predict.imr()` point predictions. Their model weights
#'   inherit the approximation described in [posterior_draws()].
#'   For log-time fits, point predictions are medians of the simulated
#'   quantities. A time-scale posterior mean need not exist under an
#'   inverse-gamma variance mixture, so a sample mean is not reported.
#' @return A list of data frames by subgroup, with `id`, `predict`, `lower`
#'   and `upper`. The point prediction is the Monte Carlo mean except for
#'   log-time survival fits, where it is the median.
#' @export
predict.imr_posterior <- function(object, newdata, platform_names = NULL,
                                  covariates = NULL, type = c("mean", "response"),
                                  level = .95, seed = 1L, ...) {
  type <- match.arg(type)
  .imr_interval_level(level)
  seed <- .imr_check_integer_scalar(seed, "seed", min = 0L)
  fit <- object$fit
  inputs <- .imr_prediction_inputs(fit, newdata, platform_names, covariates)
  if (!is.null(inputs$empty)) {
    return(lapply(inputs$empty, function(x) {
      x$lower <- numeric()
      x$upper <- numeric()
      x
    }))
  }
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(seed)
  out <- lapply(seq_along(object$beta), function(g) {
    ids <- inputs$sample_ids[[g]]
    if (!length(ids)) return(data.frame(id = ids, predict = numeric(), lower = numeric(), upper = numeric()))
    X <- .imr_posterior_design(fit, g, inputs$x_test, inputs$cova_test)
    eta <- object$beta[[g]] %*% t(X)
    v <- object$variance[[g]]
    if (fit$type_outcome == "binary") {
      values <- stats::pnorm(eta / sqrt(v))
      if (type == "response") values[] <- stats::rbinom(length(values), 1, values)
    } else {
      values <- eta
      if (type == "response") values <- eta + matrix(stats::rnorm(length(eta)), nrow(eta)) * sqrt(v)
      if (fit$type_outcome == "right.censored" && identical(fit$response_scale, "log")) {
        if (type == "mean") values <- values + v / 2
        values <- exp(values)
      }
    }
    if (any(!is.finite(values))) .imr_abort("Non-finite posterior predictions; inspect tail behavior and variance draws.")
    q <- apply(values, 2L, stats::quantile, probs = c((1-level)/2, (1+level)/2), names = FALSE)
    point <- if (fit$type_outcome == "right.censored" && identical(fit$response_scale, "log")) {
      apply(values, 2L, stats::median)
    } else colMeans(values)
    data.frame(id = ids, predict = point, lower = q[1L, ], upper = q[2L, ], row.names = NULL)
  })
  names(out) <- paste0("model:", fit$model_bitstrings)
  out
}
