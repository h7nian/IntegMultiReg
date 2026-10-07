#' Sample Regression Coefficients and Residual Variances
#'
#' Augments the retained selection models with conditional coefficient and
#' variance draws. The fitted sampler and cross-validation results are unchanged.
#'
#' @section Relation to the fitted selection model:
#' The fitted selection sampler uses a Laplace-based integration over coefficients and
#' variances while exploring selection states. This function adds conditional
#' draws from the coefficient model specified in [imr()]. For one subgroup,
#' the resulting model-averaged approximation has the form
#' \deqn{\widetilde p(b,v\mid\mathcal D)=\sum_m \widehat w_m\,
#'       p(b,v\mid\mathcal D,\gamma=m).}{Approximate posterior(b, v | data) = sum_m empirical_weight_m * conditional_posterior(b, v | data, selection_model=m).}
#' Here \eqn{\mathcal D}{data} is the observed data and \eqn{\widehat w_m}{empirical_weight_m} is the empirical
#' frequency of selection state \eqn{m} in the fitted chain. Whole joint
#' selection draws are resampled before extracting subgroup states, preserving
#' their dependence across subgroups. Inactive molecular coefficients are zero.
#' The conditional Gibbs sampler alternates active coefficients, residual
#' variance and, when needed, the augmented responses. Its coefficient update
#' has density proportional to
#' \deqn{b_j^2\phi(b_j;\mu_j,v/A_{jj}),}{b_j^2 * NormalDensity(b_j; mu_j, v/A_jj),}
#' with matrix
#' \deqn{A=Z^TZ+\mathrm{diag}(1/\tau_j)}{A = transpose(Z) * Z + diag(1/tau_j)}
#' and conditional normal location
#' \deqn{\mu_j=\frac{(Z^Ty^*)_j-\sum_{k\ne j}A_{jk}b_k}{A_{jj}}.}{mu_j = ((transpose(Z) * y*)_j - sum over k != j of A_jk * b_k) / A_jj.}
#'
#' The design, working response and pMOM scales are defined in [imr()].
#'
#' This is an additional computation. Its selection weights inherit the fit's
#' Laplace approximation and finite-chain exploration. Fits made with historical
#' updates or shifted precision boundaries must be refitted under the current
#' model before this sampler is used; converting stored labels does not repair
#' those draws. See the [methods and reproducibility guide](https://h7nian.github.io/IntegMultiReg/method-coverage.html).
#'
#' @section Output size and conditional-chain length:
#' `output_draws` is the number of model-averaged rows returned per subgroup. Suppose
#' one subgroup-model combination is assigned \eqn{r_m} of those rows and
#' `chains` is \eqn{J}. Each of its conditional chains retains
#' \deqn{n_m=\max\{d_{\min},\lceil r_m/J\rceil\},}{Retained iterations per conditional chain: n_m = max(min_draws_per_model_chain, ceiling(returned_rows_m / number_of_chains)).}
#' where \eqn{d_{\min}}{d_min} is `min_draws_per_model_chain`, after discarding `burnin`
#' iterations. Diagnostics use the conditional chains before output subsetting;
#' \eqn{r_m} rows are then sampled without replacement from the pooled draws for the
#' output. Increasing `output_draws` need not lengthen a low-frequency model's chains.
#' Increase `min_draws_per_model_chain`, and if needed `burnin`, to investigate a
#' conditional-chain warning.
#'
#' @section Conditional split R-hat:
#' Each of the \eqn{J} chains is split into two halves of length
#' \deqn{H=\lfloor n_m/2\rfloor.}{H = floor(n_m/2).}
#' An odd-length chain omits its middle draw from this diagnostic. The number
#' of split sequences is
#' \deqn{M=2J.}{M = 2J.}
#'
#' For a scalar parameter \eqn{q}, these sequences have means
#' \eqn{\bar q_c}{mean_q_c}, variances \eqn{s_c^2}
#' and overall mean \eqn{\bar q}{mean_q}. Sums below run over these \eqn{M} sequences:
#' \deqn{W=\frac{1}{M}\sum_c s_c^2,}{W = mean of the M within-sequence variances,}
#' \deqn{B=\frac{H}{M-1}\sum_c(\bar q_c-\bar q)^2,}{B = H * sum_c (mean_q_c - overall_mean_q)^2 / (M - 1),}
#' \deqn{\widehat R=\sqrt{\frac{\frac{H-1}{H}W+\frac{B}{H}}{W}}.}{Split R-hat = sqrt(((H - 1) * W / H + B / H) / W).}
#' `diagnostics` records the maximum over active coefficients and residual
#' variance for each subgroup-model combination. This is classical split
#' R-hat, not the rank-normalized or folded version. An undefined value or a
#' value above 1.05 triggers a warning. It does not diagnose the original
#' selection chain, and passing this check alone does not establish convergence.
#'
#' @param object An `imr` fit with stored data and selection draws.
#' @param output_draws Number of model-averaged draws to return (default `1000`).
#' @param burnin Conditional Gibbs burn-in iterations for each distinct
#'   subgroup selection model and each chain (default `1000`).
#' @param chains Number of conditional chains, at least two (default `2`).
#' @param min_draws_per_model_chain Minimum retained iterations per conditional chain
#'   used for sampling and split R-hat (default `200`). Increase with `burnin`
#'   if the reported conditional diagnostics are poor.
#' @param seed Integer seed. The caller's random-number state is restored.
#' @param latent Retain the augmented response draws. Only binary and
#'   right-censored outcomes have a latent response; for a continuous outcome
#'   the response is observed and nothing is retained. Storing them costs one
#'   numeric per draw per subject, so the default is `FALSE`.
#'
#' @details
#' All active coefficients, including the always-included intercept and clinical
#' effects, follow the original first-order pMOM prior. Inactive molecular
#' coefficients are exactly zero. Molecular and clinical slopes are on the
#' subgroup-standardized predictor scale; responses are not standardized.
#'
#' Whole selection draws are resampled uniformly from the fit's retained
#' selection draws,
#' preserving their empirical joint distribution across subgroups. For each
#' distinct subgroup model, a Gibbs sampler draws coefficients and variance
#' conditional on the observed data. Binary and censored responses are augmented
#' anew; posterior mean latent responses are not treated as observed data.
#' The binary variance uses the same highly concentrated inverse-gamma prior
#' as the fitted model (approximately, not exactly, one).
#'
#' This is a two-stage posterior approximation: model weights inherit the
#' fitted selection sampler's Laplace approximation and finite-chain exploration. The
#' conditional Gibbs draws do not make those weights exact. Classical split
#' R-hat is reported for each conditional model and does not diagnose the
#' original selection chain. Inspect both stages and increase simulation effort
#' before interpreting intervals. Runtime grows with the number of distinct
#' subgroup selection models, not just `output_draws`.
#'
#' Old survival fits without explicit response-scale metadata must be refitted.
#' Log-time fits return coefficients for log time; identity-scale fits are
#' retained only for historical compatibility.
#'
#' @return An `imr_posterior` object containing `beta` (one draws-by-coefficients
#'   matrix per subgroup), `variance`, source `selection_draw_index` indices, conditional
#'   `diagnostics`, and the originating `fit`. With `latent = TRUE` it also
#'   contains `latent`, one draws-by-subject matrix per subgroup holding the
#'   augmented response, paired row by row with `beta`. Coefficient column names
#'   distinguish clinical variables from platform features. Use `summary()`,
#'   `confint()`, `coef()` and `predict()` on this object.
#' @references
#' Chekouo et al. (2017). \doi{10.1111/biom.12587}, Section 3.1 and Web Appendix C.
#' [Read paper](https://academic.oup.com/biometrics/article/73/2/615/7537638) |
#' [Publisher PDF](https://academic.oup.com/biometrics/article-pdf/73/2/615/55973435/biometrics_73_2_615.pdf).
#' @export
#' @examples
#' \donttest{
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            forced_prior_scale = 1, draws = 500, burnin = 250, seed = 1)
#' draws <- sample_regression_posterior(fit, output_draws = 1000, burnin = 1000,
#'                          min_draws_per_model_chain = 1000, seed = 2)
#' summary(draws)
#' }
sample_regression_posterior <- function(object, output_draws = 1000L, burnin = 1000L,
                            chains = 2L, min_draws_per_model_chain = 200L, seed = 1L,
                            latent = FALSE) {
  validate_imr_object(object)
  .imr_require_current_updates(object)
  output_draws <- .imr_check_integer_scalar(output_draws, "output_draws", min = 2L)
  burnin <- .imr_check_integer_scalar(burnin, "burnin", min = 0L)
  chains <- .imr_check_integer_scalar(chains, "chains", min = 2L)
  min_draws_per_model_chain <- .imr_check_integer_scalar(
    min_draws_per_model_chain, "min_draws_per_model_chain", min = 4L)
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
  .imr_check_posterior_data(object)
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(seed)
  selection_draw_index <- sample.int(length(object$posterior$selection_draws), output_draws, replace = TRUE)
  beta <- variance <- diagnostics <- augmented <-
    vector("list", length(object$model$subgroup_names))
  priors <- object$control$priors
  for (g in seq_along(beta)) {
    design <- .imr_posterior_design(object, g)
    masks <- lapply(selection_draw_index, function(s) {
      c(rep(TRUE, 1L + length(object$model$covariate_names)), unlist(lapply(
        object$model$subgroup_platforms[[g]], function(p) {
          object$posterior$selection_draws[[s]][[p]][match(g, object$model$platform_subgroups[[p]]), ] == 1
        }), use.names = FALSE))
    })
    keys <- vapply(masks, function(x) paste(as.integer(x), collapse = ""), "")
    beta[[g]] <- matrix(0, output_draws, ncol(design), dimnames = list(NULL, colnames(design)))
    variance[[g]] <- numeric(output_draws)
    if (keep_latent) augmented[[g]] <- matrix(NA_real_, output_draws,
      nrow(object$preprocessing$response[[g]]),
      dimnames = list(NULL, rownames(object$preprocessing$response[[g]])))
    records <- list()
    for (key in unique(keys)) {
      positions <- which(keys == key)
      active <- masks[[positions[1L]]]
      X <- design[, active, drop = FALSE]
      h <- c(rep(priors$forced_scale, 1L + length(object$model$covariate_names)),
             rep(priors$molecular_scale, ncol(design) - 1L - length(object$model$covariate_names)))[active]
      n <- max(min_draws_per_model_chain, ceiling(length(positions) / chains))
      y <- object$preprocessing$response[[g]][, 1L]
      status <- if (object$control$outcome_type == "right.censored") object$preprocessing$response[[g]][, 2L] else NULL
      samples <- lapply(seq_len(chains), function(chain) {
        .imr_conditional_chain(X, y, h, rep(1, ncol(X)),
          priors$residual[["shape"]], priors$residual[["rate"]],
          draws = n, burnin = burnin,
          initial_beta = rep(if (chain %% 2L) -.5 else .5, ncol(X)),
          initial_variance = 1, outcome_type = object$control$outcome_type,
          status = status, keep_latent = keep_latent)
      })
      rhat <- .imr_split_rhat(samples)
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
      records[[length(records) + 1L]] <- data.frame(
        subgroup = object$model$subgroup_names[g], selection_model = key,
        output_draws = length(positions), draws_per_model_chain = n,
        max_split_rhat = max(rhat), row.names = NULL)
    }
    diagnostics[[g]] <- do.call(rbind, records)
  }
  names(beta) <- names(variance) <- object$model$subgroup_names
  if (keep_latent) names(augmented) <- object$model$subgroup_names
  diagnostics <- do.call(rbind, diagnostics)
  if (any(!is.finite(diagnostics$max_split_rhat) | diagnostics$max_split_rhat > 1.05)) {
    .imr_warn("Conditional split R-hat exceeds 1.05 or is undefined; inspect `diagnostics`, increase `min_draws_per_model_chain` and, if needed, `burnin`, then reassess before interpreting intervals.")
  }
  structure(list(beta = beta, variance = variance,
    latent = if (keep_latent) augmented else NULL,
    selection_draw_index = selection_draw_index,
    diagnostics = diagnostics, fit = object,
    control = list(output_draws = output_draws, burnin = burnin, chains = chains,
                   min_draws_per_model_chain = min_draws_per_model_chain, seed = seed),
    approximation = "Empirical selection weights from the Laplace-based fit; conditional pMOM Gibbs draws."),
    class = "imr_posterior")
}

.imr_posterior_design <- function(fit, g, xx = fit$preprocessing$features,
                                  cc = fit$preprocessing$covariates) {
  n <- nrow(cc[[g]])
  X <- cbind(rep(1, n), cc[[g]])
  nm <- c("(Intercept)", paste0("clinical:", fit$model$covariate_names))
  if (!length(fit$model$covariate_names)) nm <- "(Intercept)"
  for (p in fit$model$subgroup_platforms[[g]]) {
    X <- cbind(X, xx[[g]][[p]])
    nm <- c(nm, paste0(fit$model$platform_names[p], ":", fit$model$feature_names[[p]]))
  }
  colnames(X) <- make.unique(nm)
  X
}

.imr_check_posterior_data <- function(fit) {
  priors <- fit$control$priors
  if (any(!is.finite(c(priors$forced_scale, priors$molecular_scale,
                       priors$residual))) ||
      any(c(priors$forced_scale, priors$molecular_scale, priors$residual) <= 0)) {
    .imr_abort("The fit is missing valid posterior data or hyperparameters.")
  }
  for (g in seq_along(fit$model$subgroup_names)) {
    X <- .imr_posterior_design(fit, g)
    y <- fit$preprocessing$response[[g]]
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

#' Summarize Coefficient and Latent-Response Draws
#'
#' Extract posterior means, spreads and intervals from an `imr_posterior`
#' object without running additional sampling.
#' @section Posterior summaries:
#' For \eqn{S} stored draws of coefficient \eqn{b_j}, `coef()` returns
#' \deqn{\bar b_j=\frac{1}{S}\sum_{d=1}^{S}b_j^{(d)}.}{Posterior mean of coefficient j = sum_d beta_j[d] / S.}
#' `summary()` and `confint()` report this mean, the sample standard deviation,
#' the median and an equal-tail interval:
#' \deqn{[Q_{(1-L)/2}(b_j),Q_{(1+L)/2}(b_j)].}{Equal-tail interval: [quantile(beta_j, (1 - level)/2), quantile(beta_j, (1 + level)/2)].}
#'
#' The nonzero-probability column is
#' \deqn{\widehat{\Pr}(b_j\ne0\mid\mathcal D)=
#'       \frac{1}{S}\sum_{d=1}^{S}I(b_j^{(d)}\ne0).}{Estimated probability that coefficient j is nonzero = sum_d I(beta_j[d] != 0) / S.}
#' Here \eqn{L} is `level`, \eqn{\mathcal D}{data} denotes the observed data,
#' and quantiles use `stats::quantile(type = 7)`. The interval includes the
#' point mass at zero from inactive models. `confint()` reuses `summary()`;
#' neither method runs further MCMC. `print()` reports the stored approximation
#' and maximum conditional R-hat described in [sample_regression_posterior()].
#'
#' With `parm = "latent"`, the same mean, standard deviation and quantiles are
#' computed per subject, without a nonzero-probability column. For log-time
#' survival fits these are log-time summaries: an observed event is constant
#' across draws, whereas a censored response is sampled above its observed
#' log-time bound. For binary fits the latent normal value is below or above
#' zero according to the observed class. These latent summaries are conditional
#' on observed outcomes; they are not predictions for new subjects.
#'
#' @param object,x An `imr_posterior` object returned by [sample_regression_posterior()].
#' @param level Equal-tail credible level, between zero and one.
#' @param parm `"coefficients"` for regression coefficients, `"variance"` for
#'   residual variance, or `"latent"` for a binary/censored augmented response. `"latent"`
#'   requires draws retained with `sample_regression_posterior(latent = TRUE)` and gives one
#'   row per subject, without a `probability_nonzero` column.
#' @param ... Unused arguments are rejected.
#' @return `summary()` and `confint()` return tables for the requested parameter
#'   set, by subgroup.
#'   `coef()` returns posterior mean coefficient vectors. `print()` returns
#'   its argument invisibly. Effects are on subgroup-standardized predictor
#'   scales. Intervals include zero-mass from inactive selection models.
#' @name imr_posterior_methods
NULL

#' @rdname imr_posterior_methods
#' @export
summary.imr_posterior <- function(object, level = .95, parm = "coefficients", ...) {
  .imr_reject_dots(...)
  validate_imr_object(object)
  .imr_check_interval_level(level)
  parm <- .imr_posterior_parm(parm)
  if (parm == "latent") return(.imr_latent_summary(object, level))
  probs <- c((1-level)/2, .5, (1+level)/2)
  if (parm == "variance") return(lapply(object$variance, function(x)
    .imr_draw_interval(x, probs, quantile_type = 7L, parameter = "residual_variance")))
  lapply(object$beta, function(x) do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    out <- .imr_draw_interval(x[, j], probs, quantile_type = 7L, term = colnames(x)[j])
    out$probability_nonzero <- mean(x[, j] != 0)
    out
  })))
}

#' @rdname imr_posterior_methods
#' @export
confint.imr_posterior <- function(object, parm = "coefficients",
                                  level = .95, ...) {
  .imr_reject_dots(...)
  summary(object, level = level, parm = .imr_posterior_parm(parm))
}

# Named explicitly rather than through match.arg(), so that a bad value is
# reported against `parm` as the other methods report their arguments.
.imr_posterior_parm <- function(parm) {
  choices <- c("coefficients", "variance", "latent")
  if (length(parm) == 1L && !is.na(parm) && parm %in% choices) return(parm)
  .imr_abort("`parm` must be 'coefficients', 'variance' or 'latent'.")
}

# One row per subject, summarising the augmented response that the conditional
# chains sampled. There is no `probability_nonzero` column: a latent response is
# continuous and never exactly zero.
.imr_latent_summary <- function(object, level) {
  if (is.null(object$latent)) {
    .imr_abort(if (object$fit$control$outcome_type == "continuous")
      "A continuous outcome has no latent response; its response is observed."
      else "Latent draws were not retained; call `sample_regression_posterior(latent = TRUE)`.")
  }
  probs <- c((1 - level) / 2, .5, (1 + level) / 2)
  lapply(object$latent, function(x) do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    id <- if (is.null(colnames(x))) j else colnames(x)[j]
    .imr_draw_interval(x[, j], probs, quantile_type = 7L, id = id)
  })))
}

#' @rdname imr_posterior_methods
#' @export
coef.imr_posterior <- function(object, ...) {
  .imr_reject_dots(...)
  validate_imr_object(object)
  lapply(object$beta, colMeans)
}

#' @rdname imr_posterior_methods
#' @export
print.imr_posterior <- function(x, ...) {
  .imr_reject_dots(...)
  validate_imr_object(x)
  cat("Regression posterior for", toupper(x$fit$control$model_variant), "fit:",
      length(x$selection_draw_index), "draws;",
      length(x$beta), "availability subgroups\n")
  if (!is.null(x$latent)) cat("Latent response draws retained.\n")
  cat(x$approximation, "\n")
  cat("Maximum conditional split R-hat:", max(x$diagnostics$max_split_rhat), "\n")
  invisible(x)
}

#' Predict Using Coefficient Posterior Draws
#'
#' Compute point summaries and equal-tail intervals for new subjects from
#' stored coefficient and variance draws, optionally including outcome noise.
#' @section Posterior prediction:
#' For a new subject, let \eqn{z} be the design row transformed with the
#' training subgroup's centers, scales and formula encoding. For each stored
#' posterior draw, form
#' \deqn{\eta^{(d)}=z^Tb^{(d)}}{eta[d] = transpose(z) * beta[d]}
#' and use its paired variance \eqn{v^{(d)}}.
#'
#' For a continuous response, `quantity = "conditional_mean"` summarizes \eqn{\eta^{(d)}}{eta[d]}.
#' `quantity = "new_observation"` instead generates
#' \deqn{Y_{\mathrm{new}}^{(d)}\sim N(\eta^{(d)},v^{(d)}).}{Y_new[d] ~ N(eta[d], v[d]).}
#'
#' For a binary response, the mean draws are event probabilities
#' \deqn{p^{(d)}=\Phi\left(\frac{\eta^{(d)}}{\sqrt{v^{(d)}}}\right),}{p[d] = Phi(eta[d] / sqrt(v[d])),}
#' and response draws follow
#' \deqn{Y_{\mathrm{new}}^{(d)}\sim\mathrm{Bernoulli}(p^{(d)}).}{Y_new[d] ~ Bernoulli(p[d]).}
#'
#' For default log-time survival fits, the corresponding time-scale quantities
#' are
#' \deqn{m^{(d)}=\exp\{\eta^{(d)}+v^{(d)}/2\},}{Conditional mean time m[d] = exp(eta[d] + v[d]/2),}
#' \deqn{T_{\mathrm{new}}^{(d)}=\exp\{\eta^{(d)}+\epsilon^{(d)}\},}{New time T[d] = exp(eta[d] + error[d]),}
#' where
#' \deqn{\epsilon^{(d)}\sim N(0,v^{(d)}).}{error[d] ~ N(0, v[d]).}
#'
#' The first expression is the
#' conditional log-normal mean at a parameter draw, not the median survival
#' time at that draw. Identity-scale survival fits use the normal working
#' response without exponentiation.
#'
#' The interval endpoints are empirical equal-tail quantiles (type 7) of the
#' chosen quantity. The point summary is its sample mean for continuous and
#' binary outcomes and its sample median for log-time survival outcomes, whose
#' posterior mean may not exist. Response intervals include future outcome
#' variation; mean intervals describe parameter and model uncertainty. These
#' draw-based summaries inherit the approximation in [sample_regression_posterior()] and
#' can differ from the ranked-model plug-in predictions in [predict.imr()].
#'
#' @param object An `imr_posterior` object.
#' @param newdata,platform_names,covariates As in [predict.imr()].
#' @param quantity `"conditional_mean"` returns uncertainty in the conditional response mean
#'   (event probability for binary data). `"new_observation"` additionally generates
#'   new outcomes, including residual variability. Binary response intervals
#'   summarize 0/1 draws with interpolated quantiles; use `"conditional_mean"` for
#'   event-probability intervals.
#' @param level Equal-tail interval level.
#' @param seed Integer simulation seed; the caller's RNG state is restored.
#' @param ... Unused arguments are rejected.
#' @details Log-time survival fits return time-scale results: for `"conditional_mean"`,
#'   each draw is exp(eta + variance/2); for `"new_observation"`, each draw is
#'   exp(eta + error). Censoring times for future observations are not generated.
#'   Identity-scale fits return their historical working scale. These results
#'   integrate conditional parameter uncertainty and may differ from the
#'   existing plug-in `predict.imr()` point predictions. Their model weights
#'   inherit the approximation described in [sample_regression_posterior()].
#'   For log-time fits, point predictions are medians of the simulated
#'   quantities. A time-scale posterior mean need not exist under an
#'   inverse-gamma variance mixture, so a sample mean is not reported.
#' @return A list of data frames by subgroup, with `id`, `prediction`, `lower`
#'   and `upper`. The point prediction is the Monte Carlo mean except for
#'   log-time survival fits, where it is the median.
#' @export
predict.imr_posterior <- function(object, newdata, platform_names = NULL,
                                  covariates = NULL, quantity = c("conditional_mean", "new_observation"),
                                  level = .95, seed = 1L, ...) {
  .imr_reject_dots(...)
  validate_imr_object(object)
  .imr_require_current_updates(object$fit)
  quantity <- match.arg(quantity)
  .imr_check_interval_level(level)
  seed <- .imr_check_integer_scalar(seed, "seed", min = 0L)
  fit <- object$fit
  inputs <- .imr_prediction_inputs(fit, newdata, platform_names, covariates)
  if (!is.null(inputs$empty)) {
    return(.imr_prediction_result(lapply(inputs$empty, function(x) {
      x$lower <- numeric()
      x$upper <- numeric()
      x
    }), fit))
  }
  rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(rng), add = TRUE)
  set.seed(seed)
  out <- lapply(seq_along(object$beta), function(g) {
    ids <- inputs$sample_ids[[g]]
    if (!length(ids)) return(data.frame(id = ids, prediction = numeric(), lower = numeric(), upper = numeric()))
    X <- .imr_posterior_design(fit, g, inputs$x_test, inputs$cova_test)
    eta <- object$beta[[g]] %*% t(X)
    v <- object$variance[[g]]
    if (fit$control$outcome_type == "binary") {
      values <- stats::pnorm(eta / sqrt(v))
      if (quantity == "new_observation") values[] <- stats::rbinom(length(values), 1, values)
    } else {
      values <- eta
      if (quantity == "new_observation") values <- eta + matrix(stats::rnorm(length(eta)), nrow(eta)) * sqrt(v)
      if (fit$control$outcome_type == "right.censored" && identical(fit$control$response_scale, "log")) {
        if (quantity == "conditional_mean") values <- values + v / 2
        values <- exp(values)
      }
    }
    if (any(!is.finite(values))) .imr_abort("Non-finite posterior predictions; inspect tail behavior and variance draws.")
    q <- apply(values, 2L, stats::quantile, probs = c((1-level)/2, (1+level)/2), names = FALSE)
    point <- if (fit$control$outcome_type == "right.censored" && identical(fit$control$response_scale, "log")) {
      apply(values, 2L, stats::median)
    } else colMeans(values)
    data.frame(id = ids, prediction = point, lower = q[1L, ], upper = q[2L, ], row.names = NULL)
  })
  .imr_prediction_result(out, fit)
}
