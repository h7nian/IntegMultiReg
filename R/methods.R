.imr_mpip <- function(object, platform) {
  if (!is.null(object$posterior$selection_mean)) {
    return(object$posterior$selection_mean[[platform]])
  }
  model <- object$model
  members <- model$platform_subgroups[[platform]]
  features <- model$feature_names[[platform]]
  result <- matrix(0, length(members), length(features),
    dimnames = list(model$subgroup_names[members], features)
  )
  for (i in seq_along(members)) {
    g <- members[i]
    columns <- .imr_feature_columns(model, g, platform)
    for (j in seq_along(features)) result[i, j] <- mean(object$posterior$coefficients[[g]][, , columns[j]] != 0)
  }
  result
}

#' Extract Posterior Inclusion Probabilities
#'
#' Computes the fraction of retained draws in which a molecular feature is
#' included, separately for each platform and availability subgroup. Fits with
#' regression draws identify inclusion by nonzero coefficients; Laplace fits
#' use their stored selection indicators.
#' \deqn{\widehat\pi_{lsj}=\frac{1}{BC}\sum_{c=1}^{C}\sum_{b=1}^{B}I(\beta_{lsj}^{(b,c)}\ne0).}{mPIP_lsj = fraction of all retained chain draws in which beta_lsj is nonzero.}
#' This is a selection probability; [coef()] extracts regression effects.
#' @param object An `imr` fit.
#' @return A named list of subgroup-by-feature matrices, one per platform.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' inclusion_probabilities(fit)
#' @export
inclusion_probabilities <- function(object) {
  .imr_check_fit(object)
  stats::setNames(lapply(seq_len(object$model$n_platforms), function(l) .imr_mpip(object, l)), object$model$platform_names)
}

#' Extract Regression Coefficients
#'
#' Returns posterior means, including exclusion zeros in model averaging.
#' \deqn{\widehat\beta_j=\frac{1}{BC}\sum_{c,b}\beta_j^{(b,c)}.}{Posterior coefficient mean = average over retained iterations and chains.}
#' Coefficients use the fitted subgroup predictor scales; log-time survival
#' fits describe log time. Extraction performs no sampling.
#' @param object An `imr` fit.
#' @param ... Unused; unsupported arguments fail.
#' @return A named list of coefficient vectors by availability subgroup.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' coef(fit)
#' @export
coef.imr <- function(object, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(object)
  if (is.null(object$posterior$coefficients)) .imr_abort("This selection fit has no stored regression coefficients; use conditional posterior sampling before coefficient extraction.")
  lapply(object$posterior$coefficients, function(x) {
    values <- vapply(seq_len(dim(x)[3L]), function(j) mean(x[, , j]), 0)
    stats::setNames(values, dimnames(x)[[3L]])
  })
}

.imr_feature_ranking <- function(m, select) {
  if (!nrow(m)) {
    return(data.frame(feature = character(), max_mpip = numeric(), subgroup = character()))
  }
  maximum <- apply(m, 2L, max)
  order <- select(maximum)
  data.frame(
    feature = colnames(m)[order], max_mpip = round(maximum[order], 3),
    subgroup = rownames(m)[apply(m, 2L, which.max)][order], row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Print an IMR Fit
#'
#' Reports the model, stored chain budget, feature ranking and available MCMC
#' diagnostics. Ranking uses the maximum marginal inclusion probability over
#' subgroups, which is not the probability of selection in at least one subgroup.
#' @param x An `imr` fit.
#' @param threshold Inclusion-probability cutoff for selected-feature counts.
#' @param rank Display ranked features, including low-probability features.
#' @param top Maximum ranked features per platform.
#' @param ... Unused arguments are rejected.
#' @return `x`, invisibly.
#' @export
print.imr <- function(x, threshold = .5, rank = FALSE, top = 5L, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(x)
  threshold <- .imr_check_threshold(threshold)
  .imr_check_flag(rank, "rank")
  top <- .imr_check_integer_scalar(top, "top", min = 1)
  control <- x$control
  model <- x$model
  cat(toupper(control$model_variant), " ", .imr_sampler_description(control), "\n", sep = "")
  cat("Approximation:", if (.imr_marginalize(control) == "coefficients_and_variance") {
    "Laplace with variance-matched Gaussian moments"
  } else {
    "none"
  }, "\n")
  .imr_print_effective(control$effective)
  cat("Outcome:", control$outcome_type, "(", control$response_scale, "scale )\n")
  cat(sprintf(
    "%d chains; %d retained draws per chain after %d burn-in updates; thinning %d.\n",
    control$mcmc$chains, control$mcmc$draws, control$mcmc$burnin, control$mcmc$thin
  ))
  cat("Subgroups:", paste(paste0(model$subgroup_names, " (n=", model$sample_sizes, ")"), collapse = ", "), "\n")
  for (l in seq_len(model$n_platforms)) {
    m <- .imr_mpip(x, l)
    selected <- if (nrow(m)) sum(apply(m, 2, max) > threshold) else 0L
    cat(sprintf("  %s: %d/%d features above %.3f in a subgroup\n", model$platform_names[l], selected, ncol(m), threshold))
    if (rank && nrow(m)) print(.imr_feature_ranking(m, function(p) utils::head(order(p, decreasing = TRUE), top)), row.names = FALSE)
  }
  if (is.null(x$diagnostics)) cat("MCMC diagnostics have not been computed; use mcmc_diagnostics().\n") else .imr_print_diagnostics(x$diagnostics)
  invisible(x)
}

.imr_print_effective <- function(effective) {
  if (is.null(effective)) {
    return(invisible(NULL))
  }
  fields <- intersect(c("theta_step", "swap_rate", "variance_step"), names(effective$mcmc))
  settings <- c(effective$mcmc[fields], effective$numerical)
  if (!length(settings)) {
    return(invisible(NULL))
  }
  values <- vapply(settings, function(x) {
    if (is.null(x)) {
      return("automatic")
    }
    value <- format(x, trim = TRUE, digits = 4)
    if (length(x) > 1L && !is.null(names(x))) value <- paste0(names(x), "=", value)
    paste(value, collapse = ",")
  }, "")
  cat("Computation:", paste(paste0(names(values), "=", values), collapse = "; "), "\n")
  invisible(NULL)
}

.imr_print_diagnostics <- function(x) {
  defined <- !is.na(x$rhat)
  cat(
    "Rank R-hat:", if (any(defined)) format(max(x$rhat[defined]), digits = 4) else "undefined", "maximum;",
    sum(x$status == "constant"), "constant parameter(s) remain unassessed.\n"
  )
  if (any(x$status == "constant_in_chain")) cat(sum(x$status == "constant_in_chain"), "parameter(s) are constant in at least one chain but vary across samples.\n")
}

#' Summarize an IMR Posterior
#'
#' Summarizes coefficients, residual variances, selection indicators and
#' interactions from the stored draws. Laplace fits contain selection and
#' interaction draws without regression parameters. Stored latent responses can be
#' included with `parm = "latent"` or `"all"`.
#' @param object An `imr` fit.
#' @param parm A parameter family, or `"all"` for all stored families.
#' @param level Equal-tail probability level, strictly between zero and one.
#' @param ... Unused arguments are rejected.
#' @return A `summary.imr` object with a tidy `parameters` table and the fit's
#'   model/chain description. Quantiles use the empirical inverse CDF (type 1),
#'   respecting point masses at zero and binary support. Diagnostics are
#'   computed by the posterior package; undefined values are retained.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' summary(fit, parm = "variance")
#' @export
summary.imr <- function(object, parm = "all", level = .95, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(object)
  level <- .imr_check_level(level)
  .imr_summary(object, parm, level)
}

.imr_summary <- function(object, parm, level, conditional = FALSE) {
  structure(list(
    parameters = .imr_posterior_tables(object, parm, level, include_diagnostics = !conditional), level = level,
    draw_type = if (conditional) "conditional mixture" else if (is.null(object$posterior$coefficients)) "selection chains" else "regression chains",
    model_variant = object$control$model_variant, chains = if (conditional) NA_integer_ else object$control$mcmc$chains,
    draws_per_chain = if (conditional) NA_integer_ else object$control$mcmc$draws,
    output_draws = if (conditional) object$control$mcmc$draws else NULL
  ), class = "summary.imr")
}

#' @rdname summary.imr
#' @param x A `summary.imr` object.
#' @param digits Significant digits for printing.
#' @param max_rows Maximum displayed rows; the returned summary retains all rows.
#' @export
print.summary.imr <- function(x, digits = 4L, max_rows = 30L, ...) {
  .imr_reject_dots(...)
  max_rows <- .imr_check_integer_scalar(max_rows, "max_rows", min = 1)
  cat(sprintf("%s posterior summary (%s): %.1f%% equal-tail intervals\n", toupper(x$model_variant), x$draw_type, 100 * x$level))
  print(utils::head(x$parameters, max_rows), digits = digits, row.names = FALSE)
  if (nrow(x$parameters) > max_rows) cat(nrow(x$parameters) - max_rows, "additional rows in $parameters.\n")
  invisible(x)
}

#' Credible Intervals from an IMR Fit
#'
#' Extracts posterior means, spreads and equal-tail intervals from the stored
#' samples. The default parameter family is regression coefficients; Laplace
#' fits require an explicit selection or interaction family.
#' \deqn{[Q_{(1-L)/2},Q_{(1+L)/2}].}{Interval = empirical quantiles at (1-level)/2 and (1+level)/2.}
#' Quantiles use type 1 to preserve discrete support and exclusion point masses.
#' @param object An `imr` fit.
#' @param parm `"coefficients"` (default), `"variance"`, `"selection"`,
#'   `"interaction"`, `"latent"`, `"all"`, or coefficient names/indices.
#' @param level Equal-tail probability level.
#' @param ... Unused arguments are rejected.
#' @return A data frame with parameter family, group, parameter, posterior mean,
#'   standard deviation and interval endpoints. No additional sampling occurs.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' confint(fit)
#' @export
confint.imr <- function(object, parm = "coefficients", level = .95, ...) {
  .imr_reject_dots(...)
  .imr_check_fit(object)
  level <- .imr_check_level(level)
  .imr_intervals(object, parm, level)
}

.imr_intervals <- function(object, parm, level) {
  families <- c("coefficients", "variance", "selection", "interaction", "latent", "all")
  family <- if (length(parm) == 1L && is.character(parm) && parm %in% families) parm else "coefficients"
  result <- .imr_posterior_tables(object, family, level, include_diagnostics = FALSE)
  if (!identical(parm, family)) {
    terms <- unique(result$parameter)
    if (is.numeric(parm)) {
      if (!.imr_is_integerish(parm) || any(!parm %in% seq_along(terms))) .imr_abort("Invalid coefficient indices.")
      parm <- terms[parm]
    }
    if (!is.character(parm) || any(!parm %in% terms)) .imr_abort("Unknown coefficient name in `parm`.")
    result <- result[result$parameter %in% parm, , drop = FALSE]
  }
  rownames(result) <- NULL
  result
}

#' Compare Descriptive Summaries of Joint IMR Fits
#'
#' Reports data structure, model variant and selected-feature counts at a common
#' threshold. It supplies no predictive ranking, Bayes factor or information
#' criterion. Predictive comparisons should use matched subjects and CV folds.
#' @param ... Two or more fits, or a single named list of fits.
#' @param threshold Common marginal inclusion-probability cutoff.
#' @return One descriptive row per fit. Outcome scales and feature/subgroup
#'   definitions must agree; identical subject cohorts are the caller's responsibility.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' compare_fit_summaries(list(first = fit, repeated = fit))
#' @export
compare_fit_summaries <- function(..., threshold = .5) {
  fits <- list(...)
  if (length(fits) == 1L && !inherits(fits[[1L]], "imr")) fits <- fits[[1L]]
  if (!is.list(fits) || length(fits) < 2L) .imr_abort("Supply at least two joint fits.")
  threshold <- .imr_check_threshold(threshold)
  invisible(lapply(fits, .imr_check_fit))
  reference <- fits[[1L]]
  same <- vapply(fits, function(f) {
    identical(f$control$outcome_type, reference$control$outcome_type) &&
      identical(f$control$response_scale, reference$control$response_scale) &&
      identical(f$model$feature_names, reference$model$feature_names) &&
      identical(f$model$platform_names, reference$model$platform_names) &&
      identical(f$model$subgroup_names, reference$model$subgroup_names)
  }, TRUE)
  if (!all(same)) .imr_abort("Fits must share outcome type/scale, platforms, features and availability subgroups.")
  labels <- names(fits)
  if (is.null(labels) || any(!nzchar(labels))) labels <- paste0("fit", seq_along(fits))
  do.call(rbind, lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    data.frame(
      fit = labels[i], outcome = f$control$outcome_type,
      model_variant = f$control$model_variant, response_scale = f$control$response_scale,
      platforms = f$model$n_platforms, subgroups = length(f$model$subgroup_names),
      chains = f$control$mcmc$chains, draws_per_chain = f$control$mcmc$draws,
      selected_features = sum(vapply(inclusion_probabilities(f), function(m) {
        if (nrow(m)) sum(apply(m, 2, max) > threshold) else 0L
      }, 1L)), stringsAsFactors = FALSE
    )
  }))
}
