#' Validate a Fitted IMR Model
#'
#' Checks the public structure, dimensions and numerical ranges of a fitted IMR
#' object. This is useful after loading a saved fit or before comparing fits.
#'
#' @param object A fitted `"imr"` object.
#' @return `TRUE`, invisibly. Invalid objects fail with an informative error.
#' @export
validate_imr <- function(object) {
  if (!inherits(object, "imr") || !is.list(object)) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  required <- c(
    "gam_mean", "theta_mean", "log_posterior", "gam_sample",
    "model_bitstrings", "sample_size", "platform_names", "feature_names",
    "n_platform", "type_outcome", "method", "sample_mcmc", "theta_sample",
    "platform_models"
  )
  missing <- setdiff(required, names(object))
  if (length(missing) > 0L) {
    .imr_abort(sprintf("The fitted object is missing `%s`.", missing[1L]))
  }
  if (!is.numeric(object$n_platform) || length(object$n_platform) != 1L ||
      is.na(object$n_platform) || !is.finite(object$n_platform) ||
      object$n_platform < 1L || object$n_platform != as.integer(object$n_platform)) {
    .imr_abort("`n_platform` must be one positive integer.")
  }
  if (!is.numeric(object$sample_mcmc) || length(object$sample_mcmc) != 2L ||
      !all(c("total", "burnin") %in% names(object$sample_mcmc)) ||
      any(!is.finite(object$sample_mcmc)) ||
      object$sample_mcmc[["total"]] < 1L ||
      object$sample_mcmc[["burnin"]] < 0L ||
      any(object$sample_mcmc != as.integer(object$sample_mcmc))) {
    .imr_abort("`sample_mcmc` must contain integer `total` and `burnin` counts.")
  }
  if (!is.list(object$platform_models) ||
      length(object$platform_models) != object$n_platform) {
    .imr_abort("`platform_models` has inconsistent platform structure.")
  }
  if (length(object$gam_mean) != object$n_platform ||
      length(object$feature_names) != object$n_platform) {
    .imr_abort("Platform-level components have inconsistent lengths.")
  }
  for (l in seq_len(object$n_platform)) {
    m <- object$gam_mean[[l]]
    if (!is.matrix(m) || any(!is.finite(m)) || any(m < 0 | m > 1)) {
      .imr_abort(sprintf("`gam_mean[[%d]]` must be a finite matrix in [0, 1].", l))
    }
    if (ncol(m) != length(object$feature_names[[l]])) {
      .imr_abort(sprintf("Feature names do not match `gam_mean[[%d]]`.", l))
    }
  }
  retained <- unname(object$sample_mcmc[["total"]])
  if (length(object$gam_sample) != retained) {
    .imr_abort("`gam_sample` does not contain the recorded retained draws.")
  }
  for (s in seq_along(object$gam_sample)) {
    draw <- object$gam_sample[[s]]
    if (!is.list(draw) || length(draw) != object$n_platform) {
      .imr_abort(sprintf("`gam_sample[[%d]]` has inconsistent platforms.", s))
    }
    for (l in seq_len(object$n_platform)) {
      if (!identical(dim(draw[[l]]), dim(object$gam_mean[[l]])) ||
          any(!draw[[l]] %in% c(0, 1))) {
        .imr_abort(sprintf("`gam_sample[[%d]][[%d]]` is inconsistent.", s, l))
      }
    }
  }
  if (!is.list(object$theta_mean) ||
      length(object$theta_mean) != object$n_platform ||
      !is.list(object$theta_sample) ||
      length(object$theta_sample) != object$n_platform) {
    .imr_abort("Theta components have inconsistent platform structure.")
  }
  for (l in seq_len(object$n_platform)) {
    n_platform_models <- length(object$platform_models[[l]])
    expected_theta_columns <- choose(n_platform_models, 2L)
    if (!is.matrix(object$theta_mean[[l]]) ||
        !identical(dim(object$theta_mean[[l]]),
                   c(n_platform_models, n_platform_models)) ||
        any(!is.finite(object$theta_mean[[l]]))) {
      .imr_abort(sprintf("`theta_mean[[%d]]` must be a finite matrix.", l))
    }
    samples <- object$theta_sample[[l]]
    if (identical(object$method, "BMS")) {
      if (!is.null(samples)) {
        .imr_abort("BMS fits must not contain sampled theta interactions.")
      }
    } else if (!is.matrix(samples) || nrow(samples) != retained ||
               ncol(samples) != expected_theta_columns ||
               any(!is.finite(samples))) {
      .imr_abort(sprintf("`theta_sample[[%d]]` is inconsistent.", l))
    }
  }
  if (!is.numeric(object$log_posterior) ||
      any(!is.finite(object$log_posterior))) {
    .imr_abort("`log_posterior` must be finite and numeric.")
  }
  if (length(object$model_bitstrings) != length(object$sample_size)) {
    .imr_abort("Subgroup labels and sample sizes have inconsistent lengths.")
  }
  invisible(TRUE)
}


#' Posterior Uncertainty Summary for an IMR Fit
#'
#' Summarizes retained MCMC draws for variable-selection indicators and MRF
#' interaction parameters. The selection tables report posterior means,
#' posterior standard deviations and equal-tail credible intervals for each
#' platform-feature/subgroup indicator. The theta tables provide the same
#' summaries for each pair of linked availability subgroups.
#'
#' @param object A fitted `"imr"` object.
#' @param level Credible interval level between zero and one (default `0.95`).
#' @param ... Unused; present for future methods.
#' @return An object of class `"posterior_summary.imr"` containing `selection`
#'   and `theta` tables.
#' @export
posterior_summary <- function(object, ...) {
  UseMethod("posterior_summary")
}

#' @rdname posterior_summary
#' @export
posterior_summary.imr <- function(object, level = 0.95, ...) {
  validate_imr(object)
  level <- .imr_check_numeric_vector(
    level, "level", length = 1L, positive = TRUE
  )
  if (level >= 1) .imr_abort("`level` must be less than 1.")
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)

  selection <- lapply(seq_len(object$n_platform), function(l) {
    template <- .imr_mpip(object, l)
    rows <- vector("list", nrow(template) * ncol(template))
    at <- 0L
    for (i in seq_len(nrow(template))) {
      for (j in seq_len(ncol(template))) {
        at <- at + 1L
        draws <- vapply(
          object$gam_sample,
          function(draw) as.numeric(draw[[l]][i, j]), numeric(1L)
        )
        qs <- stats::quantile(draws, probs = probs, names = FALSE, type = 8)
        rows[[at]] <- data.frame(
          subgroup = rownames(template)[i],
          feature = colnames(template)[j],
          mean = mean(draws), sd = stats::sd(draws),
          lower = qs[1L], median = qs[2L], upper = qs[3L],
          row.names = NULL, stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  })
  names(selection) <- object$platform_names

  theta <- lapply(seq_len(object$n_platform), function(l) {
    samples <- object$theta_sample[[l]]
    if (is.null(samples) || ncol(samples) == 0L) {
      return(data.frame(
        subgroup1 = character(), subgroup2 = character(),
        mean = numeric(), sd = numeric(), lower = numeric(),
        median = numeric(), upper = numeric()
      ))
    }
    subgroup_names <- object$model_bitstrings[object$platform_models[[l]]]
    pairs <- utils::combn(seq_along(subgroup_names), 2L)
    rows <- lapply(seq_len(ncol(samples)), function(j) {
      draws <- samples[, j]
      qs <- stats::quantile(draws, probs = probs, names = FALSE, type = 8)
      data.frame(
        subgroup1 = subgroup_names[pairs[1L, j]],
        subgroup2 = subgroup_names[pairs[2L, j]],
        mean = mean(draws), sd = stats::sd(draws),
        lower = qs[1L], median = qs[2L], upper = qs[3L],
        row.names = NULL, stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  names(theta) <- object$platform_names

  out <- list(level = level, selection = selection, theta = theta)
  class(out) <- "posterior_summary.imr"
  out
}

#' @export
print.posterior_summary.imr <- function(x, ...) {
  cat(sprintf("IMR posterior summary (%.1f%% credible intervals)\n", 100 * x$level))
  for (nm in names(x$selection)) {
    cat(sprintf(
      "  %s: %d selection indicators; %d theta interaction(s)\n",
      nm, nrow(x$selection[[nm]]), nrow(x$theta[[nm]])
    ))
  }
  invisible(x)
}


#' Credible Intervals for an IMR Fit
#'
#' Standard `confint()` interface to [posterior_summary()].
#'
#' @param object A fitted `"imr"` object.
#' @param parm Either `"all"`, `"selection"` or `"theta"`.
#' @param level Credible interval level.
#' @param ... Additional arguments passed to [posterior_summary()].
#' @return A list of per-platform credible-interval tables, or a list with both
#'   selection and theta results when `parm = "all"`.
#' @export
confint.imr <- function(object, parm = c("all", "selection", "theta"),
                        level = 0.95, ...) {
  parm <- match.arg(parm)
  out <- posterior_summary(object, level = level, ...)
  if (parm == "all") return(list(selection = out$selection, theta = out$theta))
  out[[parm]]
}


#' Compare Fitted IMR Models
#'
#' Creates a compact descriptive comparison of compatible IMR fits. The table
#' deliberately does not treat raw log-posterior values as likelihood criteria;
#' it reports model structure and the number of features passing a common mPIP
#' threshold.
#'
#' @param ... Fitted `"imr"` objects, or one list of fitted objects.
#' @param threshold Common mPIP threshold used to count selected features.
#' @return A data frame with one row per fit.
#' @export
compare_imr <- function(..., threshold = 0.5) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1L]]) &&
      !inherits(fits[[1L]], "imr")) {
    fits <- fits[[1L]]
  }
  if (length(fits) < 2L) {
    .imr_abort("Supply at least two fitted `imr` objects.")
  }
  threshold <- .imr_check_numeric_vector(
    threshold, "threshold", length = 1L, nonnegative = TRUE
  )
  if (threshold > 1) .imr_abort("`threshold` must be between 0 and 1.")
  invisible(lapply(fits, validate_imr))
  reference <- fits[[1L]]
  compatible <- vapply(fits[-1L], function(fit) {
    identical(fit$type_outcome, reference$type_outcome) &&
      identical(fit$platform_names, reference$platform_names) &&
      identical(fit$feature_names, reference$feature_names) &&
      identical(fit$model_bitstrings, reference$model_bitstrings)
  }, logical(1L))
  if (any(!compatible)) {
    .imr_abort(paste0(
      "Fits must have the same outcome type, platforms, features and ",
      "availability subgroups."
    ))
  }
  fit_names <- names(fits)
  if (is.null(fit_names) || any(!nzchar(fit_names))) {
    fit_names <- paste0("fit", seq_along(fits))
  }
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    selected <- sum(vapply(
      fit$gam_mean,
      function(m) sum(apply(m, 2L, max) > threshold), integer(1L)
    ))
    data.frame(
      fit = fit_names[i], outcome = fit$type_outcome, method = fit$method,
      platforms = fit$n_platform, subgroups = length(fit$model_bitstrings),
      retained_draws = fit$sample_mcmc[["total"]],
      selected_features = selected,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, rows)
}
