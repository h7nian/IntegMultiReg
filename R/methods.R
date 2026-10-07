## S3 methods for objects of class "imr".

## Internal: marginal posterior inclusion probability (mPIP) matrix for one
## platform, with subgroup (row) and feature (column) names attached.
#' @keywords internal
#' @noRd
.imr_mpip <- function(object, platform) {
  m <- object$posterior$inclusion_probabilities[[platform]]
  rn <- object$model$subgroup_names[object$model$platform_subgroups[[platform]]]
  cn <- object$model$feature_names[[platform]]
  if (!is.null(rn) && length(rn) == nrow(m)) rownames(m) <- rn
  if (!is.null(cn) && length(cn) == ncol(m)) colnames(m) <- cn
  m
}

## Internal: generic platform labels used to decode availability bitstrings.
#' @keywords internal
#' @noRd
.imr_platform_codes <- function(n_platform) {
  paste0("P", seq_len(n_platform))
}


## Internal: translate a bitstring into generic platform labels.  The right-most
## bit corresponds to P1, matching the convention used throughout the package.
#' @keywords internal
#' @noRd
.imr_bitstring_codes <- function(bitstring, n_platform) {
  bits <- strsplit(as.character(bitstring), "", fixed = TRUE)[[1]]
  if (length(bits) < n_platform) {
    bits <- c(rep("0", n_platform - length(bits)), bits)
  }
  present <- which(rev(bits) == "1")
  if (!length(present)) return("none")
  paste(.imr_platform_codes(n_platform)[present], collapse = " + ")
}


#' Marginal Posterior Inclusion Probabilities of an IMR Fit
#'
#' Extracts the posterior mean variable-selection probabilities (the marginal
#' posterior inclusion probabilities, mPIP) of a fitted model.
#'
#' @section Statistical definition:
#' For platform \eqn{l}, subgroup \eqn{s} and feature \eqn{j}, the returned
#' entry is the retained-chain average of its binary selection indicator:
#' \deqn{\widehat\pi_{lsj} = \frac{1}{B}\sum_{b=1}^{B}\gamma_{lsj}^{(b)}.}{mPIP_lsj = sum_b gamma_lsj[b] / B.}
#' Here \eqn{B} excludes burn-in. This estimates a marginal inclusion
#' probability under the fitted sampler and its approximations. It is not a
#' regression coefficient; `coef()` on an `imr_posterior` object instead
#' returns coefficient posterior means (see [imr_posterior_methods]). Use [posterior_summary()] for
#' selection-indicator and interaction summaries, or [posterior_draws()] for
#' coefficient uncertainty.
#'
#' @param object A fitted object of class `"imr"`.
#' @param ... Unused; present for S3 compatibility.
#' @return A named list with one matrix per platform.  Rows are the subgroups
#'   containing that platform (labelled by their availability bitstrings) and
#'   columns are the platform features.
#' @seealso [imr()]
#' @export
coef.imr <- function(object, ...) {
  validate_imr(object)
  out <- lapply(seq_len(object$model$n_platforms), function(l) .imr_mpip(object, l))
  names(out) <- object$model$platform_names
  out
}


# Features ranked by their highest inclusion probability across subgroups.
# The two callers differ only in which columns they keep, so `select` receives
# those probabilities and returns the column order to report.
.imr_feature_ranking <- function(m, select) {
  maxp <- apply(m, 2, max)
  ord <- select(maxp)
  data.frame(
    feature = colnames(m)[ord],
    max_mpip = round(maxp[ord], 3),
    subgroup = rownames(m)[apply(m, 2, which.max)][ord],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Print Method for IMR Fits
#'
#' Prints a compact overview of the fit, including a platform key (`P1`,
#' `P2`, ...) that decodes the availability-subgroup bitstrings.  Optionally,
#' it can also print the top-ranked features per platform, ranked by each
#' feature's maximum mPIP over availability subgroups containing that platform.
#'
#' @section Reading the display:
#' Selected-feature counts use the strict threshold and maximum subgroup mPIP
#' defined in [summary.imr()]. With `rank = TRUE`, the display shows the `top`
#' highest-ranking features even if some are below `threshold`; the threshold
#' controls the counts, not the displayed ranking. Printing does not rerun
#' sampling or change the stored inclusion probabilities.
#'
#' @param x A fitted object of class `"imr"`.
#' @param threshold Inclusion-probability threshold used to count selected
#'   features (default `0.5`).
#' @param rank Logical; if `TRUE`, print a short ranked feature table for each
#'   platform (default `FALSE`).
#' @param top Integer; when `rank = TRUE`, the number of top-ranked features to
#'   show per platform (default `5`).
#' @param ... Unused; present for S3 compatibility.
#' @return `x`, invisibly.
#' @export
print.imr <- function(x, threshold = 0.5, rank = FALSE, top = 5, ...) {
  threshold <- .imr_check_threshold(threshold)
  .imr_check_flag(rank, "rank")
  top <- .imr_check_integer_scalar(top, "top", min = 1)
  cat("Integrative Bayesian Multi-Platform Regression (IMR)\n")
  cat("----------------------------------------------------\n")
  if (!is.null(x$control$call)) {
    cat("Call:\n  ")
    print(x$control$call)
  }
  validate_imr(x)
  cat(sprintf("\nOutcome type : %s\n", x$control$outcome_type))
  cat(sprintf("Method       : %s\n", toupper(x$control$method)))
  cat(sprintf("Platforms    : %d (%s)\n", x$model$n_platforms,
              paste(x$model$platform_names, collapse = ", ")))
  cat(sprintf("MCMC         : %d retained draws after %d burn-in\n",
              x$control$mcmc$draws, x$control$mcmc$burnin))

  codes <- .imr_platform_codes(x$model$n_platforms)
  cat("\nPlatform key:\n")
  key <- paste(sprintf("  %s = %s", codes, x$model$platform_names), collapse = "\n")
  cat(key, "\n", sep = "")

  cat("\nAvailability subgroups modelled (bitstring : platforms : size):\n")
  subgroup_codes <- vapply(
    x$model$subgroup_names, .imr_bitstring_codes, character(1),
    n_platform = x$model$n_platforms
  )
  st <- paste(sprintf("  %-*s : %-*s : %d",
                      max(nchar(x$model$subgroup_names)), x$model$subgroup_names,
                      max(nchar(subgroup_codes)), subgroup_codes,
                      as.integer(x$model$sample_sizes)),
              collapse = "\n")
  cat(st, "\n", sep = "")

  cat(sprintf("\nFeatures with mPIP > %.2f (in any subgroup):\n", threshold))
  for (l in seq_len(x$model$n_platforms)) {
    m <- x$posterior$inclusion_probabilities[[l]]
    sel <- if (nrow(m) > 0 && ncol(m) > 0) sum(apply(m, 2, max) > threshold) else 0L
    cat(sprintf("  %-12s : %d of %d\n", x$model$platform_names[l], sel, ncol(m)))
  }
  if (rank) {
    cat(sprintf("\nTop %d ranked features by maximum subgroup mPIP:\n", top))
    for (l in seq_len(x$model$n_platforms)) {
      m <- .imr_mpip(x, l)
      cat(sprintf("  %s\n", x$model$platform_names[l]))
      if (nrow(m) == 0 || ncol(m) == 0) {
        cat("    (no selectable features)\n")
        next
      }
      tab <- .imr_feature_ranking(m, function(maxp)
        utils::head(order(maxp, decreasing = TRUE), top))
      lines <- utils::capture.output(print(tab, row.names = FALSE))
      cat(paste0("    ", lines), sep = "\n")
      cat("\n")
    }
  }
  invisible(x)
}


#' Summarize an IMR Fit
#'
#' Produces a per-platform summary of the selected features (those whose
#' marginal posterior inclusion probability exceeds `threshold` in at least one
#' subgroup), ranked by their maximum inclusion probability.
#'
#' @section Selection and ranking:
#' Let \eqn{\widehat\pi_{lsj}}{mPIP_lsj} denote the subgroup mPIP defined in [coef.imr()].
#' For each platform-feature pair, the ranking score and selected set are
#' \deqn{r_{lj}=\max_{s\in\mathcal S_l}\widehat\pi_{lsj}, \qquad
#'       \mathcal A_l(t)=\{j:r_{lj}>t\},}{r_lj = maximum subgroup mPIP for feature j on platform l; A_l(t) contains features with r_lj > t.}
#' where \eqn{\mathcal S_l}{S_l} contains subgroups with platform \eqn{l}, and
#' \eqn{t} is `threshold`. The inequality is strict. Rows are sorted by
#' \eqn{r_{lj}}; the reported subgroup is the first maximizing row in the fitted
#' order. This maximum is a ranking score, not the posterior probability of
#' selection in at least one subgroup. A common feature is counted once per
#' platform, even if it exceeds the threshold in several subgroups.
#'
#' @param object A fitted object of class `"imr"`.
#' @param threshold Inclusion-probability threshold for selection (default
#'   `0.5`).
#' @param ... Unused; present for S3 compatibility.
#' @return An object of class `"summary.imr"`: a list with the run metadata and,
#'   for each platform, a data frame of selected features with their maximum
#'   mPIP and the subgroup achieving it.
#' @export
summary.imr <- function(object, threshold = 0.5, ...) {
  threshold <- .imr_check_threshold(threshold)
  validate_imr(object)
  selected <- vector("list", object$model$n_platforms)
  names(selected) <- object$model$platform_names
  for (l in seq_len(object$model$n_platforms)) {
    m <- .imr_mpip(object, l)
    if (nrow(m) == 0 || ncol(m) == 0) {
      selected[[l]] <- data.frame(feature = character(0), max_mpip = numeric(0),
                                  subgroup = character(0))
      next
    }
    selected[[l]] <- .imr_feature_ranking(m, function(maxp) {
      keep <- which(maxp > threshold)
      keep[order(maxp[keep], decreasing = TRUE)]
    })
  }
  out <- list(
    call = object$control$call,
    outcome_type = object$control$outcome_type,
    method = object$control$method,
    threshold = threshold,
    sample_sizes = object$model$sample_sizes,
    subgroup_names = object$model$subgroup_names,
    platform_names = object$model$platform_names,
    selected = selected
  )
  class(out) <- "summary.imr"
  out
}

#' @rdname summary.imr
#' @param x A `"summary.imr"` object.
#' @export
print.summary.imr <- function(x, ...) {
  cat("Integrative Bayesian Multi-Platform Regression (IMR) -- summary\n")
  cat("--------------------------------------------------------------\n")
  cat(sprintf("Outcome type : %s   Method: %s\n", x$outcome_type, toupper(x$method)))
  cat(sprintf("Selection threshold (mPIP) : %.2f\n\n", x$threshold))
  for (l in seq_along(x$selected)) {
    df <- x$selected[[l]]
    cat(sprintf("Platform '%s': %d selected feature(s)\n",
                x$platform_names[l], nrow(df)))
    if (nrow(df) > 0) {
      print(df, row.names = FALSE)
    }
    cat("\n")
  }
  invisible(x)
}
