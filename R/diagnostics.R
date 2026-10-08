#' Validate an IMR Object
#'
#' Checks object schema, data mappings, chain dimensions and numerical ranges.
#' It does not establish statistical convergence.
#' @param object An `imr` fit produced by this version of [imr()].
#' @return `TRUE`, invisibly; malformed or older selection-only objects fail.
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
#' validate_imr_object(fit)
#' @export
validate_imr_object <- function(object) {
  if (inherits(object, "imr_posterior")) {
    return(.imr_check_regression_posterior(object))
  }
  .imr_validate_fit(object, check_values = TRUE)
}

.imr_validate_fit <- function(object, check_values = FALSE) {
  if (!inherits(object, "imr") || !is.list(object) || !identical(object$schema_version, 4L)) {
    .imr_abort("An imr fit (schema 4) is required. Earlier selection-only fits need refitting; use their archived package for historical inspection.")
  }
  if (anyDuplicated(names(object)) || !setequal(names(object), c("schema_version", "control", "model", "preprocessing", "posterior", "diagnostics"))) {
    .imr_abort("The fit has missing or unsupported top-level fields.")
  }
  control <- object$control
  model <- object$model
  prep <- object$preprocessing
  post <- object$posterior
  if (!is.list(control) || !is.list(model) || !is.list(prep) || !is.list(post) ||
    length(control$inference) != 1L || !control$inference %in% c(
    "joint_pmom_mrf", "variance_marginal_pmom_mrf",
    "coefficients_marginal_pmom_mrf", "coefficients_and_variance_marginal_pmom_mrf"
  )) {
    .imr_abort("Invalid fit sections.")
  }
  if (length(control$outcome_type) != 1L || !control$outcome_type %in% c("continuous", "binary", "right.censored") ||
    length(control$model_variant) != 1L || !control$model_variant %in% c("imr", "bms")) {
    .imr_abort("Invalid outcome or model variant.")
  }
  marginalize <- .imr_marginalize(control)
  if (length(marginalize) != 1L || !marginalize %in% c("none", "variance", "coefficients", "coefficients_and_variance")) {
    .imr_abort("Invalid marginalization setting.")
  }
  expected_inference <- if (marginalize == "none") "joint_pmom_mrf" else paste0(marginalize, "_marginal_pmom_mrf")
  if (!identical(control$inference, expected_inference)) .imr_abort("The recorded target and marginalization setting disagree.")
  expected_scale <- switch(control$outcome_type,
    continuous = "identity",
    binary = "probit",
    right.censored = c("log", "identity")
  )
  if (length(control$response_scale) != 1L || !control$response_scale %in% expected_scale) {
    .imr_abort("Invalid fitted response scale.")
  }
  mcmc <- .imr_validate_specification(control$mcmc, imr_mcmc, "imr_mcmc")
  .imr_validate_specification(control$priors, imr_priors, "imr_priors")
  groups <- model$subgroup_names
  platforms <- model$platform_names
  if (!is.character(groups) || !length(groups) || anyNA(groups) || anyDuplicated(groups) ||
    !is.character(platforms) || !length(platforms) || anyNA(platforms) || anyDuplicated(platforms) ||
    !identical(model$n_platforms, as.integer(length(platforms))) || length(model$feature_names) != length(platforms) ||
    length(model$subgroup_platforms) != length(groups) || length(model$platform_subgroups) != length(platforms) ||
    length(model$sample_sizes) != length(groups) || any(model$sample_sizes < 1L)) {
    .imr_abort("Inconsistent subgroup/platform metadata.")
  }
  regression <- .imr_marginalize(control) != "coefficients_and_variance"
  if (inherits(object, "imr_selection") == regression) {
    .imr_abort("The fitted class and marginalization setting disagree.")
  }
  if (!regression) {
    for (name in c("selection_mean", "interaction_mean")) {
      if (!is.list(post[[name]]) || !identical(names(post[[name]]), platforms)) {
        .imr_abort("Invalid retained Laplace means.")
      }
    }
    if (!is.list(post$latent_mean) || !identical(names(post$latent_mean), groups)) {
      .imr_abort("Invalid retained Laplace response means.")
    }
    for (g in seq_along(groups)) {
      value <- post$latent_mean[[g]]
      if (!is.numeric(value) || length(value) != model$sample_sizes[g] ||
        (check_values && any(!is.finite(value)))) {
        .imr_abort("Invalid retained Laplace response means.")
      }
    }
    for (l in seq_along(platforms)) {
      members <- length(model$platform_subgroups[[l]])
      for (name in c("selection_mean", "interaction_mean")) {
        value <- post[[name]][[l]]
        columns <- if (name == "selection_mean") length(model$feature_names[[l]]) else members
        if (!is.numeric(value) || !identical(dim(value), as.integer(c(members, columns))) ||
          (check_values && any(!is.finite(value)))) {
          .imr_abort("Invalid retained Laplace mean dimensions or values.")
        }
      }
      if (check_values && any(post$selection_mean[[l]] < 0 | post$selection_mean[[l]] > 1)) {
        .imr_abort("Retained selection means must be probabilities.")
      }
    }
  }
  if ((regression && (!is.list(post$coefficients) || !identical(names(post$coefficients), groups))) ||
    !is.list(post$interaction) || !identical(names(post$interaction), platforms) ||
    length(prep$features) != length(groups) || length(prep$covariates) != length(groups) ||
    length(prep$response) != length(groups) || !identical(names(prep$subject_ids), groups)) {
    .imr_abort("Inconsistent posterior or preprocessing groups.")
  }
  for (seed_name in c("chain_seeds", "initial_seeds")) {
    seeds <- control[[seed_name]]
    if (length(seeds) != mcmc$chains || !.imr_is_integerish(seeds) || any(seeds < 0)) {
      .imr_abort("Invalid recorded chain seeds.")
    }
  }
  if (length(control$priors$nu) != length(platforms) || is.null(control$priors$residual)) {
    .imr_abort("Fitted prior settings must be resolved for every platform and outcome.")
  }
  check_array <- function(x, n, label, positive = FALSE) {
    if (!is.numeric(x) || !identical(dim(x), c(mcmc$draws, mcmc$chains, as.integer(n)))) {
      .imr_abort(sprintf("Invalid %s draw dimensions.", label))
    }
    if (check_values && (any(!is.finite(x)) || (positive && any(x <= 0)))) {
      .imr_abort(sprintf("Invalid %s draw values.", label))
    }
  }
  forced <- 1L + length(model$covariate_names)
  for (g in seq_along(groups)) {
    available <- model$subgroup_platforms[[g]]
    if (!.imr_is_integerish(available) || any(!available %in% seq_along(platforms)) || anyDuplicated(available)) {
      .imr_abort("Invalid subgroup-to-platform mapping.")
    }
    terms <- c(
      "(Intercept)", if (forced > 1L) paste0("clinical:", model$covariate_names),
      unlist(lapply(available, function(l) paste0(platforms[l], ":", model$feature_names[[l]])))
    )
    if (regression) {
      check_array(post$coefficients[[g]], length(terms), "coefficient")
      if (!identical(dimnames(post$coefficients[[g]])[[3L]], terms)) .imr_abort("Coefficient names do not match the design.")
      if (check_values && any(post$coefficients[[g]][, , seq_len(forced), drop = FALSE] == 0)) {
        .imr_abort("A forced coefficient cannot be an exclusion zero.")
      }
    }
    if (length(prep$subject_ids[[g]]) != model$sample_sizes[g] || anyDuplicated(prep$subject_ids[[g]]) ||
      nrow(prep$response[[g]]) != model$sample_sizes[g]) {
      .imr_abort("Inconsistent training subjects.")
    }
  }
  for (l in seq_along(platforms)) {
    expected <- which(vapply(model$subgroup_platforms, function(x) l %in% x, TRUE))
    if (!identical(as.integer(model$platform_subgroups[[l]]), as.integer(expected))) {
      .imr_abort("Platform and subgroup mappings disagree.")
    }
    count <- if (control$model_variant == "bms") 0L else choose(length(expected), 2L)
    check_array(post$interaction[[l]], count, "interaction", positive = TRUE)
    if (!regression) {
      check_array(post$selection[[l]], length(expected) * length(model$feature_names[[l]]), "selection")
      if (check_values && any(!post$selection[[l]] %in% 0:1)) .imr_abort("Selection draws must be zero or one.")
    }
  }
  if (regression) {
    check_array(post$variance, length(groups), "variance", positive = TRUE)
    if (!identical(dimnames(post$variance)[[3L]], groups)) {
      .imr_abort("Variance parameter labels do not match fitted subgroups.")
    }
  } else if (!is.null(post$coefficients) || !is.null(post$variance) ||
    !is.list(post$selection) || !identical(names(post$selection), platforms)) {
    .imr_abort("A Laplace selection fit must store selection draws without fabricated regression draws.")
  }
  if (!is.matrix(post$log_density) || !identical(dim(post$log_density), c(mcmc$draws, mcmc$chains)) ||
    (check_values && any(!is.finite(post$log_density)))) {
    .imr_abort("Invalid target-density trace.")
  }
  if (!is.null(post$latent)) {
    if (!is.list(post$latent) || !identical(names(post$latent), groups)) .imr_abort("Invalid latent-response groups.")
    for (g in seq_along(groups)) {
      x <- post$latent[[g]]
      if (control$outcome_type == "continuous") {
        if (!is.null(x)) .imr_abort("Continuous responses are observed, not augmented.")
      } else {
        check_array(x, model$sample_sizes[g], "latent response")
        if (check_values) {
          for (i in seq_len(model$sample_sizes[g])) {
            y <- prep$response[[g]][i, 1L]
            if (control$outcome_type == "binary") {
              if (any(if (y == 1) x[, , i] < 0 else x[, , i] > 0)) .imr_abort("Latent binary draws have the wrong sign.")
            } else if (prep$response[[g]][i, 2L] == 1) {
              if (any(x[, , i] != y)) .imr_abort("Observed events must remain fixed in latent draws.")
            } else if (any(x[, , i] < y)) .imr_abort("Latent survival draws fall below a censoring bound.")
          }
        }
      }
    }
  } else if (mcmc$keep_latent) .imr_abort("Requested latent storage is missing.")
  invisible(TRUE)
}

#' MCMC Convergence and Monte Carlo Precision Diagnostics
#'
#' Reports rank-normalized split/folded R-hat, bulk/tail effective sample sizes
#' and mean Monte Carlo standard errors using the posterior package. A fit's
#' cached diagnostics are returned when available; otherwise stored draws are
#' inspected without running MCMC. Large diagnostic tasks use the fit's
#' `imr_mcmc(workers = ...)` limit; short tasks run serially.
#'
#' @section Diagnostic quantities:
#' Each chain is split into equal beginning/end halves; an odd middle draw is
#' discarded. On rank-normalized split draws of length n, with within-chain
#' variance W and between-chain variance B, the basic statistic is
#' \deqn{\widehat R=\sqrt{\frac{(n-1)W/n+B/n}{W}}.}{R-hat = sqrt(((n-1)*W/n + B/n) / W).}
#' The reported value is the maximum of this rank statistic and its folded
#' counterpart. Bulk ESS uses rank-normalized draws; tail ESS assesses the
#' lower and upper tails. Mean MCSE estimates Monte Carlo error in the mean,
#' separately from the posterior standard deviation returned by [summary()].
#' Constant chains remain undefined, because these statistics cannot distinguish
#' a fixed parameter from a chain that did not explore it.
#' @param object An `imr` fit or an `imr_posterior` object. Conditional-model
#'   diagnostics for the latter are separate from the source selection chains.
#' @param type `"parameters"` returns R-hat, ESS and MCSE. `"sampler"`
#'   returns a common table of update methods and available Metropolis counts.
#'   Gibbs, integrated, fixed and observed quantities have no acceptance rate.
#' @return A data frame naming parameter family, subgroup/platform and parameter,
#'   with `rhat`, `ess_bulk`, `ess_tail`, `mcse_mean` and diagnostic status.
#'   Constant samples remain undefined. Observed, unaugmented event responses
#'   have status `observed`, not a convergence assessment.
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
#' head(mcmc_diagnostics(fit))
#' @export
mcmc_diagnostics <- function(object, type = c("parameters", "sampler")) {
  type <- match.arg(type)
  if (inherits(object, "imr_posterior")) {
    .imr_check_regression_posterior(object)
    if (type == "sampler") {
      return(object$control$sampler_diagnostics)
    }
    if (is.null(object$diagnostics)) {
      .imr_abort("Conditional-chain diagnostics were not retained; use imr_mcmc(diagnostics = TRUE) when sampling. Mixture rows cannot replace the conditional chains.")
    }
    return(object$diagnostics)
  }
  .imr_check_fit(object)
  if (type == "sampler") {
    return(object$control$acceptance)
  }
  object$diagnostics %||% .imr_compute_diagnostics(object)
}

# Parameter blocks are independent diagnostic tasks. Workers receive only the
# arrays they inspect, not another copy of the complete fitted object.
.imr_compute_diagnostics <- function(object) {
  blocks <- .imr_parameter_blocks(object, "all")
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    observed <- identical(block$family, "latent") && object$control$outcome_type == "right.censored"
    blocks[[i]]$observed <- if (observed) {
      object$preprocessing$response[[match(block$group, object$model$subgroup_names)]][, 2L] == 1
    } else {
      rep(FALSE, dim(block$draws)[3L])
    }
  }
  .imr_diagnose_blocks(blocks, object$control$mcmc$workers)
}

.imr_diagnose_blocks <- function(blocks, workers = 1L) {
  # Short chains rarely recover the cost of starting and loading workers.
  work <- sum(vapply(blocks, function(block) {
    prod(dim(block$draws)[1:2]) * sum(!block$observed)
  }, 0))
  workers <- if (work >= 1e6) workers else 1L
  result <- .imr_map_tasks(blocks, .imr_diagnostic_block, workers)
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  out
}

.imr_diagnostic_block <- function(block) {
  x <- block$draws
  if (!dim(x)[3L]) {
    return(NULL)
  }
  do.call(rbind, lapply(seq_len(dim(x)[3L]), function(j) {
    values <- matrix(x[, , j], nrow = dim(x)[1L], ncol = dim(x)[2L])
    constant <- apply(values, 2L, function(v) all(v == v[1L]))
    status <- if (block$observed[j]) "observed" else if (length(unique(as.vector(values))) == 1L) "constant" else if (ncol(values) < 2L) "single_chain" else if (any(constant)) "constant_in_chain" else "computed"
    computed <- status == "computed"
    data.frame(
      family = block$family, group = block$group,
      parameter = dimnames(x)[[3L]][j],
      rhat = if (computed) posterior::rhat(values) else NA_real_,
      ess_bulk = if (computed) posterior::ess_bulk(values) else NA_real_,
      ess_tail = if (computed) posterior::ess_tail(values) else NA_real_,
      mcse_mean = if (computed) posterior::mcse_mean(values) else NA_real_,
      status = status, stringsAsFactors = FALSE
    )
  }))
}

.imr_warn_diagnostics <- function(x) {
  bad <- sum(x$rhat > 1.01, na.rm = TRUE)
  stuck <- sum(x$status == "constant_in_chain")
  if (bad || stuck) {
    warning(structure(
      list(message = sprintf(
        "Chain diagnostics: %d parameter(s) have R-hat above 1.01; %d vary across samples but are constant within at least one chain. Inspect mcmc_diagnostics() and increase/reassess MCMC before inference.", bad, stuck
      ), call = NULL),
      class = c("imr_mcmc_warning", "warning", "condition")
    ))
  }
  invisible(NULL)
}
