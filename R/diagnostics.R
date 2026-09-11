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
  required <- c("schema_version", "control", "model", "preprocessing", "posterior")
  if (!setequal(names(object), required)) {
    missing <- setdiff(required, names(object))
    .imr_abort(if (length(missing)) sprintf("The fitted object is missing `%s`.", missing[1L]) else
      "The fitted object contains unsupported top-level fields.")
  }
  if (!identical(object$schema_version, 2L)) {
    .imr_abort("Unsupported fit schema; use `upgrade_imr_fit()` for a 0.1.x object.")
  }
  control <- object$control
  model <- object$model
  prep <- object$preprocessing
  posterior <- object$posterior
  if (!is.list(control) || !is.list(model) || !is.list(prep) || !is.list(posterior)) {
    .imr_abort("Fit schema sections must be lists.")
  }
  required_control <- c("call", "outcome_type", "response_scale", "method",
    "min_subgroup_size", "priors", "mcmc", "seed")
  required_model <- c("n_platforms", "platform_names", "feature_names",
    "covariate_names", "subgroup_names", "sample_sizes",
    "subgroup_platforms", "platform_subgroups")
  required_prep <- c("input_data", "features", "response", "covariates",
    "feature_center", "feature_scale", "covariate_center", "covariate_scale",
    "formula", "formula_data", "terms", "contrasts", "xlevels", "id")
  required_posterior <- c("inclusion_probabilities", "interaction_means",
    "latent_response_mean", "log_posterior", "selection_draws",
    "interaction_draws")
  section_requirements <- list(control = required_control, model = required_model,
    preprocessing = required_prep, posterior = required_posterior)
  sections <- list(control = control, model = model, preprocessing = prep,
                   posterior = posterior)
  for (section in names(sections)) {
    missing <- setdiff(section_requirements[[section]], names(sections[[section]]))
    if (length(missing)) {
      .imr_abort(sprintf("The `%s` section is missing `%s`.", section, missing[[1L]]))
    }
  }
  if (length(control$outcome_type) != 1L || is.na(control$outcome_type) ||
      !control$outcome_type %in% c("right.censored", "binary", "continuous") ||
      length(control$method) != 1L || is.na(control$method) ||
      !control$method %in% c("imr", "bms")) {
    .imr_abort("The fit has an invalid outcome type or method.")
  }
  if (!is.list(control$priors) ||
      !all(c("nu", "molecular_scale", "forced_scale", "residual", "interaction") %in%
           names(control$priors)) || !is.list(control$mcmc)) {
    .imr_abort("The fit has incomplete prior or MCMC controls.")
  }
  draws <- control$mcmc$draws
  burnin <- control$mcmc$burnin
  if (!.imr_is_integerish(c(draws, burnin)) || length(draws) != 1L ||
      length(burnin) != 1L || draws < 1L || burnin < 0L ||
      as.double(draws) + as.double(burnin) > .Machine$integer.max) {
    .imr_abort("`control$mcmc` must contain valid `draws` and `burnin` counts.")
  }
  n_platforms <- model$n_platforms
  if (length(n_platforms) != 1L || !.imr_is_integerish(n_platforms) || n_platforms < 1L ||
      length(model$platform_names) != n_platforms || anyDuplicated(model$platform_names) ||
      length(model$feature_names) != n_platforms ||
      length(model$platform_subgroups) != n_platforms) {
    .imr_abort("The fit has inconsistent platform metadata.")
  }
  n_subgroups <- length(model$subgroup_names)
  if (!n_subgroups || length(model$sample_sizes) != n_subgroups ||
      length(model$subgroup_platforms) != n_subgroups ||
      !.imr_is_integerish(model$sample_sizes) || any(model$sample_sizes < 1L)) {
    .imr_abort("The fit has inconsistent subgroup metadata.")
  }
  .imr_check_mrf_capacity(model$platform_subgroups)
  for (g in seq_len(n_subgroups)) {
    platforms <- model$subgroup_platforms[[g]]
    if (!.imr_is_integerish(platforms) || anyDuplicated(platforms) ||
        any(platforms < 1L | platforms > n_platforms)) {
      .imr_abort(sprintf("Subgroup %d has invalid platform indices.", g))
    }
    for (l in platforms) {
      if (!g %in% model$platform_subgroups[[l]]) {
        .imr_abort("Subgroup/platform mappings are not reciprocal.")
      }
    }
  }
  if (length(posterior$inclusion_probabilities) != n_platforms ||
      length(posterior$interaction_means) != n_platforms ||
      length(posterior$interaction_draws) != n_platforms) {
    .imr_abort("Posterior platform components have inconsistent lengths.")
  }
  for (l in seq_len(n_platforms)) {
    m <- posterior$inclusion_probabilities[[l]]
    if (!is.matrix(m) || any(!is.finite(m)) || any(m < 0 | m > 1)) {
      .imr_abort(sprintf("Inclusion probabilities for platform %d are invalid.", l))
    }
    models <- model$platform_subgroups[[l]]
    if (!.imr_is_integerish(models) || anyDuplicated(models) ||
        any(models < 1L | models > n_subgroups) ||
        nrow(m) != length(models)) {
      .imr_abort(sprintf("Platform %d has invalid subgroup indices.", l))
    }
    if (ncol(m) != length(model$feature_names[[l]])) {
      .imr_abort(sprintf("Feature names do not match platform %d posterior output.", l))
    }
  }
  if (length(posterior$selection_draws) != draws) {
    .imr_abort("Selection draws do not match `control$mcmc$draws`.")
  }
  for (s in seq_along(posterior$selection_draws)) {
    draw <- posterior$selection_draws[[s]]
    if (!is.list(draw) || length(draw) != n_platforms) {
      .imr_abort(sprintf("Selection draw %d has inconsistent platforms.", s))
    }
    for (l in seq_len(n_platforms)) {
      if (!identical(dim(draw[[l]]), dim(posterior$inclusion_probabilities[[l]])) ||
          any(!draw[[l]] %in% c(0, 1))) {
        .imr_abort(sprintf("Selection draw %d for platform %d is inconsistent.", s, l))
      }
    }
  }
  for (l in seq_len(n_platforms)) {
    n_platform_models <- length(model$platform_subgroups[[l]])
    expected_theta_columns <- choose(n_platform_models, 2L)
    if (!is.matrix(posterior$interaction_means[[l]]) ||
        !identical(dim(posterior$interaction_means[[l]]),
                   c(n_platform_models, n_platform_models)) ||
        any(!is.finite(posterior$interaction_means[[l]]))) {
      .imr_abort(sprintf("Interaction means for platform %d are invalid.", l))
    }
    samples <- posterior$interaction_draws[[l]]
    if (identical(control$method, "bms")) {
      if (!is.null(samples)) {
        .imr_abort("BMS fits must not contain sampled theta interactions.")
      }
    } else if (!is.matrix(samples) || nrow(samples) != draws ||
               ncol(samples) != expected_theta_columns ||
               any(!is.finite(samples))) {
      .imr_abort(sprintf("Interaction draws for platform %d are inconsistent.", l))
    }
  }
  if (!is.numeric(posterior$log_posterior) || any(!is.finite(posterior$log_posterior)) ||
      length(posterior$log_posterior) != draws + burnin) {
    .imr_abort("`log_posterior` must be finite and numeric with one entry per iteration.")
  }
  if (!is.list(prep$features) || length(prep$features) != n_subgroups ||
      !is.list(prep$response) || length(prep$response) != n_subgroups ||
      !is.list(prep$covariates) || length(prep$covariates) != n_subgroups) {
    .imr_abort("Preprocessed native inputs have inconsistent subgroup structure.")
  }
  if (!inherits(prep$input_data, "imr_data")) {
    .imr_abort("`preprocessing$input_data` must be a validated `imr_data` object.")
  }
  validate_imr_data(prep$input_data)
  expected_response_columns <- if (control$outcome_type == "right.censored") 2L else 1L
  for (g in seq_len(n_subgroups)) {
    if (!is.list(prep$features[[g]]) || length(prep$features[[g]]) != n_platforms ||
        !is.matrix(prep$response[[g]]) || nrow(prep$response[[g]]) != model$sample_sizes[[g]] ||
        ncol(prep$response[[g]]) != expected_response_columns ||
        !is.double(prep$response[[g]]) || any(!is.finite(prep$response[[g]])) ||
        !is.matrix(prep$covariates[[g]]) || nrow(prep$covariates[[g]]) != model$sample_sizes[[g]] ||
        ncol(prep$covariates[[g]]) != length(model$covariate_names) ||
        !is.double(prep$covariates[[g]]) || any(!is.finite(prep$covariates[[g]]))) {
      .imr_abort(sprintf("Preprocessed subgroup %d has inconsistent dimensions.", g))
    }
    for (l in seq_len(n_platforms)) {
      x <- prep$features[[g]][[l]]
      expected_rows <- if (l %in% model$subgroup_platforms[[g]]) model$sample_sizes[[g]] else 0L
      if (!is.matrix(x) || !is.double(x) || nrow(x) != expected_rows ||
          ncol(x) != length(model$feature_names[[l]]) || any(!is.finite(x))) {
        .imr_abort(sprintf("Preprocessed subgroup %d platform %d is invalid.", g, l))
      }
    }
  }
  if (!is.list(posterior$latent_response_mean) ||
      length(posterior$latent_response_mean) != n_subgroups ||
      any(vapply(seq_len(n_subgroups), function(g) {
        y <- posterior$latent_response_mean[[g]]
        !is.numeric(y) || length(y) != model$sample_sizes[[g]] || any(!is.finite(y))
      }, logical(1L)))) {
    .imr_abort("Latent-response summaries do not match the fitted subgroups.")
  }
  invisible(TRUE)
}

#' Upgrade a Legacy IMR Fit
#'
#' Convert a structurally complete 0.1.x `imr` object to the named schema used
#' by version 0.2.0. Corrupted or incomplete objects must be refitted.
#'
#' @param object A fitted `imr` object created by IntegMultiReg 0.1.x.
#' @return A validated schema-version-2 `imr` object.
#' @export
upgrade_imr_fit <- function(object) {
  if (!inherits(object, "imr") || !is.list(object)) {
    .imr_abort("`object` must be a legacy `imr` fit.")
  }
  if (identical(object$schema_version, 2L)) {
    validate_imr(object)
    return(object)
  }
  required <- c("gam_mean", "theta_mean", "estimate_latent_y", "log_posterior",
    "gam_sample", "theta_sample", "list_hyperpara", "data1", "data2",
    "type_outcome", "method", "platform_names", "feature_names",
    "covariate_names", "model_bitstrings", "sample_size", "model_platforms",
    "platform_models", "sample_mcmc", "input_data")
  missing <- setdiff(required, names(object))
  if (length(missing)) {
    .imr_abort(sprintf("Legacy fit is missing `%s`; refit the model.", missing[1L]))
  }
  if (length(object$list_hyperpara) != 7L + object$n_platform ||
      length(object$data1) != 9L || length(object$data2) != 7L) {
    .imr_abort("Legacy native payload is incomplete; refit the model.")
  }
  input_data <- object$input_data
  if (!is.null(input_data$type_outcome) && is.null(input_data$outcome_type)) {
    input_data$outcome_type <- input_data$type_outcome
    input_data$type_outcome <- NULL
  }
  h <- object$list_hyperpara
  fit <- list(
    schema_version = 2L,
    control = list(
      call = object$call, outcome_type = object$type_outcome,
      response_scale = object$response_scale,
      method = tolower(object$method), min_subgroup_size = object$ssize,
      priors = list(nu = object$nu, molecular_scale = h[[2L]],
        forced_scale = h[[1L]], residual = c(shape = h[[3L]], rate = h[[4L]]),
        interaction = c(shape = h[[5L]], rate = h[[6L]])),
      mcmc = list(draws = unname(object$sample_mcmc[["total"]]),
                  burnin = unname(object$sample_mcmc[["burnin"]])),
      seed = as.integer(h[[7L]])
    ),
    model = list(
      n_platforms = object$n_platform, platform_names = object$platform_names,
      feature_names = object$feature_names, covariate_names = object$covariate_names,
      subgroup_names = object$model_bitstrings, sample_sizes = object$sample_size,
      subgroup_platforms = object$model_platforms,
      platform_subgroups = object$platform_models
    ),
    preprocessing = list(
      input_data = input_data, features = object$data2[[1L]],
      response = object$data2[[2L]], covariates = object$data2[[3L]],
      feature_center = object$data2[[4L]], feature_scale = object$data2[[5L]],
      covariate_center = object$data2[[6L]], covariate_scale = object$data2[[7L]],
      formula = object$formula, formula_data = object$formula_data,
      terms = object$terms, contrasts = object$contrasts,
      xlevels = object$xlevels, id = object$formula_id %||% "id"
    ),
    posterior = list(
      inclusion_probabilities = object$gam_mean,
      interaction_means = object$theta_mean,
      latent_response_mean = object$estimate_latent_y,
      log_posterior = object$log_posterior,
      selection_draws = object$gam_sample,
      interaction_draws = object$theta_sample
    )
  )
  class(fit) <- "imr"
  validate_imr(fit)
  fit
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

  selection <- lapply(seq_len(object$model$n_platforms), function(l) {
    template <- .imr_mpip(object, l)
    rows <- vector("list", nrow(template) * ncol(template))
    at <- 0L
    for (i in seq_len(nrow(template))) {
      for (j in seq_len(ncol(template))) {
        at <- at + 1L
        draws <- vapply(
          object$posterior$selection_draws,
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
  names(selection) <- object$model$platform_names

  theta <- lapply(seq_len(object$model$n_platforms), function(l) {
    samples <- object$posterior$interaction_draws[[l]]
    if (is.null(samples) || ncol(samples) == 0L) {
      return(data.frame(
        subgroup1 = character(), subgroup2 = character(),
        mean = numeric(), sd = numeric(), lower = numeric(),
        median = numeric(), upper = numeric()
      ))
    }
    subgroup_names <- object$model$subgroup_names[object$model$platform_subgroups[[l]]]
    pairs <- do.call(cbind, lapply(seq.int(2L, length(subgroup_names)),
      function(i) rbind(seq_len(i - 1L), i)))
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
  names(theta) <- object$model$platform_names

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
    identical(fit$control$outcome_type, reference$control$outcome_type) &&
      identical(fit$model$platform_names, reference$model$platform_names) &&
      identical(fit$model$feature_names, reference$model$feature_names) &&
      identical(fit$model$subgroup_names, reference$model$subgroup_names)
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
      fit$posterior$inclusion_probabilities,
      function(m) sum(apply(m, 2L, max) > threshold), integer(1L)
    ))
    data.frame(
      fit = fit_names[i], outcome = fit$control$outcome_type,
      method = fit$control$method,
      platforms = fit$model$n_platforms,
      subgroups = length(fit$model$subgroup_names),
      retained_draws = fit$control$mcmc$draws,
      selected_features = selected,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, rows)
}
