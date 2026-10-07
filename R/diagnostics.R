#' Check the Structure of an IMR Object
#'
#' Checks the stored structure, dimensions and numerical ranges of an IMR fit
#' or regression-posterior draw object. This is useful after loading a saved fit or before comparing fits.
#'
#' @section Validation scope:
#' Validation checks the fitted-object schema, data dimensions, subgroup/platform
#' mappings, parameter ranges and consistency of retained-draw counts. It does
#' not rerun MCMC or establish posterior convergence. For statistical inspection,
#' use [plot.imr()] and [selection_summary()] for the selection stage and
#' [sample_regression_posterior()] for conditional coefficient sampling and its diagnostics.
#' The statistical model and computational conventions are described in [imr()].
#'
#' @param object An `imr` fit or `imr_posterior` regression-draw object.
#' @return `TRUE`, invisibly. Invalid objects fail with an informative error.
#' @examples
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            forced_prior_scale = 1, draws = 500, burnin = 250, seed = 1)
#'
#' validate_imr_object(fit)
#'
#' # A damaged object is rejected rather than used.
#' damaged <- fit
#' damaged$posterior <- NULL
#' try(validate_imr_object(damaged))
#' @export
validate_imr_object <- function(object) {
  if (inherits(object, "imr_posterior")) return(.imr_validate_regression_draws(object))
  if (!inherits(object, "imr") || !is.list(object)) {
    .imr_abort("`object` must be an `imr` object returned by `imr()`.")
  }
  required <- c("schema_version", "control", "model", "preprocessing", "posterior")
  if (!setequal(names(object), required)) {
    missing <- setdiff(required, names(object))
    .imr_abort(if (length(missing)) sprintf("The fitted object is missing `%s`.", missing[1L]) else
      "The fitted object contains unsupported top-level fields.")
  }
  if (!identical(object$schema_version, 3L)) {
    .imr_abort("Unsupported fit schema; use `upgrade_imr_object()` for a saved object from an earlier version.")
  }
  control <- object$control
  model <- object$model
  prep <- object$preprocessing
  posterior <- object$posterior
  if (!is.list(control) || !is.list(model) || !is.list(prep) || !is.list(posterior)) {
    .imr_abort("Fit schema sections must be lists.")
  }
  required_control <- c("call", "outcome_type", "response_scale", "model_variant", "selection_update",
    "min_subgroup_size", "priors", "mcmc", "seed", "numerical", "standardize")
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
      length(control$model_variant) != 1L || is.na(control$model_variant) ||
      !control$model_variant %in% c("imr", "bms")) {
    .imr_abort("The fit has an invalid outcome type or method.")
  }
  .imr_fit_numerical_control(control)
  if (!is.null(control$standardize)) .imr_check_flag(control$standardize, "standardize")
  if (!is.character(control$selection_update) || length(control$selection_update) != 1L ||
      is.na(control$selection_update) ||
      !control$selection_update %in% c("symmetric_mrf_hastings", "unadjusted_flip_swap"))
    .imr_abort("The fit has an invalid selection-update record.")
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
  if (!is.null(control$laplace_diagnostics)) {
    d <- control$laplace_diagnostics
    measures <- c("calls", "iteration_limit", "nonfinite", "factorization_failures")
    if (!is.data.frame(d) || !identical(names(d), c("subgroup", "stage", measures)) ||
        !identical(d$subgroup, rep(model$subgroup_names, each = 3L)) ||
        !identical(d$stage, rep(c("initial", "selection", "latent"), n_subgroups)) ||
        !.imr_is_integerish(unlist(d[measures], use.names = FALSE)) ||
        any(as.matrix(d[measures]) < 0) || any(d$iteration_limit > d$calls) ||
        any(d$nonfinite > d$calls))
      .imr_abort("The fit has invalid Laplace diagnostics.")
  }
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
    if (identical(control$model_variant, "bms")) {
      if (!is.null(samples)) {
        .imr_abort("BMS fits must not contain sampled theta interactions.")
      }
    } else if (!is.matrix(samples) || nrow(samples) != draws ||
               ncol(samples) != expected_theta_columns ||
               any(!is.finite(samples)) || any(samples <= 0)) {
      .imr_abort(sprintf("Interaction draws for platform %d are inconsistent.", l))
    }
    theta <- posterior$interaction_means[[l]]
    if (!isTRUE(all.equal(theta, t(theta), tolerance = 0)) ||
        any(diag(theta) != 0) ||
        (control$model_variant == "imr" && any(theta[row(theta) != col(theta)] <= 0)) ||
        (control$model_variant == "bms" && any(theta != 0))) {
      .imr_abort(sprintf("Interaction means for platform %d are invalid.", l))
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

.imr_validate_regression_draws <- function(object) {
  validate_imr_object(object$fit)
  groups <- object$fit$model$subgroup_names
  n <- object$control$output_draws
  if (!.imr_is_integerish(n) || length(n) != 1L || n < 2L ||
      !is.list(object$beta) || !identical(names(object$beta), groups) ||
      !is.list(object$variance) || !identical(names(object$variance), groups)) {
    .imr_abort("The regression posterior has inconsistent groups or output-draw controls; upgrade an older saved object explicitly.")
  }
  for (g in seq_along(groups)) {
    b <- object$beta[[g]]; v <- object$variance[[g]]
    if (!is.matrix(b) || !is.numeric(b) || nrow(b) != n ||
        any(!is.finite(b)) || !is.numeric(v) || length(v) != n ||
        any(!is.finite(v)) || any(v <= 0)) {
      .imr_abort("The regression posterior has invalid coefficient or variance draws.")
    }
    expected <- colnames(.imr_posterior_design(object$fit, g))
    if (!identical(colnames(b), expected)) {
      .imr_abort("Regression coefficient columns do not match the fitted design.")
    }
    if (!is.null(object$latent)) {
      z <- object$latent[[g]]
      if (!is.matrix(z) || !is.numeric(z) ||
          !identical(dim(z), as.integer(c(n, object$fit$model$sample_sizes[g]))) ||
          any(!is.finite(z))) .imr_abort("The regression posterior has invalid latent-response draws.")
    }
  }
  if (!.imr_is_integerish(object$selection_draw_index) || length(object$selection_draw_index) != n ||
      any(object$selection_draw_index > object$fit$control$mcmc$draws | object$selection_draw_index < 1)) {
    .imr_abort("Regression posterior source-draw indices are invalid.")
  }
  required <- c("subgroup", "selection_model", "output_draws", "draws_per_model_chain", "max_split_rhat")
  if (!is.data.frame(object$diagnostics) ||
      !all(required %in% names(object$diagnostics))) {
    .imr_abort("The regression posterior has incomplete conditional-chain diagnostics.")
  }
  invisible(TRUE)
}

# Old objects may be inspected after explicit structural conversion. Reusing
# their draws for a new inference must not silently change the target model.
.imr_require_current_updates <- function(object) {
  numerical <- .imr_fit_numerical_control(object$control)
  if (!identical(object$control$selection_update, "symmetric_mrf_hastings") ||
      !identical(numerical$prior_indexing, "coefficient_blocks")) {
    .imr_abort(paste0(
      "This saved fit uses historical updates or precision indexing. ",
      "Use its archived source for numerical replay, or refit with `imr()`. ",
      "Upgrading metadata does not correct posterior draws."
    ))
  }
  invisible(TRUE)
}

#' Upgrade Stored IMR Objects
#'
#' Convert a structurally complete 0.1.x `imr` object to the current named
#' schema, or update convention names in an earlier 0.2.0 fit. Corrupted or
#' incomplete objects must be refitted.
#'
#' @section What conversion preserves:
#' The conversion reorganizes an existing fit into the current named schema and
#' then validates it. It does not refit the model, regenerate random draws or
#' replace a historical update rule with the current model target. Stored
#' earlier convention labels become explicit update and precision records.
#' Regression-draw controls and diagnostic column names are also updated when
#' needed. The recorded original call is preserved as provenance. Historical
#' draws remain inspectable; new prediction, refitting and conditional sampling
#' require a fit produced under the current model conventions. Statistical
#' interpretation remains tied to the original fit's outcome scale and
#' computational settings. [imr()] documents those conventions; missing data or
#' insufficient response-scale information can require refitting before later
#' operations are available.
#'
#' @param object A saved `imr` fit or `imr_posterior` regression-draw object.
#' @return The validated object, with schema-version-3 fit metadata.
#' @examples
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            forced_prior_scale = 1, draws = 500, burnin = 250, seed = 1)
#'
#' # A current fit already uses the named schema, so it is returned unchanged
#' # after validation.
#' identical(upgrade_imr_object(fit), fit)
#' @export
upgrade_imr_object <- function(object) {
  if (inherits(object, "imr_posterior")) {
    object$fit <- upgrade_imr_object(object$fit)
    if (!is.null(object$model_draw) && is.null(object$selection_draw_index)) {
      object$selection_draw_index <- object$model_draw
      object$model_draw <- NULL
    }
    rename <- function(x, old, new) {
      if (old %in% names(x) && !new %in% names(x)) names(x)[names(x) == old] <- new
      x
    }
    object$control <- rename(object$control, "draws", "output_draws")
    object$control <- rename(object$control, "conditional_draws", "min_draws_per_model_chain")
    object$diagnostics <- rename(object$diagnostics, "model", "selection_model")
    object$diagnostics <- rename(object$diagnostics, "returned_draws", "output_draws")
    object$diagnostics <- rename(object$diagnostics, "conditional_draws", "draws_per_model_chain")
    validate_imr_object(object)
    return(object)
  }
  if (!inherits(object, "imr") || !is.list(object)) {
    .imr_abort("`object` must be a legacy `imr` fit.")
  }
  if (identical(object$schema_version, 3L)) {
    validate_imr_object(object)
    return(object)
  }
  if (identical(object$schema_version, 2L)) {
    if ("model_variant" %in% names(object$control)) {
      .imr_abort("Saved fit has conflicting model-variant fields; it cannot be upgraded.")
    }
    names(object$control)[names(object$control) == "method"] <- "model_variant"
    sampler <- object$control$sampler_method %||% "legacy"
    if (!is.character(sampler) || length(sampler) != 1L || is.na(sampler) ||
        !sampler %in% c("legacy", "paper", "original", "corrected")) {
      .imr_abort("Saved fit has an unknown sampler convention; it cannot be upgraded.")
    }
    if ("sampler_method" %in% names(object$control)) {
      names(object$control)[names(object$control) == "sampler_method"] <- "selection_update"
    }
    object$control$selection_update <- if (sampler %in% c("paper", "corrected"))
      "symmetric_mrf_hastings" else "unadjusted_flip_swap"
    if (is.null(object$control$standardize)) object$control$standardize <- TRUE
    if (is.null(object$control$numerical)) {
      object$control$numerical <- .imr_numerical_control()
    } else {
      indexing <- object$control$numerical$prior_indexing
      if (!is.character(indexing) || length(indexing) != 1L || is.na(indexing) ||
          !indexing %in% c("standard", "code2017", "original")) {
        .imr_abort("Saved fit has an unknown precision convention; it cannot be upgraded.")
      }
      object$control$numerical$prior_indexing <- if (indexing == "standard")
        "coefficient_blocks" else "historical_shifted_boundary"
    }
    object$schema_version <- 3L
    validate_imr_object(object)
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
    schema_version = 3L,
    control = list(
      call = object$call, outcome_type = object$type_outcome,
      response_scale = object$response_scale,
      model_variant = tolower(object$method), selection_update = "unadjusted_flip_swap", min_subgroup_size = object$ssize,
      priors = list(nu = object$nu, molecular_scale = h[[2L]],
        forced_scale = h[[1L]], residual = c(shape = h[[3L]], rate = h[[4L]]),
        interaction = c(shape = h[[5L]], rate = h[[6L]])),
      mcmc = list(draws = unname(object$sample_mcmc[["total"]]),
                  burnin = unname(object$sample_mcmc[["burnin"]])),
      seed = as.integer(h[[7L]]), numerical = .imr_numerical_control(),
      standardize = TRUE
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
  validate_imr_object(fit)
  fit
}


# One credible-interval row per parameter. The identifying columns differ
# between the selection indicators and the theta interactions; the location,
# scale and quantile summary does not.
.imr_draw_interval <- function(draws, probs, quantile_type = 8L, ...) {
  qs <- stats::quantile(draws, probs = probs, names = FALSE, type = quantile_type)
  data.frame(..., mean = mean(draws), sd = stats::sd(draws),
             lower = qs[1L], median = qs[2L], upper = qs[3L],
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Selection and Interaction Posterior Summaries
#'
#' Summarizes retained MCMC draws for variable-selection indicators and MRF
#' interaction parameters. The selection tables report posterior means,
#' posterior standard deviations and equal-tail credible intervals for each
#' platform-feature/subgroup indicator. The theta tables provide the same
#' summaries for each pair of linked availability subgroups.
#'
#' @section Summaries of retained parameters:
#' For retained values \eqn{a^{(1)},\ldots,a^{(B)}}{a[1], ..., a[B]} of a selection indicator
#' or MRF interaction, the reported mean, standard deviation and interval are
#' \deqn{\bar a=\frac{1}{B}\sum_b a^{(b)},}{mean(a) = sum_b a[b] / B,}
#' \deqn{s_a=\sqrt{\frac{1}{B-1}\sum_b(a^{(b)}-\bar a)^2},}{sd(a) = sqrt(sum_b (a[b] - mean(a))^2 / (B - 1)),}
#' \deqn{[Q_{(1-L)/2}(a),\ Q_{(1+L)/2}(a)],}{Equal-tail interval: [quantile(a, (1 - level)/2), quantile(a, (1 + level)/2)].}
#' where \eqn{L} is `level` and \eqn{Q} is the empirical quantile computed by
#' `stats::quantile(type = 8)`. The median is \eqn{Q_{0.5}}.
#' For a binary selection indicator, the mean is its mPIP. Type 8 interpolates
#' between ordered draws and can produce fractional interval endpoints even
#' though the indicator only takes values zero and one. These are interpolated
#' empirical quantiles, not exact credible sets on the indicator's support.
#' The interval is not an
#' interval for the Monte Carlo error of the estimated mPIP. A platform's theta
#' table is empty for BMS, or when that platform occurs in only one retained
#' subgroup and therefore has no subgroup pair to connect.
#'
#' This function summarizes the original fitted selection and interaction
#' draws without additional sampling. Its `sd` is posterior spread, not a
#' Monte Carlo standard error. For regression coefficient intervals, first
#' create an object with [sample_regression_posterior()] and use its `confint()` method.
#'
#' @param object A fitted `"imr"` object.
#' @param level Credible interval level between zero and one (default `0.95`).
#' @param ... Unused; present for future methods.
#' @return An object of class `"selection_summary.imr"` containing `selection`
#'   and `theta` tables.
#' @examples
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            forced_prior_scale = 1, draws = 500, burnin = 250, seed = 1)
#'
#' # Posterior mean, standard deviation and interval for every indicator,
#' # alongside the MRF interactions that link the availability subgroups.
#' summaries <- selection_summary(fit)
#' head(summaries$selection$assay)
#' @export
selection_summary <- function(object, ...) {
  UseMethod("selection_summary")
}

#' @rdname selection_summary
#' @export
selection_summary.imr <- function(object, level = 0.95, ...) {
  .imr_reject_dots(...)
  validate_imr_object(object)
  .imr_check_interval_level(level)
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
        rows[[at]] <- .imr_draw_interval(draws, probs,
          subgroup = rownames(template)[i],
          feature = colnames(template)[j])
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
      .imr_draw_interval(draws, probs,
        subgroup1 = subgroup_names[pairs[1L, j]],
        subgroup2 = subgroup_names[pairs[2L, j]])
    })
    do.call(rbind, rows)
  })
  names(theta) <- object$model$platform_names

  out <- list(level = level, selection = selection, theta = theta,
              selection_update = object$control$selection_update)
  class(out) <- "selection_summary.imr"
  out
}

#' @export
print.selection_summary.imr <- function(x, ...) {
  cat(sprintf("Selection and interaction summary (%.1f%% intervals)\n", 100 * x$level))
  if (!identical(x$selection_update, "symmetric_mrf_hastings"))
    cat("Stored historical draws; metadata conversion has not changed their target.\n")
  for (nm in names(x$selection)) {
    cat(sprintf(
      "  %s: %d selection indicators; %d theta interaction(s)\n",
      nm, nrow(x$selection[[nm]]), nrow(x$theta[[nm]])
    ))
  }
  invisible(x)
}


#' Credible Intervals for Selection and Interaction Parameters
#'
#' Standard `confint()` interface to [selection_summary()].
#'
#' @section Which parameters are summarized:
#' For an `imr` fit, this method extracts the selection-indicator and/or MRF
#' interaction tables computed by [selection_summary()]. That page defines
#' the posterior mean, standard deviation and equal-tail quantiles. Calling
#' `confint(fit)` does not sample regression coefficients. Their intervals
#' are obtained by `confint(sample_regression_posterior(fit))`; see
#' [imr_posterior_methods] for the different parameter set.
#'
#' @param object A fitted `"imr"` object.
#' @param parm Required parameter set: `"all"`, `"selection"` or `"theta"`.
#'   Naming the set explicitly distinguishes these intervals from regression
#'   coefficient intervals, available after [sample_regression_posterior()].
#' @param level Credible interval level.
#' @param ... Additional arguments passed to [selection_summary()].
#' @return A list of per-platform credible-interval tables, or a list with both
#'   selection and theta results when `parm = "all"`.
#' @export
confint.imr <- function(object, parm,
                        level = 0.95, ...) {
  if (missing(parm)) .imr_abort(paste0(
    "Specify `parm = \"selection\"`, `\"theta\"` or `\"all\"`. ",
    "For regression coefficient intervals, first use `sample_regression_posterior()`."
  ))
  parm <- match.arg(parm, c("all", "selection", "theta"))
  out <- selection_summary(object, level = level, ...)
  if (parm == "all") return(list(selection = out$selection, theta = out$theta))
  out[[parm]]
}


#' Compare Descriptive Summaries of IMR Fits
#'
#' Creates a compact descriptive comparison of compatible IMR fits. The table
#' deliberately does not treat raw log-posterior values as likelihood criteria;
#' it reports model structure and the number of features passing a common mPIP
#' threshold.
#'
#' @section Descriptive comparison:
#' The selected-feature count is
#' \deqn{N_{\mathrm{selected}}(t)=\sum_l |\mathcal A_l(t)|,}{Selected-feature count = sum over platforms of the number of features passing the common threshold.}
#' using the strict maximum-subgroup mPIP threshold set defined in
#' [summary.imr()]. Counts are by platform-feature pair. Other columns describe
#' the model variant, response scale, dimensions and retained sampling budget.
#'
#' The returned table is descriptive: it computes no Bayes factor, information
#' criterion or predictive ranking. Comparability checks require the same
#' outcome type and response scale, platforms, feature names and availability
#' subgroups. They do not check identical subject samples or CV folds. Fits on
#' log-time and identity scales cannot be combined in this table. For predictive model
#' comparison, evaluate prespecified candidates on matched subjects and folds
#' with [cv_imr()]. Choosing a specification from the data requires an outer
#' validation layer when estimating the performance of that choice.
#'
#' @param ... Fitted `"imr"` objects, or one list of fitted objects.
#' @param threshold Common mPIP threshold used to count selected features.
#' @return A data frame with one row per fit.
#' @examples
#' x <- data.frame(id = 1:40, marker = seq(-1, 1, length.out = 40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + sin(x$id) / 3)
#' fit <- imr(list(assay = x), y, outcome_type = "continuous",
#'            forced_prior_scale = 1, draws = 500, burnin = 250, seed = 1)
#' sparse <- imr(list(assay = x), y, outcome_type = "continuous",
#'               forced_prior_scale = 1, nu = -6,
#'               draws = 500, burnin = 250, seed = 1)
#'
#' # One row per fit: structure, and how many features clear the threshold.
#' compare_fit_summaries(default = fit, sparser_prior = sparse)
#' @export
compare_fit_summaries <- function(..., threshold = 0.5) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1L]]) &&
      !inherits(fits[[1L]], "imr")) {
    fits <- fits[[1L]]
  }
  if (length(fits) < 2L) {
    .imr_abort("Supply at least two fitted `imr` objects.")
  }
  threshold <- .imr_check_threshold(threshold)
  invisible(lapply(fits, validate_imr_object))
  reference <- fits[[1L]]
  compatible <- vapply(fits[-1L], function(fit) {
    identical(fit$control$outcome_type, reference$control$outcome_type) &&
      identical(fit$control$response_scale, reference$control$response_scale) &&
      identical(fit$model$platform_names, reference$model$platform_names) &&
      identical(fit$model$feature_names, reference$model$feature_names) &&
      identical(fit$model$subgroup_names, reference$model$subgroup_names)
  }, logical(1L))
  if (any(!compatible)) {
    .imr_abort(paste0(
      "Fits must have the same outcome type, response scale, platforms, features and ",
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
      model_variant = fit$control$model_variant,
      response_scale = fit$control$response_scale %||% NA_character_,
      platforms = fit$model$n_platforms,
      subgroups = length(fit$model$subgroup_names),
      retained_draws = fit$control$mcmc$draws,
      selected_features = selected,
      stringsAsFactors = FALSE, row.names = NULL
    )
  })
  do.call(rbind, rows)
}
