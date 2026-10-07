#' Fit Joint Integrative Bayesian Regressions
#'
#' Fits a regression in each availability subgroup and samples selection,
#' regression coefficients, residual variances and sharing parameters jointly.
#' Use [coef()], [confint()], [summary()] and [predict()] directly on the fit;
#' [posterior_draws()] extracts its stored samples without further MCMC.
#'
#' @section Model and computation:
#' On the working response scale, subgroup s has
#' \deqn{z_s=X_s\beta_s+\epsilon_s,\quad\epsilon_s\sim N(0,v_s I).}{z_s = X_s beta_s + error_s; error_s ~ Normal(0, v_s I).}
#' Continuous responses are observed. Binary responses constrain latent values
#' above/below zero. Survival uses log time by default; right-censored responses
#' are augmented above their observed censoring bounds.
#'
#' Every active coefficient, including the intercept and clinical effects, has
#' the first-order product-moment prior
#' \deqn{p(\beta_j\mid v_s)=\frac{\beta_j^2}{\tau_jv_s}\phi(\beta_j;0,\tau_jv_s).}{p(beta_j | v_s) = beta_j^2 / (tau_j v_s) * NormalDensity(beta_j; 0, tau_j v_s).}
#' Molecular effects are exactly zero when excluded. Clinical effects are always
#' active. `forced_scale` and `molecular_scale` specify their respective tau's.
#' Residual variance has an inverse-gamma shape/rate prior; binary fits anchor
#' its latent scale near one with shape and rate both equal to 100000.
#'
#' IMR shares selection information through the normalized symmetric MRF
#' \deqn{p(\gamma_{lj}\mid\Theta_l)\propto
#' \exp\{\nu_l\mathbf 1^T\gamma_{lj}+\gamma_{lj}^T\Theta_l\gamma_{lj}\}.}{p(gamma_lj | Theta_l) is proportional to exp(nu_l * sum(gamma_lj) + transpose(gamma_lj) * Theta_l * gamma_lj).}
#' BMS sets interactions to zero. The sampler jointly updates a feature's
#' selection indicators across its subgroups and then its active coefficients,
#' using exact scalar pMOM Bayes factors. Whole-feature exchanges aid movement
#' between correlated predictors. Variances and augmented responses use their
#' full conditionals; interactions use a Metropolis step with the exact MRF
#' normalizer. This engine does not use a Laplace model score.
#'
#' @section Scaling and interpretation:
#' When `standardize = TRUE`, a predictor in subgroup s becomes
#' \deqn{z_{sij}=(x_{sij}-\bar x_{sj})/s_{sj}.}{z_sij = (x_sij - training_mean_sj) / training_sd_sj.}
#' Constant columns retain their training mean and use scale one, so their
#' training values become zero; different new values need not become zero.
#' Prediction reuses the stored training transformations. Coefficients refer to
#' these transformed columns. With `standardize = FALSE`, supplied columns are
#' retained as given.
#'
#' Multiple chains and [mcmc_diagnostics()] assess finite-chain exploration.
#' Constant indicators have undefined R-hat and are reported separately, not
#' certified as converged. Neither a successful fit nor an R-hat cutoff alone
#' establishes adequate inference for a particular scientific analysis.
#'
#' @param x An [imr_data()] object, a formula, or a named list of platform data
#'   frames. Each frame contains unique `id` values and numeric feature columns.
#' @param outcome Outcome data frame with `id` and response; survival also has
#'   an event indicator (one for an observed event, zero for right censoring).
#' @param covariates Optional numeric clinical data frame, including `id`.
#' @param outcome_type `"continuous"`, `"binary"`, or `"right.censored"`.
#' @param model_variant `"imr"` shares information; `"bms"` fits independent
#'   subgroup regressions under the same coefficient and variance priors.
#' @param priors Statistical settings from [imr_priors()].
#' @param mcmc Sampling and storage settings from [imr_mcmc()].
#' @param min_subgroup_size Groups with at most this many subjects are excluded.
#' @param standardize Center/scale predictors within each training subgroup.
#' @param survival_scale `"log"` for an AFT model or `"identity"` for a Gaussian
#'   working-time model. Applies only to survival outcomes.
#' @param verbose Print progress for each chain.
#' @param ... Passed by the generic to its method; unsupported arguments fail.
#' @return An `imr` object with schema version 4, containing controls, model
#'   metadata, preprocessing, joint posterior arrays and MCMC diagnostics.
#'   Array dimensions preserve iteration, chain and parameter identities.
#' @seealso [posterior_draws()], [inclusion_probabilities()], [cv_imr()]
#' @references Chekouo T, Stingo FC, Doecke JD, Do K-A (2017).
#'   A Bayesian Integrative Approach for Multi-Platform Genomic Data: A Kidney
#'   Cancer Case Study. Biometrics, 73(2), 615--624. \doi{10.1111/biom.12587}.
#' @export
#' @examples
#' x <- data.frame(id = 1:40, marker = sin(1:40))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(draws = 500, burnin = 500, chains = 2, seed = 1)
#' )
#' coef(fit)
#' confint(fit)
imr <- function(x, ...) UseMethod("imr")

#' @rdname imr
#' @export
imr.default <- function(x, ...) {
  .imr_abort("`x` must be an imr_data object, formula, or platform list.")
}

#' @rdname imr
#' @export
imr.list <- function(x, outcome, covariates = NULL,
                     outcome_type = c("right.censored", "binary", "continuous"),
                     model_variant = c("imr", "bms"),
                     priors = imr_priors(), mcmc = imr_mcmc(),
                     min_subgroup_size = 30L, standardize = TRUE,
                     survival_scale = c("log", "identity"), verbose = FALSE, ...) {
  .imr_reject_dots(...)
  outcome_type <- match.arg(outcome_type)
  model_variant <- match.arg(model_variant)
  survival_scale <- match.arg(survival_scale)
  .imr_check_flag(standardize, "standardize")
  .imr_check_flag(verbose, "verbose")
  min_subgroup_size <- .imr_check_integer_scalar(min_subgroup_size, "min_subgroup_size", min = 0)
  priors <- .imr_validate_specification(priors, imr_priors, "imr_priors")
  mcmc <- .imr_validate_specification(mcmc, imr_mcmc, "imr_mcmc")
  data <- imr_data(x, outcome, covariates, outcome_type = outcome_type)
  if (length(priors$nu) == 1L) priors$nu <- rep(priors$nu, length(data$platforms))
  if (length(priors$nu) != length(data$platforms)) .imr_abort("`priors$nu` needs one value or one per platform.")
  default_residual <- if (outcome_type == "binary") c(shape = 1e5, rate = 1e5) else c(shape = .001, rate = .001)
  if (outcome_type == "binary" && !is.null(priors$residual) && !identical(priors$residual, default_residual)) {
    .imr_abort("Binary latent variance uses the anchored prior c(shape = 1e5, rate = 1e5); leave `residual = NULL` or specify that pair.")
  }
  priors$residual <- priors$residual %||% default_residual
  prepared <- .imr_prepare_fit(data, min_subgroup_size, standardize, survival_scale)
  control <- list(
    call = match.call(), outcome_type = outcome_type,
    model_variant = model_variant, response_scale = if (outcome_type == "right.censored") survival_scale else if (outcome_type == "binary") "probit" else "identity",
    priors = priors, mcmc = mcmc, rng_kind = RNGkind(),
    min_subgroup_size = min_subgroup_size, standardize = standardize,
    inference = "joint_pmom_mrf", package_version = "0.3.0"
  )
  spec <- .imr_joint_specification(prepared$model, prepared$preprocessing, control)
  parameters <- sum(vapply(spec$groups, function(g) ncol(g$design), 1L)) + length(spec$groups) + 1L
  if (model_variant == "imr") {
    parameters <- parameters + sum(vapply(
      spec$platforms,
      function(p) choose(length(p$groups), 2), 0
    ))
  }
  if (mcmc$keep_latent && outcome_type != "continuous") parameters <- parameters + sum(prepared$model$sample_sizes)
  expected_bytes <- 8 * as.double(mcmc$draws) * mcmc$chains * parameters
  if (!is.finite(expected_bytes) || expected_bytes > mcmc$max_draw_memory_mb * 1024^2) {
    .imr_abort(sprintf("Retained posterior arrays need about %.1f MiB; reduce retained `draws` or raise `mcmc$max_draw_memory_mb`. This excludes temporary fitting memory.", expected_bytes / 1024^2))
  }
  control$retained_draw_bytes <- expected_bytes
  seed <- mcmc$seed
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  saved_rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(saved_rng), add = TRUE)
  set.seed(seed)
  seeds <- matrix(sample.int(.Machine$integer.max, 2L * mcmc$chains, replace = FALSE), nrow = 2L)
  control$seed <- seed
  control$chain_seeds <- seeds[1L, ]
  control$initial_seeds <- seeds[2L, ]
  tasks <- lapply(seq_len(mcmc$chains), function(i) list(chain = i, seed = seeds[1L, i], initial_seed = seeds[2L, i]))
  chains <- .imr_map_tasks(tasks, .imr_joint_chain, mcmc$workers,
    spec = spec, model = prepared$model, control = control, verbose = verbose
  )
  control$initial <- lapply(chains, `[[`, "initial")
  control$acceptance <- data.frame(
    chain = seq_len(mcmc$chains),
    do.call(rbind, lapply(chains, `[[`, "acceptance"))
  )
  names(control$acceptance)[-1L] <- c("swap_proposals", "swap_accepts", "interaction_proposals", "interaction_accepts")
  fit <- structure(list(
    schema_version = 4L, control = control,
    model = prepared$model, preprocessing = prepared$preprocessing,
    posterior = .imr_combine_chains(chains, prepared$model, prepared$preprocessing, spec, mcmc),
    diagnostics = NULL
  ), class = "imr")
  validate_imr_object(fit)
  if (mcmc$diagnostics) {
    fit$diagnostics <- .imr_compute_diagnostics(fit)
    .imr_warn_diagnostics(fit$diagnostics)
  }
  fit
}

#' @rdname imr
#' @export
imr.imr_data <- function(x, ...) {
  validate_imr_data(x)
  if (is.null(x$outcome)) .imr_abort("Fitting requires an outcome in the imr_data object.")
  fit <- imr.list(x$platforms, x$outcome, x$covariates, outcome_type = x$outcome_type, ...)
  fit$control$call <- match.call()
  fit$preprocessing$input_data <- x
  fit
}
