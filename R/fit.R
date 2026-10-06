#' Fit the Integrative Bayesian Multi-Platform Regression Model
#'
#' @description
#' `imr()` fits the integrative multi-regression (IMR) model of Chekouo et al.
#' (2017) by Markov chain Monte Carlo (MCMC).
#' It identifies biomarkers associated with a time-to-event, binary or
#' continuous outcome while borrowing information across all subjects,
#' regardless of which platforms each subject has measured.  Subjects are
#' partitioned into availability subgroups defined by their pattern of platform
#' availability; one regression is built per subgroup and information is shared
#' across subgroups through a Markov random field (MRF) prior on the
#' variable-selection indicators together with non-local priors on the
#' regression coefficients.  The sampler is implemented in C for efficiency.
#'
#' The arguments are grouped by the orthogonal aspect of the analysis that each
#' one controls: the *data* (`x`, `outcome`, `covariates`), the
#' *likelihood* (`outcome_type`), the *model* (`method`), subgroup *filtering*
#' (`min_subgroup_size`), the *priors* (`nu`, `molecular_prior_scale`, `forced_prior_scale`, `residual_prior`, `interaction_prior`),
#' the *computation* (`draws`, `burnin`, `seed`) and the *output* (`verbose`).
#'
#' @section Outcome models and priors:
#' Let \eqn{s} index an availability subgroup. Its design matrix \eqn{Z_s}
#' contains an intercept, the forced clinical covariates and the selected
#' molecular features. On the working response scale,
#' \deqn{y_s^* = Z_s b_s + \epsilon_s, \qquad
#'       \epsilon_s \sim N(0, v_s I).}
#' For continuous outcomes, \eqn{y_s^*} is observed. For binary outcomes,
#' \eqn{Y_i = I(y_i^*>0)} and the latent normal response is sampled on the
#' appropriate side of zero. For the default survival model, an event has
#' \eqn{y_i^*=\log t_i}; a censored observation has
#' \eqn{y_i^*>\log t_i}. The binary prior concentrates \eqn{v_s} near one;
#' the implemented variance is not fixed exactly at one.
#'
#' Every active coefficient has the first-order pMOM density
#' \deqn{p(b_j\mid v_s) = \frac{b_j^2}{\tau_j v_s}
#'       \phi(b_j;0,\tau_j v_s),}
#' where \eqn{\phi(\cdot;\mu,v)} is a normal density with variance \eqn{v}.
#' The scale \eqn{\tau_j} is `forced_prior_scale` for the intercept and clinical
#' effects and `molecular_prior_scale` for molecular effects. Inactive molecular
#' coefficients are exactly zero. The residual prior is
#' \eqn{v_s\sim\mathrm{IG}(a,b)}, parameterized by a density proportional to
#' \eqn{v_s^{-a-1}\exp(-b/v_s)}.
#'
#' For feature \eqn{j} on platform \eqn{l}, collect its subgroup selection
#' indicators in \eqn{\gamma_{lj}}. The published MRF prior is
#' \deqn{p(\gamma_{lj}\mid\nu_l,\Theta_l) \propto
#'       \exp\{\nu_l\mathbf{1}^{T}\gamma_{lj}
#'                +\gamma_{lj}^{T}\Theta_l\gamma_{lj}\}.}
#' The symmetric matrix \eqn{\Theta_l} has zero diagonal and positive
#' off-diagonal interactions with the Gamma shape/rate `interaction_prior`.
#' For BMS, these interactions are zero. The prior conditional log-odds is
#' \eqn{\nu_l+2\sum_{h\ne s}\theta_{l,sh}\gamma_{l,hj}}; this factor of two
#' is used by `sampler_method = "paper"`. The legacy update differs as described
#' below. These equations describe the stated model, not a claim that the legacy
#' transition targets the same posterior.
#'
#' @section Computation and interpretation:
#' The selection sampler integrates regression coefficients and variances using
#' the package's Laplace-based marginal-likelihood calculation. For augmented
#' outcomes, this calculation uses the current latent responses. The retained
#' selection and interaction draws are summarized by [coef.imr()] and
#' [posterior_summary()]. Coefficient samples require the separate conditional
#' sampling step [posterior_draws()].
#'
#' With `standardize = TRUE`, each predictor uses its training subgroup's
#' \eqn{z_{ij}=(x_{ij}-\bar x_j)/s_j}. A constant training column is stored as
#' zero with scale one. [predict.imr()] reuses the stored centers and scales;
#' refit [cv_imr()] estimates them again within each training fold.
#'
#' @param x A list of data frames, one per platform.  Each
#'   data frame must contain an `id` column (the subject identifier, taken to be
#'   the first column); the remaining columns are finite numeric features
#'   measured on that platform.  Subject identifiers must be unique within each
#'   data frame.  The `id` column links subjects across platforms and to the
#'   outcome and covariate data.
#' @param outcome A data frame containing an `id` column and the response.  For
#'   `outcome_type = "right.censored"` it must also contain the event/censoring
#'   time and a censoring indicator (three columns in total).  For `"binary"`
#'   the response must be coded 0/1 and for `"continuous"` it is a numeric
#'   response (two columns in total).  The `id` column must be first and unique.
#' @param covariates An optional data frame of clinical covariates including an `id`
#'   column followed by finite numeric covariates.  Covariates are always
#'   included in every regression (they are not subject to selection).  Defaults
#'   to `NULL` (no covariates).
#' @param outcome_type Character string specifying the outcome type, one of
#'   `"right.censored"` (default), `"binary"` or `"continuous"`.
#' @param method Character string specifying the method, `"imr"` (default) for
#'   the integrative model that shares information across subgroups via the MRF
#'   prior, or `"bms"` for the non-integrative Bayesian multi-step variant that
#'   fits each subgroup independently (MRF interaction parameters set to zero).
#' @param min_subgroup_size Minimum availability subgroup size for a subgroup to be
#'   modelled. Subgroups with at most `min_subgroup_size` subjects are dropped. Default is
#'   `30`.
#' @param nu A numeric vector of prior log-odds of inclusion, one value per
#'   platform, controlling the prior sparsity of the selected features.
#'   Defaults to `rep(-3, length(x))`.
#' @param molecular_prior_scale Scale of the non-local (product moment) prior on the slab
#'   regression effects (default `0.087`).
#' @param forced_prior_scale pMOM prior scale for the always-included intercept and clinical
#'   coefficients (default `10000`), corresponding to tau_0 in the original
#'   paper. This is not the marginal prior variance.
#' @param residual_prior Length-2 numeric vector with the shape and rate of the
#'   inverse-gamma prior on the response error variance (default
#'   `c(0.001, 0.001)`).  Ignored when `outcome_type = "binary"`: a probit model
#'   uses a highly concentrated inverse-gamma prior with shape and rate
#'   `100000` to approximate unit residual variance for identifiability.
#' @param interaction_prior Length-2 numeric vector with the shape and rate of the
#'   gamma prior on the MRF interaction parameters theta, which borrow
#'   information across subgroups (default `c(40, 10)`).
#' @param draws Number of post-burn-in MCMC draws to retain (default `2000`).
#' @param burnin Number of initial MCMC iterations to discard (default `1000`).
#' @param seed Optional integer seed for the sampler.  If `NULL` (the default) a
#'   seed is drawn from the current R RNG state, so a run is reproducible
#'   whenever [set.seed()] is called beforehand or an explicit `seed` is passed.
#' @param verbose Logical; if `TRUE`, print the sampler's progress and
#'   diagnostics to the console.  Defaults to `FALSE` (quiet).
#' @param survival_scale Working response scale for right-censored outcomes.
#'   `"log"` (default) logs the supplied positive event/censoring times, as in
#'   the original AFT model. `"identity"` reproduces historical package analyses
#'   and does not fit a log-time AFT model. Ignored for other outcome types.
#' @param sampler_method Selection-update convention. `"legacy"` (default)
#'   retains the existing package and 2017 code updates. `"paper"` uses the
#'   symmetric MRF conditional log-odds, the boundary flip/swap Hastings
#'   correction, and the negative Gamma rate term in the log-posterior trace.
#'   This is a computational convention, separate from IMR/BMS `method`.
#' @param prior_indexing `"standard"` (default) assigns prior precision by
#'   coefficient block. `"code2017"` reproduces the strict index boundaries
#'   in the released C code; in particular, the last forced covariate receives
#'   molecular precision when covariates are present. Use for historical
#'   comparisons, not as a different scientific prior specification.
#' @param laplace_max_iter Maximum coefficient-mode iterations, either a
#'   positive integer for all stages or a named vector in the order `initial`,
#'   `selection`, `latent`, `prediction`. Defaults are 25, 40, 25, 40. The
#'   released C code used 25 for prediction; these are not MCMC draw counts.
#' @param laplace_tolerance Positive relative coefficient-mode convergence
#'   tolerance (default `0.001`). Numerical controls are stored in the fit and
#'   reused by prediction and cross-validation, including training-fold refits.
#'   `control$laplace_diagnostics` records calls, iteration-limit hits, nonfinite
#'   likelihoods and factorization failures by subgroup for the initial,
#'   selection and latent fitting stages. These counters do not assess MCMC
#'   convergence or subsequent prediction-stage optimization.
#' @param initial Optional named list with `selection` and/or `interaction`.
#'   Each component is a list of matrices named by platform. Selection matrices
#'   have subgroup rows and feature columns, with entries zero or one; interaction
#'   matrices have subgroup rows/columns, positive symmetric off-diagonals and
#'   zero diagonals (IMR only). Dimnames must match the fitted model order.
#'   NULL preserves the original random initialization and RNG consumption.
#'   Supplied selection bypasses its initialization draws. Refit CV reuses the
#'   specified starting values; missing components retain their defaults.
#' @param standardize Logical; standardize features and forced covariates
#'   within availability subgroups (default `TRUE`). Set `FALSE` for inputs
#'   already transformed by an external historical experiment. Prediction then
#'   uses those supplied units. Refit CV cannot undo external preprocessing;
#'   record how the supplied data were constructed.
#' @param ... Additional fitting arguments passed from the formula or
#'   `imr_data` method to the list method. Unused arguments are rejected.
#'
#' @details
#' All feature and covariate data are standardized internally (mean 0, standard
#' deviation 1); the centring and scaling factors are stored in the returned
#' object so that [predict.imr()] can apply the same transformation to new
#' subjects.  For right-censored outcomes the latent log-survival times of
#' censored subjects are imputed within the sampler; for binary outcomes a
#' probit data-augmentation latent variable is sampled.
#'
#' For the symmetric interaction matrix in the paper, the conditional
#' log-odds contribution is `nu + 2 * sum(theta * neighboring_indicators)`.
#' The released code used a factor of one and omitted the proposal ratio when
#' a flip moved between an empty/full model and an interior model. These
#' conventions remain available as `sampler_method = "legacy"`; they are not
#' mathematically equivalent to the stated MRF posterior. The `"paper"` option
#' corrects these updates without changing defaults. Exact historical table
#' reproduction additionally depends on data, preprocessing, initialization,
#' random-number consumption and validation settings.
#'
#' @return An object of class `"imr"` with `schema_version = 2L` and four
#'   named sections: `control` (outcome, method, priors, MCMC and seed), `model`
#'   (platform, feature and subgroup metadata), `preprocessing` (validated
#'   inputs, standardized matrices and formula metadata), and `posterior`
#'   (inclusion probabilities, interaction draws, latent-response summaries
#'   and the log-posterior trace).
#'
#' @references
#' Chekouo T, Stingo FC, Doecke JD, Do K-A (2017). "A Bayesian Integrative
#' Approach for Multi-Platform Genomic Data: A Kidney Cancer Case Study."
#' \emph{Biometrics}, \strong{73}(2), 615--624. \doi{10.1111/biom.12587}
#'
#' @seealso [predict.imr()], [cv_imr()], [summary.imr()], [plot.imr()]
#'
#' @examples
#' \donttest{
#' data("simIMR", package = "IntegMultiReg")
#' fit <- imr(
#'   x = simIMR$platforms,
#'   outcome = simIMR$outcome,
#'   covariates = simIMR$covariates,
#'   outcome_type = "binary",
#'   nu = c(-4, -3, -4),
#'   draws = 200, burnin = 100,
#'   min_subgroup_size = 5,
#'   seed = 1
#' )
#' fit
#' }
#' @export
imr <- function(x, ...) UseMethod("imr")

#' @rdname imr
#' @export
imr.default <- function(x, ...) {
  .imr_abort("`x` must be a platform list, formula, or `imr_data` object.")
}

.imr_new_fit <- function(control, model, preprocessing, posterior) {
  structure(
    list(
      schema_version = 2L,
      control = control,
      model = model,
      preprocessing = preprocessing,
      posterior = posterior
    ),
    class = "imr"
  )
}

#' @rdname imr
#' @export
imr.list <- function(x, outcome, covariates = NULL,
                     outcome_type = c("right.censored", "binary", "continuous"),
                     method = c("imr", "bms"), min_subgroup_size = 30L,
                     nu = rep(-3, length(x)), molecular_prior_scale = 0.087,
                     forced_prior_scale = 10000,
                     residual_prior = c(shape = 0.001, rate = 0.001),
                     interaction_prior = c(shape = 40, rate = 10),
                     draws = 2000L, burnin = 1000L, seed = NULL,
                     verbose = FALSE,
                     survival_scale = c("log", "identity"),
                     sampler_method = c("legacy", "paper"),
                     prior_indexing = c("standard", "code2017"),
                     laplace_max_iter = c(initial = 25L, selection = 40L,
                                          latent = 25L, prediction = 40L),
                     laplace_tolerance = 1e-3, standardize = TRUE, initial = NULL, ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    .imr_abort(sprintf("Unused argument: `%s`.", names(dots)[1L]))
  }
  call <- match.call()
  outcome_type <- match.arg(outcome_type)
  survival_scale <- match.arg(survival_scale)
  sampler_method <- match.arg(sampler_method)
  numerical <- .imr_numerical_control(match.arg(prior_indexing),
                                      laplace_max_iter, laplace_tolerance)
  method <- match.arg(method)

  .imr_check_flag(verbose, "verbose")
  .imr_check_flag(standardize, "standardize")
  validated <- imr_data(
    platforms = x,
    outcome = outcome,
    covariates = covariates,
    outcome_type = outcome_type
  )
  platforms <- validated$platforms
  outcome <- validated$outcome
  covariates <- validated$covariates
  n_platforms <- length(platforms)

  min_subgroup_size <- .imr_check_integer_scalar(
    min_subgroup_size, "min_subgroup_size", min = 0
  )
  nu <- .imr_check_numeric_vector(nu, "nu", length = n_platforms)
  forced_prior_scale <- .imr_check_numeric_vector(
    forced_prior_scale, "forced_prior_scale", length = 1, positive = TRUE
  )
  molecular_prior_scale <- .imr_check_numeric_vector(
    molecular_prior_scale, "molecular_prior_scale", length = 1, positive = TRUE
  )
  residual_prior <- .imr_check_named_pair(residual_prior, "residual_prior")
  interaction_prior <- .imr_check_named_pair(interaction_prior, "interaction_prior")
  draws <- .imr_check_integer_scalar(draws, "draws", min = 1)
  burnin <- .imr_check_integer_scalar(burnin, "burnin", min = 0)
  if (as.double(draws) + as.double(burnin) > .Machine$integer.max) {
    .imr_abort("The sum of `draws` and `burnin` must not exceed the native integer limit.")
  }
  if (!is.null(seed)) {
    seed <- .imr_check_integer_scalar(seed, "seed", min = 0)
  }

  ## The sampler draws on both R's RNG and a GSL RNG seeded by `seed`.  With an
  ## explicit seed we call set.seed() so the run is fully reproducible; when
  ## `seed` is NULL we instead draw one from the ambient R RNG, so results follow
  ## the user's own set.seed() like other R modelling functions.
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  } else {
    set.seed(seed)
  }

  ## Record human-readable platform and feature names (the first column of each
  ## platform is the 'id' and is dropped before modelling).
  platform_names <- names(platforms)
  if (is.null(platform_names) || any(platform_names == "")) {
    platform_names <- paste0("platform", seq_len(n_platforms))
  }
  feature_names <- lapply(platforms, function(platform) colnames(platform)[-1])

  dat <- .imr_subgroup_data(outcome, covariates, platforms)
  # Prepare scalar parameters.
  h0_c <- as.numeric(forced_prior_scale)
  hh_c <- as.numeric(molecular_prior_scale)
  alpha_c <- as.numeric(residual_prior[["shape"]])
  psi_c <- as.numeric(residual_prior[["rate"]])
  if (outcome_type == "binary") {
    # A probit outcome fixes the residual variance at 1 for identifiability
    # (the latent utility is z = eta + e with e ~ N(0, 1)).  The marginal
    # likelihood otherwise integrates sigma^2 out under this inverse-gamma
    # prior, which leaves the latent scale only weakly identified and makes the
    # binary chain drift and mix poorly.  Concentrating the prior at 1 pins
    # sigma^2 = 1; any large shape/rate gives an effectively fixed unit variance.
    probit_unit_variance <- 1e5
    alpha_c <- psi_c <- probit_unit_variance
  }
  alpha0_c <- as.numeric(interaction_prior[["shape"]])
  beta0_c <- as.numeric(interaction_prior[["rate"]])
  seed_c <- as.numeric(seed)

  storage.mode(h0_c) <- "double"
  storage.mode(hh_c) <- "double"
  storage.mode(alpha_c) <- "double"
  storage.mode(psi_c) <- "double"
  storage.mode(alpha0_c) <- "double"
  storage.mode(beta0_c) <- "double"
  storage.mode(seed_c) <- "double"

  #################################################################
  # Process the platform data and extract subgroup id vectors.
  #################################################################
  dat_orig <- dat # Preserve the original data (with "id" columns intact)
  # For each subgroup, extract the id vector from the first nonempty platform.
  dat_ids <- lapply(dat_orig[[3]], function(subgroup) {
    nonempty <- which(sapply(subgroup, function(platform_data) {
      if (is.data.frame(platform_data) && "id" %in% colnames(platform_data)) {
        return(nrow(platform_data))
      } else {
        return(0)
      }
    }) > 0)
    if (length(nonempty) > 0) {
      return(subgroup[[nonempty[1]]]$id)
    } else {
      return(character(0))
    }
  })

  # Now convert the platform data into a list of numeric matrices by dropping
  # the "id" column.
  dat[[3]] <- .imr_platform_matrices(dat[[3]])

  n_features <- as.integer(vapply(
    platforms, function(platform) ncol(platform) - 1L, integer(1)
  ))
  # Compute representative row counts for each subgroup.
  subgroup_rows <- sapply(dat[[3]], function(subgroup) {
    nonempty <- which(sapply(subgroup, function(mat) nrow(mat)) > 0)
    if (length(nonempty) > 0) {
      nrow(subgroup[[nonempty[1]]])
    } else {
      0
    }
  })

  # Keep only subgroups above the requested minimum size.
  model_index <- as.integer(which(subgroup_rows > min_subgroup_size))
  if (length(model_index) == 0) {
    .imr_abort(
      "No availability subgroup has more than `min_subgroup_size` subjects; lower it or check the data."
    )
  }
  sample_size <- as.integer(subgroup_rows[model_index])

  # Filter the platform data to include only subgroups meeting the threshold.
  dat_filtered <- lapply(dat, function(x) x[model_index])

  # Also filter the original id vectors accordingly.
  subgroup_ids_filtered <- dat_ids[model_index]

  #############################################################
  # Process the response and covariates using the filtered subgroup ids.
  #############################################################
  ### Match ids and drop the id column from the outcome and covariate frames.
  dat_filtered[[1]] <- lapply(
    seq_along(dat_filtered[[1]]),
    function(i) {
      .imr_match_rows(
        dat_filtered[[1]][[i]], subgroup_ids_filtered[[i]], "outcome"
      )[, -1, drop = FALSE]
    }
  )
  if (!is.null(covariates)) {
    dat_filtered[[2]] <- lapply(
      seq_along(dat_filtered[[2]]),
      function(i) {
        .imr_match_rows(
          dat_filtered[[2]][[i]], subgroup_ids_filtered[[i]], "covariates"
        )[, -1, drop = FALSE]
      }
    )
  }

  dat_normalized <- dat_filtered

  platform_preprocessing <- lapply(dat_filtered[[3]], function(platforms) {
    lapply(platforms, .imr_prepare_matrix, standardize = standardize)
  })
  platform_field <- function(field) lapply(platform_preprocessing, function(platforms) {
    lapply(platforms, `[[`, field)
  })
  mean_train <- platform_field("mean")
  sd_train <- platform_field("sd")
  dat_normalized[[3]] <- platform_field("normalized")
  if (!is.null(covariates)) {
    covariate_preprocessing <- lapply(dat_filtered[[2]], .imr_prepare_matrix,
                                     standardize = standardize)
    mean_cov_train <- lapply(covariate_preprocessing, `[[`, "mean")
    sd_cov_train <- lapply(covariate_preprocessing, `[[`, "sd")
    dat_normalized[[2]] <- lapply(covariate_preprocessing, `[[`, "normalized")
  } else {
    mean_cov_train <- NULL
    sd_cov_train <- NULL
  }

  ## outcome: force double storage (the C sampler reads it with REAL())
  dat_normalized[[1]] <- lapply(dat_filtered[[1]], function(x) {
    m <- as.matrix(x)
    storage.mode(m) <- "double"
    m
  })

  if (verbose) {
    cat("Sample sizes for modelled availability subgroups:", sample_size, "\n")
  }

  ### For each model, list the corresponding platforms involved in the model
  n_models <- length(model_index)

  model_platforms_c <- sapply(1:n_models, function(x) {
    as.integer((seq_along(dat_normalized[[3]][[x]]) - 1)[unlist(lapply(
      seq_along(dat_normalized[[3]][[x]]),
      function(i) nrow(dat_normalized[[3]][[x]][[i]])
    )) > 0])
  }, simplify = FALSE)

  ### For each platform, obtain the model indices where that platform is involved
  platform_models_c <- lapply(seq_len(n_platforms), function(x) {
    as.integer((seq(1, n_models) - 1)[unlist(lapply(
      model_platforms_c,
      function(y) (x - 1) %in% y
    ))])
  })
  .imr_check_mrf_capacity(platform_models_c)

  n_platform_c <- n_platforms
  x_filtered <- dat_normalized[[3]]
  y_list <- dat_normalized[[1]]
  if (outcome_type == "right.censored" && survival_scale == "log") {
    y_list <- lapply(y_list, function(y) {
      y[, 1L] <- log(y[, 1L])
      y
    })
  }

  if (!is.null(covariates)) {
    n_cov <- ncol(dat_normalized[[2]][[1]])
  } else {
    n_cov <- 0
  }
  if (n_cov == 0) {
    cov_list <- lapply(y_list, function(x) {
      matrix(numeric(0), nrow = nrow(x), ncol = 0)
    })
  } else {
    cov_list <- dat_normalized[[2]]
  }

  nu_c <- as.double(nu)

  ###########################################
  # Call the compiled MCMC sampler.
  ###########################################
  effective_priors <- list(
    forced_scale = h0_c, molecular_scale = hh_c,
    residual = c(shape = alpha_c, rate = psi_c),
    interaction = c(shape = alpha0_c, rate = beta0_c)
  )
  initial <- .imr_initial_state(initial, platform_names,
    lapply(platform_models_c, function(i) names(x_filtered)[i + 1L]), feature_names, method)
  results <- .imr_call_fit_native(
    priors = effective_priors, seed = seed_c, nu = nu_c, method = method,
    n_platforms = n_platform_c, platform_subgroups = platform_models_c,
    subgroup_platforms = model_platforms_c, sample_sizes = sample_size,
    n_features = n_features, n_covariates = n_cov, features = x_filtered,
    response = y_list, outcome_type = outcome_type, covariates = cov_list,
    draws = draws, burnin = burnin, verbose = verbose, sampler_method = sampler_method, numerical = numerical, initial = initial
  )

  ## Guard against tiny floating-point drift in the running averages so that
  ## the reported inclusion probabilities are exactly within [0, 1].
  results$gam_mean <- lapply(results$gam_mean, function(m) {
    m[m > 1] <- 1
    m[m < 0] <- 0
    m
  })

  subgroup_names <- names(x_filtered)
  fit <- .imr_new_fit(
    control = list(
      call = call, outcome_type = outcome_type,
      response_scale = if (outcome_type == "right.censored") survival_scale else
        if (outcome_type == "binary") "probit" else "identity",
      method = method, sampler_method = sampler_method, min_subgroup_size = min_subgroup_size,
      priors = list(nu = nu, molecular_scale = molecular_prior_scale,
        forced_scale = forced_prior_scale,
        residual = c(shape = alpha_c, rate = psi_c),
        interaction = interaction_prior),
      mcmc = list(draws = draws, burnin = burnin), seed = seed, numerical = numerical,
      standardize = standardize, initial = initial,
      laplace_diagnostics = data.frame(
        subgroup = rep(subgroup_names, each = 3L),
        stage = rep(c("initial", "selection", "latent"), length(subgroup_names)),
        stats::setNames(as.data.frame(results$laplace_diagnostics),
          c("calls", "iteration_limit", "nonfinite", "factorization_failures"))),
      rng_state = list(generator = "gsl_rand48", endian = .Platform$endian,
                       state = results$rng_state)
    ),
    model = list(
      n_platforms = n_platforms, platform_names = platform_names,
      feature_names = feature_names,
      covariate_names = if (!is.null(covariates)) colnames(covariates)[-1] else character(),
      subgroup_names = subgroup_names,
      sample_sizes = stats::setNames(sample_size, subgroup_names),
      subgroup_platforms = lapply(model_platforms_c, function(index) index + 1L),
      platform_subgroups = lapply(platform_models_c, function(index) index + 1L)
    ),
    preprocessing = list(
      input_data = validated, features = x_filtered, response = y_list,
      covariates = cov_list, feature_center = mean_train,
      feature_scale = sd_train, covariate_center = mean_cov_train,
      covariate_scale = sd_cov_train, formula = NULL, formula_data = NULL,
      terms = NULL, contrasts = NULL, xlevels = NULL, id = "id"
    ),
    posterior = list(
      inclusion_probabilities = results$gam_mean,
      interaction_means = results$theta_mean,
      latent_response_mean = results$estimate_latent_y,
      log_posterior = results$log_posterior,
      selection_draws = results$gam_sample,
      interaction_draws = results$theta_sample
    )
  )
  validate_imr(fit)
  fit
}


#' @rdname imr
#' @export
imr.imr_data <- function(x, ...) {
  validate_imr_data(x)
  if (is.null(x$outcome)) {
    .imr_abort("An `imr_data` object used for fitting must contain an outcome.")
  }
  fit <- imr.list(
    x = x$platforms,
    outcome = x$outcome,
    covariates = x$covariates,
    outcome_type = x$outcome_type,
    ...
  )
  fit$control$call <- match.call()
  fit$preprocessing$input_data <- x
  fit
}
