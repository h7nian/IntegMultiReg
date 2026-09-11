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
#' @param ... Additional fitting arguments passed from the formula or
#'   `imr_data` method to the default method. Unused arguments are rejected.
#'
#' @details
#' All feature and covariate data are standardized internally (mean 0, standard
#' deviation 1); the centring and scaling factors are stored in the returned
#' object so that [predict.imr()] can apply the same transformation to new
#' subjects.  For right-censored outcomes the latent log-survival times of
#' censored subjects are imputed within the sampler; for binary outcomes a
#' probit data-augmentation latent variable is sampled.
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
                     survival_scale = c("log", "identity"), ...) {
  dots <- list(...)
  if (length(dots) > 0L) {
    .imr_abort(sprintf("Unused argument: `%s`.", names(dots)[1L]))
  }
  call <- match.call()
  outcome_type <- match.arg(outcome_type)
  survival_scale <- match.arg(survival_scale)
  method <- match.arg(method)

  .imr_check_flag(verbose, "verbose")
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
  outcome_code <- match(outcome_type, c("right.censored", "binary", "continuous"))

  ## Record human-readable platform and feature names (the first column of each
  ## platform is the 'id' and is dropped before modelling).
  platform_names <- names(platforms)
  if (is.null(platform_names) || any(platform_names == "")) {
    platform_names <- paste0("platform", seq_len(n_platforms))
  }
  feature_names <- lapply(platforms, function(platform) colnames(platform)[-1])

  dat <- subgroup_data(outcome, covariates, platforms)
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

  method_c <- toupper(method)

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
  dat[[3]] <- lapply(dat[[3]], function(subgroup) {
    lapply(subgroup, function(platform_data) {
      if (is.data.frame(platform_data) && "id" %in% colnames(platform_data)) {
        mat <- as.matrix(platform_data[, -1, drop = FALSE])
      } else if (!is.null(colnames(platform_data)) &&
        colnames(platform_data)[1] == "id") {
        mat <- as.matrix(platform_data[, -1, drop = FALSE])
      } else {
        mat <- as.matrix(platform_data)
      }
      storage.mode(mat) <- "double"
      return(mat)
    })
  })

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

  mean_train <- mean_nested_list(dat_filtered[[3]])
  sd_train <- sd_nested_list(dat_filtered[[3]])
  if (!is.null(covariates)) {
    mean_cov_train <- lapply(dat_filtered[[2]], mean_matrix)
    sd_cov_train <- lapply(dat_filtered[[2]], sd_matrix)
  } else {
    mean_cov_train <- NULL
    sd_cov_train <- NULL
  }

  dat_normalized[[3]] <- normalize_nested_list(dat_filtered[[3]])
  ### Normalized covariates
  if (!is.null(covariates)) {
    dat_normalized[[2]] <- lapply(dat_filtered[[2]], normalize_matrix)
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
  results <- .quietly(verbose, .Call("imr_fit", h0_c, hh_c, alpha_c, psi_c,
    alpha0_c, beta0_c, seed_c, nu_c, method_c,
    n_platform_c = as.integer(n_platform_c),
    platform_models_c = platform_models_c, model_platforms_c = model_platforms_c,
    n_models = as.integer(n_models),
    sample_size = as.integer(sample_size),
    n_features = as.integer(n_features),
    n_cov = as.integer(n_cov),
    x_filtered = x_filtered, y_list = y_list,
    outcome_type = as.integer(outcome_code),
    cov_list = cov_list,
    sample = as.integer(draws),
    burnin = as.integer(burnin)
  ))

  ## Guard against tiny floating-point drift in the running averages so that
  ## the reported inclusion probabilities are exactly within [0, 1].
  results$gam_mean <- lapply(results$gam_mean, function(m) {
    m[m > 1] <- 1
    m[m < 0] <- 0
    m
  })

  subgroup_names <- names(x_filtered)
  fit <- list(
    schema_version = 2L,
    control = list(
      call = call, outcome_type = outcome_type,
      response_scale = if (outcome_type == "right.censored") survival_scale else
        if (outcome_type == "binary") "probit" else "identity",
      method = method, min_subgroup_size = min_subgroup_size,
      priors = list(nu = nu, molecular_scale = molecular_prior_scale,
        forced_scale = forced_prior_scale,
        residual = c(shape = alpha_c, rate = psi_c),
        interaction = interaction_prior),
      mcmc = list(draws = draws, burnin = burnin), seed = seed
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
  class(fit) <- "imr"
  fit
}


#' @rdname imr
#' @param formula A model formula passed as `x`.
#'   The model always includes an intercept: `0`/`-1` and `offset()` terms are
#'   rejected. The identifier column is excluded when expanding `.`.
#' @param data A data frame used with the formula interface.
#' @param platforms A list of platform data frames used with the formula
#'   interface.
#' @param id Name of the identifier column in `data` and `platforms`.
#' @export
imr.formula <- function(x, data, platforms, id = "id",
                        outcome_type = c("right.censored", "binary", "continuous"),
                        ...) {
  formula <- x
  outcome_type <- match.arg(outcome_type)
  if (!inherits(formula, "formula")) {
    .imr_abort("The first argument must be a formula.")
  }
  if (!is.data.frame(data)) {
    .imr_abort("`data` must be a data frame for the formula interface.")
  }
  if (!id %in% names(data)) {
    .imr_abort(sprintf("`data` must contain the identifier column `%s`.", id))
  }
  if (anyNA(data[[id]]) || anyDuplicated(data[[id]])) {
    .imr_abort("The identifier column in `data` must be complete and unique.")
  }

  # The identifier aligns subjects; it must not become a predictor via `.`.
  formula_data <- data[, setdiff(names(data), id), drop = FALSE]
  mf <- stats::model.frame(formula, data = formula_data,
                           na.action = stats::na.fail)
  response <- stats::model.response(mf)
  terms_object <- stats::terms(mf)
  if (attr(terms_object, "intercept") != 1L) {
    .imr_abort("IMR requires an intercept; formulas with `0` or `-1` are not supported.")
  }
  if (length(attr(terms_object, "offset")) > 0L) {
    .imr_abort("IMR does not support `offset()` terms in formulas.")
  }
  model_matrix <- stats::model.matrix(terms_object, mf)
  contrasts <- attr(model_matrix, "contrasts")
  xlevels <- stats::.getXlevels(terms_object, mf)
  keep <- attr(model_matrix, "assign") != 0L
  model_matrix <- model_matrix[, keep, drop = FALSE]

  response_matrix <- if (is.matrix(response) || is.data.frame(response)) {
    as.matrix(response)
  } else {
    matrix(response, ncol = 1L)
  }
  expected_response_columns <- if (outcome_type == "right.censored") 2L else 1L
  if (ncol(response_matrix) != expected_response_columns) {
    .imr_abort(sprintf(
      "The formula response must produce %d column(s) for `%s` outcomes.",
      expected_response_columns, outcome_type
    ))
  }
  outcome <- data.frame(
    id = data[[id]], response_matrix,
    check.names = FALSE, row.names = NULL
  )
  names(outcome)[1L] <- id
  names(outcome) <- c(
    id, if (outcome_type == "right.censored") c("time", "status") else "response"
  )
  covariates <- if (ncol(model_matrix) == 0L) {
    NULL
  } else {
    data.frame(
      id = data[[id]], model_matrix,
      check.names = FALSE, row.names = NULL
    )
  }
  if (!is.null(covariates)) names(covariates)[1L] <- id
  dat <- imr_data(
    platforms = platforms, outcome = outcome, covariates = covariates,
    outcome_type = outcome_type, id = id
  )
  fit <- imr(dat, ...)
  fit$control$call <- match.call()
  required_columns <- unique(c(id, all.vars(terms_object)))
  fit$preprocessing$formula_data <- data[, required_columns, drop = FALSE]
  fit$preprocessing$formula <- formula
  fit$preprocessing$terms <- terms_object
  # Single-bracket assignment preserves an explicit NULL in the fixed schema.
  fit$preprocessing["contrasts"] <- list(contrasts)
  fit$preprocessing$xlevels <- xlevels
  fit$preprocessing$id <- id
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


#' @keywords internal
#' @noRd
subgroup_data <- function(outcome, covariates = NULL, platforms) {
  # Collect all the input data frames into a list
  nplat <- length(platforms)
  ### Intersection of outcome and covariate ids with the union of platform ids
  ids_out_cov <- outcome$id
  if (!is.null(covariates)) {
    ids_out_cov <- intersect(outcome$id, covariates$id)
  }

  platforms <- lapply(platforms, function(x) {
    x[x$id %in% ids_out_cov, , drop = FALSE]
  })

  id_outcome <- unique(unlist(lapply(platforms, function(x) x$id)))
  if (length(id_outcome) == 0L) {
    .imr_abort(
      "No subjects have both outcome/covariate data and at least one platform."
    )
  }
  outcome1 <- outcome[outcome$id %in% id_outcome, , drop = FALSE]
  cov1 <- if (!is.null(covariates)) {
    covariates[covariates$id %in% id_outcome, , drop = FALSE]
  } else {
    NULL
  }

  # Check that each platform has at least one column named "id"
  if (!all(sapply(platforms, function(df) "id" %in% colnames(df)))) {
    .imr_abort("Each input data frame must have an `id` column.")
  }

  # Create the union of all ids across platforms
  all_ids <- unique(unlist(lapply(platforms, function(df) df$id)))

  # Build a presence matrix (rows = subjects, columns = platforms).
  presence <- do.call(cbind, lapply(platforms, function(df) all_ids %in% df$id))

  # Create a binary string for each subject.
  # Convention: the rightmost digit is presence in the first platform, etc.
  bitstrings <- apply(presence, 1, function(x) paste(as.integer(rev(x)), collapse = ""))

  # Identify the unique binary patterns (subgroups) sorted in ascending order.
  unique_patterns <- sort(unique(bitstrings))

  x1 <- list()
  sample_ids <- list()

  for (pat in unique_patterns) {
    subgroup_ids <- all_ids[bitstrings == pat]
    subgroup_list <- vector("list", nplat)
    names(subgroup_list) <- paste0("platform", seq_len(nplat))
    for (i in seq_len(nplat)) {
      # For the i-th platform, the corresponding bit is at position nplat - i + 1.
      bit <- substr(pat, nplat - i + 1, nplat - i + 1)
      if (bit == "1") {
        subgroup_list[[i]] <- .imr_match_rows(
          platforms[[i]], subgroup_ids, sprintf("x[[%d]]", i)
        )
      } else {
        subgroup_list[[i]] <- platforms[[i]][FALSE, , drop = FALSE]
      }
    }
    sample_ids[[pat]] <- subgroup_ids
    x1[[pat]] <- subgroup_list
  }
  ## Align the outcome and covariate rows to each subgroup's canonical subject
  ## order (`sample_ids`, taken from the platform data) with match(), rather than
  ## relying on the inputs sharing the same row order.  The compiled sampler
  ## pairs outcome / covariate / platform rows by position, so this keeps them
  ## correctly aligned even when the outcome or covariate frames are supplied in
  ## a different order (for id-sorted inputs it is a no-op).
  outcome2 <- lapply(sample_ids, function(x) .imr_match_rows(outcome1, x, "outcome"))
  cov2 <- if (!is.null(cov1)) {
    lapply(sample_ids, function(x) .imr_match_rows(cov1, x, "covariates"))
  } else {
    lapply(sample_ids, function(x) data.frame(id = x))
  }

  return(list(outcome = outcome2, covariate = cov2, platform_data = x1))
}

#' @keywords internal
#' @noRd
normalize_matrix <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(matrix(numeric(0), nrow = nrow(mat), ncol = ncol(mat)))
  }
  nm <- apply(mat, 2, function(col) {
    m <- mean(col, na.rm = TRUE)
    s <- sd(col, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      rep(0, length(col))
    } else {
      (col - m) / s
    }
  })
  if (is.null(dim(nm))) {
    nm <- matrix(nm, nrow = nrow(mat), ncol = ncol(mat))
  }
  return(nm)
}

#' @keywords internal
#' @noRd
mean_matrix <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(rep(0, ncol(mat)))
  }
  apply(mat, 2, function(col) mean(col, na.rm = TRUE))
}

#' @keywords internal
#' @noRd
sd_matrix <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(rep(1, ncol(mat)))
  }
  apply(mat, 2, function(col) {
    s <- sd(col, na.rm = TRUE)
    if (is.na(s) || s == 0) 1 else s
  })
}

#' @keywords internal
#' @noRd
normalize_matrix_known_mean_variance <- function(mat, mean, sd) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(matrix(numeric(0), nrow = nrow(mat), ncol = ncol(mat)))
  }
  nm <- sapply(1:NCOL(mat), function(col) {
    (mat[, col] - mean[col]) / sd[col]
  })
  if (is.null(dim(nm))) {
    nm <- matrix(nm, nrow = nrow(mat), ncol = ncol(mat))
  }
  return(nm)
}

#' @keywords internal
#' @noRd
normalize_nested_list <- function(nested_list) {
  lapply(nested_list, function(subgroup) lapply(subgroup, normalize_matrix))
}

#' @keywords internal
#' @noRd
mean_nested_list <- function(nested_list) {
  lapply(nested_list, function(subgroup) lapply(subgroup, mean_matrix))
}

#' @keywords internal
#' @noRd
sd_nested_list <- function(nested_list) {
  lapply(nested_list, function(subgroup) lapply(subgroup, sd_matrix))
}

#' @keywords internal
#' @noRd
normalize_nested_list_known_mean_sd <- function(nested_list, mean_nested_list, sd_nested_list) {
  mapply(function(x, y, z) {
    mapply(normalize_matrix_known_mean_variance, x, y, z, SIMPLIFY = FALSE)
  }, nested_list, mean_nested_list, sd_nested_list, SIMPLIFY = FALSE)
}
