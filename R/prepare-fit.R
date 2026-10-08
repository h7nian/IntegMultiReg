# Prepare observed data once; no fitting or random draws occur here.
.imr_prepare_fit <- function(validated, min_subgroup_size, standardize, survival_scale) {
  platforms <- validated$platforms
  outcome <- validated$outcome
  covariates <- validated$covariates
  outcome_type <- validated$outcome_type
  n_platforms <- length(platforms)
  platform_names <- names(platforms)
  if (is.null(platform_names) || any(platform_names == "")) {
    platform_names <- paste0("platform", seq_len(n_platforms))
  }
  feature_names <- lapply(platforms, function(platform) colnames(platform)[-1])

  dat <- .imr_subgroup_data(outcome, covariates, platforms)
  # For each subgroup, extract the id vector from the first nonempty platform.
  dat_ids <- lapply(dat[[3]], function(subgroup) {
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
  subgroup_index <- as.integer(which(subgroup_rows > min_subgroup_size))
  if (length(subgroup_index) == 0) {
    .imr_abort(
      "No availability subgroup has more than `min_subgroup_size` subjects; lower it or check the data."
    )
  }
  sample_size <- as.integer(subgroup_rows[subgroup_index])

  # Filter the platform data to include only subgroups meeting the threshold.
  dat_filtered <- lapply(dat, function(x) x[subgroup_index])

  # Also filter the original id vectors accordingly.
  subgroup_ids_filtered <- dat_ids[subgroup_index]

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
  platform_field <- function(field) {
    lapply(platform_preprocessing, function(platforms) {
      lapply(platforms, `[[`, field)
    })
  }
  mean_train <- platform_field("mean")
  sd_train <- platform_field("sd")
  dat_normalized[[3]] <- platform_field("normalized")
  if (!is.null(covariates)) {
    covariate_preprocessing <- lapply(dat_filtered[[2]], .imr_prepare_matrix,
      standardize = standardize
    )
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


  subgroup_platforms <- lapply(dat_normalized[[3]], function(group) as.integer(which(vapply(group, nrow, 1L) > 0L)))
  platform_subgroups <- lapply(seq_len(n_platforms), function(l) {
    as.integer(which(vapply(subgroup_platforms, function(index) l %in% index, TRUE)))
  })
  .imr_check_mrf_capacity(platform_subgroups)
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


  subgroup_names <- names(x_filtered)
  model <- list(
    n_platforms = n_platforms, platform_names = platform_names,
    feature_names = feature_names,
    covariate_names = if (is.null(covariates)) character() else names(covariates)[-1L],
    subgroup_names = subgroup_names, sample_sizes = stats::setNames(sample_size, subgroup_names),
    subgroup_platforms = subgroup_platforms,
    platform_subgroups = platform_subgroups
  )
  prep <- list(
    input_data = validated, features = x_filtered, response = y_list,
    covariates = cov_list, subject_ids = stats::setNames(subgroup_ids_filtered, subgroup_names),
    feature_center = mean_train, feature_scale = sd_train,
    covariate_center = mean_cov_train, covariate_scale = sd_cov_train,
    formula = NULL, formula_data = NULL, terms = NULL, contrasts = NULL,
    xlevels = NULL, id = "id"
  )
  list(model = model, preprocessing = prep)
}
