# Subject alignment and matrix preparation shared by fitting, prediction and CV.

#' @keywords internal
#' @noRd
.imr_subgroup_data <- function(outcome, covariates = NULL, platforms) {
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
.imr_prepare_matrix <- function(mat, standardize = TRUE) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(list(
      mean = rep(0, ncol(mat)), sd = rep(1, ncol(mat)),
      normalized = matrix(numeric(0), nrow(mat), ncol(mat))
    ))
  }
  mat <- as.matrix(mat)
  if (!standardize) {
    storage.mode(mat) <- "double"
    return(list(
      mean = stats::setNames(rep(0, ncol(mat)), colnames(mat)),
      sd = stats::setNames(rep(1, ncol(mat)), colnames(mat)), normalized = mat
    ))
  }
  centers <- scales <- numeric(ncol(mat))
  names(centers) <- names(scales) <- colnames(mat)
  # Extract columns directly, avoiding apply's full matrix permutation copy.
  # Keep moment arithmetic and the raw degenerate-scale decision unchanged.
  normalized_columns <- lapply(seq_len(ncol(mat)), function(column_index) {
    column <- mat[, column_index]
    center <- mean(column, na.rm = TRUE)
    scale <- sd(column, na.rm = TRUE)
    degenerate <- is.na(scale) || scale == 0
    centers[column_index] <<- center
    scales[column_index] <<- if (degenerate) 1 else scale
    if (degenerate) {
      rep(0, length(column))
    } else {
      (column - center) / scale
    }
  })
  names(normalized_columns) <- colnames(mat)
  normalized <- simplify2array(normalized_columns)
  # apply retains result row names only when every column returns the same names.
  if (nrow(mat) > 1L) {
    result_names <- names(normalized_columns[[1L]])
    if (!all(vapply(normalized_columns, function(column) {
      identical(names(column), result_names)
    }, logical(1)))) {
      result_names <- NULL
    }
    result_dimnames <- list(result_names, colnames(mat))
    dimension_names <- names(dimnames(mat))
    if (!is.null(dimension_names)) {
      names(result_dimnames) <- c(
        if (length(result_names) == nrow(mat)) dimension_names[1L] else "",
        dimension_names[2L]
      )
    }
    dimnames(normalized) <- if (is.null(dimension_names) &&
      all(vapply(result_dimnames, is.null, logical(1)))) {
      NULL
    } else {
      result_dimnames
    }
  }
  if (is.null(dim(normalized))) {
    normalized <- matrix(normalized, nrow = nrow(mat), ncol = ncol(mat))
  }
  list(mean = centers, sd = scales, normalized = normalized)
}

#' @keywords internal
#' @noRd
.imr_standardize_matrix <- function(mat, mean, sd) {
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
.imr_standardize_platforms <- function(nested_list, mean_nested_list, sd_nested_list) {
  mapply(function(x, y, z) {
    mapply(.imr_standardize_matrix, x, y, z, SIMPLIFY = FALSE)
  }, nested_list, mean_nested_list, sd_nested_list, SIMPLIFY = FALSE)
}

#' @keywords internal
#' @noRd
# Platform tables reach the sampler as plain double matrices with the leading
# "id" column dropped. Fitting and prediction take the same path, so the
# conversion lives in one place.
.imr_platform_matrices <- function(subgroups) {
  lapply(subgroups, function(subgroup) {
    lapply(subgroup, function(platform_data) {
      labelled <- (is.data.frame(platform_data) &&
        "id" %in% colnames(platform_data)) ||
        (!is.null(colnames(platform_data)) &&
          colnames(platform_data)[1] == "id")
      mat <- as.matrix(
        if (labelled) platform_data[, -1, drop = FALSE] else platform_data
      )
      storage.mode(mat) <- "double"
      mat
    })
  })
}
