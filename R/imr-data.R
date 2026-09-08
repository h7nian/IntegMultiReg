#' Validated Multi-Platform Data for IMR
#'
#' `imr_data()` creates a validated data object shared by model fitting and
#' prediction. It standardizes the subject identifier to `id`, records the
#' observed-platform pattern for every subject and catches alignment problems
#' before the MCMC sampler is called.
#'
#' @param platforms A non-empty named list of data frames, one per platform.
#'   Every data frame must contain the subject identifier and at least one
#'   finite numeric feature.
#' @param outcome Optional outcome data frame. It is required when the object is
#'   passed to [imr()] and may be omitted for prediction data. For binary and
#'   continuous outcomes it contains the identifier and response; for
#'   right-censored outcomes it contains the identifier, time and status.
#' @param covariates Optional data frame containing the identifier followed by
#'   clinical covariates.
#' @param type_outcome Optional outcome type: `"binary"`, `"continuous"` or
#'   `"right.censored"`. It is required when `outcome` is supplied.
#' @param id Name of the subject-identifier column in every supplied data frame.
#'
#' @return An object of class `"imr_data"` with components `platforms`,
#'   `outcome`, `covariates`, `type_outcome`, `id`, `availability` and
#'   `subgroup_sizes`, `n_platform_subjects` and `excluded_ids`. Availability
#'   summaries include only subjects with all required outcome/covariate rows.
#' @export
#'
#' @examples
#' data("simIMR", package = "IntegMultiReg")
#' dat <- imr_data(
#'   platforms = simIMR$platforms,
#'   outcome = simIMR$outcome.binary,
#'   covariates = simIMR$covariates,
#'   type_outcome = "binary"
#' )
#' dat
imr_data <- function(platforms, outcome = NULL, covariates = NULL,
                     type_outcome = NULL, id = "id") {
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(id)) {
    .imr_abort("`id` must be one non-empty column name.")
  }
  if (!is.list(platforms) || length(platforms) == 0L) {
    .imr_abort("`platforms` must be a non-empty list of data frames.")
  }
  platform_names <- names(platforms)
  if (is.null(platform_names)) {
    platform_names <- paste0("platform", seq_along(platforms))
  } else {
    missing_names <- is.na(platform_names) | !nzchar(platform_names)
    platform_names[missing_names] <- paste0("platform", which(missing_names))
    if (anyDuplicated(platform_names)) {
      .imr_abort("Platform names must be unique after unnamed elements are labelled.")
    }
  }
  platforms <- lapply(seq_along(platforms), function(i) {
    .imr_standardize_id_frame(
      platforms[[i]], id, sprintf("platforms[[%d]]", i),
      require_features = TRUE
    )
  })
  names(platforms) <- platform_names

  if (!is.null(outcome)) {
    if (is.null(type_outcome)) {
      .imr_abort("`type_outcome` is required when `outcome` is supplied.")
    }
    type_outcome <- match.arg(
      type_outcome, c("right.censored", "binary", "continuous")
    )
    outcome <- .imr_standardize_id_frame(
      outcome, id, "outcome", require_features = TRUE
    )
    .imr_validate_outcome(outcome, type_outcome)
  } else if (!is.null(type_outcome)) {
    type_outcome <- match.arg(
      type_outcome, c("right.censored", "binary", "continuous")
    )
  }

  if (!is.null(covariates)) {
    covariates <- .imr_standardize_id_frame(
      covariates, id, "covariates", require_features = TRUE
    )
  }

  platform_ids <- unique(unlist(lapply(platforms, `[[`, "id"), use.names = FALSE))
  all_ids <- platform_ids
  if (!is.null(outcome)) all_ids <- intersect(all_ids, outcome$id)
  if (!is.null(covariates)) all_ids <- intersect(all_ids, covariates$id)
  if (length(all_ids) == 0L) {
    .imr_abort("No platform subjects have all required outcome/covariate rows.")
  }
  presence <- do.call(cbind, lapply(platforms, function(x) all_ids %in% x$id))
  colnames(presence) <- platform_names
  bitstrings <- apply(
    presence, 1L, function(x) paste(as.integer(rev(x)), collapse = "")
  )
  availability <- data.frame(
    id = all_ids,
    presence,
    subgroup = bitstrings,
    check.names = FALSE,
    row.names = NULL
  )
  subgroup_sizes <- sort(table(bitstrings))

  out <- list(
    platforms = platforms,
    outcome = outcome,
    covariates = covariates,
    type_outcome = type_outcome,
    id = id,
    availability = availability,
    subgroup_sizes = subgroup_sizes,
    n_platform_subjects = length(platform_ids),
    excluded_ids = setdiff(platform_ids, all_ids)
  )
  class(out) <- "imr_data"
  validate_imr_data(out)
  out
}

#' Validate an IMR Data Object
#'
#' Rechecks the structure and subject alignment of an object created by
#' [imr_data()]. Invalid objects fail with an informative error.
#'
#' @param x An `"imr_data"` object.
#' @return `TRUE`, invisibly.
#' @export
validate_imr_data <- function(x) {
  if (!inherits(x, "imr_data") || !is.list(x)) {
    .imr_abort("`x` must be an `imr_data` object.")
  }
  required <- c(
    "platforms", "outcome", "covariates", "type_outcome", "id",
    "availability", "subgroup_sizes", "n_platform_subjects", "excluded_ids"
  )
  if (!all(required %in% names(x))) {
    .imr_abort("The `imr_data` object is missing required components.")
  }
  if (!is.list(x$platforms) || length(x$platforms) == 0L) {
    .imr_abort("`x$platforms` must be a non-empty list.")
  }
  if (is.null(names(x$platforms)) || anyNA(names(x$platforms)) ||
      any(!nzchar(names(x$platforms))) ||
      anyDuplicated(names(x$platforms))) {
    .imr_abort("`x$platforms` must have complete, unique names.")
  }
  for (i in seq_along(x$platforms)) {
    arg <- sprintf("x$platforms[[%d]]", i)
    .imr_check_id_frame(x$platforms[[i]], arg)
    .imr_check_numeric_columns(x$platforms[[i]], arg)
  }
  if (!is.null(x$outcome)) {
    .imr_check_id_frame(x$outcome, "x$outcome")
    .imr_validate_outcome(x$outcome, x$type_outcome)
  }
  if (!is.null(x$covariates)) {
    .imr_check_id_frame(x$covariates, "x$covariates")
    .imr_check_numeric_columns(x$covariates, "x$covariates")
  }
  if (!is.data.frame(x$availability) ||
      !all(c("id", "subgroup") %in% names(x$availability))) {
    .imr_abort("`x$availability` is not a valid availability table.")
  }
  platform_ids <- unique(unlist(lapply(x$platforms, `[[`, "id"), use.names = FALSE))
  all_ids <- platform_ids
  if (!is.null(x$outcome)) all_ids <- intersect(all_ids, x$outcome$id)
  if (!is.null(x$covariates)) all_ids <- intersect(all_ids, x$covariates$id)
  presence <- do.call(cbind, lapply(x$platforms, function(z) all_ids %in% z$id))
  colnames(presence) <- names(x$platforms)
  bitstrings <- apply(
    presence, 1L, function(z) paste(as.integer(rev(z)), collapse = "")
  )
  expected_availability <- data.frame(
    id = all_ids, presence, subgroup = bitstrings,
    check.names = FALSE, row.names = NULL
  )
  if (!identical(x$availability, expected_availability) ||
      !identical(x$subgroup_sizes, sort(table(bitstrings))) ||
      !identical(x$n_platform_subjects, length(platform_ids)) ||
      !identical(x$excluded_ids, setdiff(platform_ids, all_ids))) {
    .imr_abort("Availability metadata does not match the platform data.")
  }
  invisible(TRUE)
}

#' Convert IMR Data to an Availability Data Frame
#'
#' @param x An `"imr_data"` object.
#' @param row.names Unused; present for compatibility with [as.data.frame()].
#' @param optional Unused; present for compatibility with [as.data.frame()].
#' @param ... Unused.
#' @return The subject-by-platform availability table, including bitstrings.
#' @export
as.data.frame.imr_data <- function(x, row.names = NULL, optional = FALSE, ...) {
  validate_imr_data(x)
  x$availability
}

#' @export
print.imr_data <- function(x, ...) {
  validate_imr_data(x)
  cat("Validated IMR multi-platform data\n")
  cat("---------------------------------\n")
  cat(sprintf(
    "Platforms : %d (%s)\n", length(x$platforms),
    paste(names(x$platforms), collapse = ", ")
  ))
  cat(sprintf(
    "Subjects  : %d eligible of %d with at least one platform",
    nrow(x$availability), x$n_platform_subjects
  ))
  if (length(x$excluded_ids)) {
    cat(sprintf(" (%d missing required outcome/covariate rows)",
                length(x$excluded_ids)))
  }
  cat("\n")
  cat(sprintf("Outcome   : %s\n", if (is.null(x$outcome)) {
    "not supplied (prediction data)"
  } else {
    x$type_outcome
  }))
  cat("Availability subgroups:\n")
  for (nm in names(x$subgroup_sizes)) {
    cat(sprintf("  %s : %d\n", nm, x$subgroup_sizes[[nm]]))
  }
  invisible(x)
}

#' @export
summary.imr_data <- function(object, ...) {
  validate_imr_data(object)
  feature_counts <- vapply(
    object$platforms, function(x) ncol(x) - 1L, integer(1L)
  )
  out <- list(
    n_subjects = nrow(object$availability),
    n_platform_subjects = object$n_platform_subjects,
    n_excluded = length(object$excluded_ids),
    n_platforms = length(object$platforms),
    feature_counts = feature_counts,
    subgroup_sizes = object$subgroup_sizes,
    has_outcome = !is.null(object$outcome),
    has_covariates = !is.null(object$covariates),
    type_outcome = object$type_outcome
  )
  class(out) <- "summary.imr_data"
  out
}

#' @export
print.summary.imr_data <- function(x, ...) {
  cat("IMR data summary\n")
  cat(sprintf("  Subjects: %d; platforms: %d\n", x$n_subjects, x$n_platforms))
  cat("  Features per platform:\n")
  print(x$feature_counts)
  cat("  Availability subgroups:\n")
  print(x$subgroup_sizes)
  invisible(x)
}

#' @keywords internal
#' @noRd
.imr_standardize_id_frame <- function(x, id, arg, require_features = TRUE) {
  if (!is.data.frame(x)) {
    .imr_abort(sprintf("`%s` must be a data frame.", arg))
  }
  .imr_check_column_names(x, arg)
  if (!identical(id, "id") && "id" %in% names(x)) {
    .imr_abort("A non-identifier column named `id` conflicts with the standardized identifier.")
  }
  if (!id %in% names(x)) {
    if (identical(id, "id")) {
      .imr_abort(sprintf("`%s` must have `id` as its first column.", arg))
    }
    .imr_abort(sprintf("`%s` must contain the identifier column `%s`.", arg, id))
  }
  x <- x[, c(id, setdiff(names(x), id)), drop = FALSE]
  names(x)[1L] <- "id"
  .imr_check_id_frame(x, arg, require_features = require_features)
  .imr_check_numeric_columns(x, arg)
  x
}

#' @keywords internal
#' @noRd
.imr_validate_outcome <- function(outcome, type_outcome) {
  if (is.null(type_outcome) || length(type_outcome) != 1L) {
    .imr_abort("A valid `type_outcome` is required for outcome validation.")
  }
  if (type_outcome %in% c("binary", "continuous")) {
    if (ncol(outcome) != 2L) {
      .imr_abort("`outcome` must have exactly two columns: `id` and the response.")
    }
    .imr_check_numeric_columns(outcome, "outcome", names(outcome)[2L])
    if (type_outcome == "binary" && !all(outcome[[2L]] %in% c(0, 1))) {
      .imr_abort("For `type_outcome = \"binary\"`, the response must be coded 0/1.")
    }
  } else if (type_outcome == "right.censored") {
    if (ncol(outcome) != 3L) {
      .imr_abort(paste0(
        "`outcome` must have exactly three columns for right-censored data: ",
        "`id`, time and status."
      ))
    }
    .imr_check_numeric_columns(outcome, "outcome", names(outcome)[2:3])
    if (any(outcome[[2L]] <= 0)) {
      .imr_abort("Right-censored event times in `outcome` must be positive.")
    }
    if (!all(outcome[[3L]] %in% c(0, 1))) {
      .imr_abort("Right-censored status values in `outcome` must be coded 0/1.")
    }
  } else {
    .imr_abort("Unknown `type_outcome`.")
  }
  invisible(outcome)
}
