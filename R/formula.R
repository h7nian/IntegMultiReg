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
