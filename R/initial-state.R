# Explicit starts use the same named platform/subgroup/feature layout as the fit.
.imr_initial_state <- function(initial, platforms, subgroups, features, method) {
  if (is.null(initial)) return(NULL)
  if (!is.list(initial) || !length(initial) || is.null(names(initial)) ||
      anyDuplicated(names(initial)) || any(!names(initial) %in% c("selection", "interaction")))
    .imr_abort("`initial` must be a named list containing selection and/or interaction.")
  if (method == "bms" && !is.null(initial$interaction))
    .imr_abort("`initial$interaction` is not applicable to BMS.")
  normalize <- function(values, interaction) {
    if (is.null(values)) return(NULL)
    if (!is.list(values) || length(values) != length(platforms) ||
        is.null(names(values)) || anyDuplicated(names(values)) || !setequal(names(values), platforms))
      .imr_abort("Initial matrices must be named by every platform.")
    values <- values[platforms]
    for (p in seq_along(platforms)) {
      x <- values[[p]]
      rows <- subgroups[[p]]
      columns <- if (interaction) rows else features[[p]]
      if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || anyNA(x) ||
          any(!is.finite(x)) || !identical(dim(x), c(length(rows), length(columns))) ||
          (length(rows) && !identical(rownames(x), rows)) ||
          (length(columns) && !identical(colnames(x), columns)))
        .imr_abort("Initial matrices must have finite values and matching subgroup/feature dimnames.")
      if (interaction) {
        if (!isTRUE(all.equal(x, t(x), tolerance = 0)) || any(diag(x) != 0) ||
            any(x[row(x) != col(x)] <= 0))
          .imr_abort("Initial interactions must be symmetric with zero diagonal and positive off-diagonal entries.")
        storage.mode(x) <- "double"
      } else {
        if (any(!x %in% 0:1)) .imr_abort("Initial selection entries must be zero or one.")
        storage.mode(x) <- "integer"
      }
      values[[p]] <- x
    }
    values
  }
  list(selection = normalize(initial$selection, FALSE),
       interaction = normalize(initial$interaction, TRUE))
}
