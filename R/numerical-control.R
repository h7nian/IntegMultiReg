# These settings describe coefficient-mode approximation, not MCMC length.
.imr_numerical_control <- function(prior_indexing = "standard",
    laplace_max_iter = c(initial = 25L, selection = 40L, latent = 25L, prediction = 40L),
    laplace_tolerance = 1e-3) {
  if (!is.character(prior_indexing) || length(prior_indexing) != 1L ||
      is.na(prior_indexing) || !prior_indexing %in% c("standard", "code2017"))
    .imr_abort("`prior_indexing` must be \"standard\" or \"code2017\".")
  stages <- c("initial", "selection", "latent", "prediction")
  if (length(laplace_max_iter) == 1L && is.null(names(laplace_max_iter)))
    laplace_max_iter <- stats::setNames(rep(laplace_max_iter, 4L), stages)
  if (!.imr_is_integerish(laplace_max_iter) || length(laplace_max_iter) != 4L ||
      !identical(names(laplace_max_iter), stages) ||
      any(laplace_max_iter < 1 | laplace_max_iter > .Machine$integer.max))
    .imr_abort("`laplace_max_iter` must be a positive integer or a named integer vector: initial, selection, latent, prediction.")
  laplace_tolerance <- .imr_check_numeric_vector(laplace_tolerance,
    "laplace_tolerance", length = 1L, positive = TRUE)
  list(prior_indexing = prior_indexing,
       laplace_max_iter = stats::setNames(as.integer(laplace_max_iter), stages),
       laplace_tolerance = laplace_tolerance)
}

.imr_fit_numerical_control <- function(control) {
  if (is.null(control$numerical)) return(.imr_numerical_control())
  if (!is.list(control$numerical) ||
      !identical(names(control$numerical),
        c("prior_indexing", "laplace_max_iter", "laplace_tolerance")))
    .imr_abort("The fit has invalid numerical controls.")
  do.call(.imr_numerical_control, control$numerical)
}
