.imr_cv_folds <- function(strata, k) {
  out <- integer(length(strata))
  offset <- 0L
  for (level in unique(strata)) {
    idx <- which(strata == level)
    idx <- idx[sample.int(length(idx))]
    out[idx] <- (offset + seq_along(idx) - 1L) %% k + 1L
    offset <- offset + length(idx)
  }
  out
}

.imr_cv_refit <- function(object, train_ids, seed, verbose = FALSE) {
  control <- object$control
  model <- object$model
  preprocessing <- object$preprocessing
  dat <- preprocessing$input_data
  # Keep rows outside the eligible cohort for otherwise unused platforms;
  # .imr_subgroup_data() intersects these with the training outcome IDs before any
  # normalization or sampling. No held-out outcome is passed to the refit.
  train_platforms <- lapply(dat$platforms, function(x) {
    keep <- x$id %in% train_ids | !x$id %in% dat$availability$id
    x[keep, , drop = FALSE]
  })
  # A platform occurring only in filtered-out subgroups still needs a valid
  # input frame, but none of its rows will enter a training subgroup.
  for (p in seq_along(train_platforms)) {
    if (nrow(train_platforms[[p]]) == 0L && length(model$platform_subgroups[[p]]) == 0L) {
      train_platforms[[p]] <- dat$platforms[[p]]
    }
  }
  mcmc <- control$mcmc
  mcmc$seed <- seed
  mcmc$workers <- 1L
  mcmc$keep_latent <- FALSE
  refit_control <- list(
    outcome_type = control$outcome_type, model_variant = control$model_variant,
    standardize = control$standardize, min_subgroup_size = 0L,
    priors = control$priors, mcmc = mcmc,
    survival_scale = if (control$outcome_type == "right.censored") control$response_scale else "identity",
    verbose = verbose
  )
  if (!is.null(preprocessing$terms)) {
    id <- preprocessing$id
    train_platforms <- lapply(train_platforms, function(x) {
      names(x)[1L] <- id
      x
    })
    return(do.call(imr, c(list(
      x = preprocessing$formula,
      data = preprocessing$formula_data[preprocessing$formula_data[[id]] %in% train_ids, , drop = FALSE],
      platforms = train_platforms, id = id
    ), refit_control)))
  }
  do.call(imr, c(list(
    x = train_platforms,
    outcome = .imr_match_rows(dat$outcome, train_ids, "outcome"),
    covariates = if (!is.null(dat$covariates)) .imr_match_rows(dat$covariates, train_ids, "covariates") else NULL
  ), refit_control))
}

.imr_cv_accuracy <- function(type, prediction, outcome) {
  if (!length(prediction)) {
    return(NA_real_)
  }
  y <- outcome[[2L]]
  if (type == "continuous") {
    return(mean((prediction - y)^2))
  }
  if (type == "binary") {
    positive <- sum(y == 1)
    negative <- sum(y == 0)
    if (!positive || !negative) {
      return(NA_real_)
    }
    return((sum(rank(prediction, ties.method = "average")[y == 1]) -
      positive * (positive + 1) / 2) / (positive * negative))
  }
  .Call("imr_concordance", as.double(prediction), as.double(y),
    as.integer(outcome[[3L]]),
    PACKAGE = "IntegMultiReg"
  )
}
