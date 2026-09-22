# Resolve public choices once. Numerical engines receive explicit settings;
# the validation mode supplies defaults, not hidden overrides downstream.
.imr_cv_settings <- function(cv_method, ridge = NULL, model_set = NULL,
                              df_method = NULL, score_method = NULL, fold_rng = NULL) {
  choice <- function(value, choices, arg) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !value %in% choices) {
      .imr_abort(sprintf("`%s` must be one of %s.", arg,
                         paste(sprintf('"%s"', choices), collapse = ", ")))
    }
    value
  }
  if (cv_method == "refit") {
    supplied <- c(ridge = !is.null(ridge), model_set = !is.null(model_set),
                  df_method = !is.null(df_method), fold_rng = !is.null(fold_rng))
    if (any(supplied)) .imr_abort(sprintf(
      "`%s` applies only to post-fit CV; it cannot be used with `cv_method = \"refit\"`.",
      names(supplied)[which(supplied)[1L]]))
  } else {
    if (is.null(ridge)) ridge <- 0.001
    ridge <- .imr_check_numeric_vector(ridge, "ridge", length = 1L,
                                        nonnegative = TRUE)
    model_set <- choice(model_set %||% if (cv_method == "legacy")
      "ranked_unique" else "draws", c("draws", "ranked_unique"), "model_set")
    df_method <- choice(df_method %||% if (cv_method == "legacy")
      "legacy_integer" else "fractional",
      c("fractional", "legacy_integer"), "df_method")
    fold_rng <- choice(fold_rng %||% "reset", c("reset", "continue"), "fold_rng")
  }
  score_method <- choice(score_method %||% if (cv_method == "legacy")
    "legacy" else "standard", c("standard", "legacy"), "score_method")
  list(ridge = ridge, model_set = model_set, df_method = df_method,
       score_method = score_method, fold_rng = fold_rng)
}

.imr_cv_rng_state <- function(object, fold_rng) {
  if (!identical(fold_rng, "continue")) return(NULL)
  state <- object$control$rng_state
  if (!is.list(state) || !identical(state$generator, "gsl_rand48") ||
      !identical(state$endian, .Platform$endian) || !is.raw(state$state))
    .imr_abort("`fold_rng = \"continue\"` requires a fit with a saved GSL state from a compatible platform; refit this object.")
  state$state
}

# Fold tables are keyed by subject ID, never by their input row positions.
# row_order describes a permutation within each subgroup and round, retaining
# historical training/test summation order when a native plan is replayed.
.imr_cv_validate_folds <- function(folds, object, k, rounds) {
  if (!is.data.frame(folds) || !nrow(folds) ||
      anyDuplicated(names(folds)) || !all(c("id", "round", "fold") %in% names(folds)))
    .imr_abort("`folds` must be a nonempty data frame with unique columns `id`, `round`, and `fold`.")
  for (name in c("round", "fold")) {
    value <- folds[[name]]
    if (!.imr_is_integerish(value) || anyNA(value) ||
        any(value < 1 | value > .Machine$integer.max))
      .imr_abort(sprintf("`folds$%s` must contain positive finite integers.", name))
    folds[[name]] <- as.integer(value)
  }
  if (is.null(k)) k <- max(folds$fold)
  if (is.null(rounds)) rounds <- max(folds$round)
  k <- .imr_check_integer_scalar(k, "k", min = 2)
  rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  if (k != max(folds$fold) || rounds != max(folds$round) ||
      length(unique(folds$round)) != rounds)
    .imr_abort("`folds` must agree with `k` and `rounds`, with consecutive round and fold labels starting at 1.")
  dat <- object$preprocessing$input_data
  subjects <- dat$availability[dat$availability$subgroup %in%
    object$model$subgroup_names, c("id", "subgroup"), drop = FALSE]
  if (!is.atomic(folds$id) || anyNA(folds$id) ||
      any(!folds$id %in% subjects$id) ||
      nrow(folds) != as.double(nrow(subjects)) * rounds ||
      anyDuplicated(folds[c("id", "round")]))
    .imr_abort("`folds` must contain each modelled subject ID exactly once per round, with no missing or extra IDs.")
  subject_index <- match(folds$id, subjects$id)
  subgroup <- subjects$subgroup[subject_index]
  if (!"row_order" %in% names(folds)) {
    order_in_group <- stats::ave(seq_len(nrow(subjects)), subjects$subgroup, FUN = seq_along)
    folds$row_order <- as.integer(order_in_group[subject_index])
  }
  if (!.imr_is_integerish(folds$row_order) || anyNA(folds$row_order) ||
      any(folds$row_order < 1 | folds$row_order > .Machine$integer.max))
    .imr_abort("`folds$row_order` must be a permutation within each subgroup and round.")
  folds$row_order <- as.integer(folds$row_order)
  for (round in seq_len(rounds)) for (group in object$model$subgroup_names) {
    index <- which(folds$round == round & subgroup == group)
    if (!identical(sort(folds$row_order[index]), seq_along(index)))
      .imr_abort("`folds$row_order` must be a permutation within each subgroup and round.")
    counts <- tabulate(folds$fold[index], nbins = k)
    if (any(counts == 0L) || any(counts == length(index)))
      .imr_abort("`folds` must provide nonempty training and test rows in every subgroup and fold.")
  }
  list(k = k, rounds = rounds,
       folds = folds[c("id", "round", "fold", "row_order")])
}

.imr_cv_fold_matrices <- function(folds, ids, rounds) {
  if (is.null(folds)) return(list(folds = NULL, row_order = NULL))
  columns <- lapply(seq_len(rounds), function(round) {
    rows <- folds[folds$round == round, , drop = FALSE]
    rows[match(ids, rows$id), , drop = FALSE]
  })
  list(folds = vapply(columns, `[[`, integer(length(ids)), "fold"),
       row_order = vapply(columns, `[[`, integer(length(ids)), "row_order"))
}

.imr_cv_control <- function(settings, object, k, rounds, max_models, folds,
                             fold_source) {
  c(settings, list(k = k, rounds = rounds,
    max_models = if (identical(settings$model_set, "draws")) NULL else max_models,
    seed = object$control$seed,
    sampler_method = object$control$sampler_method %||% "legacy",
    initial = object$control$initial,
    standardize = object$control$standardize %||% TRUE,
    numerical = .imr_fit_numerical_control(object$control), folds = folds, fold_source = fold_source,
    package_version = as.character(utils::packageVersion("IntegMultiReg"))))
}
