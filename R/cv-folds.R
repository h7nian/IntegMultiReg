# Validate actual subject memberships; the same plan serves both CV algorithms.
.imr_cv_validate_folds <- function(folds, object, k, rounds) {
  if (!is.data.frame(folds) || !nrow(folds) ||
    anyDuplicated(names(folds)) || !all(c("id", "round", "fold") %in% names(folds))) {
    .imr_abort("`folds` must be a nonempty data frame with unique columns `id`, `round`, and `fold`.")
  }
  for (name in c("round", "fold")) {
    value <- folds[[name]]
    if (!.imr_is_integerish(value) || anyNA(value) ||
      any(value < 1 | value > .Machine$integer.max)) {
      .imr_abort(sprintf("`folds$%s` must contain positive finite integers.", name))
    }
    folds[[name]] <- as.integer(value)
  }
  if (is.null(k)) k <- max(folds$fold)
  if (is.null(rounds)) rounds <- max(folds$round)
  k <- .imr_check_integer_scalar(k, "k", min = 2)
  rounds <- .imr_check_integer_scalar(rounds, "rounds", min = 1)
  if (k != max(folds$fold) || rounds != max(folds$round) ||
    length(unique(folds$round)) != rounds) {
    .imr_abort("`folds` must agree with `k` and `rounds`, with consecutive round and fold labels starting at 1.")
  }
  dat <- object$preprocessing$input_data
  subjects <- dat$availability[dat$availability$subgroup %in%
    object$model$subgroup_names, c("id", "subgroup"), drop = FALSE]
  if (!is.atomic(folds$id) || anyNA(folds$id) ||
    any(!folds$id %in% subjects$id) ||
    nrow(folds) != as.double(nrow(subjects)) * rounds ||
    anyDuplicated(folds[c("id", "round")])) {
    .imr_abort("`folds` must contain each modelled subject ID exactly once per round, with no missing or extra IDs.")
  }
  subject_index <- match(folds$id, subjects$id)
  subgroup <- subjects$subgroup[subject_index]
  for (round in seq_len(rounds)) {
    for (group in object$model$subgroup_names) {
      index <- which(folds$round == round & subgroup == group)
      counts <- tabulate(folds$fold[index], nbins = k)
      if (any(counts == 0L) || any(counts == length(index))) {
        .imr_abort("`folds` must provide nonempty training and test rows in every subgroup and fold.")
      }
    }
  }
  list(
    k = k, rounds = rounds,
    folds = folds[c("id", "round", "fold")]
  )
}
