# Reproducible predictive comparison and selection of prespecified formulas.
# Rscript IntegMultiReg-covariate-comparison.R [--quick] [--out-dir DIRECTORY]
# Uses only public IntegMultiReg APIs. Synthetic covariate names are illustrative,
# not measurements from a clinical study. This is prediction, not confounder selection.
run_covariate_comparison <- function(out_dir = "covariate-comparison", quick = FALSE) {
  stopifnot(is.logical(quick), length(quick) == 1L, !is.na(quick))
  if (utils::packageVersion("IntegMultiReg") != "0.1.4") stop("Requires IntegMultiReg 0.1.4")
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) saved_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_rng) assign(".Random.seed", saved_rng, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  env <- new.env(); utils::data("simIMR", package = "IntegMultiReg", envir = env)
  dat <- env$simIMR
  clinical <- merge(dat$outcome.continuous, dat$covariates, by = "id", sort = TRUE)
  input <- IntegMultiReg::imr_data(dat$platforms, dat$outcome.continuous,
    covariates = dat$covariates, type_outcome = "continuous")
  cohort <- as.data.frame(input)
  sizes <- table(cohort$subgroup)
  cohort <- cohort[cohort$subgroup %in% names(sizes)[sizes > 30], ]
  clinical <- clinical[match(cohort$id, clinical$id), ]
  stopifnot(!anyNA(clinical), identical(clinical$id, cohort$id))
  # Both candidate formulas and all settings are fixed before evaluation.
  formulas <- list(age_only = y ~ age, age_sex_stage = y ~ age + sex + stage)
  draws <- if (quick) c(80, 40) else c(2000, 500)
  k <- if (quick) 2L else 3L
  rounds <- if (quick) 1L else 2L
  assays <- function(ids) lapply(dat$platforms, function(x) x[x$id %in% ids, , drop = FALSE])
  fit_candidate <- function(formula, ids, seed) IntegMultiReg::imr(
    formula, data = clinical[match(ids, clinical$id), , drop = FALSE],
    platforms = assays(ids), type_outcome = "continuous", ssize = 0,
    nu = c(-4, -3, -4), sample_mcmc = draws, seed = seed)
  fit_set <- function(ids, seed) lapply(formulas, fit_candidate, ids = ids, seed = seed)
  evaluate <- function(fits, rounds) {
    cv <- lapply(fits, IntegMultiReg::cv_imr, k = k, rounds = rounds)
    keys <- lapply(cv, function(x) attr(x, "predictions")[, c("round", "id", "subgroup", "fold")])
    # Identical seeds alone are not evidence of paired folds: check the actual assignments.
    stopifnot(identical(keys[[1]], keys[[2]]))
    cv
  }
  fits <- fit_set(cohort$id, 24019L)
  cv <- evaluate(fits, rounds)
  metrics <- lapply(cv, function(x) x$total_cindex[, "all"])
  paired_summary <- data.frame(candidate = names(formulas),
    formula = vapply(formulas, function(f) paste(deparse(f), collapse = ""), ""),
    mean_mse = vapply(metrics, mean, 0),
    sd_across_rounds = vapply(metrics, function(x) if (length(x)>1) stats::sd(x) else NA_real_, 0),
    row.names = NULL)
  paired_rounds <- data.frame(round = seq_len(rounds), age_only = metrics[[1]],
    age_sex_stage = metrics[[2]], expanded_minus_age_only = metrics[[2]] - metrics[[1]])
  metadata <- IntegMultiReg::compare_imr(fits)

  # Outer splits use availability groups only; held-out responses never choose a formula.
  set.seed(81043L)
  outer <- integer(nrow(cohort))
  for (group in unique(cohort$subgroup)) {
    idx <- which(cohort$subgroup == group)
    outer[idx[sample.int(length(idx))]] <- rep(seq_len(k), length.out = length(idx))
  }
  prediction <- rep(NA_real_, nrow(cohort))
  selections <- inner_records <- vector("list", k)
  for (fold in seq_len(k)) {
    heldout <- which(outer == fold); train <- which(outer != fold)
    candidate_fits <- fit_set(cohort$id[train], 92000L + fold)
    inner <- evaluate(candidate_fits, 1L)
    scores <- vapply(inner, function(x) unname(x$total_cindex[1, "all"]), 0)
    stopifnot(all(is.finite(scores)))
    chosen <- which.min(scores)  # deterministic candidate-order tie break
    records <- attr(inner[[1]], "predictions")
    stopifnot(setequal(records$id, cohort$id[train]),
      !any(records$id %in% cohort$id[heldout]))
    inner_records[[fold]] <- data.frame(outer_fold = fold,
      records[, c("round", "id", "subgroup", "fold")])
    p <- do.call(rbind, stats::predict(candidate_fits[[chosen]],
      newdata = assays(cohort$id[heldout]),
      covariates = clinical[match(cohort$id[heldout], clinical$id), , drop = FALSE]))
    stopifnot(setequal(p$id, cohort$id[heldout]), !anyDuplicated(p$id))
    prediction[heldout] <- p$predict[match(cohort$id[heldout], p$id)]
    selections[[fold]] <- data.frame(outer_fold = fold,
      selected = names(formulas)[chosen], n_train = length(train), n_test = length(heldout),
      inner_mse_age_only = scores[1], inner_mse_age_sex_stage = scores[2],
      outer_mse = mean((prediction[heldout] - clinical$y[heldout])^2), row.names = NULL)
  }
  stopifnot(all(is.finite(prediction)), all(table(outer) > 0))
  nested_summary <- data.frame(procedure = "inner-CV formula selection",
    n_subjects = nrow(cohort), outer_folds = k, inner_folds = k,
    pooled_outer_mse = mean((prediction - clinical$y)^2))
  result <- list(paired_summary = paired_summary, paired_rounds = paired_rounds,
    model_metadata = metadata, paired_folds = attr(cv[[1]], "predictions")[, c("round", "id", "subgroup", "fold")],
    nested_summary = nested_summary, selected_by_fold = do.call(rbind, selections),
    nested_inner_folds = do.call(rbind, inner_records),
    outer_predictions = data.frame(cohort, outer_fold = outer, observed = clinical$y, prediction),
    settings = list(quick = quick, sample_mcmc = draws, k = k, rounds = rounds,
      formulas = formulas, version = as.character(utils::packageVersion("IntegMultiReg"))))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result)[vapply(result, is.data.frame, TRUE)])
    utils::write.csv(result[[name]], file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  saveRDS(result, file.path(out_dir, "comparison.rds"))
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  cat("\nPrespecified candidates: paired CV (MSE; lower is better)\n")
  print(metadata, row.names = FALSE); print(paired_summary, row.names = FALSE)
  cat("\nPaired round differences (not independent inferential replicates)\n")
  print(paired_rounds, row.names = FALSE)
  cat("\nNested selection: the outer responses are used only for scoring\n")
  print(result$selected_by_fold, row.names = FALSE); print(nested_summary, row.names = FALSE)
  invisible(result)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match("--out-dir", args)
  out <- if (!is.na(hit) && hit < length(args)) args[hit + 1L] else "covariate-comparison"
  run_covariate_comparison(out, quick = "--quick" %in% args)
}
