# Post-run acceptance audit; never changes scientific results or references.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 2L))
root <- normalizePath(args[1], mustWork = TRUE)
require_original <- length(args) == 2L && identical(args[2], "--require-original")
if (length(args) == 2L && !require_original) stop("Unknown option")
s <- readRDS(file.path(root, "settings.rds"))
stopifnot(s$experiment == "table1", s$replicates == 1L)
if (require_original) stopifnot(!s$quick, s$draws == 350000L,
 s$burnin == 50000L, s$k == 10L, s$rounds == 10L)
# Verify provenance before loading any large fitted object.
stopifnot(identical(tools::md5sum(names(s$source_hashes)), s$source_hashes),
 identical(tools::md5sum(names(s$data_hash)), s$data_hash))
helper <- names(s$source_hashes)[basename(names(s$source_hashes)) == "original-experiment-helpers.R"]
stopifnot(length(helper) == 1L)
source(helper)
library(IntegMultiReg)
job <- file.path(root, "configuration-01-replicate-001")
data <- readRDS(file.path(job, "data.rds"))$data
cohort <- imr_data(data$platforms, data$outcome, covariates = data$covariates,
 outcome_type = "right.censored")$availability
labels <- c("imr-clinical-molecular", "bms-clinical-molecular",
 "imr-molecular", "bms-molecular", "cph-clinical")
status <- lapply(labels, function(label) readRDS(file.path(job, paste0(label,
  ".rds.status.rds"))))
stopifnot(all(vapply(status, function(x) identical(x$status, "completed"),
                     logical(1))))
acceptance_rows <- list()
first_folds <- NULL
id_round_key <- function(x) paste(x$id, x$round, sep = ":")
fold_cell_key <- function(x) paste(x$round, x$fold, x$subgroup, sep = ":")
for (j in seq_along(labels)) {
  label <- labels[j]
  path <- file.path(job, paste0(label, ".rds"))
  result <- read_experiment_result(path) # one at a time; legacy RDS may be large
  is_bayesian_job <- j <= 4L
  if (is_bayesian_job) {
    stopifnot(result$fit$control$mcmc$draws == s$draws,
      result$fit$control$mcmc$burnin == s$burnin,
      length(result$fit$posterior$selection_draws) == s$draws,
      identical(result$fit$control$method, sub("-.*", "", label)))
    for (name in names(s$reference_arguments$fit)) {
      actual <- if (name == "sampler_method") result$fit$control[[name]] else result$fit$control$numerical[[name]]
      stopifnot(identical(actual, s$reference_arguments$fit[[name]]))
    }
    stopifnot(result$fit$control$seed == s$seed,
      identical(result$fit$control$priors$nu, c(-4, -3, -4)),
      isTRUE(result$fit$control$standardize),
      length(result$fit$model$covariate_names) == if (j <= 2L) ncol(data$covariates) - 1L else 0L)
    stopifnot(result$cv$control$k == s$k, result$cv$control$rounds == s$rounds)
    for (name in setdiff(names(s$reference_arguments$cv), "cv_method"))
      stopifnot(identical(result$cv$control[[name]],
                          s$reference_arguments$cv[[name]]))
    stopifnot(identical(experiment_result_folds(path), result$cv$control$folds))
    pred <- result$cv$predictions
    score_method <- result$cv$control$score_method
    folds <- result$cv$control$folds
    if (is.null(first_folds)) first_folds <- folds
    result$fit <- NULL
    gc()
  } else {
    pred <- result$predictions
    score_method <- "standard"
    stopifnot(!is.null(first_folds))
    folds <- first_folds
  }
  stopifnot(!anyDuplicated(id_round_key(folds)))
  indices <- match(id_round_key(pred), id_round_key(folds))
  stopifnot(!anyNA(indices), identical(pred$fold, folds$fold[indices]))
  stopifnot(nrow(pred) == nrow(cohort) * s$rounds,
            all(is.finite(pred$prediction)),
   identical(sort(unique(pred$round)), seq_len(s$rounds)),
   identical(sort(unique(pred$fold)), seq_len(s$k)))
  for (round in seq_len(s$rounds)) {
    round_rows <- pred[pred$round == round, ]
    stopifnot(!anyDuplicated(round_rows$id), setequal(round_rows$id, cohort$id),
      identical(as.character(round_rows$subgroup),
        as.character(cohort$subgroup[match(round_rows$id, cohort$id)])))
  }
  # Recalculate with the recorded scoring convention, not another estimator.
  rescored <- fold_scores(pred, data$outcome, score_method)
  stopifnot(isTRUE(all.equal(result$scores, rescored, tolerance = 1e-12)))
  expected <- expand.grid(round = seq_len(s$rounds), fold = seq_len(s$k),
   subgroup = c(sort(unique(cohort$subgroup)), "all"), stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(fold_cell_key(rescored)),
    setequal(fold_cell_key(rescored), fold_cell_key(expected)))
  # Independent aggregation of the 100 validation sets; do not average rounds.
  # This deliberately re-derives summarize_scores() from the helpers by hand, so
  # folding it back into a call to that function would remove the independence.
  aggregated <- do.call(rbind, lapply(split(rescored, rescored$subgroup),
                                      function(rows) {
    values <- rows$value[is.finite(rows$value)]
    n <- length(values)
    data.frame(subgroup = rows$subgroup[1], mean = if (n) sum(values) / n else NA_real_,
      sd = if (n > 1) sqrt(sum((values - mean(values))^2) / (n - 1)) else NA_real_,
      se = if (n > 1) sqrt(sum((values - mean(values))^2) / (n - 1) / n) else NA_real_,
      n_valid = n, n_expected = s$k * s$rounds)
  }))
  rownames(aggregated) <- NULL
  recorded <- result$summary[match(aggregated$subgroup, result$summary$subgroup), ]
  rownames(recorded) <- NULL
  csv <- read.csv(sub("[.]rds$", "-scores.csv", path),
                  colClasses = c(subgroup = "character"))
  stopifnot(isTRUE(all.equal(aggregated, recorded, tolerance = 1e-12)),
   isTRUE(all.equal(recorded, csv, tolerance = 1e-12,
                    check.attributes = FALSE)))
  acceptance_rows[[j]] <- cbind(model = label, aggregated,
                                warnings = length(status[[j]]$warnings))
  rm(result, pred, rescored)
  gc()
  cat("PASS:", label, "\n")
}
acceptance <- do.call(rbind, acceptance_rows)
write.csv(acceptance, file.path(root, "acceptance-audit.csv"),
          row.names = FALSE)
audit_script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
saveRDS(list(settings = s, report = acceptance, session = sessionInfo(),
             verified_at = Sys.time(),
 audit_source = tools::md5sum(audit_script),
 scope = "settings, provenance, retained draw count, subject/fold coverage, scoring replay and independent aggregation; not convergence or historical digit agreement"),
 file.path(root, "acceptance-audit.rds"))
cat("PASS: five jobs audited. Non-comparable folds remain visible in n_valid.\n")
