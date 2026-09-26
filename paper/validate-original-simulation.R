# Audit the executed simulation manifest; a one-replicate run is not a full study.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 1:2)
root <- normalizePath(args[1])
full_study <- length(args) == 2L && identical(args[2], '--require-complete-study')
if (length(args) == 2L && !full_study) stop('Unknown option')
s <- readRDS(file.path(root, 'settings.rds'))
stopifnot(s$experiment %in% c('simulation', 'correlated'))
if (full_study) stopifnot(!s$quick, s$draws == 350000L, s$burnin == 50000L,
 s$replicates == if (s$experiment == 'simulation') 50L else 30L)
stopifnot(identical(tools::md5sum(names(s$source_hashes)), s$source_hashes),
 identical(tools::md5sum(names(s$data_hash)), s$data_hash))
source(names(s$source_hashes)[basename(names(s$source_hashes)) == 'original-experiment-helpers.R'])
# Deliberately recomputes selection_auc() from the helpers by a different route:
# the pairwise probability that a positive exceeds a negative, rather than the
# rank form. Keeping the second derivation is what makes the check independent.
independent_pairwise_auc <- function(score, truth) {
  a <- score[truth == 1]
  b <- score[truth == 0]
  if (!length(a) || !length(b)) return(NA_real_)
  cmp <- outer(a, b, '-')
  mean((cmp > 0) + .5 * (cmp == 0))
}
acceptance_rows <- list()
configs <- if (s$experiment == 'simulation') 3L else 6L
task_selection <- list(configurations = seq_len(configs), replicates = seq_len(s$replicates))
if (!full_study && file.exists(file.path(root, "task.rds"))) {
  task_selection <- readRDS(file.path(root, "task.rds"))
  stopifnot(identical(task_selection$kind, s$experiment),
    all(task_selection$configurations %in% seq_len(configs)),
    all(task_selection$replicates %in% seq_len(s$replicates)))
}
id_round_key <- function(x) paste(x$id, x$round, sep = ':')
for (config in task_selection$configurations) for (replicate_index in task_selection$replicates) {
  job <- file.path(root, sprintf('configuration-%02d-replicate-%03d', config, replicate_index))
  data <- readRDS(file.path(job, 'data.rds'))$data
  # The unpenalized baseline is scored on the folds the IMR fit used, so that
  # plan is captured here and must come from this replicate, not a previous one.
  first_folds <- NULL
  for (method in c('imr-molecular', 'bms-molecular', 'l1-cph')) {
    path <- file.path(job, paste0(method, '.rds'))
    status <- readRDS(paste0(path, '.status.rds'))
    stopifnot(identical(status$status, 'completed'))
    result <- read_experiment_result(path)
    if (method != 'l1-cph') {
      fit <- result$fit
      cv <- result$cv
      stopifnot(identical(cv$validation, s$reference_arguments$cv$cv_method))
      stopifnot(fit$control$mcmc$draws == s$draws, fit$control$mcmc$burnin == s$burnin,
        length(fit$posterior$selection_draws) == s$draws,
        fit$control$seed == s$seed + replicate_index - 1L,
        cv$control$k == s$k, cv$control$rounds == s$rounds,
        identical(cv$control$folds, experiment_result_folds(path)))
      for (name in names(s$reference_arguments$fit)) {
        actual <- if (name == 'sampler_method') fit$control[[name]] else fit$control$numerical[[name]]
        stopifnot(identical(actual, s$reference_arguments$fit[[name]]))
      }
      for (name in setdiff(names(s$reference_arguments$cv), 'cv_method'))
        if (!(name == 'fold_rng' && method == 'bms-molecular' &&
              !identical(s$reference_arguments$cv$fold_rng, 'continue')))
          stopifnot(identical(cv$control[[name]], s$reference_arguments$cv[[name]]))
      selection_summary <- result$selection
      for (platform in seq_along(data$truth)) {
        groups <- fit$model$subgroup_names[fit$model$platform_subgroups[[platform]]]
        truth <- data$truth[[platform]][groups, , drop = FALSE]
        prob <- fit$posterior$inclusion_probabilities[[platform]]
        recomputed <- vapply(seq_along(groups),
          function(i) independent_pairwise_auc(prob[i, ], truth[i, ]), numeric(1))
        stored <- selection_summary$auc[match(paste(fit$model$platform_names[platform], groups),
          paste(selection_summary$platform, selection_summary$subgroup))]
        stopifnot(isTRUE(all.equal(unname(stored), unname(recomputed), tolerance = 1e-12)))
      }
      pred <- cv$predictions
      folds <- cv$control$folds
      score_method <- cv$control$score_method
      if (method == 'imr-molecular') first_folds <- folds
      result$fit <- NULL
      rm(fit, cv)
      gc()
    } else {
      pred <- result$predictions
      stopifnot(!is.null(first_folds))
      folds <- first_folds
      score_method <- 'standard'
    }
    stopifnot(all(is.finite(pred$prediction)), nrow(pred) == nrow(data$outcome) * s$rounds)
    stopifnot(!anyDuplicated(id_round_key(pred)),
      setequal(id_round_key(pred), id_round_key(folds)),
      identical(pred$fold, folds$fold[match(id_round_key(pred), id_round_key(folds))]))
    for (round in seq_len(s$rounds))
      stopifnot(setequal(pred$id[pred$round == round], data$outcome$id))
    calculated <- fold_scores(pred, data$outcome, score_method, weighted_overall = method == 'l1-cph')
    stopifnot(isTRUE(all.equal(result$scores, calculated, tolerance = 1e-12)),
      isTRUE(all.equal(result$summary, summarize_scores(calculated), tolerance = 1e-12)))
    acceptance_rows[[length(acceptance_rows) + 1L]] <- data.frame(config,
      rep = replicate_index, method, draws = s$draws,
      rounds = s$rounds, k = s$k, predictions = nrow(pred), warnings = length(status$warnings))
    cat('PASS:', config, replicate_index, method, '\n')
    rm(result)
    gc()
  }
  path <- file.path(job, 'uni-cph-selection.csv')
  stopifnot(readRDS(paste0(path, '.status.rds'))$status == 'completed')
  uni <- read.csv(path, colClasses = c(subgroup = 'character'))
  stopifnot(all(is.finite(uni$auc)), all(uni$auc >= 0 & uni$auc <= 1))
}
acceptance <- do.call(rbind, acceptance_rows)
expected <- expand.grid(config = task_selection$configurations, rep = task_selection$replicates,
 method = c('imr-molecular', 'bms-molecular', 'l1-cph'), stringsAsFactors = FALSE)
acceptance_key <- function(x) paste(x$config, x$rep, x$method, sep = ':')
stopifnot(nrow(acceptance) == nrow(expected), !anyDuplicated(acceptance_key(acceptance)),
 setequal(acceptance_key(acceptance), acceptance_key(expected)))
write.csv(acceptance, file.path(root, 'simulation-acceptance.csv'), row.names = FALSE)
validator_path <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
saveRDS(list(settings = s, report = acceptance, full_study = full_study, verified_at = Sys.time(),
 validator_hash = tools::md5sum(validator_path),
 scope = 'Executed manifest, fit settings, pairwise IMR/BMS selection AUC, folds and predictive scores. Univariate/L1 scientific stability needs separate warning review.'),
 file.path(root, 'simulation-acceptance.rds'))
cat('PASS executed manifest. Complete-study requirement:', full_study, '\n')
