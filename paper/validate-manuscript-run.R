# Inspect the scientific outputs, not merely the replication completion log.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 2L))
library(IntegMultiReg)
library(survival)
id_round_key <- function(rows) paste(rows$id, rows$round, sep = ':')
validate_run <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  result <- readRDS(file.path(root, 'replication_results.rds'))
  stopifnot(!result$quick, result$sampler_method == 'paper', result$cv_method == 'refit')
  outcomes <- list(binary = simIMR$outcome.binary, continuous = simIMR$outcome.continuous,
   right.censored = simIMR$outcome.survival)
  acceptance_rows <- list()
  for (type in names(outcomes)) for (method in c('imr', 'bms')) {
    saved <- readRDS(file.path(root, paste0('validation-', type, '-', method, '.rds')))
    stopifnot(saved$fit$control$mcmc$draws == 4000L, saved$fit$control$mcmc$burnin == 1000L,
     saved$fit$control$sampler_method == 'paper', saved$fit$control$method == method,
     saved$cv$validation == 'refit', saved$cv$control$k == 5L, saved$cv$control$rounds == 10L,
     saved$cv$control$sampler_method == 'paper', saved$cv$control$score_method == 'standard')
    predictions <- saved$cv$predictions
    folds <- saved$cv$control$folds
    stopifnot(nrow(predictions) == 3000L, all(is.finite(predictions$prediction)))
    stopifnot(!anyDuplicated(id_round_key(predictions)), !anyDuplicated(id_round_key(folds)),
     setequal(id_round_key(predictions), id_round_key(folds)),
     identical(predictions$fold, folds$fold[match(id_round_key(predictions), id_round_key(folds))]))
    outcome_table <- outcomes[[type]]
    calculate <- function(scored_rows) {
      obs <- outcome_table[match(scored_rows$id, outcome_table$id), ]
      predicted <- scored_rows$prediction
      if (type == 'continuous') return(mean((obs$y - predicted)^2))
      if (type == 'binary') {
        n1 <- sum(obs$y == 1)
        n0 <- sum(obs$y == 0)
        return((sum(rank(predicted)[obs$y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0))
      }
      survival::concordance(survival::Surv(obs$time, obs$status) ~ predicted)$concordance
    }
    computed <- saved$cv$pooled
    for (round in 1:10) {
      round_rows <- predictions[predictions$round == round, ]
      stopifnot(nrow(round_rows) == 300L, setequal(round_rows$id, outcome_table$id),
        setequal(round_rows$fold, 1:5))
      for (group in colnames(computed))
        computed[round, group] <- calculate(
          if (group == 'all') round_rows else round_rows[round_rows$subgroup == group, ])
    }
    stopifnot(isTRUE(all.equal(computed, saved$cv$pooled, tolerance = 1e-12)))
    prediction_table <- result$prediction_table
    index <- which(prediction_table$outcome == type & prediction_table$method == toupper(method))
    stopifnot(length(index) == 1L, isTRUE(all.equal(unname(colMeans(computed)),
     unname(as.numeric(prediction_table[index, colnames(computed)])), tolerance = 1e-12)))
    acceptance_rows[[length(acceptance_rows) + 1L]] <- data.frame(outcome = type, method = method,
     rounds = 10L, folds = 5L, predictions = nrow(predictions),
     max_score_difference = max(abs(computed - saved$cv$pooled)))
  }
  kirc <- readRDS(file.path(root, 'kirc_reduced_fit.rds'))
  stopifnot(kirc$control$mcmc$draws == 4000L, kirc$control$mcmc$burnin == 1000L,
   kirc$control$sampler_method == 'paper', result$simulated_fit$control$mcmc$draws == 8000L,
   result$simulated_fit$control$mcmc$burnin == 2000L)
  comparison <- readRDS(file.path(root, 'covariate-comparison/comparison.rds'))
  stopifnot(!comparison$settings$quick, comparison$settings$draws == 2000, comparison$settings$burnin == 500,
   comparison$settings$k == 3, comparison$settings$rounds == 2)
  outer_predictions <- comparison$outer_predictions
  stopifnot(nrow(outer_predictions) == 300L, !anyDuplicated(outer_predictions$id),
   setequal(outer_predictions$id, simIMR$outcome.continuous$id),
   all(is.finite(outer_predictions$prediction)),
   isTRUE(all.equal(mean((outer_predictions$observed - outer_predictions$prediction)^2),
    comparison$nested_summary$pooled_outer_mse, tolerance = 1e-12)))
  for (fold in 1:3) {
    test <- outer_predictions[outer_predictions$outer_fold == fold, ]
    inner <- comparison$nested_inner_folds[comparison$nested_inner_folds$outer_fold == fold, ]
    stopifnot(nrow(inner) == 200L, !anyDuplicated(inner$id),
     setequal(inner$id, setdiff(outer_predictions$id, test$id)), setequal(inner$fold, 1:3),
     !any(inner$id %in% test$id))
    fold_choice <- comparison$selected_by_fold[comparison$selected_by_fold$outer_fold == fold, ]
    chosen <- c('age_only', 'age_sex_stage')[which.min(c(fold_choice$inner_mse_age_only,
      fold_choice$inner_mse_age_sex_stage))]
    stopifnot(fold_choice$selected == chosen, fold_choice$n_train == nrow(inner),
     fold_choice$n_test == nrow(test), isTRUE(all.equal(mean((test$observed - test$prediction)^2),
      fold_choice$outer_mse, tolerance = 1e-12)))
  }
  for (candidate in c('age_only', 'age_sex_stage')) {
    paired_row <- comparison$paired_summary[comparison$paired_summary$candidate == candidate, ]
    round_mse <- comparison$paired_rounds[[candidate]]
    stopifnot(isTRUE(all.equal(mean(round_mse), paired_row$mean_mse, tolerance = 1e-12)),
     isTRUE(all.equal(sd(round_mse), paired_row$sd_across_rounds, tolerance = 1e-12)))
  }
  acceptance <- do.call(rbind, acceptance_rows)
  write.csv(acceptance, file.path(root, 'manuscript-acceptance.csv'), row.names = FALSE)
  cat('PASS: six full refit configurations, independent scores, nested holdout separation and summary aggregation:', root, '\n')
  invisible(acceptance)
}
for (root in args) validate_run(root)
if (length(args) == 2L) {
  a <- normalizePath(args[1])
  b <- normalizePath(args[2])
  csv <- setdiff(list.files(a, pattern = '[.]csv$', recursive = TRUE), 'manuscript-acceptance.csv')
  rds <- c(paste0('validation-', rep(c('binary', 'continuous', 'right.censored'), each = 2), '-', c('imr', 'bms'), '.rds'),
   'kirc_reduced_fit.rds', 'replication_results.rds', 'covariate-comparison/comparison.rds')
  comparison_rows <- lapply(c(csv, rds), function(name) {
    read <- if (endsWith(name, '.csv')) function(path) read.csv(path, check.names = FALSE) else readRDS
    first <- read(file.path(a, name))
    second <- read(file.path(b, name))
    same <- isTRUE(all.equal(first, second, tolerance = 1e-12))
    data.frame(file = name, identical = identical(first, second), equal_1e12 = same)
  })
  comparison_report <- do.call(rbind, comparison_rows)
  write.csv(comparison_report, file.path(dirname(b), 'AB-numeric-comparison.csv'), row.names = FALSE)
  stopifnot(all(comparison_report$equal_1e12))
  cat('PASS:', nrow(comparison_report), 'CSV/RDS artifacts agree across A/B\n')
}
