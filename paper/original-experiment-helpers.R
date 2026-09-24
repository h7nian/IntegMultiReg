# Research orchestration; fitting and validation use public package functions.
parse_experiment_args <- function(args) {
  value_flags <- c("--experiment", "--reference", "--out-dir", "--draws",
    "--burnin", "--k", "--rounds", "--replicates", "--workers", "--seed", "--data", "--marker-design", "--chain", "--configuration", "--replicate",
    "--ridge", "--cv-method", "--sampler-method", "--prior-indexing",
    "--model-set", "--df-method", "--score-method", "--fold-rng")
  options <- list()
  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (!flag %in% c(value_flags, "--quick", "--audit")) stop("Unknown argument: ", flag)
    if (flag %in% names(options)) stop("Duplicate argument: ", flag)
    if (flag %in% c("--quick", "--audit")) {
      options[[flag]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args) || startsWith(args[[i + 1L]], "--") ||
          !nzchar(args[[i + 1L]])) stop("Missing value for ", flag)
      options[[flag]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  options
}

# MSI task identity is separate from global scientific settings.
record_task_selection <- function(out, selection) {
  path <- file.path(out,'task.rds')
  if(file.exists(path) && !identical(readRDS(path),selection))
    stop('Task selector changed; choose a new output directory')
  saveRDS(selection,path)
}

# Lossless research-only checkpoint container. Public fitted objects stay in
# their ordinary format; unique matrices plus per-draw indices survive RDS.
pack_experiment_fit <- function(fit) {
  if (!inherits(fit, "imr")) return(fit)
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The research checkpoint writer requires the digest package")
  history <- fit$posterior$selection_draws
  stopifnot(length(history) > 0L)
  platforms <- length(history[[1L]])
  pools <- indices <- vector("list", platforms)
  for (p in seq_len(platforms)) {
    states <- lapply(history, `[[`, p)
    # Hash short serialized-matrix keys instead of matching lists, which coerces
    # large matrices to long character strings. Verify every reconstruction so
    # even a hash collision cannot silently change a retained state.
    keys <- vapply(states, digest::digest, character(1), algo = "xxhash64")
    unique <- !duplicated(keys)
    pool <- states[unique]
    index <- match(keys, keys[unique])
    if (anyNA(index) || !all(vapply(seq_along(states), function(i)
        identical(states[[i]], pool[[index[i]]]), logical(1))))
      stop("Selection history could not be packed losslessly")
    pools[[p]] <- pool
    indices[[p]] <- index
  }
  layout <- lapply(history, function(draw) {
    draw[] <- rep(list(NULL), length(draw)); draw
  })
  attributes(layout) <- attributes(history)
  fit$posterior["selection_draws"] <- list(NULL)
  structure(list(fit = fit, layout = layout, pools = pools, indices = indices),
    class = "imr_experiment_fit_checkpoint_v1")
}

unpack_experiment_fit <- function(saved) {
  if (!inherits(saved, "imr_experiment_fit_checkpoint_v1")) return(saved)
  history <- saved$layout
  for (p in seq_along(saved$pools)) {
    index <- saved$indices[[p]]
    stopifnot(length(index) == length(history), !anyNA(index),
      all(index >= 1L & index <= length(saved$pools[[p]])))
    for (s in seq_along(history)) history[[s]][[p]] <- saved$pools[[p]][[index[s]]]
  }
  # Preserve posterior field order, too: replacement into a NULL-valued slot.
  saved$fit$posterior$selection_draws <- history
  saved$fit
}

experiment_fit_checkpoint <- function(path, code) {
  if (file.exists(path)) return(unpack_experiment_fit(readRDS(path)))
  fit <- eval(substitute(code), parent.frame())
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(pack_experiment_fit(fit), temporary)
  if (!file.rename(temporary, path)) stop("Cannot commit fit checkpoint: ", path)
  fit
}

experiment_result_folds <- function(path) {
  sidecar <- paste0(path, ".folds.rds")
  if (file.exists(sidecar)) return(readRDS(sidecar))
  # Compatibility for results written before the small metadata sidecar.
  read_experiment_result(path, include_fit = FALSE)$cv$control$folds
}

experiment_job_complete <- function(path) {
  status_path <- paste0(path, ".status.rds")
  file.exists(path) && file.exists(status_path) &&
    identical(readRDS(status_path)$status, "completed")
}

record_experiment_job <- function(path, code) {
  expression <- substitute(code)
  caller <- parent.frame()
  status_path <- paste0(path, ".status.rds")
  attempt <- 1L
  while (file.exists(paste0(path, ".attempt-", attempt, ".rds"))) attempt <- attempt + 1L
  attempt_path <- paste0(path, ".attempt-", attempt, ".rds")
  record <- list(status = "running", started = Sys.time(), ended = NULL,
    pid = Sys.getpid(), attempt = attempt, warnings = list(), error = NULL,
    memory_scope = "R heap in this process; excludes native allocations and workers")
  record$heap_start <- gc(reset = TRUE)
  started <- proc.time()
  saveRDS(record, attempt_path)
  saveRDS(record, status_path)
  on.exit({
    if (record$status == "running") record$status <- "interrupted"
    record$ended <- Sys.time()
    record$elapsed <- proc.time() - started
    record$heap_end <- gc()
    saveRDS(record, attempt_path)
    saveRDS(record, status_path)
  }, add = TRUE)
  result <- tryCatch(withCallingHandlers(eval(expression, caller), warning = function(w) {
    record$warnings[[length(record$warnings) + 1L]] <<- list(
      time = Sys.time(), message = conditionMessage(w), call = conditionCall(w))
  }), error = function(e) {
    record$status <<- "failed"
    record$error <<- list(message = conditionMessage(e), call = conditionCall(e))
    stop(e)
  })
  record$status <- "completed"
  result
}

reference_arguments <- function(reference) {
  switch(reference,
    paper = list(fit = list(sampler_method = "paper", prior_indexing = "standard"),
      cv = list(cv_method = "importance", ridge = 0, model_set = "draws",
                df_method = "fractional", score_method = "standard")),
    code2017 = list(fit = list(sampler_method = "legacy", prior_indexing = "code2017",
      laplace_max_iter = c(initial = 25L, selection = 40L, latent = 25L, prediction = 25L)),
      cv = list(cv_method = "legacy", ridge = 0, model_set = "ranked_unique",
                df_method = "legacy_integer", score_method = "legacy", fold_rng = "continue")),
    package = list(fit = list(), cv = list(cv_method = "legacy")),
    stop("reference must be paper, code2017 or package"))
}

# Reference names provide ordinary argument defaults. Explicit overrides are
# stored in the manifest and then validated by the public package entrypoints.
resolve_experiment_arguments <- function(reference, options) {
  arguments <- reference_arguments(reference)
  if (identical(options[["--cv-method"]], "refit"))
    arguments$cv <- list(cv_method = "refit")
  fit_names <- c("sampler_method", "prior_indexing")
  cv_names <- c("cv_method", "ridge", "model_set", "df_method", "score_method", "fold_rng")
  for (name in c(fit_names, cv_names)) {
    flag <- paste0("--", gsub("_", "-", name, fixed = TRUE))
    value <- options[[flag]]
    if (is.null(value)) next
    if (name == "ridge") {
      value <- suppressWarnings(as.numeric(value))
      if (length(value) != 1L || !is.finite(value) || value < 0)
        stop("--ridge must be a finite nonnegative number")
    }
    section <- if (name %in% fit_names) "fit" else "cv"
    arguments[[section]][[name]] <- value
  }
  if (identical(arguments$cv$cv_method, "refit") &&
      any(c("ridge", "model_set", "df_method", "fold_rng") %in% names(arguments$cv)))
    stop("Refit CV cannot use post-fit-only overrides")
  arguments
}

selection_auc <- function(score, truth) {
  positive <- truth == 1; n1 <- sum(positive); n0 <- sum(!positive)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

survival_score <- function(time, event, prediction, method = "standard") {
  if (method == "standard") return(as.numeric(survival::concordance(
    survival::Surv(time, event) ~ prediction)$concordance))
  # Independent pairwise transcription of the archived historical metric.
  pair <- which(upper.tri(matrix(FALSE, length(time), length(time))), arr.ind = TRUE)
  i <- pair[, 1]; j <- pair[, 2]
  denominator <- (time[j] > time[i]) * event[i] + (time[j] < time[i]) * event[j] +
    (time[j] == time[i]) * (event[i] != event[j])
  numerator <- (prediction[j] > prediction[i]) * (time[j] > time[i]) * event[i] +
    (prediction[j] < prediction[i]) * (time[j] < time[i]) * event[j] +
    .5 * ((prediction[j] == prediction[i]) | (time[j] == time[i])) * (event[i] != event[j])
  if (sum(denominator)) sum(numerator) / sum(denominator) else NA_real_
}

fold_scores <- function(predictions, outcome, score_method, weighted_overall = FALSE) {
  rows <- split(predictions, interaction(predictions$round, predictions$fold, drop = TRUE))
  do.call(rbind, lapply(rows, function(x) {
    groups <- split(x, x$subgroup)
    values <- vapply(groups, function(g) {
      y <- outcome[match(g$id, outcome$id), ]
      survival_score(y$time, y$status, g$prediction, score_method)
    }, 0)
    y <- outcome[match(x$id, outcome$id), ]
    full <- if (weighted_overall) weighted.mean(values, vapply(groups, nrow, 0L)) else
      survival_score(y$time, y$status, x$prediction, score_method)
    data.frame(round = x$round[1], fold = x$fold[1], subgroup = c(names(values), "all"),
               value = c(values, full))
  }))
}

summarize_scores <- function(scores) {
  do.call(rbind, lapply(split(scores, scores$subgroup), function(x) {
    finite <- is.finite(x$value); n <- sum(finite)
    data.frame(subgroup = x$subgroup[1], mean = if (n) mean(x$value[finite]) else NA_real_,
      sd = if (n > 1) sd(x$value[finite]) else NA_real_,
      se = if (n > 1) sd(x$value[finite]) / sqrt(n) else NA_real_,
      n_valid = n, n_expected = nrow(x))
  }))
}

biomarker_auc <- function(fit, truth) {
  values <- coef(fit)
  do.call(rbind, lapply(seq_along(values), function(p) {
    groups <- fit$model$subgroup_names[fit$model$platform_subgroups[[p]]]
    data.frame(platform = fit$model$platform_names[p], subgroup = groups,
      auc = vapply(seq_along(groups), function(g)
        selection_auc(values[[p]][g, ], truth[[p]][groups[g], ]), 0))
  }))
}

# Record the exact inner partitions and regularization path used by glmnet.
fit_cox_path <- function(design, outcome, seed, inner_k) {
  set.seed(seed)
  foldid <- sample(rep(seq_len(inner_k), length.out = nrow(design)))
  warnings <- character()
  fit <- withCallingHandlers(glmnet::cv.glmnet(design,
    survival::Surv(outcome$time, outcome$status), family = "cox",
    foldid = foldid, alpha = 1, standardize = TRUE), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
    })
  list(fit = fit, diagnostics = list(warnings = warnings,
    inner_folds = data.frame(id = outcome$id, fold = foldid),
    lambda_min = fit$lambda.min, full_fit_error = fit$glmnet.fit$jerr,
    path = data.frame(lambda = fit$lambda, cv_mean = fit$cvm,
      cv_sd = fit$cvsd, nonzero = fit$nzero)))
}

benchmark_cox <- function(data, folds, lasso, seed, inner_k = 10L) {
  cohort <- IntegMultiReg::imr_data(data$platforms, data$outcome,
    covariates = data$covariates, outcome_type = "right.censored")$availability
  record <- merge(folds, cohort[c("id", "subgroup")], by = "id", sort = FALSE)
  record$prediction <- NA_real_
  selection <- list()
  diagnostics <- list()
  for (group in unique(record$subgroup)) {
    ids <- cohort$id[cohort$subgroup == group]
    outcome <- data$outcome[match(ids, data$outcome$id), ]
    if (lasso) {
      present <- vapply(data$platforms, function(x) all(ids %in% x$id), TRUE)
      blocks <- lapply(data$platforms[present], function(x) as.matrix(x[match(ids, x$id), -1]))
      design <- do.call(cbind, blocks)
      path_fit <- fit_cox_path(design, outcome, seed, inner_k)
      full <- path_fit$fit
      diagnostics[[length(diagnostics) + 1L]] <- c(list(stage = "selection",
        subgroup = group, n_train = nrow(outcome), n_events = sum(outcome$status)),
        path_fit$diagnostics)
      coefficients <- abs(as.numeric(coef(full, s = "lambda.min")))
      if (!is.null(data$truth)) {
        offset <- 0L
        for (p in which(present)) {
          n <- ncol(data$platforms[[p]]) - 1L
          selection[[length(selection) + 1L]] <- data.frame(platform = names(data$platforms)[p],
            subgroup = group, auc = selection_auc(coefficients[offset + seq_len(n)], data$truth[[p]][group, ]))
          offset <- offset + n
        }
      }
    } else design <- as.matrix(data$covariates[match(ids, data$covariates$id), -1])
    for (round in sort(unique(record$round))) for (fold in sort(unique(record$fold))) {
      target <- which(record$subgroup == group & record$round == round & record$fold == fold)
      test <- match(record$id[target], ids); train <- setdiff(seq_along(ids), test)
      y <- outcome[train, ]
      fit_warnings <- character()
      withCallingHandlers({
      if (lasso) {
        path_fit <- fit_cox_path(design[train, , drop = FALSE], y,
          seed + 100L * round + fold, inner_k)
        fitted <- path_fit$fit
        score <- as.numeric(predict(fitted, design[test, , drop = FALSE], s = "lambda.min", type = "link"))
      } else {
        training <- data.frame(time = y$time, status = y$status, design[train, , drop = FALSE])
        fitted <- survival::coxph(survival::Surv(time, status) ~ ., training, ties = "breslow")
        score <- as.numeric(predict(fitted, data.frame(design[test, , drop = FALSE]), type = "lp"))
      }
      }, warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
      })
      diagnostics[[length(diagnostics) + 1L]] <- list(
        stage = "prediction", subgroup = group, round = round, fold = fold, n_train = length(train),
        n_events = sum(y$status), variables = colnames(design),
        warnings = fit_warnings,
        coefficients = if (lasso) NULL else stats::coef(fitted),
        regularization = if (lasso) path_fit$diagnostics else NULL)
      record$prediction[target] <- -score
    }
  }
  if (any(!is.finite(record$prediction))) stop("Cox benchmark produced non-finite predictions")
  list(predictions = record, selection = if (length(selection)) do.call(rbind, selection) else NULL,
    diagnostics = diagnostics)
}

univariate_cox_auc <- function(data) {
  do.call(rbind, lapply(seq_along(data$platforms), function(p) {
    platform <- data$platforms[[p]]
    y <- data$outcome[match(platform$id, data$outcome$id), ]
    pvalues <- vapply(platform[-1], function(marker) {
      fit <- survival::coxph(survival::Surv(y$time, y$status) ~ marker)
      summary(fit)$coefficients[1, "Pr(>|z|)"]
    }, 0)
    if (any(!is.finite(pvalues))) stop("Univariate Cox produced non-finite p-values")
    score <- 1 - p.adjust(pvalues, "hommel")
    data.frame(platform = names(data$platforms)[p], subgroup = rownames(data$truth[[p]]),
      auc = apply(data$truth[[p]], 1, function(truth) selection_auc(score, truth)))
  }))
}

# Versioned result container reuses lossless fit packing. Legacy ordinary RDS
# results remain readable; numerical values and selection-state order are kept.
write_experiment_result <- function(result, path) {
  saved <- result
  if (inherits(saved$fit, 'imr')) saved$fit <- pack_experiment_fit(saved$fit)
  container <- structure(list(result = saved), class = 'imr_experiment_result_v1')
  temporary <- paste0(path, '.tmp-', Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(container, temporary, compress = "xz")
  if (!file.rename(temporary, path)) stop('Cannot commit result: ', path)
  invisible(path)
}

read_experiment_result <- function(path, include_fit = TRUE) {
  saved <- readRDS(path)
  if (inherits(saved, 'imr_experiment_result_v1')) saved <- saved$result
  if (include_fit) saved$fit <- unpack_experiment_fit(saved$fit)
  else saved$fit <- NULL
  saved
}
