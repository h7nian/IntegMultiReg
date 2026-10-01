# Refit selection baselines from saved data and exact inner folds. This checks
# reproducibility, not the inferential validity of a warned Cox fit.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[1])
s <- readRDS(file.path(root, "settings.rds"))
stopifnot(s$experiment %in% c("simulation", "correlated"),
  identical(tools::md5sum(names(s$source_hashes)), s$source_hashes),
  identical(tools::md5sum(names(s$data_hash)), s$data_hash))
# Deliberately recomputes selection_auc() from the helpers by a different route:
# the pairwise probability that a positive exceeds a negative, rather than the
# rank form. Keeping the second derivation is what makes the check independent.
independent_pairwise_auc <- function(score, truth) {
  stopifnot(length(score) == length(truth), all(is.finite(score)))
  positive <- score[truth == 1]
  negative <- score[truth == 0]
  if (!length(positive) || !length(negative)) return(NA_real_)
  difference <- outer(positive, negative, "-")
  mean((difference > 0) + .5 * (difference == 0))
}
compare_rows <- function(actual, expected) {
  platform_subgroup_key <- function(x) paste(x$platform, x$subgroup, sep = ":")
  stopifnot(!anyDuplicated(platform_subgroup_key(actual)),
            !anyDuplicated(platform_subgroup_key(expected)),
    setequal(platform_subgroup_key(actual), platform_subgroup_key(expected)),
    isTRUE(all.equal(actual$auc[match(platform_subgroup_key(expected),
                                      platform_subgroup_key(actual))],
      expected$auc, tolerance = 1e-12, check.attributes = FALSE)))
}
reports <- list()
warning_records <- list()
configurations <- seq_len(if (s$experiment == "simulation") 3L else 6L)
replicates <- seq_len(s$replicates)
if (file.exists(file.path(root, "task.rds"))) {
  selection <- readRDS(file.path(root, "task.rds"))
  stopifnot(identical(selection$kind, s$experiment),
    all(selection$configurations %in% configurations), all(selection$replicates %in% replicates))
  configurations <- selection$configurations
  replicates <- selection$replicates
}
for (configuration in configurations) {
  for (replicate in replicates) {
    job <- file.path(root, sprintf("configuration-%02d-replicate-%03d",
                                   configuration, replicate))
    data <- readRDS(file.path(job, "data.rds"))$data
    l1 <- readRDS(file.path(job, "l1-cph.rds"))
    diagnostics <- Filter(function(x) x$stage == "selection", l1$diagnostics)
    stopifnot(length(diagnostics) > 0L,
      !anyDuplicated(vapply(diagnostics, `[[`, "", "subgroup")))
    warning_start <- length(warning_records)
    capture <- function(expr, method, platform = "", subgroup = "",
                        feature = "") {
    withCallingHandlers(expr, warning = function(w) {
      warning_records[[length(warning_records) + 1L]] <<- data.frame(
        configuration, replicate, method, platform, subgroup, feature,
        warning = conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    }
    rows <- list()
    for (d in diagnostics) {
    ids <- d$inner_folds$id
    stopifnot(!anyDuplicated(ids), all(ids %in% data$outcome$id))
    outcome <- data$outcome[match(ids, data$outcome$id), ]
    present <- vapply(data$platforms, function(p) all(ids %in% p$id), TRUE)
    design <- do.call(cbind, lapply(data$platforms[present], function(p)
      as.matrix(p[match(ids, p$id), -1])))
    fit <- capture(glmnet::cv.glmnet(design,
      survival::Surv(outcome$time, outcome$status), family = "cox", alpha = 1,
      standardize = TRUE, foldid = d$inner_folds$fold), "l1-cph", subgroup = d$subgroup)
    stopifnot(isTRUE(all.equal(fit$lambda.min, d$lambda_min, tolerance = 1e-12)),
      isTRUE(all.equal(fit$lambda, d$path$lambda, tolerance = 1e-12)),
      isTRUE(all.equal(fit$cvm, d$path$cv_mean, tolerance = 1e-12)),
      isTRUE(all.equal(fit$cvsd, d$path$cv_sd, tolerance = 1e-12)),
      fit$glmnet.fit$jerr == d$full_fit_error)
    coefficients <- abs(as.numeric(coef(fit, s = "lambda.min")))
    offset <- 0L
    for (p in which(present)) {
      n <- ncol(data$platforms[[p]]) - 1L
      rows[[length(rows) + 1L]] <- data.frame(platform = names(data$platforms)[p],
        subgroup = d$subgroup, auc = independent_pairwise_auc(coefficients[offset + seq_len(n)],
          data$truth[[p]][d$subgroup, ]))
      offset <- offset + n
    }
    stopifnot(offset == length(coefficients))
    }
    compare_rows(do.call(rbind, rows), l1$selection)
    rows <- list()
    for (p in seq_along(data$platforms)) {
    platform <- data$platforms[[p]]
    outcome <- data$outcome[match(platform$id, data$outcome$id), ]
    values <- vapply(seq.int(2L, ncol(platform)), function(j) {
      marker <- platform[[j]]
      fit <- capture(survival::coxph(survival::Surv(outcome$time, outcome$status) ~ marker),
        "uni-cph", names(data$platforms)[p], feature = names(platform)[j])
      summary(fit)$coefficients[1, "Pr(>|z|)"]
    }, numeric(1))
    stopifnot(all(is.finite(values)))
    score <- 1 - stats::p.adjust(values, method = "hommel")
    groups <- rownames(data$truth[[p]])
    rows[[p]] <- data.frame(platform = names(data$platforms)[p],
                            subgroup = groups,
      auc = vapply(groups, function(g) independent_pairwise_auc(score, data$truth[[p]][g, ]), numeric(1)))
    }
    compare_rows(do.call(rbind, rows), read.csv(file.path(job, "uni-cph-selection.csv"),
      colClasses = c(subgroup = "character")))
    reports[[length(reports) + 1L]] <- data.frame(configuration, replicate,
      l1_groups = length(diagnostics), uni_features = sum(vapply(data$platforms, ncol, 1L) - 1L),
      warnings = length(warning_records) - warning_start)
    cat("PASS baseline selection:", configuration, replicate, "\n")
  }
}
report <- do.call(rbind, reports)
warnings <- if (length(warning_records)) do.call(rbind, warning_records) else data.frame(
 configuration = integer(), replicate = integer(), method = character(), platform = character(),
 subgroup = character(), feature = character(), warning = character())
write.csv(report, file.path(root, "baseline-selection-acceptance.csv"),
          row.names = FALSE)
write.csv(warnings, file.path(root, "baseline-selection-warnings.csv"),
          row.names = FALSE)
saveRDS(list(report = report, warnings = warnings, settings = s,
             verified_at = Sys.time(),
  validator_hash = tools::md5sum(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
  session = sessionInfo(), scope = "L1 selection paths replayed with saved inner folds; L1 and Hommel-adjusted univariate selection AUC independently recomputed. Warnings do not certify scientific stability."),
  file.path(root, "baseline-selection-acceptance.rds"))
cat("PASS all executed baseline selection jobs; warning count:", nrow(warnings),
    "\n")
