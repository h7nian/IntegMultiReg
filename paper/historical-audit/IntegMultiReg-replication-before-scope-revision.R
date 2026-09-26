###############################################################################
## Replication code for:
##
## Zhang S, Wang J, Chekouo T
## "IntegMultiReg: An R Package for Integrative Bayesian Multi-Regression of
##  Multi-Platform Biomarkers"
##
## Corresponding R package:
##   https://CRAN.R-project.org/package=IntegMultiReg
##
## This single script reproduces the manuscript's empirical figures and tables.
## The default run reproduces the simulated-data results and the data-summary
## tables.  The historical TCGA-KIRC Table 1 reanalysis is included below and is
## enabled with --full-kirc because its published configuration can take many
## hours.  Use --smoke-kirc to exercise that complete code path quickly.
##
## Exact manuscript run (figures and simulation tables):
##   Rscript paper/IntegMultiReg-replication.R
##
## Quick code-path check (does not reproduce manuscript numbers):
##   Rscript paper/IntegMultiReg-replication.R --quick
##
## Add the full historical KIRC reanalysis:
##   Rscript paper/IntegMultiReg-replication.R --full-kirc
##
## Requirements:
##   * R >= 4.0.0
##   * IntegMultiReg 0.1.4, installed with its GSL system requirement
##   * survival (only for --full-kirc or --smoke-kirc)
##
## Reproducibility notes:
##   * Version 0.1.4 uses independent training-fold MCMC refits and corrected
##     binary probabilities/concordance; old 0.1.3 CV tables are not references.
##   * Version 0.1.4 simulation/reduced KIRC results use log time.
##   * Optional full/smoke KIRC paths retain explicit identity scale solely
##     for the legacy appendix; they do not validate the corrected log model.
##   * All stochastic analyses use explicit seeds.
##   * The normal run writes to paper/replication-output and never overwrites
##     the manuscript's checked-in figures.
##   * --full-kirc expects paper/data/kirc_table1_full.rda, created from the
##     Biometrics article's supplementary data as documented in data/README.md.
###############################################################################

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) any(args == flag)

arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) default else args[hit + 1L]
}

## Locate the project whether this script is started from the repository root,
## from paper/, or by an IDE whose working directory is elsewhere.
script_flag <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_flag)) sub("^--file=", "", script_flag[1L]) else ""
script_dir <- if (nzchar(script_file)) {
  dirname(normalizePath(script_file, mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

root_candidates <- unique(c(dirname(script_dir), script_dir,
                            normalizePath(getwd(), mustWork = TRUE)))
root_hit <- root_candidates[vapply(root_candidates, function(x) {
  dir.exists(file.path(x, "IntegMultiReg")) && dir.exists(file.path(x, "paper"))
}, logical(1L))]
inside_repository <- length(root_hit) > 0L
project_root <- if (inside_repository) root_hit[1L] else script_dir
materials_root <- if (inside_repository) {
  file.path(project_root, "paper")
} else {
  script_dir
}

quick <- has_flag("--quick")
run_full_kirc <- has_flag("--full-kirc")
run_smoke_kirc <- has_flag("--smoke-kirc")
run_reduced_kirc_fit <- !has_flag("--skip-reduced-kirc-fit")
update_reference <- has_flag("--update-reference")
reference_dir <- arg_value("--reference-dir", file.path(materials_root, "expected-results-0.1.4"))
if (update_reference && quick) stop("Reference tables require a full run, not --quick.")

default_output <- file.path(materials_root, "replication-output")
if (quick) default_output <- paste0(default_output, "-quick")
out_dir <- normalizePath(arg_value("--out-dir", default_output),
                         mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("IntegMultiReg", quietly = TRUE)) {
  stop("Install IntegMultiReg before running this script. For a local source ",
       "checkout, run: R CMD INSTALL IntegMultiReg")
}
if (utils::packageVersion("IntegMultiReg") != "0.1.4") {
  stop("This replication script requires exactly IntegMultiReg 0.1.4; installed: ",
       as.character(utils::packageVersion("IntegMultiReg")))
}

suppressPackageStartupMessages(library(IntegMultiReg))
data("simIMR", package = "IntegMultiReg")
data("kircIMR", package = "IntegMultiReg")

cat("Replicating the IntegMultiReg manuscript\n")
cat("Package version:", as.character(utils::packageVersion("IntegMultiReg")),
    "\nOutput directory:", out_dir, "\n")
if (quick) {
  cat("QUICK MODE: outputs validate the code path but do not match the paper.\n")
}

write_table <- function(x, name) {
  path <- file.path(out_dir, paste0(name, ".csv"))
  utils::write.csv(x, path, row.names = FALSE)
  cat("Saved", path, "\n")
  invisible(path)
}

all_ids <- function(dat) unique(unlist(lapply(dat$platforms, `[[`, "id")))

complete_ids <- function(dat) {
  Reduce(intersect, lapply(dat$platforms, `[[`, "id"))
}

subgroup_counts <- function(dat) {
  ids <- all_ids(dat)
  present <- vapply(dat$platforms, function(x) ids %in% x$id,
                    logical(length(ids)))
  bits <- apply(present, 1L, function(x) {
    paste(as.integer(rev(x)), collapse = "")
  })
  table(bits)
}

number_of_features <- function(dat) {
  sum(vapply(dat$platforms, function(x) ncol(x) - 1L, integer(1L)))
}

sim_subgroups <- subgroup_counts(simIMR)
sim_modelled_subgroups <- sim_subgroups[sim_subgroups >= 30L]

## --------------------------------------------------------------------------
## Manuscript Table: computational scale of the replication workflows
## --------------------------------------------------------------------------

full_data_path <- arg_value(
  "--kirc-data",
  file.path(materials_root, "data", "kirc_table1_full.rda")
)

computational_scale <- data.frame(
  workflow = c("simIMR figures and prediction table",
               "kircIMR reduced example",
               "Biometrics supplement KIRC"),
  subjects = c(length(all_ids(simIMR)),
               sum(kircIMR$model_subgroup_sizes),
               448L),
  features = c(number_of_features(simIMR),
               number_of_features(kircIMR),
               1598L),
  subgroups = c(length(sim_modelled_subgroups),
                length(kircIMR$model_subgroup_sizes),
                4L),
  configuration = c(
    paste("Figures: 8000 retained + 2000 burn-in; prediction:",
          "4000 retained + 1000 burn-in, 10 x 5-fold CV"),
    "4000 retained + 1000 burn-in for the illustrative fit",
    "400000 retained + 50000 burn-in, 10 x 10-fold CV"
  ),
  check.names = FALSE
)
cat("\n=== Computational scale ===\n")
print(computational_scale, row.names = FALSE)
write_table(computational_scale, "table_computational_scale")

## --------------------------------------------------------------------------
## Manuscript Table: complete-case sample loss
## --------------------------------------------------------------------------

complete_case <- data.frame(
  data_set = c("simIMR", "kircIMR reduced example"),
  any_platform = c(length(all_ids(simIMR)), length(all_ids(kircIMR))),
  complete_cases = c(length(complete_ids(simIMR)),
                     length(complete_ids(kircIMR))),
  imr_modelled = c(sum(sim_modelled_subgroups),
                   sum(kircIMR$model_subgroup_sizes)),
  check.names = FALSE
)
complete_case$complete_case_percent <- round(
  100 * complete_case$complete_cases / complete_case$any_platform, 1L
)
cat("\n=== Complete-case sample loss ===\n")
print(complete_case, row.names = FALSE)
write_table(complete_case, "table_complete_case")

## --------------------------------------------------------------------------
## Reduced TCGA-KIRC data summary printed in the manuscript
## --------------------------------------------------------------------------

kirc_dimensions <- as.data.frame(t(sapply(kircIMR$platforms, dim)))
names(kirc_dimensions) <- c("subjects", "columns_including_id")
kirc_dimensions$platform <- rownames(kirc_dimensions)
kirc_dimensions$variables <- kirc_dimensions$columns_including_id - 1L
kirc_dimensions <- kirc_dimensions[, c("platform", "subjects", "variables")]

kirc_subgroups <- data.frame(
  subgroup = names(kircIMR$model_subgroup_sizes),
  subjects = as.integer(kircIMR$model_subgroup_sizes),
  row.names = NULL
)
cat("\n=== Reduced TCGA-KIRC dimensions ===\n")
print(kirc_dimensions, row.names = FALSE)
cat("\n=== Reduced TCGA-KIRC modelled subgroups ===\n")
print(kirc_subgroups, row.names = FALSE)
write_table(kirc_dimensions, "kirc_reduced_dimensions")
write_table(kirc_subgroups, "kirc_reduced_subgroups")

## Execute the exact data-inspection commands displayed in the manuscript.
cat("\n=== Manuscript code: kircIMR dimensions ===\n")
dims <- sapply(kircIMR$platforms, dim)
rownames(dims) <- c("subjects", "variables")
print(dims)
cat("\n=== Manuscript code: kircIMR subgroup sizes ===\n")
print(kircIMR$model_subgroup_sizes)

## Execute the reduced real-data fit displayed in the manuscript. Quick mode
## uses a short chain; the normal replication uses the displayed MCMC settings.
if (run_reduced_kirc_fit) {
  cat("\n=== Fitting the reduced TCGA-KIRC illustration ===\n")
  reduced_mcmc <- if (quick) c(50L, 20L) else c(4000L, 1000L)
  fit_kirc <- imr(
    platform_data_list = kircIMR$platforms,
    outcome = kircIMR$outcome.survival,
    cov = kircIMR$covariates,
    type_outcome = "right.censored",
    nu = c(-4, -3, -4),
    sample_mcmc = reduced_mcmc,
    ssize = 30,
    seed = 1
  )
  kirc_summary <- capture.output(summary(fit_kirc, threshold = 0.5))
  cat(paste(kirc_summary, collapse = "\n"), "\n")
  writeLines(kirc_summary, file.path(out_dir, "kirc_reduced_summary.txt"))
  saveRDS(fit_kirc, file.path(out_dir, "kirc_reduced_fit.rds"))
}

## --------------------------------------------------------------------------
## Manuscript Figures and biomarker-recovery table: simulated survival data
## --------------------------------------------------------------------------

figure_mcmc <- if (quick) c(50L, 20L) else c(8000L, 2000L)
prediction_mcmc <- if (quick) c(50L, 20L) else c(4000L, 1000L)
cv_k <- if (quick) 2L else 5L
cv_rounds <- if (quick) 1L else 10L
nu <- c(-4, -3, -4)

cat("\n=== Fitting the simulated right-censored model ===\n")
analysis_data <- imr_data(
  platforms = simIMR$platforms,
  outcome = simIMR$outcome.survival,
  covariates = simIMR$covariates,
  type_outcome = "right.censored"
)
fit <- imr(
  platform_data_list = analysis_data,
  nu = nu,
  sample_mcmc = figure_mcmc,
  ssize = 30,
  seed = 1
)

## Execute and print every fitted-object command shown in the manuscript.
cat("\n=== Manuscript code: print(fit, threshold = 0.5) ===\n")
print(fit, threshold = 0.5)
cat("\n=== Manuscript code: print(fit, rank = TRUE, top = 3) ===\n")
print(fit, rank = TRUE, top = 3)
cat("\n=== Manuscript code: summary(fit, threshold = 0.5) ===\n")
print(summary(fit, threshold = 0.5))

cat("\n=== Manuscript code: predict(fit, ...) ===\n")
newX <- simIMR$platforms$genomic[1:20, ]
newP <- simIMR$platforms$proteomic[1:20, ]
pred <- predict(
  fit,
  newdata = list(newX, newP),
  covariates = simIMR$covariates
)
print(pred[["model:011"]][1:3, ], digits = 3)

## Figure: platform- and subgroup-specific marginal inclusion probabilities.
grDevices::pdf(file.path(out_dir, "selection.pdf"), width = 9, height = 3.7)
plot(fit, type = "selection")
grDevices::dev.off()

## Figure: ranked marginal inclusion probabilities for the leading biomarkers.
grDevices::pdf(file.path(out_dir, "top_features.pdf"), width = 6.5, height = 4)
plot_top_features(fit, top = 10)
grDevices::dev.off()

## Figure: log-posterior trace used for the diagnostic illustration.
grDevices::pdf(file.path(out_dir, "trace.pdf"), width = 6, height = 3.6)
plot(fit, type = "trace")
grDevices::dev.off()

truth <- simIMR$truth
mpip <- coef(fit)
recovery <- data.frame(
  platform = names(mpip),
  n_features = vapply(mpip, ncol, integer(1L)),
  n_true = vapply(names(truth), function(p) length(truth[[p]]), integer(1L)),
  max_mpip_true = round(mapply(function(m, idx) {
    max(apply(m[, idx, drop = FALSE], 2L, max))
  }, mpip, truth), 3L),
  max_mpip_null = round(mapply(function(m, idx) {
    null_index <- setdiff(seq_len(ncol(m)), idx)
    max(apply(m[, null_index, drop = FALSE], 2L, max))
  }, mpip, truth), 3L),
  row.names = NULL,
  check.names = FALSE
)
cat("\n=== Biomarker recovery ===\n")
print(recovery, row.names = FALSE)
write_table(recovery, "table_biomarker_recovery")

if (!quick) {
  expected_recovery <- c(0.055, 0.250, 0.065)
  if (max(abs(recovery$max_mpip_null - expected_recovery)) > 0.005) {
    stop("Recovery values differ from the manuscript by more than 0.005.")
  }
}

## --------------------------------------------------------------------------
## Manuscript Table: IMR versus BMS across all three outcome types
## --------------------------------------------------------------------------

outcomes <- list(
  binary = simIMR$outcome.binary,
  continuous = simIMR$outcome.continuous,
  right.censored = simIMR$outcome.survival
)
metric_name <- c(binary = "AUC", continuous = "MSE",
                 right.censored = "C-index")

compare_methods <- function(type) {
  result <- lapply(c("IMR", "BMS"), function(method) {
    cat("Starting training-fold CV:", type, method, "\n")
    model <- imr(
      platform_data_list = simIMR$platforms,
      outcome = outcomes[[type]],
      cov = simIMR$covariates,
      type_outcome = type,
      method = method,
      nu = nu,
      sample_mcmc = prediction_mcmc,
      ssize = 30,
      seed = 1
    )
    cv <- cv_imr(model, k = cv_k, rounds = cv_rounds)
    cat("Completed training-fold CV:", type, method, "\n")
    colMeans(cv$total_cindex)
  })
  names(result) <- c("IMR", "BMS")
  result
}

set.seed(1)
comparison <- lapply(names(outcomes), compare_methods)
names(comparison) <- names(outcomes)

prediction_table <- do.call(rbind, lapply(names(comparison), function(type) {
  value <- comparison[[type]]
  data.frame(
    outcome = type,
    measure = unname(metric_name[type]),
    method = c("IMR", "BMS"),
    rbind(round(value$IMR, 3L), round(value$BMS, 3L)),
    check.names = FALSE,
    row.names = NULL
  )
}))
cat("\n=== IMR versus BMS prediction ===\n")
print(prediction_table, row.names = FALSE)
write_table(prediction_table, "table_prediction_comparison")

if (!quick) {
  reference_path <- file.path(reference_dir, "table_prediction_comparison.csv")
  if (update_reference) {
    dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(prediction_table, reference_path, row.names = FALSE)
    cat("Wrote candidate reference; verify with a second run without --update-reference.\n")
  } else {
    if (!file.exists(reference_path)) {
      stop("No 0.1.4 reference exists. Run once with --update-reference, then verify independently.")
    }
    expected_prediction <- utils::read.csv(reference_path, check.names = FALSE)
    if (!identical(names(prediction_table), names(expected_prediction)) ||
        !identical(prediction_table[, 1:3], expected_prediction[, 1:3]) ||
        !isTRUE(all.equal(as.matrix(prediction_table[, -(1:3)]),
                          as.matrix(expected_prediction[, -(1:3)]), tolerance = 1e-12))) {
      stop("Prediction values differ from the 0.1.4 reference. Check source, settings and sessionInfo().")
    }
  }
}

## --------------------------------------------------------------------------
## Historical TCGA-KIRC Table 1 package-based reanalysis
## --------------------------------------------------------------------------
## This section reanalyses the supplementary KIRC data. It is not run
## by default because the manuscript configuration uses four long MCMC fits.

availability_bitstrings <- function(dat) {
  outcome <- if (!is.null(dat$outcome.survival)) dat$outcome.survival else dat$outcome
  ids <- outcome$id
  availability <- data.frame(id = ids, check.names = FALSE)
  for (name in names(dat$platforms)) {
    availability[[name]] <- ids %in% dat$platforms[[name]]$id
  }
  bits <- apply(availability[names(dat$platforms)], 1L, function(x) {
    paste(as.integer(rev(x)), collapse = "")
  })
  data.frame(id = ids, bitstring = bits, check.names = FALSE)
}

cindex_from_score <- function(time, status, score) {
  ok <- is.finite(time) & !is.na(status) & is.finite(score)
  if (sum(ok) < 3L || length(unique(status[ok])) < 2L) return(NA_real_)
  as.numeric(survival::concordance(
    survival::Surv(time[ok], status[ok]) ~ score[ok],
    reverse = TRUE
  )$concordance)
}

summarise_cv <- function(x) {
  mean <- colMeans(x, na.rm = TRUE)
  sd <- apply(x, 2L, stats::sd, na.rm = TRUE)
  if (nrow(x) == 1L) sd[] <- NA_real_
  list(mean = mean, sd = sd)
}

fold_ids <- function(n, k) sample(rep(seq_len(k), length.out = n))

run_cph_cv <- function(dat, k = 10L, rounds = 10L, seed = 1L) {
  set.seed(seed)
  outcome <- if (!is.null(dat$outcome.survival)) dat$outcome.survival else dat$outcome
  covariates <- dat$covariates
  covariate_names <- setdiff(names(covariates), "id")
  covariates <- covariates[match(outcome$id, covariates$id), , drop = FALSE]
  availability <- availability_bitstrings(dat)
  bits <- availability$bitstring[match(outcome$id, availability$id)]
  paper_bits <- c("111", "011", "101", "001")
  paper_groups <- c("E1", "E2", "E3", "E5")

  group_values <- matrix(NA_real_, nrow = rounds, ncol = length(paper_bits),
                         dimnames = list(NULL, paper_groups))
  pooled_values <- numeric(rounds)

  for (round in seq_len(rounds)) {
    pooled_time <- pooled_status <- pooled_score <- numeric(0L)
    for (group in seq_along(paper_bits)) {
      index <- which(bits == paper_bits[group])
      if (length(index) <= k) next
      folds <- fold_ids(length(index), k)
      group_score <- rep(NA_real_, length(index))
      for (fold in seq_len(k)) {
        train <- index[folds != fold]
        test <- index[folds == fold]
        train_data <- data.frame(
          outcome[train, c("time", "status")],
          covariates[train, covariate_names, drop = FALSE],
          check.names = FALSE
        )
        fit_cph <- try(survival::coxph(
          survival::Surv(time, status) ~ ., data = train_data,
          ties = "breslow"
        ), silent = TRUE)
        if (inherits(fit_cph, "try-error")) next
        group_score[folds == fold] <- as.numeric(stats::predict(
          fit_cph,
          newdata = covariates[test, covariate_names, drop = FALSE],
          type = "lp"
        ))
      }
      group_values[round, group] <- cindex_from_score(
        outcome$time[index], outcome$status[index], group_score
      )
      pooled_time <- c(pooled_time, outcome$time[index])
      pooled_status <- c(pooled_status, outcome$status[index])
      pooled_score <- c(pooled_score, group_score)
    }
    pooled_values[round] <- cindex_from_score(
      pooled_time, pooled_status, pooled_score
    )
  }
  summarise_cv(cbind(group_values, Full = pooled_values))
}

run_imr_cv <- function(dat, method = c("IMR", "BMS"), use_covariates = TRUE,
                       mcmc = c(400000L, 50000L), k = 10L, rounds = 10L,
                       seed = 1L) {
  method <- match.arg(method)
  outcome <- if (!is.null(dat$outcome.survival)) dat$outcome.survival else dat$outcome
  model <- imr(
    platform_data_list = dat$platforms,
    outcome = outcome,
    cov = if (use_covariates) dat$covariates else NULL,
    # This optional path reproduces the explicitly labelled legacy appendix.
    survival_scale = "identity",
    type_outcome = "right.censored",
    method = method,
    nu = c(-4, -3, -4),
    h0 = 10000,
    hh = 0.087,
    sig_alpha_psi = c(0.001, 0.001),
    thet_alph_bet = c(40, 10),
    sample_mcmc = mcmc,
    ssize = 30,
    seed = seed
  )
  cv <- cv_imr(model, k = k, rounds = rounds, method = method)
  paper_bits <- c("111", "011", "101", "001")
  keep <- c(intersect(paper_bits, colnames(cv$total_cindex)), "all")
  values <- cv$total_cindex[, keep, drop = FALSE]
  colnames(values) <- c(c("E1", "E2", "E3", "E5")[
    match(keep[-length(keep)], paper_bits)], "Full")
  summarise_cv(values)
}

format_mean_sd <- function(mean, sd) {
  if (is.na(sd)) sprintf("%.3f", mean) else sprintf("%.3f (%.3f)", mean, sd)
}

assemble_kirc_table <- function(results) {
  means <- do.call(rbind, lapply(results, `[[`, "mean"))
  sds <- do.call(rbind, lapply(results, `[[`, "sd"))
  cells <- matrix("", nrow = nrow(means), ncol = ncol(means),
                  dimnames = dimnames(means))
  for (i in seq_len(nrow(means))) {
    for (j in seq_len(ncol(means))) {
      cells[i, j] <- format_mean_sd(means[i, j], sds[i, j])
    }
  }
  list(cells = cells, mean = means, sd = sds)
}

load_kirc_data <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot find ", path, ". See paper/data/README.md for preparation ",
         "instructions, or omit --full-kirc/--smoke-kirc.")
  }
  environment <- new.env(parent = emptyenv())
  load(path, envir = environment)
  candidates <- c("table1_kirc", "kirc_table1", "kirc_full", "kircIMR")
  found <- candidates[candidates %in% ls(environment)]
  if (!length(found)) {
    stop("Expected one of these objects in ", path, ": ",
         paste(candidates, collapse = ", "))
  }
  get(found[1L], envir = environment)
}

kirc_table <- NULL
if (run_full_kirc || run_smoke_kirc) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for the KIRC Table 1 workflow.")
  }
  kirc_data <- load_kirc_data(full_data_path)
  kirc_mcmc <- if (run_smoke_kirc) c(50L, 20L) else c(400000L, 50000L)
  kirc_k <- if (run_smoke_kirc) 2L else 10L
  kirc_rounds <- if (run_smoke_kirc) 1L else 10L

  cat("\n=== Historical TCGA-KIRC Table 1 package-based reanalysis ===\n")
  cat("MCMC retained/burn-in:", paste(kirc_mcmc, collapse = "/"),
      "; CV:", kirc_rounds, "x", kirc_k, "folds\n")

  kirc_results <- list()
  cat("Running CPH + C\n")
  kirc_results[["CPH + C"]] <- run_cph_cv(
    kirc_data, kirc_k, kirc_rounds, seed = 1L
  )
  cat("Running IMR + C + M\n")
  kirc_results[["IMR + C + M"]] <- run_imr_cv(
    kirc_data, "IMR", TRUE, kirc_mcmc, kirc_k, kirc_rounds, 1L
  )
  cat("Running BMS + C + M\n")
  kirc_results[["BMS + C + M"]] <- run_imr_cv(
    kirc_data, "BMS", TRUE, kirc_mcmc, kirc_k, kirc_rounds, 1L
  )
  cat("Running IMR + M\n")
  kirc_results[["IMR + M"]] <- run_imr_cv(
    kirc_data, "IMR", FALSE, kirc_mcmc, kirc_k, kirc_rounds, 1L
  )
  cat("Running BMS + M\n")
  kirc_results[["BMS + M"]] <- run_imr_cv(
    kirc_data, "BMS", FALSE, kirc_mcmc, kirc_k, kirc_rounds, 1L
  )
  kirc_table <- assemble_kirc_table(kirc_results)
  print(kirc_table$cells, quote = FALSE)

  kirc_cells <- data.frame(model = rownames(kirc_table$cells),
                           kirc_table$cells, row.names = NULL,
                           check.names = FALSE)
  write_table(kirc_cells, if (run_smoke_kirc) {
    "table_kirc_supplement_smoke"
  } else {
    "table_kirc_supplement"
  })
  saveRDS(kirc_table, file.path(out_dir, "table_kirc_supplement.rds"))
}

## --------------------------------------------------------------------------
## Save one machine-readable object and the computational environment
## --------------------------------------------------------------------------

saveRDS(
  list(
    quick = quick,
    package_version = as.character(utils::packageVersion("IntegMultiReg")),
    computational_scale = computational_scale,
    complete_case = complete_case,
    kirc_dimensions = kirc_dimensions,
    kirc_subgroups = kirc_subgroups,
    simulated_fit = fit,
    recovery = recovery,
    comparison = comparison,
    prediction_table = prediction_table,
    kirc_table = kirc_table
  ),
  file.path(out_dir, "replication_results.rds")
)

writeLines(capture.output(utils::sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))

cat("\nReplication script completed successfully.\n")
cat("Results are in:", out_dir, "\n")
