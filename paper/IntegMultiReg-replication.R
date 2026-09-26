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
## Current manuscript workflow: figures, six root-level tables, the reduced
## KIRC interface example, and paired/nested covariate comparison tables.
##   Rscript IntegMultiReg-replication.R
##   Rscript IntegMultiReg-replication.R --quick
## Requires IntegMultiReg 0.2.0 and the sibling covariate-comparison script.
## Legacy raw-time KIRC calculations are archived separately; this script does
## not produce or validate the historical long-chain KIRC table.
###############################################################################

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) any(args == flag)

option_value <- function(flag, default = NULL) {
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
if (has_flag("--full-kirc") || has_flag("--smoke-kirc")) {
  stop("Legacy KIRC analysis is archived separately in historical-audit; it is not a current performance workflow.")
}
run_reduced_kirc_fit <- !has_flag("--skip-reduced-kirc-fit")
update_reference <- has_flag("--update-reference")
reference_dir <- option_value("--reference-dir", file.path(materials_root, "expected-results-0.2.0-paper"))
if (update_reference && quick) stop("Reference tables require a full run, not --quick.")

default_output <- file.path(materials_root, "replication-output")
if (quick) default_output <- paste0(default_output, "-quick")
out_dir <- normalizePath(option_value("--out-dir", default_output),
                         mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("IntegMultiReg", quietly = TRUE)) {
  stop("Install IntegMultiReg before running this script. For a local source ",
       "checkout, run: R CMD INSTALL IntegMultiReg")
}
if (utils::packageVersion("IntegMultiReg") != "0.2.0") {
  stop("This replication script requires exactly IntegMultiReg 0.2.0; installed: ",
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

computational_scale <- data.frame(
  workflow = c("simIMR figures and prediction table", "kircIMR reduced example",
               "covariate comparison and nested selection"),
  subjects = c(length(all_ids(simIMR)), sum(kircIMR$model_subgroup_sizes), 300L),
  features = c(number_of_features(simIMR), number_of_features(kircIMR), 38L),
  subgroups = c(length(sim_modelled_subgroups), length(kircIMR$model_subgroup_sizes), 4L),
  configuration = c(
    paste("Figures: 8000 retained + 2000 burn-in; prediction:",
          "4000 retained + 1000 burn-in, 10 x 5-fold CV"),
    "4000 retained + 1000 burn-in for the illustrative fit",
    "2000 retained + 500 burn-in; paired 2 x 3-fold CV; 3 outer x 3 inner folds"
  ), check.names = FALSE
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
    x = kircIMR$platforms,
    outcome = kircIMR$outcome.survival,
    covariates = kircIMR$covariates,
    outcome_type = "right.censored",
    nu = c(-4, -3, -4),
    draws = reduced_mcmc[1], burnin = reduced_mcmc[2],
    min_subgroup_size = 30,
    seed = 1, sampler_method = "paper"
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
  outcome_type = "right.censored"
)
fit <- imr(
  x = analysis_data,
  nu = nu,
  draws = figure_mcmc[1], burnin = figure_mcmc[2],
  min_subgroup_size = 30,
  seed = 1, sampler_method = "paper"
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
  result <- lapply(c("imr", "bms"), function(method) {
    cat("Starting training-fold CV:", type, method, "\n")
    model <- imr(
      x = simIMR$platforms,
      outcome = outcomes[[type]],
      covariates = simIMR$covariates,
      outcome_type = type,
      method = method,
      nu = nu,
      draws = prediction_mcmc[1], burnin = prediction_mcmc[2],
      min_subgroup_size = 30,
      seed = 1, sampler_method = "paper"
    )
    cv <- cv_imr(model, k = cv_k, rounds = cv_rounds, cv_method = "refit")
    cat("Completed training-fold CV:", type, method, "\n")
    saveRDS(list(fit = model, cv = cv), file.path(out_dir,
      paste0("validation-", type, "-", method, ".rds")))
    colMeans(cv$pooled)
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
    rbind(value$IMR, value$BMS),
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
      stop("No 0.2.0 reference exists. Run once with --update-reference, then verify independently.")
    }
    expected_prediction <- utils::read.csv(reference_path, check.names = FALSE)
    if (!identical(names(prediction_table), names(expected_prediction)) ||
        !identical(prediction_table[, 1:3], expected_prediction[, 1:3]) ||
        !isTRUE(all.equal(as.matrix(prediction_table[, -(1:3)]),
                          as.matrix(expected_prediction[, -(1:3)]), tolerance = 1e-12))) {
      stop("Prediction values differ from the 0.2.0 reference. Check source, settings and sessionInfo().")
    }
  }
}

## --------------------------------------------------------------------------
## Paired predictive comparison and nested selection of covariate formulas
## --------------------------------------------------------------------------
local({
  previous_wd <- getwd()
  on.exit(setwd(previous_wd))
  setwd(materials_root)
  source("IntegMultiReg-covariate-comparison.R", local = .GlobalEnv)
})
covariate_comparison <- run_covariate_comparison(
  out_dir = file.path(out_dir, "covariate-comparison"), quick = quick)
print(covariate_comparison$paired_summary, digits = 6)
print(covariate_comparison$nested_summary, digits = 6)

## --------------------------------------------------------------------------
## Save one machine-readable object and the computational environment
## --------------------------------------------------------------------------

saveRDS(
  list(
    quick = quick, sampler_method = "paper", cv_method = "refit",
    package_version = as.character(utils::packageVersion("IntegMultiReg")),
    computational_scale = computational_scale,
    complete_case = complete_case,
    kirc_dimensions = kirc_dimensions,
    kirc_subgroups = kirc_subgroups,
    simulated_fit = fit,
    recovery = recovery,
    comparison = comparison,
    prediction_table = prediction_table,
    covariate_comparison = covariate_comparison
  ),
  file.path(out_dir, "replication_results.rds")
)

writeLines(capture.output(utils::sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))

cat("\nReplication script completed successfully.\n")
cat("Results are in:", out_dir, "\n")
