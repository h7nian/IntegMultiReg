# Original-scale experiment runner. Run from any directory; outputs are separate
# for quick/full, experiment, reference, configuration, replicate and method.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
materials <- dirname(normalizePath(script))
 root <- dirname(materials)
source(file.path(materials, "original-reference.R"))
source(file.path(materials, "original-experiment-helpers.R"))
options <- parse_experiment_args(commandArgs(TRUE))
if (!is.null(options[["--chain"]])) stop("--chain applies only to the chain runner")
value <- function(flag, default) if (is.null(options[[flag]])) default else options[[flag]]
experiment <- match.arg(value("--experiment", "table1"), c("table1", "simulation", "correlated"))
reference <- match.arg(value("--reference", "paper"), c("paper", "code2017", "package"))
quick <- isTRUE(options[["--quick"]])
out <- value("--out-dir", file.path(root, "output", "original-experiments", if (quick) "quick" else "full", reference, experiment))
integer_value <- function(flag, default, minimum) {
  x <- as.numeric(value(flag, default))
  if (length(x) != 1L || !is.finite(x) || x != floor(x) || x < minimum || x > .Machine$integer.max)
    stop("Invalid ", flag)
  as.integer(x)
}
settings <- list(experiment = experiment, reference = reference, quick = quick,
  draws = integer_value("--draws", if (quick) 50 else 350000, 1),
  burnin = integer_value("--burnin", if (quick) 20 else 50000, 0),
  k = integer_value("--k", if (quick) 2 else 10, 2),
  rounds = integer_value("--rounds", if (quick || experiment != "table1") 1 else 10, 1),
  replicates = integer_value("--replicates", if (quick || experiment == "table1") 1 else
    if (experiment == "correlated") 30 else 50, 1),
  workers = integer_value("--workers", 1, 1),
  seed = integer_value("--seed", 100, 0),
  reference_arguments = resolve_experiment_arguments(reference, options),
  audit = isTRUE(options[["--audit"]]))
if (experiment == "table1" && !is.null(options[["--marker-design"]]))
  stop("--marker-design applies only to simulation or correlated experiments")
settings$marker_design <- if (experiment == "table1") NULL else
  match.arg(value("--marker-design", if (reference == "paper") "random" else "fixed"),
    c("random", "fixed"))
library(IntegMultiReg)
if (!requireNamespace("digest", quietly = TRUE))
  stop("The research checkpoint writer requires the digest package")
if (settings$audit) source(file.path(materials, "audit-unpenalized-cv.R"))
stopifnot(packageVersion("IntegMultiReg") == "0.2.0")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
 out <- normalizePath(out)
data_path <- value("--data", file.path(materials, "data", "kirc_table1_full.rda"))
settings$data_hash <- tools::md5sum(data_path)
settings$source_hashes <- tools::md5sum(c(script,
  file.path(materials, c("original-reference.R", "original-experiment-helpers.R", "reference/original-generator.c")),
  if (settings$audit) file.path(materials, "audit-unpenalized-cv.R"),
  system.file("R", "IntegMultiReg.rdb", package = "IntegMultiReg"),
  system.file("libs", paste0("IntegMultiReg", .Platform$dynlib.ext), package = "IntegMultiReg")))
manifest <- file.path(out, "settings.rds")
if (file.exists(manifest) && !identical(readRDS(manifest), settings))
  stop("Output directory contains different settings; choose a new directory")
saveRDS(settings, manifest)
writeLines(capture.output(dput(settings)), file.path(out, "settings.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
e <- new.env()
 load(data_path, e)
 original <- e$kirc_full
# Old prepared objects retained log time in outcome.survival; always choose raw.
original$outcome <- original$outcome.raw
original$standardize <- TRUE
stopifnot(all(original$outcome$time > 0))
writeLines(c(paste("Data MD5:", tools::md5sum(data_path)),
  paste(names(original$platforms), vapply(original$platforms, ncol, 0L) - 1L),
  "Published Table 1 parentheses are standard errors; save SD and SE separately.",
  "Published seeds/start states are unavailable. These are documented reruns, not historical digit assertions.",
  "Uni-CPH uses all available subjects per platform with Hommel adjustment.",
  "L1-CPH uses lambda.min, 10-fold inner CV (3 in quick mode), independently by subgroup.",
  "L1-CPH full prediction is sample-size-weighted subgroup concordance, as in the article.",
  "The archived code uses fixed marker indices; the article describes random marker selection.",
  paste("Marker design:", if (experiment == "table1") "not applicable" else settings$marker_design),
  "No proprietary IPA network analysis is reproduced."), file.path(out, "provenance.txt"))
configurations <- if (experiment == "table1") data.frame(scenario = 0, rho = 0, half = FALSE) else
  if (experiment == "simulation") data.frame(scenario = 1:3, rho = 0, half = FALSE) else
    expand.grid(scenario = 1, rho = c(.2, .5, .8), half = c(FALSE, TRUE))
if (experiment != "table1") generator <- load_original_generator(
  file.path(root, "Ref", "biom12587-sup-0002-suppdata_code.zip"),
  file.path(out, "reference-generator"), file.path(materials, "reference", "original-generator.c"))
selected_configurations <- if (is.null(options[['--configuration']])) seq_len(nrow(configurations)) else integer_value('--configuration', 1L, 1L)
selected_replicates <- if (is.null(options[['--replicate']])) seq_len(settings$replicates) else integer_value('--replicate', 1L, 1L)
stopifnot(all(selected_configurations <= nrow(configurations)), all(selected_replicates <= settings$replicates))
record_task_selection(out, list(kind = experiment, configurations = selected_configurations, replicates = selected_replicates))
for (configuration in selected_configurations) for (replicate in selected_replicates) {
  config <- configurations[configuration, ]
  seed <- settings$seed + replicate - 1L
  data <- if (experiment == "table1") original else original_simulation(original, generator,
    config$scenario, seed, config$rho, config$half, marker_design = settings$marker_design)
  job_dir <- file.path(out, sprintf("configuration-%02d-replicate-%03d", configuration, replicate))
  dir.create(job_dir, showWarnings = FALSE)
  saveRDS(list(configuration = config, data = data), file.path(job_dir, "data.rds"))
  runs <- if (experiment == "table1") expand.grid(method = c("imr", "bms"), clinical = c(TRUE, FALSE)) else
    data.frame(method = c("imr", "bms"), clinical = FALSE)
  shared_folds <- NULL
  for (run in seq_len(nrow(runs))) {
    method <- as.character(runs$method[run])
     clinical <- runs$clinical[run]
    label <- paste0(method, if (clinical) "-clinical-molecular" else "-molecular")
    path <- file.path(job_dir, paste0(label, ".rds"))
    if (experiment_job_complete(path)) {
      if (is.null(shared_folds)) shared_folds <- experiment_result_folds(path)
      next
    }
    run_method <- function() record_experiment_job(path, {
      cat(format(Sys.time()), configuration, replicate, label, "starting\n")
      elapsed <- system.time({
        fit <- experiment_fit_checkpoint(paste0(path, ".fit.rds"), {
          do.call(imr, c(list(x = data$platforms, outcome = data$outcome,
          covariates = if (clinical) data$covariates else NULL, outcome_type = "right.censored",
          method = method, nu = c(-4, -3, -4), draws = settings$draws, burnin = settings$burnin,
          standardize = data$standardize, seed = seed), settings$reference_arguments$fit))
        })
        cv_args <- settings$reference_arguments$cv
        # Code reproduction keeps each fitted chain's continued stream. Other
        # comparisons use paired partitions after the first model.
        if (!is.null(shared_folds) && !identical(cv_args$fold_rng, "continue")) {
          cv_args$folds <- shared_folds
           cv_args$fold_rng <- NULL
        }
        cv <- do.call(cv_imr, c(list(object = fit, k = settings$k, rounds = settings$rounds,
          workers = settings$workers), cv_args))
        if (is.null(shared_folds)) shared_folds <- cv$control$folds
        scores <- fold_scores(cv$predictions, data$outcome, cv$control$score_method)
        result <- list(fit = fit, cv = cv, scores = scores, summary = summarize_scores(scores),
          selection = if (!is.null(data$truth)) biomarker_auc(fit, data$truth) else NULL)
      })
      result$elapsed <- elapsed
      write_experiment_result(result, path)
      saveRDS(cv$control$folds, paste0(path, ".folds.rds"))
      write.csv(result$summary, sub("[.]rds$", "-scores.csv", path), row.names = FALSE)
      if (!is.null(result$selection)) write.csv(result$selection, sub("[.]rds$", "-selection.csv", path), row.names = FALSE)
      cat(format(Sys.time()), label, "completed; seconds", elapsed[["elapsed"]], "\n")
    })
    if (settings$audit) {
      tryCatch(run_method(), error = function(error) {
        audit_unpenalized_failure(path, settings, error)
      })
    } else run_method()
    if (experiment_job_complete(path)) {
      if (is.null(shared_folds)) shared_folds <- experiment_result_folds(path)
      unlink(paste0(path, ".fit.rds"))
    }
    gc()
  }
  if (settings$audit) next
  baseline_name <- if (experiment == "table1") "cph-clinical" else "l1-cph"
  path <- file.path(job_dir, paste0(baseline_name, ".rds"))
  if (!experiment_job_complete(path)) {
    record_experiment_job(path, {
      baseline <- benchmark_cox(data, shared_folds, experiment != "table1", seed,
        inner_k = if (quick) 3L else 10L)
      baseline$scores <- fold_scores(baseline$predictions, data$outcome, "standard",
        weighted_overall = experiment != "table1")
      baseline$summary <- summarize_scores(baseline$scores)
      saveRDS(baseline, path)
      write.csv(baseline$summary, sub("[.]rds$", "-scores.csv", path), row.names = FALSE)
      if (!is.null(baseline$selection)) write.csv(baseline$selection, sub("[.]rds$", "-selection.csv", path), row.names = FALSE)
    })
  }
  if (!is.null(data$truth)) {
    path <- file.path(job_dir, "uni-cph-selection.csv")
    if (!experiment_job_complete(path)) record_experiment_job(path, {
      write.csv(univariate_cox_auc(data), path, row.names = FALSE)
    })
  }
}
if (settings$audit) {
  writeLines("All requested Bayesian methods audited; inspect per-method status and rank evidence.",
             file.path(out, "AUDIT-COMPLETED"))
} else cat("All requested experiment jobs completed.\n")
