# Acceptance gates are generated from evidence, never from scheduler completion alone.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("audit", "pilot"))
study <- normalizePath(args[1], mustWork = TRUE)
source(file.path(study, "source/hpc/study-plan.R"))
plan <- read_study_plan(study)
stopifnot(identical(tools::md5sum(names(plan$hashes)), plan$hashes))
.libPaths(c(file.path(plan$candidate, "library"), plan$dependencies, .libPaths()))
source(file.path(plan$source, "paper/original-experiment-helpers.R"))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
run_validator <- function(name, arguments) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(plan$source, "paper", name)), vapply(arguments, shQuote, "")))
  if (status != 0L) stop("Pilot validation failed: ", name)
}
tasks <- plan$tasks
tasks <- tasks[if (args[2] == "audit") tasks$kind == "audit" else tasks$pilot, ]
rss <- numeric()
audit <- list()
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  directory <- file.path(study, "tasks", sprintf("%04d", task$task_id))
  stopifnot(file.exists(file.path(directory, "TASK-COMPLETED")),
            identical(readRDS(file.path(directory, "identity.rds")), task))
  settings <- readRDS(file.path(directory, "settings.rds"))
  if (task$kind == "chain") {
    path <- file.path(directory, sprintf("chain-%02d.rds", task$chain))
    stopifnot(experiment_job_complete(path))
    saved <- readRDS(path)
    stopifnot(inherits(saved, "imr_experiment_fit_checkpoint_v1"),
      length(saved$layout) == 1000000L, saved$fit$control$mcmc$burnin == 50000L,
      saved$fit$control$seed == task$seed,
      identical(saved$fit$control$sampler_method, "paper"),
      all(is.finite(saved$fit$posterior$log_posterior)))
    rm(saved)
     gc()
  } else {
    stopifnot(settings$replicates == task$replicates,
      identical(settings$reference, task$reference),
      settings$reference_arguments$cv$ridge == task$ridge)
    job <- file.path(directory, sprintf("configuration-%02d-replicate-%03d",
                                        task$configuration, task$replicate))
    labels <- if (task$experiment == "table1") c("imr-clinical-molecular",
      "bms-clinical-molecular", "imr-molecular", "bms-molecular") else
      c("imr-molecular", "bms-molecular")
    for (label in labels) {
      path <- file.path(job, paste0(label, ".rds"))
      status <- readRDS(paste0(path, ".status.rds"))
      if (identical(status$status, "failed") && task$kind == "audit") {
        rank <- readRDS(paste0(path, ".rank-audit.rds"))
        stopifnot(rank$example$rank < ncol(rank$example$design),
                  sum(rank$report$dimensionally_singular) > 0)
      } else {
        stopifnot(experiment_job_complete(path))
        result <- readRDS(path)
        if (inherits(result, "imr_experiment_result_v1")) result <- result$result
        fit <- if (inherits(result$fit, "imr_experiment_fit_checkpoint_v1")) result$fit$fit else result$fit
        stopifnot(fit$control$mcmc$draws == settings$draws,
          fit$control$mcmc$burnin == settings$burnin, fit$control$seed == task$seed,
          identical(fit$control$sampler_method, settings$reference_arguments$fit$sampler_method),
          result$cv$control$ridge == task$ridge,
          all(is.finite(result$cv$predictions$prediction)),
          all(result$summary$n_valid == result$summary$n_expected))
        rm(result, fit)
         gc()
      }
      audit[[length(audit) + 1L]] <- data.frame(task_id = task$task_id, method = label,
        status = status$status, warnings = length(status$warnings))
    }
    if (task$kind != "audit") {
      baseline <- if (task$experiment == "table1") "cph-clinical.rds" else "l1-cph.rds"
      stopifnot(experiment_job_complete(file.path(job, baseline)))
      if (task$experiment != "table1")
        stopifnot(experiment_job_complete(file.path(job, "uni-cph-selection.csv")))
      if (task$experiment == "table1") {
        run_validator("validate-original-table1.R", c(directory, "--require-original"))
      } else {
        run_validator("validate-original-simulation.R", directory)
        run_validator("validate-simulation-baselines.R", directory)
      }
    }
  }
  resources <- list.files(directory, "[.]resources$", full.names = TRUE)
  stopifnot(length(resources) > 0L)
  resources <- resources[which.max(file.info(resources)$mtime)]
  lines <- readLines(resources)
  stopifnot(any(grepl("Exit status: 0$", lines)))
  peak <- as.numeric(sub(".*: *", "", grep("Maximum resident set size", lines, value = TRUE)))
  stopifnot(length(peak) == 1L, is.finite(peak), peak > 0)
  rss <- c(rss, peak)
}
write.csv(do.call(rbind, audit), file.path(study, paste0(args[2], "-acceptance.csv")), row.names = FALSE)
if (args[2] == "pilot") {
  memory_mib <- max(4096, ceiling(1.5 * max(rss) / 1024 / 2048) * 2048)
  stopifnot(memory_mib <= plan$resource_limits$max_task_cpus * plan$resource_limits$memory_per_cpu_mib)
  writeLines(as.character(memory_mib), file.path(study, "production-memory-mib"))
}
saveRDS(list(tasks = tasks, rss_kib = rss, source_hashes = plan$hashes,
  verified = Sys.time()), file.path(study, paste0(args[2], "-acceptance.rds")))
writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE),
  file.path(study, if (args[2] == "audit") "AUDITS-ACCEPTED" else "PILOT-ACCEPTED"))
