# Rscript hpc/collect.R STUDY [--simulation-validator FILE]
# The optional, provenance-recorded validator replaces validation only, never data.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 3L))
simulation_validator <- NULL
if (length(args) == 3L) {
  stopifnot(args[2] == "--simulation-validator")
  simulation_validator <- normalizePath(args[3], mustWork = TRUE)
}
study <- normalizePath(args[1], mustWork = TRUE)
source(file.path(study, "source/hpc/study-plan.R"))
plan <- read_study_plan(study)
stopifnot(file.exists(file.path(plan$candidate, "NATIVE-VALIDATED")),
          identical(tools::md5sum(names(plan$hashes)), plan$hashes))
source(file.path(plan$source, "hpc/collection.R"))
.libPaths(c(file.path(plan$candidate, "library"), plan$dependencies, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
run_validator <- function(script, arguments) {
  validator <- if (identical(script, "paper/validate-original-simulation.R") &&
                   !is.null(simulation_validator)) simulation_validator else file.path(plan$source, script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(validator), vapply(arguments, shQuote, "")))
  if (status != 0L) stop("Validation failed: ", script)
}
tasks <- plan$tasks[plan$tasks$kind != "audit", ]
stopifnot(all(file.exists(file.path(study, "tasks", sprintf("%04d", tasks$task_id), "TASK-COMPLETED"))))
key <- paste(tasks$experiment, tasks$reference, ifelse(is.na(tasks$ridge), "none", tasks$ridge), sep = "-")
for (group in split(tasks, key)) {
  paths <- file.path(study, "tasks", sprintf("%04d", group$task_id))
  experiment <- group$experiment[1]
  if (experiment == "table1") {
    run_validator("paper/validate-original-table1.R", c(paths, "--require-original"))
    next
  }
  out <- file.path(study, "combined", paste(experiment, group$reference[1],
    ifelse(is.na(group$ridge[1]), "none", group$ridge[1]), sep = "-"))
  collect_tasks(experiment, paths, out)
  if (experiment == "chains") {
    run_validator("hpc/diagnose-chains.R", out)
    run_validator("paper/rank-chain-diagnostics.R", c(out, file.path(out, "diagnostics")))
    run_validator("paper/chain-ranking-stability.R", c(out, file.path(out, "diagnostics")))
    run_validator("hpc/check-halves.R", c(out, file.path(out, "half-drift")))
  } else {
    run_validator("paper/validate-original-simulation.R", c(out, "--require-complete-study"))
    run_validator("paper/validate-simulation-baselines.R", out)
    run_validator("paper/summarize-original-simulation.R", out)
  }
}
writeLines("Program checks complete; convergence, warnings and interpretation require scientific review.",
           file.path(study, "VALIDATION-COMPLETED"))
