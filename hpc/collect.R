# Rscript hpc/collect.R STUDY -- executes existing scientific validators after collection.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
study <- normalizePath(args[1], mustWork = TRUE)
source(file.path(study, "source/hpc/study-plan.R"))
plan <- read_study_plan(study)
stopifnot(identical(tools::md5sum(names(plan$hashes)), plan$hashes))
source(file.path(plan$source, "hpc/collection.R"))
.libPaths(c(file.path(plan$candidate, "library"), plan$dependencies, .libPaths()))
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
run <- function(script, arguments) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(plan$source, script)), vapply(arguments, shQuote, "")))
  if (status != 0L) stop("Validation failed: ", script)
}
tasks <- plan$tasks[plan$tasks$kind != "audit", ]
stopifnot(all(file.exists(file.path(study, "tasks", sprintf("%04d", tasks$task_id), "TASK-COMPLETED"))))
key <- paste(tasks$experiment, tasks$reference, ifelse(is.na(tasks$ridge), "none", tasks$ridge), sep = "-")
for (group in split(tasks, key)) {
  paths <- file.path(study, "tasks", sprintf("%04d", group$task_id))
  experiment <- group$experiment[1]
  if (experiment == "table1") {
    run("paper/validate-original-table1.R", c(paths, "--require-original"))
    next
  }
  out <- file.path(study, "combined", paste(experiment, group$reference[1],
    ifelse(is.na(group$ridge[1]), "none", group$ridge[1]), sep = "-"))
  collect_tasks(experiment, paths, out)
  if (experiment == "chains") {
    run("hpc/diagnose-chains.R", out)
    run("paper/rank-chain-diagnostics.R", c(out, file.path(out, "diagnostics")))
    run("paper/chain-ranking-stability.R", c(out, file.path(out, "diagnostics")))
    run("hpc/check-halves.R", c(out, file.path(out, "half-drift")))
  } else {
    run("paper/validate-original-simulation.R", c(out, "--require-complete-study"))
    run("paper/validate-simulation-baselines.R", out)
    run("paper/summarize-original-simulation.R", out)
  }
}
writeLines("Program checks complete; convergence, warnings and interpretation require scientific review.",
           file.path(study, "VALIDATION-COMPLETED"))
