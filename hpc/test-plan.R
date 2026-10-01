script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                   value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "study-plan.R"))
tasks <- make_study_plan()
stopifnot(nrow(tasks) == 679L, identical(tasks$task_id, seq_len(679L)),
  sum(tasks$kind == "audit") == 10L, sum(tasks$pilot) == 20L,
  sum(tasks$kind == "chain") == 8L, sum(tasks$kind == "table1") == 1L)
simulation <- tasks[tasks$kind == "simulation", ]
stopifnot(!anyDuplicated(simulation[c("experiment", "reference", "configuration", "replicate")]))
for (reference in c("paper", "code2017")) {
  rows <- simulation[simulation$reference == reference, ]
  stopifnot(sum(rows$experiment == "simulation") == 150L,
            sum(rows$experiment == "correlated") == 180L,
            all(rows$seed == 99L + rows$replicate))
  stopifnot(all(rows$ridge == if (reference == "paper") .001 else 0),
            all(rows$marker_design == if (reference == "paper") "random" else "fixed"))
}
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  args <- study_task_arguments(task, "directory with spaces")$args
  stopifnot(tail(args, 1) == "directory with spaces" || "--out-dir" %in% args,
            identical("--audit" %in% args, task$kind == "audit"))
  if (task$kind != "chain") {
    stopifnot(as.integer(args[match("--replicate",
                                    args) + 1L]) == task$replicate,
              as.integer(args[match("--replicates",
                                    args) + 1L]) == task$replicates,
              !"--chain" %in% args)
  }
}
source(file.path(dirname(dirname(normalizePath(script))),
                 "paper/original-experiment-helpers.R"))
strict <- resolve_experiment_arguments("paper", list())
paper_ridge_arguments <- resolve_experiment_arguments("paper", list("--ridge" = "0.001"))
stopifnot(strict$cv$ridge == 0, paper_ridge_arguments$cv$ridge == .001,
          identical(strict$fit, paper_ridge_arguments$fit),
          identical(reference_arguments("code2017"),
                    resolve_experiment_arguments("code2017", list())))
for (options in list(list("--ridge" = "NaN"), list("--ridge" = "-1"),
                    list("--cv-method" = "refit", "--ridge" = "0"))) {
  failure <- tryCatch(resolve_experiment_arguments("paper", options),
                      error = identity)
  stopifnot(inherits(failure, "error"))
}
cat("PASS: complete task coverage, selectors, seed independence and argument overrides.\n")
