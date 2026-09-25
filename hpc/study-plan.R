# One row is one independently restartable Slurm array task.
read_study_plan <- function(study) {
  path <- file.path(study, "study.rds")
  stopifnot(identical(unname(tools::md5sum(path)),
                      readLines(file.path(study, "MANIFEST.md5"))))
  readRDS(path)
}

# Audits and the pilot may overlap native checks; full production must wait.
native_validation_dependency <- function(plan) {
  if (file.exists(file.path(plan$candidate, "NATIVE-VALIDATED"))) return("")
  job <- plan$native_validation_job
  stopifnot(is.character(job), length(job) == 1L, !is.na(job),
            grepl("^[1-9][0-9]*$", job))
  paste0("afterok:", job)
}

make_study_plan <- function() {
  rows <- list()
  add <- function(kind, experiment, reference, ridge, marker_design,
                  configuration = 1L, replicate = 1L, replicates = 1L,
                  chain = NA_integer_, pilot = FALSE) {
    rows[[length(rows) + 1L]] <<- data.frame(task_id = length(rows) + 1L,
      kind, experiment, reference, ridge, marker_design,
      configuration, replicate, replicates, chain, pilot,
      seed = if (kind == "chain") 99L + chain else 99L + replicate,
      stringsAsFactors = FALSE)
  }
  for (chain in 1:8) add("chain", "chains", "paper", NA_real_, NA_character_,
    chain = chain, pilot = chain == 1L)
  for (reference in c("code2017", "paper")) {
    for (experiment in c("simulation", "correlated")) {
      replicates <- if (experiment == "simulation") 50L else 30L
      configurations <- if (experiment == "simulation") 3L else 6L
      for (configuration in seq_len(configurations)) for (replicate in seq_len(replicates))
        add("simulation", experiment, reference,
          if (reference == "paper") .001 else 0,
          if (reference == "paper") "random" else "fixed",
          configuration, replicate, replicates, pilot = replicate == 1L)
    }
  }
  add("table1", "table1", "code2017", 0, NA_character_, pilot = TRUE)
  for (experiment in c("simulation", "correlated", "table1")) {
    configurations <- switch(experiment, simulation = 3L, correlated = 6L, table1 = 1L)
    for (configuration in seq_len(configurations))
      add("audit", experiment, "paper", 0,
        if (experiment == "table1") NA_character_ else "random", configuration)
  }
  do.call(rbind, rows)
}

study_task_arguments <- function(task, directory) {
  if (task$kind == "chain") return(list(script = "run-appendix-chains.R",
    args = c("--chain", task$chain, "--draws", 1000000L, "--burnin", 50000L,
             "--seed", 100L, "--out-dir", directory)))
  arguments <- c("--experiment", task$experiment, "--reference", task$reference,
    "--ridge", format(task$ridge, scientific = FALSE),
    "--draws", if (task$kind == "audit") 20000L else 350000L,
    "--burnin", if (task$kind == "audit") 5000L else 50000L,
    "--k", 10L, "--rounds", if (task$kind == "table1") 10L else 1L,
    "--replicates", task$replicates, "--configuration", task$configuration,
    "--replicate", task$replicate, "--workers", 1L, "--seed", 100L,
    "--out-dir", directory)
  if (!is.na(task$marker_design)) arguments <- c(arguments, "--marker-design", task$marker_design)
  if (task$kind == "audit") arguments <- c(arguments, "--audit")
  list(script = "original-experiments.R", args = arguments)
}
