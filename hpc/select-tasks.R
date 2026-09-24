# Selection is separate from the immutable study settings.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% 2:3, args[2] %in% c("audit", "pilot", "full"))
study <- normalizePath(args[1], mustWork = TRUE)
source(file.path(study, "source/hpc/study-plan.R"))
tasks <- read_study_plan(study)$tasks
if (args[2] != "audit") stopifnot(file.exists(file.path(study, "AUDITS-ACCEPTED")))
if (args[2] == "full") stopifnot(file.exists(file.path(study, "PILOT-ACCEPTED")))
keep <- switch(args[2], audit = tasks$kind == "audit",
  pilot = tasks$kind != "audit" & tasks$pilot,
  full = tasks$kind != "audit")
complete <- file.exists(file.path(study, "tasks", sprintf("%04d", tasks$task_id), "TASK-COMPLETED"))
ids <- tasks$task_id[keep & !complete]
if (length(args) == 3L) {
  values <- strsplit(args[3], ",", fixed = TRUE)[[1]]
  stopifnot(length(values) > 0L, all(grepl("^[1-9][0-9]*$", values)))
  requested <- suppressWarnings(as.integer(values))
  stopifnot(!anyNA(requested), !anyDuplicated(requested), all(requested %in% ids))
  ids <- requested
}
cat(paste(ids, collapse = ","), "\n", sep = "")
