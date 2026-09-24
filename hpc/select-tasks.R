# Selection is separate from the immutable study settings.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, args[2] %in% c("audit", "pilot", "full"))
study <- normalizePath(args[1], mustWork = TRUE)
tasks <- readRDS(file.path(study, "study.rds"))$tasks
if (args[2] != "audit") stopifnot(file.exists(file.path(study, "AUDITS-ACCEPTED")))
if (args[2] == "full") stopifnot(file.exists(file.path(study, "PILOT-ACCEPTED")))
keep <- switch(args[2], audit = tasks$kind == "audit",
  pilot = tasks$kind != "audit" & tasks$pilot,
  full = tasks$kind != "audit")
complete <- file.exists(file.path(study, "tasks", sprintf("%04d", tasks$task_id), "TASK-COMPLETED"))
cat(paste(tasks$task_id[keep & !complete], collapse = ","))
