# This file is invoked from the frozen study source, never from a live runner.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
study <- normalizePath(args[1], mustWork = TRUE)
id <- as.integer(args[2])
plan <- readRDS(file.path(study, "study.rds"))
stopifnot(length(id) == 1L, !is.na(id), id >= 1L, id <= nrow(plan$tasks),
  identical(tools::md5sum(names(plan$hashes)), plan$hashes),
  identical(R.version, plan$R_version), identical(RNGkind(), plan$rng_kind),
  identical(system2("gsl-config", "--version", stdout = TRUE), plan$gsl_version))
.libPaths(c(file.path(plan$candidate, "library"), plan$dependencies, .libPaths()))
stopifnot(normalizePath(find.package("IntegMultiReg")) ==
          normalizePath(file.path(plan$candidate, "library/IntegMultiReg")))
source(file.path(plan$source, "hpc/study-plan.R"))
task <- plan$tasks[id, , drop = FALSE]
directory <- file.path(study, "tasks", sprintf("%04d", id))
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
identity_path <- file.path(directory, "identity.rds")
if (file.exists(identity_path)) {
  stopifnot(identical(readRDS(identity_path), task))
} else saveRDS(task, identity_path)
command <- study_task_arguments(task, directory)
Sys.setenv(R_LIBS = paste(c(file.path(plan$candidate, "library"), plan$dependencies),
                         collapse = .Platform$path.sep),
           R_LIBS_USER = file.path(plan$candidate, "library"))
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(file.path(plan$source, "paper", command$script)),
    vapply(command$args, shQuote, "")))
if (status != 0L) quit(status = status)
if (task$kind == "audit") stopifnot(file.exists(file.path(directory, "AUDIT-COMPLETED")))
writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(directory, "TASK-COMPLETED"))
