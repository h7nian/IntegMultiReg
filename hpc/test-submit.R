# Exercise the real submitter with mock Slurm commands: no jobs are submitted.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
directory <- dirname(script)
source(file.path(directory, "study-plan.R"))
study <- tempfile("submit-check-")
dir.create(file.path(study, "source/hpc"), recursive = TRUE)
dir.create(file.path(study, "bin"))
for (name in c("study-plan.R", "select-tasks.R", "submit.sh"))
  stopifnot(file.copy(file.path(directory, name), file.path(study, "source/hpc")))
writeLines(c(paste0("IMR_ROOT=", shQuote(file.path(study, "source"))),
             paste("source", shQuote(file.path(directory, "config.sh")))),
           file.path(study, "source/hpc/environment.sh"))
candidate <- file.path(study, "candidate")
dir.create(candidate)
saveRDS(list(tasks = make_study_plan(), candidate = candidate,
             native_validation_job = "99999"), file.path(study, "study.rds"))
writeLines(unname(tools::md5sum(file.path(study, "study.rds"))), file.path(study, "MANIFEST.md5"))
writeLines(c("#!/bin/sh", 'printf "%s\\n" "$*" >> "$IMR_SUBMIT_TEST_LOG"',
             'printf "12345\\n"'), file.path(study, "bin/sbatch"))
writeLines(c("#!/bin/sh",
             'if [ "${IMR_MOCK_ACTIVE:-0}" = 1 ]; then printf "12345\\n"; fi'),
           file.path(study, "bin/squeue"))
Sys.chmod(file.path(study, "bin", c("sbatch", "squeue")), "0755")
calls <- file.path(study, "calls.log")
submit <- function(success, ...) {
  overrides <- list(...)
  environment <- c(paste0("PATH=", shQuote(paste(file.path(study, "bin"), Sys.getenv("PATH"), sep = ":"))),
    paste0("IMR_SUBMIT_TEST_LOG=", shQuote(calls)),
    vapply(names(overrides), function(name) paste0(name, "=", shQuote(overrides[[name]])), ""))
  log <- tempfile(tmpdir = study)
  status <- system2("bash", c(shQuote(file.path(directory, "submit.sh")),
                              shQuote(study), "pilot"),
                    env = environment, stdout = log, stderr = log)
  if (!identical(status == 0L, success)) stop(paste(readLines(log), collapse = "\n"))
}
submit(FALSE)
stopifnot(!file.exists(calls))
writeLines("test fixture", file.path(study, "AUDITS-ACCEPTED"))
submit(TRUE)
arguments <- readLines(calls)
stopifnot(length(arguments) == 1L, grepl("%16", arguments, fixed = TRUE),
          grepl("--mem=16384M", arguments, fixed = TRUE),
          grepl("--cpus-per-task=9", arguments, fixed = TRUE),
          grepl("--time=48:00:00", arguments, fixed = TRUE))
submit(FALSE, IMR_MOCK_ACTIVE = "1")
submit(FALSE, IMR_TASK_IDS = "670")
submit(FALSE, IMR_MEMORY_MIB = "300000")
submit(FALSE, IMR_WALL_HOURS = "97")
stopifnot(length(readLines(calls)) == 1L)
submit(TRUE, IMR_TASK_IDS = "1", IMR_MEMORY_MIB = "32768",
       IMR_WALL_HOURS = "96")
arguments <- tail(readLines(calls), 1L)
stopifnot(grepl("--array=1%16", arguments, fixed = TRUE),
          grepl("--cpus-per-task=18", arguments, fixed = TRUE))
# Run the actual controller/submitter, replacing only phase acceptance and Slurm.
# The experiment acceptance code is covered separately by its scientific tests.
writeLines(c("#!/bin/sh", 'case "$2" in */accept-phase.R) exit 0 ;; esac',
             paste("exec", shQuote(file.path(R.home("bin"), "Rscript")),
                   '"$@"')),
           file.path(study, "bin/Rscript"))
Sys.chmod(file.path(study, "bin/Rscript"), "0755")
advance <- function() {
  log <- tempfile(tmpdir = study)
  status <- system2("bash", c(shQuote(file.path(directory, "advance.sh")),
                              shQuote(study), "audit"),
    env = c("SLURM_JOB_ID=67890", paste0("IMR_SUBMIT_TEST_LOG=", shQuote(calls)),
      paste0("PATH=", shQuote(paste(file.path(study, "bin"), Sys.getenv("PATH"),
                                    sep = ":")))),
    stdout = log, stderr = log)
  if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
  tail(readLines(calls), 1L)
}
stopifnot(grepl("--dependency=afterany:12345,afterok:99999", advance(),
                fixed = TRUE))
writeLines("test fixture", file.path(candidate, "NATIVE-VALIDATED"))
arguments <- advance()
stopifnot(grepl("--dependency=afterany:12345", arguments, fixed = TRUE),
          !grepl("afterok:", arguments, fixed = TRUE))
unlink(study, recursive = TRUE)
cat("PASS: submit gates, active-array refusal, global throttle, resource bounds and pending-native controller dependencies.\n")
