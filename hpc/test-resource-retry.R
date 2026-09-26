# Accounting is mocked; only the retry plan is exercised, never sbatch.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
directory <- dirname(script)
source(file.path(directory, "study-plan.R"))
study <- tempfile("retry-check-")
dir.create(file.path(study, "source/hpc"), recursive = TRUE)
dir.create(file.path(study, "bin"))
stopifnot(file.copy(file.path(directory, "study-plan.R"), file.path(study, "source/hpc")))
saveRDS(list(tasks = make_study_plan(), resource_limits = list(
  max_task_cpus = 128L, memory_per_cpu_mib = 1896L, max_wall_hours = 96L)),
  file.path(study, "study.rds"))
writeLines(unname(tools::md5sum(file.path(study, "study.rds"))), file.path(study, "MANIFEST.md5"))
writeLines("12345\tpilot\t16384\t1,9\t48", file.path(study, "submissions.tsv"))
writeLines(c("#!/bin/sh", 'cat "$IMR_ACCOUNTING_FIXTURE"'), file.path(study, "bin/sacct"))
Sys.chmod(file.path(study, "bin/sacct"), "0755")
fixture <- file.path(study, "accounting.txt")
run_with_accounting <- function(rows, success) {
  writeLines(rows, fixture)
  log <- tempfile(tmpdir = study)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(directory, "resource-retry.R")), shQuote(study), "pilot"),
    env = c(paste0("PATH=", shQuote(paste(file.path(study, "bin"), Sys.getenv("PATH"), sep = ":"))),
            paste0("IMR_ACCOUNTING_FIXTURE=", shQuote(fixture))), stdout = log, stderr = log)
  stopifnot(identical(status == 0L, success))
  if (success) readLines(log)
}
observed <- run_with_accounting(c("12345_1|OUT_OF_MEMORY|", "12345_9|FAILED|"), TRUE)
if (!identical(observed, c("1", "32768", "48"))) stop(paste(observed, collapse = " | "))
run_with_accounting(c("12345_1|OUT_OF_MEMORY|", "12345_9|FAILED|"), FALSE)
stopifnot(identical(run_with_accounting(c("12345_1|COMPLETED|", "12345_9|TIMEOUT|"), TRUE),
                    c("9", "16384", "96")))
run_with_accounting(c("12345_1|CANCELLED|", "12345_9|TIMEOUT|"), FALSE)
unlink(study, recursive = TRUE)
cat("PASS: bounded OOM/timeout retries; no numerical-failure or duplicate retries.\n")
