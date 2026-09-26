# Exercise preparation gates on disposable source/data/validation fixtures.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
directory <- dirname(script)
fixture <- tempfile("concurrent-start-")
root <- file.path(fixture, "root")
candidate <- file.path(fixture, "candidate")
dependencies <- file.path(fixture, "dependencies")
for (path in c("hpc", "paper/reference", "paper/data", "Ref"))
  dir.create(file.path(root, path), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(candidate, "evidence"), recursive = TRUE)
dir.create(dependencies)
for (name in c("study-plan.R", "config.sh"))
  stopifnot(file.copy(file.path(directory, name), file.path(root, "hpc")))
for (name in c("paper/reference/original-generator.c", "paper/data/kirc_table1_full.rda",
               "Ref/biom12587-sup-0002-suppdata_code.zip"))
  writeLines("preparation test fixture, never used for fitting", file.path(root, name))
for (name in c("UNIT-VALIDATED", "SANITIZER-VALIDATED"))
  writeLines("test fixture", file.path(candidate, name))
environment <- c(paste0("IMR_ROOT=", shQuote(root)),
  paste0("IMR_DEPENDENCIES=", shQuote(dependencies)), "IMR_MEMORY_PER_CPU_MIB=1896",
  "IMR_MAX_TASK_CPUS=128", "IMR_MAX_WALL_HOURS=96")
run_fixture_script <- function(name, arguments, success) {
  log <- tempfile(tmpdir = fixture)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(directory, name)), vapply(arguments, shQuote, "")),
    env = environment, stdout = log, stderr = log)
  if (!identical(status == 0L, success)) stop(paste(readLines(log), collapse = "\n"))
}
run_fixture_script("check-research-source.R", c("capture", file.path(candidate, "evidence/research-source-hashes.rds")), TRUE)
prepare <- function(success, options = character()) {
  out <- tempfile("study-", tmpdir = fixture)
  run_fixture_script("prepare-study.R", c(candidate, out, options), success)
  if (success) readRDS(file.path(out, "study.rds")) else stopifnot(!dir.exists(out))
}
prepare(FALSE)
prepare(FALSE, c("--pending-native-job", "invalid"))
prepare(FALSE, c("--unknown-option", "99999"))
plan <- prepare(TRUE, c("--pending-native-job", "99999"))
stopifnot(!plan$native_validated_at_start, identical(plan$native_validation_job, "99999"),
          nrow(plan$tasks) == 679L, identical(tools::md5sum(names(plan$hashes)), plan$hashes),
          !file.exists(file.path(candidate, "NATIVE-VALIDATED")))
unlink(file.path(candidate, "SANITIZER-VALIDATED"))
prepare(FALSE, c("--pending-native-job", "99999"))
writeLines("test fixture", file.path(candidate, "SANITIZER-VALIDATED"))
writeLines("test fixture", file.path(candidate, "NATIVE-VALIDATED"))
plan <- prepare(TRUE)
stopifnot(plan$native_validated_at_start, is.null(plan$native_validation_job))
cat("changed", file = file.path(root, "paper/reference/original-generator.c"), append = TRUE)
prepare(FALSE, c("--pending-native-job", "99999"))
unlink(fixture, recursive = TRUE)
cat("PASS: explicit concurrent preparation, unchanged default gate, required sanitizer evidence and source-integrity refusal.\n")
