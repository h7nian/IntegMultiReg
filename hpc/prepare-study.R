# Rscript hpc/prepare-study.R CANDIDATE_DIRECTORY NEW_STUDY_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(Sys.getenv("IMR_ROOT"), mustWork = TRUE)
candidate <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
stopifnot(file.exists(file.path(candidate, "UNIT-VALIDATED")),
          file.exists(file.path(candidate, "NATIVE-VALIDATED")), !dir.exists(out))
tested <- readRDS(file.path(candidate, "evidence/research-source-hashes.rds"))
stopifnot(identical(tools::md5sum(names(tested)), tested))
source(file.path(root, "hpc/study-plan.R"))
dir.create(out, recursive = TRUE)
out <- normalizePath(out)
frozen <- file.path(out, "source")
files <- c(list.files(file.path(root, "paper"), "[.]R$", full.names = TRUE),
  list.files(file.path(root, "hpc"), "[.](R|sh|py)$", full.names = TRUE),
  file.path(root, c("paper/reference/original-generator.c", "paper/data/kirc_table1_full.rda",
                  "Ref/biom12587-sup-0002-suppdata_code.zip")))
destinations <- file.path(frozen, substring(files, nchar(root) + 2L))
for (i in seq_along(files)) {
  dir.create(dirname(destinations[i]), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(files[i], destinations[i], copy.mode = TRUE))
}
dependencies <- normalizePath(Sys.getenv("IMR_DEPENDENCIES"), mustWork = TRUE)
cat("\nIMR_DEPENDENCIES=", shQuote(dependencies), "\nIMR_RUNS=",
    shQuote(dirname(out)), "\n", sep = "", append = TRUE,
    file = file.path(frozen, "hpc/config.sh"))
runtime_files <- c(list.files(file.path(candidate, "library"),
  "[.](rdb|rdx|so|dll)$", recursive = TRUE, full.names = TRUE),
  list.files(dependencies, "[.](rdb|rdx|so|dll)$", recursive = TRUE, full.names = TRUE))
tasks <- make_study_plan()
plan <- list(tasks = tasks, candidate = candidate, dependencies = dependencies,
  source = frozen, hashes = tools::md5sum(c(destinations, runtime_files)),
  R_version = R.version, rng_kind = RNGkind(),
  gsl_version = system2("gsl-config", "--version", stdout = TRUE),
  concurrency = 16L,
  resource_limits = list(memory_per_cpu_mib = as.integer(Sys.getenv("IMR_MEMORY_PER_CPU_MIB")),
    max_task_cpus = as.integer(Sys.getenv("IMR_MAX_TASK_CPUS")),
    max_wall_hours = as.integer(Sys.getenv("IMR_MAX_WALL_HOURS"))), created = Sys.time())
stopifnot(all(is.finite(unlist(plan$resource_limits))), all(unlist(plan$resource_limits) > 0))
saveRDS(plan, file.path(out, "study.rds"))
writeLines(unname(tools::md5sum(file.path(out, "study.rds"))), file.path(out, "MANIFEST.md5"))
write.csv(tasks, file.path(out, "tasks.csv"), row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Frozen", nrow(tasks), "tasks:", sum(tasks$kind == "audit"), "audits and",
    sum(tasks$kind != "audit"), "scientific tasks.\n")
