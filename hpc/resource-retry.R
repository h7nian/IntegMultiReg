# Produce one bounded retry request from terminal Slurm resource failures.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
study <- normalizePath(args[1], mustWork = TRUE)
phase <- args[2]
source(file.path(study, "source/hpc/study-plan.R"))
plan <- read_study_plan(study)
submissions <- read.delim(file.path(study, "submissions.tsv"), header = FALSE,
                          colClasses = "character", sep = "\t")
last <- submissions[nrow(submissions), ]
stopifnot(last[[2]] == phase, grepl("^[0-9]+$", last[[1]]))
output <- system2("sacct", c("-nP", "-j", last[[1]],
                           "--format=JobID%40,State%40"), stdout = TRUE)
stopifnot(is.null(attr(output, "status")))
fields <- strsplit(output, "|", fixed = TRUE)
eligible <- list()
attempts_path <- file.path(study, "resource-retries.rds")
attempts <- if (file.exists(attempts_path)) readRDS(attempts_path) else integer()
for (row in fields) {
  if (length(row) < 2L || !grepl(paste0("^", last[[1]], "_[0-9]+$"), row[1])) next
  id <- as.integer(sub(".*_", "", row[1]))
  stopifnot(id %in% as.integer(strsplit(last[[4]], ",", fixed = TRUE)[[1]]))
  state <- trimws(row[2])
  if (!state %in% c("OUT_OF_MEMORY", "TIMEOUT", "NODE_FAIL", "PREEMPTED")) next
  if (as.character(id) %in% names(attempts)) next
  stopifnot(id %in% plan$tasks$task_id)
  eligible[[length(eligible) + 1L]] <- data.frame(id, state)
}
if (!length(eligible)) stop("No untried terminal resource failures; inspect failed task logs.")
retry <- do.call(rbind, eligible)
memory_mib <- as.integer(last[[3]]) * if (any(retry$state == "OUT_OF_MEMORY")) 2L else 1L
hours <- as.integer(last[[5]]) * if (any(retry$state == "TIMEOUT")) 2L else 1L
stopifnot(memory_mib <= plan$resource_limits$max_task_cpus * plan$resource_limits$memory_per_cpu_mib,
          hours <= plan$resource_limits$max_wall_hours)
attempts[as.character(retry$id)] <- 1L
# Reserve the retry before submitting: an interrupted controller cannot duplicate it.
temporary <- paste0(attempts_path, ".tmp-", Sys.getpid())
saveRDS(attempts, temporary)
stopifnot(file.rename(temporary, attempts_path))
write.csv(retry, file.path(study, paste0("resource-retry-", last[[1]], ".csv")),
          row.names = FALSE)
cat(paste(c(paste(retry$id, collapse = ","), memory_mib, hours),
          collapse = "\n"),
    "\n", sep = "")
