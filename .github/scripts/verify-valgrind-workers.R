# Focused child-process gate, supplementary to the full parent-process suite.
# No private data, custom CV folds, or production instrumentation branches.
library(IntegMultiReg)
options(warn = 2)
launcher <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts",
                      "valgrind-worker-rscript.sh")
stopifnot(identical(Sys.info()[["sysname"]], "Linux"), file.exists(launcher),
          dir.exists(Sys.getenv("IMR_VALGRIND_LOG_DIR")))
cluster_options <- get("defaultClusterOptions", envir = asNamespace("parallel"))
stopifnot(is.environment(cluster_options),
          exists("rscript", envir = cluster_options, inherits = FALSE),
          identical(cluster_options$rscript_args, character()))
cluster_options$rscript <- normalizePath(launcher)
cluster_options$outfile <- ""

platform <- data.frame(id = 1:12, marker = sin(1:12))
outcomes <- list(binary = data.frame(id = 1:12, y = rep(0:1, 6)),
                 continuous = data.frame(id = 1:12, y = cos(1:12)),
                 right.censored = data.frame(id = 1:12, time = 2:13,
                                              status = rep(0:1, 6)))
for (type in names(outcomes)) {
  for (method in c("imr", "bms")) {
    fit <- imr(list(assay = platform), outcomes[[type]], outcome_type = type,
               method = method, draws = 5, burnin = 2,
               min_subgroup_size = 0, seed = 3)
    for (mode in c("legacy", "refit", "importance")) {
      set.seed(912)
      rng <- .Random.seed
      connections <- showConnections(all = TRUE)
      reference <- cv_imr(fit, k = 2, rounds = 1, cv_method = mode)
      result <- cv_imr(fit, k = 2, rounds = 1, cv_method = mode, workers = 2L)
      stopifnot(identical(result, reference), identical(.Random.seed, rng),
                identical(showConnections(all = TRUE), connections))
      cat(type, method, mode, "EXACT\n")
    }
  }
}

# stopCluster sends DONE but does not join the OS processes. Require every
# child's final Memcheck report, so a late error cannot escape the CI gate.
for (attempt in seq_len(60L)) {
  logs <- list.files(Sys.getenv("IMR_VALGRIND_LOG_DIR"),
                     pattern = "^worker-[0-9]+[.]log$", full.names = TRUE)
  contents <- lapply(logs, readLines, warn = FALSE)
  complete <- vapply(contents, function(lines) any(grepl("ERROR SUMMARY:", lines,
                                                        fixed = TRUE)), logical(1L))
  if (length(logs) == 36L && all(complete)) break
  Sys.sleep(1)
}
stopifnot(length(logs) == 36L, all(complete))
for (lines in contents) {
  stopifnot(any(grepl("ERROR SUMMARY: 0 errors from 0 contexts", lines, fixed = TRUE)),
            any(grepl("definitely lost: 0 bytes in 0 blocks", lines, fixed = TRUE)) ||
              any(grepl("All heap blocks were freed -- no leaks are possible", lines,
                         fixed = TRUE)))
}
cat("Verified 36 instrumented worker processes: zero Memcheck errors and definite leaks.\n")
