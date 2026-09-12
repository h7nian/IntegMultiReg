# stopCluster sends DONE without joining OS processes. Require final reports.
# R startup helpers inherit logging; matching PIDs select actual PSOCK workers.
verify_valgrind_worker_logs <- function(expected_workers = NULL) {
  for (attempt in seq_len(60L)) {
    logs <- list.files(Sys.getenv("IMR_VALGRIND_LOG_DIR"),
                       pattern = "^worker-[0-9]+-[0-9]+[.]log$", full.names = TRUE)
    logs <- logs[grepl("^worker-([0-9]+)-\\1[.]log$", basename(logs))]
    contents <- lapply(logs, readLines, warn = FALSE)
    complete <- vapply(contents, function(lines) any(grepl("ERROR SUMMARY:", lines,
                                                          fixed = TRUE)), logical(1L))
    count_matches <- length(logs) > 0L &&
      (is.null(expected_workers) || length(logs) == expected_workers)
    if (count_matches && all(complete)) break
    Sys.sleep(1)
  }
  stopifnot(count_matches, all(complete))
  for (lines in contents) {
    stopifnot(any(grepl("ERROR SUMMARY: 0 errors from 0 contexts", lines, fixed = TRUE)),
              any(grepl("definitely lost: 0 bytes in 0 blocks", lines, fixed = TRUE)) ||
                any(grepl("All heap blocks were freed -- no leaks are possible", lines,
                           fixed = TRUE)))
  }
  cat("Verified", length(logs),
       "instrumented worker processes: zero Memcheck errors and definite leaks.\n")
  invisible(length(logs))
}
