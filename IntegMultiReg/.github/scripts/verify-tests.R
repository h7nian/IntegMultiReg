# Run the installed package, not a source-loaded substitute, on every CI host.
library(IntegMultiReg)
library(testthat)

cat("Commit:", Sys.getenv("GITHUB_SHA", "local"), "\n")
print(packageVersion("IntegMultiReg"))
print(sessionInfo())
results <- test_package("IntegMultiReg", reporter = "summary",
                        stop_on_failure = FALSE, stop_on_warning = FALSE)
summary <- as.data.frame(results)
evidence <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "test-evidence")
dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
write.csv(summary[setdiff(names(summary), "result")],
          file.path(evidence, "tests.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(evidence, "sessionInfo.txt"))
writeLines(Sys.getenv("GITHUB_SHA", "local"), file.path(evidence, "commit.txt"))
stopifnot(nrow(summary) > 0L, sum(summary$passed) > 0L,
          !any(summary$failed), !any(summary$error),
          !any(summary$warning), !any(summary$skipped))
cat("Verified:", sum(summary$passed), "passed; zero failures, warnings or skips.\n")
