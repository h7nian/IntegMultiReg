# Run from the package root. Report executed source expressions, not branch or
# statistical coverage; preserve the native/R split and the public API audit.
out <- Sys.getenv("IMR_COVERAGE_DIR", file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "imr-coverage"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
coverage <- covr::package_coverage(type = "tests", quiet = FALSE)
saveRDS(coverage, file.path(out, "coverage.rds"))
lines <- as.data.frame(coverage)
write.csv(lines, file.path(out, "expressions.csv"), row.names = FALSE)
ns <- readLines("NAMESPACE")
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
methods <- sub("^S3method\\((.*),(.*)\\)$", "\\1.\\2",
               grep("^S3method\\(", ns, value = TRUE))
api <- data.frame(kind = c(rep("export", length(exports)), rep("S3", length(methods))),
                  function_name = c(exports, methods))
api$executed <- vapply(api$function_name, function(fn) {
  # Unattributed expressions must not turn a missing API hit into NA.
  any(lines$functions == fn & lines$value > 0, na.rm = TRUE)
}, logical(1))
write.csv(api, file.path(out, "public-api.csv"), row.names = FALSE)
summary <- vapply(list(overall = lines, R = lines[grepl("^R/", lines$filename), ],
                      C = lines[grepl("^src/", lines$filename), ]),
                  function(x) 100 * mean(x$value > 0), numeric(1))
# covr merges multiple expressions on the same line for percent_coverage().
capture.output(print(coverage), file = file.path(out, "summary.txt"), type = "message")
write.csv(data.frame(scope = names(summary), expression_percent = summary),
          file.path(out, "expression-summary.csv"), row.names = FALSE)
print(coverage)
print(api)
if (any(!api$executed)) stop("An exported function or registered S3 method was not exercised.")
