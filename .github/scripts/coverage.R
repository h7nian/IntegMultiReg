# Run from the package root. Report executed source expressions, not branch or
# statistical coverage; preserve the native/R split and the public API audit.
out <- Sys.getenv("IMR_COVERAGE_DIR", file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "imr-coverage"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/")
probe <- normalizePath(".github/scripts/verify-predict-dispatch.R", winslash = "/")
coverage <- covr::package_coverage(type = "tests", quiet = FALSE,
  code = sprintf("source(%s); verify_predict_dispatch(%s)",
    encodeString(probe, quote = '"'), encodeString(out, quote = '"')))
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
# Only literal aliases in the source qualify for shared-body attribution.
# Each name must still have its own observed generic dispatch; a hit on a
# different alias alone is insufficient evidence of execution.
links <- list()
for (file in list.files("R", full.names = TRUE, pattern = "[.]R$")) {
  for (node in parse(file, keep.source = FALSE)) {
    if (is.call(node) && length(node) == 3L &&
      identical(node[[1L]], as.name("<-")) &&
      is.name(node[[2L]]) && is.name(node[[3L]])) {
      alias_names <- vapply(as.list(node)[2:3], as.character, "")
      if (all(alias_names %in% api$function_name)) links[[length(links) + 1L]] <- alias_names
    }
  }
}
dispatch <- read.csv(file.path(out, "alias-dispatch.csv"), stringsAsFactors = FALSE)
api$evidence <- "source expressions"
for (i in seq_len(nrow(api))) {
  group <- api$function_name[i]
  repeat {
    expanded <- unique(c(group, unlist(links[vapply(links, function(pair) any(pair %in% group), TRUE)])))
    if (setequal(group, expanded)) break
    group <- expanded
  }
  if (length(group) == 1L) next
  shared_body <- any(lines$functions %in% group & lines$value > 0, na.rm = TRUE)
  own_dispatch <- dispatch$dispatch_hits[match(api$function_name[i], dispatch$function_name)]
  api$executed[i] <- shared_body && !is.na(own_dispatch) && own_dispatch > 0L
  api$evidence[i] <- "shared source body and observed S3 dispatch"
}
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
