# Execute every displayed CodeInput block, retain code/output mapping, and
# refresh or add CodeOutput blocks only when --refresh-output is requested.
args <- commandArgs(TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
materials <- dirname(normalizePath(script))
path <- file.path(materials, "IntegMultiReg.tex")
text <- paste(readLines(path, warn = FALSE), collapse = "\n")
pattern <- "(?s)\\\\begin\\{CodeInput\\}(.*?)\\\\end\\{CodeInput\\}"
matches <- gregexpr(pattern, text, perl = TRUE)[[1]]
blocks <- regmatches(text, list(matches))[[1]]
code <- gsub("(?m)^(R> |\\+ +)", "", sub("\\\\end\\{CodeInput\\}$", "",
  sub("^\\\\begin\\{CodeInput\\}", "", blocks)), perl = TRUE)
hit <- match("--out-dir", args)
out <- if (!is.na(hit)) args[hit + 1L] else file.path(materials, "command-audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
# Execute unchanged displayed commands in a fresh workspace: the companion
# writes its default relative outputs there, never over historical paper files.
workspace <- file.path(out, "workspace")
if (dir.exists(workspace)) stop("Audit workspace already exists; choose a new --out-dir")
dir.create(workspace)
companion <- file.path(materials, "IntegMultiReg-covariate-comparison.R")
stopifnot(file.copy(companion, workspace))
saveRDS(list(manuscript = normalizePath(path),
  hashes = tools::md5sum(c(path, script, companion,
    system.file("R", "IntegMultiReg.rdb", package = "IntegMultiReg"),
    system.file("libs", paste0("IntegMultiReg", .Platform$dynlib.ext), package = "IntegMultiReg"))),
  started = Sys.time()), file.path(out, "provenance.rds"))
setwd(workspace)
environment <- new.env(parent = globalenv())
outputs <- vector("list", length(code))
for (i in seq_along(code)) {
  cat("Executing manuscript block", i, "\n")
  writeLines(code[i], file.path(out, sprintf("block-%02d.R", i)))
  outputs[[i]] <- capture.output({
    expressions <- parse(text = code[i])
    for (expression in expressions) {
      result <- withVisible(eval(expression, environment))
      if (result$visible) print(result$value)
    }
    invisible(NULL)
  })
  writeLines(outputs[[i]], file.path(out, sprintf("block-%02d.txt", i)))
}
escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
html <- c('<!doctype html><meta charset="utf-8"><title>Manuscript command audit</title>',
  '<style>body{max-width:1000px;margin:3rem auto;font-family:system-ui}pre{white-space:pre-wrap;background:#f5f5f5;padding:1rem}</style>',
  '<h1>Executed manuscript commands</h1>')
for (i in seq_along(code)) html <- c(html, sprintf('<h2>Block %d</h2><pre>%s</pre><pre>%s</pre>',
  i, escape(code[i]), escape(paste(outputs[[i]], collapse = "\n"))))
writeLines(c(html, paste0('<pre>', escape(paste(capture.output(sessionInfo()), collapse = "\n")), '</pre>')),
  file.path(out, "code.html"))
saveRDS(list(code = code, output = outputs, session = sessionInfo()), file.path(out, "commands.rds"))
if ("--refresh-output" %in% args) {
  for (i in rev(seq_along(matches))) {
    end <- matches[i] + attr(matches, "match.length")[i] - 1L
    suffix <- substring(text, end + 1L)
    # Match across lines explicitly for the adjacent output block.
    output_match <- regexpr("(?s)^\\s*\\\\begin\\{CodeOutput\\}(.*?)\\\\end\\{CodeOutput\\}", suffix, perl = TRUE)
    if (output_match[1] == 1L) {
      replacement <- paste0("\n\\begin{CodeOutput}\n", paste(outputs[[i]], collapse = "\n"), "\n\\end{CodeOutput}")
      text <- paste0(substr(text, 1L, end), replacement,
        substring(suffix, attr(output_match, "match.length") + 1L))
    } else if (length(outputs[[i]])) {
      replacement <- paste0("\n\\begin{CodeOutput}\n", paste(outputs[[i]], collapse = "\n"),
        "\n\\end{CodeOutput}")
      text <- paste0(substr(text, 1L, end), replacement, suffix)
    }
  }
  writeLines(text, path)
}
cat("Executed", length(code), "blocks; audit saved to", out, "\n")
