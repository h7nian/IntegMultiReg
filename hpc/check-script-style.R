# Keep the research scripts spaced the way the package sources are. Checks are
# token-based, so a comma inside a string or a semicolon inside a comment is
# never reported. Indentation is deliberately not checked: continuation lines
# are aligned by judgement rather than by a rule.
args <- commandArgs(TRUE)
root <- normalizePath(Sys.getenv("IMR_ROOT"), mustWork = TRUE)
files <- if (length(args)) normalizePath(args, mustWork = TRUE) else
  c(list.files(file.path(root, "paper"), "[.]R$", full.names = TRUE),
    list.files(file.path(root, "hpc"), "[.]R$", full.names = TRUE))
spaced <- c("LEFT_ASSIGN", "RIGHT_ASSIGN", "EQ_ASSIGN", "EQ_SUB", "EQ_FORMALS",
  "EQ", "NE", "LT", "GT", "LE", "GE", "AND", "AND2", "OR", "OR2")
heads <- c("IF", "FOR", "WHILE")
faults <- character(0)
# A single-quoted literal is only accepted when double quoting it would mean
# escaping a quote of its own.
single_quoted <- function(text)
  startsWith(text, "'") & !grepl('"', text, fixed = TRUE)
report <- function(path, line, message)
  faults <<- c(faults, sprintf("%s:%d: %s", path, line, message))
for (path in files) {
  lines <- readLines(path, warn = FALSE)
  for (line in grep("\t", lines)) report(path, line, "tab indentation")
  for (line in grep("[ \t]+$", lines)) report(path, line, "trailing whitespace")
  data <- utils::getParseData(parse(path, keep.source = TRUE))
  strings <- data[data$token == "STR_CONST", , drop = FALSE]
  for (i in which(single_quoted(strings$text)))
    report(path, strings$line1[i], "single-quoted string; the package uses `\"`")
  terminal <- data[data$terminal, , drop = FALSE]
  terminal <- terminal[order(terminal$line1, terminal$col1), , drop = FALSE]
  for (i in seq_len(nrow(terminal) - 1L)) {
    left <- terminal[i, ]
    right <- terminal[i + 1L, ]
    if (left$line2 != right$line1) next
    gap <- substr(lines[left$line2], left$col2 + 1L, right$col1 - 1L)
    # An omitted subscript, as in `x[i, , drop = FALSE]`, is written with the
    # space that separates the two commas.
    if (right$token == "','" && nzchar(gap) &&
        !left$token %in% c("','", "'['", "LBB"))
      report(path, left$line2, "space before `,`")
    else if (left$token == "';'" && !identical(gap, " "))
      report(path, left$line2, "`;` needs one following space, or its own line")
    else if (left$token == "','" && !identical(gap, " "))
      report(path, left$line2, "`,` needs exactly one following space")
    else if ((left$token %in% spaced || right$token %in% spaced) && !identical(gap, " "))
      report(path, left$line2, sprintf("`%s` needs one space on each side",
        if (left$token %in% spaced) left$text else right$text))
    else if (left$token %in% heads && right$token == "'('" && !identical(gap, " "))
      report(path, left$line2, sprintf("`%s` needs a space before its condition", left$text))
  }
}
if (length(faults)) {
  writeLines(faults)
  stop(sprintf("%d script style violations; see IntegMultiReg/R for the convention.",
    length(faults)), call. = FALSE)
}
cat(sprintf("Script style checked: %d files.\n", length(files)))
