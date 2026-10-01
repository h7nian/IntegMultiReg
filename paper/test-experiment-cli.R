# Run from any directory with Rscript; requires no installed candidate package.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                   value = TRUE)[1])
source(file.path(dirname(normalizePath(script)),
                 "original-experiment-helpers.R"))
stopifnot(identical(parse_experiment_args(character()), list()))
stopifnot(identical(parse_experiment_args(c("--marker-design", "random")),
  list("--marker-design" = "random")))
stopifnot(identical(
  parse_experiment_args(c("--quick", "--draws", "50", "--out-dir", "a path")),
  list("--quick" = TRUE, "--draws" = "50", "--out-dir" = "a path")))
cases <- list(
  list(c("--drawws", "50"), "Unknown argument"),
  list(c("--draws", "50", "--draws", "60"), "Duplicate argument"),
  list(c("--quick", "--quick"), "Duplicate argument"),
  list("--draws", "Missing value"),
  list(c("--draws", "--quick"), "Missing value"),
  list(c("--quick", "true"), "Unknown argument"),
  list(c("--data", ""), "Missing value"))
for (case in cases) {
  error <- tryCatch(parse_experiment_args(case[[1]]), error = identity)
  stopifnot(inherits(error, "error"),
    grepl(case[[2]], conditionMessage(error), fixed = TRUE))
}
cat("PASS: valid options and 7 rejected malformed cases\n")
