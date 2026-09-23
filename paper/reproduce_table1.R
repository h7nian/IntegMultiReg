# Compatibility entry point for the maintained original-experiment runner.
# New usage: Rscript paper/original-experiments.R --experiment table1 --reference paper
args <- commandArgs(TRUE)
value <- function(flag, default = NULL) {
  at <- match(flag, args); if (is.na(at)) default else args[at + 1L]
}
mode <- value("--mode", "paper")
if (mode != "paper") stop("Use the manuscript replication script for reduced KIRC examples; this entry point uses original Table 1 data.")
forward <- args
for (flag in c("--mode", "--mcmc", "--prefix")) {
  at <- match(flag, forward)
  if (!is.na(at)) forward <- forward[-c(at, at + 1L)]
}
forward[forward == "--smoke"] <- "--quick"
mcmc <- value("--mcmc")
if (!is.null(mcmc)) {
  parts <- strsplit(mcmc, ",", fixed = TRUE)[[1]]
  if (length(parts) != 2L) stop("--mcmc is retained,burnin; the paper uses 350000,50000")
  forward <- c(forward, "--draws", parts[1], "--burnin", parts[2])
}
if (!is.null(value("--prefix"))) stop("Use --out-dir to identify a run; --prefix is no longer supported.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
runner <- file.path(dirname(normalizePath(script)), "original-experiments.R")
status <- system2(file.path(R.home("bin"), "Rscript"),
  c(shQuote(runner), "--experiment", "table1", vapply(forward, shQuote, "")))
quit(status = status)
