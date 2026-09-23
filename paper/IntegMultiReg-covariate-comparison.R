# Thin replication entry point; computation lives in the installed package example.
source(system.file("examples", "compare-covariates.R", package = "IntegMultiReg"))
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  hit <- match("--out-dir", args)
  out <- if (!is.na(hit) && hit < length(args)) args[hit + 1L] else "covariate-comparison"
  run_covariate_comparison(out, quick = "--quick" %in% args)
}
