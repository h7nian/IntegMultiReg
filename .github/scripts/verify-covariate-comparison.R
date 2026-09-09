source("inst/examples/compare-covariates.R")
out <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "covariate-comparison-evidence")
set.seed(12679)
rng <- .Random.seed
first <- run_covariate_comparison(file.path(out, "first"))
stopifnot(identical(.Random.seed, rng))
second <- run_covariate_comparison(file.path(out, "repeat"))
stopifnot(identical(.Random.seed, rng))
files <- list.files(file.path(out, "first"), pattern = "[.]csv$")
stopifnot(length(files) == 8L)
for (file in files) {
  a <- file.path(out, "first", file); b <- file.path(out, "repeat", file)
  stopifnot(identical(readBin(a, "raw", file.info(a)$size),
                      readBin(b, "raw", file.info(b)$size)))
}
p <- first$outer_predictions
stopifnot(nrow(p) == 300L, !anyDuplicated(p$id),
  isTRUE(all.equal(mean((p$observed - p$prediction)^2),
                   first$nested_summary$pooled_outer_mse)))
for (fold in unique(p$outer_fold)) {
  train <- first$nested_inner_folds$id[first$nested_inner_folds$outer_fold == fold]
  stopifnot(setequal(train, p$id[p$outer_fold != fold]),
            !any(train %in% p$id[p$outer_fold == fold]))
}
cat("Verified: eight identical CSVs; preserved RNG; unique held-out predictions;\n",
    "outer MSE reconstructed independently; no outer test IDs in inner folds.\n")
