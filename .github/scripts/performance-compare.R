# Compare outputs of two isolated performance-run.R processes.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
read <- function(directory, stage) readRDS(file.path(directory, paste0(stage, ".rds")))
baseline <- read(args[1L], "fit")
candidate <- read(args[2L], "fit")
stopifnot(identical(baseline$posterior, candidate$posterior),
          identical(baseline$model, candidate$model),
          identical(baseline$preprocessing, candidate$preprocessing))
close <- function(x, y, absolute = 1e-10, relative = 1e-8) {
  stopifnot(identical(dim(x), dim(y)), identical(is.na(x), is.na(y)),
            identical(is.finite(x), is.finite(y)))
  keep <- is.finite(x)
  stopifnot(all(abs(x[keep] - y[keep]) <= absolute + relative * abs(x[keep])))
}
for (method in c("legacy", "refit", "importance")) {
  a <- read(args[1L], paste0("cv-", method))
  b <- read(args[2L], paste0("cv-", method))
  stopifnot(identical(names(a), names(b)), identical(a$validation, b$validation),
    identical(a$metric, b$metric),
    identical(a$predictions[setdiff(names(a$predictions), "prediction")],
              b$predictions[setdiff(names(b$predictions), "prediction")]))
  close(a$predictions$prediction, b$predictions$prediction)
  for (field in c("pooled", "fold_mean")) {
    if (a$metric == "MSE") close(a[[field]], b[[field]])
    else close(a[[field]], b[[field]], 1e-12, 0)
  }
}
for (index in 1:2) {
  a <- read(args[1L], paste0("predict-", index))
  b <- read(args[2L], paste0("predict-", index))
  stopifnot(identical(names(a), names(b)), length(a) == length(b))
  for (group in seq_along(a)) {
    close(a[[group]]$prediction, b[[group]]$prediction)
    a[[group]]$prediction <- b[[group]]$prediction <- NULL
    stopifnot(identical(a[[group]], b[[group]]))
  }
}
a <- read(args[1L], "posterior")
b <- read(args[2L], "posterior")
a$fit <- b$fit <- NULL
stopifnot(identical(a, b))
cat("Exact sampler/posterior and tolerance-bounded prediction/CV checks passed.\n")
