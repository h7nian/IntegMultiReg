# Optional coefficient/predictive validation for the manuscript simulation.
# First run IntegMultiReg-replication.R, then:
# Rscript IntegMultiReg-uncertainty.R replication-output
# This longer validation is deliberately separate from fast CRAN tests.
library(IntegMultiReg)
stopifnot(packageVersion('IntegMultiReg') == '0.2.0')
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[1L] else 'replication-output'
out <- file.path(input, 'uncertainty')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
f <- readRDS(file.path(input, 'replication_results.rds'))$simulated_fit
d <- posterior_draws(f, draws = 2000, burnin = 5000, conditional_draws = 4000, seed = 609)
print(d)
print(lapply(confint(d), head))
write.csv(d$diagnostics, file.path(out, 'diagnostics.csv'), row.names = FALSE)
saveRDS(d, file.path(out, 'draws.rds'))
for (type in c('mean', 'response')) {
  p <- predict(d, f$preprocessing$input_data$platforms, covariates = f$preprocessing$input_data$covariates, type = type)
  stopifnot(sum(vapply(p, nrow, 1L)) == 300L,
    all(is.finite(unlist(lapply(p, function(x)x[, 2:4])))))
  write.csv(do.call(rbind, p), file.path(out, paste0(type, '.csv')), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
