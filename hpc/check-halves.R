# Descriptive within-chain drift; not a stationarity or convergence test.
# Rscript check-halves.R RESULTS OUTPUT
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1])
out <- args[2]
paths <- file.path(root, sprintf("chain-%02d-diagnostic.rds", 1:8))
rows <- list()
for (ch in 1:8) {
  stopifnot(identical(readRDS(file.path(root,
    sprintf("chain-%02d.rds.status.rds", ch)))$status, "completed"))
  f <- readRDS(paths[ch])
  for (p in seq_along(f$posterior$interaction_draws)) {
    groups <- f$model$subgroup_names[f$model$platform_subgroups[[p]]]
    if (length(groups) < 2L) next
    pairs <- do.call(cbind, lapply(2:length(groups),
                                   function(j) rbind(seq_len(j - 1), j)))
    mat <- f$posterior$interaction_draws[[p]]
    stopifnot(ncol(mat) == ncol(pairs), nrow(mat) == f$control$mcmc$draws,
              all(is.finite(mat)), nrow(mat) >= 4L)
    cut <- floor(nrow(mat) / 2)
    for (j in seq_len(ncol(mat))) {
      first <- mat[seq_len(cut), j]
      last <- mat[seq.int(cut + 1L, nrow(mat)), j]
      q1 <- quantile(first, c(.05, .5, .95), names = FALSE)
      q2 <- quantile(last, c(.05, .5, .95), names = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(chain = ch,
        platform = f$model$platform_names[p], subgroup1 = groups[pairs[1, j]],
        subgroup2 = groups[pairs[2, j]], first_n = length(first), last_n = length(last),
        first_mean = mean(first), last_mean = mean(last),
        mean_difference = mean(last) - mean(first),
        first_q05 = q1[1], last_q05 = q2[1], first_median = q1[2], last_median = q2[2],
        first_q95 = q1[3], last_q95 = q2[3])
    }
  }
}
result <- do.call(rbind, rows)
stopifnot(nrow(result) == 64L, !anyNA(result))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(out, "theta-half-drift.csv"), row.names = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                   value = TRUE))
saveRDS(list(result = result, input_hashes = tools::md5sum(paths),
  script_hash = tools::md5sum(script), session = sessionInfo(),
  interpretation = "Contiguous halves of every retained chain, no extra burn-in or thinning. Differences are descriptive; autocorrelation prevents interpreting raw draws as independent replicates. No pass/fail threshold or convergence claim."),
  file.path(out, "theta-half-drift.rds"))
print(result[order(abs(result$mean_difference), decreasing = TRUE), ][1:8, ],
      row.names = FALSE)
