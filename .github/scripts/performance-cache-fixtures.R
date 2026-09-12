# Artificial selection-state workloads, NOT posterior samples for inference.
# Keep the data, shapes and draw count fixed; vary only cache reuse patterns.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2L]))
library(IntegMultiReg)
fit <- readRDS(args[1L])
dir.create(args[2L], recursive = TRUE)
states <- fit$posterior$selection_draws
fit$posterior$selection_draws <- rep(states[1L], length(states))
validate_imr(fit)
saveRDS(fit, file.path(args[2L], "repeated.rds"))
bits <- ceiling(log2(length(states)))
stopifnot(bits <= length(states[[1L]][[1L]]), bits <= 30L)
fit$posterior$selection_draws <- lapply(seq_along(states), function(index) {
  state <- states[[1L]]
  state[[1L]][seq_len(bits)] <- as.integer(intToBits(index - 1L))[seq_len(bits)]
  state
})
validate_imr(fit)
saveRDS(fit, file.path(args[2L], "unique.rds"))
