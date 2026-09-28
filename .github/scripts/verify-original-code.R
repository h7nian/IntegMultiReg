# Run from the package root after installing the candidate package.
# Rscript .github/scripts/verify-original-code.R ../Ref/biom12587-sup-0002-suppdata_code.zip OUTPUT
# The reference functions come from the archived release, never the candidate.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, file.exists(args[1]))
archive <- normalizePath(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[2])
work <- file.path(out, "reference-source")
dir.create(work, showWarnings = FALSE)
files <- c("PredictionCV.c", "SampleGam.c", "Logposterior.c", "LogLikNonLocal.c",
           "utils.c", "utils.h", "myheader.h")
utils::unzip(archive, files = paste0("CcodeBiometrics/", files), exdir = work)
work <- file.path(work, "CcodeBiometrics")
source_hashes <- tools::md5sum(file.path(work, files))
# malloc.h is unavailable on macOS. This header-only portability change is
# recorded; it does not change the original allocation or numerical code.
utility <- file.path(work, "utils.c")
lines <- readLines(utility, warn = FALSE)
writeLines(sub("#include <malloc.h>", "#include <stdlib.h>", lines, fixed = TRUE), utility)
file.copy(".github/scripts/original-code-probe.c", file.path(work, "probe.c"), overwrite = TRUE)
old <- setwd(work)
gsl_flags <- system2("gsl-config", "--cflags", stdout = TRUE)
gsl_libs <- system2("gsl-config", "--libs", stdout = TRUE)
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "SHLIB", "-o", paste0("original", .Platform$dynlib.ext),
    "probe.c", files[grepl("[.]c$", files)]),
  env = c(paste0("PKG_CPPFLAGS=", shQuote(gsl_flags)),
          paste0("PKG_LIBS=", shQuote(gsl_libs)), "PKG_CFLAGS=-std=gnu11"),
  stdout = file.path(out, "compile.log"), stderr = file.path(out, "compile.log"))
setwd(old)
if (status != 0L) stop("Original reference compilation failed; see compile.log")
library(IntegMultiReg)
dll <- dyn.load(file.path(work, paste0("original", .Platform$dynlib.ext)))
id <- seq_len(40)
fit <- imr(list(assay = data.frame(id, a = sin(id), b = cos(id / 3))),
  data.frame(id, y = sin(id / 2) + cos(id / 5)), method = "bms", outcome_type = "continuous",
  draws = 6, burnin = 2, min_subgroup_size = 0, seed = 17,
  residual_prior = c(shape = .37, rate = .21))
states <- rbind(c(0L, 0L), c(1L, 0L), c(1L, 0L), c(1L, 1L), c(0L, 0L), c(1L, 0L))
fit$posterior$selection_draws <- lapply(seq_len(nrow(states)), function(i) list(matrix(states[i, ], 1)))
result <- cv_imr(fit, k = 2, rounds = 2, cv_method = "importance", ridge = 0,
                 df_method = "legacy_integer", score_method = "legacy")
comparisons <- list()
for (round in 1:2) for (fold in 1:2) {
  rows <- result$control$folds[result$control$folds$round == round, ]
  rows <- rows[order(rows$row_order), ]
  test <- match(rows$id[rows$fold == fold], id)
  train <- match(rows$id[rows$fold != fold], id)
  expected <- .Call("original_cv", fit$preprocessing$features[[1]][[1]],
    fit$posterior$latent_response_mean[[1]], states, as.integer(train - 1L),
    as.integer(test - 1L), PACKAGE = dll[["name"]])
  actual <- result$predictions[result$predictions$round == round, ]
  actual <- actual$prediction[match(id[test], actual$id)]
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))
  comparisons[[length(comparisons) + 1L]] <- data.frame(round, fold,
    max_abs_difference = max(abs(actual - expected)))
}
write.csv(do.call(rbind, comparisons), file.path(out, "predictions.csv"), row.names = FALSE)
writeLines(c("Reference: original pred_aftcv conditional on supplied states, latent response, and ordered folds.",
  "This is component equivalence, not a reproduction of published Table 1 or the full selection sampler.",
  "Portability patch: utils.c malloc.h include changed to stdlib.h; algorithms unchanged.",
  paste(names(source_hashes), source_hashes), capture.output(sessionInfo())),
  file.path(out, "manifest.txt"))
print(do.call(rbind, comparisons))
