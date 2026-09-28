# Focused child-process gate, supplementary to the full parent-process suite.
# No private data, custom CV folds, or production instrumentation branches.
library(IntegMultiReg)
options(warn = 2)
scripts <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts")
source(file.path(scripts, "configure-valgrind-workers.R"))
source(file.path(scripts, "verify-valgrind-logs.R"))

platform <- data.frame(id = 1:12, marker = sin(1:12))
outcomes <- list(binary = data.frame(id = 1:12, y = rep(0:1, 6)),
                 continuous = data.frame(id = 1:12, y = cos(1:12)),
                 right.censored = data.frame(id = 1:12, time = 2:13,
                                              status = rep(0:1, 6)))
for (type in names(outcomes)) {
  for (method in c("imr", "bms")) {
    fit <- imr(list(assay = platform), outcomes[[type]], outcome_type = type,
               method = method, draws = 5, burnin = 2,
               min_subgroup_size = 0, seed = 3)
    for (mode in c("legacy", "refit", "importance")) {
      set.seed(912)
      rng <- .Random.seed
      connections <- showConnections(all = TRUE)
      reference <- cv_imr(fit, k = 2, rounds = 1, cv_method = mode)
      result <- cv_imr(fit, k = 2, rounds = 1, cv_method = mode, workers = 2L)
      stopifnot(identical(result, reference), identical(.Random.seed, rng),
                identical(showConnections(all = TRUE), connections))
      cat(type, method, mode, "EXACT\n")
    }
  }
}

verify_valgrind_worker_logs(36L)
