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
configurations <- list(refit = list(cv_method = "refit"),
  reweight_all_draws = list(cv_method = "reweight", model_set = "all_draws"),
  reweight_top_unique = list(cv_method = "reweight", model_set = "top_unique"))
for (type in names(outcomes)) {
  for (method in c("imr", "bms")) {
    fit <- imr(list(assay = platform), outcomes[[type]], outcome_type = type,
               model_variant = method, draws = 5, burnin = 2,
               min_subgroup_size = 0, seed = 3)
    for (configuration in names(configurations)) {
      set.seed(912)
      rng <- .Random.seed
      connections <- showConnections(all = TRUE)
      arguments <- c(list(object = fit, k = 2, rounds = 1), configurations[[configuration]])
      reference <- do.call(cv_imr, arguments)
      result <- do.call(cv_imr, c(arguments, list(workers = 2L)))
      stopifnot(identical(result, reference), identical(.Random.seed, rng),
                identical(showConnections(all = TRUE), connections))
      cat(type, method, configuration, "EXACT\n")
    }
  }
}

verify_valgrind_worker_logs(36L)
