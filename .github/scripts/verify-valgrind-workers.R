# Instrument chain workers and both CV algorithms for every outcome/variant.
library(IntegMultiReg)
scripts <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts")
source(file.path(scripts, "configure-valgrind-workers.R"))
source(file.path(scripts, "verify-valgrind-logs.R"))
platform <- data.frame(id = 1:12, marker = sin(1:12))
outcomes <- list(binary = data.frame(id = 1:12, y = rep(0:1, 6)),
  continuous = data.frame(id = 1:12, y = cos(1:12)),
  right.censored = data.frame(id = 1:12, time = 2:13, status = rep(0:1, 6)))
# Tiny budgets deliberately make PSIS unreliable; require that warning rather
# than hiding any unrelated warning from the instrumented process.
allow_psis_warning <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (!grepl("Some PSIS fold estimates are unreliable", conditionMessage(w), fixed = TRUE)) stop(w)
  invokeRestart("muffleWarning")
})
for (type in names(outcomes)) for (variant in c("imr", "bms")) {
  fit_args <- list(x = list(assay = platform), outcome = outcomes[[type]],
    outcome_type = type, model_variant = variant, min_subgroup_size = 0)
  settings <- imr_mcmc(draws = 20, burnin = 4, chains = 2, seed = 3, diagnostics = FALSE)
  fit <- do.call(imr, c(fit_args, list(mcmc = settings)))
  set.seed(912); rng <- .Random.seed; connections <- showConnections(all = TRUE)
  settings$workers <- 2L
  parallel_fit <- do.call(imr, c(fit_args, list(mcmc = settings)))
  stopifnot(identical(fit$posterior, parallel_fit$posterior), identical(.Random.seed, rng),
    identical(showConnections(all = TRUE), connections))
  for (method in c("refit", "reweight")) {
    reference <- allow_psis_warning(cv_imr(fit, k = 2, rounds = 1, cv_method = method))
    result <- allow_psis_warning(cv_imr(fit, k = 2, rounds = 1, cv_method = method, workers = 2))
    if (!identical(result, reference)) {
      saveRDS(list(serial = reference, parallel = result), file.path(
        Sys.getenv("IMR_VALGRIND_LOG_DIR"), paste(type, variant, method, "difference.rds", sep = "-")))
      print(all.equal(result, reference, tolerance = 0))
    }
    stopifnot(identical(result, reference), identical(.Random.seed, rng),
      identical(showConnections(all = TRUE), connections))
    cat(type, variant, method, "EXACT\n")
  }
}
# 3 outcomes * 2 variants * (chain fit + refit CV + reweight CV) * 2 workers.
verify_valgrind_worker_logs(36L)
