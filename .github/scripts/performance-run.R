# One isolated measurement process. Use a separate output directory per repeat.
# Example: R_LIBS=/path/to/library Rscript .github/scripts/performance-run.R 
#   /absolute/output sim continuous imr 200 100
# KIRC: replace 'sim' with the absolute converted supplement .rda path.
# Run once as warm-up and five times for timings; /usr/bin/time measures peak RSS.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 6L, length(args) <= 8L)
output <- args[1L]
if (dir.exists(output)) stop("Refusing to overwrite an existing measurement")
dir.create(output, recursive = TRUE)
library(IntegMultiReg)
outcome_type <- match.arg(args[3L], c("continuous", "binary", "right.censored"))
method <- match.arg(args[4L], c("imr", "bms"))
draws <- as.integer(args[5L])
burnin <- as.integer(args[6L])
interface <- if (length(args) >= 7L) match.arg(args[7L], c("list", "formula")) else "list"
workers <- if (length(args) >= 8L) as.integer(args[8L]) else 1L
stopifnot(!is.na(draws), draws > 0L, !is.na(burnin), burnin >= 0L,
          !is.na(workers), workers >= 1L)
cat("Measurement PID:", Sys.getpid(), "\n")
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
saveRDS(list(arguments = args, library = find.package("IntegMultiReg"),
             source_sha = Sys.getenv("IMR_SOURCE_SHA", unset = NA_character_),
             threads = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                                    "VECLIB_MAXIMUM_THREADS"))),
        file.path(output, "provenance.rds"))
timings <- list()
measure <- function(name, expression) {
  gc()
  started <- proc.time()
  value <- force(expression)
  elapsed <- proc.time() - started
  timings[[name]] <<- data.frame(stage = name, user = elapsed[[1L]],
                                system = elapsed[[2L]], elapsed = elapsed[[3L]])
  write.csv(do.call(rbind, timings), file.path(output, "timings.csv"), row.names = FALSE)
  saveRDS(value, file.path(output, paste0(name, ".rds")))
  value
}
data <- measure("data", {
  if (args[2L] == "sim") {
    data("simIMR", package = "IntegMultiReg")
    list(platforms = simIMR$platforms, covariates = simIMR$covariates,
         outcome = simIMR[[paste0("outcome.", switch(outcome_type,
           right.censored = "survival", outcome_type))]])
  } else {
    stopifnot(outcome_type == "right.censored")
    environment <- new.env(parent = emptyenv())
    load(args[2L], envir = environment)
    dat <- environment$kirc_full
    stopifnot(!is.null(dat$outcome.raw), nrow(dat$outcome.raw) == 448L)
    writeLines(unname(tools::md5sum(args[2L])), file.path(output, "input-md5.txt"))
    list(platforms = dat$platforms, covariates = dat$covariates,
         outcome = dat$outcome.raw)
  }
})
fit_args <- list(x = data$platforms, outcome = data$outcome,
                 covariates = data$covariates, outcome_type = outcome_type,
                 method = method, draws = draws, burnin = burnin,
                 seed = 100, nu = c(-4, -3, -4), min_subgroup_size = 30)
if (args[2L] != "sim") {
  fit_args$molecular_prior_scale <- 0.087
  fit_args$forced_prior_scale <- 10000
  fit_args$residual_prior <- c(shape = .001, rate = .001)
  fit_args$interaction_prior <- c(shape = 40, rate = 10)
}
if (interface == "formula") {
  stopifnot(outcome_type != "right.censored")
  clinical <- data$covariates
  clinical$.response <- data$outcome[match(clinical$id, data$outcome$id), 2L]
  fit_args$x <- reformulate(setdiff(names(data$covariates), "id"), ".response")
  fit_args$outcome <- fit_args$covariates <- NULL
  fit_args$data <- clinical
  fit_args$platforms <- data$platforms
}
fit <- measure("fit", do.call(imr, fit_args))
state_keys <- unlist(lapply(fit$posterior$selection_draws, function(state) {
  paste(unlist(state, use.names = FALSE), collapse = "")
}), use.names = FALSE)
saveRDS(list(draws = length(state_keys), unique_states = length(unique(state_keys))),
        file.path(output, "state-counts.rds"))
for (repeat_index in 1:2) {
  measure(paste0("predict-", repeat_index), predict(fit, data$platforms,
    covariates = if (interface == "formula") clinical else data$covariates))
}
for (cv_method in c("legacy", "refit", "importance")) {
  cv_args <- list(object = fit, k = 5L, rounds = 2L, cv_method = cv_method)
  if ("workers" %in% names(formals(cv_imr))) cv_args$workers <- workers
  else if (workers != 1L) stop("Installed version does not implement workers")
  measure(paste0("cv-", cv_method), do.call(cv_imr, cv_args))
}
# These short conditional chains measure execution, not convergence. Keep warnings.
measure("posterior", posterior_draws(fit, draws = 20L, burnin = 20L,
                                     conditional_draws = 20L, seed = 101L))
