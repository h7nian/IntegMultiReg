# Run from the package root with the candidate installed.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[1])
work <- file.path(out, "probe-source")
dir.create(work, showWarnings = FALSE)
files <- c("log_lik_non_local.c", "log_posterior.c", "utils.c", "sample_gamma.c",
           "sampler_math.c")
file.copy(c(file.path("src", files), list.files("src", "[.]h$", full.names = TRUE),
            ".github/scripts/sampler-probe.c"), work, overwrite = TRUE)
old <- setwd(work)
status <- system2(file.path(R.home("bin"), "R"),
  c("CMD", "SHLIB", "-o", paste0("probe", .Platform$dynlib.ext), "sampler-probe.c", files),
  env = c(paste0("PKG_CPPFLAGS=", shQuote(system2("gsl-config", "--cflags", stdout = TRUE))),
          paste0("PKG_LIBS=", shQuote(system2("gsl-config", "--libs", stdout = TRUE)))),
  stdout = file.path(out, "compile.log"), stderr = file.path(out, "compile.log"))
setwd(old)
stopifnot(status == 0L)
library(IntegMultiReg)
dll <- dyn.load(file.path(work, paste0("probe", .Platform$dynlib.ext)))
probe <- function(name, ...) .Call(name, ..., PACKAGE = dll[["name"]])
stopifnot(isTRUE(all.equal(probe("review_precision", 0L), c(1.2, 1/3, 1/7, 1/2))),
          isTRUE(all.equal(probe("review_precision", 1L), c(1.2, 1/7, 1/2, 1/2))))
log_z <- function(t) log(1 + 2 * exp(-3) + exp(-6 + 2 * t))
for (pair in list(c(1, 2), c(.0001, .0002))) {
  actual <- probe("review_prior", pair[2], c(0L, 0L), 1L) -
    probe("review_prior", pair[1], c(0L, 0L), 1L)
  expected <- diff(dgamma(pair, shape = 40, rate = 10, log = TRUE)) - diff(log_z(pair))
  stopifnot(abs(actual - expected) < 1e-11)
}
stopifnot(abs(probe("review_prior", 2, c(1L, 1L), 1L) -
              probe("review_prior", 2, c(0L, 1L), 1L) - 1) < 1e-12)
set.seed(14)
n <- 30L
x <- data.frame(id = 1:n, a = rnorm(n), b = rnorm(n))
y <- data.frame(id = 1:n, y = rnorm(n))
states <- expand.grid(a = 0:1, b = 0:1)
results <- list()
for (method in c("legacy", "paper")) {
  fit <- imr(list(A = x), y, outcome_type = "continuous", method = "bms",
    nu = 0, forced_prior_scale = 1, molecular_prior_scale = 1,
    residual_prior = c(shape = 1, rate = 1), min_subgroup_size = 0,
    draws = 100000L, burnin = 1000L, seed = 12, sampler_method = method)
  ll <- apply(states, 1, function(mask) probe("review_loglik",
    fit$preprocessing$features[[1]][[1]][, mask == 1, drop = FALSE], y$y, 1, 1, 1, 1))
  normalize <- function(z) exp(z - max(z)) / sum(exp(z - max(z)))
  target <- normalize(ll)
  legacy_target <- normalize(ll - log(ifelse(rowSums(states) %in% c(0, 2), 1, .5)))
  sampled <- vapply(fit$posterior$selection_draws, function(s) sum(s[[1]][1, ] * c(1L, 2L)), 0)
  empirical <- tabulate(sampled + 1L, 4) / length(sampled)
  expected <- if (method == "paper") target else legacy_target
  stopifnot(max(abs(empirical - expected)) < .006)
  results[[method]] <- data.frame(method, states, target, legacy_target, empirical)
}
write.csv(do.call(rbind, results), file.path(out, "stationary-probabilities.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(do.call(rbind, results), digits = 6)
cat("Native Gamma density, prior indexing and stationary-distribution checks passed.\n")
