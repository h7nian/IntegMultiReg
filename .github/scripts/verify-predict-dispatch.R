# covr attributes a shared closure to one binding. Independently count actual
# generic dispatches, while executing the instrumented shared function body.
verify_predict_dispatch <- function(out = Sys.getenv("IMR_COVERAGE_DIR",
    file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "imr-coverage"))) {
  library(IntegMultiReg)
  classes <- c("imr", "imr_selection", "imr_posterior")
  originals <- lapply(classes, function(class) getS3method("predict", class))
  names(originals) <- classes
  generic_environment <- asNamespace("stats")
  on.exit(for (class in classes) {
    registerS3method("predict", class, originals[[class]], envir = generic_environment)
  }, add = TRUE)
  hits <- new.env(parent = emptyenv())
  for (class in classes) {
    wrapper <- local({
      name <- class
      original <- originals[[class]]
      hits[[name]] <- 0L
      function(...) {
        hits[[name]] <- hits[[name]] + 1L
        original(...)
      }
    })
    registerS3method("predict", class, wrapper, envir = generic_environment)
  }
  x <- list(assay = data.frame(id = 1:12, marker = sin(1:12)))
  y <- data.frame(id = 1:12, y = 1 + cos(1:12))
  settings <- imr_mcmc(draws = 16, burnin = 8, chains = 1, seed = 27, diagnostics = FALSE)
  joint <- imr(x, y, outcome_type = "continuous", min_subgroup_size = 0,
    priors = imr_priors(forced_scale = 1), mcmc = settings)
  selection <- imr(x, y, outcome_type = "continuous", min_subgroup_size = 0,
    marginalize = "coefficients_and_variance",
    priors = imr_priors(forced_scale = 1), mcmc = settings)
  conditional <- sample_regression_posterior(selection, output_draws = 12, mcmc = settings)
  objects <- list(joint, selection, conditional)
  expected <- c("posterior_draws", "ranked_model_modes", "conditional_posterior_draws")
  for (i in seq_along(objects)) {
    prediction <- stats::predict(objects[[i]], x)
    stopifnot(identical(attr(prediction, "prediction")$method, expected[[i]]),
      all(is.finite(prediction[[1]]$prediction)))
  }
  result <- data.frame(function_name = paste0("predict.", classes),
    dispatch_hits = vapply(classes, function(class) hits[[class]], 1L))
  stopifnot(all(result$dispatch_hits > 0L))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(out, "alias-dispatch.csv"), row.names = FALSE)
}
if (sys.nframe() == 0L) verify_predict_dispatch()
