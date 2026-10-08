# Compare fixed-seed public Laplace outputs across object-schema changes.
# Only renamed labels/fields are normalized; numerical equality is exact.
args <- commandArgs(TRUE)
if (length(args) >= 3L) .libPaths(c(args[3L], .libPaths()))
library(IntegMultiReg)
root <- args[1L]; phase <- args[2L]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
id <- 1:60
platforms <- list(A = data.frame(id = id[1:40], marker = sin(id[1:40])),
  B = data.frame(id = id[21:60], marker = cos(id[21:60] / 2)))
covariates <- data.frame(id = id, age = cos(id / 5))
reports <- list()
for (outcome in c('continuous', 'binary', 'right.censored')) for (variant in c('imr', 'bms')) {
  y <- switch(outcome,
    continuous = data.frame(id = id, y = 1 + sin(id) + covariates$age / 3),
    binary = data.frame(id = id, y = rep(0:1, 30)),
    right.censored = data.frame(id = id, time = exp(1 + sin(id)), status = as.integer(id %% 3 != 0)))
  arguments <- list(x = platforms, outcome = y, covariates = covariates,
    outcome_type = outcome, min_subgroup_size = 0)
  if (phase == 'baseline') {
    fit_formals <- names(formals(getS3method("imr", "list")))
    arguments[[if ("model_variant" %in% fit_formals) "model_variant" else "method"]] <- variant
    arguments$forced_prior_scale <- 1
    arguments$molecular_prior_scale <- .5
    arguments$interaction_prior <- c(shape = 2, rate = 2)
    arguments$draws <- 80L; arguments$burnin <- 40L; arguments$seed <- 37L
    if ("sampler_method" %in% fit_formals) arguments$sampler_method <- "paper"
  } else {
    arguments$model_variant <- variant
    arguments$marginalize <- 'coefficients_and_variance'
    arguments$priors <- imr_priors(forced_scale = 1, molecular_scale = .5, interaction = c(shape = 2, rate = 2))
    arguments$mcmc <- imr_mcmc(draws = 80, burnin = 40, chains = 1,
      seed = 37, initial = NULL, diagnostics = FALSE)
  }
  fit <- do.call(imr, arguments)
  prediction <- if (phase == 'baseline') predict(fit, platforms, covariates = covariates) else {
    predict(fit, platforms, covariates = covariates,
      type = if (outcome == 'right.censored') 'link' else 'response')
  }
  if (phase == 'baseline') {
    result <- list(probability = lapply(fit$posterior$inclusion_probabilities, unname),
      selection = unlist(fit$posterior$selection_draws, use.names = FALSE),
      interaction = lapply(fit$posterior$interaction_draws, unname),
      log_density = unname(tail(fit$posterior$log_posterior, 80)),
      latent_mean = unname(fit$posterior$latent_response_mean))
    exported <- getNamespaceExports('IntegMultiReg')
    sampler <- getExportedValue('IntegMultiReg', if ('sample_regression_posterior' %in% exported) 'sample_regression_posterior' else 'posterior_draws')
    post_args <- list(object = fit, burnin = 40L, chains = 2L, seed = 81L, latent = TRUE)
    post_args[[if ('output_draws' %in% names(formals(sampler))) 'output_draws' else 'draws']] <- 80L
    post_args[[if ('min_draws_per_model_chain' %in% names(formals(sampler))) 'min_draws_per_model_chain' else 'conditional_draws']] <- 40L
    post <- suppressWarnings(do.call(sampler, post_args))
    coefficients <- if (is.null(post$coefficients)) post$beta else post$coefficients
  } else {
    history <- IntegMultiReg:::.imr_selection_history(fit)
    result <- list(probability = lapply(inclusion_probabilities(fit), unname),
      selection = unlist(history, use.names = FALSE),
      interaction = lapply(fit$posterior$interaction, function(x) {
        if (variant == 'bms') NULL else matrix(x, nrow = 80)
      }),
      log_density = unname(as.vector(fit$posterior$log_density)),
      latent_mean = unname(fit$posterior$latent_mean))
    post <- suppressWarnings(sample_regression_posterior(fit, output_draws = 80,
      mcmc = imr_mcmc(draws = 40, burnin = 40, chains = 2, seed = 81,
        keep_latent = TRUE, initial = NULL)))
    coefficients <- post$coefficients
  }
  result$prediction <- unname(lapply(prediction, function(x) x$prediction))
  result$probability <- unname(result$probability)
  result$interaction <- unname(result$interaction)
  result$conditional <- list(coefficients = lapply(coefficients, unname),
    variance = unname(post$variance), latent = if (is.null(post$latent)) NULL else lapply(post$latent, unname),
    selection_index = if (is.null(post$selection_draw_index)) post$model_draw else post$selection_draw_index)
  name <- paste(outcome, variant, sep = '-')
  path <- file.path(root, paste0(name, '.rds'))
  if (phase == 'baseline') saveRDS(result, path)
  same <- identical(result, readRDS(path), num.eq = FALSE)
  if (!same) print(all.equal(result, readRDS(path), tolerance = 0))
  row <- data.frame(outcome, variant, phase, identical = same)
  print(row); reports[[name]] <- row
  stopifnot(same)
  write.csv(do.call(rbind, reports), file.path(root, paste0(phase, '.csv')), row.names = FALSE)
}
