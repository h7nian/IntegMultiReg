.imr_feature_columns <- function(model, subgroup, platform) {
  available <- model$subgroup_platforms[[subgroup]]
  preceding <- available[seq_len(match(platform, available) - 1L)]
  offset <- 1L + length(model$covariate_names) + sum(lengths(model$feature_names[preceding]))
  offset + seq_along(model$feature_names[[platform]])
}

# The compiled kernel receives subgroup designs and explicit column mappings.
# This keeps preprocessing independent of the inference engine and avoids a
# preliminary fit merely to obtain standardized inputs.
.imr_joint_design <- function(model, prep, subgroup) {
  platforms <- model$subgroup_platforms[[subgroup]]
  columns <- c(
    list(rep(1, model$sample_sizes[subgroup]), prep$covariates[[subgroup]]),
    prep$features[[subgroup]][platforms]
  )
  design <- do.call(cbind, columns)
  terms <- c(
    "(Intercept)", paste0("clinical:", model$covariate_names),
    unlist(lapply(platforms, function(l) paste0(model$platform_names[l], ":", model$feature_names[[l]])))
  )
  if (!length(model$covariate_names)) terms <- terms[terms != "clinical:"]
  if (anyDuplicated(terms)) .imr_abort("Clinical and platform coefficient labels collide; rename the conflicting platform or variable.")
  colnames(design) <- terms
  storage.mode(design) <- "double"
  design
}

.imr_joint_specification <- function(model, prep, control) {
  forced <- 1L + length(model$covariate_names)
  groups <- lapply(seq_along(model$subgroup_names), function(g) {
    design <- .imr_joint_design(model, prep, g)
    list(
      design = design, response = as.double(prep$response[[g]][, 1L]),
      status = if (control$outcome_type == "right.censored") as.integer(prep$response[[g]][, 2L]) else integer(),
      prior_scale = c(
        rep(control$priors$forced_scale, forced),
        rep(control$priors$molecular_scale, ncol(design) - forced)
      ),
      forced = forced,
      outcome = as.integer(match(control$outcome_type, c("continuous", "binary", "right.censored")) - 1L),
      residual_prior = as.double(control$priors$residual)
    )
  })
  names(groups) <- model$subgroup_names
  platforms <- lapply(seq_len(model$n_platforms), function(l) {
    members <- model$platform_subgroups[[l]]
    columns <- matrix(0L, nrow = length(members), ncol = length(model$feature_names[[l]]))
    for (i in seq_along(members)) {
      columns[i, ] <- .imr_feature_columns(model, members[i], l) - 1L
    }
    storage.mode(columns) <- "integer"
    list(groups = as.integer(members - 1L), columns = columns, nu = as.double(control$priors$nu[l]))
  })
  names(platforms) <- model$platform_names
  list(groups = groups, platforms = platforms)
}

.imr_chain_start <- function(chain, spec, model, control) {
  mcmc <- control$mcmc
  if (is.list(mcmc$initial)) {
    initial <- mcmc$initial[[chain]]
    if (!is.list(initial) || anyDuplicated(names(initial)) ||
      !setequal(names(initial), c("coefficients", "variance", "interaction"))) {
      .imr_abort("Each initial state needs `coefficients`, `variance` and `interaction`.")
    }
    if (!identical(names(initial$coefficients), model$subgroup_names) ||
      !identical(names(initial$variance), model$subgroup_names) ||
      !identical(names(initial$interaction), model$platform_names)) {
      .imr_abort("Initial state components must be named in fitted subgroup/platform order.")
    }
    for (g in seq_along(spec$groups)) {
      if (!identical(names(initial$coefficients[[g]]), colnames(spec$groups[[g]]$design))) {
        .imr_abort("Initial coefficient names must match the subgroup design columns.")
      }
      .imr_check_numeric_vector(initial$coefficients[[g]], "initial coefficients")
      initial$coefficients[[g]] <- stats::setNames(as.double(initial$coefficients[[g]]), colnames(spec$groups[[g]]$design))
    }
    .imr_check_numeric_vector(initial$variance, "initial variance", positive = TRUE)
    initial$variance <- stats::setNames(as.double(initial$variance), model$subgroup_names)
    for (l in seq_along(initial$interaction)) {
      x <- initial$interaction[[l]]
      labels <- model$subgroup_names[model$platform_subgroups[[l]]]
      if (!length(labels)) labels <- NULL
      if (!is.matrix(x) || !identical(dim(x), rep(length(labels), 2L)) ||
        !identical(dimnames(x), list(labels, labels))) {
        .imr_abort("Initial interaction matrices must name their fitted subgroup rows and columns in order.")
      }
      .imr_check_numeric_vector(x, "initial interaction")
      storage.mode(x) <- "double"
      initial$interaction[[l]] <- x
    }
    return(initial)
  }
  mode <- if (mcmc$initial == "dispersed") c("empty", "full", "alternating", "random")[(chain - 1L) %% 4L + 1L] else mcmc$initial
  coefficients <- lapply(seq_along(spec$groups), function(g) {
    dat <- spec$groups[[g]]
    molecular <- ncol(dat$design) - dat$forced
    selected <- switch(mode,
      empty = rep(0, molecular),
      full = rep(1, molecular),
      alternating = as.integer((seq_len(molecular) + g) %% 2 == 0),
      random = stats::rbinom(molecular, 1, .5)
    )
    sign <- if (chain %% 2) 1 else -1
    c(rep(.1 * sign, dat$forced), .1 * sign * selected)
  })
  variance <- vapply(spec$groups, function(g) {
    if (g$outcome == 1L) {
      return(1)
    }
    observed_variance <- if (length(g$response) > 1L) stats::var(g$response) else 0
    max(observed_variance, g$residual_prior[2L] / (g$residual_prior[1L] + 1), .Machine$double.eps)
  }, 0) * if (mcmc$initial == "dispersed") c(.5, 1, 2, 1.5)[(chain - 1L) %% 4L + 1L] else 1
  probability <- c(.1, .5, .9, .75)[(chain - 1L) %% 4L + 1L]
  theta <- if (control$model_variant == "bms") {
    0
  } else {
    stats::qgamma(probability,
      shape = control$priors$interaction[1L], rate = control$priors$interaction[2L]
    )
  }
  interaction <- lapply(spec$platforms, function(p) {
    x <- matrix(theta, length(p$groups), length(p$groups))
    diag(x) <- 0
    x
  })
  names(coefficients) <- model$subgroup_names
  for (g in seq_along(coefficients)) names(coefficients[[g]]) <- colnames(spec$groups[[g]]$design)
  names(interaction) <- model$platform_names
  for (l in seq_along(interaction)) {
    labels <- model$subgroup_names[model$platform_subgroups[[l]]]
    dimnames(interaction[[l]]) <- list(labels, labels)
  }
  list(coefficients = coefficients, variance = stats::setNames(as.double(variance), model$subgroup_names), interaction = interaction)
}

.imr_joint_chain <- function(task, spec, model, control, verbose = FALSE) {
  saved_rng <- .imr_save_rng()
  on.exit(.imr_restore_rng(saved_rng), add = TRUE)
  set.seed(task$initial_seed)
  initial <- .imr_chain_start(task$chain, spec, model, control)
  set.seed(task$seed)
  settings <- list(
    draws = control$mcmc$draws, burnin = control$mcmc$burnin,
    thin = control$mcmc$thin, sharing = as.integer(control$model_variant == "imr"),
    keep_latent = as.integer(control$mcmc$keep_latent), verbose = as.integer(verbose),
    theta_step = as.double(control$mcmc$theta_step), swap_rate = as.double(control$mcmc$swap_rate),
    interaction_prior = as.double(control$priors$interaction)
  )
  value <- .Call("imr_joint_sample", spec$groups, spec$platforms, settings, initial,
    PACKAGE = "IntegMultiReg"
  )
  value$initial <- initial
  value
}

.imr_combine_chains <- function(chains, model, prep, spec, mcmc) {
  combine <- function(matrices, variables) {
    out <- array(unlist(matrices, use.names = FALSE),
      dim = c(mcmc$draws, length(variables), mcmc$chains)
    )
    out <- aperm(out, c(1L, 3L, 2L))
    dimnames(out) <- list(iteration = NULL, chain = paste0("chain", seq_len(mcmc$chains)), parameter = variables)
    out
  }
  coefficients <- lapply(seq_along(spec$groups), function(g) {
    combine(lapply(chains, function(x) x$coefficients[[g]]), colnames(spec$groups[[g]]$design))
  })
  names(coefficients) <- model$subgroup_names
  interaction <- lapply(seq_len(model$n_platforms), function(l) {
    groups <- model$subgroup_names[model$platform_subgroups[[l]]]
    pairs <- if (length(groups) < 2L || dim(chains[[1L]]$interaction[[l]])[2L] == 0L) {
      character()
    } else {
      unlist(lapply(seq.int(2L, length(groups)), function(j) paste(groups[seq_len(j - 1L)], groups[j], sep = ":")))
    }
    combine(lapply(chains, function(x) x$interaction[[l]]), pairs)
  })
  names(interaction) <- model$platform_names
  latent <- if (mcmc$keep_latent) {
    lapply(seq_along(model$subgroup_names), function(g) {
      if (is.null(chains[[1L]]$latent[[g]])) {
        return(NULL)
      }
      combine(lapply(chains, function(x) x$latent[[g]]), as.character(prep$subject_ids[[g]]))
    })
  } else {
    NULL
  }
  if (!is.null(latent)) names(latent) <- model$subgroup_names
  list(
    coefficients = coefficients,
    variance = combine(lapply(chains, `[[`, "variance"), model$subgroup_names),
    interaction = interaction, latent = latent,
    log_density = do.call(cbind, lapply(chains, `[[`, "log_density"))
  )
}
