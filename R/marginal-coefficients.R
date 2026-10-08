# Coefficient-marginal updates using polynomial-exact Gaussian quadrature.
.imr_coefficients_marginal_sample <- function(groups, feature_platform, nu, draws = 30000L,
                                              burnin = 2000L, seed = 1L, model_variant = "imr",
                                              initial = "empty", keep_latent = TRUE,
                                              interaction_prior = c(shape = 2, rate = 2),
                                              thin = 1L, theta_step = .4, swap_rate = .5,
                                              max_nodes = 200000L, initial_state = NULL, variance_step = NULL) {
  set.seed(seed)
  count <- length(groups)
  blocks <- .imr_marginal_blocks(groups, feature_platform, nu, model_variant,
    coefficient_columns = FALSE
  )

  selection_state <- lapply(seq_along(groups), function(s) {
    g <- groups[[s]]
    switch(initial,
      empty = rep(0L, length(g$feature_index)),
      full = rep(1L, length(g$feature_index)),
      alternating = as.integer((g$feature_index + s) %% 2 == 0)
    )
  })
  sigma2 <- vapply(groups, function(g) g$residual_prior[["rate"]] / (g$residual_prior[["shape"]] + 1), 0)
  if (is.null(variance_step)) variance_step <- vapply(groups, function(g) min(.4, 1.5 / sqrt(g$residual_prior[["shape"]] + nrow(g$X) / 2)), 0)
  if (!is.null(initial_state)) {
    selection_state <- lapply(seq_along(groups), function(s) {
      as.integer(initial_state$coefficients[[s]][groups[[s]]$n_forced + seq_along(groups[[s]]$feature_index)] != 0)
    })
    sigma2 <- initial_state$variance
    for (l in seq_along(blocks)) blocks[[l]]$theta <- initial_state$interaction[[l]]
  }
  if (length(variance_step) == 1L) variance_step <- rep(variance_step, count)
  z <- lapply(groups, function(g) {
    if (g$outcome_type == "binary") ifelse(g$y == 1, .5, -.5) else as.double(g$y)
  })
  model_cache <- lapply(groups, function(g) new.env(parent = emptyenv()))
  quadrature_cache <- new.env(parent = emptyenv())
  gaussian_rule <- function(d) {
    key <- as.character(d)
    if (exists(key, quadrature_cache, inherits = FALSE)) {
      return(get(key, quadrature_cache))
    }
    orders <- seq.int(d + 1L, 2L)
    if (sum(log(orders)) > log(max_nodes)) stop("Exact reference cubature exceeds its node budget; no approximation or model truncation was used.")
    rules <- lapply(orders, function(q) {
      jacobi <- matrix(0, q, q)
      jacobi[cbind(1:(q - 1L), 2:q)] <- sqrt(1:(q - 1L))
      jacobi <- jacobi + t(jacobi)
      result <- eigen(jacobi, symmetric = TRUE)
      list(nodes = result$values, weights = result$vectors[1, ]^2)
    })
    index <- as.matrix(expand.grid(lapply(orders, seq_len)))
    nodes <- matrix(0, nrow(index), d)
    log_weight <- numeric(nrow(index))
    for (j in seq_len(d)) {
      nodes[, j] <- rules[[j]]$nodes[index[, j]]
      log_weight <- log_weight + log(rules[[j]]$weights[index[, j]])
    }
    result <- list(nodes = nodes, log_weight = log_weight)
    assign(key, result, quadrature_cache)
    result
  }
  statistics <- function(s, selected = selection_state[[s]], response = z[[s]]) {
    g <- groups[[s]]
    key <- paste0(selected, collapse = "")
    cache <- model_cache[[s]]
    if (exists(key, cache, inherits = FALSE)) {
      result <- get(key, cache)
      if (identical(result$response, response)) {
        return(result)
      }
    } else {
      active <- c(seq_len(g$n_forced), g$n_forced + which(selected != 0))
      X <- g$X[, active, drop = FALSE]
      tau <- g$prior_scale[active]
      d <- length(active)
      factor <- chol(crossprod(X) + diag(1 / tau, nrow = d))
      covariance <- chol2inv(factor)
      normal_factor <- chol(covariance)
      rule <- gaussian_rule(d)
      a <- g$residual_prior[["shape"]]
      b <- g$residual_prior[["rate"]]
      result <- list(
        active = active, X = X, tau = tau, covariance = covariance,
        normal_factor = normal_factor, nodes = rule$nodes %*% normal_factor,
        log_weight = rule$log_weight, H = a + nrow(X) / 2 + d,
        log_D = a * log(b) - lgamma(a) - nrow(X) / 2 * log(2 * pi) -
          1.5 * sum(log(tau)) - sum(log(diag(factor))),
        latent_precision = if (g$outcome_type == "continuous") {
          NULL
        } else {
          chol2inv(chol(diag(nrow(X)) + tcrossprod(sweep(X, 2, sqrt(tau), "*"))))
        }
      )
    }
    result$response <- response
    result$mean <- drop(result$covariance %*% crossprod(result$X, response))
    result$B <- g$residual_prior[["rate"]] + .5 *
      (sum((response - drop(result$X %*% result$mean))^2) + sum(result$mean^2 / result$tau))
    assign(key, result, cache)
    result
  }
  log_kernel <- function(model, variance) {
    values <- sweep(sqrt(variance) * model$nodes, 2, model$mean, "+")
    log_moment <- .imr_log_sum_exp(model$log_weight + 2 * rowSums(log(abs(values))))
    value <- model$log_D - (model$H + 1) * log(variance) - model$B / variance + log_moment
    if (!is.finite(value)) stop("Non-finite exact coefficient-integrated density.")
    value
  }
  restore_coefficients <- function(model, variance) {
    d <- length(model$active)
    covariance <- variance * model$covariance
    scales <- sqrt(model$mean^2 + diag(covariance))
    scaled_mean_norm <- sqrt(sum((model$mean / scales)^2))
    scaled_covariance <- covariance / outer(scales, scales)
    largest <- eigen(scaled_covariance, symmetric = TRUE, only.values = TRUE)$values[1] *
      (1 + 64 * .Machine$double.eps * d)
    # q = Normal(mean, 2 C). Bound the product by the squared Euclidean norm,
    # then maximize its radial polynomial times exp(-r^2/4).
    a <- .25
    b <- a * scaled_mean_norm / sqrt(largest)
    radius <- 2 * d / (b + sqrt(b^2 + 4 * a * d))
    log_bound <- 2 * d * log(scaled_mean_norm + sqrt(largest) * radius) -
      d * log(d) - a * radius^2 + 64 * .Machine$double.eps * d
    repeat {
      normal <- matrix(stats::rnorm(16L * d), 16L, d)
      proposed <- sweep(sqrt(2 * variance) * (normal %*% model$normal_factor), 2, model$mean, "+")
      log_acceptance <- 2 * rowSums(log(abs(sweep(proposed, 2, scales, "/")))) -
        .5 * rowSums(normal^2) - log_bound
      if (any(log_acceptance > 1e-10) || anyNA(log_acceptance)) stop("Invalid coefficient rejection envelope.")
      accepted <- which(log(stats::runif(16)) < pmin(0, log_acceptance))
      if (length(accepted)) {
        return(proposed[accepted[1], ])
      }
    }
  }
  beta_draws <- lapply(groups, function(g) matrix(0, draws, ncol(g$X)))
  sigma2_draws <- matrix(0, draws, count)
  theta_draws <- lapply(blocks, function(b) {
    matrix(
      0, draws,
      if (model_variant == "imr") choose(length(b$members), 2) else 0L
    )
  })
  latent_draws <- if (keep_latent) {
    lapply(groups, function(g) {
      if (g$outcome_type == "continuous") NULL else matrix(0, draws, nrow(g$X))
    })
  } else {
    NULL
  }
  log_density <- numeric(draws)
  acceptance <- c(
    swap_proposals = 0, swap_accepts = 0, interaction_proposals = 0,
    interaction_accepts = 0, variance_proposals = 0, variance_accepts = 0, latent_proposals = 0, latent_accepts = 0
  )
  for (iteration in seq_len(burnin + draws * thin)) {
    for (s in seq_along(groups)) {
      g <- groups[[s]]
      augmented <- if (g$outcome_type == "binary") seq_len(nrow(g$X)) else if (g$outcome_type == "right.censored") which(g$status == 0) else integer()
      for (i in augmented) {
        model <- statistics(s)
        precision <- model$latent_precision
        location <- -sum(precision[i, -i] * z[[s]][-i]) / precision[i, i]
        scale <- sqrt(sigma2[s] / precision[i, i])
        sign <- if (g$outcome_type == "binary" && g$y[i] == 0) -1 else 1
        lower <- if (g$outcome_type == "binary") 0 else g$y[i]
        log_tail <- stats::pnorm((lower - sign * location) / scale, lower.tail = FALSE, log.p = TRUE)
        value <- sign * max(lower, sign * location + scale * stats::qnorm(log(stats::runif(1)) + log_tail,
          lower.tail = FALSE, log.p = TRUE
        ))
        trial <- z[[s]]
        trial[i] <- value
        candidate <- statistics(s, response = trial)
        ratio <- log_kernel(candidate, sigma2[s]) - log_kernel(model, sigma2[s]) +
          stats::dnorm(z[[s]][i], location, scale, log = TRUE) - stats::dnorm(value, location, scale, log = TRUE)
        acceptance[7] <- acceptance[7] + 1
        if (log(stats::runif(1)) < min(0, ratio)) {
          z[[s]] <- trial
          acceptance[8] <- acceptance[8] + 1
        }
      }
    }
    for (l in seq_along(blocks)) {
      block <- blocks[[l]]
      for (f in seq_along(block$features)) {
        log_bf <- vapply(seq_along(block$members), function(i) {
          s <- block$members[i]
          j <- block$columns[i, f]
          excluded <- included <- selection_state[[s]]
          excluded[j] <- 0L
          included[j] <- 1L
          log_kernel(statistics(s, included), sigma2[s]) - log_kernel(statistics(s, excluded), sigma2[s])
        }, 0)
        weights <- .imr_marginal_energy(block, nu[l]) + drop(block$patterns %*% log_bf)
        pattern <- block$patterns[sample.int(nrow(block$patterns), 1, prob = exp(weights - max(weights))), ]
        for (i in seq_along(block$members)) selection_state[[block$members[i]]][block$columns[i, f]] <- pattern[i]
      }
      if (length(block$features) >= 2L) {
        for (attempt in seq_len(ceiling(swap_rate * length(block$features)))) {
          pair <- sample.int(length(block$features), 2)
          proposed <- list()
          ratio <- 0
          for (i in seq_along(block$members)) {
            s <- block$members[i]
            columns <- block$columns[i, pair]
            if (selection_state[[s]][columns[1]] == selection_state[[s]][columns[2]]) next
            trial <- selection_state[[s]]
            trial[columns] <- rev(trial[columns])
            ratio <- ratio + log_kernel(statistics(s, trial), sigma2[s]) - log_kernel(statistics(s), sigma2[s])
            proposed[[length(proposed) + 1L]] <- list(s = s, selection = trial)
          }
          if (length(proposed)) {
            acceptance[1] <- acceptance[1] + 1
            if (log(stats::runif(1)) < min(0, ratio)) {
              for (move in proposed) selection_state[[move$s]] <- move$selection
              acceptance[2] <- acceptance[2] + 1
            }
          }
        }
      }
      selected <- block$columns
      for (i in seq_along(block$members)) selected[i, ] <- selection_state[[block$members[i]]][block$columns[i, ]]
      update <- .imr_update_marginal_interaction(
        block, selected, nu[l],
        interaction_prior, theta_step, model_variant == "imr"
      )
      block <- update$block
      acceptance[3:4] <- acceptance[3:4] + update$counts
      blocks[[l]] <- block
    }
    for (s in seq_along(groups)) {
      model <- statistics(s)
      proposal <- exp(log(sigma2[s]) + stats::rnorm(1, 0, variance_step[s]))
      acceptance[5] <- acceptance[5] + 1
      if (!is.finite(proposal) || proposal <= 0) next
      ratio <- log_kernel(model, proposal) - log_kernel(model, sigma2[s]) + log(proposal) - log(sigma2[s])
      if (log(stats::runif(1)) < min(0, ratio)) {
        sigma2[s] <- proposal
        acceptance[6] <- acceptance[6] + 1
      }
    }
    # Conditional coefficient recovery is ancillary to the marginal chain.
    # Running it on every sweep preserves the thinning/storage replay property.
    recovered_coefficients <- lapply(seq_along(groups), function(s) {
      restore_coefficients(statistics(s), sigma2[s])
    })
    if (iteration > burnin && (iteration - burnin) %% thin == 0) {
      row <- (iteration - burnin) %/% thin
      score <- 0
      for (s in seq_along(groups)) {
        model <- statistics(s)
        beta_draws[[s]][row, model$active] <- recovered_coefficients[[s]]
        sigma2_draws[row, s] <- sigma2[s]
        if (!is.null(latent_draws[[s]])) latent_draws[[s]][row, ] <- z[[s]]
        score <- score + log_kernel(model, sigma2[s])
      }
      for (l in seq_along(blocks)) {
        block <- blocks[[l]]
        patterns <- block$columns
        for (i in seq_along(block$members)) patterns[i, ] <- selection_state[[block$members[i]]][block$columns[i, ]]
        score <- score + sum(colSums(patterns) * nu[l]) + sum(patterns * (block$theta %*% patterns)) -
          length(block$features) * .imr_log_sum_exp(.imr_marginal_energy(block, nu[l]))
        if (model_variant == "imr") {
          values <- block$theta[upper.tri(block$theta)]
          theta_draws[[l]][row, ] <- values
          score <- score + sum(stats::dgamma(values, interaction_prior[1], rate = interaction_prior[2], log = TRUE))
        }
      }
      log_density[row] <- score
    }
  }
  names(beta_draws) <- names(groups)
  list(
    coefficients = beta_draws, variance = sigma2_draws, interaction = theta_draws,
    latent = latent_draws, log_density = log_density, acceptance = acceptance
  )
}
