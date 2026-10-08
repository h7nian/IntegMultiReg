# Variance-marginal posterior updates and exact conditional recovery.
.imr_variance_marginal_sample <- function(groups, feature_platform, nu, draws = 30000L,
                                          burnin = 2000L, seed = 1L, model_variant = "imr",
                                          initial = "empty", keep_latent = TRUE,
                                          interaction_prior = c(shape = 2, rate = 2),
                                          thin = 1L, theta_step = .4, swap_rate = .5, initial_state = NULL) {
  set.seed(seed)
  count <- length(groups)
  blocks <- .imr_marginal_blocks(groups, feature_platform, nu, model_variant,
    coefficient_columns = TRUE
  )

  beta <- lapply(seq_along(groups), function(s) {
    g <- groups[[s]]
    selected <- switch(initial,
      empty = rep(0, length(g$feature_index)),
      full = rep(1, length(g$feature_index)),
      alternating = as.integer((g$feature_index + s) %% 2 == 0)
    )
    c(rep(.1, g$n_forced), .1 * selected)
  })
  if (!is.null(initial_state)) {
    beta <- initial_state$coefficients
    for (l in seq_along(blocks)) blocks[[l]]$theta <- initial_state$interaction[[l]]
  }
  z <- lapply(groups, function(g) {
    if (g$outcome_type == "binary") ifelse(g$y == 1, .5, -.5) else as.double(g$y)
  })
  square_t <- function(location, scale, df) {
    normalizer_scale <- max(abs(location), scale * sqrt(df / (df - 2)))
    a <- location / normalizer_scale
    b <- scale / normalizer_scale
    mixture <- a^2 / (a^2 + b^2 * df / (df - 2))
    repeat {
      proposal <- if (stats::runif(1) < mixture) {
        stats::rt(1, df)
      } else {
        sign <- if (stats::runif(1) < .5) -1 else 1
        sign * sqrt(df * stats::rchisq(1, 3) / stats::rchisq(1, df - 2))
      }
      acceptance <- (a + b * proposal)^2 / (2 * (a^2 + b^2 * proposal^2))
      if (!is.finite(acceptance)) stop("Non-finite weighted-t proposal.")
      if (stats::runif(1) < acceptance) {
        value <- location + scale * proposal
        if (!is.finite(value)) stop("Non-finite coefficient draw.")
        if (value != 0) {
          return(value)
        }
      }
    }
  }
  conditional <- function(s, j, omit = j) {
    g <- groups[[s]]
    other <- beta[[s]]
    other[omit] <- 0
    residual <- z[[s]] - drop(g$X %*% other)
    x <- g$X[, j]
    tau <- g$prior_scale[j]
    precision <- sum(x^2) + 1 / tau
    location <- sum(x * residual) / precision
    shape0 <- g$residual_prior[["shape"]] + nrow(g$X) / 2 + 1.5 * sum(other != 0)
    # A sum of positive terms avoids subtracting a fitted sum of squares.
    rate1 <- g$residual_prior[["rate"]] + .5 * (sum((residual - x * location)^2) +
      sum(other^2 / g$prior_scale) + location^2 / tau)
    gain <- .5 * precision * location^2
    log_bf <- -1.5 * log(tau) - .5 * log(precision) + shape0 * log1p(gain / rate1) +
      log(shape0 * location^2 / rate1 + 1 / precision)
    list(
      location = location, scale = sqrt(rate1 / (precision * (shape0 + 1))),
      df = 2 * shape0 + 2, log_bf = log_bf
    )
  }
  draw_coefficient <- function(parameters) {
    square_t(parameters$location, parameters$scale, parameters$df)
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
  acceptance <- c(swap_proposals = 0, swap_accepts = 0, interaction_proposals = 0, interaction_accepts = 0)
  for (iteration in seq_len(burnin + draws * thin)) {
    for (s in seq_along(groups)) {
      g <- groups[[s]]
      mean <- drop(g$X %*% beta[[s]])
      shape <- g$residual_prior[["shape"]] + nrow(g$X) / 2 + 1.5 * sum(beta[[s]] != 0)
      prior_squares <- sum(beta[[s]]^2 / g$prior_scale)
      augmented <- if (g$outcome_type == "binary") seq_len(nrow(g$X)) else if (g$outcome_type == "right.censored") which(g$status == 0) else integer()
      for (i in augmented) {
        rate <- g$residual_prior[["rate"]] + .5 * (sum((z[[s]][-i] - mean[-i])^2) + prior_squares)
        scale <- sqrt(rate / (shape - .5))
        df <- 2 * shape - 1
        sign <- if (g$outcome_type == "binary" && g$y[i] == 0) -1 else 1
        lower <- if (g$outcome_type == "binary") 0 else g$y[i]
        location <- sign * mean[i]
        log_tail <- stats::pt((lower - location) / scale, df, lower.tail = FALSE, log.p = TRUE)
        value <- location + scale * stats::qt(log(stats::runif(1)) + log_tail, df, lower.tail = FALSE, log.p = TRUE)
        if (!is.finite(value)) stop("Non-finite latent draw.")
        z[[s]][i] <- sign * max(lower, value)
      }
      for (j in seq_len(g$n_forced)) beta[[s]][j] <- draw_coefficient(conditional(s, j))
    }
    for (l in seq_along(blocks)) {
      block <- blocks[[l]]
      for (f in seq_along(block$features)) {
        parameters <- lapply(seq_along(block$members), function(i) conditional(block$members[i], block$columns[i, f]))
        log_bf <- vapply(parameters, `[[`, 0, "log_bf")
        log_weights <- .imr_marginal_energy(block, nu[l]) + drop(block$patterns %*% log_bf)
        pattern <- block$patterns[sample.int(nrow(block$patterns), 1, prob = exp(log_weights - max(log_weights))), ]
        for (i in seq_along(block$members)) {
          beta[[block$members[i]]][block$columns[i, f]] <- if (pattern[i]) draw_coefficient(parameters[[i]]) else 0
        }
      }
      if (length(block$features) >= 2L) {
        for (attempt in seq_len(ceiling(swap_rate * length(block$features)))) {
          pair <- sample.int(length(block$features), 2)
          proposed <- list()
          ratio <- 0
          for (i in seq_along(block$members)) {
            s <- block$members[i]
            columns <- block$columns[i, pair]
            selected <- beta[[s]][columns] != 0
            if (selected[1] == selected[2]) next
            old <- columns[selected]
            new <- columns[!selected]
            from <- conditional(s, old, columns)
            to <- conditional(s, new, columns)
            proposed[[length(proposed) + 1L]] <- list(s = s, old = old, new = new, parameters = to)
            ratio <- ratio + to$log_bf - from$log_bf
          }
          if (length(proposed)) {
            acceptance[1] <- acceptance[1] + 1
            if (log(stats::runif(1)) < min(0, ratio)) {
              for (move in proposed) {
                beta[[move$s]][move$old] <- 0
                beta[[move$s]][move$new] <- draw_coefficient(move$parameters)
              }
              acceptance[2] <- acceptance[2] + 1
            }
          }
        }
      }
      selected <- block$columns
      for (i in seq_along(block$members)) selected[i, ] <- beta[[block$members[i]]][block$columns[i, ]] != 0
      update <- .imr_update_marginal_interaction(
        block, selected, nu[l],
        interaction_prior, theta_step, model_variant == "imr"
      )
      block <- update$block
      acceptance[3:4] <- acceptance[3:4] + update$counts
      blocks[[l]] <- block
    }
    # Recovery never conditions a transition. Do it on every sweep so storage
    # thinning does not change the random-number stream of the primary chain.
    recovered_variance <- vapply(seq_along(groups), function(s) {
      g <- groups[[s]]
      shape <- g$residual_prior[["shape"]] + nrow(g$X) / 2 + 1.5 * sum(beta[[s]] != 0)
      rate <- g$residual_prior[["rate"]] + .5 * (sum((z[[s]] - drop(g$X %*% beta[[s]]))^2) +
        sum(beta[[s]]^2 / g$prior_scale))
      1 / stats::rgamma(1, shape, rate = rate)
    }, 0)
    if (iteration > burnin && (iteration - burnin) %% thin == 0) {
      row <- (iteration - burnin) %/% thin
      score <- 0
      for (s in seq_along(groups)) {
        g <- groups[[s]]
        active <- beta[[s]] != 0
        d <- sum(active)
        shape <- g$residual_prior[["shape"]] + nrow(g$X) / 2 + 1.5 * d
        rate <- g$residual_prior[["rate"]] + .5 * (sum((z[[s]] - drop(g$X %*% beta[[s]]))^2) + sum(beta[[s]]^2 / g$prior_scale))
        beta_draws[[s]][row, ] <- beta[[s]]
        sigma2_draws[row, s] <- recovered_variance[s]
        if (!is.null(latent_draws[[s]])) latent_draws[[s]][row, ] <- z[[s]]
        score <- score + g$residual_prior[["shape"]] * log(g$residual_prior[["rate"]]) -
          lgamma(g$residual_prior[["shape"]]) - (nrow(g$X) + d) / 2 * log(2 * pi) -
          1.5 * sum(log(g$prior_scale[active])) + 2 * sum(log(abs(beta[[s]][active]))) + lgamma(shape) - shape * log(rate)
      }
      for (l in seq_along(blocks)) {
        block <- blocks[[l]]
        patterns <- block$columns
        for (i in seq_along(block$members)) patterns[i, ] <- beta[[block$members[i]]][block$columns[i, ]] != 0
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
