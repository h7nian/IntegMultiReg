# Data preparation is shared; C owns transitions, recovery and sample storage.
.imr_marginal_chain <- function(groups, feature_platform, nu, draws, burnin, seed,
                                model_variant, initial, keep_latent, interaction_prior,
                                thin, theta_step, swap_rate, initial_state, verbose,
                                marginalize, max_nodes = 200000L, variance_step = NULL) {
  set.seed(seed)
  blocks <- .imr_marginal_blocks(groups, feature_platform, nu, model_variant)
  groups <- lapply(groups, function(g) {
    storage.mode(g$X) <- "double"
    g$y <- as.double(g$y)
    g$status <- as.integer(g$status)
    g$n_forced <- .imr_check_integer_scalar(g$n_forced, "n_forced", min = 1, max = ncol(g$X))
    g$prior_scale <- as.double(g$prior_scale)
    g$residual_prior <- as.double(g$residual_prior[c("shape", "rate")])
    g
  })
  if (is.null(initial_state)) {
    coefficients <- lapply(seq_along(groups), function(s) {
      g <- groups[[s]]
      selected <- switch(initial,
        empty = rep(0, length(g$feature_index)),
        full = rep(1, length(g$feature_index)),
        alternating = as.integer((g$feature_index + s) %% 2 == 0)
      )
      c(rep(.1, g$n_forced), .1 * selected)
    })
    initial_state <- list(
      coefficients = coefficients,
      variance = vapply(groups, function(g) g$residual_prior[2] / (g$residual_prior[1] + 1), 0),
      interaction = lapply(blocks, `[[`, "theta")
    )
  }
  initial_state$coefficients <- lapply(initial_state$coefficients, as.double)
  initial_state$variance <- as.double(initial_state$variance)
  initial_state$interaction <- lapply(initial_state$interaction, function(x) {
    storage.mode(x) <- "double"
    x
  })
  blocks <- lapply(blocks, function(block) {
    storage.mode(block$patterns) <- "double"
    block
  })
  if (is.null(variance_step)) {
    variance_step <- vapply(groups, function(g) {
      min(.4, 1.5 / sqrt(g$residual_prior[1] + nrow(g$X) / 2))
    }, 0)
  }
  if (length(variance_step) == 1L) variance_step <- rep(variance_step, length(groups))
  settings <- list(
    draws = draws, burnin = burnin, thin = thin,
    sharing = as.integer(model_variant == "imr"), verbose = as.integer(verbose),
    keep_latent = as.integer(keep_latent),
    coefficients_marginalized = as.integer(marginalize == "coefficients"),
    theta_step = theta_step, swap_rate = swap_rate,
    interaction_shape = interaction_prior[1], interaction_rate = interaction_prior[2],
    variance_step = as.double(variance_step), nu = as.double(nu)
  )
  algebra <- list(
    multiply = get("%*%", baseenv()), crossprod = base::crossprod,
    largest_eigenvalue = function(x) eigen(x, symmetric = TRUE, only.values = TRUE)$values[1],
    prepare_model = if (marginalize == "coefficients") .imr_marginal_model_factory(groups, max_nodes) else NULL
  )
  .Call("imr_marginal_sample", groups, blocks, settings, initial_state, algebra,
    PACKAGE = "IntegMultiReg"
  )
}

# Static factorizations are prepared lazily and retained by the native cache.
# R's existing BLAS/LAPACK entry points preserve the reference numerical policy.
.imr_marginal_model_factory <- function(groups, max_nodes) {
  quadrature_cache <- new.env(parent = emptyenv())
  gaussian_rule <- function(d) {
    key <- as.character(d)
    if (exists(key, quadrature_cache, inherits = FALSE)) {
      return(get(key, quadrature_cache))
    }
    orders <- seq.int(d + 1L, 2L)
    if (sum(log(orders)) > log(max_nodes)) {
      .imr_abort("Exact coefficient integration exceeds its node budget; no approximation or model truncation was used.")
    }
    rules <- lapply(orders, function(q) {
      jacobi <- matrix(0, q, q)
      jacobi[cbind(1:(q - 1L), 2:q)] <- sqrt(1:(q - 1L))
      result <- eigen(jacobi + t(jacobi), symmetric = TRUE)
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
  function(s, selected) {
    g <- groups[[s]]
    active <- c(seq_len(g$n_forced), g$n_forced + which(selected != 0))
    X <- g$X[, active, drop = FALSE]
    tau <- g$prior_scale[active]
    d <- length(active)
    factor <- chol(crossprod(X) + diag(1 / tau, nrow = d))
    covariance <- chol2inv(factor)
    normal_factor <- chol(covariance)
    rule <- gaussian_rule(d)
    a <- g$residual_prior[1]
    b <- g$residual_prior[2]
    list(
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
}
