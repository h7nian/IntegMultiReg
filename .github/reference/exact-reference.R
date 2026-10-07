# Independent small-model reference: integrate all coefficients and variances.
# Gaussian raw moments are polynomials in v when beta|v has covariance v*A^-1.
normal_moment_polynomial <- function(index, location, covariance) {
  if (!length(index)) return(1)
  first <- index[1]; rest <- index[-1]
  value <- location[first] * normal_moment_polynomial(rest, location, covariance)
  if (length(rest)) for (j in seq_along(rest)) {
    term <- covariance[first, rest[j]] * c(0, normal_moment_polynomial(rest[-j], location, covariance))
    length(value) <- max(length(value), length(term)); value[is.na(value)] <- 0
    value[seq_along(term)] <- value[seq_along(term)] + term
  }
  value
}

integrated_moment <- function(polynomial, shape, rate) {
  keep <- which(polynomial != 0)
  if (!length(keep)) return(list(sign = 0, log_absolute = -Inf))
  powers <- keep - 1L
  stopifnot(all(shape - powers > 0))
  log_terms <- log(abs(polynomial[keep])) + lgamma(shape - powers) -
    (shape - powers) * log(rate)
  scale <- max(log_terms)
  total <- sum(sign(polynomial[keep]) * exp(log_terms - scale))
  list(sign = sign(total), log_absolute = scale + log(abs(total)))
}

gaussian_pmom_evidence <- function(X, y, prior_scale, residual_prior) {
  n <- nrow(X); d <- ncol(X)
  precision <- crossprod(X) + diag(1 / prior_scale, nrow = d)
  covariance <- solve(precision)
  location <- drop(covariance %*% crossprod(X, y))
  shape <- residual_prior[["shape"]] + n / 2 + d
  rate <- residual_prior[["rate"]] +
    (sum(y^2) - sum(drop(crossprod(X, y)) * location)) / 2
  stopifnot(rate > 0)
  index <- rep(seq_len(d), each = 2L)
  polynomial <- normal_moment_polynomial(index, location, covariance)
  denominator <- integrated_moment(polynomial, shape, rate)
  stopifnot(denominator$sign > 0)
  log_constant <- -n / 2 * log(2 * pi) -
    (as.numeric(determinant(precision, logarithm = TRUE)$modulus) + sum(log(prior_scale))) / 2 -
    sum(log(prior_scale)) + residual_prior[["shape"]] * log(residual_prior[["rate"]]) -
    lgamma(residual_prior[["shape"]])
  ratio <- function(index_extra = integer(), variance_power = 0L) {
    numerator <- integrated_moment(normal_moment_polynomial(c(index, index_extra), location, covariance),
                                   shape - variance_power, rate)
    numerator$sign * exp(numerator$log_absolute - denominator$log_absolute)
  }
  list(log_evidence = log_constant + denominator$log_absolute,
       beta_mean = vapply(seq_len(d), function(j) ratio(j), 0),
       beta_second = vapply(seq_len(d), function(j) ratio(c(j, j)), 0),
       variance_mean = ratio(variance_power = 1L))
}

exact_two_group_posterior <- function(groups, nu = -1, interaction_prior = c(shape = 2, rate = 2),
                                      model_variant = c("imr", "bms")) {
  model_variant <- match.arg(model_variant)
  stopifnot(length(groups) == 2L, all(vapply(groups, function(g) g$outcome_type == "continuous", TRUE)))
  p <- length(groups[[1]]$feature_index)
  stopifnot(identical(groups[[1]]$feature_index, groups[[2]]$feature_index))
  states <- as.matrix(expand.grid(rep(list(0:1), 2 * p)))
  log_weight <- theta_mean <- numeric(nrow(states))
  coefficient_means <- lapply(groups, function(g) matrix(0, nrow(states), ncol(g$X)))
  coefficient_seconds <- coefficient_means
  variance_means <- matrix(0, nrow(states), 2)
  for (row in seq_len(nrow(states))) {
    selection <- matrix(states[row, ], nrow = 2, byrow = TRUE)
    if (model_variant == "bms") {
      prior_mass <- exp(nu * sum(selection) - 2 * p * log1p(exp(nu)))
      theta_mean[row] <- 0
    } else {
      # Closed two-group normalizer; no sampler normalizer routine is reused.
      log_prior <- function(theta) {
        logs <- cbind(0, nu, nu, 2 * nu + 2 * theta)
        largest <- apply(logs, 1, max)
        log_z <- largest + log(rowSums(exp(logs - largest)))
        nu * sum(selection) + 2 * theta * sum(selection[1, ] * selection[2, ]) - p * log_z +
          dgamma(theta, shape = interaction_prior[["shape"]], rate = interaction_prior[["rate"]], log = TRUE)
      }
      prior_mass <- integrate(function(theta) exp(log_prior(theta)), 0, Inf,
                              rel.tol = 1e-10, subdivisions = 500L)$value
      theta_mean[row] <- integrate(function(theta) theta * exp(log_prior(theta)), 0, Inf,
                                  rel.tol = 1e-10, subdivisions = 500L)$value / prior_mass
    }
    log_weight[row] <- log(prior_mass)
    for (g in 1:2) {
      dat <- groups[[g]]
      active <- c(seq_len(dat$n_forced), dat$n_forced + which(selection[g, ] == 1))
      evidence <- gaussian_pmom_evidence(dat$X[, active, drop = FALSE], dat$y,
                                         dat$prior_scale[active], dat$residual_prior)
      log_weight[row] <- log_weight[row] + evidence$log_evidence
      coefficient_means[[g]][row, active] <- evidence$beta_mean
      coefficient_seconds[[g]][row, active] <- evidence$beta_second
      variance_means[row, g] <- evidence$variance_mean
    }
  }
  probability <- exp(log_weight - max(log_weight)); probability <- probability / sum(probability)
  list(states = states, probability = probability,
       inclusion = drop(crossprod(probability, states)),
       beta_mean = lapply(coefficient_means, function(x) drop(crossprod(probability, x))),
       beta_second = lapply(coefficient_seconds, function(x) drop(crossprod(probability, x))),
       variance_mean = drop(crossprod(probability, variance_means)),
       theta_mean = sum(probability * theta_mean))
}
