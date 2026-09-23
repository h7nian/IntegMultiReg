# Draw from density proportional to x^2 Normal(x; mu, sigma^2).
# The proposal uses (mu^2 + sigma^2*z^2) phi(z), an equally tight envelope
# in either sign mode: (mu + sigma*z)^2 <= 2*(mu^2 + sigma^2*z^2).
# This avoids a local random walk becoming trapped on one side of zero.
.imr_pmom_normal <- function(mu, sigma) {
  stopifnot(length(mu) == 1L, is.finite(mu), length(sigma) == 1L,
            is.finite(sigma), sigma > 0)
  scale <- max(abs(mu), sigma)
  a <- mu / scale
  b <- sigma / scale
  mixture <- a^2 / (a^2 + b^2)
  z <- .Call("imr_pmom_standardized_draw", a, b, mixture, PACKAGE = "IntegMultiReg")
  mu + sigma * z
}

# Conditional continuous-outcome model:
# y | beta, v ~ N(X beta, v I), v ~ IG(alpha, psi), and
# p(beta_j | v) = [beta_j^2/(h_j*v)]^moment_j N(0, h_j*v).
# IG uses a shape/rate parameterization: 1/v ~ Gamma(shape, rate).
.imr_conditional_chain <- function(X, y, prior_scale, moment, alpha, psi,
                             draws = 2000L, burnin = 1000L,
                             initial_beta = NULL, initial_variance = NULL,
                             outcome_type = "continuous", status = NULL) {
  stopifnot(is.matrix(X), is.numeric(X), all(is.finite(X)),
            is.numeric(y), length(y) == nrow(X), all(is.finite(y)),
            nrow(X) > 0L, ncol(X) > 0L,
            length(prior_scale) == ncol(X), all(is.finite(prior_scale)),
            all(prior_scale > 0), length(moment) == ncol(X),
            all(moment %in% c(0, 1)),
            length(alpha) == 1L, is.finite(alpha), alpha > 0,
            length(psi) == 1L, is.finite(psi), psi > 0,
            length(draws) == 1L, draws >= 1L, draws == as.integer(draws),
            length(burnin) == 1L, burnin >= 0L, burnin == as.integer(burnin))
  if (!outcome_type %in% c("continuous", "binary", "right.censored")) {
    stop("Unsupported conditional outcome type.")
  }
  observed <- y
  if (outcome_type == "binary") {
    stopifnot(all(y %in% c(0, 1)))
    y <- ifelse(y == 1, .5, -.5)
  }
  if (outcome_type == "right.censored") {
    stopifnot(length(status) == length(y), all(status %in% c(0, 1)))
  }
  p <- ncol(X)
  precision <- crossprod(X) + diag(1 / prior_scale, nrow = p)
  rhs <- drop(crossprod(X, y))
  beta <- if (is.null(initial_beta)) drop(solve(precision, rhs)) else initial_beta
  v <- if (is.null(initial_variance)) {
    (psi + sum((y - X %*% beta)^2) / 2) / (alpha + length(y) / 2)
  } else initial_variance
  stopifnot(length(beta) == p, all(is.finite(beta)),
            length(v) == 1L, is.finite(v), v > 0)
  out <- matrix(NA_real_, draws, p + 1L)
  colnames(out) <- c(if (is.null(colnames(X))) paste0("beta", seq_len(p)) else
                    colnames(X), "variance")
  shape <- alpha + (length(y) + p) / 2 + sum(moment)
  for (iter in seq_len(burnin + draws)) {
    if (outcome_type != "continuous") {
      eta <- drop(X %*% beta)
      if (outcome_type == "binary") {
        direction <- ifelse(observed == 1, 1, -1)
        y <- direction * .imr_lower_normal(direction * eta, sqrt(v), 0)
      } else {
        censored <- which(status == 0)
        y[censored] <- .imr_lower_normal(eta[censored], sqrt(v), observed[censored])
      }
      rhs <- drop(crossprod(X, y))
    }
    for (j in seq_len(p)) {
      mu <- (rhs[j] - sum(precision[j, ] * beta) +
               precision[j, j] * beta[j]) / precision[j, j]
      sd <- sqrt(v / precision[j, j])
      beta[j] <- if (moment[j] == 1L) .imr_pmom_normal(mu, sd) else stats::rnorm(1L, mu, sd)
    }
    rate <- psi + (sum((y - drop(X %*% beta))^2) +
                    sum(beta^2 / prior_scale)) / 2
    v <- 1 / stats::rgamma(1L, shape = shape, rate = rate)
    if (any(!is.finite(beta)) || !is.finite(v) || v <= 0) {
      stop("Non-finite posterior state; check the data and prior scales.")
    }
    if (iter > burnin) out[iter - burnin, ] <- c(beta, v)
  }
  out
}

# Classical split R-hat. This is a diagnostic, not proof of convergence.
.imr_split_rhat <- function(chains) {
  stopifnot(length(chains) >= 2L,
            all(vapply(chains, nrow, integer(1)) == nrow(chains[[1L]])))
  n <- floor(nrow(chains[[1L]]) / 2)
  stopifnot(n >= 2L)
  halves <- unlist(lapply(chains, function(x) list(
    x[seq_len(n), , drop = FALSE],
    x[nrow(x) - n + seq_len(n), , drop = FALSE]
  )), recursive = FALSE)
  within <- Reduce(`+`, lapply(halves, function(x) apply(x, 2, stats::var))) / length(halves)
  means <- do.call(rbind, lapply(halves, colMeans))
  between <- n * apply(means, 2, stats::var)
  sqrt(((n - 1) / n * within + between / n) / within)
}

# Inverse survival-CDF sampling remains stable far into a Normal tail.
.imr_lower_normal <- function(mu, sd, lower) {
  if (!length(mu)) return(numeric())
  a <- (lower - mu) / sd
  log_tail <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
  z <- stats::qnorm(log(pmax(stats::runif(length(mu)), .Machine$double.xmin)) +
                    log_tail, lower.tail = FALSE, log.p = TRUE)
  mu + sd * z
}
