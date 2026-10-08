# Independent observed-data quadrature, without augmented responses or MCMC.
normal_quadrature <- function(n) {
  J <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) J[i, i+1L] <- J[i+1L, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1, ]^2)
}
gamma_quadrature <- function(n, shape, rate) {
  J <- diag(2 * (seq_len(n)-1) + shape)
  for (i in seq_len(n-1L)) J[i, i+1L] <- J[i+1L, i] <- sqrt(i * (i + shape - 1))
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values / rate, weights = e$vectors[1, ]^2)
}
observed_pmom_evidence <- function(dat, active, normal_order = 35L, gamma_order = 60L) {
  X <- dat$X[, active, drop = FALSE]; tau <- dat$prior_scale[active]
  a <- dat$residual_prior[['shape']]; b <- dat$residual_prior[['rate']]
  d <- length(active)
  q <- normal_quadrature(normal_order)
  indices <- as.matrix(expand.grid(rep(list(seq_len(normal_order)), d)))
  z <- matrix(q$nodes[indices], nrow(indices), d)
  u <- sweep(z, 2, sqrt(tau), '*')
  log_weight <- rowSums(log(matrix(q$weights[indices], nrow(indices), d)) + log(z^2))
  eta <- X %*% t(u)
  if (dat$outcome_type == 'binary') {
    direction <- ifelse(dat$y == 1, 1, -1)
    log_w <- log_weight + colSums(pnorm(direction * eta, log.p = TRUE))
    scale <- max(log_w); weight <- exp(log_w - scale); total <- sum(weight)
    probability <- weight / total
    root_v <- exp(.5 * log(b) + lgamma(a - .5) - lgamma(a))
    mean_v <- b / (a - 1)
    latent <- eta + direction * exp(dnorm(eta, log = TRUE) - pnorm(direction * eta, log.p = TRUE))
    return(list(log_evidence = scale + log(total),
      beta_mean = root_v * drop(crossprod(probability, u)),
      beta_second = mean_v * drop(crossprod(probability, u^2)),
      variance_mean = mean_v,
      latent_mean = root_v * drop(latent %*% probability)))
  }
  stopifnot(dat$outcome_type == 'right.censored')
  observed <- which(dat$status == 1L); censored <- which(dat$status == 0L)
  shape <- a + length(observed)/2
  rate <- b + sum(dat$y[observed]^2)/2
  g <- gamma_quadrature(gamma_order, shape, rate)
  log_constant <- -length(observed)/2 * log(2*pi) + a*log(b) - lgamma(a) + lgamma(shape) - shape*log(rate)
  linear <- if(length(observed)) drop(crossprod(dat$y[observed], eta[observed, , drop=FALSE])) else rep(0,nrow(u))
  square <- if(length(observed)) colSums(eta[observed, , drop=FALSE]^2)/2 else rep(0,nrow(u))
  # Each precision node has its own scale; accumulate positive weighted moments
  # under a common scale chosen from a preliminary pass.
  log_component <- function(k) {
    root <- sqrt(g$nodes[k])
    likelihood <- root * linear - square
    if(length(censored)) likelihood <- likelihood + colSums(pnorm(
      eta[censored, , drop=FALSE] - dat$y[censored]*root, log.p=TRUE))
    log_weight + log(g$weights[k]) + likelihood
  }
  scale <- max(vapply(seq_along(g$nodes),function(k)max(log_component(k)),0))
  total <- 0; beta_mean <- beta_second <- numeric(d); variance_mean <- 0
  latent_mean <- numeric(nrow(X))
  for(k in seq_along(g$nodes)) {
    precision <- g$nodes[k]; root <- sqrt(precision)
    weight <- exp(log_component(k) - scale)
    mass <- sum(weight); total <- total + mass
    beta_mean <- beta_mean + drop(crossprod(weight,u))/root
    beta_second <- beta_second + drop(crossprod(weight,u^2))/precision
    variance_mean <- variance_mean + mass/precision
    if(length(observed)) latent_mean[observed] <- latent_mean[observed] + dat$y[observed]*mass
    if(length(censored)) {
      standardized <- eta[censored, , drop=FALSE] - dat$y[censored]*root
      conditional <- (eta[censored, , drop=FALSE] +
        exp(dnorm(standardized,log=TRUE)-pnorm(standardized,log.p=TRUE)))/root
      latent_mean[censored] <- latent_mean[censored] + drop(conditional %*% weight)
    }
  }
  list(log_evidence=log_constant+scale+log(total), beta_mean=beta_mean/total,
       beta_second=beta_second/total,variance_mean=variance_mean/total,
       latent_mean=latent_mean/total)
}

observed_two_group_reference <- function(groups, nu=-1,
 interaction_prior=c(shape=2,rate=2), normal_order=35L,gamma_order=60L) {
  stopifnot(length(groups)==2L)
  p <- length(groups[[1]]$feature_index)
  local_states <- as.matrix(expand.grid(rep(list(0:1),p)))
  cache <- lapply(groups,function(dat) lapply(seq_len(nrow(local_states)),function(row) {
    active <- c(seq_len(dat$n_forced),dat$n_forced+which(local_states[row,]==1))
    list(active=active, moments=observed_pmom_evidence(dat,active,normal_order,gamma_order))
  }))
  states <- as.matrix(expand.grid(rep(list(0:1),2*p)))
  log_weight <- theta_mean <- numeric(nrow(states))
  beta_mean <- beta_second <- lapply(groups,function(g)matrix(0,nrow(states),ncol(g$X)))
  variance_mean <- matrix(0,nrow(states),2)
  latent_mean <- lapply(groups,function(g)matrix(0,nrow(states),nrow(g$X)))
  for(row in seq_len(nrow(states))) {
    selection <- matrix(states[row,],2,byrow=TRUE)
    density <- function(theta) {
      terms <- cbind(0,nu,nu,2*nu+2*theta)
      maximum <- apply(terms,1,max)
      logz <- maximum+log(rowSums(exp(terms-maximum)))
      exp(nu*sum(selection)+2*theta*sum(selection[1,]*selection[2,])-p*logz+
        dgamma(theta,interaction_prior[['shape']],rate=interaction_prior[['rate']],log=TRUE))
    }
    prior <- integrate(density,0,Inf,rel.tol=1e-10,subdivisions=500)$value
    theta_mean[row] <- integrate(function(x)x*density(x),0,Inf,rel.tol=1e-10,subdivisions=500)$value/prior
    log_weight[row] <- log(prior)
    for(g in 1:2) {
      at <- 1L+sum(selection[g,]*2^(seq_len(p)-1L))
      item <- cache[[g]][[at]]; m <- item$moments
      log_weight[row] <- log_weight[row]+m$log_evidence
      beta_mean[[g]][row,item$active] <- m$beta_mean
      beta_second[[g]][row,item$active] <- m$beta_second
      variance_mean[row,g] <- m$variance_mean
      latent_mean[[g]][row,] <- m$latent_mean
    }
  }
  probability <- exp(log_weight-max(log_weight));probability<-probability/sum(probability)
  list(states=states,probability=probability,inclusion=drop(crossprod(probability,states)),
    beta_mean=lapply(beta_mean,function(x)drop(crossprod(probability,x))),
    beta_second=lapply(beta_second,function(x)drop(crossprod(probability,x))),
    variance_mean=drop(crossprod(probability,variance_mean)),
    theta_mean=sum(probability*theta_mean),
    latent_mean=lapply(latent_mean,function(x)drop(crossprod(probability,x))),
    quadrature=list(normal_order=normal_order,gamma_order=gamma_order))
}
