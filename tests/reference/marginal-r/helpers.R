.imr_marginal_blocks <- function(groups, feature_platform, nu, model_variant,
                                 coefficient_columns = FALSE) {
  lapply(seq_along(nu), function(l) {
    features <- which(feature_platform == l)
    members <- which(vapply(groups, function(g) all(features %in% g$feature_index), TRUE))
    patterns <- if (length(members)) {
      as.matrix(expand.grid(rep(list(0:1), length(members))))
    } else {
      matrix(integer(), 1, 0)
    }
    columns <- vapply(features, function(f) {
      vapply(members, function(s) {
        offset <- if (coefficient_columns) groups[[s]]$n_forced else 0L
        as.integer(offset + match(f, groups[[s]]$feature_index))
      }, 1L)
    }, integer(length(members)))
    dim(columns) <- c(length(members), length(features))
    theta <- matrix(if (model_variant == "bms") 0 else .25, length(members), length(members))
    diag(theta) <- 0
    list(members = members, features = features, columns = columns, patterns = patterns, theta = theta)
  })
}

.imr_marginal_energy <- function(block, nu) {
  rowSums(block$patterns) * nu + rowSums((block$patterns %*% block$theta) * block$patterns)
}

.imr_log_sum_exp <- function(x) {
  largest <- max(x)
  largest + log(sum(exp(x - largest)))
}

.imr_update_marginal_interaction <- function(block, selected, nu, prior, step, sharing) {
  counts <- c(proposed = 0, accepted = 0)
  if (sharing && length(block$members) > 1L) {
    for (j in 2:length(block$members)) {
      for (i in seq_len(j - 1L)) {
        old <- block$theta[i, j]
        proposal <- exp(log(old) + stats::rnorm(1, 0, step))
        counts[1L] <- counts[1L] + 1
        if (!is.finite(proposal) || proposal <= 0) next
        candidate <- block
        candidate$theta[i, j] <- candidate$theta[j, i] <- proposal
        shared <- sum(selected[i, ] * selected[j, ])
        ratio <- 2 * (proposal - old) * shared - length(block$features) *
          (.imr_log_sum_exp(.imr_marginal_energy(candidate, nu)) -
            .imr_log_sum_exp(.imr_marginal_energy(block, nu))) +
          stats::dgamma(proposal, prior[1L], rate = prior[2L], log = TRUE) -
          stats::dgamma(old, prior[1L], rate = prior[2L], log = TRUE) + log(proposal) - log(old)
        if (log(stats::runif(1)) < min(0, ratio)) {
          block <- candidate
          counts[2L] <- counts[2L] + 1
        }
      }
    }
  }
  list(block = block, counts = counts)
}
