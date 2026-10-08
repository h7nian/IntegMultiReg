.imr_marginal_blocks <- function(groups, feature_platform, nu, model_variant) {
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
        offset <- groups[[s]]$n_forced
        as.integer(offset + match(f, groups[[s]]$feature_index))
      }, 1L)
    }, integer(length(members)))
    dim(columns) <- c(length(members), length(features))
    theta <- matrix(if (model_variant == "bms") 0 else .25, length(members), length(members))
    diag(theta) <- 0
    list(members = members, features = features, columns = columns, patterns = patterns, theta = theta)
  })
}
