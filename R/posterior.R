# Canonical draw blocks. Arrays preserve iteration/chain/parameter identity;
# inclusion indicators are derived from the exact exclusion zeros in beta.
.imr_parameter_blocks <- function(object, parm = "all") {
  families <- c("coefficients", "variance", "selection", "interaction", "latent")
  if (identical(parm, "all")) parm <- families
  if (!is.character(parm) || !length(parm) || any(!parm %in% families) || anyDuplicated(parm)) {
    .imr_abort("`parm` must name coefficients, variance, selection, interaction, latent, or all.")
  }
  output <- list()
  add <- function(family, group, draws) {
    output[[length(output) + 1L]] <<- list(family = family, group = group, draws = draws)
  }
  model <- object$model
  post <- object$posterior
  if ("coefficients" %in% parm) {
    for (g in model$subgroup_names) {
      add("coefficients", g, post$coefficients[[g]])
    }
  }
  if ("variance" %in% parm) {
    for (g in seq_along(model$subgroup_names)) {
      value <- post$variance[, , g, drop = FALSE]
      dimnames(value)[[3L]] <- "residual_variance"
      add("variance", model$subgroup_names[g], value)
    }
  }
  if ("selection" %in% parm) {
    for (l in seq_len(model$n_platforms)) {
      members <- model$platform_subgroups[[l]]
      features <- model$feature_names[[l]]
      value <- array(0L, c(
        object$control$mcmc$draws, object$control$mcmc$chains,
        length(members) * length(features)
      ))
      labels <- character(dim(value)[3L])
      for (i in seq_along(members)) {
        g <- members[i]
        columns <- .imr_feature_columns(model, g, l)
        target <- (i - 1L) * length(features) + seq_along(features)
        value[, , target] <- post$coefficients[[g]][, , columns, drop = FALSE] != 0
        labels[target] <- paste(model$subgroup_names[g], features, sep = ":")
      }
      dimnames(value) <- list(iteration = NULL, chain = paste0("chain", seq_len(dim(value)[2L])), parameter = labels)
      add("selection", model$platform_names[l], value)
    }
  }
  if ("interaction" %in% parm) {
    for (l in seq_len(model$n_platforms)) {
      add("interaction", model$platform_names[l], post$interaction[[l]])
    }
  }
  if ("latent" %in% parm && !is.null(post$latent)) {
    for (g in model$subgroup_names) {
      if (!is.null(post$latent[[g]])) add("latent", g, post$latent[[g]])
    }
  }
  output
}

#' Extract Stored Joint Posterior Draws
#'
#' Returns retained samples from a fitted joint chain. This function performs
#' no sampling and does not alter the random-number state.
#' @param object A joint `imr` fit.
#' @param parm `"coefficients"`, `"variance"`, `"selection"`, `"interaction"`,
#'   `"latent"` or `"all"`. The default is regression coefficients.
#' @return Named arrays with dimensions iteration by chain by parameter,
#'   grouped by availability subgroup or platform. `parm = "all"` returns a
#'   list of these families. Selection indicators are derived as `beta != 0`;
#'   they occupy no duplicate history in the fitted object. Latent responses
#'   require `imr_mcmc(keep_latent = TRUE)` for a binary or censored fit.
#' @examples
#' # Short interface example; increase the budget and inspect diagnostics for inference.
#' x <- data.frame(id = 1:20, marker = sin(1:20))
#' y <- data.frame(id = x$id, y = 1 + x$marker + cos(x$id) / 3)
#' fit <- imr(list(assay = x), y,
#'   outcome_type = "continuous",
#'   min_subgroup_size = 0, priors = imr_priors(forced_scale = 1),
#'   mcmc = imr_mcmc(
#'     draws = 100, burnin = 100, chains = 2,
#'     seed = 1, diagnostics = FALSE
#'   )
#' )
#' dim(posterior_draws(fit)[[1]])
#' @export
posterior_draws <- function(object, parm = "coefficients") {
  .imr_check_fit(object)
  parm <- match.arg(parm, c("coefficients", "variance", "selection", "interaction", "latent", "all"))
  if (parm == "latent" && (is.null(object$posterior$latent) || object$control$outcome_type == "continuous")) {
    .imr_abort("This fit has no stored augmented responses; fit binary/censored data with `imr_mcmc(keep_latent = TRUE)`.")
  }
  blocks <- .imr_parameter_blocks(object, parm)
  by_family <- split(blocks, vapply(blocks, `[[`, "", "family"))
  result <- lapply(by_family, function(x) {
    stats::setNames(
      lapply(x, `[[`, "draws"),
      vapply(x, `[[`, "", "group")
    )
  })
  if (parm == "all") result else result[[parm]]
}

# Empirical inverse CDF, with the same floating-point boundary tolerance for
# every parameter family and prediction. This prevents (1-.95)/2 from selecting
# a different order statistic than .025 at an exact empirical-CDF boundary.
.imr_quantile <- function(x, probability, weights = NULL) {
  if (!length(x) || any(!is.finite(x))) .imr_abort("Posterior quantiles need finite draws.")
  if (is.null(weights)) {
    position <- length(x) * probability
    index <- pmax(1L, pmin(length(x), ceiling(position - 4 * .Machine$double.eps * pmax(1, abs(position)))))
    return(sort(x, partial = unique(index))[index])
  }
  keep <- weights > 0
  x <- x[keep]
  weights <- weights[keep]
  order <- order(x)
  cumulative <- cumsum(weights[order]) / sum(weights)
  cumulative[length(cumulative)] <- 1
  vapply(probability, function(p) x[order[which(cumulative >= p - 4 * .Machine$double.eps)[1L]]], 0)
}

.imr_draw_table <- function(draws, level) {
  parameters <- dimnames(draws)[[3L]]
  if (!length(parameters)) {
    return(data.frame(
      parameter = character(), mean = numeric(), sd = numeric(),
      lower = numeric(), median = numeric(), upper = numeric()
    ))
  }
  alpha <- (1 - level) / 2
  rows <- lapply(seq_along(parameters), function(j) {
    x <- as.vector(draws[, , j])
    interval <- .imr_quantile(x, c(alpha, .5, 1 - alpha))
    data.frame(
      parameter = parameters[j], mean = mean(x), sd = stats::sd(x),
      lower = interval[1L], median = interval[2L], upper = interval[3L], stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.imr_check_level <- function(level) {
  level <- .imr_check_numeric_vector(level, "level", 1, positive = TRUE)
  if (level >= 1) .imr_abort("`level` must be less than one.")
  level
}

.imr_posterior_tables <- function(object, parm, level, include_diagnostics = TRUE) {
  blocks <- .imr_parameter_blocks(object, parm)
  diagnostics <- if (include_diagnostics) mcmc_diagnostics(object) else NULL
  tables <- lapply(blocks, function(block) {
    x <- .imr_draw_table(block$draws, level)
    if (!nrow(x)) {
      return(cbind(family = character(), group = character(), x))
    }
    if (!is.null(diagnostics)) {
      candidates <- diagnostics[diagnostics$family == block$family & diagnostics$group == block$group, , drop = FALSE]
      x <- cbind(x, candidates[match(x$parameter, candidates$parameter),
        c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "status"),
        drop = FALSE
      ])
    }
    cbind(family = block$family, group = block$group, x, stringsAsFactors = FALSE)
  })
  # Empty platform tables may have no diagnostics columns. Keep their absence
  # explicit rather than manufacture rows for parameters that do not exist.
  tables <- tables[vapply(tables, nrow, 1L) > 0L]
  if (!length(tables)) {
    return(data.frame(
      family = character(), group = character(), parameter = character(),
      mean = numeric(), sd = numeric(), lower = numeric(), median = numeric(), upper = numeric()
    ))
  }
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}
