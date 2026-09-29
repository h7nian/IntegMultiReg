# Test-only R reference frozen from bd4c6b2; never loaded by package users.
pmom_r_reference <- function(mu, sigma) {
  scale <- max(abs(mu), sigma)
  a <- mu / scale
  b <- sigma / scale
  mixture <- a^2 / (a^2 + b^2)
  repeat {
    z <- if (stats::runif(1L) < mixture) stats::rnorm(1L) else {
      (if (stats::runif(1L) < 0.5) -1 else 1) * sqrt(stats::rchisq(1L, df = 3))
    }
    if (stats::runif(1L) < (a + b * z)^2 / (2 * (a^2 + b^2 * z^2)))
      return(mu + sigma * z)
  }
}

test_that("native pMOM draws preserve R RNG kinds and follow-on distributions", {
  original_kind <- RNGkind()
  original_seed <- IntegMultiReg:::.imr_save_rng()
  on.exit({
    do.call(RNGkind, as.list(original_kind))
    IntegMultiReg:::.imr_restore_rng(original_seed)
  }, add = TRUE)
  uniform_kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
    "Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  normal_kinds <- c("Inversion", "Box-Muller", "Kinderman-Ramage", "Ahrens-Dieter",
                    "Buggy Kinderman-Ramage")
  parameters <- list(c(0, 1), c(1, 1), c(-1, 1), c(100, .1), c(-100, .1),
    c(1e300, 1e-300), c(-1e300, 1e300), c(1e-150, 1e150),
    c(-0, .Machine$double.xmin * .Machine$double.eps))
  sample_state <- function(sampler, parameter) {
    value <- replicate(8L, sampler(parameter[1L], parameter[2L]))
    rng <- .Random.seed
    tail <- c(runif(8L), rnorm(8L), rchisq(8L, 3))
    list(value = value, rng = rng, tail = tail, tail_rng = .Random.seed)
  }
  for (uniform in uniform_kinds) for (normal in normal_kinds) {
    messages <- character()
    # Old RNG kinds intentionally warn. Assert their exact expected count and
    # recognized messages, rather than hiding an unexpected sampler warning.
    withCallingHandlers(RNGkind(uniform, normal), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expected_warnings <- as.integer(normal == "Buggy Kinderman-Ramage") +
      as.integer(uniform == "Marsaglia-Multicarry") +
      as.integer(uniform == "Marsaglia-Multicarry" &&
                   normal %in% c("Kinderman-Ramage", "Ahrens-Dieter"))
    expect_length(messages, expected_warnings)
    expect_true(all(grepl("buggy version of Kinderman-Ramage|poor statistical properties|deviations from normality",
                         messages)))
    for (parameter in parameters) {
      set.seed(729)
      expected <- sample_state(pmom_r_reference, parameter)
      set.seed(729)
      actual <- sample_state(IntegMultiReg:::.imr_pmom_normal, parameter)
      expect_true(identical(actual, expected, num.eq = FALSE),
                  info = paste(uniform, normal, paste(parameter, collapse = "/")))
    }
  }
})

test_that("native pMOM input failures do not advance the RNG", {
  native <- function(a, b, mixture) .Call("imr_pmom_standardized_draw",
    a, b, mixture, PACKAGE = "IntegMultiReg")
  invalid <- list(list(NULL, 1, 0), list(1L, 1, .5), list(c(1, 1), 1, .5),
    list(NA_real_, 1, .5), list(Inf, 1, .5), list(0, 0, .5),
    list(.5, .5, .5), list(2, 1, .5), list(1, -1, .5),
    list(1, 1, NA_real_), list(1, 1, -1), list(1, 1, 2))
  set.seed(73)
  rng <- .Random.seed
  for (arguments in invalid) {
    expect_error(do.call(native, arguments), "Invalid normalized pMOM inputs")
    expect_identical(.Random.seed, rng)
  }
  for (arguments in list(list(NA_real_, 1), list(1, 0), list(c(1, 2), 1))) {
    expect_error(do.call(IntegMultiReg:::.imr_pmom_normal, arguments))
    expect_identical(.Random.seed, rng)
  }
})
