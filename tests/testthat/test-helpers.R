test_that("matrix preparation standardizes columns to mean 0, sd 1", {
  m <- matrix(rnorm(50, mean = 5, sd = 3), ncol = 5)
  z <- IntegMultiReg:::.imr_prepare_matrix(m)$normalized
  expect_equal(unname(colMeans(z)), rep(0, 5), tolerance = 1e-8)
  expect_equal(unname(apply(z, 2, stats::sd)), rep(1, 5), tolerance = 1e-8)
})

test_that("matrix preparation handles constant and empty columns gracefully", {
  m <- cbind(rnorm(10), rep(2, 10))
  z <- IntegMultiReg:::.imr_prepare_matrix(m)$normalized
  expect_true(all(z[, 2] == 0))          # zero-variance column -> all zeros
  empty <- IntegMultiReg:::.imr_prepare_matrix(matrix(numeric(0), nrow = 0, ncol = 3))$normalized
  expect_equal(dim(empty), c(0L, 3L))
})

test_that("known-mean/sd normalization inverts the training transform", {
  m <- matrix(rnorm(40), ncol = 4)
  prepared <- IntegMultiReg:::.imr_prepare_matrix(m)
  z <- IntegMultiReg:::normalize_matrix_known_mean_variance(m, prepared$mean, prepared$sd)
  expect_equal(z, prepared$normalized, tolerance = 1e-8)
})

test_that("subgroup_data partitions subjects by availability pattern", {
  d <- IntegMultiReg:::subgroup_data(simIMR$outcome.binary, simIMR$covariates,
                                     simIMR$platforms)
  expect_named(d, c("outcome", "covariate", "platform_data"))
  expect_setequal(names(d$platform_data), c("011", "100", "101", "111"))
})

test_that("subgroup_data aligns outcome, covariate and platform rows by id", {
  set.seed(99)
  platforms <- lapply(simIMR$platforms, function(x) x[sample(nrow(x)), ])
  outcome <- simIMR$outcome.binary[sample(nrow(simIMR$outcome.binary)), ]
  cov <- simIMR$covariates[sample(nrow(simIMR$covariates)), ]

  d <- IntegMultiReg:::subgroup_data(outcome, cov, platforms)

  for (pat in names(d$platform_data)) {
    ids <- d$outcome[[pat]]$id
    expect_identical(d$covariate[[pat]]$id, ids)
    for (platform in d$platform_data[[pat]]) {
      if (nrow(platform) > 0L) {
        expect_identical(platform$id, ids)
      }
    }
  }
})
