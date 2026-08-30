test_that("imr_data validates and summarizes multi-platform inputs", {
  dat <- imr_data(
    simIMR$platforms,
    outcome = simIMR$outcome.binary,
    covariates = simIMR$covariates,
    type_outcome = "binary"
  )
  expect_s3_class(dat, "imr_data")
  expect_true(validate_imr_data(dat))
  expect_setequal(names(dat$subgroup_sizes), c("011", "100", "101", "111"))
  expect_equal(sum(dat$subgroup_sizes), 300L)
  expect_s3_class(summary(dat), "summary.imr_data")
  expect_identical(as.data.frame(dat), dat$availability)
  expect_output(print(dat), "Validated IMR")
})

test_that("imr_data detects corrupted availability metadata", {
  dat <- imr_data(simIMR$platforms)
  dat$availability$subgroup[1L] <- "000"
  expect_error(validate_imr_data(dat), "does not match")
})

test_that("imr_data reports subjects lacking required aligned rows", {
  outcome <- simIMR$outcome.binary[-1L, ]
  dat <- imr_data(
    simIMR$platforms, outcome = outcome, type_outcome = "binary"
  )
  expect_equal(dat$n_platform_subjects, 300L)
  expect_equal(nrow(dat$availability), 299L)
  expect_equal(length(dat$excluded_ids), 1L)
  expect_output(print(dat), "299 eligible of 300")
})

test_that("imr_data supports a non-standard identifier name", {
  platforms <- lapply(simIMR$platforms, function(x) {
    names(x)[1L] <- "subject"
    x
  })
  outcome <- simIMR$outcome.continuous
  names(outcome)[1L] <- "subject"
  dat <- imr_data(
    platforms, outcome = outcome,
    type_outcome = "continuous", id = "subject"
  )
  expect_identical(names(dat$platforms[[1L]])[1L], "id")
  expect_identical(names(dat$outcome)[1L], "id")
})

test_that("an imr_data object can be fitted directly", {
  dat <- imr_data(
    simIMR$platforms,
    outcome = simIMR$outcome.binary,
    covariates = simIMR$covariates,
    type_outcome = "binary"
  )
  fit <- imr(
    dat, nu = c(-4, -3, -4), sample_mcmc = c(60, 30),
    ssize = 30, seed = 71
  )
  expect_s3_class(fit, "imr")
  expect_s3_class(fit$input_data, "imr_data")
})

test_that("the formula/data interface builds outcome and model matrix", {
  clinical <- merge(
    simIMR$outcome.binary, simIMR$covariates,
    by = "id", sort = FALSE
  )
  fit <- imr(
    y ~ age + sex + stage,
    data = clinical,
    platforms = simIMR$platforms,
    type_outcome = "binary",
    nu = c(-4, -3, -4), sample_mcmc = c(60, 30),
    ssize = 30, seed = 72
  )
  expect_s3_class(fit, "imr")
  expect_equal(fit$covariate_names, c("age", "sex", "stage"))
  expect_true(inherits(fit$formula, "formula"))
})

test_that("the formula interface accepts the standard named formula argument", {
  clinical <- merge(
    simIMR$outcome.binary, simIMR$covariates,
    by = "id", sort = FALSE
  )
  clinical$age_group <- factor(clinical$age > stats::median(clinical$age))
  fit <- imr(
    formula = y ~ sex + age_group,
    data = clinical, platforms = simIMR$platforms,
    type_outcome = "binary", nu = c(-4, -3, -4),
    sample_mcmc = c(30, 10), ssize = 30, seed = 73
  )
  expect_s3_class(fit, "imr")
  expect_true(any(grepl("^age_group", fit$covariate_names)))
})

test_that("imr_data can carry prediction platforms and covariates", {
  pred_data <- imr_data(
    platforms = list(
      genomic = simIMR$platforms$genomic[1:12, ],
      proteomic = simIMR$platforms$proteomic[1:12, ]
    ),
    covariates = simIMR$covariates
  )
  out <- predict(fit_bin, newdata = pred_data)
  expect_equal(nrow(out[["model:011"]]), 12L)
})

test_that("prediction imr_data maps named platforms to the fitted model", {
  pred_data <- imr_data(platforms = list(
    metabolomic = simIMR$platforms$metabolomic[1:12, ]
  ), covariates = simIMR$covariates)
  out <- predict(fit_bin, newdata = pred_data)
  expect_equal(nrow(out[["model:100"]]), 12L)

  wrong <- imr_data(platforms = list(
    unknown = simIMR$platforms$proteomic[1:12, ]
  ), covariates = simIMR$covariates)
  expect_error(
    predict(fit_bin, newdata = wrong),
    "uniquely match"
  )
})

test_that("predict returns full precision", {
  out <- predict(
    fit_bin,
    newdata = list(
      simIMR$platforms$genomic[1:20, ],
      simIMR$platforms$proteomic[1:20, ]
    ),
    covariates = simIMR$covariates
  )
  values <- unlist(lapply(out, `[[`, "predict"), use.names = FALSE)
  expect_true(any(abs(values - round(values, 3L)) > 1e-8))
})
