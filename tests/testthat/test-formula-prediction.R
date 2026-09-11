formula_fixture <- function(formula, clinical = NULL, id = "id") {
  if (is.null(clinical)) {
    clinical <- merge(simIMR$outcome.continuous, simIMR$covariates,
                      by = "id", sort = FALSE)
    clinical$group <- factor(ifelse(clinical$age > 0, "older", "younger"))
  }
  platforms <- simIMR$platforms
  if (id != "id") {
    names(clinical)[names(clinical) == "id"] <- id
    platforms <- lapply(platforms, function(x) {
      names(x)[1L] <- id
      x
    })
  }
  imr(formula, data = clinical, platforms = platforms, id = id,
      outcome_type = "continuous", draws = 30, burnin = 10,
      min_subgroup_size = 30, seed = 73)
}

test_that("unsupported formula semantics fail before sampling", {
  expect_error(formula_fixture(y ~ 0 + age), "requires an intercept")
  expect_error(formula_fixture(y ~ age - 1), "requires an intercept")
  expect_error(formula_fixture(y ~ age + offset(sex)), "does not support.*offset")
})

test_that("formula predictions preserve training transformations and contrasts", {
  clinical <- merge(simIMR$outcome.continuous, simIMR$covariates,
                    by = "id", sort = FALSE)
  clinical$group <- factor(ifelse(clinical$age > 0, "older", "younger"))
  contrasts(clinical$group) <- stats::contr.sum(2)
  fit <- formula_fixture(y ~ poly(age, 2) + group * sex, clinical)
  expect_equal(fit$preprocessing$contrasts$group, contrasts(clinical$group))
  expect_identical(fit$preprocessing$xlevels$group, levels(clinical$group))

  # A subset with reversed rows and a single observed factor level must still
  # use the original polynomial basis and the original factor coding.
  test <- clinical[clinical$group == "older", , drop = FALSE]
  test <- test[rev(seq_len(nrow(test))), , drop = FALSE]
  test$group <- droplevels(test$group)
  tt <- stats::delete.response(fit$preprocessing$terms)
  mf <- stats::model.frame(tt, test, xlev = fit$preprocessing$xlevels)
  mm <- stats::model.matrix(tt, mf, contrasts.arg = fit$preprocessing$contrasts)[, -1, drop = FALSE]
  encoded <- data.frame(id = test$id, mm, check.names = FALSE)
  expect_equal(
    predict(fit, simIMR$platforms, covariates = test),
    predict(fit, simIMR$platforms, covariates = encoded)
  )
  test$group <- "unseen"
  expect_error(predict(fit, simIMR$platforms, covariates = test), "new level")
  expect_error(predict(fit, simIMR$platforms,
                       covariates = clinical[c("id", "age")]), "missing formula")
  test <- clinical
  test$age[1] <- NA_real_
  expect_error(predict(fit, simIMR$platforms, covariates = test), "missing values")
})

test_that("dot expansion excludes identifiers and custom IDs work for prediction", {
  clinical <- merge(simIMR$outcome.continuous, simIMR$covariates,
                    by = "id", sort = FALSE)
  explicit <- formula_fixture(y ~ age + sex + stage, clinical)
  dotted <- formula_fixture(y ~ ., clinical)
  expect_identical(dotted$model$covariate_names, explicit$model$covariate_names)
  expect_identical(dotted$posterior$log_posterior, explicit$posterior$log_posterior)
  custom <- formula_fixture(y ~ I(age^2) + sex, clinical, id = "subject")
  names(clinical)[1L] <- "subject"
  out <- predict(custom, simIMR$platforms, covariates = clinical)
  expect_equal(sum(vapply(out, nrow, integer(1L))), nrow(clinical))
  expect_true(all(is.finite(unlist(lapply(out, `[[`, "prediction")))))
})
