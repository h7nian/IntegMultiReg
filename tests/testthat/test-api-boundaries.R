test_that("integer controls reject nonfinite, nonscalar and overflowing values cleanly", {
  fit_args <- list(platform_data_list = simIMR$platforms,
    outcome = simIMR$outcome.binary, type_outcome = "binary",
    sample_mcmc = c(2, 0), ssize = 30)
  invalid <- list(NA_real_, NaN, Inf, -Inf, numeric(), c(1, 2), 1.5,
                  .Machine$integer.max + 1, 1e100)
  for (arg in c("ssize", "seed")) for (value in invalid) {
    args <- fit_args; args[arg] <- list(value)
    expect_error(expect_warning(do.call(imr, args), NA), paste0("`", arg, "`"))
  }
  for (arg in c("k", "rounds", "max_models")) for (value in invalid) {
    args <- list(object = fit_bin); args[arg] <- list(value)
    expect_error(expect_warning(do.call(cv_imr, args), NA), paste0("`", arg, "`"))
  }
  for (value in list(c(1e100, 0), c(.Machine$integer.max, 1), c(1, NA))) {
    args <- fit_args; args$sample_mcmc <- value
    expect_error(expect_warning(do.call(imr, args), NA), "sample_mcmc")
  }
})

test_that("ambiguous frame names are rejected before data can be dropped", {
  for (nm in list(c("id", "x", "x"), c("id", "", "y"),
                  c("id", NA, "y"), c("id", "id", "y"))) {
    x <- data.frame(id = 1:4, x = 1:4, y = 4:1); names(x) <- nm
    expect_error(imr_data(list(assay = x)), "unique column names")
  }
  x <- data.frame(subject = 1:4, id = 1:4, marker = 4:1)
  expect_error(imr_data(list(assay = x), id = "subject"), "conflicts")
  dat <- imr_data(simIMR$platforms)
  names(dat$platforms)[1] <- NA_character_
  expect_error(validate_imr_data(dat), "complete, unique names")
  dat <- imr_data(simIMR$platforms)
  names(dat$platforms[[1]])[3] <- names(dat$platforms[[1]])[2]
  expect_error(validate_imr_data(dat), "unique column names")
})

test_that("saved fit validation catches corrupted counts, mappings and traces", {
  for (field in c("n_platform", "sample_mcmc")) {
    bad <- fit_bin; bad[[field]][1] <- 1e100
    expect_error(expect_warning(validate_imr(bad), NA), field)
  }
  for (indices in list(0L, 100L, c(1L, 1L), NA_integer_)) {
    bad <- fit_bin; bad$platform_models[[1]] <- indices
    expect_error(validate_imr(bad), "platform_models")
  }
  bad <- fit_bin; bad$log_posterior <- bad$log_posterior[-1]
  expect_error(validate_imr(bad), "one entry per iteration")
})

test_that("one retained draw and zero burn-in remain supported", {
  x <- data.frame(id = 1:12, marker = sin(1:12))
  for (type in c("binary", "continuous", "right.censored")) {
    y <- switch(type, binary = data.frame(id = x$id, y = rep(0:1, 6)),
      continuous = data.frame(id = x$id, y = cos(x$id)),
      right.censored = data.frame(id = x$id, time = x$id + 1, status = rep(0:1, 6)))
    f <- imr(list(assay = x), y, type_outcome = type, ssize = 0,
             sample_mcmc = c(1, 0), seed = 3)
    expect_true(validate_imr(f))
    expect_length(f$gam_sample, 1)
    expect_length(f$log_posterior, 1)
    p <- predict(f, list(assay = x))[[1]]$predict
    expect_true(all(is.finite(p)))
    if (type == "binary") expect_true(all(p >= 0 & p <= 1))
  }
})

test_that("undefined validation scores do not report artificial performance", {
  accuracy <- IntegMultiReg:::.imr_cv_accuracy
  expect_true(is.na(accuracy("binary", c(.2, .8), data.frame(id = 1:2, y = c(1, 1)))))
  expect_true(is.na(accuracy("right.censored", 1:2,
    data.frame(id = 1:2, time = 1:2, status = c(0, 0)))))
  expect_true(is.na(accuracy("continuous", numeric(), data.frame(id = integer(), y = numeric()))))
  expect_error(cv_imr(fit_bin, k = min(fit_bin$sample_size) + 1), "sample size")
  legacy <- fit_bin; legacy$input_data <- NULL
  expect_error(cv_imr(legacy), "raw inputs")
})

test_that("data summary printing preserves the summary and reports actual dimensions", {
  dat <- imr_data(simIMR$platforms)
  s <- summary(dat)
  printed <- capture.output(result <- withVisible(print(s)))
  expect_identical(result$value, s)
  expect_false(result$visible)
  expect_true(any(grepl(sprintf("Subjects: %d; platforms: %d", nrow(as.data.frame(dat)),
    length(simIMR$platforms)), printed, fixed = TRUE)))
  expect_identical(s$feature_counts, vapply(simIMR$platforms, ncol, integer(1)) - 1L)
})
