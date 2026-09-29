test_that("integer controls reject nonfinite, nonscalar and overflowing values cleanly", {
  fit_args <- list(x = simIMR$platforms,
    outcome = simIMR$outcome.binary, outcome_type = "binary",
    draws = 2, burnin = 0, min_subgroup_size = 30)
  invalid <- list(NA_real_, NaN, Inf, -Inf, numeric(), c(1, 2), 1.5,
                  .Machine$integer.max + 1, 1e100)
  for (arg in c("min_subgroup_size", "seed")) for (value in invalid) {
    args <- fit_args; args[arg] <- list(value)
    expect_error(expect_warning(do.call(imr, args), NA), paste0("`", arg, "`"))
  }
  for (arg in c("k", "rounds", "max_models")) for (value in invalid) {
    args <- list(object = fit_bin); args[arg] <- list(value)
    expect_error(expect_warning(do.call(cv_imr, args), NA), paste0("`", arg, "`"))
  }
  for (value in list(1e100, NA_real_, 1.5, numeric(), c(1, 2))) {
    args <- fit_args; args$draws <- value
    expect_error(expect_warning(do.call(imr, args), NA), "`draws`")
  }
  args <- fit_args
  args$draws <- .Machine$integer.max
  args$burnin <- 1L
  expect_error(do.call(imr, args), "sum of `draws` and `burnin`")
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
  bad <- fit_bin; bad$model$n_platforms <- 1e100
  expect_error(expect_warning(validate_imr(bad), NA), "platform metadata")
  bad <- fit_bin; bad$control$mcmc$draws <- 1e100
  expect_error(expect_warning(validate_imr(bad), NA), "control\\$mcmc")
  for (indices in list(0L, 100L, c(1L, 1L), NA_integer_)) {
    bad <- fit_bin; bad$model$platform_subgroups[[1]] <- indices
    expect_error(validate_imr(bad), "subgroup indices|not reciprocal")
  }
  bad <- fit_bin; bad$posterior$log_posterior <- bad$posterior$log_posterior[-1]
  expect_error(validate_imr(bad), "one entry per iteration")
})

test_that("one retained draw and zero burn-in remain supported", {
  x <- data.frame(id = 1:12, marker = sin(1:12))
  for (type in c("binary", "continuous", "right.censored")) {
    y <- switch(type, binary = data.frame(id = x$id, y = rep(0:1, 6)),
      continuous = data.frame(id = x$id, y = cos(x$id)),
      right.censored = data.frame(id = x$id, time = x$id + 1, status = rep(0:1, 6)))
    f <- imr(list(assay = x), y, outcome_type = type, min_subgroup_size = 0,
             draws = 1, burnin = 0, seed = 3)
    expect_true(validate_imr(f))
    expect_length(f$posterior$selection_draws, 1)
    expect_length(f$posterior$log_posterior, 1)
    p <- predict(f, list(assay = x))[[1]]$prediction
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
  expect_error(cv_imr(cv_method = "refit", fit_bin, k = min(fit_bin$model$sample_sizes) + 1), "sample size")
  damaged <- fit_bin; damaged$preprocessing$input_data <- NULL
  expect_error(cv_imr(cv_method = "refit", damaged), "raw inputs|Preprocessed|missing `input_data`")
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

test_that("the default fit method rejects unsupported input classes", {
  expect_error(imr(1), "`x` must be a platform list, formula, or `imr_data` object",
               fixed = TRUE)
})
