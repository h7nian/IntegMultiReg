test_that("print returns its argument invisibly and prints key fields", {
  out <- capture.output(res <- print(fit_bin))
  expect_identical(res, fit_bin)
  expect_true(any(grepl("Outcome type", out)))
  expect_true(any(grepl("Platform key", out)))
  expect_true(any(grepl("P1 = genomic", out, fixed = TRUE)))
  expect_true(any(grepl("Availability subgroups modelled", out)))
  expect_true(any(grepl("011 : P1 + P2", out, fixed = TRUE)))
})

test_that("summary produces a per-platform selection table", {
  s <- summary(fit_bin, threshold = 0.5)
  expect_s3_class(s, "summary.imr")
  expect_named(s$selected, c("genomic", "proteomic", "metabolomic"))
  for (df in s$selected) {
    expect_true(all(c("feature", "max_mpip", "subgroup") %in% names(df)))
    # Selection uses the unrounded max mPIP, but the reported max_mpip is
    # rounded to 3 dp, so a feature just above the threshold displays as 0.500.
    if (nrow(df) > 0) expect_true(all(df$max_mpip >= 0.5))
  }
  expect_output(print(s), "summary")
})

test_that("standard S3 generics dispatch without duplicate wrappers", {
  expect_s3_class(summary(fit_bin, threshold = 0.5), "summary.imr")
  expect_named(inclusion_probabilities(fit_bin), fit_bin$model$platform_names)
  expect_false(any(c("summary_imr", "coef_imr", "plot_imr", "predict_imr") %in%
                   getNamespaceExports("IntegMultiReg")))
})

test_that("inclusion_probabilities returns named per-platform mPIP matrices", {
  mp <- inclusion_probabilities(fit_bin)
  expect_type(mp, "list")
  expect_named(mp, c("genomic", "proteomic", "metabolomic"))
  expect_equal(rownames(mp$genomic),
               fit_bin$model$subgroup_names[fit_bin$model$platform_subgroups[[1]]])
})

test_that("plot methods run without error for each type", {
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  expect_silent(plot(fit_bin, type = "selection"))
  expect_silent(plot(fit_bin, type = "theta"))
  expect_silent(plot(fit_bin, type = "trace"))
  expect_silent(plot(fit_bin, type = "selection", platform = 1))
})


test_that("ranking displays small probabilities while threshold summaries filter them", {
  # Isolate the display contract from Monte Carlo variation in the fixture.
  display <- fit_bin
  for (platform in seq_along(display$posterior$inclusion_probabilities))
    display$posterior$inclusion_probabilities[[platform]][] <- 0
  display$posterior$inclusion_probabilities[[1L]][1L, 1:3] <- c(.9, .029, .024)
  output <- capture.output(print(display, rank = TRUE, top = 3))
  expect_true(any(grepl("0.029", output, fixed = TRUE)))
  expect_true(any(grepl("0.024", output, fixed = TRUE)))
  selected <- summary(display, threshold = .5)$selected
  expect_identical(selected[[1L]]$feature,
                   display$model$feature_names[[1L]][1L])
  expect_true(all(vapply(selected[-1L], nrow, integer(1L)) == 0L))
})
