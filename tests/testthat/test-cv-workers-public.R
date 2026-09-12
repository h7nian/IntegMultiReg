test_that("the public workers argument preserves all CV modes and validates inputs", {
  for (mode in c("legacy", "refit", "importance")) {
    reference <- cv_imr(fit_bin, k = 2, rounds = 1, max_models = 4, cv_method = mode)
    expect_identical(cv_imr(fit_bin, k = 2, rounds = 1, max_models = 4,
                            cv_method = mode, workers = 2L), reference)
    for (invalid in list(0, -1, 1.5, NA, Inf, TRUE, "2", c(1, 2))) {
      expect_error(cv_imr(fit_bin, cv_method = mode, workers = invalid), "workers")
    }
  }
  expect_identical(formals(cv_imr)$workers, 1L)
  expect_error(cv_imr(fit_bin, k = 2, rounds = .Machine$integer.max), "task.*index limit")
})

test_that("parallel formula refits can serialize a local custom transformation", {
  clinical <- merge(simIMR$outcome.continuous, simIMR$covariates,
                    by = "id", sort = FALSE)
  fit <- local({
    square <- function(x) x^2
    imr(y ~ square(age) + sex, data = clinical, platforms = simIMR$platforms,
         outcome_type = "continuous", draws = 12, burnin = 6,
         min_subgroup_size = 30, seed = 53)
  })
  expect_identical(cv_imr(fit, k = 2, rounds = 1, cv_method = "refit", workers = 2L),
                   cv_imr(fit, k = 2, rounds = 1, cv_method = "refit"))
})

test_that("worker warnings are relayed in task order without changing results", {
  warnings <- character()
  connections <- showConnections(all = TRUE)
  result <- withCallingHandlers(IntegMultiReg:::.imr_cv_map(as.list(1:3), function(task) {
    warning(paste("task", task), call. = FALSE)
    task
  }, 2L), warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  })
  expect_identical(result, as.list(1:3))
  expect_identical(warnings, paste("task", 1:3))
  expect_identical(showConnections(all = TRUE), connections)
})

test_that("a deterministic worker startup failure restores parent process state", {
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else
    Sys.setenv(`_R_CHECK_LIMIT_CORES_` = limit))
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "true")
  settings <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                           "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                           "BLIS_NUM_THREADS"), unset = NA_character_)
  connections <- showConnections(all = TRUE)
  set.seed(812)
  rng <- .Random.seed
  expect_error(IntegMultiReg:::.imr_cv_map(as.list(1:3), identity, 3L),
                "simultaneous processes")
  expect_identical(.Random.seed, rng)
  expect_identical(Sys.getenv(names(settings), unset = NA_character_), settings)
  expect_identical(showConnections(all = TRUE), connections)
})

test_that("one-draw zero-burnin fits preserve all parallel CV modes", {
  platform <- data.frame(id = 1:12, marker = sin(1:12))
  outcomes <- list(binary = data.frame(id = 1:12, y = rep(0:1, 6)),
                   continuous = data.frame(id = 1:12, y = cos(1:12)),
                   right.censored = data.frame(id = 1:12, time = 2:13,
                                                status = rep(0:1, 6)))
  for (type in names(outcomes)) {
    fit <- imr(list(assay = platform), outcomes[[type]], outcome_type = type,
               draws = 1, burnin = 0, min_subgroup_size = 0, seed = 3)
    for (mode in c("legacy", "refit", "importance")) {
      expect_identical(cv_imr(fit, k = 2, rounds = 1, cv_method = mode, workers = 2L),
                       cv_imr(fit, k = 2, rounds = 1, cv_method = mode),
                       info = paste(type, mode))
    }
  }
})

test_that("worker bootstrap finds a library not inherited through R_LIBS", {
  variable <- Sys.getenv("R_LIBS", unset = NA_character_)
  on.exit(if (is.na(variable)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = variable))
  Sys.unsetenv("R_LIBS")
  expected <- normalizePath(getNamespaceInfo(asNamespace("IntegMultiReg"), "path"))
  paths <- IntegMultiReg:::.imr_cv_map(list(1L, 2L), function(task) {
    normalizePath(getNamespaceInfo(asNamespace("IntegMultiReg"), "path"))
  }, 2L)
  expect_identical(paths, rep(list(expected), 2L))
})
