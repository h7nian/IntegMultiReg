# Run only against an isolated library built with cv-row-fault-injection.patch.
# All fixtures are synthetic. No fault control is compiled into release code.
library(IntegMultiReg)
options(warn = 2)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  source(file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts",
                   "macos-sanitizer-workers.R"))
}

verify_failures <- function() {
  original <- Sys.getenv("IMR_TEST_ROW_FAIL_AT", unset = NA_character_)
  on.exit(if (is.na(original)) Sys.unsetenv("IMR_TEST_ROW_FAIL_AT") else
    Sys.setenv(IMR_TEST_ROW_FAIL_AT = original), add = TRUE)
  socket_count <- function() {
    connections <- showConnections(all = TRUE)
    sum(connections[, "class"] == "sockconn")
  }
  platform <- data.frame(id = 1:12, marker = sin(1:12))
  outcomes <- list(continuous = data.frame(id = 1:12, y = cos(1:12)),
                   binary = data.frame(id = 1:12, y = rep(0:1, 6)),
                   right.censored = data.frame(id = 1:12, time = 2:13,
                                                status = rep(0:1, 6)))
  count <- 0L
  for (outcome_type in names(outcomes)) {
    for (method in c("imr", "bms")) {
      Sys.unsetenv("IMR_TEST_ROW_FAIL_AT")
      fit <- imr(list(assay = platform), outcomes[[outcome_type]],
        outcome_type = outcome_type, method = method, draws = 8, burnin = 2,
        min_subgroup_size = 0, seed = 3)
      # Two guaranteed distinct, valid states. Reverse traversal starts A/A/B;
      # failure at allocation 3 occurs after a live alias in importance mode.
      first <- fit$posterior$selection_draws[[1L]]
      first[[1L]][] <- 0L
      second <- first
      second[[1L]][1L] <- 1L
      fit$posterior$selection_draws <- rep(list(first), 8L)
      fit$posterior$selection_draws[6L] <- list(second)
      for (cv_method in c("legacy", "importance")) {
        for (workers in 1:2) {
          run <- function() cv_imr(fit, k = 2L, rounds = 1L,
            cv_method = cv_method, workers = workers)
          Sys.unsetenv("IMR_TEST_ROW_FAIL_AT")
          set.seed(739)
          rng <- .Random.seed
          before <- socket_count()
          reference <- run()
          stopifnot(identical(.Random.seed, rng), identical(socket_count(), before))
          for (fail_at in 1:3) {
            Sys.setenv(IMR_TEST_ROW_FAIL_AT = fail_at)
            failure <- tryCatch(run(), error = identity)
            stopifnot(inherits(failure, "error"),
              grepl("Unable to allocate CV predictions \\(round [0-9]+, fold [0-9]+\\)",
                    conditionMessage(failure)),
              identical(.Random.seed, rng), identical(socket_count(), before))
            gc()
            Sys.unsetenv("IMR_TEST_ROW_FAIL_AT")
            stopifnot(identical(run(), reference), identical(.Random.seed, rng),
                      identical(socket_count(), before))
            count <- count + 1L
            cat(outcome_type, method, cv_method, workers, fail_at,
                "failure/recovery EXACT\n")
            flush.console()
          }
        }
      }
    }
  }
  stopifnot(count == 72L)
  cat("Verified:", count, "prediction allocation failure/recovery cases.\n")
}
verify_failures()
