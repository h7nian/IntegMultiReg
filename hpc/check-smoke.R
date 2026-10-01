args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
chain <- normalizePath(args[1])
simulation <- normalizePath(args[2])
check <- function(path) stopifnot(file.exists(path),
  identical(readRDS(paste0(path, ".status.rds"))$status, "completed"))
check(file.path(chain, "chain-01.rds"))
d <- readRDS(file.path(chain, "chain-01-diagnostic.rds"))
stopifnot(d$control$seed == 100L, d$control$mcmc$draws == 100L,
  all(is.finite(d$posterior$log_posterior)))
job <- file.path(simulation, "configuration-02-replicate-001")
for (name in c("imr-molecular.rds", "bms-molecular.rds", "l1-cph.rds",
               "uni-cph-selection.csv"))
  check(file.path(job, name))
stopifnot(identical(readRDS(file.path(simulation, "task.rds")),
  list(kind = "simulation", configurations = 2L, replicates = 1L)))
stopifnot(capabilities("png"))
withCallingHandlers(png(file.path(chain, "graphics-smoke.png"),
  type = if (capabilities("aqua")) "quartz" else "cairo", width = 600, height = 400),
  warning = function(w)stop(conditionMessage(w)))
plot(1:3, 1:3, main = "Headless graphics smoke check")
dev.off()
stopifnot(file.info(file.path(chain, "graphics-smoke.png"))$size > 1000)
cat("Chain, generator, IMR/BMS CV and Cox baselines completed. Smoke only.\n")
