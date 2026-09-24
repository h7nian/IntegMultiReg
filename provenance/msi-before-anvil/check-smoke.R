root <- normalizePath(commandArgs(TRUE)[1])
check <- function(path) stopifnot(file.exists(path),
  identical(readRDS(paste0(path,'.status.rds'))$status, 'completed'))
check(file.path(root,'chain/chain-01.rds'))
d <- readRDS(file.path(root,'chain/chain-01-diagnostic.rds'))
stopifnot(d$control$seed == 100L, d$control$mcmc$draws == 100L,
  all(is.finite(d$posterior$log_posterior)))
job <- file.path(root,'simulation/configuration-02-replicate-001')
for(name in c('imr-molecular.rds','bms-molecular.rds','l1-cph.rds','uni-cph-selection.csv'))
  check(file.path(job,name))
stopifnot(identical(readRDS(file.path(root,'simulation/task.rds')),
  list(kind='simulation', configurations=2L, replicates=1L)))
stopifnot(capabilities('png'))
withCallingHandlers(png(file.path(root,'graphics-smoke.png'),
  type=if(capabilities('aqua')) 'quartz' else 'cairo',width=600,height=400),
  warning=function(w)stop(conditionMessage(w)))
plot(1:3,1:3,main='MSI headless graphics smoke check')
dev.off()
stopifnot(file.info(file.path(root,'graphics-smoke.png'))$size>1000)
cat('Chain, generator, IMR/BMS CV and Cox baselines completed. Smoke only.\n')
