# R bytecode databases and native shared libraries of every locked package.
script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
root <- dirname(dirname(script)); lib <- normalizePath(Sys.getenv('R_LIBS_USER'))
lock <- read.delim(file.path(root,'dependencies/dependencies.tsv'))
files <- unlist(lapply(c(lock$Package,'IntegMultiReg'),function(p)
  list.files(file.path(lib,p),pattern='[.](rdb|rdx|so|dylib|dll)$',recursive=TRUE,full.names=TRUE)))
stopifnot(length(files)>0L, !file.exists(file.path(lib,'runtime-hashes.rds')))
saveRDS(tools::md5sum(sort(files)),file.path(lib,'runtime-hashes.rds'))
