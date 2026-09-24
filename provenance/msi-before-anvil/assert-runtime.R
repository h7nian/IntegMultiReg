script <- normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
root <- dirname(dirname(script))
lib <- normalizePath(Sys.getenv('R_LIBS_USER'), mustWork=TRUE)
stopifnot(getRversion() >= '4.4.0')
lock <- read.delim(file.path(root,'dependencies/dependencies.tsv'), stringsAsFactors=FALSE)
lock <- rbind(lock, data.frame(Package='IntegMultiReg',Version='0.2.0'))
for(i in seq_len(nrow(lock))) {
  p <- lock$Package[i]
  stopifnot(as.character(packageVersion(p, lib.loc=lib)) == lock$Version[i],
            dirname(normalizePath(find.package(p))) == lib)
}
version <- file.path(lib,'R-version.rds')
if(file.exists(version)) stopifnot(identical(readRDS(version),R.version))
gsl <- file.path(lib,'gsl-version.txt')
if(file.exists(gsl)) stopifnot(identical(readLines(gsl),system2('gsl-config','--version',stdout=TRUE)))
runtime <- file.path(lib,'runtime-hashes.rds')
if(file.exists(runtime)) {
  hashes <- readRDS(runtime)
  stopifnot(identical(tools::md5sum(names(hashes)),hashes))
}
cat('Locked package versions and runtime verified.\n')
