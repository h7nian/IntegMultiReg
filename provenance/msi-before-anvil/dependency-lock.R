# Export the dependency versions used for the accepted local calculations.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
ip <- installed.packages()
roots <- c('digest','coda','posterior','survival','glmnet')
deps <- tools::package_dependencies(roots, db=ip,
  which=c('Depends','Imports','LinkingTo'), recursive=TRUE)
packages <- unique(c(roots, unlist(deps)))
packages <- setdiff(packages, rownames(ip)[!is.na(ip[,'Priority']) & ip[,'Priority']=='base'])
direct <- tools::package_dependencies(packages, db=ip,
  which=c('Depends','Imports','LinkingTo'), recursive=FALSE)
ordered <- character()
while(length(setdiff(packages, ordered))) {
  ready <- setdiff(packages, ordered)
  ready <- ready[vapply(ready, function(p) all(intersect(direct[[p]],packages) %in% ordered), TRUE)]
  if(!length(ready)) stop('Dependency cycle')
  ordered <- c(ordered, ready)
}
write.table(ip[ordered,c('Package','Version')], args[1], sep='\t',
  row.names=FALSE, quote=FALSE)
