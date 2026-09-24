root <- Sys.getenv("IMR_ROOT")
lib <- normalizePath(Sys.getenv("IMR_DEPENDENCIES"), mustWork = TRUE)
stopifnot(getRversion() == "4.4.1")
lock <- read.delim(file.path(root, "dependencies/dependencies.tsv"))
for (i in seq_len(nrow(lock))) {
  package <- lock$Package[i]
  stopifnot(as.character(packageVersion(package, lib.loc = lib)) == lock$Version[i],
            dirname(normalizePath(find.package(package))) == lib)
}
cat("Verified all", nrow(lock), "locked dependency versions and load paths.\n")
