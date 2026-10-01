root <- Sys.getenv("IMR_ROOT")
lib <- normalizePath(Sys.getenv("IMR_DEPENDENCIES"), mustWork = TRUE)
stopifnot(getRversion() == "4.4.1")
lock <- read.delim(file.path(root, "dependencies/dependencies.tsv"))
for (i in seq_len(nrow(lock))) {
  package <- lock$Package[i]
  # packageVersion() normalizes '-' to '.', while the source lock stores the
  # literal DESCRIPTION version. Compare the original metadata exactly.
  actual <- packageDescription(package, lib.loc = lib, fields = "Version")
  if (!identical(actual, lock$Version[i]))
    stop(sprintf("%s: expected DESCRIPTION version %s, found %s", package,
                 lock$Version[i], actual))
  namespace <- loadNamespace(package, lib.loc = lib)
  stopifnot(dirname(normalizePath(getNamespaceInfo(namespace, "path"))) == lib)
}
cat("Verified all", nrow(lock), "locked dependency versions and load paths.\n")
