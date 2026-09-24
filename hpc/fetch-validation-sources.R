# Test/documentation build dependencies are separate from the scientific lock.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, !dir.exists(args[1]))
dir.create(args[1], recursive = TRUE)
repository <- "https://cloud.r-project.org"
db <- available.packages(repos = repository, type = "source")
roots <- c("testthat", "knitr", "rmarkdown", "survival")
dependencies <- tools::package_dependencies(roots, db,
  which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
packages <- intersect(unique(c(roots, unlist(dependencies))), rownames(db))
direct <- tools::package_dependencies(packages, db,
  which = c("Depends", "Imports", "LinkingTo"))
ordered <- character()
while (length(remaining <- setdiff(packages, ordered))) {
  ready <- remaining[vapply(remaining, function(package)
    all(intersect(direct[[package]], packages) %in% ordered), TRUE)]
  stopifnot(length(ready) > 0L)
  ordered <- c(ordered, ready)
}
downloads <- download.packages(ordered, destdir = args[1],
  repos = repository, available = db, type = "source")
stopifnot(setequal(downloads[, 1], ordered))
files <- downloads[match(ordered, downloads[, 1]), 2]
expected <- db[ordered, "MD5sum"]
stopifnot(all(unname(tools::md5sum(files)) == expected))
write.table(data.frame(Package = ordered, Version = db[ordered, "Version"],
  file = basename(files), md5 = expected), file.path(args[1], "dependencies.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(repository = repository, metadata = db[ordered, , drop = FALSE],
  downloaded = Sys.time()), file.path(args[1], "repository.rds"))
cat("Verified", length(ordered), "validation dependency sources.\n")
