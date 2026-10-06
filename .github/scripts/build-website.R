# Build the guides from their maintained sources rather than keeping copies.
build_website <- function() {
  guides <- c("method-coverage.md" = "inst/METHOD-COVERAGE.md",
              "migration.md" = "inst/MIGRATION.md")
  if (any(file.exists(names(guides)))) {
    stop("Temporary website guide paths already exist.", call. = FALSE)
  }
  on.exit(unlink(names(guides)), add = TRUE)
  for (destination in names(guides)) {
    stopifnot(file.copy(guides[[destination]], destination))
  }
  pkgdown::build_site_github_pages(install = FALSE, new_process = FALSE,
                                  quiet = FALSE)
  required <- c("index.html", "reference/index.html", "reference/cv_imr.html",
                "articles/IntegMultiReg.html", "method-coverage.html",
                "migration.html", "news/index.html", "authors.html")
  stopifnot(all(file.exists(file.path("docs", required))))
  cat("Package homepage, tutorial, reference, guides and citation pages built.\n")
}
build_website()
