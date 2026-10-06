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
  # pkgdown sees temporary Markdown names; expose the maintained source paths.
  for (temporary in names(guides)) {
    page <- file.path("docs", sub("[.]md$", ".html", temporary))
    html <- xml2::read_html(page)
    links <- xml2::xml_find_all(html, "//small[contains(@class, 'dont-index')]//a")
    href <- xml2::xml_attr(links, "href")
    source <- links[endsWith(href, paste0("/", temporary))]
    stopifnot(length(source) == 1L)
    old <- xml2::xml_attr(source, "href")
    xml2::xml_attr(source, "href") <- paste0(
      substr(old, 1L, nchar(old) - nchar(temporary)), guides[[temporary]])
    xml2::xml_text(xml2::xml_find_first(source, ".//code")) <- guides[[temporary]]
    xml2::write_html(html, page)
  }
  required <- c("index.html", "reference/index.html", "reference/cv_imr.html",
                "articles/IntegMultiReg.html", "method-coverage.html",
                "migration.html", "news/index.html", "authors.html")
  stopifnot(all(file.exists(file.path("docs", required))))
  cat("Package homepage, tutorial, reference, guides and citation pages built.\n")
}
build_website()
