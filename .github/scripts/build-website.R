options(repos = c(CRAN = "https://cloud.r-project.org"))
stopifnot(requireNamespace("pkgdown", quietly = TRUE))
pkgdown::build_site(
  pkg = ".", examples = TRUE, preview = FALSE,
  new_process = FALSE, install = TRUE, quiet = FALSE
)
required <- c(
  "index.html", "reference/index.html", "reference/imr.html",
  "reference/posterior_draws.html", "reference/cv_imr.html",
  "articles/IntegMultiReg.html", "articles/uncertainty.html",
  "articles/interpretation.html", "news/index.html", "authors.html"
)
if (!all(file.exists(file.path("docs", required)))) {
  stop("The documentation build did not create all required pages.")
}
