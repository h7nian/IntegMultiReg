# Compare source annotations with the recorded namespace without editing it.
args <- commandArgs(TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("Install roxygen2 before checking generated documentation.", call. = FALSE)
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this script from the package root or supply that path.", call. = FALSE)
}

work <- tempfile("imr-documentation-")
dir.create(work)
# Keep cleanup local to a function so errors remove the temporary source too.
check_documentation <- function() {
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  for (name in c("DESCRIPTION", "NAMESPACE", "R", "man", "data")) {
    source <- file.path(root, name)
    if (file.exists(source)) {
      stopifnot(file.copy(source, work, recursive = TRUE))
    }
  }
  load_source_data <- function(path) {
    env <- roxygen2::load_source(path)
    for (file in list.files(file.path(path, "data"),
                           pattern = "[.](rda|RData)$", full.names = TRUE)) {
      load(file, envir = env)
    }
    env
  }
  roxygen2::roxygenize(work, roclets = c("namespace", "rd"),
                      load_code = load_source_data)
  entries <- function(path) {
    declarations <- as.list(parse(path, keep.source = FALSE))
    entries <- lapply(declarations, function(declaration) {
      directive <- as.character(declaration[[1L]])
      arguments <- as.list(declaration)[-1L]
      # Roxygen versions may group several imports or exports in one call.
      # Compare their bindings, independently of grouping and line wrapping.
      if (directive == "importFrom") {
        package <- as.character(arguments[[1L]])
        return(vapply(arguments[-1L], function(name) {
          paste0("importFrom(", package, ",", as.character(name), ")")
        }, character(1)))
      }
      if (directive == "export") {
        return(vapply(arguments, function(name) {
          paste0("export(", as.character(name), ")")
        }, character(1)))
      }
      paste(deparse(declaration, width.cutoff = 500L), collapse = " ")
    })
    sort(unique(unlist(entries, use.names = FALSE)))
  }
  recorded <- entries(file.path(root, "NAMESPACE"))
  generated <- entries(file.path(work, "NAMESPACE"))
  if (!identical(recorded, generated)) {
    stop(paste(c("Source annotations and NAMESPACE differ.",
                 paste("Only recorded:", setdiff(recorded, generated)),
                 paste("Only generated:", setdiff(generated, recorded))),
               collapse = "\n"), call. = FALSE)
  }
  if (any(grepl("^export\\(\\.imr_", generated))) {
    stop("An internal .imr_ helper is exported.", call. = FALSE)
  }
  cat("Documentation exports match source annotations; internal helpers remain private.\n")
}
check_documentation()
