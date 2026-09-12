# Match the parent's execution environment for exact serial/parallel checks.
launcher <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts",
                      "valgrind-worker-rscript.sh")
stopifnot(identical(Sys.info()[["sysname"]], "Linux"), file.exists(launcher),
          dir.exists(Sys.getenv("IMR_VALGRIND_LOG_DIR")))
cluster_options <- get("defaultClusterOptions", envir = asNamespace("parallel"))
stopifnot(is.environment(cluster_options),
          exists("rscript", envir = cluster_options, inherits = FALSE),
          identical(cluster_options$rscript_args, character()))
cluster_options$rscript <- normalizePath(launcher)
cluster_options$outfile <- ""
