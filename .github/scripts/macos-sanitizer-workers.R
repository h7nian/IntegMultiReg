# CI-only launcher configuration: production code has no sanitizer special case.
# Keep base parallel's scheduling and worker protocol; replace only its launcher
# so each child loads ASAN before the instrumented package, despite macOS SIP.
stopifnot(identical(Sys.info()[["sysname"]], "Darwin"))
launcher <- file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts",
                      "macos-sanitizer-rscript.sh")
stopifnot(file.exists(launcher), nzchar(Sys.getenv("SANITIZER_RUNTIME")))
cluster_options <- get("defaultClusterOptions", envir = asNamespace("parallel"))
stopifnot(is.environment(cluster_options),
          exists("rscript", envir = cluster_options, inherits = FALSE),
          identical(cluster_options$rscript_args, character()))
cluster_options$rscript <- normalizePath(launcher)
# Keep child diagnostics visible to the workflow's sanitizer-log gate.
cluster_options$outfile <- ""
