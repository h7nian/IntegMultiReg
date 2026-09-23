# Compile the archived data generator independently of the package.
load_original_generator <- function(archive, out_dir, adapter) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir)
  files <- c("ReadData.c", "utils.c", "utils.h", "myheader.h")
  utils::unzip(archive, files = paste0("CcodeBiometrics/", files), exdir = out_dir)
  work <- file.path(out_dir, "CcodeBiometrics")
  hashes <- tools::md5sum(file.path(work, files))
  path <- file.path(work, "utils.c")
  writeLines(sub("#include <malloc.h>", "#include <stdlib.h>",
    readLines(path, warn = FALSE), fixed = TRUE), path)
  # The archived generator leaks its auxiliary RNG. Free it at function exit;
  # no random values or numerical operations are changed.
  path <- file.path(work, "ReadData.c")
  lines <- readLines(path, warn = FALSE)
  at <- tail(which(trimws(lines) == "}"), 1L)
  lines <- append(lines, "gsl_rng_free(r1);", after = at - 1L)
  writeLines(lines, path)
  file.copy(adapter, file.path(work, "generator.c"), overwrite = TRUE)
  old <- setwd(work); on.exit(setwd(old))
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "-o",
    paste0("generator", .Platform$dynlib.ext), "generator.c", "ReadData.c", "utils.c"),
    env = c(paste0("PKG_CPPFLAGS=", shQuote(system2("gsl-config", "--cflags", stdout = TRUE))),
            paste0("PKG_LIBS=", shQuote(system2("gsl-config", "--libs", stdout = TRUE))),
            "PKG_CFLAGS=-std=gnu11"),
    stdout = file.path(out_dir, "compile.log"), stderr = file.path(out_dir, "compile.log"))
  if (status != 0) stop("Cannot compile original generator; see compile.log")
  writeLines(c("Original source hashes before portability/resource patches:",
    paste(names(hashes), hashes), "utils.c: malloc.h -> stdlib.h",
    "ReadData.c: free auxiliary r1 RNG at generdata exit."), file.path(out_dir, "manifest.txt"))
  dyn.load(file.path(work, paste0("generator", .Platform$dynlib.ext)))[["name"]]
}

original_simulation <- function(data, dll, scenario, seed, rho = 0, half = FALSE,
                                marker_design = c("fixed", "random")) {
  marker_design <- match.arg(marker_design)
  stopifnot(scenario %in% 1:3, seed >= 0, rho >= 0, rho < 1,
            identical(vapply(data$platforms, ncol, 0L) - 1L, c(mrna = 778L, mirna = 91L, methylation = 729L)))
  bits <- c("111", "011", "101", "001")
  groups <- list(1:4, 1:2, c(1, 3))
  availability <- data$platform_availability
  ids <- lapply(c("E1", "E2", "E3", "E5"), function(g) availability$id[availability$paper_subgroup == g])
  permutation <- lapply(data$platforms, function(x) seq_len(ncol(x) - 1L))
  if (marker_design == "random") {
    # A common permutation per platform randomizes identities while preserving
    # the archived generator's exact within/across-group overlap pattern.
    # Keep this RNG separate from the archived GSL response generator.
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    permutation <- lapply(permutation, sample)
  }
  matrices <- lapply(1:3, function(p) lapply(groups[[p]], function(g) {
    value <- as.matrix(data$platforms[[p]][match(ids[[g]], data$platforms[[p]]$id), -1])
    value <- value[, permutation[[p]], drop = FALSE]
    storage.mode(value) <- "double"; value
  }))
  output <- .Call("original_generate", matrices, as.integer(scenario), as.double(seed),
                   as.double(rho), as.integer(half), PACKAGE = dll)
  platforms <- lapply(1:3, function(p) {
    value <- do.call(rbind, output[[2]][[p]])
    value <- value[, order(permutation[[p]]), drop = FALSE]
    colnames(value) <- names(data$platforms[[p]])[-1]
    data.frame(id = unlist(ids[groups[[p]]], use.names = FALSE), value, check.names = FALSE)
  })
  names(platforms) <- names(data$platforms)
  y <- do.call(rbind, output[[1]])
  outcome <- data.frame(id = unlist(ids, use.names = FALSE), time = exp(y[, 1]), status = as.integer(y[, 2]))
  truth <- output[[3]]; names(truth) <- names(platforms)
  for (p in 1:3) truth[[p]] <- truth[[p]][, order(permutation[[p]]), drop = FALSE]
  for (p in 1:3) dimnames(truth[[p]]) <- list(bits[groups[[p]]], names(platforms[[p]])[-1])
  list(platforms = platforms, outcome = outcome, covariates = NULL, truth = truth,
       standardize = FALSE, generator = list(scenario = scenario, seed = seed, rho = rho, half = half,
         marker_design = marker_design, marker_permutation = permutation,
         marker_rng = if (marker_design == "random") c("Mersenne-Twister", "Inversion", "Rejection") else NULL,
         source = "archived generdata", censoring_fraction = mean(outcome$status == 0)))
}
