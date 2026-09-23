## Generates the `kircIMR` real-data example shipped with IntegMultiReg.
##
## The data are a reduced TCGA-KIRC multi-platform survival example aligned
## with Chekouo, Stingo, Doecke and Do (2017, Biometrics): mRNA expression,
## miRNA expression, DNA methylation, clinical covariates and right-censored
## survival outcome.  The complete source files are public UCSC Xena TCGA data,
## not controlled-access TCGA/GDC files.  The shipped data use package-internal
## patient identifiers rather than TCGA barcodes.
##
## Run from the package root:
##   Rscript data-raw/make_kircIMR.R
##
## Note: the full HM450 methylation file is large (~423 MB compressed).  The
## script streams it and keeps only a small screened feature panel; the package
## data object itself is small.

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required to build kircIMR.")
}

options(stringsAsFactors = FALSE)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop("Environment variable ", name, " must be an integer.")
  out
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") return(default)
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) stop("Environment variable ", name, " must be numeric.")
  out
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") default else value
}

cfg <- list(
  cache_dir = env_chr("IMR_CACHE_DIR", file.path("data-raw", "cache")),
  n_mrna = env_int("IMR_N_MRNA", 50L),
  n_mirna = env_int("IMR_N_MIRNA", 30L),
  n_methylation = env_int("IMR_N_METHYLATION", 50L),
  methylation_candidates = env_int("IMR_METHYLATION_CANDIDATES", 1000L),
  max_missing = env_num("IMR_MAX_MISSING", 0.20),
  analysis_ssize = env_int("IMR_ANALYSIS_SSIZE", 30L),
  output_file = env_chr("IMR_OUTPUT_FILE",
                        file.path("data", "kircIMR.rda"))
)

urls <- c(
  clinical = "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/KIRC_clinicalMatrix",
  mrna = "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HiSeqV2.gz",
  mirna = "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/miRNA_HiSeq_gene.gz",
  methylation = "https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/HumanMethylation450.gz"
)

dir.create(cfg$cache_dir, showWarnings = FALSE, recursive = TRUE)

paths <- c(
  clinical = file.path(cfg$cache_dir, "KIRC_clinicalMatrix.tsv"),
  mrna = file.path(cfg$cache_dir, "KIRC_HiSeqV2.tsv.gz"),
  mirna = file.path(cfg$cache_dir, "KIRC_miRNA_HiSeq_gene.tsv.gz"),
  methylation = file.path(cfg$cache_dir, "KIRC_HumanMethylation450.tsv.gz")
)

download_if_missing <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 0) {
    message("Using cached file: ", path)
    return(invisible(path))
  }
  message("Downloading ", url)
  utils::download.file(url, path, mode = "wb", quiet = FALSE)
  invisible(path)
}

tcga_patient_id <- function(x) substr(x, 1L, 12L)

is_primary_tumor_sample <- function(x) {
  nchar(x) >= 15L & substr(x, 14L, 15L) == "01"
}

as_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out
}

stage_to_num <- function(x) {
  x <- toupper(trimws(x))
  out <- rep(NA_real_, length(x))
  out[grepl("STAGE I([^I]|$)|STAGE IA|STAGE IB|^I$", x)] <- 1
  out[grepl("STAGE II([^I]|$)|STAGE IIA|STAGE IIB|^II$", x)] <- 2
  out[grepl("STAGE III|^III$", x)] <- 3
  out[grepl("STAGE IV|^IV$", x)] <- 4
  out
}

grade_to_num <- function(x) {
  x <- toupper(trimws(x))
  out <- suppressWarnings(as.numeric(sub("^G", "", x)))
  out[!out %in% 1:4] <- NA_real_
  out
}

make_clinical <- function(path) {
  clin <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  clin <- clin[clin$sample_type == "Primary Tumor", , drop = FALSE]

  patient <- clin$bcr_patient_barcode
  status <- ifelse(toupper(clin$vital_status) == "DECEASED", 1L, 0L)
  death <- as_num(clin$days_to_death)
  followup <- as_num(clin$days_to_last_followup)
  time <- ifelse(status == 1L & !is.na(death), death, followup)

  covariates <- data.frame(
    id = patient,
    age = as_num(clin$age_at_initial_pathologic_diagnosis),
    female = as.numeric(toupper(clin$gender) == "FEMALE"),
    stage = stage_to_num(clin$pathologic_stage),
    grade = grade_to_num(clin$neoplasm_histologic_grade),
    check.names = FALSE
  )
  outcome <- data.frame(id = patient, time = time, status = status,
                        check.names = FALSE)

  keep <- !duplicated(patient) &
    !is.na(outcome$time) & outcome$time > 0 &
    stats::complete.cases(covariates)

  covariates <- covariates[keep, , drop = FALSE]
  outcome <- outcome[keep, , drop = FALSE]
  covariates <- covariates[order(covariates$id), , drop = FALSE]
  outcome <- outcome[match(covariates$id, outcome$id), , drop = FALSE]

  list(outcome = outcome, covariates = covariates)
}

collapse_duplicate_patients <- function(mat, samples) {
  patient <- tcga_patient_id(samples)
  keep <- !duplicated(patient)
  mat <- mat[, keep, drop = FALSE]
  colnames(mat) <- patient[keep]
  mat
}

row_missing <- function(mat) rowMeans(is.na(mat))

row_var <- function(mat) {
  apply(mat, 1L, function(x) stats::var(x, na.rm = TRUE))
}

impute_by_row_median <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    miss <- is.na(mat[i, ])
    if (any(miss)) {
      med <- stats::median(mat[i, ], na.rm = TRUE)
      if (is.na(med)) med <- 0
      mat[i, miss] <- med
    }
  }
  mat
}

cox_scores <- function(mat, outcome) {
  common <- intersect(colnames(mat), outcome$id)
  y <- outcome[match(common, outcome$id), , drop = FALSE]
  surv_y <- survival::Surv(y$time, y$status)
  mat <- mat[, common, drop = FALSE]

  score_one <- function(x) {
    if (length(unique(x[!is.na(x)])) < 3L) return(NA_real_)
    fit <- try(survival::coxph(surv_y ~ x, ties = "breslow"),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    p <- try(summary(fit)$coefficients[1L, "Pr(>|z|)"], silent = TRUE)
    if (inherits(p, "try-error") || is.na(p) || p <= 0) return(NA_real_)
    -log10(p)
  }

  scores <- apply(mat, 1L, score_one)
  scores[is.na(scores)] <- -Inf
  scores
}

sanitize_feature_names <- function(x, prefix) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  empty <- is.na(x) | x == ""
  x[empty] <- paste0(prefix, seq_len(sum(empty)))
  make.unique(x)
}

matrix_to_platform <- function(mat, prefix) {
  mat <- t(mat)
  colnames(mat) <- sanitize_feature_names(colnames(mat), prefix)
  data.frame(id = rownames(mat), mat, check.names = FALSE, row.names = NULL)
}

make_internal_id_map <- function(ids, prefix = "KIRC") {
  stats::setNames(sprintf("%s%03d", prefix, seq_along(ids)), ids)
}

replace_ids <- function(dat, id_map) {
  dat$id <- unname(id_map[dat$id])
  if (anyNA(dat$id)) stop("Internal ID mapping is incomplete.")
  dat
}

screen_loaded_matrix <- function(path, outcome, n_features, prefix,
                                 max_missing = 0.20) {
  dat <- utils::read.delim(gzfile(path), check.names = FALSE,
                           stringsAsFactors = FALSE)
  feature <- dat[[1L]]
  mat <- as.matrix(dat[, -1L, drop = FALSE])
  storage.mode(mat) <- "double"

  samples <- colnames(mat)
  keep <- is_primary_tumor_sample(samples) &
    tcga_patient_id(samples) %in% outcome$id
  mat <- collapse_duplicate_patients(mat[, keep, drop = FALSE], samples[keep])
  rownames(mat) <- feature

  keep_rows <- row_missing(mat) <= max_missing
  mat <- mat[keep_rows, , drop = FALSE]

  scores <- cox_scores(mat, outcome)
  ord <- order(scores, decreasing = TRUE, na.last = NA)
  ord <- ord[seq_len(min(n_features, length(ord)))]
  selected <- impute_by_row_median(mat[ord, , drop = FALSE])

  list(
    platform = matrix_to_platform(selected, prefix),
    scores = data.frame(
      feature = rownames(selected),
      cox_score = round(scores[ord], 4),
      stringsAsFactors = FALSE
    )
  )
}

update_top <- function(features, scores, values, candidate_n) {
  score_vec <- c(values$score, scores)
  feature_vec <- c(values$feature, features)
  ord <- order(score_vec, decreasing = TRUE, na.last = NA)
  ord <- ord[seq_len(min(candidate_n, length(ord)))]
  list(feature = feature_vec[ord], score = score_vec[ord])
}

stream_methylation_candidates <- function(path, outcome, candidate_n,
                                          max_missing = 0.20,
                                          progress_every = 50000L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1L]]
  samples <- header[-1L]
  keep <- is_primary_tumor_sample(samples) &
    tcga_patient_id(samples) %in% outcome$id
  keep_idx_all <- which(keep) + 1L
  keep_samples <- samples[keep]
  keep_patient <- tcga_patient_id(keep_samples)
  keep_first <- !duplicated(keep_patient)
  keep_idx <- keep_idx_all[keep_first]

  top <- list(feature = character(0), score = numeric(0))
  n_read <- 0L
  repeat {
    lines <- readLines(con, n = 1000L)
    if (length(lines) == 0L) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    features <- vapply(parts, `[`, character(1), 1L)
    scores <- vapply(parts, function(z) {
      x <- suppressWarnings(as.numeric(z[keep_idx]))
      if (mean(is.na(x)) > max_missing) return(NA_real_)
      stats::var(x, na.rm = TRUE)
    }, numeric(1))
    scores[is.na(scores)] <- -Inf
    top <- update_top(features, scores, top, candidate_n)
    n_read <- n_read + length(lines)
    if (n_read %% progress_every < length(lines)) {
      message("Scanned methylation rows: ", n_read)
    }
  }
  top$feature
}

read_methylation_features <- function(path, features, outcome) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1L]]
  samples <- header[-1L]
  keep <- is_primary_tumor_sample(samples) &
    tcga_patient_id(samples) %in% outcome$id
  keep_idx_all <- which(keep) + 1L
  keep_samples <- samples[keep]
  keep_patient <- tcga_patient_id(keep_samples)
  keep_first <- !duplicated(keep_patient)
  keep_idx <- keep_idx_all[keep_first]
  keep_patient <- keep_patient[keep_first]

  wanted <- setNames(seq_along(features), features)
  mat <- matrix(NA_real_, nrow = length(features), ncol = length(keep_idx))
  rownames(mat) <- features
  colnames(mat) <- keep_patient

  repeat {
    lines <- readLines(con, n = 5000L)
    if (length(lines) == 0L) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    got <- vapply(parts, `[`, character(1), 1L)
    hit <- got %in% features
    if (any(hit)) {
      for (z in parts[hit]) {
        mat[wanted[[z[1L]]], ] <- suppressWarnings(as.numeric(z[keep_idx]))
      }
    }
    if (all(rowSums(!is.na(mat)) > 0L)) break
  }
  mat
}

screen_methylation <- function(path, outcome, n_features, candidate_n,
                               max_missing = 0.20) {
  cand_cache <- file.path(
    dirname(path),
    sprintf("KIRC_methylation_candidates_n%d_miss%s.rds",
            candidate_n, gsub("[^0-9]+", "", as.character(max_missing)))
  )
  if (file.exists(cand_cache)) {
    candidates <- readRDS(cand_cache)
    message("Using cached methylation candidates: ", cand_cache)
  } else {
    message("Streaming HM450 methylation to select high-variance candidates.")
    candidates <- stream_methylation_candidates(
      path, outcome, candidate_n = candidate_n, max_missing = max_missing)
    saveRDS(candidates, cand_cache)
  }

  mat <- read_methylation_features(path, candidates, outcome)
  scores <- cox_scores(mat, outcome)
  ord <- order(scores, decreasing = TRUE, na.last = NA)
  ord <- ord[seq_len(min(n_features, length(ord)))]
  selected <- impute_by_row_median(mat[ord, , drop = FALSE])

  list(
    platform = matrix_to_platform(selected, "cg"),
    scores = data.frame(
      feature = rownames(selected),
      cox_score = round(scores[ord], 4),
      stringsAsFactors = FALSE
    )
  )
}

message("Downloading public TCGA-KIRC source files if needed.")
invisible(Map(download_if_missing, urls, paths))

message("Preparing clinical covariates and survival outcome.")
clin <- make_clinical(paths["clinical"])

message("Screening mRNA expression features.")
mrna <- screen_loaded_matrix(
  paths["mrna"], clin$outcome, n_features = cfg$n_mrna,
  prefix = "G", max_missing = cfg$max_missing)

message("Screening miRNA expression features.")
mirna <- screen_loaded_matrix(
  paths["mirna"], clin$outcome, n_features = cfg$n_mirna,
  prefix = "miR", max_missing = cfg$max_missing)

message("Screening DNA methylation features.")
methylation <- screen_methylation(
  paths["methylation"], clin$outcome, n_features = cfg$n_methylation,
  candidate_n = cfg$methylation_candidates, max_missing = cfg$max_missing)

platforms <- list(
  mrna = mrna$platform,
  mirna = mirna$platform,
  methylation = methylation$platform
)

all_platform_ids <- unique(unlist(lapply(platforms, `[[`, "id")))
outcome <- clin$outcome[clin$outcome$id %in% all_platform_ids, , drop = FALSE]
covariates <- clin$covariates[clin$covariates$id %in% all_platform_ids, ,
                              drop = FALSE]
outcome <- outcome[match(covariates$id, outcome$id), , drop = FALSE]

id_map <- make_internal_id_map(covariates$id)
platforms <- lapply(platforms, replace_ids, id_map = id_map)
covariates <- replace_ids(covariates, id_map)
outcome <- replace_ids(outcome, id_map)

availability <- data.frame(id = covariates$id, check.names = FALSE)
for (nm in names(platforms)) {
  availability[[nm]] <- availability$id %in% platforms[[nm]]$id
}
bitstrings <- apply(availability[names(platforms)], 1L, function(x) {
  paste(as.integer(rev(x)), collapse = "")
})
subgroup_sizes <- sort(table(bitstrings))
model_subgroup_sizes <- subgroup_sizes[subgroup_sizes > cfg$analysis_ssize]

kircIMR <- list(
  platforms = platforms,
  covariates = covariates,
  outcome = outcome,
  outcome.survival = outcome,
  feature_screening = list(
    mrna = mrna$scores,
    mirna = mirna$scores,
    methylation = methylation$scores,
    max_missing = cfg$max_missing,
    methylation_candidates = cfg$methylation_candidates,
    method = paste(
      "mRNA and miRNA features ranked by univariable Cox score;",
      "HM450 methylation first reduced to high-variance candidates,",
      "then ranked by univariable Cox score."
    )
  ),
  platform_availability = availability,
  subgroup_sizes = subgroup_sizes,
  model_subgroup_sizes = model_subgroup_sizes,
  paper_alignment = list(
    reference_platforms = c("mRNA expression", "miRNA expression",
                            "DNA methylation"),
    reference_covariates = c("age", "sex", "pathologic stage",
                             "histologic grade"),
    reference_outcome = "right-censored survival time",
    reference_screened_features = c(mrna = 776L, mirna = 91L,
                                    methylation = 729L),
    package_screened_features = c(mrna = cfg$n_mrna, mirna = cfg$n_mirna,
                                  methylation = cfg$n_methylation),
    analysis_ssize = cfg$analysis_ssize,
    note = paste(
      "This is a reduced public TCGA-KIRC example for package examples.",
      "The package object was generated from public UCSC Xena TCGA-KIRC",
      "sampleMap files, not from controlled-access TCGA/GDC data.",
      "It mirrors the platform/covariate/outcome structure of the Biometrics",
      "case study but does not reproduce the original full-scale analysis."
    )
  ),
  source = list(
    project = "TCGA-KIRC",
    hub = "UCSC Xena TCGA KIRC sampleMap",
    urls = urls,
    data_access = paste(
      "Derived from public UCSC Xena TCGA-KIRC files.",
      "No controlled-access TCGA/GDC files are distributed in this object."
    ),
    id_policy = paste(
      "Participant ids are package-internal labels (KIRC001, KIRC002, ...).",
      "No TCGA barcode mapping is distributed with the package; users should",
      "not attempt participant re-identification or external linkage."
    ),
    citations = c(
      tcga = "The Cancer Genome Atlas Research Network.",
      gdc = "National Cancer Institute Genomic Data Commons.",
      xena = "UCSC Xena platform.",
      methodology = paste(
        "Chekouo T, Stingo FC, Doecke JD, Do K-A (2017).",
        "A Bayesian Integrative Approach for Multi-Platform Genomic Data:",
        "A Kidney Cancer Case Study. Biometrics 73(2), 615-624."
      )
    ),
    reference = paste(
      "Chekouo T, Stingo FC, Doecke JD, Do K-A (2017).",
      "A Bayesian Integrative Approach for Multi-Platform Genomic Data:",
      "A Kidney Cancer Case Study. Biometrics 73(2), 615-624."
    )
  )
)

dir.create(dirname(cfg$output_file), showWarnings = FALSE, recursive = TRUE)
save(kircIMR, file = cfg$output_file, compress = "xz")

message("Wrote ", cfg$output_file)
message("Platform dimensions:")
print(lapply(kircIMR$platforms, dim))
message("Availability subgroup sizes:")
print(kircIMR$subgroup_sizes)
message("Modelled subgroup sizes with ssize = ", cfg$analysis_ssize, ":")
print(kircIMR$model_subgroup_sizes)
