## Convert the Biometrics supplement data into the IntegMultiReg table format.
##
## Source:
##   Chekouo et al. (2017), Biometrics supplement code/data zip
##   Ref/biom12587-sup-0002-suppdata_code.zip
##
## Output:
##   paper/data/kirc_table1_full.rda

args <- commandArgs(trailingOnly = TRUE)

option_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) default else args[hit + 1L]
}

data_dir <- option_value(
  "--data-dir",
  "paper/data/biom12587_supplement/CcodeBiometrics/Data"
)
zip_path <- option_value(
  "--zip",
  "Ref/biom12587-sup-0002-suppdata_code.zip"
)
out_path <- option_value(
  "--out",
  "paper/data/kirc_table1_full.rda"
)

if (!dir.exists(data_dir)) {
  if (!file.exists(zip_path)) {
    stop("Could not find data_dir or zip_path: ", data_dir, " / ", zip_path)
  }
  message("Extracting supplement data from ", zip_path)
  utils::unzip(zip_path, files = grep(
    "^CcodeBiometrics/Data/",
    utils::unzip(zip_path, list = TRUE)$Name,
    value = TRUE
  ), exdir = dirname(dirname(data_dir)))
}

read_matrix <- function(file) {
  read.table(file.path(data_dir, file), header = FALSE, check.names = FALSE)
}

read_feature_names <- function(file, expected, prefix) {
  path <- file.path(data_dir, file)
  value <- readLines(path, warn = FALSE)
  if (length(value) == 1L && grepl("[[:space:]]", value)) {
    value <- strsplit(value, "[[:space:]]+")[[1L]]
    value <- value[nzchar(value)]
  }
  if (length(value) != expected) {
    warning(file, " has ", length(value), " names; expected ", expected,
            ". Falling back to generated names.")
    value <- paste0(prefix, seq_len(expected))
  }
  value
}

safe_names <- function(value, prefix) {
  value <- ifelse(nzchar(value), value, paste0(prefix, seq_along(value)))
  make.unique(make.names(value, allow_ = TRUE))
}

ids_from_group_labels <- function(platform_group, all_group, ids) {
  all_by_group <- split(seq_along(all_group), all_group)
  counters <- setNames(integer(length(all_by_group)), names(all_by_group))
  out <- character(length(platform_group))

  for (i in seq_along(platform_group)) {
    g <- as.character(platform_group[i])
    counters[g] <- counters[g] + 1L
    src <- all_by_group[[g]]
    if (is.null(src) || counters[g] > length(src)) {
      stop("Could not map platform row ", i, " for group ", g)
    }
    out[i] <- ids[src[counters[g]]]
  }

  out
}

clinical <- read_matrix("Clinical4Modelc.txt")
gene <- read_matrix("GeneExp4Modelc.txt")
mirna <- read_matrix("MirExpModel12c.txt")
methylation <- read_matrix("MethyExpModel13c.txt")
covariates <- read_matrix("PrognosticFactors.txt")

group_all <- read_matrix("Group4Modelc.txt")[[1L]]
group_mirna <- read_matrix("GroupModel12c.txt")[[1L]]
group_methylation <- read_matrix("GroupModel13c.txt")[[1L]]

stopifnot(
  nrow(clinical) == nrow(gene),
  nrow(clinical) == nrow(covariates),
  length(group_all) == nrow(clinical),
  length(group_mirna) == nrow(mirna),
  length(group_methylation) == nrow(methylation)
)

ids <- sprintf("KIRC%03d", seq_len(nrow(clinical)))

gene_names <- read_feature_names(
  "NamesGeneExp4Modelc.txt", ncol(gene), "gene"
)
mirna_names <- read_feature_names(
  "NamesMirExpModel12c.txt", ncol(mirna), "mirna"
)
methylation_names <- read_feature_names(
  "NamesMethyExpModel13c.txt", ncol(methylation), "probe"
)
covariate_names <- read_feature_names(
  "NamesPrognosticFactors.txt", ncol(covariates), "cov"
)

colnames(gene) <- safe_names(gene_names, "gene")
colnames(mirna) <- safe_names(mirna_names, "mirna")
colnames(methylation) <- safe_names(methylation_names, "probe")
colnames(covariates) <- safe_names(covariate_names, "cov")

mrna_ids <- ids
mirna_ids <- ids_from_group_labels(group_mirna, group_all, ids)
methylation_ids <- ids_from_group_labels(group_methylation, group_all, ids)

platforms <- list(
  mrna = data.frame(id = mrna_ids, gene, check.names = FALSE),
  mirna = data.frame(id = mirna_ids, mirna, check.names = FALSE),
  methylation = data.frame(id = methylation_ids, methylation,
                           check.names = FALSE)
)

covariates <- data.frame(id = ids, covariates, check.names = FALSE)

## The original C reader applies log() to the first clinical column before the
## survival sampler sees it. The package now applies that transformation itself;
## store raw positive months as the public outcome and log times only for audit.
outcome_raw <- data.frame(
  id = ids,
  time = as.numeric(clinical[[1L]]),
  status = as.integer(clinical[[2L]]),
  check.names = FALSE
)
outcome_survival <- data.frame(
  id = ids,
  time = log(outcome_raw$time),
  status = outcome_raw$status,
  check.names = FALSE
)

group_label <- c("1" = "E1", "2" = "E2", "3" = "E3", "4" = "E5")
paper_subgroup <- unname(group_label[as.character(group_all)])

platform_availability <- data.frame(
  id = ids,
  mrna = ids %in% platforms$mrna$id,
  mirna = ids %in% platforms$mirna$id,
  methylation = ids %in% platforms$methylation$id,
  paper_group_code = group_all,
  paper_subgroup = paper_subgroup,
  check.names = FALSE
)

metadata <- list(
  source_zip = zip_path,
  source_hash = tools::md5sum(zip_path),
  source_file_hashes = tools::md5sum(list.files(data_dir, pattern = "[.]txt$", full.names = TRUE)),
  source_data_dir = data_dir,
  original_readme_command = paste(
    "./mainIMR",
    "ClinOut=Clinical4Modelc.txt",
    "GeneExp=GeneExp4Modelc.txt",
    "GrpSampleGenEx=Group4Modelc.txt",
    "miRNAExp=MirExpModel12c.txt",
    "GrpSamplemiR=GroupModel12c.txt",
    "MethyExp=MethyExpModel13c.txt",
    "GrpSampleMethy=GroupModel13c.txt",
    "ClinicInfo=PrognosticFactors.txt",
    "sampleafterburnin=200000",
    "burnin=20000",
    "CVfoldRound=2",
    "TypaData=Kid",
    "Method=IMR"
  ),
  dimensions = list(
    clinical = dim(clinical),
    mrna = dim(platforms$mrna),
    mirna = dim(platforms$mirna),
    methylation = dim(platforms$methylation),
    covariates = dim(covariates)
  ),
  group_counts = table(paper_subgroup),
  event_counts = table(outcome_raw$status),
  outcome_time_scale = paste(
    "outcome.survival$time and outcome.raw$time are positive original months;",
    "fit with survival_scale=log. outcome.log is for auditing ImportData1()."
  ),
  original_feature_names = list(
    mrna = gene_names,
    mirna = mirna_names,
    methylation = methylation_names,
    covariates = covariate_names
  )
)

kirc_full <- list(
  platforms = platforms,
  covariates = covariates,
  outcome = outcome_raw,
  outcome.survival = outcome_raw,
  outcome.log = outcome_survival,
  outcome.raw = outcome_raw,
  platform_availability = platform_availability,
  subgroup_sizes = table(paper_subgroup),
  paper_alignment = metadata
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
save(kirc_full, file = out_path, compress = "xz")

message("Saved ", out_path)
message("Platform dimensions:")
print(sapply(kirc_full$platforms, dim))
message("Subgroup sizes:")
print(kirc_full$subgroup_sizes)
message("Outcome status counts, 0=censored and 1=event:")
print(table(kirc_full$outcome.raw$status))
