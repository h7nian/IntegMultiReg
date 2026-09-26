## Reproduce the Table 1 workflow from Chekouo et al. (2017).
##
## The paper-aligned run requires the supplementary KIRC data distributed with
## the Biometrics article: the authors' C supplement uses preprocessed 778 mRNA,
## 91 miRNA and 729 methylation columns on the original E1/E2/E3/E5 groups.
## That data object is expected to
## contain a list with components:
##   platforms = list(mrna = ..., mirna = ..., methylation = ...)
##   outcome.survival or outcome = data.frame(id, time, status)
##   covariates = data.frame(id, age, sex/female, stage, grade, ...)
##
## Examples:
##   Rscript paper/prepare_biom12587_table1_data.R
##   R_LIBS=.rlib Rscript paper/reproduce_table1.R --mode reduced --smoke
##   R_LIBS=.rlib Rscript paper/reproduce_table1.R --mode paper \
##     --data paper/data/kirc_table1_full.rda --mcmc 400000,50000 \
##     --k 10 --rounds 10

suppressPackageStartupMessages({
  library(IntegMultiReg)
  library(survival)
})

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit) || hit == length(args)) default else args[hit + 1L]
}

has_flag <- function(flag) any(args == flag)

mode <- arg_value("--mode", "paper")
data_path <- arg_value("--data", "paper/data/kirc_table1_full.rda")
out_dir <- arg_value("--out-dir", "paper/figures")
out_prefix <- arg_value("--prefix", NULL)
smoke <- has_flag("--smoke")

if (smoke) {
  mcmc <- c(50L, 20L)
  k_cv <- 2L
  rounds <- 1L
} else {
  mcmc <- as.integer(strsplit(arg_value("--mcmc", "400000,50000"), ",")[[1L]])
  k_cv <- as.integer(arg_value("--k", "10"))
  rounds <- as.integer(arg_value("--rounds", "10"))
}

if (length(mcmc) != 2L || any(is.na(mcmc))) {
  stop("--mcmc must be formatted as retained,burnin, e.g. 400000,50000")
}

target_mean <- rbind(
  "CPH + C"     = c(E1 = 0.72, E2 = 0.76, E3 = 0.79, E5 = 0.71, Full = 0.71),
  "IMR + C + M" = c(E1 = 0.89, E2 = 0.82, E3 = 0.89, E5 = 0.80, Full = 0.86),
  "BMS + C + M" = c(E1 = 0.87, E2 = 0.82, E3 = 0.89, E5 = 0.79, Full = 0.85),
  "IMR + M"     = c(E1 = 0.88, E2 = 0.71, E3 = 0.84, E5 = 0.72, Full = 0.82),
  "BMS + M"     = c(E1 = 0.88, E2 = 0.67, E3 = 0.84, E5 = 0.72, Full = 0.81)
)

target_sd <- rbind(
  "CPH + C"     = c(E1 = 0.01, E2 = 0.01, E3 = 0.01, E5 = 0.01, Full = 0.01),
  "IMR + C + M" = c(E1 = 0.02, E2 = 0.04, E3 = 0.02, E5 = 0.04, Full = 0.01),
  "BMS + C + M" = c(E1 = 0.03, E2 = 0.04, E3 = 0.02, E5 = 0.04, Full = 0.01),
  "IMR + M"     = c(E1 = 0.02, E2 = 0.06, E3 = 0.02, E5 = 0.05, Full = 0.01),
  "BMS + M"     = c(E1 = 0.02, E2 = 0.06, E3 = 0.02, E5 = 0.05, Full = 0.01)
)

paper_order <- c("111", "011", "101", "001")
paper_names <- c("E1", "E2", "E3", "E5")

load_table1_data <- function(mode, path) {
  if (mode == "reduced") {
    data("kircIMR", package = "IntegMultiReg")
    message("Using reduced package example kircIMR. ",
            "This validates the workflow but is not the original Table 1 data.")
    return(kircIMR)
  }

  if (!file.exists(path)) {
    stop(
      "The original Table 1 full-data object is not available at ", path, ".\n",
      "Place the Biometrics supplementary KIRC data there, or run with ",
      "--mode reduced --smoke to validate the workflow on kircIMR."
    )
  }

  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  candidates <- c("table1_kirc", "kirc_table1", "kirc_full", "kircIMR")
  found <- candidates[candidates %in% ls(env)]
  if (!length(found)) {
    stop("Could not find one of these objects in ", path, ": ",
         paste(candidates, collapse = ", "))
  }
  get(found[1L], envir = env)
}

get_outcome <- function(dat) {
  if (!is.null(dat$outcome.survival)) dat$outcome.survival else dat$outcome
}

availability_bitstrings <- function(dat) {
  ids <- get_outcome(dat)$id
  avail <- data.frame(id = ids, check.names = FALSE)
  for (nm in names(dat$platforms)) {
    avail[[nm]] <- ids %in% dat$platforms[[nm]]$id
  }
  bits <- apply(avail[names(dat$platforms)], 1L, function(x) {
    paste(as.integer(rev(x)), collapse = "")
  })
  data.frame(id = ids, bitstring = bits, check.names = FALSE)
}

cindex_from_score <- function(time, status, score) {
  ok <- is.finite(time) & !is.na(status) & is.finite(score)
  if (sum(ok) < 3L || length(unique(status[ok])) < 2L) return(NA_real_)
  as.numeric(survival::concordance(
    survival::Surv(time[ok], status[ok]) ~ score[ok],
    reverse = TRUE)$concordance)
}

summarise_cv <- function(mat) {
  mean <- colMeans(mat, na.rm = TRUE)
  sd <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  if (nrow(mat) == 1L) sd[] <- NA_real_
  list(mean = mean, sd = sd)
}

fold_ids <- function(n, k) sample(rep(seq_len(k), length.out = n))

run_cph_cv <- function(dat, k = 10L, rounds = 10L, seed = 1L) {
  set.seed(seed)
  outcome <- get_outcome(dat)
  cov <- dat$covariates
  keep_cov <- setdiff(names(cov), "id")
  cov <- cov[match(outcome$id, cov$id), , drop = FALSE]
  av <- availability_bitstrings(dat)
  bits <- av$bitstring[match(outcome$id, av$id)]

  group_values <- matrix(NA_real_, nrow = rounds, ncol = length(paper_order))
  colnames(group_values) <- paper_names
  pooled_values <- numeric(rounds)

  for (r in seq_len(rounds)) {
    pooled_time <- pooled_status <- pooled_score <- numeric(0)
    for (j in seq_along(paper_order)) {
      bit <- paper_order[j]
      idx <- which(bits == bit)
      if (length(idx) <= k) next
      folds <- fold_ids(length(idx), k)
      group_score <- rep(NA_real_, length(idx))
      for (fold in seq_len(k)) {
        train <- idx[folds != fold]
        test <- idx[folds == fold]
        train_dat <- data.frame(outcome[train, c("time", "status")],
                                cov[train, keep_cov, drop = FALSE],
                                check.names = FALSE)
        test_dat <- cov[test, keep_cov, drop = FALSE]
        fit <- try(survival::coxph(
          survival::Surv(time, status) ~ ., data = train_dat,
          ties = "breslow"), silent = TRUE)
        if (inherits(fit, "try-error")) next
        group_score[folds == fold] <- as.numeric(stats::predict(
          fit, newdata = test_dat, type = "lp"))
      }
      group_values[r, j] <- cindex_from_score(
        outcome$time[idx], outcome$status[idx], group_score)
      pooled_time <- c(pooled_time, outcome$time[idx])
      pooled_status <- c(pooled_status, outcome$status[idx])
      pooled_score <- c(pooled_score, group_score)
    }
    pooled_values[r] <- cindex_from_score(pooled_time, pooled_status,
                                          pooled_score)
  }
  out <- cbind(group_values, Full = pooled_values)
  summarise_cv(out)
}

run_imr_cv <- function(dat, method = c("IMR", "BMS"), use_cov = TRUE,
                       mcmc = c(400000L, 50000L), k = 10L, rounds = 10L,
                       seed = 1L) {
  method <- match.arg(method)
  cov <- if (use_cov) dat$covariates else NULL
  fit <- imr(
    platform_data_list = dat$platforms,
    outcome = get_outcome(dat),
    cov = cov,
    type_outcome = "right.censored",
    method = method,
    nu = c(-4, -3, -4),
    h0 = 10000,
    hh = 0.087,
    sig_alpha_psi = c(0.001, 0.001),
    thet_alph_bet = c(40, 10),
    sample_mcmc = mcmc,
    ssize = 30,
    seed = seed)

  cv <- cv_imr(fit, k = k, rounds = rounds, method = method)
  mat <- cv$total_cindex
  keep <- c(intersect(paper_order, colnames(mat)), "all")
  mat <- mat[, keep, drop = FALSE]
  colnames(mat) <- c(paper_names[match(keep[-length(keep)], paper_order)],
                     "Full")
  summarise_cv(mat)
}

format_cell <- function(mean, sd) {
  if (is.na(sd)) sprintf("%.3f", mean) else sprintf("%.3f (%.3f)", mean, sd)
}

assemble_table <- function(results) {
  means <- do.call(rbind, lapply(results, `[[`, "mean"))
  sds <- do.call(rbind, lapply(results, `[[`, "sd"))
  out <- matrix("", nrow = nrow(means), ncol = ncol(means),
                dimnames = dimnames(means))
  for (i in seq_len(nrow(means))) {
    for (j in seq_len(ncol(means))) out[i, j] <- format_cell(means[i, j],
                                                             sds[i, j])
  }
  list(cells = out, mean = means, sd = sds)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dat <- load_table1_data(mode, data_path)

prefix <- if (!is.null(out_prefix)) {
  out_prefix
} else if (mode == "paper") {
  "table1_reproduction"
} else {
  "table1_reduced"
}
if (smoke) prefix <- paste0(prefix, "_smoke")

message("Platform dimensions:")
print(sapply(dat$platforms, dim))
message("MCMC retained/burn-in: ", paste(mcmc, collapse = "/"),
        "; CV: ", rounds, " x ", k_cv, "-fold")

results <- list()

checkpoint_results <- function(stage) {
  if (!length(results)) return(invisible(NULL))
  tab_partial <- assemble_table(results)
  payload <- list(
    mode = mode,
    smoke = smoke,
    mcmc = mcmc,
    k = k_cv,
    rounds = rounds,
    stage = stage,
    completed = names(results),
    table = tab_partial,
    target_mean = target_mean,
    target_sd = target_sd
  )
  saveRDS(payload, file = file.path(out_dir, paste0(prefix, "_partial.rds")))
  utils::write.csv(tab_partial$cells,
                   file = file.path(out_dir, paste0(prefix, "_partial.csv")))
  message("Checkpoint saved after ", stage)
}

message("Running CPH + C")
results[["CPH + C"]] <- run_cph_cv(dat, k = k_cv, rounds = rounds, seed = 1)
checkpoint_results("CPH + C")
message("Running IMR + C + M")
results[["IMR + C + M"]] <- run_imr_cv(dat, "IMR", TRUE, mcmc, k_cv, rounds, 1)
checkpoint_results("IMR + C + M")
message("Running BMS + C + M")
results[["BMS + C + M"]] <- run_imr_cv(dat, "BMS", TRUE, mcmc, k_cv, rounds, 1)
checkpoint_results("BMS + C + M")
message("Running IMR + M")
results[["IMR + M"]] <- run_imr_cv(dat, "IMR", FALSE, mcmc, k_cv, rounds, 1)
checkpoint_results("IMR + M")
message("Running BMS + M")
results[["BMS + M"]] <- run_imr_cv(dat, "BMS", FALSE, mcmc, k_cv, rounds, 1)
checkpoint_results("BMS + M")

tab <- assemble_table(results)

print(tab$cells, quote = FALSE)

saveRDS(list(mode = mode, smoke = smoke, mcmc = mcmc, k = k_cv,
             rounds = rounds, table = tab, target_mean = target_mean,
             target_sd = target_sd),
        file = file.path(out_dir, paste0(prefix, ".rds")))
utils::write.csv(tab$cells, file = file.path(out_dir, paste0(prefix, ".csv")))
message("Saved results to ", out_dir)
