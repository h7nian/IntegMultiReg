## Reproduces every figure and number reported in the IntegMultiReg JSS paper.
## Run from the repository root:  Rscript paper/make_analysis.R
suppressMessages(library(IntegMultiReg))
data("simIMR", package = "IntegMultiReg")
data("kircIMR", package = "IntegMultiReg")

figdir <- "paper/figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

FIG_MCMC <- c(8000, 2000)  # (retained draws, burn-in) for reported figures
CV_MCMC  <- c(4000, 1000)  # kept shorter for repeated prediction summaries
NU       <- c(-4, -3, -4)

## ---------------------------------------------------------------------------
## 1. Summarise the reduced TCGA-KIRC real-data example
## ---------------------------------------------------------------------------
kirc_dims <- sapply(kircIMR$platforms, dim)
rownames(kirc_dims) <- c("subjects", "variables")
cat("\n=== reduced TCGA-KIRC example ===\n")
print(kirc_dims)
cat("\n--- KIRC subgroup sizes modelled with ssize =",
    kircIMR$paper_alignment$analysis_ssize, "---\n")
print(kircIMR$model_subgroup_sizes)

## ---------------------------------------------------------------------------
## 2. Fit the right-censored model used for the illustrative figures
## ---------------------------------------------------------------------------
fit <- imr(
  platform_data_list = simIMR$platforms, outcome = simIMR$outcome.survival,
  cov = simIMR$covariates, type_outcome = "right.censored",
  nu = NU, sample_mcmc = FIG_MCMC, ssize = 30, seed = 1)

print(fit)
cat("\n--- selected biomarkers (mPIP > 0.5) ---\n")
print(summary(fit))

## Figure: subjects per availability subgroup (data overview)
pdf(file.path(figdir, "subgroups.pdf"), width = 6, height = 3.6)
plot_subgroup_sizes(fit)
dev.off()

## Figure: mPIP heatmaps for the three platforms
pdf(file.path(figdir, "selection.pdf"), width = 9, height = 3.7)
plot(fit, type = "selection")
dev.off()

## Figure: ranked bar chart of the top selected biomarkers across platforms
pdf(file.path(figdir, "top_features.pdf"), width = 6.5, height = 4)
plot_top_features(fit, top = 10)
dev.off()

## Figure: MRF interaction (theta) heatmaps
pdf(file.path(figdir, "theta.pdf"), width = 9, height = 3.2)
plot(fit, type = "theta")
dev.off()

## Figure: log-posterior trace
pdf(file.path(figdir, "trace.pdf"), width = 6, height = 3.6)
plot(fit, type = "trace")
dev.off()

## ---------------------------------------------------------------------------
## 3. Selection-recovery summary against the planted truth
## ---------------------------------------------------------------------------
truth <- simIMR$truth
mp <- coef(fit)
recovery <- data.frame(
  platform = names(mp),
  n.features = sapply(mp, ncol),
  n.true = sapply(names(truth), function(p) length(truth[[p]])),
  max_mpip_true = round(mapply(function(m, idx) max(apply(m[, idx, drop = FALSE], 2, max)),
                               mp, truth), 3),
  max_mpip_null = round(mapply(function(m, idx) {
    nullidx <- setdiff(seq_len(ncol(m)), idx)
    max(apply(m[, nullidx, drop = FALSE], 2, max))
  }, mp, truth), 3)
)
cat("\n--- selection recovery (true vs null max mPIP) ---\n")
print(recovery)

## ---------------------------------------------------------------------------
## 4. IMR vs BMS predictive comparison across outcome types
## ---------------------------------------------------------------------------
outcomes <- list(
  binary         = simIMR$outcome.binary,
  continuous     = simIMR$outcome.continuous,
  right.censored = simIMR$outcome.survival
)
metric.name <- c(binary = "AUC", continuous = "MSE", right.censored = "C-index")

compare <- function(type) {
  res <- lapply(c("IMR", "BMS"), function(meth) {
    f <- imr(
      platform_data_list = simIMR$platforms, outcome = outcomes[[type]],
      cov = simIMR$covariates, type_outcome = type, method = meth,
      nu = NU, sample_mcmc = CV_MCMC, ssize = 30, seed = 1)
    cv <- cv_imr(f, k = 5, rounds = 10)
    colMeans(cv$total_cindex)   # per-subgroup + overall
  })
  names(res) <- c("IMR", "BMS")
  res
}

set.seed(1)
comparison <- lapply(names(outcomes), compare)
names(comparison) <- names(outcomes)

cat("\n=== IMR vs BMS cross-validated accuracy (mean over 10 rounds) ===\n")
tab <- do.call(rbind, lapply(names(comparison), function(ty) {
  m <- comparison[[ty]]
  data.frame(outcome = ty, metric = unname(metric.name[ty]), method = c("IMR", "BMS"),
             rbind(round(m$IMR, 3), round(m$BMS, 3)),
             check.names = FALSE, row.names = NULL)
}))
print(tab, row.names = FALSE)

saveRDS(list(kirc_dims = kirc_dims,
             kirc_subgroups = kircIMR$model_subgroup_sizes,
             fit = fit,
             recovery = recovery,
             comparison = comparison,
             table = tab),
        file = file.path(figdir, "analysis_results.rds"))
cat("\nSaved figures and results to", figdir, "\n")
