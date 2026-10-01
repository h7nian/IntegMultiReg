# Plot retained diagnostic projections; never reconstruct or thin saved chains.
plot_appendix_chains <- function(run_dir, out_dir = file.path(run_dir, "figures")) {
  settings <- readRDS(file.path(run_dir, "settings.rds"))
  fits <- lapply(1:8, function(i) readRDS(file.path(run_dir, sprintf("chain-%02d-diagnostic.rds", i))))
  stopifnot(length(settings$seeds) == 8L)
  correlations <- read.csv(file.path(run_dir, "diagnostics/mpip_correlations.csv"), colClasses = c(subgroup = "character"))
  acfs <- read.csv(file.path(run_dir, "diagnostics/theta_acf.csv"), colClasses = c(subgroup1 = "character", subgroup2 = "character"))
  n <- settings$draws
   burn <- settings$burnin
  stopifnot(all(vapply(fits, function(f) length(f$posterior$log_posterior) == n + burn, logical(1))))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  colours <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
               "#555555", "#000000")
  note <- sprintf("%s | %s retained + %s burn-in per chain | seeds %s-%s",
    if (isTRUE(settings$quick)) "SHORT RUN: workflow check only" else "Diagnostic review; no automatic convergence claim",
    format(n, big.mark = ","), format(burn, big.mark = ","), min(settings$seeds), max(settings$seeds))
  device <- function(name, code) {
    png(file.path(out_dir, paste0(name, ".png")), width = 1800, height = 1500,
        res = 170)
    on.exit(dev.off())
    par(mfrow = c(4, 2), mar = c(3.4, 5.5, 2.1, 1), oma = c(3.5, 0, 3, 0),
        mgp = c(2.6, .65, 0), tcl = -.25, las = 1, family = "sans", cex = .85)
    force(code)
  }
  for (retained in c(FALSE, TRUE)) {
    keep <- if (retained) seq.int(burn + 1L, burn + n) else seq_len(burn + n)
    ys <- lapply(fits, function(f) f$posterior$log_posterior[keep])
    stopifnot(all(vapply(ys, function(y) all(is.finite(y)), logical(1))))
    ylim <- range(unlist(ys, use.names = FALSE))
    device(if (retained) "log-posterior-retained" else "log-posterior-full", {
      for (i in 1:8) {
        plot(keep, ys[[i]], type = "n", ylim = ylim, xlab = "MCMC iteration",
             ylab = "", main = sprintf("Chain %d (seed %d)", i,
                                       settings$seeds[i]))
        mtext("Log posterior", side = 2, line = 4.1, las = 0)
        if (!retained && burn > 0) {
          rect(0, par("usr")[3], burn, par("usr")[4], col = "#EEEEEE",
               border = NA)
          abline(v = burn, lty = 2, col = "#888888")
        }
        lines(keep, ys[[i]], col = colours[i], lwd = .65)
      }
      mtext(if (retained) "Log-posterior traces after burn-in" else "Full log-posterior traces (shading marks burn-in)",
            outer = TRUE, side = 3, line = 1, font = 2, cex = 1.05)
      mtext(note, outer = TRUE, side = 1, line = 1.8, cex = .8)
    })
  }
  keys <- unique(acfs[c("platform", "subgroup1", "subgroup2")])
  stopifnot(nrow(keys) == 8L)
  for (k in seq_len(nrow(keys))) {
    p <- match(keys$platform[k], fits[[1]]$model$platform_names)
    groups <- fits[[1]]$model$subgroup_names[fits[[1]]$model$platform_subgroups[[p]]]
    pairs <- do.call(cbind, lapply(2:length(groups),
                                   function(j) rbind(seq_len(j - 1L), j)))
    column <- which(groups[pairs[1, ]] == keys$subgroup1[k] &
                    groups[pairs[2, ]] == keys$subgroup2[k])
    stopifnot(length(column) == 1L)
    ys <- lapply(fits, function(f) f$posterior$interaction_draws[[p]][, column])
    stopifnot(all(vapply(ys, function(y) length(y) == n && all(is.finite(y)),
                         logical(1))))
    ylim <- range(unlist(ys, use.names = FALSE))
    name <- sprintf("theta-trace-%s-%s-%s", keys$platform[k], keys$subgroup1[k],
                    keys$subgroup2[k])
    device(name, {
      for (i in 1:8) {
        plot(seq_len(n), ys[[i]], type = "l", ylim = ylim, col = colours[i],
             lwd = .65,
             xlab = "Retained iteration", ylab = "",
             main = sprintf("Chain %d (seed %d)", i, settings$seeds[i]))
        mtext(expression(theta), side = 2, line = 4.1, las = 0)
        abline(h = mean(ys[[i]]), col = "#555555", lty = 2, lwd = .8)
      }
      mtext(sprintf("%s interaction: subgroups %s / %s", keys$platform[k],
                    keys$subgroup1[k], keys$subgroup2[k]),
            outer = TRUE, side = 3, line = 1, font = 2, cex = 1.05)
      mtext("Shared axes across chains; dashed line marks each chain mean.",
            outer = TRUE, side = 1, line = .4, cex = .75)
      mtext(note, outer = TRUE, side = 1, line = 2, cex = .8)
    })
  }
  device("theta-acf", {
    for (k in seq_len(nrow(keys))) {
      x <- subset(acfs, platform == keys$platform[k] & subgroup1 == keys$subgroup1[k] & subgroup2 == keys$subgroup2[k])
      plot(range(x$lag), c(-1, 1), type = "n",
           xlab = "Lag (retained iterations)", ylab = "Autocorrelation",
           main = sprintf("%s: %s / %s", keys$platform[k], keys$subgroup1[k],
                          keys$subgroup2[k]))
      abline(h = 0, col = "#CCCCCC")
      for (i in 1:8) {
        z <- x[x$chain == i, ]
        lines(z$lag, z$acf, col = colours[i], lty = if (i <= 4) 1 else 2, lwd = 1)
      }
    }
    mtext("Interaction-parameter autocorrelation by chain", outer = TRUE,
          side = 3, line = 1, font = 2, cex = 1.05)
    for (i in 1:8) mtext(sprintf("Chain %d (%s)", i, if (i <= 4) "solid" else "dashed"),
      outer = TRUE, side = 1, line = .4, at = seq(.06, .94, length.out = 8)[i],
      col = colours[i], cex = .65)
    mtext(note, outer = TRUE, side = 1, line = 2.4, cex = .8)
  })
  groups <- unique(correlations[c("platform", "subgroup")])
  stopifnot(nrow(groups) == 8L)
  device("mpip-correlations", {
    for (k in seq_len(nrow(groups))) {
      x <- subset(correlations, platform == groups$platform[k] & subgroup == groups$subgroup[k])
      p <- match(groups$platform[k], fits[[1]]$model$platform_names)
      subgroup_names <- fits[[1]]$model$subgroup_names[fits[[1]]$model$platform_subgroups[[p]]]
      s <- match(groups$subgroup[k], subgroup_names)
      stopifnot(!is.na(p), !is.na(s))
      diagonal <- vapply(fits, function(f) {
        v <- f$posterior$inclusion_probabilities[[p]][s, ]
        if (length(v) > 1L && sd(v) > 0) 1 else NA_real_
      }, numeric(1))
      mat <- diag(diagonal)
      mat[row(mat) != col(mat)] <- NA_real_
      for (j in seq_len(nrow(x))) mat[x$chain1[j], x$chain2[j]] <- mat[x$chain2[j], x$chain1[j]] <- x$pearson[j]
      image(1:8, 1:8, mat, zlim = c(-1, 1), col = gray.colors(101, start = 1, end = .15),
            xlab = "Chain", ylab = "Chain", axes = FALSE,
            main = sprintf("%s: subgroup %s", groups$platform[k],
                           groups$subgroup[k]))
      axis(1, at = 1:8)
       axis(2, at = 1:8, cex.axis = .65, gap.axis = 0)
       box()
      for (i in 1:8) for (j in 1:8) text(i, j,
        if (is.na(mat[i, j])) "NA" else sprintf("%.2f", mat[i, j]),
        col = if (!is.na(mat[i, j]) && mat[i, j] > .2) "white" else "black", cex = .68)
    }
    mtext("Feature-mPIP Pearson correlations between chains", outer = TRUE,
          side = 3, line = 1, font = 2, cex = 1.05)
    mtext("Shared scale: white = -1; dark = +1. Labels give values; NA denotes an undefined correlation.",
          outer = TRUE, side = 1, line = .4, cex = .75)
    mtext(note, outer = TRUE, side = 1, line = 2, cex = .8)
  })
  writeLines(c("All panels of each trace figure share x and y limits; full and post-burn-in views are separate.",
    "Every saved log-posterior draw is plotted; no numerical thinning is performed.",
    "Each interaction trace figure uses every retained draw and shares axes across all eight chains; dashed lines mark chain means.",
    "ACF panels share [-1,1] limits and chain colours; line style supplements colour.",
    "mPIP correlations share [-1,1] shading and numeric cell labels.", note), file.path(out_dir, "FIGURE-NOTES.txt"))
  invisible(out_dir)
}
