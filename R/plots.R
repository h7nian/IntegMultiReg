#' Plot Method for IMR Fits
#'
#' @description
#' Visualizes a fitted `"imr"` object. Five plot types are available:
#' \describe{
#'   \item{`"selection"`}{Heatmap of the marginal posterior inclusion
#'     probabilities (mPIP), one panel per platform, with features on the
#'     horizontal axis, availability subgroups on the vertical axis and a small
#'     intensity legend showing that darker values are closer to 1.}
#'   \item{`"theta"`}{Heatmap of the posterior mean MRF interaction parameters
#'     between availability subgroups, one panel per platform.}
#'   \item{`"trace"`}{Trace plot of the log-posterior across MCMC iterations.}
#'   \item{`"theta_trace"`}{Trace plot for one retained MRF interaction
#'     parameter, selected by `platform` and `parameter`.}
#'   \item{`"selection_trace"`}{Trace plot for one retained selection
#'     indicator, selected by `platform`, `subgroup` and `feature`.}
#' }
#' See [plot_top_features()] and [plot_subgroup_sizes()] for two further ready
#' made displays.
#'
#' @section Statistical interpretation:
#' The selection heatmap displays \eqn{\widehat\pi_{lsj}}{mPIP_lsj}, the retained-chain
#' indicator average defined in [coef.imr()]. The theta heatmap displays
#' \eqn{B^{-1}\sum_{b=1}^{B}\theta_{l,sh}^{(b)}}{mean_b theta_l,sh[b]}. An MRF interaction measures
#' prior coupling of selection indicators; it is not a correlation between
#' measured biomarkers or a regression effect.
#'
#' `theta_trace` and `selection_trace` display retained parameter draws in
#' iteration order. The log-posterior trace includes burn-in and reflects the
#' fitted sampler's density convention. These plots can reveal mixing or drift,
#' but a flat trace alone does not establish convergence. They display the
#' original selection stage; the separate conditional coefficient-stage
#' diagnostic is defined in [posterior_draws()].
#'
#' @param x A fitted object of class `"imr"`.
#' @param type Character; one of `"selection"` (default), `"theta"`,
#'   `"trace"`, `"theta_trace"` or `"selection_trace"`.
#' @param platform Optional integer vector selecting which platforms to display
#'   for the `"selection"` and `"theta"` plots; defaults to all platforms.
#' @param parameter Positive integer selecting a theta-pair column for
#'   `type = "theta_trace"`.
#' @param subgroup Positive integer selecting a platform-specific subgroup row
#'   for `type = "selection_trace"`.
#' @param feature Positive integer selecting a feature column for
#'   `type = "selection_trace"`.
#' @param base_cex Overall text-size multiplier. The `cex_*` arguments default
#'   to values derived from this multiplier (default `1`).
#' @param cex_axis Axis-label size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_lab Axis-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_main Main-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param col Optional colours. For heatmaps this is the colour scale; for
#'   trace plots this is the line colour.
#' @param palette Optional heatmap palette name used when `col = NULL`.
#'   Selection plots default to `"grey"` to preserve the mPIP intensity scale;
#'   theta plots default to the muted grey-blue `"heatmap"` palette.
#' @param legend Logical; for `"selection"` plots, should the mPIP intensity
#'   legend be drawn (default `TRUE`)?
#' @param legend_width Relative width of the intensity-legend panel for
#'   `"selection"` plots (default `0.28`).
#' @param mar,mgp Optional graphical margin and axis-title placement vectors
#'   passed to [graphics::par()] for finer layout control.
#' @param ... Further graphical parameters passed to the underlying plotting
#'   functions.
#' @return `NULL`, invisibly; called for the side effect of producing a plot.
#' @seealso [imr()], [plot_top_features()], [plot_subgroup_sizes()]
#' @export
plot.imr <- function(x, type = c("selection", "theta", "trace",
                                 "theta_trace", "selection_trace"),
                     platform = NULL, base_cex = 1, cex_axis = NULL,
                     cex_lab = NULL, cex_main = NULL, col = NULL,
                     palette = NULL,
                     legend = TRUE, legend_width = 0.28,
                     mar = NULL, mgp = NULL, parameter = 1L,
                     subgroup = 1L, feature = 1L, ...) {
  if (!inherits(x, "imr")) {
    .imr_abort("`x` must be an `imr` object returned by `imr()`.")
  }
  type <- match.arg(type)
  dots <- list(...)
  cex_defaults <- if (type == "selection") {
    list(axis = 1.05, lab = 1.3, main = 1.4,
         names = 1, legend = 1, values = 1)
  } else {
    list(axis = 0.95, lab = 1.1, main = 1.15,
         names = 1, legend = 1, values = 1)
  }
  sz <- .imr_plot_cex(
    dots, base_cex = base_cex, cex_axis = cex_axis,
    cex_lab = cex_lab, cex_main = cex_main,
    defaults = cex_defaults
  )
  cex_axis <- sz$axis
  cex_lab <- sz$lab
  cex_main <- sz$main
  dots <- sz$dots
  .imr_check_flag(legend, "legend")
  legend_width <- .imr_check_numeric_vector(
    legend_width, "legend_width", length = 1, positive = TRUE
  )
  if (!is.null(platform)) {
    if (!is.numeric(platform) || length(platform) == 0L ||
        any(!is.finite(platform)) || any(platform != as.integer(platform)) ||
        any(platform < 1L) || any(platform > x$model$n_platforms)) {
      .imr_abort(sprintf(
        "`platform` must contain whole-number indices between 1 and %d.",
        x$model$n_platforms
      ))
    }
    platform <- as.integer(platform)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    if (type == "selection") try(graphics::layout(1), silent = TRUE)
    try(graphics::par(op), silent = TRUE)
  }, add = TRUE)

  if (type == "trace") {
    lp <- x$posterior$log_posterior
    pp <- .imr_plot_par(mar, mgp, default_mar = c(4.8, 4.8, 3, 1))
    graphics::par(mar = pp$mar, mgp = pp$mgp)
    trace_col <- if (is.null(col)) .imr_plot_trace_colour() else col
    do.call(graphics::plot, c(list(
      x = seq_along(lp), y = lp, type = "l",
      xlab = "MCMC iteration", ylab = "Log-posterior",
      main = "Log-posterior trace", cex.axis = cex_axis,
      cex.lab = cex_lab, cex.main = cex_main, col = trace_col
    ), dots))
    graphics::abline(v = x$control$mcmc$burnin, lty = 2, col = "grey50")
    return(invisible(NULL))
  }

  if (type %in% c("theta_trace", "selection_trace")) {
    if (is.null(platform)) platform <- 1L
    platform <- .imr_check_integer_scalar(
      platform, "platform", min = 1L, max = x$model$n_platforms
    )
    pp <- .imr_plot_par(mar, mgp, default_mar = c(4.8, 4.8, 3, 1))
    graphics::par(mar = pp$mar, mgp = pp$mgp)
    trace_col <- if (is.null(col)) .imr_plot_trace_colour() else col
    if (type == "theta_trace") {
      samples <- x$posterior$interaction_draws[[platform]]
      if (is.null(samples) || ncol(samples) == 0L) {
        .imr_abort("The selected platform has no sampled theta interactions.")
      }
      parameter <- .imr_check_integer_scalar(
        parameter, "parameter", min = 1L, max = ncol(samples)
      )
      values <- samples[, parameter]
      title <- sprintf("Theta trace: %s, pair %d",
                       x$model$platform_names[platform], parameter)
      ylab <- "Theta"
    } else {
      template <- .imr_mpip(x, platform)
      subgroup <- .imr_check_integer_scalar(
        subgroup, "subgroup", min = 1L, max = nrow(template)
      )
      feature <- .imr_check_integer_scalar(
        feature, "feature", min = 1L, max = ncol(template)
      )
      values <- vapply(
        x$posterior$selection_draws,
        function(draw) as.numeric(draw[[platform]][subgroup, feature]),
        numeric(1L)
      )
      title <- sprintf(
        "Selection trace: %s / %s / %s", x$model$platform_names[platform],
        rownames(template)[subgroup], colnames(template)[feature]
      )
      ylab <- "Selection indicator"
    }
    do.call(graphics::plot, c(list(
      x = seq_along(values), y = values, type = "l",
      xlab = "Retained MCMC draw", ylab = ylab, main = title,
      cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_main,
      col = trace_col
    ), dots))
    return(invisible(NULL))
  }

  plats <- if (is.null(platform)) seq_len(x$model$n_platforms) else platform
  heat_palette <- if (is.null(palette)) {
    if (type == "selection") "grey" else "heatmap"
  } else {
    palette
  }
  heat_col <- if (is.null(col)) {
    .imr_plot_palette(64, palette = heat_palette)
  } else {
    .imr_plot_palette(length(col), col = col)
  }

  if (type == "selection") {
    layout_widths <- rep(1, length(plats))
    if (legend) layout_widths <- c(layout_widths, legend_width)
    graphics::layout(
      matrix(seq_along(layout_widths), nrow = 1L),
      widths = layout_widths
    )
  } else if (length(plats) > 1) {
    graphics::par(mfrow = c(1, length(plats)))
  }

  panel_par <- .imr_plot_par(
    mar, mgp,
    default_mar = if (type == "selection") {
      c(5.7, 4.8, 2.9, 0.5)
    } else {
      c(4.8, 4.8, 3, 1)
    },
    default_mgp = if (type == "selection") c(3.6, 0.9, 0) else c(2.7, 0.8, 0)
  )
  for (l in plats) {
    if (type == "selection") {
      m <- .imr_mpip(x, l)
      main <- sprintf("mPIP: %s", x$model$platform_names[l])
      xlab <- "Features"; ylab <- "Availability subgroups"
      rlab <- rownames(m)
    } else {
      m <- x$posterior$interaction_means[[l]]
      rlab <- x$model$subgroup_names[x$model$platform_subgroups[[l]]]
      if (!is.null(rlab) && length(rlab) == nrow(m)) {
        rownames(m) <- colnames(m) <- rlab
      }
      main <- sprintf("Theta: %s", x$model$platform_names[l])
      xlab <- "Availability subgroups"; ylab <- "Availability subgroups"
    }
    graphics::par(mar = panel_par$mar, mgp = panel_par$mgp)
    if (nrow(m) == 0 || ncol(m) == 0) {
      graphics::plot.new()
      graphics::title(main = paste(main, "(empty)"), cex.main = cex_main)
      next
    }
    zlim <- if (type == "selection") c(0, 1) else range(m, na.rm = TRUE)
    do.call(graphics::image, c(list(
      x = seq_len(ncol(m)), y = seq_len(nrow(m)), z = t(m),
      col = heat_col, zlim = zlim, axes = FALSE,
      xlab = xlab, ylab = ylab, main = main,
      cex.lab = cex_lab, cex.main = cex_main
    ), dots))
    if (!is.null(colnames(m)) && ncol(m) <= 40) {
      graphics::axis(1, at = seq_len(ncol(m)), labels = colnames(m),
                     las = 2, cex.axis = cex_axis)
    } else {
      graphics::axis(1, cex.axis = cex_axis)
    }
    graphics::axis(2, at = seq_len(nrow(m)), labels = rlab, las = 2,
                   cex.axis = cex_axis)
    graphics::box()
  }
  if (type == "selection" && legend) {
    legend_par <- .imr_plot_par(
      NULL, mgp, default_mar = c(5.7, 0.1, 2.9, 3.0),
      default_mgp = c(3.6, 0.9, 0)
    )
    graphics::par(mar = legend_par$mar, mgp = legend_par$mgp)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    yb <- seq(0, 1, length.out = length(heat_col) + 1L)
    graphics::rect(0.18, yb[-length(yb)], 0.52, yb[-1L],
                   col = heat_col, border = NA)
    graphics::axis(4, at = c(0, 0.25, 0.5, 0.75, 1),
                   labels = c("0", "0.25", "0.5", "0.75", "1"),
                   las = 1, cex.axis = cex_axis, tck = -0.18)
    graphics::mtext("mPIP", side = 3, line = 0.1, at = 0.35,
                    cex = 0.75 * cex_main)
    graphics::box(bty = "n")
  }
  invisible(NULL)
}

## Stand-alone plotting helpers for objects of class "imr".

#' Plot the Top Selected Features of an IMR Fit
#'
#' @description
#' Draws a horizontal bar chart of the features with the highest marginal
#' posterior inclusion probability (mPIP).  For each platform-feature pair, the
#' plotted score is the highest mPIP attained among the availability subgroups
#' containing that platform; the displayed bars are the largest such scores
#' across all platforms.  This complements the per-platform heatmap of
#' `plot(fit, type = "selection")`.
#'
#' @section Plotted quantity:
#' Each bar is the maximum subgroup mPIP for one platform-feature pair,
#' \eqn{r_{lj}=\max_{s\in\mathcal S_l}\widehat\pi_{lsj}}{r_lj = maximum subgroup mPIP for platform l and feature j} as defined in
#' [summary.imr()]. Bars are ranked jointly across platforms. This maximum
#' is the inclusion probability of a maximizing subgroup; it is not a
#' coefficient effect size or the posterior probability of selection in at
#' least one subgroup. `reference` adds a
#' visual guide and does not filter the selected `top` bars.
#'
#' @param object A fitted object of class `"imr"` returned by [imr()].
#' @param top Integer; the number of highest-mPIP features to display
#'   (default `10`).
#' @param base_cex Overall text-size multiplier. The `cex_*` arguments default
#'   to values derived from this multiplier (default `1`).
#' @param cex_names Feature-label size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_axis Axis-label size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_lab Axis-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_main Main-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_legend Legend text size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param col Optional bar colours. If `NULL`, colours are generated from
#'   `palette`.
#' @param palette Palette name used when `col = NULL`; `"platform"` gives the
#'   package's standard muted platform colours (default).
#' @param reference Optional vertical reference line; use `NULL` to suppress it
#'   (default `0.5`).
#' @param show_source Logical; should each bar be annotated with the
#'   availability subgroup (bitstring) in which the feature attained its maximum
#'   mPIP (default `TRUE`)?
#' @param legend Logical; should the platform legend be drawn (default `TRUE`)?
#' @param xlim Numeric vector of length two giving the horizontal axis limits
#'   (default `c(0, 1)`).
#' @param mar,mgp Optional graphical margin and axis-title placement vectors
#'   passed to [graphics::par()] for finer layout control.
#' @param ... Further graphical parameters passed to [graphics::barplot()].
#'
#' @return Invisibly, a data frame of the displayed features with columns
#'   `platform`, `feature`, `mpip` and `subgroup`, where `mpip` is the maximum
#'   subgroup mPIP for that platform-feature pair and `subgroup` is the
#'   subgroup where the maximum is attained. Rows are ordered by decreasing
#'   mPIP.
#' @seealso [imr()], [plot.imr()], [plot_subgroup_sizes()]
#' @examples
#' \donttest{
#' data("simIMR", package = "IntegMultiReg")
#' fit <- imr(
#'   x = simIMR$platforms, outcome = simIMR$outcome,
#'   covariates = simIMR$covariates, outcome_type = "binary",
#'   nu = c(-4, -3, -4), draws = 200, burnin = 100,
#'   min_subgroup_size = 5, seed = 1
#' )
#' plot_top_features(fit, top = 8)
#' }
#' @export
plot_top_features <- function(object, top = 10, base_cex = 1,
                              cex_names = NULL, cex_axis = NULL,
                              cex_lab = NULL, cex_main = NULL,
                              cex_legend = NULL, col = NULL,
                              palette = "platform", reference = 0.5,
                              show_source = TRUE, legend = TRUE, xlim = c(0, 1),
                              mar = NULL, mgp = NULL, ...) {
  .imr_check_fit(object)
  dots <- list(...)
  top <- .imr_check_integer_scalar(top, "top", min = 1)
  sz <- .imr_plot_cex(
    dots, base_cex = base_cex, cex_names = cex_names,
    cex_axis = cex_axis, cex_lab = cex_lab, cex_main = cex_main,
    cex_legend = cex_legend,
    defaults = list(axis = 0.95, lab = 1.05, main = 1.1,
                    names = 0.95, legend = 0.9, values = 0.9)
  )
  dots <- sz$dots
  .imr_check_flag(show_source, "show_source")
  .imr_check_flag(legend, "legend")
  xlim <- .imr_check_numeric_vector(xlim, "xlim", length = 2)
  if (xlim[1] >= xlim[2]) {
    .imr_abort("`xlim` must be increasing.")
  }
  if (!is.null(reference)) {
    reference <- .imr_check_numeric_vector(
      reference, "reference", length = 1, nonnegative = TRUE
    )
    if (reference > 1) {
      .imr_abort("`reference` must be between 0 and 1, or NULL.")
    }
  }
  ## Collect every feature's best mPIP (and the subgroup achieving it).
  rows <- list()
  validate_imr(object)
  for (l in seq_len(object$model$n_platforms)) {
    m <- .imr_mpip(object, l)
    if (nrow(m) == 0 || ncol(m) == 0) next
    maxp <- apply(m, 2, max)
    sg <- rownames(m)[apply(m, 2, which.max)]
    rows[[l]] <- data.frame(
      platform = object$model$platform_names[l],
      feature = colnames(m),
      mpip = as.numeric(maxp),
      subgroup = sg,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab) || nrow(tab) == 0) {
    .imr_abort("The fit contains no selectable features to plot.")
  }
  tab <- tab[order(tab$mpip, decreasing = TRUE), , drop = FALSE]
  top <- min(top, nrow(tab))
  tab <- tab[seq_len(top), , drop = FALSE]

  ## One colour per platform, drawn in increasing order so the largest bar is
  ## at the top of the horizontal chart.
  platforms <- object$model$platform_names
  pal <- .imr_plot_palette(length(platforms), col = col, palette = palette)
  ord <- rev(seq_len(top))
  labels <- paste0(tab$feature[ord], " (", tab$platform[ord], ")")

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  pp <- .imr_plot_par(mar, mgp, default_mar = c(4.8, 8.8, 3.2, 1))
  graphics::par(mar = pp$mar, mgp = pp$mgp)
  if (is.null(dots$border)) dots$border <- NA
  bp <- do.call(graphics::barplot, c(list(
    height = tab$mpip[ord], names.arg = labels, horiz = TRUE, las = 1,
    xlim = xlim, xlab = "Maximum subgroup mPIP",
    main = sprintf("Top %d features by max subgroup mPIP", top),
    col = pal[match(tab$platform[ord], platforms)],
    cex.names = sz$names, cex.axis = sz$axis,
    cex.lab = sz$lab, cex.main = sz$main
  ), dots))
  if (!is.null(reference)) {
    graphics::abline(v = reference, lty = 2,
                     col = .imr_plot_reference_colour())
  }
  ## Annotate each bar with the availability subgroup where its maximum mPIP was
  ## attained.  Long bars carry the label inside the coloured region (light
  ## text); short bars carry it just past the tip (dark text) so it never
  ## overplots the reference line or spills past the axis.
  if (show_source) {
    xend <- tab$mpip[ord]
    src <- tab$subgroup[ord]
    inside <- xend >= xlim[1] + 0.5 * (xlim[2] - xlim[1])
    if (any(inside)) {
      graphics::text(xend[inside], bp[inside], labels = src[inside],
                     pos = 2, offset = 0.35, col = "white",
                     font = 2, cex = sz$values)
    }
    if (any(!inside)) {
      graphics::text(xend[!inside], bp[!inside], labels = src[!inside],
                     pos = 4, offset = 0.35,
                     col = .imr_plot_reference_colour(),
                     font = 2, cex = sz$values)
    }
  }
  if (legend) {
    graphics::legend("bottomright", legend = platforms, fill = pal,
                     bty = "n", cex = sz$legend)
  }
  invisible(tab[, c("platform", "feature", "mpip", "subgroup")])
}


#' Plot the Availability Subgroup Sizes of an IMR Fit
#'
#' @description
#' Draws a bar chart of the number of subjects in each modelled
#' availability subgroup (the non-empty regions of the Venn diagram), giving a
#' quick picture of how the sample is distributed across subgroups.
#'
#' @section Counted subjects:
#' For a retained availability pattern \eqn{s}, the bar height is
#' \eqn{n_s=\sum_i I(A_i=s)}{n_s = count of subjects with availability pattern s}, where \eqn{A_i} is the subject's platform
#' availability pattern defined in [imr_data()]. Counts are the observed sample
#' sizes stored in the fit, after eligibility checks and subgroup filtering.
#' They are not posterior estimates or imputed numbers of complete cases.
#'
#' @param object A fitted object of class `"imr"` returned by [imr()].
#' @param base_cex Overall text-size multiplier. The `cex_*` arguments default
#'   to values derived from this multiplier (default `1`).
#' @param cex_axis Axis-label size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_lab Axis-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_main Main-title size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param cex_values Bar-value label size multiplier. If `NULL`, a plot-specific
#'   default derived from `base_cex` is used.
#' @param col Optional bar colours. If `NULL`, colours are generated from
#'   `palette`.
#' @param palette Palette name used when `col = NULL`; `"subgroup"` gives the
#'   package's standard neutral grey-blue subgroup colours (default).
#' @param show_values Logical; should sample sizes be printed above the bars
#'   (default `TRUE`)?
#' @param ylim Optional numeric vector of length two giving the vertical axis
#'   limits. If `NULL`, limits are chosen from the subgroup sizes.
#' @param mar,mgp Optional graphical margin and axis-title placement vectors
#'   passed to [graphics::par()] for finer layout control.
#' @param ... Further graphical parameters passed to [graphics::barplot()].
#'
#' @return Invisibly, the named integer vector of subgroup sizes.
#' @seealso [imr()], [plot.imr()], [plot_top_features()]
#' @examples
#' \donttest{
#' data("simIMR", package = "IntegMultiReg")
#' fit <- imr(
#'   x = simIMR$platforms, outcome = simIMR$outcome,
#'   covariates = simIMR$covariates, outcome_type = "binary",
#'   nu = c(-4, -3, -4), draws = 200, burnin = 100,
#'   min_subgroup_size = 5, seed = 1
#' )
#' plot_subgroup_sizes(fit)
#' }
#' @export
plot_subgroup_sizes <- function(object, base_cex = 1, cex_axis = NULL,
                                cex_lab = NULL, cex_main = NULL,
                                cex_values = NULL, col = NULL,
                                palette = "subgroup", show_values = TRUE,
                                ylim = NULL, mar = NULL, mgp = NULL, ...) {
  .imr_check_fit(object)
  dots <- list(...)
  sz <- .imr_plot_cex(
    dots, base_cex = base_cex, cex_axis = cex_axis,
    cex_lab = cex_lab, cex_main = cex_main, cex_values = cex_values,
    defaults = list(axis = 0.95, lab = 1.05, main = 1.1,
                    names = 1, legend = 0.9, values = 0.95)
  )
  dots <- sz$dots
  .imr_check_flag(show_values, "show_values")
  validate_imr(object)
  sizes <- as.integer(object$model$sample_sizes)
  names(sizes) <- object$model$subgroup_names
  if (is.null(ylim)) {
    ylim <- c(0, max(sizes) * 1.15)
  } else {
    ylim <- .imr_check_numeric_vector(ylim, "ylim", length = 2)
    if (ylim[1] >= ylim[2]) {
      .imr_abort("`ylim` must be increasing.")
    }
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  pp <- .imr_plot_par(mar, mgp, default_mar = c(4.8, 4.8, 3.2, 1))
  graphics::par(mar = pp$mar, mgp = pp$mgp)
  if (is.null(dots$border)) dots$border <- NA
  bp <- do.call(graphics::barplot, c(list(
    height = sizes, xlab = "Availability subgroup (bitstring)",
    ylab = "Number of subjects", main = "Subjects per availability subgroup",
    col = .imr_plot_palette(length(sizes), col = col, palette = palette),
    ylim = ylim, cex.axis = sz$axis,
    cex.lab = sz$lab, cex.main = sz$main
  ), dots))
  if (show_values) {
    graphics::text(bp, sizes, labels = sizes, pos = 3, cex = sz$values)
  }
  invisible(sizes)
}
