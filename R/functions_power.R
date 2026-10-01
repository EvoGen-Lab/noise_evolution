###############################################################
### functions_power.R
### Functions for the power analysis (Section 10), sourced by
### power_analysis.R and by the SLURM cluster script power_grid.R
### after functions.R, whose .fit_one(), .fit_split() and
### perm_pval() they build on.
###
### Outline (function name - purpose):
###     open_grid_pdf() - Opens a PDF sized to the panel grid it is about to hold, with tight margins throughout.
###     plot_lines() - Plots one line per column of a matrix against a shared x vector, with reference lines at power 0.05 and 0.9.
###     legend_page() - One legend page mapping each line color to the value it represents.
###     power_line_figure() - Multi-panel PDF of power curves from the POWER array: one panel per (outer, inner) pair, one line per series value, with a legend page.
###     power_grid_row() - power of one grid row i (mean reads x cell count x burst-frequency size) over the SIZE.RATIO sweep.
###############################################################

## ============================================================
## Panel layout and power-curve figures (power_analysis.R)
## ============================================================
## Every power figure shares one compact panel layout, so the grids drawn
## over different pairs of axes read the same way.

## Opens a PDF sized to the nr x nc panel grid it is about to hold.
## Margins and label spacing (mgp) are set tight throughout, so every
## power-analysis figure shares the same compact panel layout.
open_grid_pdf <- function(file, nr, nc, panel_w = 2.3, panel_h = 2.3, mar = c(3, 3, 2, 1), mgp = c(1.6, 0.5, 0)) {
  pdf(file, width = nc * panel_w, height = nr * panel_h, useDingbats = FALSE)
  par(mfrow = c(nr, nc), mar = mar, mgp = mgp)
}

## Plots one line per column of `mat` against `xv`, using a color ramp
## sized to ncol(mat). Reference lines at power 0.05 (nominal false
## positive rate) and 0.9 (a common target) give every panel the same
## visual anchor regardless of which quantity is on the x-axis.
plot_lines <- function(xv, mat, xlab = "", ylab = "Power", main = "", show_axes = TRUE, lwd = 1.4, x_at = NULL, x_labels = NULL) {
  cols <- colorRampPalette(POWER.COLOR)(ncol(mat))
  plot(xv, mat[, 1], type = "l", lwd = lwd, col = cols[1], ylim = c(0, 1), xlab = if (show_axes && is.null(x_at)) xlab else "", ylab = if (show_axes) ylab else "", xaxt = if (!is.null(x_at) || !show_axes) "n" else "s", yaxt = if (show_axes) "s" else "n", main = main, cex.main = 0.8)
  abline(h = 0.05, col = COLOR.ACCENT, lty = 2); abline(h = 0.9, col = COLOR.ACCENT, lty = 3)
  if (ncol(mat) > 1) for (k in 2:ncol(mat)) lines(xv, mat[,k], lwd = lwd, col = cols[k])
  #x_at/x_labels give exact tick positions and real-valued labels for a
  #log-scaled axis (e.g. log2(SIZE.RATIO) with SIZE.RATIO as labels)
  if (!is.null(x_at)) {
    axis(1, at = x_at, labels = x_labels, cex.axis = 0.55, las = 2)
    title(xlab = xlab, mgp = c(2.6,1,0), cex.lab = 0.75)
  }
  invisible(cols)
}

## One legend page appended after a grid of panels, mapping each line
## color to the value it represents.
legend_page <- function(labels, title) {
  cols <- colorRampPalette(POWER.COLOR)(length(labels))
  plot.new()
  legend("center", legend = labels, col = cols, lwd = 2, title = title, ncol = ceiling(length(labels)/15), bty = "n", cex = 0.9)
}

## power_line_figure: one multi-panel PDF of power curves from the POWER array (dimensions M mean reads, N cell
## count, P burst-frequency SIZE, Q SIZE.RATIO; grid holds the value of each dimension). Panels run over the
## outer and inner dimension; every panel draws one line per value of the series dimension against the x
## dimension, with the remaining two dimensions fixed at the panel's values. main_fmt is a sprintf format whose
## %1$s is the outer and %2$s the inner value. The figure ends with a legend page (legend_labels, legend_title).
## Plain style (styled = NULL) lets plot_lines() draw the axes (xlab; ticks = TRUE puts exact tick labels on the
## log-scaled axis). Styled panels draw percent labels on the y axis and a custom x axis described by
## styled = list(at, labels, cex_axis, mgp, xlab, xlab_mgp, xlab_cex[, las, tcl]). The PDF has one row per outer
## value and one column per inner value.
power_line_figure <- function(file, power, grid, outer, inner, x, series, main_fmt, legend_labels, legend_title,
                              log_x = TRUE, xlab = "", ticks = FALSE, styled = NULL, mar = c(3, 3, 2, 1)) {
  dim_ix <- c(M = 1, N = 2, P = 3, Q = 4)
  open_grid_pdf(file.path(FIGURE.DIR, "extended", file), nr = length(grid[[outer]]), nc = length(grid[[inner]]), mar = mar)
  xv <- if (log_x) log2(grid[[x]]) else grid[[x]]
  for (o in seq_along(grid[[outer]])) {
    for (i in seq_along(grid[[inner]])) {
      idx <- list(M = TRUE, N = TRUE, P = TRUE, Q = TRUE)
      idx[[outer]] <- o; idx[[inner]] <- i
      sl  <- do.call(`[`, c(list(power), unname(idx)))
      mat <- if (dim_ix[[x]] > dim_ix[[series]]) t(sl) else sl
      main <- sprintf(main_fmt, grid[[outer]][o], grid[[inner]][i])
      if (is.null(styled)) {
        if (ticks) plot_lines(xv, mat, xlab = xlab, ylab = "Power", main = main, x_at = xv, x_labels = grid[[x]])
        else plot_lines(xv, mat, xlab = xlab, ylab = "Power", main = main)
      } else {
        plot_lines(xv, mat, show_axes = FALSE)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1, labels = c("0%", "20%", "40%", "60%", "80%", "100%"), cex.axis = 0.6, mgp = c(3, 0.6, 0))
        title(ylab = "Power", mgp = c(2, 1, 0))
        do.call(axis, c(list(1, at = styled$at, labels = styled$labels, cex.axis = styled$cex_axis, mgp = styled$mgp), styled[intersect(c("las", "tcl"), names(styled))]))
        title(xlab = styled$xlab, mgp = styled$xlab_mgp, cex.lab = styled$xlab_cex)
        title(main = main, cex.main = 0.8, line = 0.2)
      }
    }
  }
  legend_page(legend_labels, legend_title)
  dev.off()
}

## ============================================================
## Simulation of one grid row (the cluster script power_grid.R)
## ============================================================

## power_grid_row: power of one grid row i (mean reads x cell count x burst-frequency size) over the SIZE.RATIO sweep. The reference group's
## counts and fit are drawn once per replicate and reused across the sweep, so the sweep isolates the effect of SIZE.RATIO.
power_grid_row <- function(i, ALPHA, CELL.RATIO, EXPOSURE.CV, GRID, MEAN.READS, N.CELLS, N.MIX, NI, NJ, PI1, SEED.BASE, SIZE, SIZE.RATIO) {
  ## size_log2_ratio: log2 ratio of the NB size (disp, the burst-frequency axis) between two
  ## .fit_one() results, NA if either side is non-positive or non-finite.
  size_log2_ratio <- function(fit_a, fit_b) {
    a <- fit_a[["disp"]]; b <- fit_b[["disp"]]
    if (is.finite(a) && a > 0 && is.finite(b) && b > 0) log2(a) - log2(b) else NA_real_
  }

  row <- GRID[i, ]
  MEAN.READS.X <- MEAN.READS[row$m]
  N.SC.X       <- N.CELLS[row$n]
  N.SE.X       <- round(N.SC.X * CELL.RATIO)
  SIZE.1.X     <- SIZE[row$p]

  set.seed(SEED.BASE + 1e6 + row$n)  #shared draws, keyed only by cell count
  sdlog <- sqrt(log(1 + EXPOSURE.CV^2))
  EXPO.X.LIST <- replicate(NJ, rlnorm(N.SC.X, meanlog = -0.5*sdlog^2, sdlog = sdlog), simplify = FALSE)
  EXPO.Y.LIST <- replicate(NJ, rlnorm(N.SE.X, meanlog = -0.5*sdlog^2, sdlog = sdlog), simplify = FALSE)
  PERM.LIST   <- replicate(NI, sample.int(N.SC.X + N.SE.X), simplify = FALSE)

  set.seed(SEED.BASE + i)  #row-specific draws (X depends on m, p, n)
  X.LIST  <- vector("list", NJ)
  FA.LIST <- vector("list", NJ)
  for (j in seq_len(NJ)) {
    X.LIST[[j]]  <- rnbinom(n = N.SC.X, size = SIZE.1.X, mu = MEAN.READS.X*EXPO.X.LIST[[j]])
    FA.LIST[[j]] <- .fit_one(X.LIST[[j]], EXPO.X.LIST[[j]])  #shared across SIZE.RATIO below
  }

  P.ALL <- matrix(NA_real_, NJ, length(SIZE.RATIO))
  for (qi in seq_along(SIZE.RATIO)) {
    SIZE.2.X <- SIZE.1.X / SIZE.RATIO[qi]
    for (j in seq_len(NJ)) {
      EXPO.Y <- EXPO.Y.LIST[[j]]
      Y  <- rnbinom(n = N.SE.X, size = SIZE.2.X, mu = MEAN.READS.X*EXPO.Y)
      fb <- .fit_one(Y, EXPO.Y)   #same estimator as the null below
      OBS <- size_log2_ratio(FA.LIST[[j]], fb)

      XY <- c(X.LIST[[j]], Y); EXPO.XY <- c(EXPO.X.LIST[[j]], EXPO.Y)
      NULL.DIST <- numeric(NI)
      for (k in seq_len(NI)) {
        sp <- .fit_split(XY, EXPO.XY, PERM.LIST[[k]], N.SC.X)  #same estimator as the observed contrast
        NULL.DIST[k] <- size_log2_ratio(sp$a, sp$b)
      }
      P.ALL[j, qi] <- perm_pval(OBS, NULL.DIST)
    }
  }

  i0     <- which(SIZE.RATIO == 1)
  P.NULL <- P.ALL[, i0]
  n_alt  <- min(NJ, round(NJ * PI1 / (1 - PI1)))
  vapply(seq_along(SIZE.RATIO), function(qi) {
    if (qi == i0) return(mean(P.NULL < ALPHA, na.rm = TRUE))
    mean(replicate(N.MIX, {
      q <- p.adjust(c(P.NULL, P.ALL[sample.int(NJ, n_alt), qi]), method = "BH")
      mean(q[NJ + seq_len(n_alt)] < ALPHA, na.rm = TRUE)
    }), na.rm = TRUE)
  }, numeric(1))
}
