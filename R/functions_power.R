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
###     gamma_unit() - Gamma draws with mean 1 and a given CV^2 (1 when the CV^2 is 0).
###     sim_allele_counts() - Counts of one hybrid allele: Poisson with rate mu * expo * u * v, u the factor shared by both alleles of a cell and v the allele's own.
###     paired_alleles_check() - Checks the paired-allele generator: recovered NB size against the target and the allele correlation against phi.
###     bh_power() - Power at a Benjamini-Hochberg FDR from the p-values of null and alternative datasets (the rule power_grid_row() uses).
###     power_modes_row() - power of the cis and trans contrasts (trans under two nulls) for one grid row, on paired-allele hybrid data, with the pipeline's own nulls.
###     mde_interp() - Smallest log2 effect at which a power curve reaches a target, by linear interpolation.
###     modes_long() - Stacks the per-row results of power_modes_row() into long tables with the grid values attached.
###     modes_mde() - Minimum detectable log2 ratio at a target power for every grid point, mode and axis.
###     modes_f_cross() - The hybrid split fraction at which cis and trans reach equal power, interpolated over f.
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

## ============================================================
## Power of the cis and trans contrasts on paired-allele hybrid data (power_modes_grid.R)
## ============================================================
## power_grid_row() simulates one two-group contrast. The functions below simulate the contrasts the
## pipeline actually tests. cis compares the two alleles of the same hybrid cells and is nulled by the
## within-cell allele swap. trans is the parental ratio minus the hybrid allele ratio and is nulled by
## pooling parent and hybrid cells (two versions: independent allele shuffles, as in permute_contrasts_one(),
## and a pairing-preserving one). The nulls come from the helpers the pipeline calls (null_cis_axes(),
## null_trans_axes(), perm_trans_pool(), perm_trans_pool_paired(), fit_ratio_axes()).

## gamma_unit: n Gamma draws with mean 1 and CV^2 cv2 (shape = rate = 1 / cv2); all ones when cv2 is 0.
gamma_unit <- function(n, cv2) if (cv2 > 0) rgamma(n, shape = 1 / cv2, rate = 1 / cv2) else rep(1, n)

## sim_allele_counts: counts of one hybrid allele in cells with exposure expo. Each count is Poisson with rate
## mu * expo * u * v. u is the extrinsic factor, drawn once per cell by the caller and shared by both alleles
## of the cell; v is the allele's own (intrinsic) gamma factor. The extrinsic CV^2 is cv2_ext. The allele's
## total latent CV^2 is set to 1 / size_a, so (1 + cv2_ext)(1 + cv2_int) - 1 = 1 / size_a fixes the intrinsic
## CV^2. The count variance is mean + mean^2 / size_a exactly, and the marginal law is close to negative
## binomial with that size (a product of gammas is not exactly gamma; paired_alleles_check() tests the fit).
## The latent correlation between two alleles of equal size is cv2_ext * size_a, which is phi when
## cv2_ext = phi / size_a.
sim_allele_counts <- function(mu, expo, u, size_a, cv2_ext) {
  stopifnot(is.numeric(mu), mu > 0, size_a > 0, cv2_ext >= 0, length(u) == length(expo))
  cv2_int <- (1 + 1 / size_a) / (1 + cv2_ext) - 1
  stopifnot(cv2_int >= -1e-9)
  rpois(length(expo), mu * expo * u * gamma_unit(length(expo), max(cv2_int, 0)))
}

## paired_alleles_check: validates the generator for one (mu, size, phi) on n cells. The extrinsic CV^2 is
## phi / size and the intrinsic CV^2 solves the total latent CV^2 = 1 / size. Fits each allele with .fit_one()
## and stops unless the recovered NB size is within tol_size (relative) of `size` and the latent allele
## correlation is within tol_cor of phi. Returns one row of the recovered size, the correlation and their targets.
paired_alleles_check <- function(mu, size, phi, n = 40000, expo_cv = 0.8, tol_size = 0.1, tol_cor = 0.03) {
  stopifnot(phi >= 0, phi <= 1)
  sdlog <- sqrt(log(1 + expo_cv^2))
  expo  <- rlnorm(n, -0.5 * sdlog^2, sdlog)
  cv2_ext <- phi / size
  u <- gamma_unit(n, cv2_ext)
  cv2_int <- (1 + 1 / size) / (1 + cv2_ext) - 1
  v1 <- gamma_unit(n, cv2_int); v2 <- gamma_unit(n, cv2_int)
  y1 <- rpois(n, mu * expo * u * v1); y2 <- rpois(n, mu * expo * u * v2)
  rec <- c(.fit_one(y1, expo)[["disp"]], .fit_one(y2, expo)[["disp"]])
  cor_lat <- cor(u * v1, u * v2)
  stopifnot("recovered NB size is outside the tolerance" = all(abs(rec / size - 1) <= tol_size),
            "allele correlation is outside the tolerance of phi" = abs(cor_lat - phi) <= tol_cor)
  data.frame(mu = mu, size = size, phi = phi, size_rec = mean(rec), cor_latent = cor_lat,
             cor_count = cor(y1 / expo, y2 / expo))
}

## bh_power: power at a Benjamini-Hochberg FDR. p_null: p-values of datasets simulated under the null
## (ratio 1); p_alt: p-values of the same number of datasets under an alternative. A random subset of
## the alternatives (a fraction pi1 of the pooled set) is pooled with the nulls, BH is applied, and power is the
## fraction of the alternatives called at q < alpha; the mean over n_mix random subsets uses every dataset.
## This is the rule power_grid_row() applies to the total contrast.
bh_power <- function(p_null, p_alt, alpha, pi1, n_mix) {
  nj <- length(p_null); stopifnot(length(p_alt) == nj, pi1 > 0, pi1 < 1)
  n_alt <- min(nj, round(nj * pi1 / (1 - pi1)))
  mean(replicate(n_mix, {
    q <- p.adjust(c(p_null, p_alt[sample.int(nj, n_alt)]), method = "BH")
    mean(q[nj + seq_len(n_alt)] < alpha, na.rm = TRUE)
  }), na.rm = TRUE)
}

## mde_interp: the smallest log2 effect at which power reaches `target`, by linear interpolation of
## power against log2(ratios) between the two bracketing points (ratios ascending, the first being 1). NA
## when the curve never reaches the target, 0 when it starts above it.
mde_interp <- function(power, ratios, target = 0.8) {
  stopifnot(length(power) == length(ratios), !is.unsorted(ratios))
  x <- log2(ratios)
  ok <- is.finite(power)
  x <- x[ok]; power <- power[ok]
  if (length(power) == 0 || max(power) < target) return(NA_real_)
  k <- which(power >= target)[1]
  if (k == 1) return(x[1])
  x[k - 1] + (target - power[k - 1]) / (power[k] - power[k - 1]) * (x[k] - x[k - 1])
}

## power_modes_row: power of the cis and trans contrasts for one row of GRID (indices m, p, h, f, n, d into
## MEAN.READS, SIZE, PHI, FRAC, DESIGN and DEPTH.RATIO) over the SIZE.RATIO sweep. SIZE.RATIO is the effect on
## the ratio of NB sizes (the burst-frequency axis), applied as
##   null          ratio 1 everywhere;
##   cis effect    hybrid allele ratio = ratio, parental ratio = ratio, so trans = 0;
##   trans effect  hybrid allele ratio = 1, parental ratio = ratio, so cis = 0.
## Sc-side groups have NB size SIZE and Se-side groups SIZE / ratio. Hybrid cells come with both alleles
## (sim_allele_counts()): the shared fraction phi sets the extrinsic CV^2 to phi / SIZE and the alleles' own CV^2
## makes the total latent CV^2 1 / size. phi = 0 gives unpaired alleles. The N.H hybrid cells of a design are split
## as in the pipeline: round(f * N.H) cells (HYC) test cis and the rest (HYT) test trans. Hybrid allele mean reads
## are MEAN.READS * DEPTH.RATIO, the hybrid-to-parent per-allele depth ratio. Exposures are lognormal with CV
## EXPOSURE.CV. Parental exposures are keyed by the design (as in power_grid_row()), hybrid exposures and the
## null relabelings by design and f, so rows that share them draw the same values; the counts use their own
## seed (SEED.BASE + i).
## Each dataset's contrasts are tested on every axis (mu, disp, cv2) against the pipeline's nulls with
## perm_pval(): cis by the allele swap (null_cis_axes()), trans by pooled relabelings (null_trans_axes()) in two
## versions, "trans_indep" (independent Sc and Se shuffles, perm_trans_pool()) and "trans_paired" (the same hybrid
## cells move in both alleles, perm_trans_pool_paired()), which share the observed contrast. Power is bh_power()
## at ratios above 1. At ratio 1 the entry is the raw false-positive rate at p < ALPHA (type I error), and the
## raw rejection rate at p < ALPHA is returned for every ratio. Returns list(power = data.frame of
## mode, ratio, axis, power_bh, power_raw; sd_null = SD of the observed contrast at ratio 1, by mode and axis).
power_modes_row <- function(i, ALPHA, CELL.RATIO, DEPTH.RATIO, DESIGN, EXPOSURE.CV, FRAC, GRID, MEAN.READS, N.MIX, NI, NJ,
                            PHI, PI1, SEED.BASE, SIZE, SIZE.RATIO) {
  row <- GRID[i, ]
  mu_p <- MEAN.READS[row$m]; size <- SIZE[row$p]; phi <- PHI[row$h]; frac <- FRAC[row$f]
  n_p <- DESIGN$N.P[row$n]; n_h <- DESIGN$N.H[row$n]; mu_h <- mu_p * DEPTH.RATIO[row$d]
  n_se <- round(n_p * CELL.RATIO)
  n_cis <- max(2, round(frac * n_h)); n_trans <- n_h - n_cis
  stopifnot(n_trans >= 2, any(SIZE.RATIO == 1), phi >= 0, phi <= 1, phi <= min(SIZE.RATIO))
  modes <- c("cis", "trans_indep", "trans_paired"); axes <- c("mu", "disp", "cv2")
  sdlog <- sqrt(log(1 + EXPOSURE.CV^2))
  expo_draw <- function(n) rlnorm(n, meanlog = -0.5 * sdlog^2, sdlog = sdlog)
  cv2_ext <- phi / size

  set.seed(SEED.BASE + 1e6 + row$n)                      # parental exposures, keyed by design
  EXPO.PSC <- replicate(NJ, expo_draw(n_p), simplify = FALSE)
  EXPO.PSE <- replicate(NJ, expo_draw(n_se), simplify = FALSE)
  set.seed(SEED.BASE + 2e6 + 1000 * row$n + row$f)       # hybrid exposures, keyed by design and f
  EXPO.C <- replicate(NJ, expo_draw(n_cis), simplify = FALSE)
  EXPO.T <- replicate(NJ, expo_draw(n_trans), simplify = FALSE)
  set.seed(SEED.BASE + 3e6 + 1000 * row$n + row$f)       # null relabelings, keyed by design and f
  SWAP     <- replicate(NI, runif(n_cis) < 0.5, simplify = FALSE)
  PERM.IND <- replicate(NI, list(sc = perm_trans_pool(n_p, seq_len(n_trans), seq_len(n_trans))$m,
                                 se = perm_trans_pool(n_se, seq_len(n_trans), seq_len(n_trans))$m), simplify = FALSE)
  PERM.PAI <- replicate(NI, perm_trans_pool_paired(n_p, n_se, n_trans), simplify = FALSE)

  set.seed(SEED.BASE + i)
  nq <- length(SIZE.RATIO)
  P   <- array(NA_real_, c(length(modes), NJ, nq, 3), dimnames = list(modes, NULL, NULL, axes))
  OBS <- array(NA_real_, c(2, NJ, nq, 3), dimnames = list(c("cis", "trans"), NULL, NULL, axes))
  pvals <- function(obs, null) vapply(axes, function(a) perm_pval(obs[[a]], null[a, ]), numeric(1))
  for (j in seq_len(NJ)) {
    ## Draws shared across the ratio sweep: the Sc parent, the trans-effect hybrid (allele ratio 1), and the
    ## Sc allele and extrinsic factor of the cis-effect hybrid
    par_sc <- rnbinom(n_p, size = size, mu = mu_p * EXPO.PSC[[j]]); fit_psc <- .fit_one(par_sc, EXPO.PSC[[j]])
    u_t <- gamma_unit(n_trans, cv2_ext)
    t_sc <- sim_allele_counts(mu_h, EXPO.T[[j]], u_t, size, cv2_ext)
    t_se <- sim_allele_counts(mu_h, EXPO.T[[j]], u_t, size, cv2_ext)
    fit_tsc <- .fit_one(t_sc, EXPO.T[[j]]); fit_tse <- .fit_one(t_se, EXPO.T[[j]])
    u_c <- gamma_unit(n_cis, cv2_ext)
    c_sc <- sim_allele_counts(mu_h, EXPO.C[[j]], u_c, size, cv2_ext); fit_csc <- .fit_one(c_sc, EXPO.C[[j]])
    for (qi in seq_len(nq)) {
      ratio <- SIZE.RATIO[qi]
      par_se <- rnbinom(n_se, size = size / ratio, mu = mu_p * EXPO.PSE[[j]]); fit_pse <- .fit_one(par_se, EXPO.PSE[[j]])
      c_se <- sim_allele_counts(mu_h, EXPO.C[[j]], u_c, size / ratio, cv2_ext); fit_cse <- .fit_one(c_se, EXPO.C[[j]])
      obs_cis   <- fit_ratio_axes(fit_csc, fit_cse)
      obs_trans <- fit_ratio_axes(fit_psc, fit_pse) - fit_ratio_axes(fit_tsc, fit_tse)
      null_cis   <- vapply(SWAP, function(s) null_cis_axes(c_sc, c_se, EXPO.C[[j]], s), numeric(3))
      null_ind   <- vapply(PERM.IND, function(p) null_trans_axes(par_sc, par_se, t_sc, t_se, EXPO.PSC[[j]], EXPO.PSE[[j]], EXPO.T[[j]], p$sc, p$se), numeric(3))
      null_pai   <- vapply(PERM.PAI, function(p) null_trans_axes(par_sc, par_se, t_sc, t_se, EXPO.PSC[[j]], EXPO.PSE[[j]], EXPO.T[[j]], p$sc, p$se), numeric(3))
      P[1, j, qi, ] <- pvals(obs_cis, null_cis)
      P[2, j, qi, ] <- pvals(obs_trans, null_ind)
      P[3, j, qi, ] <- pvals(obs_trans, null_pai)
      OBS[1, j, qi, ] <- obs_cis; OBS[2, j, qi, ] <- obs_trans
    }
  }
  i0 <- which(SIZE.RATIO == 1)[1]
  power <- do.call(rbind, lapply(modes, function(md) do.call(rbind, lapply(axes, function(ax) {
    pn <- P[md, , i0, ax]
    data.frame(mode = md, ratio = SIZE.RATIO, axis = ax,
               power_bh = vapply(seq_len(nq), function(qi) if (qi == i0) mean(pn < ALPHA, na.rm = TRUE)
                                 else bh_power(pn, P[md, , qi, ax], ALPHA, PI1, N.MIX), numeric(1)),
               power_raw = vapply(seq_len(nq), function(qi) mean(P[md, , qi, ax] < ALPHA, na.rm = TRUE), numeric(1)))
  }))))
  sd_null <- do.call(rbind, lapply(c("cis", "trans"), function(md)
    data.frame(mode = md, axis = axes, sd = vapply(axes, function(ax) sd(OBS[md, , i0, ax], na.rm = TRUE), numeric(1)),
               mean = vapply(axes, function(ax) mean(OBS[md, , i0, ax], na.rm = TRUE), numeric(1)))))
  list(power = power, sd_null = sd_null)
}

## modes_long: stacks the list `res` of power_modes_row() results (one per row of GRID) into two data.frames, power
## (mode, ratio, axis, power_bh, power_raw) and sd_null (mode, axis, sd, mean of the observed null contrast), each with
## the row index i and the grid values (mean reads, size, phi, f, N.P, N.H, depth ratio) of its row.
modes_long <- function(res, GRID, MEAN.READS, SIZE, PHI, FRAC, DESIGN, DEPTH.RATIO) {
  vals <- data.frame(i = seq_len(nrow(GRID)), mean_reads = MEAN.READS[GRID$m], size = SIZE[GRID$p], phi = PHI[GRID$h],
                     f = FRAC[GRID$f], N.P = DESIGN$N.P[GRID$n], N.H = DESIGN$N.H[GRID$n], depth = DEPTH.RATIO[GRID$d])
  stack <- function(part) do.call(rbind, lapply(seq_along(res), function(i) cbind(vals[i, ], res[[i]][[part]], row.names = NULL)))
  list(power = stack("power"), sd_null = stack("sd_null"))
}

## modes_mde: the minimum detectable log2 ratio (mde_interp()) at power `target` on the BH-power curve of every
## (grid row, mode, axis) in the long power table `pw` from modes_long(). Returns the grid values of each row with a mde_log2 column.
modes_mde <- function(pw, target = 0.8) {
  key <- interaction(pw$i, pw$mode, pw$axis, drop = TRUE)
  do.call(rbind, lapply(split(pw, key), function(d) {
    d <- d[order(d$ratio), ]
    out <- d[1, setdiff(names(d), c("ratio", "power_bh", "power_raw"))]
    out$mde_log2 <- mde_interp(d$power_bh, d$ratio, target)
    out
  }))
}

## modes_f_cross: the split fraction f at which cis and trans have equal power. For each combination of the grid
## values other than f (and each axis), the power of cis minus the power of the trans version `trans_mode` at SIZE.RATIO ==
## `ratio` is interpolated linearly over f, and the first sign change gives f_cross (NA when the two curves do not
## cross). pw: long power table from modes_long().
modes_f_cross <- function(pw, ratio, trans_mode = "trans_indep", axis = "disp") {
  d <- pw[pw$ratio == ratio & pw$axis == axis & pw$mode %in% c("cis", trans_mode), ]
  wide <- reshape(d[, c("mean_reads", "size", "phi", "f", "N.P", "N.H", "depth", "mode", "power_bh")],
                  idvar = c("mean_reads", "size", "phi", "f", "N.P", "N.H", "depth"), timevar = "mode", direction = "wide")
  wide$diff <- wide[["power_bh.cis"]] - wide[[paste0("power_bh.", trans_mode)]]
  keys <- c("mean_reads", "size", "phi", "N.P", "N.H", "depth")
  do.call(rbind, lapply(split(wide, interaction(wide[keys], drop = TRUE)), function(w) {
    w <- w[order(w$f), ]
    k <- which(w$diff[-nrow(w)] * w$diff[-1] < 0)[1]
    fc <- if (is.na(k)) NA_real_ else w$f[k] + (0 - w$diff[k]) / (w$diff[k + 1] - w$diff[k]) * (w$f[k + 1] - w$f[k])
    cbind(w[1, keys], f_cross = fc, ratio = ratio, trans_mode = trans_mode, axis = axis)
  }))
}
