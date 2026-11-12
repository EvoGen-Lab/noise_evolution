###############################################################
### estimator_mse.R
### Bias, variance and mean squared error of the two dispersion estimators on the experiment's own sequencing depths: the pipeline's Cox-Reid adjusted
### NB likelihood (fit_dispersion(method = "NB")) and the Appendix A weighted method of moments (method = "weighted"), both truncated at zero as in the
### pipeline. The comparison bears on the choice of the primary estimator (Notes for revision: Cox-Reid adjustment): Cox-Reid is the more precise when the
### counts are NB, and the moment estimator the less biased when the expression factor is not gamma.
### Needs only the Section 1 checkpoint. The depths are the Section 2 offsets of the observed cells: a parent's library over the median parental library
### (Pc, Pe) and half the hybrid depth for a hybrid allele (H); every simulated dataset uses all the observed cells of its group. The mu and phi values span
### the data: mu (the expected count of a cell at the reference depth) runs over MSE.N.MU log-spaced values between the MSE.MU.Q quantiles of the genes'
### fitted rates, and phi takes the MSE.PHI.Q quantiles of the positive Cox-Reid fits of MSE.FIT.GENES genes per group (share corrected in the hybrid).
### Each cell of the grid (group x factor family x mu x phi) simulates MSE.R datasets (mse_cell()); the bias, variance and MSE are averaged over
### MSE.NBATCH batches of the replicates with their SE. The simulated hybrid alleles carry no transcriptome share variation (eta about 7e-5 in the data).
### Results: results/checks/M1.rds and M1_estimator_mse.csv (raw and relative scales), a row in summary.csv, and one PDF per measure in
### figures/extra/check_M1_<measure>.pdf: x = mu (log scale), one colored line per phi, the mean as a solid line and mean +- 1 SE as light lines, panels by
### group (rows) and estimator (columns), one page per factor family. MSE.PLOT.SCALE "relative" plots bias / phi, variance / phi^2 and MSE / phi^2 so the
### phi values share an axis; "raw" plots the raw quantities.
### M2 repeats the comparison on the log2 total, cis and trans contrasts the pipeline tests, with a shared extrinsic factor for the two alleles of a hybrid
### cell, the hybrid cells split into cis and trans cells, and three designs (null, a cis effect, a trans effect); results in M2.rds, M2_contrast_mse.csv,
### M2_contrast_mse_ratio.csv and figures/extra/check_M2_<measure>.pdf. M3 measures the false-positive rate of the pipeline's dispersion test (shifted
### bootstrap null, resolution gate) when a rare cell state occurs at the same frequency in every group; results in M3.rds, M3_rare_state_type1.csv and
### figures/extra/check_M3_type1.pdf.
### Runtime: M1 about MSE.R x 2 fits per cell; at MSE.R = 500 on one core expect roughly 30 to 60 minutes for the default grid. M2 fits six groups with
### both estimators per dataset; at MSE.C.R = 300 and 5 mu values it takes about as long as M1.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT   <- "analysis/checks/estimator_mse.R"
CORES          <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.MSE       <- 9026
MSE.R          <- 500                     # simulated datasets per grid cell
MSE.NBATCH     <- 20                      # batches of the replicates for the SEs
MSE.N.MU       <- 10                      # mu values, log-spaced
MSE.MU.Q       <- c(0.05, 0.95)           # quantiles of the genes' fitted rates that bound the mu range
MSE.PHI.Q      <- c(0.10, 0.25, 0.50, 0.75, 0.90)   # quantiles of the positive fitted dispersions that give the phi values
MSE.FIT.GENES  <- 1000                    # genes per group fitted to find the phi quantiles
MSE.MIN.MEAN   <- 0.02                    # a gene enters the ranges with this mean count per cell and
MSE.MIN.CELLS  <- 50                      # this many cells with a count in its group
MSE.FAMILIES   <- c("gamma", "lognormal", "bimodal")
MSE.GROUPS     <- c("Pc", "Pe", "H")
MSE.PLOT.SCALE <- "relative"              # "relative" or "raw"
MSE.RUN        <- c("M1", "M2", "M3")     # which parts to run; M2 and M3 need only the ranges, not M1
set.seed(SEED.SECTION + 25)

stopifnot("estimator_mse.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE (cells after QC), C.BAR, R.C, R.E

## Offsets of Section 2: D is the median parental library; a hybrid allele is offset by half the cell's hybrid depth.
DEPTH.REF.M <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
MSE.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF.M, Pe = colSums(MIX.SE) / DEPTH.REF.M, H = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF.M / 2)
cat(sprintf("Cells: Pc %d, Pe %d, hybrid %d; offset CV %.2f, %.2f, %.2f\n", length(MSE.EXPOS$Pc), length(MSE.EXPOS$Pe), length(MSE.EXPOS$H),
            sd(MSE.EXPOS$Pc) / mean(MSE.EXPOS$Pc), sd(MSE.EXPOS$Pe) / mean(MSE.EXPOS$Pe), sd(MSE.EXPOS$H) / mean(MSE.EXPOS$H)))

## The ranges in the data. Each group's genes with MSE.MIN.MEAN counts per cell and MSE.MIN.CELLS cells with a count give the fitted rates
## sum(y) / sum(O) (the two hybrid alleles pooled, each with its share correction for phi); a random MSE.FIT.GENES of them per group are fitted
## with the Cox-Reid estimator for the dispersion quantiles.
MSE.DATA <- list(Pc = list(MIX.SC), Pe = list(MIX.SE), H = list(HYB.SC, HYB.SE))
MSE.ETA  <- list(Pc = 0, Pe = 0, H = c(R.C, R.E))
MU.OBS <- list(); PHI.OBS <- list()
for (k in MSE.GROUPS) for (a in seq_along(MSE.DATA[[k]])) {
  m <- MSE.DATA[[k]][[a]]; e <- MSE.EXPOS[[k]]
  ok <- rowMeans(m) >= MSE.MIN.MEAN & rowSums(m > 0) >= MSE.MIN.CELLS
  MU.OBS[[paste(k, a)]] <- rowSums(m[ok, , drop = FALSE]) / sum(e)
  gs <- with_local_seed(SEED.MSE + a, sample(which(ok), min(MSE.FIT.GENES, sum(ok))))
  ph <- vapply(gs, function(i) fit_dispersion(as.numeric(m[i, ]), e, MSE.ETA[[k]][[a]], "NB")[["disp"]], numeric(1))
  PHI.OBS[[paste(k, a)]] <- ph[is.finite(ph) & ph > 0]
}
MU.RANGE <- quantile(unlist(MU.OBS), MSE.MU.Q, names = FALSE)
MSE.MU   <- signif(exp(seq(log(MU.RANGE[1]), log(MU.RANGE[2]), length.out = MSE.N.MU)), 3)
MSE.PHI  <- signif(quantile(unlist(PHI.OBS), MSE.PHI.Q, names = FALSE), 3)
cat("mu grid (expected count at the reference depth):", MSE.MU, "\n")
cat("phi grid (quantiles", MSE.PHI.Q, "of the positive Cox-Reid fits):", MSE.PHI, "\n")
print(round(t(vapply(PHI.OBS, quantile, numeric(length(MSE.PHI.Q)), probs = MSE.PHI.Q)), 3))   # per group and allele, for reference

## Line colors and labels of the figures: one color per phi value in a single-hue ramp (light = small phi, dark = large phi).
MSE.COLS <- colorRampPalette(c("#6BAED6", "#08306B"))(length(MSE.PHI))
MSE.EST.LABEL <- c(NB = "Cox-Reid NB", weighted = "Weighted moments (Appendix A)")
MSE.GROUP.LABEL <- c(Pc = "Sc parent", Pe = "Se parent", H = "Hybrid allele")

##############################################################################
## M1. Bias, variance and MSE of phi in one group
if ("M1" %in% MSE.RUN) {
  ## The simulation grid, one task per cell, each on its own L'Ecuyer stream so the results do not depend on the core count.
  MSE.GRID <- expand.grid(group = MSE.GROUPS, family = MSE.FAMILIES, mu = MSE.MU, phi = MSE.PHI, stringsAsFactors = FALSE)
  MSE.STREAMS <- check_streams(nrow(MSE.GRID), SEED.MSE)
  t0 <- Sys.time()
  M1 <- do.call(rbind, mclapply(seq_len(nrow(MSE.GRID)), mse_cell, grid = MSE.GRID, expos = MSE.EXPOS, R = MSE.R, nbatch = MSE.NBATCH,
                                streams = MSE.STREAMS, mc.cores = CORES, mc.preschedule = FALSE))
  cat(sprintf("%d grid cells x %d datasets in %.1f min\n", nrow(MSE.GRID), MSE.R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ## Which estimator has the smaller MSE in each cell, and by how much (the MSE ratio Cox-Reid / moment, with the cell's z of the difference).
  M1.WIDE <- merge(M1[M1$estimator == "NB", c("group", "family", "mu", "phi", "mse", "mse_se", "bias", "var")],
                   M1[M1$estimator == "weighted", c("group", "family", "mu", "phi", "mse", "mse_se", "bias", "var")], by = c("group", "family", "mu", "phi"), suffixes = c("_cr", "_mom"))
  M1.WIDE$mse_ratio_cr_over_mom <- M1.WIDE$mse_cr / M1.WIDE$mse_mom
  M1.WIDE$z_diff <- (M1.WIDE$mse_cr - M1.WIDE$mse_mom) / sqrt(M1.WIDE$mse_se_cr^2 + M1.WIDE$mse_se_mom^2)
  M1.SUMMARY <- do.call(rbind, lapply(split(M1.WIDE, list(M1.WIDE$group, M1.WIDE$family), drop = TRUE), function(w) data.frame(group = w$group[1], family = w$family[1], cells = nrow(w), cr_better = sum(w$mse_ratio_cr_over_mom < 1 & w$z_diff < -2), mom_better = sum(w$mse_ratio_cr_over_mom > 1 & w$z_diff > 2), median_ratio = median(w$mse_ratio_cr_over_mom))))
  rownames(M1.SUMMARY) <- NULL
  print(M1.SUMMARY, digits = 3, row.names = FALSE)
  write.csv(M1, file.path(CHECK.DIR, "M1_estimator_mse.csv"), row.names = FALSE)
  write.csv(M1.WIDE, file.path(CHECK.DIR, "M1_estimator_mse_ratio.csv"), row.names = FALSE)
  check_save("M1", list(cells = M1, wide = M1.WIDE, summary = M1.SUMMARY, mu = MSE.MU, phi = MSE.PHI, mu_obs = MU.OBS, phi_obs = PHI.OBS), "run",
             headline = sprintf("%s %s: Cox-Reid MSE below the moment estimator's by more than 2 SE in %d of %d cells, above in %d; median MSE ratio CR / MoM %.3f",
                                M1.SUMMARY$group, M1.SUMMARY$family, M1.SUMMARY$cr_better, M1.SUMMARY$cells, M1.SUMMARY$mom_better, M1.SUMMARY$median_ratio),
             mc_error = sprintf("SE over %d batches of %d datasets per cell", MSE.NBATCH, MSE.R / MSE.NBATCH), script = CHECK.SCRIPT)

  ## Figures: one PDF per measure, one page per factor family, panels by group (rows) and estimator (columns). x is mu on a log scale; each phi value is one
  ## line in a single-hue ramp (light = small phi, dark = large phi), its mean solid and mean +- 1 SE in the same color, thinner and lighter. Panels of a row
  ## share the y range, so the two estimators of a group read on one scale; MSE and variance use a log y axis, bias a linear one with a zero line.
  for (ms in c("mse", "var", "bias")) {
    col_m <- if (MSE.PLOT.SCALE == "relative") paste0("rel_", ms) else ms; col_s <- paste0(col_m, "_se")
    ylab <- switch(ms, mse = if (MSE.PLOT.SCALE == "relative") "MSE / phi^2" else "MSE", var = if (MSE.PLOT.SCALE == "relative") "variance / phi^2" else "variance",
                   bias = if (MSE.PLOT.SCALE == "relative") "bias / phi" else "bias")
    fig_pdf(sprintf("extra/check_M1_%s.pdf", ms), 9, 2.6 * length(MSE.GROUPS) + 0.8)
    for (fam in MSE.FAMILIES) {
      layout(rbind(matrix(seq_len(2 * length(MSE.GROUPS)), ncol = 2, byrow = TRUE), 2 * length(MSE.GROUPS) + 1), heights = c(rep(1, length(MSE.GROUPS)), 0.25))
      par(mar = c(3.6, 4, 2.2, 1), mgp = c(2.2, 0.6, 0), cex.axis = 0.8)
      for (k in MSE.GROUPS) {
        d <- M1[M1$family == fam & M1$group == k, ]
        lo <- d[[col_m]] - d[[col_s]]; hi <- d[[col_m]] + d[[col_s]]
        yr <- if (ms == "bias") range(c(lo, hi, 0), na.rm = TRUE) else range(c(lo[lo > 0], hi), na.rm = TRUE)
        for (est in c("NB", "weighted")) {
          plot(NA, xlim = range(MSE.MU), ylim = yr, log = if (ms == "bias") "x" else "xy", xlab = "mu (expected count at the reference depth)", ylab = ylab,
               main = sprintf("%s, %s (%s factor)", MSE.GROUP.LABEL[[k]], MSE.EST.LABEL[[est]], fam), cex.main = 0.85, las = 1)
          abline(v = MSE.MU, col = "grey92", lwd = 0.5)
          if (ms == "bias") abline(h = 0, col = "grey55", lty = 2)
          for (j in seq_along(MSE.PHI)) {
            x <- d[d$estimator == est & d$phi == MSE.PHI[j], ]; x <- x[order(x$mu), ]
            m <- x[[col_m]]; s <- x[[col_s]]; low <- m - s; if (ms != "bias") low[low <= 0] <- NA
            lines(x$mu, low, col = adjustcolor(MSE.COLS[j], 0.35), lwd = 1)
            lines(x$mu, m + s, col = adjustcolor(MSE.COLS[j], 0.35), lwd = 1)
            lines(x$mu, m, col = MSE.COLS[j], lwd = 2)
          }
        }
      }
      par(mar = c(0, 0, 0, 0)); plot.new()
      legend("center", legend = sprintf("phi = %s", MSE.PHI), col = MSE.COLS, lwd = 2, horiz = TRUE, bty = "n", cex = 0.85,
             title = "solid: mean over batches; light: mean +- 1 SE")
    }
    dev.off()
  }
}

##############################################################################
## M2. The same comparison on the log2 contrasts the pipeline tests (Methods: Contrasts between groups)
# A bias that both groups of a contrast share cancels in their log2 ratio, so the estimator error that matters for the tests is the error of the contrast.
# Each grid cell simulates the two parents and the two alleles of every hybrid cell on the observed offsets (mse_contrast_cell()), with the alleles sharing
# an extrinsic factor of variance MSE.C.EXT.SHARE x the smaller allele dispersion. The hybrid cells are ordered by depth and split depth-balanced into cis
# and trans cells at MSE.C.F.CIS (split_indices_by_depth(), as Section 2.2). Three designs set the true log2 dispersion ratios: null (all equal), cis (the
# alleles and the parents differ by MSE.C.EFFECT, so cis = total = effect and trans = 0) and trans (only the parents differ, so trans = total = effect and
# cis = 0); in each, the reference dispersion is the grid phi on the Sc side. mu takes MSE.C.N.MU values of the M1 grid and phi the M1 values. The bias,
# variance and MSE of the log2 total, cis and trans contrasts are computed over the datasets in which the contrast is defined (every group phi > 0), with
# that share (defined); the pipeline additionally requires every group to be resolved in its bootstrap, which this check does not simulate.
if ("M2" %in% MSE.RUN) {
  MSE.C.R          <- 300                   # simulated datasets per cell
  MSE.C.N.MU       <- 5                     # mu values taken from MSE.MU, evenly spread over its range
  MSE.C.EFFECT     <- 1                     # log2 dispersion ratio of the cis and trans designs
  MSE.C.F.CIS      <- 0.5                   # fraction of the hybrid cells in the cis partition
  MSE.C.EXT.SHARE  <- 0.5                   # extrinsic variance shared by the two alleles of a cell, as a fraction of the smaller allele dispersion
  MSE.C.MU <- MSE.MU[unique(round(seq(1, length(MSE.MU), length.out = MSE.C.N.MU)))]
  MSE.C.HYB.SPLIT <- local({ ord <- order(MSE.EXPOS$H, decreasing = TRUE); sp <- split_indices_by_depth(length(ord), MSE.C.F.CIS); list(c = ord[sp$c], t = ord[sp$t]) })
  MSE.C.GRID <- rbind(expand.grid(family = MSE.FAMILIES, mu = MSE.C.MU, phi = MSE.PHI, design = "null", effect = 0, stringsAsFactors = FALSE),
                      expand.grid(family = MSE.FAMILIES, mu = MSE.C.MU, phi = MSE.PHI, design = c("cis", "trans"), effect = MSE.C.EFFECT, stringsAsFactors = FALSE))
  MSE.C.STREAMS <- check_streams(nrow(MSE.C.GRID), SEED.MSE + 1)
  t0 <- Sys.time()
  M2 <- do.call(rbind, mclapply(seq_len(nrow(MSE.C.GRID)), mse_contrast_cell, grid = MSE.C.GRID, expos = MSE.EXPOS, hyb_split = MSE.C.HYB.SPLIT, R = MSE.C.R,
                                nbatch = MSE.NBATCH, ext_share = MSE.C.EXT.SHARE, streams = MSE.C.STREAMS, mc.cores = CORES, mc.preschedule = FALSE))
  cat(sprintf("M2: %d grid cells x %d datasets in %.1f min\n", nrow(MSE.C.GRID), MSE.C.R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  ## Which estimator has the smaller contrast MSE in each cell (z of the difference), by contrast, family and design.
  M2.KEYS <- c("family", "mu", "phi", "design", "effect", "contrast")
  M2.WIDE <- merge(M2[M2$estimator == "NB", c(M2.KEYS, "mse", "mse_se", "bias", "var", "defined")], M2[M2$estimator == "weighted", c(M2.KEYS, "mse", "mse_se", "bias", "var", "defined")],
                   by = M2.KEYS, suffixes = c("_cr", "_mom"))
  M2.WIDE$mse_ratio_cr_over_mom <- M2.WIDE$mse_cr / M2.WIDE$mse_mom
  M2.WIDE$z_diff <- (M2.WIDE$mse_cr - M2.WIDE$mse_mom) / sqrt(M2.WIDE$mse_se_cr^2 + M2.WIDE$mse_se_mom^2)
  M2.SUMMARY <- do.call(rbind, lapply(split(M2.WIDE, list(M2.WIDE$contrast, M2.WIDE$family, M2.WIDE$design), drop = TRUE), function(w) data.frame(contrast = w$contrast[1], family = w$family[1], design = w$design[1], cells = nrow(w), cr_better = sum(w$z_diff < -2, na.rm = TRUE), mom_better = sum(w$z_diff > 2, na.rm = TRUE), median_ratio = median(w$mse_ratio_cr_over_mom, na.rm = TRUE))))
  M2.SUMMARY <- M2.SUMMARY[order(M2.SUMMARY$contrast, M2.SUMMARY$design, M2.SUMMARY$family), ]; rownames(M2.SUMMARY) <- NULL
  print(M2.SUMMARY, digits = 3, row.names = FALSE)
  write.csv(M2, file.path(CHECK.DIR, "M2_contrast_mse.csv"), row.names = FALSE)
  write.csv(M2.WIDE, file.path(CHECK.DIR, "M2_contrast_mse_ratio.csv"), row.names = FALSE)
  check_save("M2", list(cells = M2, wide = M2.WIDE, summary = M2.SUMMARY, mu = MSE.C.MU, phi = MSE.PHI, hyb_split = MSE.C.HYB.SPLIT,
                        settings = list(R = MSE.C.R, effect = MSE.C.EFFECT, f_cis = MSE.C.F.CIS, ext_share = MSE.C.EXT.SHARE)), "run",
             headline = sprintf("%s %s design, %s: Cox-Reid contrast MSE lower by > 2 SE in %d of %d cells, higher in %d; median MSE ratio CR / MoM %.3f",
                                M2.SUMMARY$contrast, M2.SUMMARY$design, M2.SUMMARY$family, M2.SUMMARY$cr_better, M2.SUMMARY$cells, M2.SUMMARY$mom_better, M2.SUMMARY$median_ratio),
             mc_error = sprintf("SE over %d batches of %d datasets per cell", MSE.NBATCH, MSE.C.R / MSE.NBATCH), script = CHECK.SCRIPT)

  ## Figures: one PDF per measure (figures/extra/check_M2_<measure>.pdf), one page per factor family and design, panels by contrast (rows: total, cis, trans) and
  ## estimator (columns). x is mu (log scale); each phi value is one line in the M1 ramp, the mean solid and mean +- 1 SE light. The measures are on the log2
  ## contrast scale: bias in log2 units, variance and MSE in squared log2 units. Panels of a row share the y range.
  for (ms in c("mse", "var", "bias")) {
    ylab <- switch(ms, mse = "MSE of the log2 contrast", var = "variance of the log2 contrast", bias = "bias of the log2 contrast")
    fig_pdf(sprintf("extra/check_M2_%s.pdf", ms), 9, 2.6 * 3 + 0.8)
    for (fam in MSE.FAMILIES) for (dz in c("null", "cis", "trans")) {
      layout(rbind(matrix(1:6, ncol = 2, byrow = TRUE), 7), heights = c(1, 1, 1, 0.25))
      par(mar = c(3.6, 4, 2.2, 1), mgp = c(2.2, 0.6, 0), cex.axis = 0.8)
      for (cn in c("total", "cis", "trans")) {
        d <- M2[M2$family == fam & M2$design == dz & M2$contrast == cn, ]
        lo <- d[[ms]] - d[[paste0(ms, "_se")]]; hi <- d[[ms]] + d[[paste0(ms, "_se")]]
        yr <- if (ms == "bias") range(c(lo, hi, 0), na.rm = TRUE) else range(c(lo[lo > 0], hi), na.rm = TRUE)
        for (est in c("NB", "weighted")) {
          plot(NA, xlim = range(MSE.C.MU), ylim = yr, log = if (ms == "bias") "x" else "xy", xlab = "mu (expected count at the reference depth)", ylab = ylab,
               main = sprintf("%s contrast, %s (%s factor, %s design)", cn, MSE.EST.LABEL[[est]], fam, dz), cex.main = 0.8, las = 1)
          abline(v = MSE.C.MU, col = "grey92", lwd = 0.5)
          if (ms == "bias") abline(h = 0, col = "grey55", lty = 2)
          for (j in seq_along(MSE.PHI)) {
            x <- d[d$estimator == est & d$phi == MSE.PHI[j], ]; x <- x[order(x$mu), ]
            m <- x[[ms]]; s <- x[[paste0(ms, "_se")]]; low <- m - s; if (ms != "bias") low[low <= 0] <- NA
            lines(x$mu, low, col = adjustcolor(MSE.COLS[j], 0.35), lwd = 1)
            lines(x$mu, m + s, col = adjustcolor(MSE.COLS[j], 0.35), lwd = 1)
            lines(x$mu, m, col = MSE.COLS[j], lwd = 2)
          }
        }
      }
      par(mar = c(0, 0, 0, 0)); plot.new()
      legend("center", legend = sprintf("phi = %s (Sc side)", MSE.PHI), col = MSE.COLS, lwd = 2, horiz = TRUE, bty = "n", cex = 0.8,
             title = "solid: mean over batches; light: mean +- 1 SE")
    }
    dev.off()
  }
}

##############################################################################
## M3. Rare cell states and the false-positive rate of the dispersion test (Methods: Null for the log2 dispersion contrasts, Informative estimates)
# A rare cell state (a few cells with expression J times the bulk) at the same frequency p in every group adds binomial noise in how many such cells each group
# happens to contain, for example 3 in one parent and 0 in the other. That noise is part of the true sampling variance of a dispersion contrast. The pipeline's
# shifted bootstrap resamples each group's own cells, so a group that sampled none of the rare cells never produces one in its resamples, and the null may be
# too narrow for total and trans. cis is expected to be protected, because the state belongs to the cell and lifts both alleles together. M3 measures it:
# every group shares mu, phi, p and J (all true contrasts zero), each dataset is tested with the pipeline's rule (shifted bootstrap null, resolution gate),
# and the rejection rate among tested datasets is the type I error (mse_rare_cell()). p = 0 is the baseline. sd_ratio is the SD of the observed contrast over
# datasets divided by the median SD of its null; above 1 means the null is too narrow. The rate among datasets whose parents' rare-cell counts differ by two
# or more isolates the case of concern. Runtime grows with M3.R x M3.B x 6 fits per estimator; on 8 cores the defaults take roughly half an hour.
if ("M3" %in% MSE.RUN) {
  M3.R     <- 200                         # datasets per cell
  M3.B     <- 199                         # shifted bootstrap draws per dataset
  M3.P     <- c(0, 0.001, 0.003, 0.01)    # rare-state frequency per cell
  M3.J     <- c(10, 50)                   # fold of expression in the rare state
  M3.MU    <- MSE.MU[unique(round(quantile(seq_along(MSE.MU), c(0.3, 0.7))))]   # two mu values of the M1 grid
  M3.PHI   <- MSE.PHI[unique(round(quantile(seq_along(MSE.PHI), c(0.5, 0.75))))] # bulk phi: the median and upper-quartile M1 values
  M3.EST   <- "NB"                        # estimators to test; add "weighted" for the moment estimator (doubles the runtime)
  M3.GRID <- rbind(expand.grid(mu = M3.MU, phi = M3.PHI, p = 0, J = 1), expand.grid(mu = M3.MU, phi = M3.PHI, p = M3.P[M3.P > 0], J = M3.J))
  M3.STREAMS <- check_streams(nrow(M3.GRID), SEED.MSE + 2)
  M3.HYB.SPLIT <- local({ ord <- order(MSE.EXPOS$H, decreasing = TRUE); sp <- split_indices_by_depth(length(ord), 0.5); list(c = ord[sp$c], t = ord[sp$t]) })
  t0 <- Sys.time()
  M3 <- do.call(rbind, mclapply(seq_len(nrow(M3.GRID)), mse_rare_cell, grid = M3.GRID, expos = MSE.EXPOS, hyb_split = M3.HYB.SPLIT, R = M3.R, B = M3.B,
                                ext_share = 0.5, estimators = M3.EST, streams = M3.STREAMS, mc.cores = CORES, mc.preschedule = FALSE))
  cat(sprintf("M3: %d grid cells x %d datasets x %d null draws in %.1f min\n", nrow(M3.GRID), M3.R, M3.B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(M3[, c("mu", "phi", "p", "J", "estimator", "contrast", "tested", "type1_05", "type1_05_se", "sd_ratio", "rare_split_share", "type1_05_rare_split")], digits = 3, row.names = FALSE)
  write.csv(M3, file.path(CHECK.DIR, "M3_rare_state_type1.csv"), row.names = FALSE)
  M3.FLAG <- M3[M3$type1_05 - 0.05 > 2 * M3$type1_05_se & is.finite(M3$type1_05_se), ]
  check_save("M3", list(cells = M3, grid = M3.GRID, settings = list(R = M3.R, B = M3.B, ext_share = 0.5, f_cis = 0.5)), "run",
             headline = c(sprintf("%s %s, p %g, J %g, mu %s, phi %s: tested %.2f, type I %.3f (SE %.3f), SD ratio %.2f", M3$estimator, M3$contrast, M3$p, M3$J, M3$mu, M3$phi,
                                  M3$tested, M3$type1_05, M3$type1_05_se, M3$sd_ratio),
                          sprintf("type I above 0.05 by more than 2 SE in %d of %d cells", nrow(M3.FLAG), nrow(M3))),
             mc_error = sprintf("binomial SE over the tested datasets of %d per cell; %d null draws per dataset", M3.R, M3.B), script = CHECK.SCRIPT)

  ## Figure (figures/extra/check_M3_type1.pdf): page 1 the type I error at 0.05 against the rare-state frequency p, one line per fold J (single-hue ramp), the
  ## mean solid and +- 1 SE light, the nominal 0.05 dashed; panels by contrast (rows) and (mu, phi) (columns). p is drawn at equal spacing with p = 0 first.
  ## Page 2 the SD ratio of the observed contrast to its null on the same layout (1 dashed).
  M3.COLS <- colorRampPalette(c("#6BAED6", "#08306B"))(length(M3.J))
  M3.CELLS <- unique(M3.GRID[, c("mu", "phi")])
  fig_pdf("extra/check_M3_type1.pdf", 3 * nrow(M3.CELLS) + 0.5, 2.6 * 3 + 0.8)
  for (ms in c("type1_05", "sd_ratio")) for (est in M3.EST) {
    layout(rbind(matrix(seq_len(3 * nrow(M3.CELLS)), nrow = 3, byrow = TRUE), 3 * nrow(M3.CELLS) + 1), heights = c(1, 1, 1, 0.25))
    par(mar = c(3.6, 4, 2.2, 1), mgp = c(2.2, 0.6, 0), cex.axis = 0.8)
    for (cn in c("total", "cis", "trans")) for (q in seq_len(nrow(M3.CELLS))) {
      d <- M3[M3$estimator == est & M3$contrast == cn & M3$mu == M3.CELLS$mu[q] & M3$phi == M3.CELLS$phi[q], ]
      base <- d[d$p == 0, ]
      yr <- if (ms == "type1_05") range(c(0, 0.1, d$type1_05 + d$type1_05_se), na.rm = TRUE) else range(c(0.8, 1.2, d$sd_ratio), na.rm = TRUE)
      plot(NA, xlim = c(1, length(M3.P)), ylim = yr, xaxt = "n", xlab = "rare-state frequency p", ylab = if (ms == "type1_05") "type I error at 0.05" else "SD(observed) / SD(null)",
           main = sprintf("%s, mu %s, phi %s (%s)", cn, M3.CELLS$mu[q], M3.CELLS$phi[q], est), cex.main = 0.85, las = 1)
      axis(1, at = seq_along(M3.P), labels = M3.P)
      abline(h = if (ms == "type1_05") 0.05 else 1, col = "grey55", lty = 2)
      for (j in seq_along(M3.J)) {
        x <- rbind(base, d[d$J == M3.J[j], ]); x <- x[order(x$p), ]; xp <- match(x$p, M3.P)
        if (ms == "type1_05") {
          lines(xp, x$type1_05 - x$type1_05_se, col = adjustcolor(M3.COLS[j], 0.35), lwd = 1)
          lines(xp, x$type1_05 + x$type1_05_se, col = adjustcolor(M3.COLS[j], 0.35), lwd = 1)
        }
        lines(xp, x[[ms]], col = M3.COLS[j], lwd = 2); points(xp, x[[ms]], col = M3.COLS[j], pch = 19, cex = 0.6)
      }
    }
    par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", legend = sprintf("J = %s", M3.J), col = M3.COLS, lwd = 2, horiz = TRUE, bty = "n", cex = 0.85,
           title = if (ms == "type1_05") "solid: rejection rate among tested datasets; light: +- 1 binomial SE" else "SD ratio of the observed contrast to its shifted bootstrap null")
  }
  dev.off()
}
