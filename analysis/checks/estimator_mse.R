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
### Runtime: about MSE.R x 2 fits per cell; at MSE.R = 500 on one core expect roughly 30 to 60 minutes for the default grid.
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
MSE.COLS <- colorRampPalette(c("#6BAED6", "#08306B"))(length(MSE.PHI))
MSE.EST.LABEL <- c(NB = "Cox-Reid NB", weighted = "Weighted moments (Appendix A)")
MSE.GROUP.LABEL <- c(Pc = "Sc parent", Pe = "Se parent", H = "Hybrid allele")
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
