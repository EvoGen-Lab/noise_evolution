###############################################################
### power_analysis.R
### Section 10 of the pipeline: power of the permutation test for
### a burst-frequency (NB size) difference between two groups,
### simulated over mean expression, cell count, burst frequency and
### size ratio.
###
### Run after analysis.R has written the Section 2 checkpoint
### (section2_checkpoint.rda, which supplies PR and N.PERM).
### Functions live in R/functions_power.R (this analysis and the
### cluster script power_grid.R) and R/functions.R (shared NB fit).
### Cluster round trip: this script writes power_inputs.rda (10.1),
### power.sub runs power_grid.R, and 10.2 onward loads the
### power_output_<k>of<K>.rda files. Sections 10.9 to 10.14 add the cis
### and trans contrasts on paired-allele hybrid data (power_modes_inputs.rda,
### power_modes.sub, power_modes_grid.R, power_modes_output_<k>of<K>.rda).
###############################################################

# here() anchors every path to the project root, so R/setup.R loads from any working
# directory. setup.R checks renv, defines the project paths, class names and color
# palettes, and sources R/functions.R.
library(here)
source(here("R", "setup.R"))
source(here("R", "functions_power.R"))

# PR (per-gene permutation results with BH q-values) and N.PERM come from the
# Section 2 checkpoint
SEC2 <- new.env()
load(ckpt_path(2), envir = SEC2)
if (!all(c("PR", "N.PERM") %in% ls(SEC2)))
  stop("The Section 2 checkpoint lacks PR or N.PERM. Checkpoints written before the power ",
       "analysis moved out of analysis.R do not hold N.PERM; add it with ",
       "load(ckpt_path(2)); N.PERM <- <value used in analysis.R Section 2>; save(list = ls(), file = ckpt_path(2)) ",
       "or rerun the Section 2 checkpoint save.")
PR     <- SEC2$PR
N.PERM <- SEC2$N.PERM
# SEC2 stays loaded for Section 10.9, which reads the real design parameters

console_start(10)
set.seed(SEED.SECTION + 10)
##############################################################################
## 10. POWER ANALYSIS                                                       ##
##############################################################################
# Run through 10.1 to build power_inputs.rda, then submit
# power.sub (a SLURM job array) to the cluster. Once every power_output_<k>of<K>.rda
# exists, run 10.2 onward to reshape them back into POWER and reproduce every plot.
# PR (permutation results with BH q-values) and N.PERM come from the Section 2
# checkpoint: PR sets PI1, the fraction of truly different genes in the simulated
# mixture, and N.PERM matches the permutation count to the gene-level analysis.


# Significance level for power analysis. Power is judged at this Benjamini-Hochberg
# FDR (q < ALPHA), matching the gene-level analysis.
ALPHA <- 0.05
# Average number of reads per cell, before per-cell capture-depth scaling.
MEAN.READS <- c(0.25,0.5,1,2,4,8,16,32,64,128)
# Burst frequency (NB SIZE / dispersion parameter, theta) for the
# reference group.
SIZE <- c(0.5,1,2,4,8,16,32,64,128,256)
# Ratio of SIZE between the two groups. 1 confirms the false positive
# rate is controlled.
SIZE.RATIO <- c(1,1.1,1.25,1.5,1.75,2,3,4,6,8,12)

# Reference-group (Sc) cell count, evenly spaced.
N.CELLS <- c(400, 800, 1600, 2400)
# Se:Sc cell-count ratio. Set to 1, equal group sizes.
CELL.RATIO <- 1

# Per-cell capture-depth coefficient of variation.
EXPOSURE.CV <- 0.8

NJ <- 1000 # Simulated datasets per row
NI <- N.PERM # Permutations per dataset, matched to the gene-level analysis
PI1 <- min(0.5, max(0.02, mean(PR$bfreq_total_q < ALPHA, na.rm = TRUE)))
cat(sprintf("PI1 (fraction of genes called by BH at q < %.2f, parental burst frequency) = %.3f\n", ALPHA, PI1))
N.MIX <- 50 # Random mixtures averaged per grid cell
SEED.BASE <- 1

## 10.1 Grid definition and cluster submission
# Each row covers one (MEAN.READS, N.CELLS, SIZE) combination and returns
# power for every SIZE.RATIO value at once
M <- 1:length(MEAN.READS); N <- 1:length(N.CELLS)
P <- 1:length(SIZE);       Q <- 1:length(SIZE.RATIO)

GRID <- expand.grid(m = M, n = N, p = P)

## ---- Cluster round trip: SLURM job array power.sub (Rscript power_grid.R, one task per array index) ----
## Reads power_inputs.rda (saved below), writes one power_output_<k>of<K>.rda
## per array task (ROWS.PART, POWER.PART), reshaped into POWER in 10.2
save(GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE,
     file = file.path(INPUT.DIR, "power_inputs.rda"))

cat(sprintf("power grid inputs written: %d rows x %d SIZE.RATIO values\n", nrow(GRID), length(SIZE.RATIO)))
cat("Submit power.sub (job array) now; resume below once every power_output_<k>of<K>.rda exists.\n")
## ---- end cluster round trip (power_output_<k>of<K>.rda files are loaded in 10.2) ----

## 10.2 Reshape cluster output into the POWER array
# Each array task writes the GRID rows it computed (ROWS.PART) and their power values
# (POWER.PART). Rows are placed back in GRID order, so the result does not depend on
# how many array tasks were used.
PART.FILES <- list.files(OUTPUT.DIR, pattern = "^power_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(PART.FILES) > 0)
POWER.MAT <- matrix(NA_real_, nrow(GRID), length(SIZE.RATIO))
ROW.DONE  <- logical(nrow(GRID))
for (f in PART.FILES) {
  e <- new.env(); load(f, envir = e)
  POWER.MAT[e$ROWS.PART, ] <- e$POWER.PART
  ROW.DONE[e$ROWS.PART]    <- TRUE
}
stopifnot(all(ROW.DONE))   # every array task must have finished before reshaping

POWER <- array(NA_real_, dim = c(length(M), length(N), length(P), length(Q)))
for (i in seq_len(nrow(GRID))) {
  POWER[GRID$m[i], GRID$n[i], GRID$p[i], ] <- POWER.MAT[i, ]
}
save(POWER, file = file.path(OUTPUT.DIR, "power.rda"))

## Value of each POWER dimension: mean reads (M), cell count (N), burst-frequency SIZE (P), SIZE.RATIO (Q)
POWER.GRID <- list(M = MEAN.READS, N = N.CELLS, P = SIZE, Q = SIZE.RATIO)
LOG2.READS.AXIS <- list(at = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), labels = -2:10, mgp = c(3, 0.4, 0))

## 10.3 Figure 1: power vs mean expression, panels over N.CELLS x SIZE.RATIO
power_line_figure("Power.Analysis.grid_mean_by_ncells.pdf", POWER, POWER.GRID, "N", "M", x = "Q", series = "P",
  main_fmt = "Mean=%2$s Sc N.Cells=%1$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  xlab = "SIZE Ratio", ticks = TRUE)
power_line_figure("Power.Analysis.1.pdf", POWER, POWER.GRID, "N", "Q", x = "M", series = "P",
  main_fmt = "SIZE Ratio = %2$s\nSc # Cells = %1$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  styled = c(LOG2.READS.AXIS, list(cex_axis = 0.6, xlab = expression('log'[2]*'(Average Reads per Cell)'), xlab_mgp = c(1.2, 1, 0), xlab_cex = 0.6)))

## 10.4 Supplementary diagnostic grids: every other pairing of the four
# grid axes (mean, burst frequency, SIZE.RATIO, N.CELLS)
power_line_figure("Power.Analysis.grid_ncells_by_mean.pdf", POWER, POWER.GRID, "Q", "M", x = "N", series = "P",
  main_fmt = "Ratio=%1$s Mean=%2$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  log_x = FALSE, xlab = "Sc Cell Count")
power_line_figure("Power.Analysis.grid_logmean_by_size_ratio.pdf", POWER, POWER.GRID, "P", "Q", x = "M", series = "N",
  main_fmt = "SIZE=%1$s Ratio=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = "log2(Mean Reads)")
power_line_figure("Power.Analysis.grid_size_by_ratio.pdf", POWER, POWER.GRID, "M", "Q", x = "P", series = "N",
  main_fmt = "MEAN=%1$s Ratio=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = expression('log'[2]*'(SIZE)'))
power_line_figure("Power.Analysis.grid_sizeratio_by_mean.pdf", POWER, POWER.GRID, "P", "M", x = "Q", series = "N",
  main_fmt = "SIZE=%1$s Mean=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = "SIZE.RATIO", ticks = TRUE)
power_line_figure("Power.Analysis.grid_sizeratio_by_ncells_full.pdf", POWER, POWER.GRID, "P", "N", x = "Q", series = "M",
  main_fmt = "SIZE=%1$s Sc Cell Number=%2$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  xlab = "SIZE.RATIO", ticks = TRUE)
power_line_figure("Power.Analysis.grid_ncells_by_size_ratio_full.pdf", POWER, POWER.GRID, "P", "Q", x = "N", series = "M",
  main_fmt = "SIZE=%1$s Ratio=%2$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  log_x = FALSE, xlab = "Sc Cell Count")

## 10.5 Figure 2: power vs burst frequency, panels over N.CELLS x SIZE.RATIO
power_line_figure("Power.Analysis.2.pdf", POWER, POWER.GRID, "N", "Q", x = "P", series = "M", mar = c(4.5, 3, 2, 1),
  main_fmt = "SIZE Ratio = %2$s\nSc # Cells = %1$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  styled = list(at = log2(SIZE), labels = SIZE, cex_axis = 0.55, mgp = c(3, 0.7, 0), las = 2, tcl = -0.3,
                xlab = "Burst Frequency (SIZE)", xlab_mgp = c(3.3, 1, 0), xlab_cex = 0.6))
power_line_figure("Power.Analysis.grid_size_by_ncells_mean.pdf", POWER, POWER.GRID, "N", "M", x = "P", series = "Q",
  main_fmt = "Sc Cell Number=%1$s Mean=%2$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  xlab = expression('log'[2]*'(SIZE)'))
power_line_figure("Power.Analysis.grid_ncells_by_size_mean.pdf", POWER, POWER.GRID, "P", "M", x = "N", series = "Q",
  main_fmt = "SIZE=%1$s Mean=%2$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  log_x = FALSE, xlab = "Sc Cell Count")

## 10.6 Figure 3: power vs mean expression, panels over N.CELLS x SIZE
power_line_figure("Power.Analysis.3.pdf", POWER, POWER.GRID, "N", "P", x = "M", series = "Q",
  main_fmt = "SIZE = %2$s\nSc # Cells = %1$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  styled = c(LOG2.READS.AXIS, list(cex_axis = 0.55, xlab = expression('log'[2]*'(Average Reads per Cell)'), xlab_mgp = c(1.4, 1, 0), xlab_cex = 0.8)))

## 10.7 Heatmap A: power over mean expression x cell count, faceted by burst frequency
HM.Q <- which.min(abs(SIZE.RATIO - 2))
HM.COLS <- COLOR.SEQ

open_grid_pdf(file.path(FIGURE.DIR, "extended", "Power.Analysis.heatmap_mean_by_ncells.pdf"), ceiling(length(P)/2), 2, panel_w = 3.4, panel_h = 3.2)
for (p in P) {
  mat <- POWER[,,p,HM.Q]  # length(M) x length(N)
  image(log2(MEAN.READS), N.CELLS, mat, zlim = c(0,1), col = HM.COLS,
        xaxt = "n", yaxt = "n", xlab = "Mean Reads per Cell", ylab = "Sc Cell Count",
        main = paste0("SIZE = ", SIZE[p]), cex.main = 0.85)
  axis(1, at = log2(MEAN.READS), labels = MEAN.READS, las = 2, cex.axis = 0.55)
  axis(2, at = N.CELLS, labels = N.CELLS, cex.axis = 0.65)
  if (any(mat >= 0.8, na.rm = TRUE) && any(mat < 0.8, na.rm = TRUE)) {
    contour(log2(MEAN.READS), N.CELLS, mat, levels = 0.8, add = TRUE, col = "white", lwd = 1.5, labcex = 0.6)
  }
}
dev.off()

fig_pdf(file.path("extended", "Power.Analysis.heatmap_colorbar.pdf"), 4, 1.3, useDingbats = TRUE)
par(mar = c(2.2,1,1,1))
image(seq(0,1,length.out=100), 1, matrix(seq(0,1,length.out=100), ncol=1), col = HM.COLS,
      axes = FALSE, xlab = paste0("Power to detect SIZE.RATIO = ", round(SIZE.RATIO[HM.Q],2)), ylab = "")
axis(1, at = c(0,0.2,0.4,0.6,0.8,1.0), labels = c("0%","20%","40%","60%","80%","100%"))
dev.off()

## 10.8 Heatmap B: minimum detectable SIZE.RATIO, faceted by burst frequency

## Minimum detectable SIZE.RATIO: the smallest SIZE.RATIO (ratios in ascending order) at which power
## reaches 80%, one number per (MEAN.READS, N.CELLS, SIZE) in place of the full power-vs-ratio curve.
## Lower means more sensitive. NA when even the largest tested ratio misses the target, which image()
## renders as blank.
MDR <- array(NA_real_, dim = c(length(M), length(N), length(P)))
for (m in M) for (n in N) for (p in P) MDR[m,n,p] <- local({
  power_vec <- POWER[m, n, p, ]
  ratios <- SIZE.RATIO
  target <- 0.8
  hit <- which(power_vec >= target)
  if (length(hit) == 0) return(NA_real_)
  ratios[min(hit)]
})

MDR.COLS <- rev(COLOR.SEQ)
MDR.ZLIM <- log2(range(SIZE.RATIO))

open_grid_pdf(file.path(FIGURE.DIR, "extended", "Power.Analysis.min_detectable_ratio.pdf"), ceiling(length(P)/2), 2, panel_w = 3.4, panel_h = 3.2)
for (p in P) {
  mat <- log2(MDR[,,p])  # NA stays NA
  image(log2(MEAN.READS), N.CELLS, mat, zlim = MDR.ZLIM, col = MDR.COLS,
        xaxt = "n", yaxt = "n", xlab = "Mean Reads per Cell", ylab = "Sc Cell Count",
        main = paste0("SIZE = ", SIZE[p], "\n(blank = not detectable by ratio=12)"), cex.main = 0.7)
  axis(1, at = log2(MEAN.READS), labels = MEAN.READS, las = 2, cex.axis = 0.55)
  axis(2, at = N.CELLS, labels = N.CELLS, cex.axis = 0.65)
}
dev.off()

fig_pdf(file.path("extended", "Power.Analysis.min_detectable_ratio_colorbar.pdf"), 4, 1.3, useDingbats = TRUE)
par(mar = c(2.2,1,1,1))
zseq <- seq(MDR.ZLIM[1], MDR.ZLIM[2], length.out = 100)
image(zseq, 1, matrix(zseq, ncol=1), col = MDR.COLS,
      axes = FALSE, xlab = "Minimum detectable SIZE.RATIO (80% power)", ylab = "")
RATIO.TICKS <- SIZE.RATIO
axis(1, at = log2(RATIO.TICKS), labels = RATIO.TICKS, cex.axis = 0.7)
dev.off()

# Checkpoint
save(GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE,
     POWER, MDR, file = ckpt_path(10))

console_stop()

##############################################################################
## 10.9 CIS AND TRANS POWER ON PAIRED-ALLELE DATA                           ##
##############################################################################
# Sections 10.1 to 10.8 simulate one two-group contrast. The pipeline tests two others:
# cis (the alleles of the same hybrid cells, nulled by the within-cell allele swap) and
# trans (the parental ratio minus the hybrid allele ratio, nulled by pooling parent and
# hybrid cells). Sections 10.9 to 10.14 simulate each with the data structure and the null the
# pipeline uses (power_modes_row() in R/functions_power.R, null helpers in R/functions.R), so the
# modes are directly comparable at the same effect and depth. The effect stays SIZE.RATIO, a
# ratio of NB sizes (burst frequency), in log2 units:
#   null          ratio 1 everywhere
#   cis effect    hybrid allele ratio = ratio, parental ratio = ratio (trans = 0)
#   trans effect  hybrid allele ratio = 1, parental ratio = ratio (cis = 0)
# Hybrid cells carry two alleles that share an extrinsic gamma factor u (CV^2 = phi / SIZE) and each have
# an intrinsic gamma factor v, so the alleles' latent correlation is phi and each allele has latent CV^2
# 1 / SIZE; phi = 0 gives unpaired alleles and, as a control, cis and trans should then follow the
# total result. The counts are Poisson, so the marginal law is close to negative binomial; the check below
# confirms the recovered size. The trans null runs in two versions on the same data: the independent
# shuffles that permute_contrasts_one() uses and a pairing-preserving swap (perm_trans_pool_paired()),
# so the SIZE.RATIO == 1 type I error validates the pairing fix.

# Real design parameters from the Section 2 checkpoint: hybrid and parent cell counts, the split
# fraction f* the pipeline chose, the hybrid-to-parent per-allele depth ratio, and the allele
# correlation across genes that sets the phi grid.
MODES.REAL <- local({
  hyb_sc <- SEC2$HYB.SC; hyb_se <- SEC2$HYB.SE; genes <- SEC2$GENES; fits <- SEC2$CONTRAST.FITS
  ## Per-allele depth ratio: the hybrid allele's share of the hybrid library times the hybrid-to-parent
  ## library ratio, averaged over the two species
  share_sc <- sum(hyb_sc) / (sum(hyb_sc) + sum(hyb_se))
  depth <- mean(c(share_sc * mean(SEC2$EXPO.HYB) / mean(SEC2$EXPO.MIX.SC),
                  (1 - share_sc) * mean(SEC2$EXPO.HYB) / mean(SEC2$EXPO.MIX.SE)))
  ## Allele-pair residual correlation per gene, disattenuated by the alleles' reliabilities as in Section 5.1
  ## (the latent correlation phi of the generator is the quantity this estimates)
  r_sc <- nb_residuals(hyb_sc[genes, ], SEC2$EXPO.HYB, fits$HYB.SC[genes, ])
  r_se <- nb_residuals(hyb_se[genes, ], SEC2$EXPO.HYB, fits$HYB.SE[genes, ])
  attn <- sqrt(gene_reliability(fits["HYB.SC"], genes) * gene_reliability(fits["HYB.SE"], genes))
  rho_true <- row_cor(r_sc, r_se) / attn
  rho_true <- rho_true[is.finite(rho_true) & attn > 0.2]
  list(N.H = ncol(hyb_sc), N.P = round(mean(c(length(SEC2$EXPO.MIX.SC), length(SEC2$EXPO.MIX.SE)))),
       depth = depth, rho_true = rho_true,
       f_star = if (is.null(SEC2$SPLIT.FRAC)) c(f_mean = NA_real_, f_disp = NA_real_) else unlist(SEC2$SPLIT.FRAC[c("f_mean", "f_disp")]))
})
rm(SEC2)
cat(sprintf("real design: N.H = %d hybrid cells, N.P = %d parent cells per species, hybrid-to-parent per-allele depth ratio %.2f\n",
            MODES.REAL$N.H, MODES.REAL$N.P, MODES.REAL$depth))
cat("disattenuated allele-pair correlation across genes (quantiles):\n"); print(round(quantile(MODES.REAL$rho_true, c(0.1, 0.25, 0.5, 0.75, 0.9)), 3))
cat("pipeline split fractions (f_mean, f_disp):", round(MODES.REAL$f_star, 3), "\n")

# Grid. The reduced grid below is a first pass: a subset of mean reads, cell count and size with
# the full phi and f sweeps; extend MEAN.READS.M, SIZE.M, MODES.N.IDX, NJ.M and NI.M for the full grid.
# The phi grid is the control (0) and the 10th, 50th and 90th percentiles of the disattenuated allele
# correlation across genes. The map from that correlation to phi is approximate: the Pearson-residual
# correlation is attenuated by the alleles' reliabilities and the product-of-gammas model is only
# approximately negative binomial, so phi is a latent-correlation scale, not a fitted parameter.
# The f grid spans the split fractions and adds the pipeline's own f_mean and f_disp. The depth grid is
# the observed hybrid-to-parent per-allele depth ratio.
MEAN.READS.M <- MEAN.READS[c(5, 7)]                 # 4 and 16 reads per cell
SIZE.M       <- SIZE[c(3, 5)]                       # 2 and 8
MODES.N.IDX  <- 1                                   # first N.CELLS value (400 parent cells)
DESIGN.M     <- data.frame(N.P = N.CELLS[MODES.N.IDX],
                           N.H = round(N.CELLS[MODES.N.IDX] * MODES.REAL$N.H / MODES.REAL$N.P))
PHI.M        <- round(c(0, pmin(0.95, pmax(0, quantile(MODES.REAL$rho_true, c(0.1, 0.5, 0.9))))), 3)
FRAC.M       <- sort(unique(round(c(0.1, 0.2, 0.3, 0.4, 0.5, MODES.REAL$f_star[is.finite(MODES.REAL$f_star)]), 3)))
DEPTH.RATIO.M <- round(MODES.REAL$depth, 3)
SIZE.RATIO.M <- c(1, 1.25, 1.5, 2, 3, 4, 6, 8)
NJ.M <- 200    # simulated datasets per row
NI.M <- 2000   # permutations per dataset (the pipeline uses N.PERM)
# BH power needs the permutation p-values to resolve small values: the smallest p-value is 1 / (NI.M + 1), and BH
# over the pooled null and alternative datasets rejects only when enough alternatives sit at or below
# ALPHA * rank / total. A small NI.M therefore makes the BH power curve a step that stays at zero until the raw power
# is high. The number of alternatives that must reach the p-value floor is checked here, and the tables fall back
# to the raw rejection rate (p < ALPHA) when it is large.
BH.NEED <- ceiling((NJ.M + round(NJ.M * PI1 / (1 - PI1))) * (1 / (NI.M + 1)) / ALPHA)
MODES.PW.COL <- if (BH.NEED <= 0.1 * round(NJ.M * PI1 / (1 - PI1))) "power_bh" else "power_raw"
cat(sprintf("BH needs %d alternatives at the p-value floor (of %d); tables and curves use %s\n", BH.NEED, round(NJ.M * PI1 / (1 - PI1)), MODES.PW.COL))
GRID.M <- expand.grid(m = seq_along(MEAN.READS.M), p = seq_along(SIZE.M), h = seq_along(PHI.M),
                      f = seq_along(FRAC.M), n = seq_len(nrow(DESIGN.M)), d = seq_along(DEPTH.RATIO.M))

# Generator check on the grid's (size, phi) values at a high mean: stops unless the recovered NB size is
# within 10% of the target and the latent allele correlation within 0.03 of phi
set.seed(SEED.BASE)
GEN.GRID  <- expand.grid(size = SIZE.M, phi = PHI.M)
GEN.CHECK <- do.call(rbind, Map(paired_alleles_check, size = GEN.GRID$size, phi = GEN.GRID$phi, MoreArgs = list(mu = max(MEAN.READS.M))))
print(round(GEN.CHECK, 3))

MODES.INPUTS <- list(GRID = GRID.M, MEAN.READS = MEAN.READS.M, SIZE = SIZE.M, PHI = PHI.M, FRAC = FRAC.M,
                     DESIGN = DESIGN.M, DEPTH.RATIO = DEPTH.RATIO.M, SIZE.RATIO = SIZE.RATIO.M,
                     EXPOSURE.CV = EXPOSURE.CV, CELL.RATIO = CELL.RATIO, ALPHA = ALPHA, NJ = NJ.M, NI = NI.M,
                     PI1 = PI1, N.MIX = N.MIX, SEED.BASE = SEED.BASE)
## ---- Cluster round trip: SLURM job array power_modes.sub (Rscript power_modes_grid.R, one task per array index) ----
## Reads power_modes_inputs.rda (saved below), writes one power_modes_output_<k>of<K>.rda per array task
## (ROWS.PART, RESULTS.PART), reassembled in 10.10
save(MODES.INPUTS, file = file.path(INPUT.DIR, "power_modes_inputs.rda"))
cat(sprintf("power modes inputs written: %d rows x %d SIZE.RATIO values x 3 modes\n", nrow(GRID.M), length(SIZE.RATIO.M)))
cat("Submit power_modes.sub (job array) now; resume below once every power_modes_output_<k>of<K>.rda exists.\n")
## ---- end cluster round trip ----

## 10.10 Reassemble the cluster output
MODES.FILES <- list.files(OUTPUT.DIR, pattern = "^power_modes_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(MODES.FILES) > 0)
MODES.RES <- vector("list", nrow(GRID.M))
for (fl in MODES.FILES) {
  e <- new.env(); load(fl, envir = e)
  MODES.RES[e$ROWS.PART] <- e$RESULTS.PART
}
stopifnot(!any(vapply(MODES.RES, is.null, logical(1))))   # every array task must have finished
MODES <- modes_long(MODES.RES, GRID.M, MEAN.READS.M, SIZE.M, PHI.M, FRAC.M, DESIGN.M, DEPTH.RATIO.M)
MODES.POWER <- MODES$power
MODES.SD    <- MODES$sd_null
MODES.COL   <- c(cis = SPECIES.COLOR[["Sc"]], trans_indep = SPECIES.COLOR[["Se"]], trans_paired = COLOR.ACCENT)
MODES.LTY   <- c(cis = 1, trans_indep = 1, trans_paired = 2)

## 10.11 Power against effect size, by mode, at every grid point (size axis)
# One PDF per (mean reads, size): panels over f (rows) and phi (columns), one curve per mode, power at
# BH q < ALPHA (or the raw rate when MODES.PW.COL is power_raw) against the log2 SIZE.RATIO, with the 80% line.
for (mr in MEAN.READS.M) for (sz in SIZE.M) {
  open_grid_pdf(file.path(FIGURE.DIR, "extended", sprintf("Power.Modes.curves_mean%s_size%s.pdf", mr, sz)),
                length(FRAC.M), length(PHI.M), panel_w = 2.6, panel_h = 2.2)
  for (fr in FRAC.M) for (ph in PHI.M) {
    d <- MODES.POWER[MODES.POWER$mean_reads == mr & MODES.POWER$size == sz & MODES.POWER$f == fr & MODES.POWER$phi == ph & MODES.POWER$axis == "disp", ]
    plot(NA, xlim = range(log2(SIZE.RATIO.M)), ylim = c(0, 1), xlab = "log2 SIZE.RATIO", ylab = "Power",
         main = sprintf("f = %s, phi = %s", fr, ph), cex.main = 0.8)
    abline(h = 0.8, col = COLOR.GREY[["mid"]], lty = 3)
    for (md in names(MODES.COL)) with(d[d$mode == md, ], lines(log2(ratio), get(MODES.PW.COL), col = MODES.COL[[md]], lty = MODES.LTY[[md]], lwd = 1.6))
  }
  dev.off()
}
fig_pdf("extended/Power.Modes.legend.pdf", 3, 2)
plot.new(); legend("center", legend = c("cis", "trans, independent shuffles", "trans, pairing-preserving"), col = MODES.COL, lty = MODES.LTY, lwd = 2, bty = "n")
dev.off()

## 10.12 Minimum detectable log2 ratio at 80% power, and the trans-to-cis ratio (open item 4)
# mde_interp() interpolates each power curve (column MODES.PW.COL) over log2(SIZE.RATIO); NA means the curve never reaches 80% by the
# largest ratio tested. MDE ratio > 1 means trans needs a larger effect than cis at the same depth.
MODES.MDE <- modes_mde(MODES.POWER[MODES.POWER$axis == "disp", ], target = 0.8, col = MODES.PW.COL)
MODES.MDE.WIDE <- reshape(MODES.MDE[, c("i", "mean_reads", "size", "phi", "f", "N.P", "N.H", "depth", "mode", "mde_log2")],
                          idvar = c("i", "mean_reads", "size", "phi", "f", "N.P", "N.H", "depth"), timevar = "mode", direction = "wide")
write.csv(MODES.MDE, file.path(TABLE.DIR, "power_modes_mde.csv"), row.names = FALSE)
MDE.RATIO <- do.call(rbind, lapply(c("trans_indep", "trans_paired"), function(tm) {
  w <- MODES.MDE.WIDE
  data.frame(w[c("mean_reads", "size", "phi", "f", "N.P", "N.H", "depth")], trans_mode = tm,
             mde_cis = w[["mde_log2.cis"]], mde_trans = w[[paste0("mde_log2.", tm)]],
             mde_ratio = w[[paste0("mde_log2.", tm)]] / w[["mde_log2.cis"]])
}))
write.csv(MDE.RATIO, file.path(TABLE.DIR, "power_modes_mde_ratio.csv"), row.names = FALSE)
cat("MDE_trans / MDE_cis (independent-shuffle trans null), median over the grid by phi and f:\n")
print(round(tapply(MDE.RATIO$mde_ratio[MDE.RATIO$trans_mode == "trans_indep"],
                   list(phi = MDE.RATIO$phi[MDE.RATIO$trans_mode == "trans_indep"], f = MDE.RATIO$f[MDE.RATIO$trans_mode == "trans_indep"]),
                   median, na.rm = TRUE), 2))
fig_pdf("extended/Power.Modes.mde_ratio_by_f.pdf", 3.2 * length(SIZE.M), 3 * length(MEAN.READS.M))
par(mfrow = c(length(MEAN.READS.M), length(SIZE.M)), mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
for (mr in MEAN.READS.M) for (sz in SIZE.M) {
  d <- MDE.RATIO[MDE.RATIO$mean_reads == mr & MDE.RATIO$size == sz & MDE.RATIO$trans_mode == "trans_indep", ]
  plot(NA, xlim = range(FRAC.M), ylim = c(0.5, 3), xlab = "hybrid split fraction f (cis cells)", ylab = "MDE trans / MDE cis",
       main = sprintf("mean %s, size %s", mr, sz), cex.main = 0.8)
  abline(h = 1, col = COLOR.GREY[["mid"]], lty = 3)
  abline(v = MODES.REAL$f_star[is.finite(MODES.REAL$f_star)], col = COLOR.ACCENT, lty = 2)
  for (k in seq_along(PHI.M)) with(d[d$phi == PHI.M[k], ], lines(f, mde_ratio, type = "b", pch = 19, cex = 0.5, col = colorRampPalette(POWER.COLOR)(length(PHI.M))[k]))
}
plot.new(); legend("center", legend = paste("phi =", PHI.M), col = colorRampPalette(POWER.COLOR)(length(PHI.M)), pch = 19, lty = 1, bty = "n")
dev.off()

## 10.13 Type I error at ratio 1, by mode and null version, and the split fraction at which cis and trans cross
# The SIZE.RATIO == 1 entry of power_bh is the raw false-positive rate at p < ALPHA. The binomial SE with NJ.M
# datasets is sqrt(ALPHA * (1 - ALPHA) / NJ.M). The size axis carries the effect; the mean axis has no effect in
# this parametrization, so its rate is a second control.
MODES.TYPE1 <- aggregate(power_bh ~ mode + axis + phi, data = MODES.POWER[MODES.POWER$ratio == 1, ], FUN = function(x) c(mean = mean(x), max = max(x)))
MODES.TYPE1 <- do.call(data.frame, MODES.TYPE1)
names(MODES.TYPE1) <- c("mode", "axis", "phi", "type1_mean", "type1_max")
write.csv(MODES.TYPE1, file.path(TABLE.DIR, "power_modes_type1.csv"), row.names = FALSE)
cat(sprintf("type I error at SIZE.RATIO = 1 (nominal %.2f; binomial SE per row %.3f), mean and max over the grid:\n", ALPHA, sqrt(ALPHA * (1 - ALPHA) / NJ.M)))
print(MODES.TYPE1, digits = 3, row.names = FALSE)
CROSS.RATIO <- SIZE.RATIO.M[which.min(abs(SIZE.RATIO.M - 2))]
F.CROSS <- do.call(rbind, lapply(c("trans_indep", "trans_paired"), function(tm) modes_f_cross(MODES.POWER, CROSS.RATIO, trans_mode = tm, axis = "disp", col = MODES.PW.COL)))
write.csv(F.CROSS, file.path(TABLE.DIR, "power_modes_f_cross.csv"), row.names = FALSE)
cat(sprintf("split fraction at which cis and trans power are equal (SIZE.RATIO = %s, size axis); pipeline f_disp = %.3f:\n", CROSS.RATIO, MODES.REAL$f_star[["f_disp"]]))
print(F.CROSS[F.CROSS$trans_mode == "trans_indep", c("mean_reads", "size", "phi", "f_cross")], digits = 3, row.names = FALSE)

## 10.14 Observed null contrast: simulated SD against the variance statement in the project notes
# project_context.md states var_trans = var_cis + 4 / N_p. At ratio 1 the SD of the observed size contrast is
# simulated for cis (paired alleles) and trans (parental ratio minus hybrid ratio). The statement compares the two
# on the same hybrid cells, which holds at f = 0.5 only (cis and trans then use N.H / 2 cells each), so the
# excess variance times N.P, the constant the statement puts at 4, is summarized at f = 0.5. It is reported
# in log2 and natural-log units because the notes do not name the scale, and it varies with mean, size and phi.
SD.WIDE <- reshape(MODES.SD[MODES.SD$axis == "disp", c("i", "mean_reads", "size", "phi", "f", "N.P", "N.H", "depth", "mode", "sd")],
                   idvar = c("i", "mean_reads", "size", "phi", "f", "N.P", "N.H", "depth"), timevar = "mode", direction = "wide")
SD.WIDE$excess_log2   <- SD.WIDE$sd.trans^2 - SD.WIDE$sd.cis^2
SD.WIDE$const_log2    <- SD.WIDE$excess_log2 * SD.WIDE$N.P
SD.WIDE$const_ln      <- SD.WIDE$const_log2 * log(2)^2
write.csv(SD.WIDE, file.path(TABLE.DIR, "power_modes_null_sd.csv"), row.names = FALSE)
cat("simulated SD of the observed size contrast under the null (ratio 1), median by f and phi:\n")
print(round(tapply(SD.WIDE$sd.cis, list(phi = SD.WIDE$phi, f = SD.WIDE$f), median), 3))
print(round(tapply(SD.WIDE$sd.trans, list(phi = SD.WIDE$phi, f = SD.WIDE$f), median), 3))
SD.HALF <- SD.WIDE[SD.WIDE$f == 0.5, ]
if (nrow(SD.HALF) > 0) {
  cat("(var_trans - var_cis) * N.P at f = 0.5, the notes' constant of 4: log2 scale, then natural-log scale (median by phi):\n")
  print(round(tapply(SD.HALF$const_log2, SD.HALF$phi, median, na.rm = TRUE), 2))
  print(round(tapply(SD.HALF$const_ln, SD.HALF$phi, median, na.rm = TRUE), 2))
}

save(MODES.INPUTS, MODES.REAL, GEN.CHECK, MODES.POWER, MODES.SD, MODES.MDE, MDE.RATIO, MODES.TYPE1, F.CROSS, SD.WIDE,
     file = file.path(CHECKPOINT.DIR, "section10_modes_checkpoint.rda"))
