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
### power_output_<k>of<K>.rda files.
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
rm(SEC2)

console_start(10)
##############################################################################
## 10. POWER ANALYSIS                                                       ##
##############################################################################
# Run through 10.1 to build power_inputs.rda, then submit
# power.sub (a SLURM job array) to the cluster. Once every power_output_<k>of<K>.rda
# exists, run 10.2 onward to reshape them back into POWER and reproduce every plot.
# PR (permutation results with BH q-values) and N.PERM come from the Section 2
# checkpoint: PR sets PI1, the fraction of truly different genes in the simulated
# mixture, and N.PERM matches the permutation count to the gene-level analysis.

set.seed(1)

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
