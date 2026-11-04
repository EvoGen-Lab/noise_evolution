###############################################################
### estimator_leak.R
### Two checks of the Cox-Reid NB dispersion estimator against the Appendix A moment estimator on the data, run locally from the
### Section 1 checkpoint alone (no Section 2 rerun and no cluster output):
###   L1  estimator comparison: log2(phi Cox-Reid / phi moment) per gene and group, summarized in bins of the mean count per cell and of phi;
###   L2  leak test: across genes, Delta log2 phi (Cox-Reid) - Delta log2 phi (moment) regressed on Delta log2 mu for the cis and trans
###       contrasts. A slope near zero means the shape assumption of the NB fit does not carry mean differences into the dispersion.
### The offsets and share corrections are those of analysis.R Section 2: a parent's offset is its library over the median parental library,
### each hybrid allele's is half the cell's hybrid depth, the allele dispersions carry the share corrections R.C and R.E, and the allele means
### are divided by their abundance factors 2 c-bar and 2 (1 - c-bar). The contrasts use all hybrid cells: cis compares the two alleles and
### trans is the parental log2 ratio minus the allele log2 ratio.
### Mean and dispersion fitted on the same cells have correlated errors unless the counts follow the NB (docs/validation_checks.md, check A10),
### and the two estimators carry that covariance to different degrees, which alone can give the leak regression a slope. The test is therefore
### run in two arms: "same" (all cells for both) and "disjoint" (each group's cells split at random into halves, the mean from one half and both
### dispersions from the other). Only a slope in the disjoint arm reads as leakage of the mean into the dispersion.
###   L3  the L1 gap and the L2 slopes with the top 0.5% of cells by scaled count left out of every fit (influence of a few cells);
###   L4  each group's expression-factor third moment against the gamma (check A11 on the Section 1 data);
###   L5  the leak test on raw dispersion contrasts (Delta phi);
###   L6  a variance term linear in the mean that the NB lacks, Var(Y) = (1 + a) M + phi M^2, per gene and group, which makes the NB phi fall with the mean;
###   L7  the cells that carry the extreme counts: per cell the number of genes for which it is among the top 0.5% (against the same tally on NB data
###       simulated from the fits), a description of the cells far above that null, and L1, L2 and L6 refitted without them;
###   L8  how each gene's counts scale with cell depth (slope beta of log E[y] on log O; the library-size offset assumes beta = 1), and L6 and the
###       Delta phi / Delta mu slopes refitted with each gene's own scaling O^beta.
### Results go to results/checks/L1.rds to L8.rds with rows in summary.csv; the figure goes to figures/extra/check_L2_leak.pdf.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT    <- "analysis/checks/estimator_leak.R"
CORES           <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.LEAK       <- 5026
LEAK.MIN.MEAN   <- 0.02     # mean count per cell required in every group
LEAK.MIN.CELLS  <- 50       # cells with a count required in every group
LEAK.MAX.GENES  <- Inf      # a random subset of the eligible genes when finite (for a quick run)
LEAK.BOOT       <- 1000     # gene-resampling draws for the slope interval
MBAR.BREAKS     <- c(0, 0.1, 0.3, 1, 3, 10, Inf)
PHI.BREAKS      <- c(0, 0.05, 0.2, 0.5, 1, Inf)
set.seed(SEED.SECTION + 23)

stopifnot("estimator_leak.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE (cells after QC, genes after the union filter), C.BAR, R.C, R.E

## Offsets and share corrections of Section 2. The hybrid cells keep their Section 1 order; nothing here depends on the depth order Section 2 uses.
DEPTH.REF.L <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
LEAK.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
LEAK.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF.L, Pe = colSums(MIX.SE) / DEPTH.REF.L,
                   Hc = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF.L / 2, He = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF.L / 2)
LEAK.ETA   <- c(Pc = 0, Pe = 0, Hc = R.C, He = R.E)

## Genes with data in all four groups: a mean count of LEAK.MIN.MEAN per cell and LEAK.MIN.CELLS cells with a count in every group.
LEAK.GENES <- rownames(MIX.SC)[Reduce(`&`, lapply(LEAK.MATS, function(m) rowMeans(m) >= LEAK.MIN.MEAN & rowSums(m > 0) >= LEAK.MIN.CELLS))]
if (is.finite(LEAK.MAX.GENES) && length(LEAK.GENES) > LEAK.MAX.GENES) LEAK.GENES <- with_local_seed(SEED.LEAK, sample(LEAK.GENES, LEAK.MAX.GENES))
cat(sprintf("Estimator comparison: %d genes; cells Pc %d, Pe %d, hybrid %d; c-bar %.4f, eta_Hc %.3g, eta_He %.3g\n", length(LEAK.GENES),
            ncol(MIX.SC), ncol(MIX.SE), ncol(HYB.SC), C.BAR, R.C, R.E))

## Random halves of each group's cells for the disjoint arm, drawn once; the two hybrid alleles share their halves because they share cells.
LEAK.HALVES <- with_local_seed(SEED.LEAK + 1, local({
  h <- sample.int(ncol(HYB.SC), floor(ncol(HYB.SC) / 2))
  list(Pc = sample.int(ncol(MIX.SC), floor(ncol(MIX.SC) / 2)), Pe = sample.int(ncol(MIX.SE), floor(ncol(MIX.SE) / 2)), Hc = h, He = h)
}))

LEAK.FITS <- do.call(rbind, mclapply(LEAK.GENES, estimator_pair_gene, mats = LEAK.MATS, expos = LEAK.EXPOS, eta = LEAK.ETA, halves = LEAK.HALVES, mc.cores = CORES))
rownames(LEAK.FITS) <- LEAK.GENES

##############################################################################
## L1. log2(phi Cox-Reid / phi moment) by bins of the mean count per cell and of phi (Cox-Reid), per group, from all cells
# A gene enters a bin when both estimates are positive; the shares with either estimate at zero are reported beside the ratio.
L1 <- list()
for (k in names(LEAK.MATS)) {
  cr <- LEAK.FITS[, paste0("same.", k, ".phi_cr")]; mo <- LEAK.FITS[, paste0("same.", k, ".phi_mom")]; mb <- LEAK.FITS[, paste0("same.", k, ".mbar")]
  bm <- cut(mb, MBAR.BREAKS, right = FALSE); bp <- cut(cr, PHI.BREAKS, right = FALSE)
  for (lm_ in levels(bm)) for (lp in levels(bp)) {
    ix <- which(bm == lm_ & bp == lp); pos <- ix[cr[ix] > 0 & mo[ix] > 0]
    if (!length(ix)) next
    r <- log2(cr[pos] / mo[pos])
    L1[[length(L1) + 1]] <- data.frame(group = k, mean_count = lm_, phi_cr = lp, genes = length(ix), both_positive = length(pos),
                                       moment_zero = sum(mo[ix] == 0), median_log2_ratio = if (length(r)) median(r) else NA_real_,
                                       q25 = if (length(r)) unname(quantile(r, 0.25)) else NA_real_, q75 = if (length(r)) unname(quantile(r, 0.75)) else NA_real_,
                                       mean_log2_ratio = if (length(r)) mean(r) else NA_real_, se_mean = if (length(r) > 1) sd(r) / sqrt(length(r)) else NA_real_)
  }
}
L1 <- do.call(rbind, L1)
L1.PHI0 <- data.frame(group = names(LEAK.MATS),
                      cr_zero = vapply(names(LEAK.MATS), function(k) mean(LEAK.FITS[, paste0("same.", k, ".phi_cr")] == 0), numeric(1)),
                      moment_zero = vapply(names(LEAK.MATS), function(k) mean(LEAK.FITS[, paste0("same.", k, ".phi_mom")] == 0), numeric(1)), row.names = NULL)
print(L1[L1$both_positive >= 20, c("group", "mean_count", "phi_cr", "both_positive", "median_log2_ratio", "q25", "q75")], digits = 3, row.names = FALSE)
x <- L1[L1$both_positive >= 20, ]
check_save("L1", list(bins = L1, zero_shares = L1.PHI0, fits = LEAK.FITS), "run",
           headline = c(sprintf("%s mean count %s, phi %s: median log2(CR / moment) %.3f (IQR %.3f to %.3f), %d genes", x$group, x$mean_count, x$phi_cr,
                                x$median_log2_ratio, x$q25, x$q75, x$both_positive),
                        sprintf("%s: phi = 0 in %.3f (Cox-Reid) and %.3f (moment) of genes", L1.PHI0$group, L1.PHI0$cr_zero, L1.PHI0$moment_zero)),
           mc_error = "bins with at least 20 genes in the headline; se_mean per bin in the rds", script = CHECK.SCRIPT)

##############################################################################
## L2. Leak test: [Delta log2 phi (Cox-Reid) - Delta log2 phi (moment)] on Delta log2 mu across genes, for cis and trans, in both arms
# cis: Hc against He, means divided by their abundance factors. trans: the parental log2 ratio minus the allele log2 ratio. A gene enters a contrast
# when all of its groups have both dispersion estimates positive. The fit is repeated within terciles of the genes' mean count per cell.
L2.DATA <- rbind(leak_contrasts(LEAK.FITS, "same", C.BAR), leak_contrasts(LEAK.FITS, "disjoint", C.BAR))
L2 <- list()
for (arm in c("same", "disjoint")) for (md in c("cis", "trans")) {
  d <- L2.DATA[L2.DATA$arm == arm & L2.DATA$contrast == md, ]
  L2[[paste(arm, md, "all")]] <- data.frame(arm = arm, contrast = md, genes_in = "all", leak_regression(d$dlog2mu, d$diff, d$dphi_cr, d$dphi_mom, LEAK.BOOT, SEED.LEAK + 2))
  for (t in levels(d$tercile)) { i <- d$tercile == t
    L2[[paste(arm, md, t)]] <- data.frame(arm = arm, contrast = md, genes_in = paste("mean count tercile", t), leak_regression(d$dlog2mu[i], d$diff[i], d$dphi_cr[i], d$dphi_mom[i], LEAK.BOOT, SEED.LEAK + 3)) }
}
L2 <- do.call(rbind, L2); rownames(L2) <- NULL
print(L2, digits = 3, row.names = FALSE)
x <- L2[L2$genes_in == "all", ]
check_save("L2", list(slopes = L2, data = L2.DATA), "run",
           headline = c(sprintf("%s %s: slope %.4f (SE %.4f; gene bootstrap 95%% %.4f to %.4f), Spearman %.3f, %d genes; own slopes Cox-Reid %.4f, moment %.4f",
                                x$arm, x$contrast, x$slope, x$slope_se, x$boot_lo, x$boot_hi, x$spearman, x$genes, x$slope_cr, x$slope_mom)),
           mc_error = sprintf("gene bootstrap, %d draws", LEAK.BOOT), script = CHECK.SCRIPT)

fig_pdf("extra/check_L2_leak.pdf", 9, 9)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (arm in c("same", "disjoint")) for (md in c("cis", "trans")) {
  d <- L2.DATA[L2.DATA$arm == arm & L2.DATA$contrast == md & is.finite(L2.DATA$dlog2mu) & is.finite(L2.DATA$diff), ]
  s <- L2[L2$arm == arm & L2$contrast == md & L2$genes_in == "all", ]
  plot(d$dlog2mu, d$diff, pch = 19, cex = 0.3, col = adjustcolor(COLOR.GREY[["dark"]], 0.4), xlab = "Delta log2 mu", ylab = "Delta log2 phi: Cox-Reid - moment",
       main = sprintf("%s, %s cells: slope %.3f [%.3f, %.3f]", md, arm, s$slope, s$boot_lo, s$boot_hi))
  abline(h = 0, lty = 2, col = COLOR.ACCENT); abline(lm(diff ~ dlog2mu, data = d), col = COLOR.ACCENT, lwd = 2)
}
dev.off()

##############################################################################
## L3. Influence of the cells with the largest counts: every fit repeated without the top LEAK.TRIM fraction of cells by scaled count y / O
# The moment estimator weights the squared deviation of every cell, so a few cells with very large counts move it more than the NB likelihood. If the
# Cox-Reid / moment gap of L1 and the leak slopes of L2 shrink once those cells are left out, a few cells carry them. Trimming lowers both estimates, so
# the comparison is between the estimators on the same trimmed cells, next to the untrimmed values.
LEAK.TRIM <- 0.005
LEAK.FITS.TRIM <- do.call(rbind, mclapply(LEAK.GENES, estimator_pair_gene, mats = LEAK.MATS, expos = LEAK.EXPOS, eta = LEAK.ETA, halves = LEAK.HALVES, trim = LEAK.TRIM, mc.cores = CORES))
rownames(LEAK.FITS.TRIM) <- LEAK.GENES
L3.GAP <- list()
for (k in names(LEAK.MATS)) {
  bp <- cut(LEAK.FITS[, paste0("same.", k, ".phi_cr")], PHI.BREAKS, right = FALSE)
  for (fits in c("all cells", "trimmed")) {
    F <- if (fits == "all cells") LEAK.FITS else LEAK.FITS.TRIM
    cr <- F[, paste0("same.", k, ".phi_cr")]; mo <- F[, paste0("same.", k, ".phi_mom")]; ok <- cr > 0 & mo > 0
    L3.GAP[[paste(k, fits)]] <- data.frame(group = k, fits = fits, phi_cr_bin = levels(bp), genes = as.vector(table(bp[ok])),
                                           median_log2_ratio = vapply(levels(bp), function(l) median(log2(cr / mo)[ok & bp == l]), numeric(1)),
                                           median_log2_trim_over_all_cr = vapply(levels(bp), function(l) median(log2(LEAK.FITS.TRIM[, paste0("same.", k, ".phi_cr")] / LEAK.FITS[, paste0("same.", k, ".phi_cr")])[ok & bp == l]), numeric(1)),
                                           median_log2_trim_over_all_mom = vapply(levels(bp), function(l) median(log2(LEAK.FITS.TRIM[, paste0("same.", k, ".phi_mom")] / LEAK.FITS[, paste0("same.", k, ".phi_mom")])[ok & bp == l]), numeric(1)),
                                           row.names = NULL)
  }
}
L3.GAP <- do.call(rbind, L3.GAP)
L3.DATA <- leak_contrasts(LEAK.FITS.TRIM, "disjoint", C.BAR)
L3.SLOPES <- list()
for (md in c("cis", "trans")) {
  d <- L3.DATA[L3.DATA$contrast == md, ]
  L3.SLOPES[[md]] <- data.frame(arm = "disjoint, trimmed", contrast = md, leak_regression(d$dlog2mu, d$diff, d$dphi_cr, d$dphi_mom, LEAK.BOOT, SEED.LEAK + 4))
}
L3.SLOPES <- do.call(rbind, L3.SLOPES); rownames(L3.SLOPES) <- NULL
print(L3.GAP[L3.GAP$genes >= 20, ], digits = 3, row.names = FALSE); print(L3.SLOPES, digits = 3, row.names = FALSE)
x <- L3.GAP[L3.GAP$genes >= 20, ]
check_save("L3", list(gap = L3.GAP, slopes = L3.SLOPES, fits = LEAK.FITS.TRIM, trim = LEAK.TRIM), "run",
           headline = c(sprintf("%s phi %s, %s: median log2(CR / moment) %.3f", x$group, x$phi_cr_bin, x$fits, x$median_log2_ratio),
                        sprintf("trimmed %s: leak slope %.4f (95%% %.4f to %.4f); own slopes Cox-Reid %.3f, moment %.3f", L3.SLOPES$contrast, L3.SLOPES$slope, L3.SLOPES$boot_lo,
                                L3.SLOPES$boot_hi, L3.SLOPES$slope_cr, L3.SLOPES$slope_mom)),
           mc_error = sprintf("top %.1f%% of cells by scaled count left out of every fit; gene bootstrap %d draws", 100 * LEAK.TRIM, LEAK.BOOT), script = CHECK.SCRIPT)

##############################################################################
## L4. Third moment of the expression factor against the gamma, per group (check A11 on the Section 1 data)
# The NB fit is exact for a gamma factor; the excess E[(a - 1)^3] - 2 phi^2 measures the departure. Each tested gene is compared with NB data simulated at its
# own mean and variance on the same offsets (factor_excess_test()); the group statistic is the mean standardized deviation z over the tested genes (0 for a
# gamma factor, positive for a factor more right-skewed than the gamma). Genes are a random sample of those with phi > 0 in the group.
L4.GENES.PER.GROUP <- 400; L4.SIM <- 200
L4 <- list()
for (k in names(LEAK.MATS)) {
  gk <- LEAK.GENES[LEAK.FITS[, paste0("same.", k, ".phi_cr")] > 0]
  gk <- with_local_seed(SEED.LEAK + 5, sample(gk, min(L4.GENES.PER.GROUP, length(gk))))
  tst <- with_local_seed(SEED.LEAK + 6, t(vapply(gk, factor_excess_test, numeric(6), mat = LEAK.MATS[[k]], expo = LEAK.EXPOS[[k]], B = L4.SIM)))
  mb <- cut(LEAK.FITS[gk, paste0("same.", k, ".mbar")], MBAR.BREAKS, right = FALSE); ok <- is.finite(tst[, "z"])
  for (l in c("all", levels(mb))) { i <- ok & (l == "all" | mb == l)
    if (sum(i) < 10) next
    L4[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), mean_z = mean(tst[i, "z"]), mean_z_se = 1 / sqrt(sum(i)), median_ratio = median(tst[i, "ratio"]),
                                    median_null_ratio = median(tst[i, "null_ratio"]), share_p05 = mean(tst[i, "p"] < 0.05), row.names = NULL) }
}
L4 <- do.call(rbind, L4)
print(L4, digits = 3, row.names = FALSE)
check_save("L4", L4, "run", headline = sprintf("%s mean count %s: mean z %.2f (SE %.2f), k3 / (2 phi^2) median %.2f against %.2f under a gamma factor, p < 0.05 in %.3f of %d genes",
                                               L4$group, L4$mean_count, L4$mean_z, L4$mean_z_se, L4$median_ratio, L4$median_null_ratio, L4$share_p05, L4$genes),
           mc_error = sprintf("%d null draws per gene; 0.05 expected at p < 0.05 under a gamma factor", L4.SIM), script = CHECK.SCRIPT)

##############################################################################
## L5. Leak test on the raw scale: the difference of the two estimators' raw dispersion contrasts (Delta phi) on Delta log2 mu, both arms
# A slope on the log2 scale that vanishes on the raw scale comes from small phi, where a log2 contrast magnifies small differences.
L5.DATA <- rbind(leak_contrasts(LEAK.FITS, "same", C.BAR, "raw"), leak_contrasts(LEAK.FITS, "disjoint", C.BAR, "raw"))
L5 <- list()
for (arm in c("same", "disjoint")) for (md in c("cis", "trans")) {
  d <- L5.DATA[L5.DATA$arm == arm & L5.DATA$contrast == md, ]
  L5[[paste(arm, md)]] <- data.frame(arm = arm, contrast = md, leak_regression(d$dlog2mu, d$diff, d$dphi_cr, d$dphi_mom, LEAK.BOOT, SEED.LEAK + 7))
}
L5 <- do.call(rbind, L5); rownames(L5) <- NULL
print(L5, digits = 3, row.names = FALSE)
check_save("L5", list(slopes = L5, data = L5.DATA), "run",
           headline = sprintf("%s %s, raw Delta phi: slope %.4f (95%% %.4f to %.4f), %d genes; own slopes Cox-Reid %.4f, moment %.4f", L5$arm, L5$contrast, L5$slope, L5$boot_lo, L5$boot_hi,
                              L5$genes, L5$slope_cr, L5$slope_mom),
           mc_error = sprintf("gene bootstrap, %d draws", LEAK.BOOT), script = CHECK.SCRIPT)

##############################################################################
## L6. A variance term linear in the mean that the NB lacks: Var(Y) = (1 + a) M + phi M^2 (linear_excess_gene())
# An NB gene has a = 0. A technical source of Poisson-like extra variance (for example counts duplicated in a fraction of molecules) gives a > 0, and the NB
# fit then returns about phi + a / M, a dispersion that falls with the mean count. That alone makes Delta log2 phi fall with Delta log2 mu across genes
# with no change in the expression factor, and it moves the moment estimator, whose weights favor low-count cells, more than the NB fit. a and phi are
# separated within each gene by the spread of the offsets. Reported per group and mean-count bin: the median and mean of a with the share of genes beyond
# two SEs on each side (0.025 each expected for NB genes), and the trend of log2 phi (Cox-Reid) on log2 mu across genes, with and without the a / M term.
L6 <- list(); L6.GENE <- list()
for (k in names(LEAK.MATS)) {
  r <- t(vapply(LEAK.GENES, linear_excess_gene, numeric(4), mat = LEAK.MATS[[k]], expo = LEAK.EXPOS[[k]]))
  z <- r[, "a"] / r[, "a_se"]; mb <- cut(r[, "mbar"], MBAR.BREAKS, right = FALSE)
  L6.GENE[[k]] <- data.frame(group = k, gene = LEAK.GENES, r, z = z, phi_cr = LEAK.FITS[, paste0("same.", k, ".phi_cr")], row.names = NULL)
  ph <- LEAK.FITS[, paste0("same.", k, ".phi_cr")]; mu <- LEAK.FITS[, paste0("same.", k, ".mu")]; ok <- ph > 0 & r[, "phi_lin"] > 0 & is.finite(r[, "phi_lin"])
  trend_cr <- unname(coef(lm(log2(ph[ok]) ~ log2(mu[ok])))[2]); trend_lin <- unname(coef(lm(log2(r[ok, "phi_lin"]) ~ log2(mu[ok])))[2])
  for (l in c("all", levels(mb))) { i <- is.finite(z) & (l == "all" | mb == l)
    if (sum(i) < 10) next
    L6[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), median_a = median(r[i, "a"]), mean_a = mean(r[i, "a"]), mean_a_se = sd(r[i, "a"]) / sqrt(sum(i)),
                                    share_above = mean(z[i] > 2), share_below = mean(z[i] < -2), trend_log2phi_cr_on_log2mu = trend_cr, trend_log2phi_lin_on_log2mu = trend_lin, row.names = NULL) }
}
L6 <- do.call(rbind, L6); L6.GENE <- do.call(rbind, L6.GENE)
print(L6, digits = 3, row.names = FALSE)
x <- L6[L6$mean_count == "all", ]
check_save("L6", list(summary = L6, genes = L6.GENE), "run",
           headline = c(sprintf("%s: a median %.3f, mean %.3f (SE %.3f); a > 2 SE in %.3f and < -2 SE in %.3f of %d genes; trend of log2 phi on log2 mu across genes: Cox-Reid %.3f, with the linear term removed %.3f",
                                x$group, x$median_a, x$mean_a, x$mean_a_se, x$share_above, x$share_below, x$genes, x$trend_log2phi_cr_on_log2mu, x$trend_log2phi_lin_on_log2mu)),
           mc_error = "sandwich SE per gene; SE of the mean a over genes", script = CHECK.SCRIPT)

##############################################################################
## L7. Cells that carry the extreme counts, and the fits without them
# For every cell, the number of genes for which it is among the top LEAK.TRIM fraction of cells, ranked by the upper-tail probability of its count under the
# gene's NB fit (top_cell_tally()), which puts cells of every depth on one scale. The tally is read against the same tally on counts simulated from each
# gene's NB fit (mean and Cox-Reid phi) on the same cells. Cells above the largest simulated tally of their group (one simulation) are described by library size, detected genes, allele share (hybrid)
# and the share of their library in the group's 50 most expressed genes, and removed; the fits of L1, L2 and L6 are then repeated on the remaining cells.
L7.TALLY <- list(); L7.SUMMARY <- list(); L7.OUT <- list()
TOP50 <- lapply(LEAK.MATS, function(m) order(rowSums(m), decreasing = TRUE)[1:50])
for (k in names(LEAK.MATS)) {
  m <- LEAK.MATS[[k]]; e <- LEAK.EXPOS[[k]]; mu <- LEAK.FITS[, paste0("same.", k, ".mu")]; ph <- LEAK.FITS[, paste0("same.", k, ".phi_cr")]
  sim <- with_local_seed(SEED.LEAK + 8, t(vapply(seq_along(LEAK.GENES), function(j) rnbinom(length(e), mu = mu[j] * e, size = if (ph[j] > 0) 1 / ph[j] else 1e8), numeric(length(e)))))
  rownames(sim) <- LEAK.GENES
  names(mu) <- names(ph) <- LEAK.GENES
  obs <- top_cell_tally(m, e, LEAK.GENES, LEAK.TRIM, mu, ph); nul <- top_cell_tally(sim, e, LEAK.GENES, LEAK.TRIM, mu, ph)
  lib <- colSums(m); det <- colSums(m > 0); top50 <- colSums(m[TOP50[[k]], , drop = FALSE]) / lib
  out <- obs > max(nul); L7.OUT[[k]] <- which(out)
  L7.TALLY[[k]] <- data.frame(group = k, cell = if (is.null(colnames(m))) seq_len(ncol(m)) else colnames(m), tally = obs, tally_null_sorted = sort(nul), library = lib, detected = det, top50_share = top50, outlier = out, row.names = NULL)
  k1 <- max(1, ceiling(0.01 * length(e)))
  L7.SUMMARY[[k]] <- data.frame(group = k, cells = length(e), genes = length(LEAK.GENES), expected_tally = length(LEAK.GENES) * max(1, ceiling(LEAK.TRIM * length(e))) / length(e),
                                max_null = max(nul), outlier_cells = sum(out), share_trims_top1pct_cells = sum(sort(obs, TRUE)[1:k1]) / sum(obs),
                                share_trims_top1pct_cells_null = sum(sort(nul, TRUE)[1:k1]) / sum(nul),
                                spearman_tally_library = cor(obs, lib, method = "spearman"), spearman_tally_library_null = cor(nul, lib, method = "spearman"),
                                outlier_library_quantile = if (any(out)) median(rank(lib)[out] / length(lib)) else NA_real_,
                                outlier_detected_quantile = if (any(out)) median(rank(det)[out] / length(det)) else NA_real_,
                                outlier_top50_share = if (any(out)) median(top50[out]) else NA_real_, other_top50_share = median(top50[!out]), outlier_chat_z = NA_real_, row.names = NULL)
}
## The hybrid alleles share their cells, so a cell flagged in either allele leaves both.
L7.OUT$Hc <- L7.OUT$He <- sort(union(L7.OUT$Hc, L7.OUT$He))
## The allele share of the flagged hybrid cells, as a z-score against all hybrid cells (C.HAT is in the cell order of the Section 1 matrices).
if (length(C.HAT) == ncol(HYB.SC) && length(L7.OUT$Hc)) L7.SUMMARY$Hc$outlier_chat_z <- L7.SUMMARY$He$outlier_chat_z <- median((C.HAT[L7.OUT$Hc] - mean(C.HAT)) / sd(C.HAT))
L7.SUMMARY <- do.call(rbind, L7.SUMMARY)
print(L7.SUMMARY, digits = 3, row.names = FALSE)
## Refit without the outlier cells
L7.MATS <- lapply(setNames(names(LEAK.MATS), names(LEAK.MATS)), function(k) if (length(L7.OUT[[k]])) LEAK.MATS[[k]][, -L7.OUT[[k]], drop = FALSE] else LEAK.MATS[[k]])
L7.EXPOS <- lapply(setNames(names(LEAK.EXPOS), names(LEAK.EXPOS)), function(k) if (length(L7.OUT[[k]])) LEAK.EXPOS[[k]][-L7.OUT[[k]]] else LEAK.EXPOS[[k]])
L7.HALVES <- with_local_seed(SEED.LEAK + 9, local({
  h <- sample.int(ncol(L7.MATS$Hc), floor(ncol(L7.MATS$Hc) / 2))
  list(Pc = sample.int(ncol(L7.MATS$Pc), floor(ncol(L7.MATS$Pc) / 2)), Pe = sample.int(ncol(L7.MATS$Pe), floor(ncol(L7.MATS$Pe) / 2)), Hc = h, He = h)
}))
L7.FITS <- do.call(rbind, mclapply(LEAK.GENES, estimator_pair_gene, mats = L7.MATS, expos = L7.EXPOS, eta = LEAK.ETA, halves = L7.HALVES, mc.cores = CORES))
rownames(L7.FITS) <- LEAK.GENES
L7.GAP <- list(); L7.SLOPES <- list(); L7.A <- list()
for (k in names(LEAK.MATS)) {
  bp <- cut(LEAK.FITS[, paste0("same.", k, ".phi_cr")], PHI.BREAKS, right = FALSE)
  for (fits in c("all cells", "outlier cells removed")) {
    F <- if (fits == "all cells") LEAK.FITS else L7.FITS
    cr <- F[, paste0("same.", k, ".phi_cr")]; mo <- F[, paste0("same.", k, ".phi_mom")]; ok <- cr > 0 & mo > 0
    L7.GAP[[paste(k, fits)]] <- data.frame(group = k, fits = fits, phi_cr_bin = levels(bp), genes = as.vector(table(bp[ok])),
                                           median_log2_ratio = vapply(levels(bp), function(l) median(log2(cr / mo)[ok & bp == l]), numeric(1)),
                                           median_log2_phi_cr_change = vapply(levels(bp), function(l) median(log2(L7.FITS[, paste0("same.", k, ".phi_cr")] / LEAK.FITS[, paste0("same.", k, ".phi_cr")])[ok & bp == l]), numeric(1)),
                                           row.names = NULL)
  }
  r <- t(vapply(LEAK.GENES, linear_excess_gene, numeric(4), mat = L7.MATS[[k]], expo = L7.EXPOS[[k]])); z <- r[, "a"] / r[, "a_se"]
  L7.A[[k]] <- data.frame(group = k, removed_cells = length(L7.OUT[[k]]), median_a = median(r[, "a"], na.rm = TRUE), share_above = mean(z > 2, na.rm = TRUE), share_below = mean(z < -2, na.rm = TRUE), row.names = NULL)
}
L7.DATA <- leak_contrasts(L7.FITS, "disjoint", C.BAR)
for (md in c("cis", "trans")) { d <- L7.DATA[L7.DATA$contrast == md, ]
  L7.SLOPES[[md]] <- data.frame(arm = "disjoint, outlier cells removed", contrast = md, leak_regression(d$dlog2mu, d$diff, d$dphi_cr, d$dphi_mom, LEAK.BOOT, SEED.LEAK + 10)) }
L7.GAP <- do.call(rbind, L7.GAP); L7.SLOPES <- do.call(rbind, L7.SLOPES); L7.A <- do.call(rbind, L7.A); rownames(L7.SLOPES) <- NULL
print(L7.GAP[L7.GAP$genes >= 20, ], digits = 3, row.names = FALSE); print(L7.SLOPES, digits = 3, row.names = FALSE); print(L7.A, digits = 3, row.names = FALSE)
check_save("L7", list(summary = L7.SUMMARY, cells = do.call(rbind, L7.TALLY), removed = L7.OUT, gap = L7.GAP, slopes = L7.SLOPES, linear_excess = L7.A, fits = L7.FITS), "run",
           headline = c(sprintf("%s: %d cells above the largest simulated tally (%d; expected tally %.1f); top 1%% of cells hold %.3f of the trims (NB %.3f); Spearman(tally, library) %.2f (NB %.2f)",
                                L7.SUMMARY$group, L7.SUMMARY$outlier_cells, L7.SUMMARY$max_null, L7.SUMMARY$expected_tally, L7.SUMMARY$share_trims_top1pct_cells,
                                L7.SUMMARY$share_trims_top1pct_cells_null, L7.SUMMARY$spearman_tally_library, L7.SUMMARY$spearman_tally_library_null),
                        sprintf("without them, %s: leak slope %.4f (95%% %.4f to %.4f); own slopes Cox-Reid %.3f, moment %.3f", L7.SLOPES$contrast, L7.SLOPES$slope, L7.SLOPES$boot_lo,
                                L7.SLOPES$boot_hi, L7.SLOPES$slope_cr, L7.SLOPES$slope_mom),
                        sprintf("without them, %s: a median %.3f, a > 2 SE in %.3f of genes", L7.A$group, L7.A$median_a, L7.A$share_above)),
           mc_error = "one NB simulation per group for the tally null; gene bootstrap for the slopes", script = CHECK.SCRIPT)

##############################################################################
## L8. Depth scaling of each gene, and L6 and the Delta phi / Delta mu slopes with each gene's own scaling
# The offset assumes every gene keeps a fixed share of the library (beta = 1). A gene whose share changes with depth (beta != 1) leaves a misfit that grows
# with its mean, which L6 would read as a linear term rising with the mean. beta is reported by group and mean-count bin; L6 is then repeated with offsets
# O^beta, and the own slopes of Delta log2 phi (Cox-Reid) on Delta log2 mu are refitted with those offsets (all cells).
L8 <- list(); L8.BETA <- list(); L8.A <- list(); L8.FIT <- list()
for (k in names(LEAK.MATS)) {
  b <- t(vapply(LEAK.GENES, depth_slope_gene, numeric(3), mat = LEAK.MATS[[k]], expo = LEAK.EXPOS[[k]])); L8.BETA[[k]] <- b
  zb <- (b[, "beta"] - 1) / b[, "beta_se"]; mb <- cut(LEAK.FITS[, paste0("same.", k, ".mbar")], MBAR.BREAKS, right = FALSE)
  for (l in c("all", levels(mb))) { i <- is.finite(zb) & (l == "all" | mb == l)
    if (sum(i) < 10) next
    L8[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), median_beta = median(b[i, "beta"]), q25 = unname(quantile(b[i, "beta"], 0.25)), q75 = unname(quantile(b[i, "beta"], 0.75)),
                                    share_below_1 = mean(zb[i] < -2), share_above_1 = mean(zb[i] > 2), row.names = NULL) }
  beta <- setNames(ifelse(is.finite(b[, "beta"]), b[, "beta"], 1), LEAK.GENES)
  r <- t(vapply(LEAK.GENES, linear_excess_gene, numeric(4), mat = LEAK.MATS[[k]], expo = LEAK.EXPOS[[k]], beta = beta)); z <- r[, "a"] / r[, "a_se"]
  L8.FIT[[k]] <- t(vapply(LEAK.GENES, function(g) fit_dispersion(as.numeric(LEAK.MATS[[k]][g, ]), LEAK.EXPOS[[k]]^beta[[g]], LEAK.ETA[[k]])[c("mu", "disp")], numeric(2)))
  ok <- L8.FIT[[k]][, "disp"] > 0
  L8.A[[k]] <- data.frame(group = k, median_a = median(r[, "a"], na.rm = TRUE), share_above = mean(z > 2, na.rm = TRUE), share_below = mean(z < -2, na.rm = TRUE),
                          trend_log2phi_on_log2mu = unname(coef(lm(log2(L8.FIT[[k]][ok, "disp"]) ~ log2(L8.FIT[[k]][ok, "mu"])))[2]), row.names = NULL)
}
L8 <- do.call(rbind, L8); L8.A <- do.call(rbind, L8.A)
L8.OWN <- local({
  lmu <- list(Pc = log2(L8.FIT$Pc[, "mu"]), Pe = log2(L8.FIT$Pe[, "mu"]), Hc = log2(L8.FIT$Hc[, "mu"] / (2 * C.BAR)), He = log2(L8.FIT$He[, "mu"] / (2 * (1 - C.BAR))))
  lph <- lapply(L8.FIT, function(f) ifelse(f[, "disp"] > 0, log2(f[, "disp"]), NA_real_))
  dx_cis <- lmu$Hc - lmu$He; dy_cis <- lph$Hc - lph$He; dx_tr <- (lmu$Pc - lmu$Pe) - dx_cis; dy_tr <- (lph$Pc - lph$Pe) - dy_cis     # lm() drops genes with an NA (phi = 0)
  data.frame(contrast = c("cis", "trans"), own_slope_cox_reid_scaled = c(unname(coef(lm(dy_cis ~ dx_cis))[2]), unname(coef(lm(dy_tr ~ dx_tr))[2])),
             own_slope_cox_reid_library_offset = L2$slope_cr[L2$arm == "same" & L2$genes_in == "all"])
})
print(L8, digits = 3, row.names = FALSE); print(L8.A, digits = 3, row.names = FALSE); print(L8.OWN, digits = 3, row.names = FALSE)
x <- L8[L8$mean_count == "all", ]
check_save("L8", list(beta = L8, beta_genes = L8.BETA, linear_excess_scaled = L8.A, own_slopes = L8.OWN, fits_scaled = L8.FIT), "run",
           headline = c(sprintf("%s: beta median %.3f (IQR %.3f to %.3f); beta < 1 by 2 SE in %.3f, > 1 in %.3f of %d genes", x$group, x$median_beta, x$q25, x$q75, x$share_below_1, x$share_above_1, x$genes),
                        sprintf("with O^beta, %s: a median %.3f, a > 2 SE in %.3f; trend of log2 phi on log2 mu %.3f", L8.A$group, L8.A$median_a, L8.A$share_above, L8.A$trend_log2phi_on_log2mu),
                        sprintf("%s own slope of Delta log2 phi (Cox-Reid) on Delta log2 mu: %.3f with O^beta, %.3f with the library offset", L8.OWN$contrast, L8.OWN$own_slope_cox_reid_scaled,
                                L8.OWN$own_slope_cox_reid_library_offset)),
           mc_error = "quasi-Poisson SE of beta per gene", script = CHECK.SCRIPT)
