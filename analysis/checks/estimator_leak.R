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
### Results go to results/checks/L1.rds and L2.rds with rows in summary.csv; the figure goes to figures/extra/check_L2_leak.pdf.
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
L2.DATA <- list(); L2 <- list()
for (arm in c("same", "disjoint")) {
  F <- LEAK.FITS[, startsWith(colnames(LEAK.FITS), paste0(arm, ".")), drop = FALSE]; colnames(F) <- sub(paste0("^", arm, "\\."), "", colnames(F))
  lmu <- list(Pc = log2(F[, "Pc.mu"]), Pe = log2(F[, "Pe.mu"]), Hc = log2(F[, "Hc.mu"] / (2 * C.BAR)), He = log2(F[, "He.mu"] / (2 * (1 - C.BAR))))
  lcr <- lapply(setNames(names(LEAK.MATS), names(LEAK.MATS)), function(k) ifelse(F[, paste0(k, ".phi_cr")] > 0, log2(F[, paste0(k, ".phi_cr")]), NA_real_))
  lmo <- lapply(setNames(names(LEAK.MATS), names(LEAK.MATS)), function(k) ifelse(F[, paste0(k, ".phi_mom")] > 0, log2(F[, paste0(k, ".phi_mom")]), NA_real_))
  mbar <- rowMeans(F[, paste0(names(LEAK.MATS), ".mbar")])
  terc <- cut(rank(mbar, ties.method = "first"), 3, labels = c("low", "mid", "high"))
  for (md in c("cis", "trans")) {
    if (md == "cis") { x <- lmu$Hc - lmu$He; ycr <- lcr$Hc - lcr$He; ymo <- lmo$Hc - lmo$He }
    else { x <- (lmu$Pc - lmu$Pe) - (lmu$Hc - lmu$He); ycr <- (lcr$Pc - lcr$Pe) - (lcr$Hc - lcr$He); ymo <- (lmo$Pc - lmo$Pe) - (lmo$Hc - lmo$He) }
    y <- ycr - ymo
    L2.DATA[[paste(arm, md)]] <- data.frame(gene = LEAK.GENES, arm = arm, contrast = md, dlog2mu = x, diff = y, dlog2phi_cr = ycr, dlog2phi_mom = ymo, tercile = terc)
    L2[[paste(arm, md, "all")]] <- data.frame(arm = arm, contrast = md, genes_in = "all", leak_regression(x, y, ycr, ymo, LEAK.BOOT, SEED.LEAK + 2))
    for (t in levels(terc)) { i <- terc == t
      L2[[paste(arm, md, t)]] <- data.frame(arm = arm, contrast = md, genes_in = paste("mean count tercile", t), leak_regression(x[i], y[i], ycr[i], ymo[i], LEAK.BOOT, SEED.LEAK + 3)) }
  }
}
L2 <- do.call(rbind, L2); rownames(L2) <- NULL; L2.DATA <- do.call(rbind, L2.DATA)
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
