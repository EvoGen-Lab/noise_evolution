###############################################################
### state_conditioned.R
### State-conditioned dispersion: how much of each gene's dispersion, and of the cis and trans noise contrasts, comes from cell states that other genes
### mark. Runs from the Section 1 checkpoint alone (no Section 2 rerun and no cluster output), with the offsets and share corrections of Section 2 as in
### estimator_leak.R.
### No cells are removed and no cells are clustered. Each gene is fitted twice by the same Cox-Reid NB fit (state_cr_glm()): with a constant mean (phi0,
### the total dispersion) and with continuous state scores in the mean (phiZ, the dispersion given the states). Clustering cells on expression and then
### measuring expression variance inside the clusters removes variance by construction; here a gene's own counts never enter its state scores, so the
### scores carry only the variation the gene shares with its module. phiZ is a sensitivity estimate reported beside the total, not a replacement: both
### describe real variation, at two levels of conditioning.
### State scores: each module's summed counts per cell over the offset, log2, standardized within the group, the gene left out of its own module. In the
### hybrid cells the score sums both alleles, so the state is the cell's and the two alleles share it. Modules (GO through org.Sc.sgd.db, as in Section 7):
### glycolysis GO:0006096, oxidative phosphorylation GO:0006119, ribosome biogenesis GO:0042254, non-fermentative carbon use (glyoxylate cycle GO:0006097
### and fatty-acid beta-oxidation GO:0006635) and S phase (nucleosome GO:0000786, histone genes).
###   SC1  the state effects: per group and module the median coefficient, the share of genes whose likelihood gain exceeds the chi-square 95% point, and
###        the median log2(phiZ / phi0) of the module's own genes;
###   SC2  log2(phiZ / phi0) per group, by mean-count bin, and the genes whose phi falls to zero once the states are in the mean;
###   SC3  cis and trans contrasts of log2 phi, total against state-conditioned: Spearman correlation, median absolute change, sign changes among
###        contrasts above 0.5 in absolute value, and the own slopes of Delta log2 phi on Delta log2 mu (the coupling) for both, with the slope of their
###        difference and its gene-bootstrap interval.
### A score measured with sampling error conditions only partly, so phiZ stays above the dispersion given the true state (errors-in-variables); low-depth
### cells measure their state least well.
### Results go to results/checks/SC1.rds to SC3.rds with rows in summary.csv; the figure goes to figures/extra/check_SC3_state_conditioned.pdf.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel); library(org.Sc.sgd.db)

CHECK.SCRIPT     <- "analysis/checks/state_conditioned.R"
CORES            <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.STATE       <- 8026
STATE.MIN.MEAN   <- 0.02     # mean count per cell required in every group (as estimator_leak.R)
STATE.MIN.CELLS  <- 50       # cells with a count required in every group
STATE.MAX.GENES  <- Inf      # a random subset of the eligible genes when finite (for a quick run)
STATE.BOOT       <- 1000     # gene-resampling draws for the slope intervals
STATE.SETS       <- list(Glycolysis = "GO:0006096", OXPHOS = "GO:0006119", RiBi = "GO:0042254", NonFerm = c("GO:0006097", "GO:0006635"), Sphase = "GO:0000786")
MBAR.BREAKS      <- c(0, 0.1, 0.3, 1, 3, 10, Inf)
set.seed(SEED.SECTION + 24)

stopifnot("state_conditioned.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE (cells after QC, genes after the union filter), C.BAR, R.C, R.E

## Offsets and share corrections of Section 2 (as estimator_leak.R).
DEPTH.REF.S <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
STATE.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
STATE.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF.S, Pe = colSums(MIX.SE) / DEPTH.REF.S,
                    Hc = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF.S / 2, He = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF.S / 2)
STATE.ETA   <- c(Pc = 0, Pe = 0, Hc = R.C, He = R.E)

## Module members present in the data, and each module's summed counts per cell (both alleles together in the hybrid cells).
STATE.MEMBERS <- lapply(STATE.SETS, function(ids) intersect(unique(unlist(lapply(ids, go_gene_set))), rownames(MIX.SC)))
stopifnot("every state module needs at least 5 genes in the data" = all(lengths(STATE.MEMBERS) >= 5))
HYB.TOT <- vapply(STATE.MEMBERS, function(m) colSums(HYB.SC[m, , drop = FALSE]) + colSums(HYB.SE[m, , drop = FALSE]), numeric(ncol(HYB.SC)))
STATE.TOTALS <- list(Pc = vapply(STATE.MEMBERS, function(m) colSums(MIX.SC[m, , drop = FALSE]), numeric(ncol(MIX.SC))),
                     Pe = vapply(STATE.MEMBERS, function(m) colSums(MIX.SE[m, , drop = FALSE]), numeric(ncol(MIX.SE))), Hc = HYB.TOT, He = HYB.TOT)
cat("State modules (genes in the data):", paste(sprintf("%s %d", names(STATE.MEMBERS), lengths(STATE.MEMBERS)), collapse = ", "), "\n")
print(round(cor(STATE.TOTALS$Pc / STATE.EXPOS$Pc, method = "spearman"), 2))   # how distinct the module scores are in one group

## Genes with data in all four groups (as estimator_leak.R).
STATE.GENES <- rownames(MIX.SC)[Reduce(`&`, lapply(STATE.MATS, function(m) rowMeans(m) >= STATE.MIN.MEAN & rowSums(m > 0) >= STATE.MIN.CELLS))]
if (is.finite(STATE.MAX.GENES) && length(STATE.GENES) > STATE.MAX.GENES) STATE.GENES <- with_local_seed(SEED.STATE, sample(STATE.GENES, STATE.MAX.GENES))
cat(sprintf("State-conditioned fits: %d genes; cells Pc %d, Pe %d, hybrid %d\n", length(STATE.GENES), ncol(MIX.SC), ncol(MIX.SE), ncol(HYB.SC)))

STATE.FITS <- do.call(rbind, mclapply(STATE.GENES, state_gene, mats = STATE.MATS, expos = STATE.EXPOS, eta = STATE.ETA, totals = STATE.TOTALS,
                                      members = STATE.MEMBERS, mc.cores = CORES, mc.preschedule = FALSE))
rownames(STATE.FITS) <- STATE.GENES
GROUPS <- names(STATE.MATS)
STATE.COLS <- c("mu", "mbar", "phi_pipe", "phi0", "phiZ", "lr", paste0("beta_", names(STATE.MEMBERS)))
FIT <- lapply(setNames(GROUPS, GROUPS), function(k) setNames(as.data.frame(STATE.FITS[, paste(k, STATE.COLS, sep = "."), drop = FALSE]), STATE.COLS))   # per-group columns

## Agreement of the intercept-only fit with the pipeline fit (they differ only in the mean: NB maximum likelihood against the ratio estimate).
SC.AGREE <- vapply(FIT, function(f) median(abs(log2(f$phi0 / f$phi_pipe))[f$phi0 > 0 & f$phi_pipe > 0], na.rm = TRUE), numeric(1))
cat("Median |log2(phi0 / phi pipeline)| per group:", paste(sprintf("%s %.4f", GROUPS, SC.AGREE), collapse = ", "), "\n")

##############################################################################
## SC1. State effects: per group and module the median coefficient (change in log mean per SD of the score) and the share of genes whose likelihood
## gain from all modules together exceeds the chi-square 95% point with one degree of freedom per module.
SC1 <- list()
for (k in GROUPS) {
  B <- FIT[[k]][paste0("beta_", names(STATE.MEMBERS))]
  SC1[[k]] <- data.frame(group = k, module = names(STATE.MEMBERS), genes = sum(is.finite(FIT[[k]]$lr)), median_beta = vapply(B, median, numeric(1), na.rm = TRUE),
                         median_abs_beta = vapply(abs(B), median, numeric(1), na.rm = TRUE),
                         share_lr_above_crit = mean(FIT[[k]]$lr > qchisq(0.95, length(STATE.MEMBERS)), na.rm = TRUE),
                         member_genes = vapply(STATE.MEMBERS, function(m) sum(STATE.GENES %in% m), numeric(1)),
                         members_median_log2_phiZ_over_phi0 = vapply(STATE.MEMBERS, function(m) median(log2(FIT[[k]]$phiZ / FIT[[k]]$phi0)[STATE.GENES %in% m & FIT[[k]]$phiZ > 0 & FIT[[k]]$phi0 > 0]), numeric(1)),
                         row.names = NULL)
}
SC1 <- do.call(rbind, SC1)
print(SC1, digits = 3, row.names = FALSE)
check_save("SC1", list(table = SC1, agreement = SC.AGREE, members = STATE.MEMBERS), "run",
           headline = c(sprintf("%s %s: median coefficient %.3f (median |coefficient| %.3f); likelihood gain above chi-square 95%% in %.3f of %d genes; module's own %d genes: median log2(phiZ / phi0) %.3f",
                                SC1$group, SC1$module, SC1$median_beta, SC1$median_abs_beta, SC1$share_lr_above_crit, SC1$genes, SC1$member_genes, SC1$members_median_log2_phiZ_over_phi0),
                        sprintf("median |log2(phi0 / phi pipeline)|: %s", paste(sprintf("%s %.4f", GROUPS, SC.AGREE), collapse = ", "))),
           mc_error = "none (deterministic fits)", script = CHECK.SCRIPT)

##############################################################################
## SC2. log2(phiZ / phi0) per group and mean-count bin, and the genes whose phi falls to zero once the states are in the mean
SC2 <- list()
for (k in GROUPS) {
  p0 <- FIT[[k]]$phi0; pz <- FIT[[k]]$phiZ; bins <- as.character(cut(FIT[[k]]$mbar, MBAR.BREAKS, right = FALSE)); r <- log2(pz / p0)
  for (b in c("all", levels(cut(0, MBAR.BREAKS, right = FALSE)))) {
    i <- which(p0 > 0 & (b == "all" | bins == b))
    if (length(i) < 10) next
    SC2[[paste(k, b)]] <- data.frame(group = k, mean_count = b, genes = length(i), median_log2_phiZ_over_phi0 = median(r[i][pz[i] > 0], na.rm = TRUE),
                                     q10 = unname(quantile(r[i][pz[i] > 0], 0.1, na.rm = TRUE)), share_phiZ_zero = mean(pz[i] == 0, na.rm = TRUE), row.names = NULL)
  }
}
SC2 <- do.call(rbind, SC2); rownames(SC2) <- NULL
print(SC2, digits = 3, row.names = FALSE)
check_save("SC2", SC2, "run",
           headline = sprintf("%s mean count %s: median log2(phiZ / phi0) %.3f (10th percentile %.3f); phi to zero in %.3f of %d genes",
                              SC2$group, SC2$mean_count, SC2$median_log2_phiZ_over_phi0, SC2$q10, SC2$share_phiZ_zero, SC2$genes),
           mc_error = "none (deterministic fits)", script = CHECK.SCRIPT)

##############################################################################
## SC3. Cis and trans contrasts of log2 phi, total against state-conditioned. Allele means are divided by their abundance factors (as leak_contrasts()).
## The own slope is the slope of Delta log2 phi on Delta log2 mu; leak_regression() takes the total and the conditioned contrasts as its two estimators,
## so its slope is that of (total - conditioned) on Delta log2 mu, with a gene-bootstrap interval.
LMU <- list(Pc = log2(FIT[["Pc"]][["mu"]]), Pe = log2(FIT[["Pe"]][["mu"]]), Hc = log2(FIT[["Hc"]][["mu"]] / (2 * C.BAR)), He = log2(FIT[["He"]][["mu"]] / (2 * (1 - C.BAR))))
LPH <- lapply(c(total = "phi0", state = "phiZ"), function(q) lapply(setNames(GROUPS, GROUPS), function(k) ifelse(FIT[[k]][[q]] > 0, log2(FIT[[k]][[q]]), NA_real_)))
CONTRAST.W <- list(cis = c(Hc = 1, He = -1), trans = c(Pc = 1, Pe = -1, Hc = -1, He = 1))   # cis: allele log2 ratio; trans: parental minus allele log2 ratio
SC3 <- list(); SC3.PTS <- list()
for (md in c("cis", "trans")) {
  w <- CONTRAST.W[[md]]
  x <- Reduce(`+`, Map(`*`, LMU[names(w)], w)); a <- Reduce(`+`, Map(`*`, LPH$total[names(w)], w)); b <- Reduce(`+`, Map(`*`, LPH$state[names(w)], w))
  ok <- is.finite(x) & is.finite(a) & is.finite(b)
  big <- ok & abs(a) > 0.5
  reg <- leak_regression(x[ok], a[ok] - b[ok], a[ok], b[ok], STATE.BOOT, SEED.STATE + 1)
  SC3[[md]] <- data.frame(contrast = md, genes = sum(ok), spearman_total_state = cor(a[ok], b[ok], method = "spearman"),
                          median_abs_change = median(abs(a[ok] - b[ok])), large_contrasts = sum(big), sign_changes_among_large = mean(sign(a[big]) != sign(b[big])),
                          own_slope_total = reg$slope_cr, own_slope_state = reg$slope_mom, slope_difference = reg$slope, boot_lo = reg$boot_lo, boot_hi = reg$boot_hi, row.names = NULL)
  SC3.PTS[[md]] <- data.frame(gene = STATE.GENES[ok], contrast = md, dlog2mu = x[ok], total = a[ok], state = b[ok])
}
SC3 <- do.call(rbind, SC3); SC3.PTS <- do.call(rbind, SC3.PTS); rownames(SC3) <- rownames(SC3.PTS) <- NULL
print(SC3, digits = 3, row.names = FALSE)
check_save("SC3", list(summary = SC3, points = SC3.PTS), "run",
           headline = sprintf("%s: total vs state-conditioned Delta log2 phi Spearman %.3f, median |change| %.3f, sign changes in %.3f of %d contrasts above 0.5; own slopes %.3f (total) and %.3f (state), difference %.3f (%.3f to %.3f), %d genes",
                              SC3$contrast, SC3$spearman_total_state, SC3$median_abs_change, SC3$sign_changes_among_large, SC3$large_contrasts,
                              SC3$own_slope_total, SC3$own_slope_state, SC3$slope_difference, SC3$boot_lo, SC3$boot_hi, SC3$genes),
           mc_error = sprintf("slope interval from %d gene-resampling draws", STATE.BOOT), script = CHECK.SCRIPT)

fig_pdf("extra/check_SC3_state_conditioned.pdf", 9, 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (md in c("cis", "trans")) { p <- SC3.PTS[SC3.PTS$contrast == md, ]
  plot(p$total, p$state, pch = 19, cex = 0.3, col = adjustcolor(COLOR.GREY[["dark"]], 0.4), xlab = "Delta log2 phi, total", ylab = "Delta log2 phi, state-conditioned", main = md)
  abline(0, 1, col = COLOR.ACCENT) }
dev.off()
