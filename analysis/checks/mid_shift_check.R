###############################################################
### mid_shift_check.R
### MID1: the size of a design shift in the midparent contrast d (the per-genome hybrid against the mean of the two parents, mode "mid"), on the log2 mean
### and the log2 dispersion. The parents grew in one culture and the hybrid in another, and the hybrid carries both genomes, so culture and ploidy can move d
### for every gene alike. Genes with no detectable trans divergence and no detectable total divergence are an internal calibrator: their two parental trans
### environments act on the gene alike, so their d should sit near zero unless the design moves it. The median d over those genes estimates the design shift
### with little of the real trans dominance mixed in.
### Per quantity (mean, disp): the calibrator genes (matched set GENES, finite d, BH q >= MID.SIG for both the trans and the total contrast of that quantity),
### their median d with two 95% intervals, and the scale to read it against:
###   gene bootstrap  MID.B.GENE resamples of the calibrator genes, the median recomputed in each (the spread a new set of such genes would give);
###   cell bootstrap  MID.B.CELL resamples of the parental and undivided hybrid cells, with the four groups refitted (the split correction on the
###                   resampled cells), the median of a
###                   random subsample of MID.CELL.GENES calibrator genes recomputed in each replicate (mid_shift_gene()); this is the replicate-by-replicate
###                   median an adjustment would subtract;
###   scale           the genome-wide median d, the median |d| of the genes whose midparent contrast is significant (q < MID.SIG), the IQR of d over the
###                   calibrators and their median bootstrap SE of d.
### The calibrator median is also given by tercile of mean expression, since a capture or culture shift on the log2 mean should not depend on expression.
### Decision read from the table: a calibrator median small against the typical significant d and with an interval that covers zero leaves d unadjusted;
### otherwise subtracting it (recomputed in each replicate) removes the shared shift.
### d is recomputed here from the observed allele fits of all hybrid cells (CONTRAST.FITS) with the current contrast definition (geometric midpoints on the
### log2 scale); its PR q-values (trans, total, mid) and SEs come from the checkpoint as written (the pipeline's mid contrast uses the trans cells).
### Needs the Section 2 checkpoint (GENES, BURST.CONTRASTS, PR, CONTRAST.FITS, CONTRAST.MATS, CONTRAST.EXPOS, SPLIT.FIT.ETAS). Runs locally; no cluster output.
### Results: results/checks/MID1.rds, results/checks/MID1_mid_shift.csv, a row in summary.csv, figures/extra/check_MID1_mid_shift.pdf.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/mid_shift_check.R"
CORES        <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.MID     <- 4026
if (!exists("MID.SIG"))        MID.SIG        <- 0.05   # BH q threshold: calibrators have q >= MID.SIG for trans and total; "significant d" has q < MID.SIG
if (!exists("MID.B.GENE"))     MID.B.GENE     <- 2000   # gene-bootstrap resamples (fast)
if (!exists("MID.B.CELL"))     MID.B.CELL     <- 200    # cell-bootstrap replicates (each refits 4 groups per subsample gene)
if (!exists("MID.CELL.GENES")) MID.CELL.GENES <- 200    # calibrator genes in the cell bootstrap; runtime grows with MID.B.CELL x MID.CELL.GENES
set.seed(SEED.SECTION + 30)

stopifnot("mid_shift_check.R needs the Section 2 checkpoint" = file.exists(ckpt_path(2)))
load(ckpt_path(2))
need <- c("GENES", "BURST.CONTRASTS", "PR", "CONTRAST.FITS", "CONTRAST.MATS", "CONTRAST.EXPOS", "SPLIT.FIT.ETAS")
stopifnot("the Section 2 checkpoint lacks an object this check needs" = all(vapply(need, exists, logical(1))))
MID.Q <- c("mean", "disp")
stopifnot(all(paste0(rep(MID.Q, each = 3), c("_trans_q", "_total_q", "_mid_q")) %in% names(PR)), identical(PR$gene, BURST.CONTRASTS$gene))

## Observed d of every matched-set gene from its fits, with the current definition; the allele fits of all hybrid cells take the place of the trans cells
MID.OBS <- t(vapply(GENES, function(g) { fv <- fit_group_values(setNames(CONTRAST.FITS[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], c("MIX.SC", "MIX.SE", "HYT.SC", "HYT.SE")), g)
  .mode_contrasts(fv$mu, fv$disp, modes = "mid")["mid", c("mean", "disp")] }, numeric(2)))

## Calibrator genes and the scale of d, per quantity
MID.SETS <- lapply(setNames(MID.Q, MID.Q), function(q) {
  d <- unname(MID.OBS[match(BURST.CONTRASTS$gene, GENES), q])
  in_g <- BURST.CONTRASTS$gene %in% GENES & is.finite(d)
  cal <- in_g & (PR[[paste0(q, "_trans_q")]] >= MID.SIG) %in% TRUE & (PR[[paste0(q, "_total_q")]] >= MID.SIG) %in% TRUE
  sig <- in_g & (PR[[paste0(q, "_mid_q")]] < MID.SIG) %in% TRUE
  list(d = setNames(d, BURST.CONTRASTS$gene), se = setNames(BURST.CONTRASTS[[paste0(q, "_mid_se")]], BURST.CONTRASTS$gene),
       genes = BURST.CONTRASTS$gene[cal], sig = BURST.CONTRASTS$gene[sig], all = BURST.CONTRASTS$gene[in_g])
})
for (q in MID.Q) cat(sprintf("%s: %d calibrator genes (trans and total q >= %.2f) of %d matched-set genes with a finite d; %d with a significant d\n",
                             q, length(MID.SETS[[q]]$genes), MID.SIG, length(MID.SETS[[q]]$all), length(MID.SETS[[q]]$sig)))

## Gene bootstrap of the calibrator median (indices drawn first, on their own stream)
MID.GENE.BOOT <- with_local_seed(SEED.MID, lapply(MID.Q, function(q) { x <- MID.SETS[[q]]$d[MID.SETS[[q]]$genes]
  replicate(MID.B.GENE, median(x[sample.int(length(x), length(x), replace = TRUE)])) }))
names(MID.GENE.BOOT) <- MID.Q

## Cell bootstrap of the median over a random subsample of calibrator genes. Every replicate's indices are drawn before any fit, and all genes share them,
## so the median over genes is taken within each replicate. The subsample is drawn from the genes that calibrate both quantities when there are enough,
## otherwise from the union.
MID.CELL.SET <- with_local_seed(SEED.MID + 1, {
  both <- intersect(MID.SETS$mean$genes, MID.SETS$disp$genes)
  pool <- if (length(both) >= MID.CELL.GENES) both else union(MID.SETS$mean$genes, MID.SETS$disp$genes)
  sort(sample(pool, min(MID.CELL.GENES, length(pool))))
})
MID.DRAWS <- with_local_seed(SEED.MID + 2, {
  n <- c(Pc = ncol(CONTRAST.MATS$MIX.SC), Pe = ncol(CONTRAST.MATS$MIX.SE), H = ncol(CONTRAST.MATS$HYB.SC))
  lapply(seq_len(MID.B.CELL), function(b) lapply(n, function(k) sample.int(k, k, replace = TRUE)))
})
t0 <- Sys.time()
MID.CELL.REPS <- mclapply(MID.CELL.SET, mid_shift_gene, mats = CONTRAST.MATS[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], expos = CONTRAST.EXPOS,
                          etas = SPLIT.FIT.ETAS[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], draws = MID.DRAWS, mc.cores = CORES)
cat(sprintf("cell bootstrap: %d genes x %d replicates in %.1f min\n", length(MID.CELL.SET), MID.B.CELL, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
## median over the subsample genes in each replicate, for each quantity; a gene enters a quantity's median only if it calibrates that quantity
MID.CELL.MED <- sapply(MID.Q, function(q) { k <- MID.CELL.SET %in% MID.SETS[[q]]$genes
  apply(vapply(MID.CELL.REPS[k], function(m) m[, q], numeric(MID.B.CELL)), 1, median, na.rm = TRUE) })

## The table
MID1 <- do.call(rbind, lapply(MID.Q, function(q) {
  s <- MID.SETS[[q]]; dc <- s$d[s$genes]; ci_g <- quantile(MID.GENE.BOOT[[q]], c(0.025, 0.975), names = FALSE)
  sub <- intersect(MID.CELL.SET, s$genes); ci_c <- quantile(MID.CELL.MED[, q], c(0.025, 0.975), type = 1, na.rm = TRUE, names = FALSE)
  typ <- median(abs(s$d[s$sig]))
  data.frame(quantity = q, calibrators = length(dc), median_d = median(dc), gene_lo = ci_g[1], gene_hi = ci_g[2],
             cell_genes = length(sub), cell_median_d = median(s$d[sub]), cell_lo = ci_c[1], cell_hi = ci_c[2],
             iqr_d = IQR(dc), median_se_d = median(s$se[s$genes], na.rm = TRUE), genomewide_median_d = median(s$d[s$all]),
             significant = length(s$sig), median_abs_sig_d = typ, shift_over_typical_sig = abs(median(dc)) / typ,
             gene_ci_covers_0 = ci_g[1] <= 0 && ci_g[2] >= 0, cell_ci_covers_0 = ci_c[1] <= 0 && ci_c[2] >= 0, row.names = NULL)
}))
cat("\nMidparent contrast d of the calibrator genes (log2):\n"); print(MID1, digits = 3, row.names = FALSE)

## The calibrator median by tercile of mean expression (the mean of the two parental fitted means)
MID.TERC <- do.call(rbind, lapply(MID.Q, function(q) {
  g <- MID.SETS[[q]]$genes; m <- 0.5 * (CONTRAST.FITS$MIX.SC[g, "MU"] + CONTRAST.FITS$MIX.SE[g, "MU"])
  tc <- cut(m, quantile(m, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("low", "mid", "high"))
  do.call(rbind, lapply(levels(tc), function(t) { x <- MID.SETS[[q]]$d[g[tc %in% t]]
    data.frame(quantity = q, tercile = t, genes = length(x), median_d = median(x), iqr_d = IQR(x), row.names = NULL) }))
}))
cat("\nCalibrator median d by tercile of mean expression:\n"); print(MID.TERC, digits = 3, row.names = FALSE)

## Figure: d of the calibrator genes, of all matched-set genes and of the significant genes, with the calibrator median and its gene-bootstrap interval
dir.create(CHECK.FIG.DIR, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(CHECK.FIG.DIR, "check_MID1_mid_shift.pdf"), width = 9, height = 4, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
for (q in MID.Q) {
  s <- MID.SETS[[q]]; r <- MID1[MID1$quantity == q, ]
  br <- pretty(range(s$d[s$all], na.rm = TRUE), 60)
  h_all <- hist(s$d[s$all], breaks = br, plot = FALSE); h_cal <- hist(s$d[s$genes], breaks = br, plot = FALSE)
  plot(h_all, col = COLOR.GREY[["light"]], border = NA, main = sprintf("%s: median d of calibrators %.3f", q, r$median_d), xlab = sprintf("d (log2 %s, midparent)", q))
  plot(h_cal, col = adjustcolor(COLOR.ACCENT, 0.5), border = NA, add = TRUE)
  abline(v = 0, lty = 3); abline(v = r$median_d, col = COLOR.ACCENT, lwd = 2); abline(v = c(r$gene_lo, r$gene_hi), col = COLOR.ACCENT, lty = 2)
  if (length(s$sig)) rug(s$d[s$sig], col = "black")
  legend("topright", legend = c("all matched-set genes", "calibrator genes", "significant d (rug)"), fill = c(COLOR.GREY[["light"]], adjustcolor(COLOR.ACCENT, 0.5), NA),
         border = NA, bty = "n", cex = 0.7)
}
dev.off()

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(MID1, file.path(CHECK.DIR, "MID1_mid_shift.csv"), row.names = FALSE)
check_save("MID1", list(table = MID1, terciles = MID.TERC, gene_boot = MID.GENE.BOOT, cell_median = MID.CELL.MED, cell_genes = MID.CELL.SET,
                        settings = list(sig = MID.SIG, B_gene = MID.B.GENE, B_cell = MID.B.CELL, cell_genes = MID.CELL.GENES)), "run",
           headline = sprintf("%s: %d calibrators, median d %.3f (gene 95%% %.3f to %.3f; cell 95%% %.3f to %.3f), typical significant |d| %.3f (ratio %.2f)",
                              MID1$quantity, MID1$calibrators, MID1$median_d, MID1$gene_lo, MID1$gene_hi, MID1$cell_lo, MID1$cell_hi,
                              MID1$median_abs_sig_d, MID1$shift_over_typical_sig),
           mc_error = sprintf("gene bootstrap %d resamples; cell bootstrap %d replicates on %d genes", MID.B.GENE, MID.B.CELL, length(MID.CELL.SET)),
           script = CHECK.SCRIPT)
