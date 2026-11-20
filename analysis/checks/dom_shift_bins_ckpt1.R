###############################################################
### dom_shift_bins_ckpt1.R
### MID5: is one median per expression decile enough to describe the hybrid's genome-wide rise in dispersion? The dominance contrasts of log2 phi,
### Delta_dom = mid ((log2 phi_Hc + log2 phi_He) / 2 - (log2 phi_Pc + log2 phi_Pe) / 2) and par ((log2 phi_Hc + log2 phi_He) / 2 - log2 phi_P), are
### plotted against log2 mu_H, the per-genome hybrid mean on the log2 scale ((log2 mu_Hc + log2 mu_He) / 2, allele means on their own library offsets),
### with the median of each equal-count decile of log2 mu_H overlaid and a loess fit for comparison. Equal counts keep every median stable;
### the end deciles then span wide ranges of expression, so each decile is also split at its median log2 mu_H, and the medians of its two halves are
### compared: a difference beyond 2 SE (SE of a median about 1.2533 MAD / sqrt(n)) marks a decile whose shift slopes, which one median would under- or
### over-correct at its edges (split it, or use a smooth fit). With 30 deciles tested at 2 SE, about 1 to 2 flags arise by chance; a slope that matters is
### large against its SE and agrees with the loess fit in the figure.
### Fits: the pipeline's Cox-Reid NB on all cells of each group, offsets of Section 2, the library split correction (SPLIT.ETA). Genes need phi > 0 and a mean count
### of at least MID5.MIN.COUNT per cell in all four groups. Section 1 checkpoint only. Point estimates.
### Writes results/checks/MID5_dom_bins.csv, MID5.rds, a row in summary.csv and figures/extra/check_MID5_dom_shift.pdf.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/dom_shift_bins_ckpt1.R"
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
if (!exists("MID5.MIN.COUNT")) MID5.MIN.COUNT <- 0.1
if (!exists("MID5.BINS"))      MID5.BINS      <- 10    # equal-count bins of log2 mu_H

stopifnot("dom_shift_bins_ckpt1.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

## Fits of the four groups (as MID1): per-genome mean rate and split-corrected phi
DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
M5.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
M5.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF, Pe = colSums(MIX.SE) / DEPTH.REF, Hc = colSums(HYB.SC) / DEPTH.REF, He = colSums(HYB.SE) / DEPTH.REF)
M5.ETA   <- list(Pc = 0, Pe = 0, Hc = SPLIT.ETA$Hc, He = SPLIT.ETA$He)
M5.GENE.IDS <- rownames(MIX.SC)[Reduce(`&`, lapply(M5.MATS, function(m) rowMeans(m) >= MID5.MIN.COUNT))]
M5.FIT <- lapply(setNames(names(M5.MATS), names(M5.MATS)), function(k) {
  f <- do.call(rbind, mclapply(M5.GENE.IDS, function(g) fit_dispersion(as.numeric(M5.MATS[[k]][g, ]), M5.EXPOS[[k]], M5.ETA[[k]], "NB")[c("mu", "disp")], mc.cores = CORES))
  cbind(mu = f[, "mu"], phi = f[, "disp"])
})

## Per gene: log2 mu_H and the three dominance contrasts of log2 phi, on the genes with phi > 0 in every group
M5.GENES <- local({
  lm <- sapply(M5.FIT, function(f) log2(f[, "mu"])); lp <- suppressWarnings(sapply(M5.FIT, function(f) log2(f[, "phi"])))
  ok <- rowSums(is.finite(lm) & is.finite(lp)) == 4
  hp <- 0.5 * (lp[, "Hc"] + lp[, "He"])
  data.frame(gene = M5.GENE.IDS, log2_mu_H = 0.5 * (lm[, "Hc"] + lm[, "He"]), mid = hp - 0.5 * (lp[, "Pc"] + lp[, "Pe"]), par_sc = hp - lp[, "Pc"], par_se = hp - lp[, "Pe"],
             stringsAsFactors = FALSE)[ok, ]
})
cat(sprintf("%d genes with phi > 0 in every group (of %d with a mean count of at least %g)\n", nrow(M5.GENES), length(M5.GENE.IDS), MID5.MIN.COUNT))

## Equal-count bins of log2 mu_H; per bin the median shift, its range of log2 mu_H, and the medians of the two halves of the bin with the slope flag
M5.BIN <- cut(M5.GENES$log2_mu_H, quantile(M5.GENES$log2_mu_H, seq(0, 1, length.out = MID5.BINS + 1)), include.lowest = TRUE, labels = FALSE)
MID5.BINS.TAB <- do.call(rbind, lapply(c("mid", "par_sc", "par_se"), function(q) do.call(rbind, lapply(seq_len(MID5.BINS), function(b) {
  i <- M5.BIN == b; x <- M5.GENES$log2_mu_H[i]; y <- M5.GENES[[q]][i]; lo <- x <= median(x)
  d <- median(y[!lo]) - median(y[lo]); se <- sqrt(median_se(y[lo])^2 + median_se(y[!lo])^2)
  data.frame(contrast = q, bin = b, genes = sum(i), log2_mu_H_from = min(x), log2_mu_H_to = max(x), median_shift = median(y), se_median = median_se(y),
             median_lower_half = median(y[lo]), median_upper_half = median(y[!lo]), half_difference = d, half_difference_se = se,
             sloped = abs(d) > 2 * se, row.names = NULL)
}))))
cat("\nShift (log2 phi) per equal-count decile of log2 mu_H, with the medians of each decile's lower and upper half:\n")
print(MID5.BINS.TAB, digits = 3, row.names = FALSE)

## Figure: each dominance contrast against log2 mu_H, decile medians as horizontal segments over each decile's range, the half medians of every decile,
## a loess fit, and the deciles whose halves differ by more than 2 SE marked in red
dir.create(CHECK.FIG.DIR, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(CHECK.FIG.DIR, "check_MID5_dom_shift.pdf"), width = 12, height = 4.2, useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.6, 0))
for (q in c("mid", "par_sc", "par_se")) {
  x <- M5.GENES$log2_mu_H; y <- M5.GENES[[q]]; tb <- MID5.BINS.TAB[MID5.BINS.TAB$contrast == q, ]
  plot(x, y, pch = 16, cex = 0.3, col = adjustcolor("grey40", 0.35), xlab = expression(log[2] ~ mu[H]), ylab = expression(Delta[dom] ~ (log[2] ~ phi)),
       main = q, ylim = quantile(y, c(0.01, 0.99)))
  abline(h = 0, lty = 3)
  lf <- loess(y ~ x, span = 0.3); xs <- seq(min(x), max(x), length.out = 200); lines(xs, predict(lf, xs), col = "steelblue", lwd = 2)
  for (j in seq_len(nrow(tb))) {
    cl <- if (tb$sloped[j]) "firebrick" else "black"
    segments(tb$log2_mu_H_from[j], tb$median_shift[j], tb$log2_mu_H_to[j], tb$median_shift[j], col = cl, lwd = 2.5)
    xm <- median(x[M5.BIN == j])
    points(c((tb$log2_mu_H_from[j] + xm) / 2, (xm + tb$log2_mu_H_to[j]) / 2), c(tb$median_lower_half[j], tb$median_upper_half[j]), pch = 1, col = cl, cex = 0.9)
  }
  legend("topright", legend = c("decile median", "half-decile medians", "decile with sloped halves", "loess"), col = c("black", "black", "firebrick", "steelblue"),
         lty = c(1, NA, 1, 1), pch = c(NA, 1, NA, NA), lwd = c(2.5, NA, 2.5, 2), bty = "n", cex = 0.75)
}
dev.off()

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(MID5.BINS.TAB, file.path(CHECK.DIR, "MID5_dom_bins.csv"), row.names = FALSE)
check_save("MID5", list(bins = MID5.BINS.TAB, genes = M5.GENES, settings = list(min_count = MID5.MIN.COUNT, bins = MID5.BINS)), "run",
           headline = c(sprintf("%s: decile medians %.2f to %.2f; deciles with sloped halves: %s", c("mid", "par_sc", "par_se"),
                                tapply(MID5.BINS.TAB$median_shift, MID5.BINS.TAB$contrast, min)[c("mid", "par_sc", "par_se")],
                                tapply(MID5.BINS.TAB$median_shift, MID5.BINS.TAB$contrast, max)[c("mid", "par_sc", "par_se")],
                                vapply(c("mid", "par_sc", "par_se"), function(q) { s <- MID5.BINS.TAB$bin[MID5.BINS.TAB$contrast == q & MID5.BINS.TAB$sloped]
                                  if (length(s)) paste(s, collapse = ", ") else "none" }, ""))),
           mc_error = "point estimates; SE of a median 1.2533 MAD / sqrt(n)", script = CHECK.SCRIPT)
