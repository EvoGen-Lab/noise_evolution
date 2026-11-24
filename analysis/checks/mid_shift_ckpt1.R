###############################################################
### mid_shift_ckpt1.R
### MID1 point estimates from the Section 1 checkpoint alone: the median midparent contrast d (the per-genome hybrid against the midparent, both geometric
### means on the log2 scale) of the log2 mean and the log2 dispersion over calibrator genes, the genes with no trans and no total divergence, whose d
### should sit near zero unless culture or ploidy moves every gene alike.
### Without the Section 2 tests there are no q-values, so a calibrator is chosen by its point estimates: |total| < tau and |trans| < tau on the log2
### scale of the quantity, for every tau in MID1.TAU, so the dependence of the median on the choice shows. Point estimates only; no intervals.
### Fits: the pipeline's Cox-Reid NB estimator (fit_dispersion()) on all cells of each group, with the offsets of analysis.R Section 2 (a parent's library
### over the median parental library; the library of its own genome over the same reference for each hybrid allele, so the allele means are per-genome
### rates), and the library split correction of the allele dispersions with the split variances SPLIT.ETA of each cell. Contrasts on all hybrid cells (no cis and trans
### partition): total = log2 Pc - log2 Pe, cis = log2 Hc - log2 He, trans = total - cis, d = (log2 Hc + log2 He) / 2 - (log2 Pc + log2 Pe) / 2.
### A gene enters a quantity when all four groups have a finite positive value (mean, or phi > 0) and a mean count of at least MID1.MIN.COUNT per cell.
### Writes results/checks/MID1_ckpt1.csv (summary), MID1_ckpt1_genes.csv (per gene), MID1_ckpt1.rds and a row in summary.csv.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/mid_shift_ckpt1.R"
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
if (!exists("MID1.TAU"))       MID1.TAU       <- c(0.1, 0.25, 0.5)   # log2 thresholds on |total| and |trans| that define a calibrator
if (!exists("MID1.MIN.COUNT")) MID1.MIN.COUNT <- 0.1                 # minimum mean count per cell in every group

stopifnot("mid_shift_ckpt1.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

## Offsets and split variances of Section 2
DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
M1.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
M1.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF, Pe = colSums(MIX.SE) / DEPTH.REF, Hc = colSums(HYB.SC) / DEPTH.REF, He = colSums(HYB.SE) / DEPTH.REF)
M1.ETA   <- list(Pc = 0, Pe = 0, Hc = SPLIT.ETA$Hc, He = SPLIT.ETA$He)

## Fits of every gene in the four groups: per-genome mean rate mu, dispersion phi and mean count per cell
M1.FIT <- lapply(names(M1.MATS), function(k) {
  f <- do.call(rbind, mclapply(seq_len(nrow(M1.MATS[[k]])), function(i) fit_dispersion(M1.MATS[[k]][i, ], M1.EXPOS[[k]], M1.ETA[[k]], "NB")[c("mu", "disp")],
                               mc.cores = CORES))
  cbind(mu = f[, "mu"], phi = f[, "disp"], count = rowMeans(M1.MATS[[k]]))
})
names(M1.FIT) <- names(M1.MATS)
cat(sprintf("fitted %d genes in each of Pc, Pe, Hc, He\n", nrow(M1.FIT$Pc)))

## Per-gene log2 contrasts of the mean and the dispersion
M1.GENES <- local({
  mu <- sapply(M1.FIT, function(f) f[, "mu"]); ph <- sapply(M1.FIT, function(f) f[, "phi"]); ct <- sapply(M1.FIT, function(f) f[, "count"])
  enough <- rowSums(ct >= MID1.MIN.COUNT) == 4
  ok_mu <- enough & rowSums(is.finite(mu) & mu > 0) == 4
  ok_ph <- enough & rowSums(is.finite(ph) & ph > 0) == 4
  lm <- suppressWarnings(log2(mu)); lp <- suppressWarnings(log2(ph))
  out <- data.frame(gene = rownames(MIX.SC), mean_count = rowMeans(ct[, c("Pc", "Pe")]), ok_mean = ok_mu, ok_disp = ok_ph, stringsAsFactors = FALSE)
  for (q in c("mean", "disp")) {
    l <- if (q == "mean") lm else lp; ok <- if (q == "mean") ok_mu else ok_ph
    tot <- l[, "Pc"] - l[, "Pe"]; cis <- l[, "Hc"] - l[, "He"]
    out[[paste0(q, "_total")]] <- ifelse(ok, tot, NA); out[[paste0(q, "_cis")]] <- ifelse(ok, cis, NA)
    out[[paste0(q, "_trans")]] <- ifelse(ok, tot - cis, NA)
    out[[paste0(q, "_mid")]]   <- ifelse(ok, 0.5 * (l[, "Hc"] + l[, "He"]) - 0.5 * (l[, "Pc"] + l[, "Pe"]), NA)
  }
  out
})

## Median d of the calibrators at each threshold, against the spread of d over all genes of the quantity
MID1.CKPT1 <- do.call(rbind, lapply(c("mean", "disp"), function(q) {
  d <- M1.GENES[[paste0(q, "_mid")]]; ok <- is.finite(d)
  do.call(rbind, lapply(MID1.TAU, function(tau) {
    cal <- ok & abs(M1.GENES[[paste0(q, "_total")]]) < tau & abs(M1.GENES[[paste0(q, "_trans")]]) < tau
    data.frame(quantity = q, tau = tau, genes = sum(ok), calibrators = sum(cal), median_d_calibrators = median(d[cal]), iqr_d_calibrators = IQR(d[cal]),
               median_d_all = median(d[ok]), median_abs_d_all = median(abs(d[ok])), q90_abs_d_all = unname(quantile(abs(d[ok]), 0.9)), row.names = NULL)
  }))
}))
cat("\nMidparent contrast d (log2), calibrators with |total| and |trans| below tau:\n"); print(MID1.CKPT1, digits = 3, row.names = FALSE)

## The calibrator median by tercile of mean expression, at the middle threshold
MID1.TERC <- do.call(rbind, lapply(c("mean", "disp"), function(q) {
  tau <- MID1.TAU[ceiling(length(MID1.TAU) / 2)]
  cal <- is.finite(M1.GENES[[paste0(q, "_mid")]]) & abs(M1.GENES[[paste0(q, "_total")]]) < tau & abs(M1.GENES[[paste0(q, "_trans")]]) < tau
  m <- M1.GENES$mean_count[cal]; d <- M1.GENES[[paste0(q, "_mid")]][cal]
  tc <- cut(m, quantile(m, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))
  do.call(rbind, lapply(levels(tc), function(t) data.frame(quantity = q, tau = tau, tercile = t, calibrators = sum(tc == t),
                                                           median_d = median(d[tc == t]), row.names = NULL)))
}))
cat("\nCalibrator median d by tercile of mean count:\n"); print(MID1.TERC, digits = 3, row.names = FALSE)

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(MID1.CKPT1, file.path(CHECK.DIR, "MID1_ckpt1.csv"), row.names = FALSE)
write.csv(M1.GENES, file.path(CHECK.DIR, "MID1_ckpt1_genes.csv"), row.names = FALSE)
check_save("MID1_ckpt1", list(summary = MID1.CKPT1, terciles = MID1.TERC, genes = M1.GENES, settings = list(tau = MID1.TAU, min_count = MID1.MIN.COUNT)),
           "reduced", headline = sprintf("%s, tau %g: %d calibrators, median d %.3f (all genes: median d %.3f, median |d| %.3f)", MID1.CKPT1$quantity,
                                         MID1.CKPT1$tau, MID1.CKPT1$calibrators, MID1.CKPT1$median_d_calibrators, MID1.CKPT1$median_d_all, MID1.CKPT1$median_abs_d_all),
           mc_error = "point estimates only; calibrators by |total| and |trans| below tau, no tests", script = CHECK.SCRIPT)
