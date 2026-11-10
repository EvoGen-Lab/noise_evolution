###############################################################
### checks_data.R
### Data checks of the methods document (docs/methods.md), one block per check ID of the validation handoff: A3, A5, A7 (data arm), A8, B1, B2 (data
### arm), B5, B6, C1, C2, C3, C6, D9, E2, E3 (data arm), E4, F1(b, c), F5 and G1, and a local bootstrap of a gene subsample for A4, A6, C7 and D8.
### Reads the Section 1 and Section 2 checkpoints (results/checkpoints/section1_checkpoint.rda, section2_checkpoint.rda); a check whose inputs are
### missing is recorded as "not run" with the reason. The checks that read the bootstrap and null outputs (C3, D9, E2 significance counts, E3) use the
### cluster results already merged into the Section 2 checkpoint (BOOT.CONTRASTS, PR); the local bootstrap needs no cluster output.
### Each result goes to results/checks/<ID>.rds with a row in results/checks/summary.csv; figures go to figures/extra/.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT  <- "analysis/checks/checks_data.R"
CORES         <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)
SEED.CHECKS   <- 3026
N.CHECK.GENES <- 300      # genes in the local bootstrap and refit subsamples (random genes of the matched set, resolved in every dataset)
N.CHECK.BOOT  <- 200      # bootstrap replicates per gene in the local bootstrap
set.seed(SEED.SECTION + 21)

have1 <- file.exists(ckpt_path(1)); have2 <- file.exists(ckpt_path(2))
if (have1) { E1.ENV <- new.env(); load(ckpt_path(1), envir = E1.ENV) }
if (have2) load(ckpt_path(2))
if (!have2) {
  for (id in c("A3", "A5", "A7_data", "A8", "B1", "B2_data", "B5", "B6", "C1", "C2", "C3", "C6", "D9", "E2", "E3_data", "E4", "F1_data", "F5", "G1", "A4_data", "A6", "C7_data", "D8"))
    check_save(id, NULL, "not run", headline = "the Section 2 checkpoint is missing; run analysis.R Sections 1 and 2 first", script = CHECK.SCRIPT)
  stop("checks_data.R needs the Section 2 checkpoint")
}
FIT.GROUPS <- names(SPLIT.FITS)
SUB.GENES  <- with_local_seed(SEED.CHECKS, sample(GENES, min(N.CHECK.GENES, length(GENES))))
dir.create(CHECK.FIG.DIR, recursive = TRUE, showWarnings = FALSE)

##############################################################################
## A3. Fits whose maximum reaches either end of the search range (Methods: Estimating phi)
# LOG.THETA.RANGE c(-4, 15) spans phi from exp(-15) = 3.06e-7 to exp(4) = 54.6. A fit at the lower end has phi-hat near 3e-7 (the boundary rule sends
# phi = 0 fits to the boundary instead, with HIT_RANGE 0), one at the upper end phi-hat near 55.
A3 <- list()
for (k in FIT.GROUPS) {
  f <- SPLIT.FITS[[k]]; hit <- f$HIT_RANGE %in% 1
  A3[[k]] <- data.frame(dataset = k, genes = nrow(f), hit_range = sum(hit), lower_end = sum(hit & f$DISP_APP < 1e-5, na.rm = TRUE), upper_end = sum(hit & f$DISP_APP > 1, na.rm = TRUE),
                        phi_zero = sum(f$DISP_APP %in% 0), row.names = NULL)
}
A3 <- do.call(rbind, A3)
check_save("A3", A3, "run", headline = sprintf("%s: %d at the lower end, %d at the upper end (of %d; %d at phi = 0 by the slope rule)", A3$dataset, A3$lower_end, A3$upper_end, A3$genes, A3$phi_zero),
           mc_error = "exact counts", script = CHECK.SCRIPT)

##############################################################################
## A5. Variance ratio of the ratio mean to the NB maximum likelihood mean, against the cap 1 + CV^2(O) (Methods: Estimating mu)
A5.GROUPS <- c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE", "HYB.COMB")
A5 <- list()
for (k in A5.GROUPS) {
  e <- CONTRAST.EXPOS[[k]]; cap <- 1 + mean((e - mean(e))^2) / mean(e)^2           # population CV of the offsets, as in the cap's derivation
  r <- vapply(rownames(SPLIT.FITS[[k]]), var_ratio_gene, numeric(1), expo = e, fit = SPLIT.FITS[[k]])
  A5[[k]] <- data.frame(dataset = k, cv_offset = sqrt(cap - 1), cap = cap, genes = sum(is.finite(r)), q50 = median(r, na.rm = TRUE), q90 = quantile(r, 0.9, na.rm = TRUE),
                        max = max(r, na.rm = TRUE), max_over_cap = max(r / cap, na.rm = TRUE), below_1 = sum(r < 1 - 1e-12, na.rm = TRUE), above_cap = sum(r > cap * (1 + 1e-12), na.rm = TRUE), row.names = NULL)
}
A5 <- do.call(rbind, A5)
check_save("A5", A5, "run", headline = sprintf("%s: CV(O) %.2f, cap %.3f, ratio median %.3f, max %.3f (max/cap %.3f); below 1: %d, above cap: %d", A5$dataset, A5$cv_offset, A5$cap,
                                               A5$q50, A5$max, A5$max_over_cap, A5$below_1, A5$above_cap), mc_error = "exact", script = CHECK.SCRIPT)

##############################################################################
## A7 (data arm). Offset adequacy by depth tercile, and the largest library share of each gene
A7.GROUPS <- c(MIX.SC = "MIX.SC", MIX.SE = "MIX.SE", HYB.SC = "HYB.SC", HYB.SE = "HYB.SE")
A7 <- list()
for (k in names(A7.GROUPS)) {
  m <- if (k == "HYB.SC") HYB.SC else if (k == "HYB.SE") HYB.SE else if (have1) get(k, E1.ENV) else NULL
  if (is.null(m)) next
  e <- CONTRAST.EXPOS[[k]]
  stopifnot(ncol(m) == length(e))
  terc <- cut(rank(e, ties.method = "first"), 3, labels = FALSE)
  ph <- t(vapply(SUB.GENES, tercile_phi_gene, numeric(6), mat = m, expo = e, terc = terc))
  lib <- if (k %in% c("HYB.SC", "HYB.SE")) SHARE.CELLS$lib else colSums(m)
  share <- apply(sweep(as.matrix(m[SUB.GENES, , drop = FALSE]), 2, lib, "/"), 1, max)
  A7[[k]] <- data.frame(dataset = k, tercile = 1:3, median_phi_offset = apply(ph[, 1:3], 2, median, na.rm = TRUE), median_phi_constant = apply(ph[, 4:6], 2, median, na.rm = TRUE),
                        median_max_share = median(share), max_max_share = max(share), row.names = NULL)
}
A7 <- if (length(A7)) do.call(rbind, A7) else NULL
if (is.null(A7)) check_save("A7_data", NULL, "not run", headline = "count matrices missing", script = CHECK.SCRIPT) else
  check_save("A7_data", A7, "run", headline = c(sprintf("%s tercile %d: median phi with offset %.3f, constant offset %.3f", A7$dataset, A7$tercile, A7$median_phi_offset, A7$median_phi_constant),
                                                 sprintf("largest library share of a gene in a cell: median over genes %.2g, max %.2g", median(A7$median_max_share), max(A7$max_max_share))),
             mc_error = sprintf("%d genes per group; tercile medians", length(SUB.GENES)), script = CHECK.SCRIPT)

##############################################################################
## A8. Binomial thinning invariance of the parental fits (UMI capture): keep each UMI with probability 0.5, recompute libraries and offsets, refit
if (have1) {
  # Thinning halves every library and the reference depth D with it, so the offsets are unchanged while the counts halve: mu-hat halves and phi-hat stays.
  A8 <- list(); set.seed(SEED.CHECKS + 8)
  D0 <- median(c(colSums(get("MIX.SC", E1.ENV)), colSums(get("MIX.SE", E1.ENV))))
  TH <- lapply(c(MIX.SC = "MIX.SC", MIX.SE = "MIX.SE"), function(k) { m <- as.matrix(get(k, E1.ENV)); matrix(rbinom(length(m), m, 0.5), nrow(m), dimnames = dimnames(m)) })
  D1 <- median(c(colSums(TH$MIX.SC), colSums(TH$MIX.SE)))
  for (k in c("MIX.SC", "MIX.SE")) {
    m <- as.matrix(get(k, E1.ENV)); e0 <- colSums(m) / D0; e1 <- colSums(TH[[k]]) / D1
    f0 <- t(vapply(SUB.GENES, function(g) fit_dispersion(m[g, ], e0)[c("mu", "disp")], numeric(2)))
    f1 <- t(vapply(SUB.GENES, function(g) fit_dispersion(TH[[k]][g, ], e1)[c("mu", "disp")], numeric(2)))
    strat <- cut(log10(f0[, 1] * mean(e0)), c(-Inf, -0.5, 0, 0.5, Inf), labels = c("<0.3", "0.3-1", "1-3", ">3"))
    for (s in levels(strat)) {
      ix <- strat == s & f0[, 2] > 0
      A8[[paste(k, s)]] <- data.frame(dataset = k, mean_count_stratum = s, genes = sum(ix), median_mu_ratio = median(f1[ix, 1] / f0[ix, 1]), median_phi_ratio = median(f1[ix, 2] / f0[ix, 2]),
                                      reference_depth_ratio = D1 / D0, row.names = NULL)
    }
  }
  A8 <- do.call(rbind, A8)
  check_save("A8", A8, "run", headline = sprintf("%s mean count %s: mu ratio %.3f (expect 0.5), phi ratio %.3f (expect 1), %d genes", A8$dataset, A8$mean_count_stratum, A8$median_mu_ratio, A8$median_phi_ratio, A8$genes),
             mc_error = "one thinning draw; medians over genes", script = CHECK.SCRIPT)
} else check_save("A8", NULL, "not run", headline = "the Section 1 checkpoint (parental counts) is missing", script = CHECK.SCRIPT)

##############################################################################
## A11. Third moment of the expression factor against the gamma (Methods: Separate partitions; Notes for revision: non-gamma expression factors)
# The mean and dispersion errors of one group are uncorrelated on shared cells only when the factor's third central moment equals the gamma's 2 phi^2;
# the excess E[(a - 1)^3] - 2 phi^2 sets their covariance and the bias of the NB fit under a non-gamma factor (B4). Factorial moments estimate it for
# every gene resolved in the dataset (factor_moment_table()). Per-gene third moments are noisy and skewed at these cell numbers, so the subsample genes are tested
# against NB data simulated at their own fitted mean and variance on the same offsets (factor_excess_test(), A11.SIM draws), which carries the
# estimator's finite-sample bias into the null; the share of genes with p < 0.05 is read against 0.05, and the observed k3 / (2 phi^2) against its
# null mean. One gene's test has little power (simulated lognormal factors with phi 0.5 and 2000 cells: p < 0.05 in about 13% of genes), so the group
# statistic is the mean over the tested genes of z = (excess - null mean) / null SD, with SE 1 / sqrt(genes), which is near 0 for gamma factors and
# positive (more skewed than gamma) or negative (less skewed) otherwise. The hybrid alleles' factor includes the transcriptome share multiplier.
A11.GROUPS <- c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE"); A11.SIM <- 200
A11 <- list(table = list(), tests = list(), summary = list())
for (k in A11.GROUPS) {
  m <- if (k == "HYB.SC") HYB.SC else if (k == "HYB.SE") HYB.SE else if (have1) get(k, E1.ENV) else NULL
  if (is.null(m)) next
  e <- CONTRAST.EXPOS[[k]]; cls1 <- GENE.TABLE$gene[GENE.TABLE[[paste0(k, ".resolved")]] %in% TRUE]     # genes whose dispersion is resolved in this dataset
  tab <- factor_moment_table(m, e, genes = cls1)
  tg <- intersect(SUB.GENES, cls1)
  tst <- with_local_seed(SEED.CHECKS + 9, t(vapply(tg, factor_excess_test, numeric(6), mat = m, expo = e, B = A11.SIM)))
  A11$table[[k]] <- data.frame(dataset = k, tab, row.names = NULL); A11$tests[[k]] <- data.frame(dataset = k, gene = tg, tst, row.names = NULL)
  ok <- is.finite(tst[, "p"])
  A11$summary[[k]] <- data.frame(dataset = k, genes = nrow(tab), pooled_k3_over_gamma = sum(tab$k3, na.rm = TRUE) / sum(tab$k3_gamma[is.finite(tab$k3)], na.rm = TRUE),
                                tested = sum(ok), median_ratio = median(tst[ok, "ratio"]), median_null_ratio = median(tst[ok, "null_ratio"]),
                                mean_z = mean(tst[ok, "z"]), mean_z_se = 1 / sqrt(sum(ok)), share_p05 = mean(tst[ok, "p"] < 0.05), share_p05_se = binom_se(0.05, sum(ok)),
                                share_above = mean(tst[ok, "p"] < 0.05 & tst[ok, "excess"] > tst[ok, "null_excess"]),
                                share_below = mean(tst[ok, "p"] < 0.05 & tst[ok, "excess"] < tst[ok, "null_excess"]), row.names = NULL)
}
A11 <- lapply(A11, function(x) do.call(rbind, x))
if (is.null(A11$summary)) check_save("A11", NULL, "not run", headline = "count matrices missing", script = CHECK.SCRIPT) else
  check_save("A11", A11, "run", headline = sprintf("%s: k3 / (2 phi^2) pooled over %d resolved genes %.2f; tested %d genes: mean z %.2f (SE %.2f), median ratio %.2f against %.2f under a gamma factor, p < 0.05 in %.3f (above %.3f, below %.3f)",
                                                 A11$summary$dataset, A11$summary$genes, A11$summary$pooled_k3_over_gamma, A11$summary$tested, A11$summary$mean_z, A11$summary$mean_z_se, A11$summary$median_ratio,
                                                 A11$summary$median_null_ratio, A11$summary$share_p05, A11$summary$share_above, A11$summary$share_below),
             mc_error = sprintf("%d null draws per gene; binomial SE of a 5%% share %.3f", A11.SIM, A11$summary$share_p05_se[1]), script = CHECK.SCRIPT)

##############################################################################
## B1. Binomial sampling and independence of the observed share (Methods: Transcriptome share estimate)
# Within a library bin Var(c-hat) = sigma2_c + E[c (1 - c) / L], so regressing the bin variance on 1 / L-bar gives slope c-bar (1 - c-bar) and intercept sigma2_c.
B1 <- local({
  L <- SHARE.CELLS$lib; bin <- cut(rank(L, ties.method = "first"), 10, labels = FALSE)
  tab <- data.frame(bin = 1:10, inv_L = tapply(1 / L, bin, mean), var_chat = tapply(C.HAT, bin, var), n = as.vector(table(bin)))
  fit <- lm(var_chat ~ inv_L, data = tab)
  cc <- if (file.exists(ckpt_path(7))) { e7 <- new.env(); load(ckpt_path(7), envir = e7); if (exists("CC.AXIS", envir = e7)) get("CC.AXIS", e7) else NULL } else NULL
  list(table = tab, coef = summary(fit)$coefficients, cbar = C.BAR, expected_slope = C.BAR * (1 - C.BAR), var_c = VAR.C,
       spearman_chat_L = cor.test(C.HAT, L, method = "spearman", exact = FALSE)[c("estimate", "p.value")], cell_cycle = cc)
})
check_save("B1", B1, "run",
           headline = c(sprintf("slope %.4f (SE %.4f) against c-bar(1 - c-bar) = %.4f; intercept %.3g (SE %.2g) against sigma2_c-hat %.3g", B1$coef[2, 1], B1$coef[2, 2], B1$expected_slope,
                                B1$coef[1, 1], B1$coef[1, 2], B1$var_c),
                        sprintf("Spearman(c-hat, L) %.3f (p %.2g)", B1$spearman_chat_L$estimate, B1$spearman_chat_L$p.value),
                        if (is.null(B1$cell_cycle)) "cell-cycle correlation not run: CC.AXIS (Section 7 checkpoint) is missing" else "cell-cycle axis loaded; see rds"),
           mc_error = "regression SEs", script = CHECK.SCRIPT)

##############################################################################
## B2 (data arm). Interval of eta-hat and the split-half estimate of sigma2_c
# Splitting each cell's UMIs at random into two halves gives two shares that, given c_i, are independent binomial estimates of it, so their covariance
# across cells estimates sigma2_c without the binomial correction.
B2D <- with_local_seed(SEED.CHECKS + 2, local({
  L <- SHARE.CELLS$lib; Lc <- round(SHARE.CELLS$chat * L); n <- length(L)          # the hybrid libraries and their Sc UMIs, as Section 1.2 measured them
  bs <- t(vapply(seq_len(2000), function(b) { ix <- sample.int(n, n, TRUE); st <- share_terms(C.HAT[ix], L[ix]); c(eta_c = st$var_obs - st$var_binom, cbar = st$cbar) }, numeric(2)))
  eta_c <- bs[, 1] / bs[, 2]^2; eta_e <- bs[, 1] / (1 - bs[, 2])^2
  half <- numeric(200)
  for (b in seq_along(half)) {
    c1 <- rbinom(n, Lc, 0.5); e1 <- rbinom(n, L - Lc, 0.5); ok <- c1 + e1 > 0 & (L - c1 - e1) > 0
    half[b] <- cov((c1 / (c1 + e1))[ok], ((Lc - c1) / (L - c1 - e1))[ok])
  }
  list(eta_Hc = R.C, eta_He = R.E, eta_Hc_ci = quantile(eta_c, c(0.025, 0.975)), eta_He_ci = quantile(eta_e, c(0.025, 0.975)), var_c_ltv = VAR.C,
       var_c_split_half = mean(half), var_c_split_half_mc_se = sd(half) / sqrt(length(half)), var_c_split_half_sd_one_split = sd(half))
}))
check_save("B2_data", B2D, "run",
           headline = c(sprintf("eta_Hc %.3g (95%% bootstrap interval %.3g to %.3g), eta_He %.3g (%.3g to %.3g)", B2D$eta_Hc, B2D$eta_Hc_ci[1], B2D$eta_Hc_ci[2], B2D$eta_He, B2D$eta_He_ci[1], B2D$eta_He_ci[2]),
                        sprintf("sigma2_c: law of total variance %.3g, split-half covariance %.3g (SD of one split %.2g)", B2D$var_c_ltv, B2D$var_c_split_half, B2D$var_c_split_half_sd_one_split)),
           mc_error = sprintf("2000 bootstrap replicates; split-half mean of 200 splits, MC SE %.2g", B2D$var_c_split_half_mc_se), script = CHECK.SCRIPT)

##############################################################################
## B5. Per-genome hybrid mean from the summed rescaled counts against the average of the corrected allele rates (Methods: Per-genome hybrid estimates)
B5 <- local({
  z <- HYB.SC[GENES, , drop = FALSE] / (2 * C.BAR) + HYB.SE[GENES, , drop = FALSE] / (2 * (1 - C.BAR))
  summed <- rowSums(z) / sum(EXPO.HYB)
  avg <- 0.5 * (SPLIT.FITS$HYB.SC[GENES, "MU"] / (2 * C.BAR) + SPLIT.FITS$HYB.SE[GENES, "MU"] / (2 * (1 - C.BAR)))
  d <- log2(summed) - log2(avg)
  data.frame(gene = GENES, log2_summed = log2(summed), log2_allele_average = log2(avg), diff = d)
})
check_save("B5", B5, "run", headline = sprintf("%d genes: max |log2 difference| %.2g, mean difference %.2g (Bland-Altman limits %.2g to %.2g)", nrow(B5), max(abs(B5$diff), na.rm = TRUE),
                                               mean(B5$diff, na.rm = TRUE), mean(B5$diff, na.rm = TRUE) - 1.96 * sd(B5$diff, na.rm = TRUE), mean(B5$diff, na.rm = TRUE) + 1.96 * sd(B5$diff, na.rm = TRUE)),
           mc_error = "exact", script = CHECK.SCRIPT)

##############################################################################
## B6. Genes whose fit the share correction truncates at zero
B6 <- list()
for (k in FIT.GROUPS) {
  f <- SPLIT.FITS[[k]]
  B6[[k]] <- data.frame(dataset = k, truncated_by_correction = sum(f$DISP %in% 0 & f$DISP_APP > 0, na.rm = TRUE), row.names = NULL)
}
B6 <- do.call(rbind, B6)
check_save("B6", B6, "run", headline = sprintf("%s: truncated to zero by the correction %d", B6$dataset, B6$truncated_by_correction),
           mc_error = "exact counts", script = CHECK.SCRIPT)

##############################################################################
## C1. Partition fractions and sizes (Methods: Estimating partition fractions)
check_save("C1", list(report = PARTITION.REPORT, frac = SPLIT.FRAC), "run",
           headline = sprintf("%s: gamma %.3f, f %.3f, cis %d, trans %d cells", PARTITION.REPORT$quantity, PARTITION.REPORT$gamma, PARTITION.REPORT$f, PARTITION.REPORT$n_cis, PARTITION.REPORT$n_trans),
           mc_error = "pilot bootstrap; see PILOT.LOGPHI.SE", script = CHECK.SCRIPT)

##############################################################################
## C2. Subsampling: the cis contrast variance falls as 1 / n (slope -1 of log variance on log n)
C2 <- with_local_seed(SEED.CHECKS + 12, local({
  n <- ncol(HYB.SC); fr <- c(0.25, 0.5, 0.75, 1)
  subs <- lapply(fr, function(f) sort(sample.int(n, round(f * n)))); boot <- lapply(subs, function(s) lapply(seq_len(100), function(b) sample.int(length(s), length(s), TRUE)))
  V <- mclapply(SUB.GENES[seq_len(min(100, length(SUB.GENES)))], subsample_cis_var_gene, yc = HYB.SC, ye = HYB.SE, expo = EXPO.HYB.SC, eta = SHARE.R, subs = subs, boot = boot, mc.cores = CORES)
  slopes <- matrix(NA_real_, length(V), 2, dimnames = list(NULL, c("mean", "disp")))
  for (j in seq_along(V)) for (q in c("mean", "disp")) {
    v <- V[[j]][, paste0("var_", q)]
    if (all(is.finite(v) & v > 0)) slopes[j, q] <- unname(coef(lm(log(v) ~ log(V[[j]][, "n"])))[2])
  }
  list(variances = V, slopes = slopes)
}))
check_save("C2", C2, "run", headline = sprintf("slope of log variance on log n: mean median %.3f (IQR %.3f to %.3f), dispersion median %.3f (IQR %.3f to %.3f)",
                                               median(C2$slopes[, "mean"], na.rm = TRUE), quantile(C2$slopes[, "mean"], 0.25, na.rm = TRUE), quantile(C2$slopes[, "mean"], 0.75, na.rm = TRUE),
                                               median(C2$slopes[, "disp"], na.rm = TRUE), quantile(C2$slopes[, "disp"], 0.25, na.rm = TRUE), quantile(C2$slopes[, "disp"], 0.75, na.rm = TRUE)),
           mc_error = "100 bootstrap replicates per fraction, 100 genes", script = CHECK.SCRIPT)

##############################################################################
## C3. Equal precision of cis and trans at the pooled f: per-gene var_cis / var_trans from the bootstrap SEs, and gene-specific against pooled gamma
C3 <- list()
for (q in c("mean", "disp", "kbal", "bsize")) {
  r <- (BOOT.CONTRASTS[[paste0(q, "_cis_se")]] / BOOT.CONTRASTS[[paste0(q, "_trans_se")]])^2; r <- r[is.finite(r) & r > 0]
  C3[[q]] <- data.frame(quantity = q, genes = length(r), median_ratio = median(r), q25 = quantile(r, 0.25), q75 = quantile(r, 0.75), median_log2_ratio = median(log2(r)), row.names = NULL)
}
C3 <- do.call(rbind, C3)
check_save("C3", C3, "run", headline = sprintf("%s: var_cis/var_trans median %.3f (IQR %.3f to %.3f), %d genes", C3$quantity, C3$median_ratio, C3$q25, C3$q75, C3$genes),
           mc_error = sprintf("from BOOT.CONTRASTS (%s bootstrap SEs)", "gene_boot"), script = CHECK.SCRIPT)

##############################################################################
## C6. Depth balance of the cis and trans subsets, and stability of calls over 50 random partitions
C6 <- local({
  ks <- data.frame(partition = c("mean", "dispersion", "kinetic balance", "burst size"),
                   ks_p = vapply(list(SPLIT.IDX.MEAN, SPLIT.IDX.DISP, SPLIT.IDX.KBAL, SPLIT.IDX.BSIZE), function(s) suppressWarnings(ks.test(EXPO.HYB[s$c], EXPO.HYB[s$t])$p.value), numeric(1)))
  g <- intersect(GENES, BOOT.CONTRASTS$gene); bc <- BOOT.CONTRASTS[match(g, BOOT.CONTRASTS$gene), ]
  S <- check_streams(50, SEED.CHECKS + 16)
  calls <- mclapply(seq_len(50), split_mean_calls, genes = g, yc = HYB.SC, ye = HYB.SE, expo = EXPO.HYB.SC, depth = EXPO.HYB, f = SPLIT.FRAC$f_mean,
                    mu_pc = SPLIT.FITS$MIX.SC[g, "MU"], mu_pe = SPLIT.FITS$MIX.SE[g, "MU"], cbar = C.BAR, se_cis = bc$mean_cis_se, se_trans = bc$mean_trans_se, streams = S, mc.cores = CORES)
  mc <- sapply(calls, function(x) x[, "cis"]); mt <- sapply(calls, function(x) x[, "trans"])
  list(ks = ks, changed_cis = mean(rowSums(mc) > 0 & rowSums(mc) < ncol(mc)), changed_trans = mean(rowSums(mt) > 0 & rowSums(mt) < ncol(mt)), genes = length(g))
})
check_save("C6", C6, "run", headline = c(sprintf("KS p (depth, cis vs trans) %s: %.3g", C6$ks$partition, C6$ks$ks_p),
                                         sprintf("genes whose cis call changes over 50 random partitions %.3f, trans %.3f (Wald approximation, %d genes)", C6$changed_cis, C6$changed_trans, C6$genes)),
           mc_error = "50 partitions", script = CHECK.SCRIPT)

##############################################################################
## D9. Genes at the p-value floor 1 / (B + 1), and whether BH at q = 0.05 can reach it (Methods: p-values)
# B is N.PERM for the permutation nulls (mean contrasts except mid) and N.SIM for the simulated nulls. With m tests, k genes at the floor are all rejected
# by BH when floor <= 0.05 k / m.
D9 <- list()
NSIM.D9 <- local({ f <- file.path(INPUT.DIR, "gene_perm_inputs.rda"); e <- new.env(); if (file.exists(f)) { load(f, envir = e); get("N.SIM", e) } else 2000 })
for (cn in grep("_p$", names(PR), value = TRUE)) {
  p <- PR[[cn]]; p <- p[is.finite(p)]
  if (!length(p)) next
  B <- if (grepl("^mean_", cn) && cn != "mean_mid_p") N.PERM else NSIM.D9; fl <- 1 / (B + 1); k <- sum(p <= fl * (1 + 1e-9))
  D9[[cn]] <- data.frame(test = cn, m = length(p), B = B, floor = fl, at_floor = k, min_p = min(p), bh_reaches_floor = k > 0 && fl <= 0.05 * k / length(p),
                         bh_one_gene = fl <= 0.05 / length(p), row.names = NULL)
}
D9 <- do.call(rbind, D9)
check_save("D9", D9, "run", headline = sprintf("%s: %d of %d at the floor %.2g; BH reaches the floor: %s", D9$test, D9$at_floor, D9$m, D9$floor, D9$bh_reaches_floor),
           mc_error = "exact counts", script = CHECK.SCRIPT)

##############################################################################
## E2. Tested and significant dispersion contrasts across resolution thresholds (Methods: Informative estimates)
# A group is resolved when at most t of its bootstrap replicates refit to zero (RESOLVED.ZERO.MAX = 0.025 in the analysis). The zero fractions are stored per gene
# and dataset in GENE.TABLE, and the unmasked null p-values in PERM.RESULTS, so every threshold is recomputed exactly: a contrast is tested when its fits define it
# and every group is resolved at t, and significant at BH q < 0.05 over the genes tested at t.
E2 <- list(counts = list(), significant = list())
E2.T <- c(0.01, 0.025, 0.05, 0.10)
PR.ALL <- PERM.RESULTS[match(GENE.TABLE$gene, PERM.RESULTS$gene), ]
for (t in E2.T) {
  res <- sapply(.CONTRAST.GROUPS, function(k) (GENE.TABLE[[paste0(k, ".zero_frac")]] <= t) %in% TRUE & (GENE.TABLE[[paste0(k, ".disp")]] > 0) %in% TRUE)
  E2$counts[[as.character(t)]] <- data.frame(zero_frac_max = t, dataset = .CONTRAST.GROUPS, resolved = colSums(res), row.names = NULL)
  for (md in .OUT_MODES) {
    gr <- .MODES[[.DISP_SOURCE[[md]]]]; ok <- rowSums(res[, gr, drop = FALSE]) == length(gr)
    pl <- PR.ALL[[paste0("disp_", md, "_p")]]; pr <- PR.ALL[[paste0("rdisp_", md, "_p")]]
    tl <- ok & GENE.TABLE[[paste0("log_", md, "_valid")]] & is.finite(pl); tr <- ok & GENE.TABLE[[paste0("raw_", md, "_valid")]] & is.finite(pr)
    E2$significant[[paste(t, md)]] <- data.frame(zero_frac_max = t, mode = md, log2_tested = sum(tl), log2_sig = sum(p.adjust(pl[tl], "BH") < 0.05),
                                                 raw_tested = sum(tr), raw_sig = sum(p.adjust(pr[tr], "BH") < 0.05), row.names = NULL)
  }
}
E2 <- lapply(E2, function(x) do.call(rbind, x))
x <- aggregate(cbind(log2_tested, log2_sig, raw_sig) ~ zero_frac_max, E2$significant, sum)
check_save("E2", E2, "run", headline = sprintf("zero fraction at most %g: tested log2 dispersion contrasts %d, significant log2 %d, raw %d (all modes)", x$zero_frac_max, x$log2_tested, x$log2_sig, x$raw_sig),
           mc_error = "exact counts from the stored zero fractions and null p-values", script = CHECK.SCRIPT)

##############################################################################
## E3 (data arm) and E4. Defined contrasts left untested because a group is unresolved, likelihood-bound calls for one resolved group; gene sets and exclusions
## per contrast and scale
E3D <- list(unresolved = EXCLUDED.COUNTS[EXCLUDED.COUNTS$scale != "mean", c("contrast", "scale", "eligible", "tested", "unresolved")], bounds = BOUNDS.SEPARATION)
check_save("E3_data", E3D, "run", headline = c(sprintf("%s %s: %d defined, %d tested, %d unresolved", E3D$unresolved$contrast, E3D$unresolved$scale, E3D$unresolved$eligible,
                                                       E3D$unresolved$tested, E3D$unresolved$unresolved),
                                               sprintf("%s: one group resolved in %d genes, likelihood bounds separate in %d", BOUNDS.SEPARATION$contrast,
                                                       BOUNDS.SEPARATION$one_group_resolved, BOUNDS.SEPARATION$bounds_separate)),
           mc_error = "from gene_boot", script = CHECK.SCRIPT)
write.csv(EXCLUDED.COUNTS, file.path(CHECK.DIR, "E4_excluded_counts.csv"), row.names = FALSE)
write.csv(GENE.SET.SIZES, file.path(CHECK.DIR, "E4_gene_set_sizes.csv"), row.names = FALSE)
write.csv(GENE.TABLE[, c("gene", grep("_valid$|^resolved_|^in_matched_set$|^inference_set$", names(GENE.TABLE), value = TRUE))], file.path(CHECK.DIR, "E4_gene_sets.csv"), row.names = FALSE)
check_save("E4", list(excluded = EXCLUDED.COUNTS, sizes = GENE.SET.SIZES), "code", headline = sprintf("tables written: E4_excluded_counts.csv, E4_gene_set_sizes.csv, E4_gene_sets.csv; matched set %d genes",
                                                                                                     sum(GENE.TABLE$in_matched_set)), script = CHECK.SCRIPT)

##############################################################################
## F1 (b, c). The Appendix A estimator on the data: convergence, boundary share, and agreement with the Cox-Reid estimate
F1D <- list()
for (k in FIT.GROUPS) {
  w <- SPLIT.FITS.W[[k]]; f <- SPLIT.FITS[[k]]; ok <- w$DISP > 0 & f$DISP > 0 & is.finite(w$DISP) & is.finite(f$DISP); d <- log2(w$DISP[ok]) - log2(f$DISP[ok])
  F1D[[k]] <- data.frame(dataset = k, converged = mean(w$CONVERGED %in% 1), boundary = mean(w$BOUNDARY %in% 1), genes_both_positive = sum(ok), mean_log2_diff = mean(d),
                         loa_low = mean(d) - 1.96 * sd(d), loa_high = mean(d) + 1.96 * sd(d), row.names = NULL)
}
F1D <- do.call(rbind, F1D)
# The weighted estimator of the code is not the Appendix A estimator as documented (mismatches.md); F1_sim reports both.
check_save("F1_data", F1D, "run", headline = sprintf("%s: converged %.4f, boundary %.3f; log2(weighted / Cox-Reid) mean %.3f (limits %.2f to %.2f)", F1D$dataset, F1D$converged, F1D$boundary,
                                                     F1D$mean_log2_diff, F1D$loa_low, F1D$loa_high),
           mc_error = "exact", script = CHECK.SCRIPT)

##############################################################################
## F5. Predicted shift delta-hat against the observed log2(phi_H / phi_Sigma) (Appendix B)
F5 <- local({
  phiH <- sqrt(SPLIT.FITS$HYB.SC[GENES, "DISP"] * SPLIT.FITS$HYB.SE[GENES, "DISP"]); obs <- log2(phiH / SPLIT.FITS$HYB.COMB[GENES, "DISP"]); pred <- PLOIDY.COMP[GENES, "delta"]
  ok <- is.finite(obs) & is.finite(pred); fit <- lm(obs[ok] ~ pred[ok])
  list(coef = summary(fit)$coefficients, genes = sum(ok), spearman = cor(obs[ok], pred[ok], method = "spearman"), route_check = ROUTE.CHECK)
})
check_save("F5", F5, "run", headline = sprintf("observed on predicted: slope %.3f (SE %.3f), intercept %.3f (SE %.3f), Spearman %.3f, %d genes", F5$coef[2, 1], F5$coef[2, 2], F5$coef[1, 1], F5$coef[1, 2], F5$spearman, F5$genes),
           mc_error = "regression SEs", script = CHECK.SCRIPT)

##############################################################################
## G1. Every number the text quotes, regenerated from the pipeline
G1 <- local({
  med_c <- median(GENE.TABLE$HYB.SC.disp[GENE.TABLE$HYB.SC.resolved]); med_e <- median(GENE.TABLE$HYB.SE.disp[GENE.TABLE$HYB.SE.resolved])
  med_c_all <- median(SPLIT.FITS$HYB.SC$DISP, na.rm = TRUE); med_e_all <- median(SPLIT.FITS$HYB.SE$DISP, na.rm = TRUE)
  nmin <- min(c(vapply(CONTRAST.EXPOS[c("MIX.SC", "MIX.SE", "HYB.SC")], length, integer(1)), PARTITION.REPORT$n_cis, PARTITION.REPORT$n_trans))
  data.frame(quantity = c("c-bar-hat", "eta_Hc", "eta_He", "median phi Hc (resolved)", "median phi He (resolved)", "median phi Hc (all fits)", "median phi He (all fits)",
                          "correction term Hc (1 + phi) eta", "correction term He (1 + phi) eta", "phi below which (1 + phi) eta / phi > 1% (Hc)", "smallest group or partition (cells)"),
             text = c(0.5024, 7.4e-5, 7.5e-5, 0.150, 0.152, 0.150, 0.152, 8.5e-5, 8.6e-5, 0.0075, 800),
             pipeline = c(C.BAR, R.C, R.E, med_c, med_e, med_c_all, med_e_all, (1 + med_c) * R.C, (1 + med_e) * R.E, R.C / (0.01 - R.C), nmin), row.names = NULL)
})
G1$relative_difference <- G1$pipeline / G1$text - 1
write.csv(G1, file.path(CHECK.DIR, "G1_numbers.csv"), row.names = FALSE)
check_save("G1", G1, "run", headline = sprintf("%s: text %g, pipeline %.4g", G1$quantity, G1$text, G1$pipeline), mc_error = "exact", script = CHECK.SCRIPT)

##############################################################################
## Local bootstrap of a gene subsample: A4 (data), A6, C7 (data) and D8
# The indices are drawn first for every replicate (make_draws(), as in gene_boot: the undivided hybrid cells resampled once, each copy keeping its partition
# memberships), so the hybrid partitions share their resampled cells.
# D8 refits each replicate twice from the same draws: with the share terms re-estimated on the resampled cells (the pipeline) and with the observed ones.
NCELLS.CK <- vapply(CONTRAST.EXPOS[names(.DRAW.KEY)], length, integer(1))
DRAWS.CK <- make_draws(NCELLS.CK, N.CHECK.BOOT, SEED.CHECKS + 30, splits = list(HYC = SPLIT.IDX.MEAN$c, HYC.N = SPLIT.IDX.DISP$c, HYC.K = SPLIT.IDX.KBAL$c, HYC.Z = SPLIT.IDX.BSIZE$c))
CK.GROUPS <- c(.CONTRAST.GROUPS, "HYB.COMB")
MATS.CK <- lapply(setNames(CK.GROUPS, CK.GROUPS), function(k) as.matrix(CONTRAST.MATS[[k]][SUB.GENES, , drop = FALSE]))
BOOT.CK <- mclapply(SUB.GENES, data_boot_gene, mats = MATS.CK, expos = CONTRAST.EXPOS, eta = SHARE.R, draws = DRAWS.CK, refit_share = TRUE, share_cells = SHARE.CELLS,
                    groups = CK.GROUPS, share_sub = setNames(PLOIDY.COMP$share_sub, rownames(PLOIDY.COMP)), mc.cores = CORES)
BOOT.CK.FIXED <- mclapply(SUB.GENES, data_boot_gene, mats = MATS.CK, expos = CONTRAST.EXPOS, eta = SHARE.R, draws = DRAWS.CK, refit_share = FALSE, mc.cores = CORES)
STATS.CK <- do.call(rbind, mclapply(BOOT.CK, boot_gene_stats, cbar = C.BAR, mc.cores = CORES)); rownames(STATS.CK) <- SUB.GENES
STATS.CK.FIXED <- do.call(rbind, mclapply(BOOT.CK.FIXED, boot_gene_stats, cbar = C.BAR, mc.cores = CORES)); rownames(STATS.CK.FIXED) <- SUB.GENES
A4D <- list(within = col_quartiles(STATS.CK, grep("^cor_mu_phi_(raw|log2)\\.", colnames(STATS.CK))), contrasts = col_quartiles(STATS.CK, grep("^cor_dmu_dphi_", colnames(STATS.CK))))
check_save("A4_data", A4D, "run",
           headline = c(sprintf("within-group cor(mu, phi) raw: median over datasets of gene medians %.3f (range %.3f to %.3f)", median(A4D$within[grepl("raw", rownames(A4D$within)), 2]),
                                min(A4D$within[grepl("raw", rownames(A4D$within)), 2]), max(A4D$within[grepl("raw", rownames(A4D$within)), 2])),
                        sprintf("%s: median %.3f (IQR %.3f to %.3f)", sub("^cor_dmu_dphi_", "", rownames(A4D$contrasts)), A4D$contrasts[, 2], A4D$contrasts[, 1], A4D$contrasts[, 3])),
           mc_error = sprintf("%d genes, %d replicates; SE of one gene's correlation about %.3f", length(SUB.GENES), N.CHECK.BOOT, 1 / sqrt(N.CHECK.BOOT - 3)), script = CHECK.SCRIPT)
A6 <- col_quartiles(STATS.CK, grep("^(cor_zeta_theta|cor_kappa_mu|se_ratio_kappa_theta|se_ratio_logphi_logmu)\\.(MIX|HYB)\\.", colnames(STATS.CK)))
check_save("A6", A6, "run", headline = sprintf("%s: median %.3f (IQR %.3f to %.3f)", rownames(A6), A6[, 2], A6[, 1], A6[, 3]),
           mc_error = sprintf("%d genes, %d replicates", length(SUB.GENES), N.CHECK.BOOT), script = CHECK.SCRIPT)
C7D <- col_quartiles(STATS.CK, c("cor_mu_Hc_phi_He_kbal", "cor_mu_He_phi_Hc_kbal", "cor_mu_Hc_phi_He_bsize", "cor_mu_He_phi_Hc_bsize"))
check_save("C7_data", C7D, "run", headline = sprintf("%s: median %.3f (IQR %.3f to %.3f)", rownames(C7D), C7D[, 2], C7D[, 1], C7D[, 3]),
           mc_error = sprintf("SE of one gene's correlation about %.3f", 1 / sqrt(N.CHECK.BOOT - 3)), script = CHECK.SCRIPT)
D8 <- local({
  se1 <- STATS.CK[, grep("^se_disp\\.", colnames(STATS.CK))]; se0 <- STATS.CK.FIXED[, grep("^se_disp\\.", colnames(STATS.CK.FIXED))]
  mbar <- SPLIT.FITS$HYB.SC[SUB.GENES, "MBAR"]; strat <- cut(mbar, quantile(mbar, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE, labels = c("low", "mid", "high"))
  r <- se1 / se0
  data.frame(stratum = rep(levels(strat), each = ncol(r)), mode = sub("^se_disp\\.", "", colnames(r)),
             median_se_ratio = as.vector(vapply(levels(strat), function(s) apply(r[strat %in% s, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(r)))))
})
check_save("D8", D8, "run", headline = sprintf("%s mean count %s: SE ratio (refit / fixed share) %.4f", D8$mode, D8$stratum, D8$median_se_ratio),
           mc_error = sprintf("same draws in both arms; %d genes, %d replicates", length(SUB.GENES), N.CHECK.BOOT), script = CHECK.SCRIPT)

##############################################################################
## F5 and B5, the two routes to the per-genome hybrid dispersion on the local bootstrap (Cox-Reid fits)
# Route A, used by every dominance contrast: log2 of the geometric mean of the two corrected allele dispersions. The summed route: log2 of the corrected
# summed-count dispersion plus the reported shift delta, which is fixed per gene, so its bootstrap SD is that of log2 phi_Sigma. Both come from the same
# resampled hybrid cells. Reported: the median over genes of SD(summed route) / SD(route A), and the Spearman correlation of the two point estimates
# (ROUTE.CHECK, all genes of the matched set). The earlier moment-based comparison (compare_phiG_routes.R) gave 0.965 and 0.97.
ROUTE.SD <- t(vapply(BOOT.CK, function(D) c(allele = sd(0.5 * (log2(D[, "phi.HYB.SC"]) + log2(D[, "phi.HYB.SE"]))[is.finite(log2(D[, "phi.HYB.SC"]) + log2(D[, "phi.HYB.SE"]))]),
                                             summed = sd(log2(D[, "phi.HYB.COMB"])[is.finite(log2(D[, "phi.HYB.COMB"]))])), numeric(2)))
check_save("F5_routes", list(sd = ROUTE.SD, route_check = ROUTE.CHECK), "run",
           headline = sprintf("bootstrap SD of log2 phi_H, summed route over allele route: median %.3f (IQR %.3f to %.3f), %d genes; Spearman of the point estimates %.3f (%d genes)",
                              median(ROUTE.SD[, "summed"] / ROUTE.SD[, "allele"], na.rm = TRUE), quantile(ROUTE.SD[, "summed"] / ROUTE.SD[, "allele"], 0.25, na.rm = TRUE),
                              quantile(ROUTE.SD[, "summed"] / ROUTE.SD[, "allele"], 0.75, na.rm = TRUE), sum(is.finite(ROUTE.SD[, "summed"] / ROUTE.SD[, "allele"])), ROUTE.CHECK$spearman, ROUTE.CHECK$n_genes),
           mc_error = sprintf("%d bootstrap replicates per gene", N.CHECK.BOOT), script = CHECK.SCRIPT)
