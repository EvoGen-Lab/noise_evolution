###############################################################
### contrast_functions_check.R
### Local end-to-end check of the contrast machinery on a small synthetic design: 17 datasets (parents, the hybrid alleles, the summed hybrid and the
### mean, noise and kinetic balance splits of the hybrid cells), observed share terms, the Cox-Reid fits, and boot_contrasts_one() and
### permute_contrasts_one() with few draws. Checks the share offset on the cis and trans mean contrasts, the ploidy shift in mid and par, the
### +-Inf mapping of a zero dispersion, percentile standard errors, the split plumbing (make_draws(), perm_label_draw()), and the gate.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(11)
G <- 6; nH <- 360; nP <- 220; cbar0 <- 0.5024
gn <- paste0("g", seq_len(G)); mu <- c(2, 5, 1, 8, 3, 0.4); phi <- c(0.3, 0.5, 0.2, 0.15, 0.05, 0.6)
lib <- rlnorm(nH, log(2e4), 0.35); S <- lib / 2e4; ci <- rbeta(nH, cbar0 * 2e3, (1 - cbar0) * 2e3)
nb <- function(m, ph) matrix(rnbinom(length(m), mu = m, size = 1 / rep(ph, ncol(m))), nrow(m), dimnames = list(gn, NULL))
Hsc <- nb(outer(mu, S * ci), phi); Hse <- nb(outer(mu, S * (1 - ci)), phi)       # no cis effect: each allele keeps its genome's share of the library
eP <- rlnorm(nP, 0, 0.3); Msc <- nb(outer(mu, eP), phi); Mse <- nb(outer(mu, eP), phi)
colnames(Hsc) <- colnames(Hse) <- paste0("h", seq_len(nH))
libh <- round(lib); chat <- rbinom(nH, libh, ci) / libh              # the observed share comes from the whole library, not from these few genes
st <- share_terms(chat, libh); SR <- share_r_groups(st$r_c, st$r_e)
eH <- S / 2; Hcomb <- Hsc + Hse
pc <- ploidy_components(Hsc, Hse, eH, st)
## three depth-balanced splits at different fractions
ix <- list(MEAN = split_indices_by_depth(nH, 0.4), DISP = split_indices_by_depth(nH, 0.35), KBAL = split_indices_by_depth(nH, 0.3))
M <- list(MIX.SC = Msc, MIX.SE = Mse, HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hcomb)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = eH, HYB.SE = eH, HYB.COMB = S)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K"))
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, ix[[sp]][[side]]]; E[[nm]] <- eH[ix[[sp]][[side]]]
}
stopifnot(setequal(names(M), names(.DRAW.KEY)))
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]], share_sub = if (k == "HYB.COMB") pc$share_sub else 0)
                                       data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
CL <- lapply(names(M), function(k) fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]], share_sub = if (k == "HYB.COMB") pc$share_sub else 0)$CLASS); names(CL) <- names(M)
stopifnot(all(vapply(CL, function(x) all(x %in% 1:3), logical(1))))
share_obs <- list(cbar = st$cbar, r = SR, sub = setNames(pc$share_sub, gn), shift = setNames(pc$delta, gn), shift_raw = setNames(pc$delta_raw, gn))
share_cells <- list(chat = chat, lib = libh)
NC <- c(MIX.SC = nP, MIX.SE = nP, HYB.COMB = nH)
B <- 12
draws <- make_draws(NC, B, 1, hyc = ix$MEAN$c, hyc_n = ix$DISP$c, hyc_k = ix$KBAL$c)
## the draws: every group's indices have the group's size and index its own columns
d1 <- draws[[1]]
stopifnot(length(d1$HYB) == nH, length(d1$HYC.K) == length(ix$KBAL$c), length(d1$HYT.K) == length(ix$KBAL$t), all(d1$HYC.K <= length(ix$KBAL$c)),
          length(d1$HYC) == length(ix$MEAN$c), length(d1$HYT.N) == length(ix$DISP$t))
gate_all <- setNames(rep(TRUE, 4 * length(.OUT_MODES)), as.vector(outer(c("mean", "raw", "log", "kbal"), .OUT_MODES, paste, sep = "_")))
GATES <- matrix(TRUE, G, length(gate_all), dimnames = list(gn, names(gate_all)))
bo <- do.call(rbind, lapply(gn, boot_contrasts_one, expos = E, fits = FITS, mats = M, draws = draws, share_cells = share_cells, share_obs = share_obs, gates = GATES, methods = c("NB", "weighted")))
stopifnot(nrow(bo) == G, all(c("mean_cis_est", "disp_cis_se", "disp_cis_lo", "disp_cis_zero", "disp_cis_zeroflag", "kbal_cis_est", "kbalsup_cis_est", "wkbal_trans_est",
                               "rkbal_total_est", "cor_cis", "ploidy_shift", "ploidy_shift_raw") %in% names(bo)))
## the share offset: cis subtracts log2[c-bar / (1 - c-bar)] from the mean contrast and trans inherits it with the opposite sign
fvg <- fit_group_values(FITS, "g4")
stopifnot(abs(.contrast_value("cis", fvg$mu, "MU", st$cbar) - (.contrast_value("cis", fvg$mu, "MU", 0.5) - log2(st$cbar / (1 - st$cbar)))) < 1e-12,
          abs(.contrast_value("trans", fvg$mu, "MU", st$cbar) - (.contrast_value("trans", fvg$mu, "MU", 0.5) + log2(st$cbar / (1 - st$cbar)))) < 1e-12)
## mid and par add the ploidy shift on the noise axis and nothing on the mean axis
stopifnot(abs(.contrast_value("par_sc", fvg$disp, "DISP", shift = 0.3) - .contrast_value("par_sc", fvg$disp, "DISP") - 0.3) < 1e-12,
          abs(.contrast_value("mid", fvg$disp, "DISP", shift = 0.3) - .contrast_value("mid", fvg$disp, "DISP") - 0.3) < 1e-12,
          .contrast_value("mid", fvg$mu, "MU", st$cbar, shift = 0.3) == .contrast_value("mid", fvg$mu, "MU", st$cbar),
          abs(.contrast_raw("par_se", fvg$disp, 0.02) - .contrast_raw("par_se", fvg$disp) - 0.02) < 1e-12)
## the per-genome rescaling of the summed hybrid's mean: equal per-genome rates give a zero mid/par mean contrast at any c-bar
g_eq <- list(MIX.SC = 3, MIX.SE = 3, HYB.SC = 3 * 2 * 0.6, HYB.SE = 3 * 2 * 0.4)
stopifnot(abs(.contrast_value("mid", g_eq, "MU", 0.6)) < 1e-12, abs(.contrast_value("par_sc", g_eq, "MU", 0.6)) < 1e-12, abs(.contrast_value("cis", c(g_eq, list(HYC.SC = 6 * 0.6, HYC.SE = 6 * 0.4)), "MU", 0.6)) < 1e-12)
## a zero dispersion maps to -Inf, and percentile SEs survive it
x <- c(-Inf, rnorm(99)); s <- .boot_summary(x); stopifnot(is.finite(s[["se"]]), s[["zero"]] == 0.01)
stopifnot(is.infinite(.boot_summary(c(rep(-Inf, 10), rnorm(90)))[["lo"]]))
stopifnot(.contrast_value("total", list(MIX.SC = 0, MIX.SE = 1), "DISP", zero_inf = TRUE) == -Inf, is.nan(.contrast_value("total", list(MIX.SC = 0, MIX.SE = 0), "DISP", zero_inf = TRUE)),
          is.na(.disp_contrasts("total", list(MIX.SC = 0, MIX.SE = 0), zero_inf = TRUE)[["log"]]), is.na(.disp_contrasts("total", list(MIX.SC = 0, MIX.SE = 1))[["log"]]))
stopifnot(perm_pval(1, c(-Inf, Inf, 0.2, NA)) == (1 + 2) / (1 + 3))
## the permutations
perms <- lapply(1:15, perm_label_draw, hyc = ix$MEAN$c, hyc_n = ix$DISP$c, hyc_k = ix$KBAL$c, hyt = ix$MEAN$t, hyt_n = ix$DISP$t, hyt_k = ix$KBAL$t, nHYB = nH, nSC = nP, nSE = nP)
p1 <- permute_contrasts_one("g1", E, FITS, M, perms, share_obs, seed = 1, gates = GATES, resamples = draws, methods = c("NB", "weighted"))
stopifnot(all(c("mean_cis_p", "disp_trans_p", "disp_trans_p_pooled", "kbal_cis_p", "kbal_trans_p", "kbalsup_cis_p", "disp_mid_p", "disp_par_sc_p", "disp_par_sc_p_ploidy",
                "disp_par_sc_p_ind", "kbal_par_se_p_ploidy", "wkbal_cis_p", "rkbal_mid_p", "disp_cis_infrac", "ploidy_shift") %in% names(p1)))
pv <- unlist(p1[grep("_p(_ploidy|_ind|_pooled)?$", names(p1))]); stopifnot(all(is.na(pv) | (pv > 0 & pv <= 1)))
## the gate NAs a contrast's estimates and p-values
G2 <- GATES; G2["g1", "log_cis"] <- FALSE; G2["g1", "kbal_trans"] <- FALSE
b2 <- boot_contrasts_one("g1", E, FITS, M, draws, share_cells, share_obs, gates = G2); p2 <- permute_contrasts_one("g1", E, FITS, M, perms, share_obs, gates = G2)
stopifnot(is.na(b2$disp_cis_est), is.na(b2$bsize_cis_se), is.na(b2$kbal_trans_est), !is.na(b2$mean_cis_est), !is.na(b2$disp_trans_est), is.na(p2$disp_cis_p), is.na(p2$kbal_trans_obs))
cat("contrast function checks passed\n")
