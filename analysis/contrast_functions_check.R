###############################################################
### contrast_functions_check.R
### Local checks of the contrast machinery (Methods: Contrasts between groups, Null distributions) on small synthetic designs: the weighted contrasts and the
### per-genome hybrid (corrected allele means, geometric and arithmetic means), the null dispersions by projection (zero contrast for every contrast type on both
### scales), the kinetic projection, the simulated alleles with a shared extrinsic factor (target dispersions and allele covariance on many cells), and a run of
### boot_contrasts_one() and null_contrasts_one() on 21 datasets with few draws, including the gate, the group zero fractions and the shifted bootstrap nulls, the
### resampling of make_draws() (undivided hybrid cells, memberships kept) and the likelihood bounds of dispersion_bound_stats().
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(11)

## ---- contrasts ----
## equal per-genome rates give zero mean contrasts at any c-bar: the allele fits carry the abundance factors 2 c-bar and 2 (1 - c-bar)
cb <- 0.6; mu0 <- 3
gmu <- list(MIX.SC = mu0, MIX.SE = mu0, HYB.SC = 2 * cb * mu0, HYB.SE = 2 * (1 - cb) * mu0, HYC.SC = 2 * cb * mu0, HYC.SE = 2 * (1 - cb) * mu0, HYT.SC = 2 * cb * mu0, HYT.SE = 2 * (1 - cb) * mu0)
cm <- .correct_means(gmu, cb)
for (md in c("total", "cis", "trans", "mid", "par_sc", "par_se")) stopifnot(abs(.contrast_value(md, cm, "MU")) < 1e-12)
## the cis mean contrast equals the old offset form log2 mu*_Hc - log2 mu*_He - log2[c-bar / (1 - c-bar)]
g2 <- gmu; g2$HYC.SC <- 5; g2$HYC.SE <- 2
stopifnot(abs(.contrast_value("cis", .correct_means(g2, cb), "MU") - (log2(5 / 2) - log2(cb / (1 - cb)))) < 1e-12)
## dispersion contrasts: geometric means on the log2 scale, arithmetic on the raw scale
ph <- list(MIX.SC = 0.2, MIX.SE = 0.4, HYB.SC = 0.1, HYB.SE = 0.3)
stopifnot(abs(.contrast_value("mid", ph, "DISP") - (0.5 * log2(0.1 * 0.3) - 0.5 * log2(0.2 * 0.4))) < 1e-12,
          abs(.contrast_raw("mid", ph) - (0.5 * (0.1 + 0.3) - 0.5 * (0.2 + 0.4))) < 1e-12,
          abs(.contrast_value("par_se", ph, "DISP") - (0.5 * log2(0.1 * 0.3) - log2(0.4))) < 1e-12)
stopifnot(.contrast_value("total", list(MIX.SC = 0, MIX.SE = 1), "DISP", zero_inf = TRUE) == -Inf,
          is.na(.disp_contrasts("total", list(MIX.SC = 0, MIX.SE = 0), zero_inf = TRUE)[["log"]]))
## the per-genome hybrid dispersion is the geometric mean of the allele dispersions on the log2 scale, the arithmetic mean on the raw scale
mc <- .mode_contrasts(gmu[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], ph, cb, modes = "par_sc")
stopifnot(colnames(mc) == c("mean", "disp", "raw"), abs(mc["par_sc", "disp"] - (0.5 * log2(0.1 * 0.3) - log2(0.2))) < 1e-12,
          abs(mc["par_sc", "raw"] - (0.5 * (0.1 + 0.3) - 0.2)) < 1e-12)

## ---- null dispersions: the projection zeroes every contrast type on both scales ----
for (md in c("total", "cis", "trans", "mid", "par_sc")) {
  lam <- .MODE.WEIGHTS[[md]]; p0 <- runif(length(lam), 0.05, 1)
  stopifnot(abs(sum(lam * log2(null_phi_log2(p0, lam)))) < 1e-12, abs(sum(lam * null_phi_raw(p0, lam))) < 1e-12, all(null_phi_raw(p0, lam) >= 0))
  dmu <- rnorm(1)
  pb <- null_phi_kinetic(p0, lam, dmu, 1, 1); pk <- null_phi_kinetic(p0, lam, dmu, -1, -2)
  stopifnot(abs(dmu + sum(lam * log2(pb))) < 1e-12, abs(-dmu - 2 * sum(lam * log2(pk))) < 1e-12)
}
## the raw projection with truncation: a group pushed below zero is held at zero and the rest re-projected
r0 <- null_phi_raw(c(0.01, 0.9, 0.02, 0.05), .MODE.WEIGHTS$trans); stopifnot(all(r0 >= 0))

## ---- simulated alleles: target dispersions and the shared factor ----
n <- 400000; m <- rep(50, n)
s <- sim_alleles(m, m, 0.3, 0.15, 0.08)
fa <- function(y) (var(y) - mean(y)) / mean(y)^2                                      # moment dispersion of a constant-mean count
stopifnot(abs(fa(s$hc) - 0.3) < 0.01, abs(fa(s$he) - 0.15) < 0.01, abs(cov(s$hc / 50, s$he / 50) - 0.08) < 0.005)
s0 <- sim_alleles(m, m, 0.02, 0.02, 0.08); stopifnot(abs(fa(s0$hc) - 0.08) < 0.01)   # a private variance below zero is held at zero

## ---- likelihood bounds: lower < phi-hat < upper for an overdispersed gene, lower = 0 for a Poisson gene ----
eb <- rlnorm(2000, 0, 0.3)
bs <- dispersion_bound_stats(rnbinom(2000, mu = 4 * eb, size = 1 / 0.4), eb)
stopifnot(bs[["disp_lower"]] > 0, bs[["disp_lower"]] < bs[["disp_upper"]], bs[["lrt"]] > BOUND.CRIT)
bp <- dispersion_bound_stats(rpois(2000, 4 * eb), eb); stopifnot(bp[["disp_upper"]] > 0, bp[["disp_lower"]] >= 0)
stopifnot(abs(BOUND.CRIT - 2.7055) < 1e-4)

## ---- boot_contrasts_one() and null_contrasts_one() on 21 synthetic datasets ----
G <- 4; nH <- 400; nP <- 260; cbar0 <- 0.5024
gn <- paste0("g", seq_len(G)); mu <- c(2, 5, 3, 8); phi <- c(0.3, 0.5, 0.2, 0.15)
lib <- rlnorm(nH, log(2e4), 0.35); S <- lib / 2e4; ci <- rbeta(nH, cbar0 * 2e3, (1 - cbar0) * 2e3)
nb <- function(m, ph) matrix(rnbinom(length(m), mu = m, size = 1 / rep(ph, ncol(m))), nrow(m), dimnames = list(gn, NULL))
Hsc <- nb(outer(mu, S * ci), phi); Hse <- nb(outer(mu, S * (1 - ci)), phi)
eP <- rlnorm(nP, 0, 0.3); Msc <- nb(outer(mu, eP), phi); Mse <- nb(outer(mu, eP), phi)
libh <- round(lib); chat <- rbinom(nH, libh, ci) / libh
st <- share_terms(chat, libh); SR <- share_r_groups(st$r_c, st$r_e); eH <- S / 2
pc <- ploidy_components(Hsc, Hse, eH, st)
ix <- list(MEAN = split_indices_by_depth(nH, 0.4), DISP = split_indices_by_depth(nH, 0.35), KBAL = split_indices_by_depth(nH, 0.3), BSIZE = split_indices_by_depth(nH, 0.33))
M <- list(MIX.SC = Msc, MIX.SE = Mse, HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = eH, HYB.SE = eH, HYB.COMB = S)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K", BSIZE = ".Z"))
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, ix[[sp]][[side]]]; E[[nm]] <- eH[ix[[sp]][[side]]]
}
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
share_obs <- list(cbar = st$cbar, r = SR, sigma_d = setNames(pc$sigma_d, gn))
stopifnot(setequal(names(M), names(.DRAW.KEY)))
gate_all <- setNames(rep(TRUE, 6 * length(.OUT_MODES)), as.vector(outer(c("mean", "raw", "log", "bsize", "rbsize", "kbal"), .OUT_MODES, paste, sep = "_")))
GATES <- matrix(TRUE, G, length(gate_all), dimnames = list(gn, names(gate_all)))
splits <- list(HYC = ix$MEAN$c, HYC.N = ix$DISP$c, HYC.K = ix$KBAL$c, HYC.Z = ix$BSIZE$c)
draws <- make_draws(c(MIX.SC = nP, MIX.SE = nP, HYB.COMB = nH), 10, 1, splits = splits)
## every resampled hybrid cell lands on the side of each partition its original was on, and the sides together hold every copy
for (d in draws) for (k in names(splits)) {
  kt <- sub("^HYC", "HYT", k)
  stopifnot(length(d[[k]]) + length(d[[kt]]) == nH, length(d[[k]]) == sum(d$HYB %in% splits[[k]]),
            all(sort(ix[[c(HYC = "MEAN", HYC.N = "DISP", HYC.K = "KBAL", HYC.Z = "BSIZE")[[k]]]]$c[d[[k]]]) == sort(d$HYB[d$HYB %in% splits[[k]]])))
}
bo <- do.call(rbind, lapply(gn, boot_contrasts_one, expos = E, fits = FITS, mats = M, draws = draws, share_cells = list(chat = chat, lib = libh), share_obs = share_obs,
                            gates = GATES, methods = c("NB", "weighted")))
stopifnot(nrow(bo) == G, all(c("mean_mid_est", "disp_par_sc_se", "disp_cis_zero", "kbal_cis_est", "bsize_trans_est", "wkbal_mid_est", "rkbal_total_est", "rbsize_cis_est",
                               "cor_cis", paste0("zero_", .CONTRAST.GROUPS), paste0("wzero_", .CONTRAST.GROUPS)) %in% names(bo)),
          !any(grepl("kbalsup|zeroflag", names(bo))), all(bo[, paste0("zero_", .CONTRAST.GROUPS)] >= 0 & bo[, paste0("zero_", .CONTRAST.GROUPS)] <= 1))
## burst size comes from its own partition, with the mean and the dispersion of the same cells
fv1 <- fit_group_values(FITS[.CONTRAST.GROUPS], "g1"); oc1 <- .mode_contrasts(fv1$mu, fv1$disp, share_obs$cbar)
stopifnot(abs(bo$bsize_cis_est[1] - (oc1["cis_z", "mean"] + oc1["cis_z", "disp"])) < 1e-12, abs(bo$kbal_trans_est[1] - (-oc1["trans_k", "mean"] - 2 * oc1["trans_k", "disp"])) < 1e-12)
perms <- lapply(1:30, perm_label_draw, hyc = ix$MEAN$c, hyt = ix$MEAN$t, nHYB = nH, nSC = nP, nSE = nP)
p1 <- null_contrasts_one("g1", E, FITS, M, perms, draws, list(chat = chat, lib = libh), share_obs, n_sim = 8, seed = 1, gates = GATES, methods = c("NB", "weighted"))
need <- c("mean_total_p", "mean_cis_p", "mean_trans_p", "mean_mid_p", "mean_par_sc_p", "disp_total_p", "disp_cis_p", "disp_trans_p", "disp_mid_p", "disp_par_se_p",
          "rdisp_cis_p", "rdisp_trans_p", "rdisp_mid_p", "bsize_cis_p", "kbal_trans_p", "wdisp_mid_p", "wkbal_par_sc_p", "wrdisp_total_p", "disp_cis_infrac", "k_trans")
stopifnot(all(need %in% names(p1)))
pv <- unlist(p1[grep("_p$", names(p1))]); stopifnot(all(is.na(pv) | (pv > 0 & pv <= 1)), all(!is.na(unlist(p1[need[need != "k_trans"]]))))
stopifnot(!any(grepl("cv2", c(names(p1), names(bo)))))
stopifnot(!any(grepl("kbalsup", names(p1))))
## the shifted bootstrap null of a log2 dispersion contrast is the resampled contrast minus the observed one: with the same draws, the null draws of disp_total
## are the bootstrap draws shifted by the estimate, so the p-value equals the share of bootstrap draws at least |estimate| away from it
bt <- vapply(draws, function(d) { a <- .fit_one(M$MIX.SC["g1", d$MIX.SC], E$MIX.SC[d$MIX.SC]); b <- .fit_one(M$MIX.SE["g1", d$MIX.SE], E$MIX.SE[d$MIX.SE])
  .disp_contrasts("total", list(MIX.SC = a[["disp"]], MIX.SE = b[["disp"]]), zero_inf = TRUE)[["log"]] }, numeric(1))
obs <- p1$disp_total_obs; stopifnot(abs(p1$disp_total_p - perm_pval(obs, bt - obs)) < 1e-12)
## the gate NAs a contrast's estimates and p-values and skips its nulls
G2 <- GATES; G2["g1", "log_cis"] <- FALSE; G2["g1", "kbal_trans"] <- FALSE; G2["g1", "bsize_cis"] <- FALSE
b2 <- boot_contrasts_one("g1", E, FITS, M, draws, list(chat = chat, lib = libh), share_obs, gates = G2)
p2 <- null_contrasts_one("g1", E, FITS, M, perms, draws[1:4], list(chat = chat, lib = libh), share_obs, n_sim = 4, gates = G2)
stopifnot(is.na(b2$disp_cis_est), is.na(b2$bsize_cis_se), is.na(b2$kbal_trans_est), !is.na(b2$mean_cis_est), !is.na(b2$disp_trans_est),
          is.na(p2$disp_cis_p), is.na(p2$kbal_trans_obs), !is.na(p2$disp_trans_p), !is.na(p2$bsize_trans_p))
cat("contrast function checks passed\n")
