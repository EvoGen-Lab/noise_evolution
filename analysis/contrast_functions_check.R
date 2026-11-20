###############################################################
### contrast_functions_check.R
### Local checks of the contrast machinery (Methods: Contrasts, Tests of contrasts) on small synthetic designs: the weighted contrasts and the per-genome
### hybrid (corrected allele means, geometric and arithmetic means; averages of the allele values for raw theta, zeta and kappa), the bootstrap sign test
### (sign_boot_p(): null side, ties, undefined replicates, p floor), and a run of boot_contrasts_one() on 21 datasets with few draws, including the gate, the
### p-values, the equality of the log2 and raw dispersion tests for total and cis, and the resampling of make_draws() (undivided hybrid cells, memberships kept).
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
## the mean midpoints are geometric on the log2 scale too: mid and par compare the mean log2 allele rate with the mean log2 parental rate or one parent
mu5 <- list(MIX.SC = 2, MIX.SE = 8, HYB.SC = 1, HYB.SE = 9)
stopifnot(abs(.contrast_value("mid", mu5, "MU") - (0.5 * log2(1 * 9) - 0.5 * log2(2 * 8))) < 1e-12,
          abs(.contrast_value("par_se", mu5, "MU") - (0.5 * log2(1 * 9) - log2(8))) < 1e-12)
## the per-genome hybrid dispersion is the geometric mean of the allele dispersions on the log2 scale, the arithmetic mean on the raw scale
mc <- .mode_contrasts(gmu[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], ph, cb, modes = "par_sc")
stopifnot(colnames(mc) == c("mean", "disp", "raw"), abs(mc["par_sc", "disp"] - (0.5 * log2(0.1 * 0.3) - log2(0.2))) < 1e-12,
          abs(mc["par_sc", "raw"] - (0.5 * (0.1 + 0.3) - 0.2)) < 1e-12)

## ---- raw kinetic quantities: hybrid values are averages of the allele values, theta and kappa infinite at phi = 0 ----
mu4 <- list(MIX.SC = 2, MIX.SE = 3, HYB.SC = 2 * cb * 4, HYB.SE = 2 * (1 - cb) * 1)
stopifnot(abs(.contrast_burst_raw("mid", mu4, ph, "theta", cb) - (0.5 * (1 / 0.1 + 1 / 0.3) - 0.5 * (1 / 0.2 + 1 / 0.4))) < 1e-12,
          abs(.contrast_burst_raw("par_sc", mu4, ph, "bsize", cb) - (0.5 * (4 * 0.1 + 1 * 0.3) - 2 * 0.2)) < 1e-12,
          abs(.contrast_burst_raw("par_se", mu4, ph, "kbal", cb) - (0.5 * (1 / (4 * 0.1^2) + 1 / (1 * 0.3^2)) - 1 / (3 * 0.4^2))) < 1e-12,
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0.4), "theta") == Inf,
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0.2, MIX.SE = 0), "kbal") == -Inf,
          is.nan(.contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0), "theta")),
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0.4), "bsize") == -3 * 0.4)

## ---- the bootstrap sign test ----
## replicates on the null side: at or below zero for a positive estimate, at or above for a negative one; ties and undefined replicates (NA, NaN) on the null side
x <- c(-1, 0, NaN, NA, 2, 3, Inf, -Inf, 0.5, 1)
stopifnot(sign_boot_p(x, 1) == min(1, 2 * (1 + 5) / 11), sign_boot_p(x, -1) == min(1, 2 * (1 + 8) / 11), sign_boot_p(x, 0) == 1,
          sign_boot_p(rep(1, 999), 1) == 2 / 1000, sign_boot_p(rep(-1, 9), 1) == 1)

## ---- boot_contrasts_one() and null_contrasts_one() on 21 synthetic datasets ----
G <- 4; nH <- 400; nP <- 260; cbar0 <- 0.5024
gn <- paste0("g", seq_len(G)); mu <- c(2, 5, 3, 8); phi <- c(0.3, 0.5, 0.2, 0.15)
lib <- rlnorm(nH, log(2e4), 0.35); S <- lib / 2e4; ci <- rbeta(nH, cbar0 * 2e3, (1 - cbar0) * 2e3)
nb <- function(m, ph) matrix(rnbinom(length(m), mu = m, size = 1 / rep(ph, ncol(m))), nrow(m), dimnames = list(gn, NULL))
Hsc <- nb(outer(mu, S * ci), phi); Hse <- nb(outer(mu, S * (1 - ci)), phi)
eP <- rlnorm(nP, 0, 0.3); Msc <- nb(outer(mu, eP), phi); Mse <- nb(outer(mu, eP), phi)
libh <- round(lib); chat <- rbinom(nH, libh, ci) / libh
st <- share_terms(chat, libh); SR <- share_r_groups(st$r_c, st$r_e); eH <- S / 2
ix <- list(MEAN = split_indices_by_depth(nH, 0.4), DISP = split_indices_by_depth(nH, 0.35), KBAL = split_indices_by_depth(nH, 0.3), BSIZE = split_indices_by_depth(nH, 0.33))
M <- list(MIX.SC = Msc, MIX.SE = Mse, HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = eH, HYB.SE = eH, HYB.COMB = S)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K", BSIZE = ".Z"))
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, ix[[sp]][[side]]]; E[[nm]] <- eH[ix[[sp]][[side]]]
}
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
share_obs <- list(cbar = st$cbar, r = SR)
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
pcols <- as.vector(outer(c("mean", "disp", "rdisp", "rtheta", "bsize", "rbsize", "kbal", "rkbal"), .OUT_MODES, function(q, m) paste0(q, "_", m, "_p")))
stopifnot(nrow(bo) == G, all(c("mean_mid_est", "disp_par_sc_se", "kbal_cis_est", "bsize_trans_est", "wkbal_mid_est", "rkbal_total_est", "rbsize_cis_est", "rtheta_mid_lo",
                               "wrtheta_cis_est", "cor_cis", pcols) %in% names(bo)),
          !any(grepl("kbalsup|zeroflag|_zero$|^zero_|^wzero_|_nfin$|^w.*_p$", names(bo))),
          all(is.na(unlist(bo[pcols])) | (unlist(bo[pcols]) >= 2 / (length(draws) + 1) & unlist(bo[pcols]) <= 1)))
## for total and cis the log2 and raw dispersion tests give the same p-value, and raw theta (1 / phi) the same as well
for (md in c("total", "cis")) stopifnot(isTRUE(all.equal(bo[[paste0("disp_", md, "_p")]], bo[[paste0("rdisp_", md, "_p")]])),
                                        isTRUE(all.equal(bo[[paste0("disp_", md, "_p")]], bo[[paste0("rtheta_", md, "_p")]])))
## burst size comes from its own partition, with the mean and the dispersion of the same cells
fv1 <- fit_group_values(FITS[.CONTRAST.GROUPS], "g1"); oc1 <- .mode_contrasts(fv1$mu, fv1$disp, share_obs$cbar)
stopifnot(abs(bo$bsize_cis_est[1] - (oc1["cis_z", "mean"] + oc1["cis_z", "disp"])) < 1e-12, abs(bo$kbal_trans_est[1] - (-oc1["trans_k", "mean"] - 2 * oc1["trans_k", "disp"])) < 1e-12)
## the p-value is the sign test on the bootstrap replicates, which for disp_total are the log2 ratios of the refitted parental dispersions, a failed refit as zero
bt <- vapply(draws, function(d) { a <- .fit_one(M$MIX.SC["g1", d$MIX.SC], E$MIX.SC[d$MIX.SC]); b <- .fit_one(M$MIX.SE["g1", d$MIX.SE], E$MIX.SE[d$MIX.SE])
  f0 <- function(v) if (is.finite(v)) v else 0
  .disp_contrasts("total", list(MIX.SC = f0(a[["disp"]]), MIX.SE = f0(b[["disp"]])), zero_inf = TRUE)[["log"]] }, numeric(1))
stopifnot(abs(bo$disp_total_p[1] - sign_boot_p(bt, sign(bo$disp_total_est[1]))) < 1e-12)
## kept replicate values (boot_contrasts_one(keep = )): one column per kept contrast, the replicates the sign test reads
bk <- boot_contrasts_one("g1", E, FITS, M, draws, list(chat = chat, lib = libh), share_obs, gates = GATES, keep = c("disp_mid", "disp_par_sc"))
kd <- attr(bk, "draws")
stopifnot(identical(dim(kd), c(length(draws), 2L)), abs(sign_boot_p(kd[, "disp_mid"], sign(bk$disp_mid_est)) - bk$disp_mid_p) < 1e-12,
          abs(sign_boot_p(kd[, "disp_par_sc"], sign(bk$disp_par_sc_est)) - bk$disp_par_sc_p) < 1e-12)

## the shared-rise curve: a planted rise 2 - 0.25 x on calibrators is recovered, extrapolated flat beyond their range, and one curve subtracted from mid and
## both parental contrasts leaves their difference (minus the total contrast) unchanged in the estimate and in every replicate
set.seed(3); ng <- 600; xs <- runif(ng, -1, 8); calg <- runif(ng) < 0.5 & xs < 7
dmid <- 2 - 0.25 * xs + rnorm(ng, 0, 0.3); f0 <- dom_shift_curve(dmid, xs, calg, 0.3)
stopifnot(max(abs((f0 - (2 - 0.25 * pmin(xs, max(xs[calg]))))[calg])) < 0.25, abs(median(dmid[calg] - f0[calg])) < 0.05,
          isTRUE(all.equal(f0[xs > max(xs[calg])], rep(f0[which.max(replace(xs, !calg, -Inf))], sum(xs > max(xs[calg]))))))
tot <- rnorm(ng, 0, 0.5); dsc <- dmid - tot / 2; dse <- dmid + tot / 2
rep_ <- function(v) v + matrix(rnorm(ng * 30, 0, 0.2), ng, 30)
Dm_ <- rep_(dmid); Dsc_ <- Dm_ - tot / 2; Dse_ <- Dm_ + tot / 2
fb <- vapply(1:30, function(b) dom_shift_curve(Dm_[, b], xs, calg, 0.3), numeric(ng))
a_sc <- dom_adjust_contrast("par_sc", dsc, Dsc_, f0, fb); a_se <- dom_adjust_contrast("par_se", dse, Dse_, f0, fb)
stopifnot(max(abs((a_sc$adisp_par_sc_est - a_se$adisp_par_se_est) - (dsc - dse))) < 1e-12, all(c("adisp_par_sc_p", "adisp_par_sc_lo") %in% names(a_sc)))

## the gate NAs a contrast's estimates, intervals and p-values
G2 <- GATES; G2["g1", "log_cis"] <- FALSE; G2["g1", "kbal_trans"] <- FALSE; G2["g1", "bsize_cis"] <- FALSE
b2 <- boot_contrasts_one("g1", E, FITS, M, draws, list(chat = chat, lib = libh), share_obs, gates = G2)
stopifnot(is.na(b2$disp_cis_est), is.na(b2$bsize_cis_se), is.na(b2$kbal_trans_est), !is.na(b2$mean_cis_est), !is.na(b2$disp_trans_est),
          is.na(b2$disp_cis_p), is.na(b2$rtheta_cis_p), is.na(b2$kbal_trans_p), is.na(b2$rkbal_trans_p), !is.na(b2$rdisp_cis_p), !is.na(b2$disp_trans_p), !is.na(b2$bsize_trans_p))
## the pilot statistic counts a failed (non-finite) refit as a zero, over all replicates
stopifnot(zero_fail_frac(c(0, 0.1, NA, 0.2)) == 0.5, zero_fail_frac(c(0.1, NaN)) == 0.5, zero_fail_frac(c(0.1, 0.2)) == 0, is.na(zero_fail_frac(numeric(0))))

cat("contrast function checks passed\n")
