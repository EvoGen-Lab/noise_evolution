###############################################################
### contrast_functions_check.R
### Local checks of the contrast machinery (Methods: Contrasts, Tests of contrasts, Dominance, Shared shift) on small synthetic designs: the weighted contrasts
### and the per-genome hybrid (geometric means on the log2 scale, arithmetic on the raw scale; averages of the allele values for raw theta, beta and kappa), the
### dominance and inheritance contrasts on the trans cells, the bootstrap sign test (sign_boot_p(): null side, ties, undefined replicates, p floor), a run of
### boot_contrasts_one() on 21 datasets with few draws (gate, p-values, the equality of the log2 and raw dispersion tests for total and cis, the moment flag,
### the resampling of make_draws()), the shared shift (boot_contrasts_one(shift = ): the log2 dominance contrasts move by -Xi, -2 Xi of kappa reversed, cis and
### trans log2 unchanged, raw hybrid values divided by 2^Xi), the degree of dominance k and its distribution (.ratio_cdf()), and the inheritance classes
### (classify_dom()).
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(11)

## ---- contrasts ----
## equal per-genome rates give zero mean contrasts: each allele is offset by its own genome's library, so its fitted rate is per genome with no rescaling
mu0 <- 3
gmu <- setNames(as.list(rep(mu0, length(.CONTRAST.GROUPS))), .CONTRAST.GROUPS)
for (md in .OUT_MODES) stopifnot(abs(.contrast_value(md, gmu, "MU")) < 1e-12)
## dispersion contrasts: geometric means on the log2 scale, arithmetic on the raw scale; mid, par and inh use the trans cells
ph <- list(MIX.SC = 0.2, MIX.SE = 0.4, HYT.SC = 0.1, HYT.SE = 0.3, HYC.SC = 0.25, HYC.SE = 0.5)
stopifnot(abs(.contrast_value("mid", ph, "DISP") - (0.5 * log2(0.1 * 0.3) - 0.5 * log2(0.2 * 0.4))) < 1e-12,
          abs(.contrast_raw("mid", ph) - (0.5 * (0.1 + 0.3) - 0.5 * (0.2 + 0.4))) < 1e-12,
          abs(.contrast_value("par_se", ph, "DISP") - (0.5 * log2(0.1 * 0.3) - log2(0.4))) < 1e-12,
          abs(.contrast_value("inh_c", ph, "DISP") - log2(0.1 / 0.2)) < 1e-12, abs(.contrast_value("inh_e", ph, "DISP") - log2(0.3 / 0.4)) < 1e-12,
          abs(.contrast_value("cis", ph, "DISP") - log2(0.25 / 0.5)) < 1e-12)
## the identities of Derived comparisons on the log2 scale: par = dom -+ total / 2, inh = dom -+ trans / 2
v <- function(md) .contrast_value(md, ph, "DISP")
stopifnot(abs(v("par_sc") - (v("mid") - v("total") / 2)) < 1e-12, abs(v("par_se") - (v("mid") + v("total") / 2)) < 1e-12,
          abs(v("inh_c") - (v("mid") - v("trans") / 2)) < 1e-12, abs(v("inh_e") - (v("mid") + v("trans") / 2)) < 1e-12)
stopifnot(.contrast_value("total", list(MIX.SC = 0, MIX.SE = 1), "DISP", zero_inf = TRUE) == -Inf,
          is.na(.disp_contrasts("total", list(MIX.SC = 0, MIX.SE = 0), zero_inf = TRUE)[["log"]]))
## no contrast uses the allele fits of all hybrid cells or the summed counts
stopifnot(!any(c("HYB.SC", "HYB.SE", "HYB.COMB") %in% unlist(.MODES)), all(unlist(.MODES) %in% .CONTRAST.GROUPS))

## ---- raw kinetic quantities: hybrid values are averages of the allele values, theta and kappa infinite at phi = 0 ----
mu4 <- list(MIX.SC = 2, MIX.SE = 3, HYT.SC.N = 4, HYT.SE.N = 1, HYT.SC.Z = 4, HYT.SE.Z = 1, HYT.SC.K = 4, HYT.SE.K = 1)
ph4 <- list(MIX.SC = 0.2, MIX.SE = 0.4, HYT.SC.N = 0.1, HYT.SE.N = 0.3, HYT.SC.Z = 0.1, HYT.SE.Z = 0.3, HYT.SC.K = 0.1, HYT.SE.K = 0.3)
stopifnot(abs(.contrast_burst_raw("mid", mu4, ph4, "theta") - (0.5 * (1 / 0.1 + 1 / 0.3) - 0.5 * (1 / 0.2 + 1 / 0.4))) < 1e-12,
          abs(.contrast_burst_raw("par_sc", mu4, ph4, "bsize") - (0.5 * (4 * 0.1 + 1 * 0.3) - 2 * 0.2)) < 1e-12,
          abs(.contrast_burst_raw("par_se", mu4, ph4, "kbal") - (0.5 * (1 / (4 * 0.1^2) + 1 / (1 * 0.3^2)) - 1 / (3 * 0.4^2))) < 1e-12,
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0.4), "theta") == Inf,
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0.2, MIX.SE = 0), "kbal") == -Inf,
          is.nan(.contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0), "theta")),
          .contrast_burst_raw("total", mu4, list(MIX.SC = 0, MIX.SE = 0.4), "bsize") == -3 * 0.4)

## ---- the bootstrap sign test ----
x <- c(-1, 0, NaN, NA, 2, 3, Inf, -Inf, 0.5, 1)
stopifnot(sign_boot_p(x, 1) == min(1, 2 * (1 + 5) / 11), sign_boot_p(x, -1) == min(1, 2 * (1 + 8) / 11), sign_boot_p(x, 0) == 1,
          sign_boot_p(rep(1, 999), 1) == 2 / 1000, sign_boot_p(rep(-1, 9), 1) == 1)

## ---- the split variances and the split bias ----
sh <- c(0.3, 0.5, 0.7); L <- c(1000, 2000, 4000); et <- split_eta(sh, L)
stopifnot(all(abs(et$Hc - (1 - sh) / (sh * L)) < 1e-15), all(abs(et$He - sh / ((1 - sh) * L)) < 1e-15), all(abs(et$x - (et$Hc + et$He - 1 / L)) < 1e-15))
## the interpolated information of one cell matches the exact sum within 1%
for (M_ in c(0.05, 1, 30)) for (p_ in c(1e-3, 0.1, 2)) stopifnot(abs(nb_cell_info(M_, p_) / .nb_info_exact(M_, p_) - 1) < 0.01)
## the split bias is subtracted from the fitted phi*, and is zero without a split
y0 <- rnbinom(500, mu = 4, size = 1 / 0.2); e0 <- rep(1, 500)
f0 <- fit_dispersion(y0, e0, 0, "NB"); f1 <- fit_dispersion(y0, e0, rep(0.01, 500), "NB")
stopifnot(f0[["split_bias"]] == 0, f1[["split_bias"]] > 0, abs(f1[["disp"]] - (f1[["disp_app"]] - f1[["split_bias"]])) < 1e-12, f1[["disp_app"]] == f0[["disp_app"]])

## ---- boot_contrasts_one() on 21 synthetic datasets ----
G <- 4; nH <- 400; nP <- 260; sbar0 <- 0.5024
gn <- paste0("g", seq_len(G)); mu <- c(2, 5, 3, 8); phi <- c(0.3, 0.5, 0.2, 0.15)
lib <- rlnorm(nH, log(2e4), 0.35); si <- rbeta(nH, sbar0 * 2e3, (1 - sbar0) * 2e3); libh <- round(lib)
nb <- function(m, ph) matrix(rnbinom(length(m), mu = m, size = 1 / rep(ph, ncol(m))), nrow(m), dimnames = list(gn, NULL))
Hsc <- nb(outer(mu, lib * si / 2e4), phi); Hse <- nb(outer(mu, lib * (1 - si) / 2e4), phi)
eP <- rlnorm(nP, 0, 0.3); Msc <- nb(outer(mu, eP), phi); Mse <- nb(outer(mu, eP), phi)
Lc <- rbinom(nH, libh, si); shat <- Lc / libh; eta <- split_eta(shat, libh)
oc <- Lc / 2e4; oe <- (libh - Lc) / 2e4                                 # allele offsets L_Hc / D and L_He / D
ix <- list(MEAN = split_indices_by_depth(nH, 0.4), DISP = split_indices_by_depth(nH, 0.35), KBAL = split_indices_by_depth(nH, 0.3), BSIZE = split_indices_by_depth(nH, 0.33))
M <- list(MIX.SC = Msc, MIX.SE = Mse, HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = oc, HYB.SE = oe, HYB.COMB = libh / 2e4)
ET <- list(MIX.SC = 0, MIX.SE = 0, HYB.SC = eta$Hc, HYB.SE = eta$He, HYB.COMB = 0)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K", BSIZE = ".Z"))
  cells <- ix[[sp]][[side]]
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, cells]; E[[nm]] <- (if (al == "SC") oc else oe)[cells]; ET[[nm]] <- (if (al == "SC") eta$Hc else eta$He)[cells]
}
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], eta = ET[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
stopifnot(setequal(names(M), names(.DRAW.KEY)))
## the allele means estimate the per-genome rates with no rescaling
stopifnot(all(abs(log2(FITS$HYB.SC$MU / mu)) < 0.15), all(abs(log2(FITS$HYB.SE$MU / mu)) < 0.15))
gate_all <- setNames(rep(TRUE, 6 * length(.OUT_MODES)), as.vector(outer(c("mean", "raw", "log", "bsize", "rbsize", "kbal"), .OUT_MODES, paste, sep = "_")))
GATES <- matrix(TRUE, G, length(gate_all), dimnames = list(gn, names(gate_all)))
splits <- list(HYC = ix$MEAN$c, HYC.N = ix$DISP$c, HYC.K = ix$KBAL$c, HYC.Z = ix$BSIZE$c)
draws <- make_draws(c(MIX.SC = nP, MIX.SE = nP, HYB.COMB = nH), 10, 1, splits = splits)
for (d in draws) for (k in names(splits)) {
  kt <- sub("^HYC", "HYT", k)
  stopifnot(length(d[[k]]) + length(d[[kt]]) == nH, length(d[[k]]) == sum(d$HYB %in% splits[[k]]),
            all(sort(ix[[c(HYC = "MEAN", HYC.N = "DISP", HYC.K = "KBAL", HYC.Z = "BSIZE")[[k]]]]$c[d[[k]]]) == sort(d$HYB[d$HYB %in% splits[[k]]])))
}
bo <- do.call(rbind, lapply(gn, boot_contrasts_one, expos = E, fits = FITS, mats = M, draws = draws, etas = ET, gates = GATES, methods = c("NB", "weighted")))
pcols <- as.vector(outer(names(.NB.QUANTITIES), .OUT_MODES, function(q, m) paste0(q, "_", m, "_p")))
stopifnot(nrow(bo) == G, all(c("mean_mid_est", "disp_par_sc_se", "kbal_cis_est", "bsize_trans_est", "wkbal_mid_est", "rkbal_total_est", "rbsize_cis_est", "rtheta_mid_lo",
                               "wrtheta_cis_est", "disp_inh_c_est", "mean_inh_e_p", "wdisp_mid_flag", "cor_cis", pcols) %in% names(bo)),
          !any(grepl("^w.*_p$|^a|^k_|^kt_", names(bo))),
          all(is.na(unlist(bo[pcols])) | (unlist(bo[pcols]) >= 2 / (length(draws) + 1) & unlist(bo[pcols]) <= 1)))
## the moment flag: the weighted estimate outside the Cox-Reid percentile interval
stopifnot(identical(bo$wdisp_total_flag, bo$wdisp_total_est < bo$disp_total_lo | bo$wdisp_total_est > bo$disp_total_hi))
for (md in c("total", "cis")) stopifnot(isTRUE(all.equal(bo[[paste0("disp_", md, "_p")]], bo[[paste0("rdisp_", md, "_p")]])),
                                        isTRUE(all.equal(bo[[paste0("disp_", md, "_p")]], bo[[paste0("rtheta_", md, "_p")]])))
## burst size comes from its own partition, with the mean and the dispersion of the same cells
## (a gene with phi > 0 in every dataset; a zero observed phi leaves no log2 estimate)
gp <- which(vapply(gn, function(g) all(unlist(fit_group_values(FITS[.CONTRAST.GROUPS], g)$disp) > 0), logical(1)))[1]
stopifnot(!is.na(gp))
fvp <- fit_group_values(FITS[.CONTRAST.GROUPS], gn[gp]); ocp <- .mode_contrasts(fvp$mu, fvp$disp)
stopifnot(abs(bo$bsize_cis_est[gp] - (ocp["cis_z", "mean"] + ocp["cis_z", "disp"])) < 1e-12, abs(bo$kbal_trans_est[gp] - (-ocp["trans_k", "mean"] - 2 * ocp["trans_k", "disp"])) < 1e-12,
          abs(bo$disp_mid_est[gp] - ocp["mid_n", "disp"]) < 1e-12)
fv1 <- fit_group_values(FITS[.CONTRAST.GROUPS], "g1")
## the p-value is the sign test on the bootstrap replicates, which for disp_total are the log2 ratios of the refitted parental dispersions, a failed refit as zero
bt <- vapply(draws, function(d) { a <- .fit_one(M$MIX.SC["g1", d$MIX.SC], E$MIX.SC[d$MIX.SC]); b <- .fit_one(M$MIX.SE["g1", d$MIX.SE], E$MIX.SE[d$MIX.SE])
  f0 <- function(v) if (is.finite(v)) v else 0
  .disp_contrasts("total", list(MIX.SC = f0(a[["disp"]]), MIX.SE = f0(b[["disp"]])), zero_inf = TRUE)[["log"]] }, numeric(1))
stopifnot(abs(bo$disp_total_p[1] - sign_boot_p(bt, sign(bo$disp_total_est[1]))) < 1e-12)
## percentile intervals: an undefined replicate counts as -Inf for the lower percentile and +Inf for the upper
bs <- .boot_summary(c(NaN, NaN, 1:38)); bs2 <- .boot_summary(c(rep(NaN, 2), 1:118))
stopifnot(bs[["lo"]] == -Inf, bs[["hi"]] == Inf, bs[["se"]] == Inf, bs2[["lo"]] == 1, bs2[["hi"]] == 117)

## ---- the shared shift ----
## the curve: a planted rise 2 - 0.25 x on calibrators is recovered and extrapolated flat beyond their range
set.seed(3); ng <- 600; xs <- runif(ng, -1, 8); calg <- runif(ng) < 0.5 & xs < 7
dmid <- 2 - 0.25 * xs + rnorm(ng, 0, 0.3); cf <- dom_shift_curve(dmid, xs, calg, 0.3)
stopifnot(max(abs((cf - (2 - 0.25 * pmin(xs, max(xs[calg]))))[calg])) < 0.25, abs(median(dmid[calg] - cf[calg])) < 0.05,
          isTRUE(all.equal(cf[xs > max(xs[calg])], rep(cf[which.max(replace(xs, !calg, -Inf))], sum(xs > max(xs[calg]))))))
## the second stage: every hybrid-allele dispersion divided by 2^Xi, with Xi fixed at its observed value. The log2 contrasts of the hybrid with the parents move
## by -Xi (disp, bsize) and +2 Xi (kbal), in the estimate and in every replicate; cis and trans log2 do not move; the raw hybrid values are divided by 2^Xi
XI <- 0.8; SH <- setNames(rep(XI, G), gn)
ba <- boot_contrasts_one("g1", E, FITS, M, draws, ET, gates = GATES, shift = SH)
b0 <- boot_contrasts_one("g1", E, FITS, M, draws, ET, gates = GATES)
for (md in .DOM.MODES) stopifnot(abs(ba[[paste0("adisp_", md, "_est")]] - (b0[[paste0("disp_", md, "_est")]] - XI)) < 1e-12,
                                 abs(ba[[paste0("absize_", md, "_est")]] - (b0[[paste0("bsize_", md, "_est")]] - XI)) < 1e-12,
                                 abs(ba[[paste0("akbal_", md, "_est")]] - (b0[[paste0("kbal_", md, "_est")]] + 2 * XI)) < 1e-12,
                                 abs(ba[[paste0("adisp_", md, "_lo")]] - (b0[[paste0("disp_", md, "_lo")]] - XI)) < 1e-12)
stopifnot(!any(c("adisp_cis_est", "adisp_trans_est", "adisp_total_est", "ardisp_total_est", "amean_mid_est") %in% names(ba)),
          all(c("ardisp_cis_est", "ardisp_trans_est", "artheta_mid_p") %in% names(ba)))
fvn <- fv1$disp; fvn[startsWith(names(fvn), "HY")] <- lapply(fvn[startsWith(names(fvn), "HY")], function(z) z / 2^XI)
stopifnot(abs(ba$ardisp_cis_est - .contrast_raw("cis_n", fvn)) < 1e-12, abs(ba$ardisp_trans_est - .contrast_raw("trans_n", fvn)) < 1e-12,
          abs(ba$ardisp_cis_est - b0$rdisp_cis_est / 2^XI) < 1e-12)
## the degree of dominance k = 2 Delta_dom / Delta_total from the shifted Delta_dom, with its distribution at -1, 0 and 1
stopifnot(abs(ba$k_disp_est - 2 * ba$adisp_mid_est / ba$disp_total_est) < 1e-12, abs(ba$kt_disp_est - 2 * ba$adisp_mid_est / ba$disp_trans_est) < 1e-12,
          abs(ba$k_mean_est - 2 * ba$mean_mid_est / ba$mean_total_est) < 1e-12, ba$k_disp_F0 >= ba$k_disp_G0, ba$k_disp_F1 >= ba$k_disp_F0)
## the ratio rules: finite over infinite gives 0, infinite over finite an infinite k by the sign of the quotient, two infinities, a zero denominator or an
## undefined contrast an undefined replicate, counted in the lower tail (F) and not in the upper (G)
stopifnot(identical(.ratio_k(c(1, Inf, -Inf, Inf, 1, 0, NaN, 1), c(Inf, 2, 2, -Inf, 0, 0, 1, 2)), c(0, Inf, -Inf, NA, NA, NA, NA, 1)))
rc <- .ratio_cdf(1, 2, c(1, 2, -1, 0, 3, NaN), c(2, 1, 1, 1, 0, 1))      # k = 1, 4, -2, 0, undefined, undefined
stopifnot(rc[["est"]] == 1, rc[["n"]] == 6, rc[["n_undef"]] == 2, rc[["F1"]] == 5 / 6, rc[["G1"]] == 2 / 6, rc[["F0"]] == 4 / 6, rc[["G0"]] == 1 / 6,
          rc[["Fm1"]] == 3 / 6, rc[["Gm1"]] == 1 / 6, is.na(.ratio_cdf(1, 0, 1:5, 1:5)[["est"]]))

## ---- inheritance classes ----
## a family of 100 genes with 20 selected by both tests sets the interval level 1 - 0.05 * 20 / 100, so a = 0.01, a / 2 = 0.005
cdf <- function(Fm1, F0, F1) data.frame(Fm1 = Fm1, F0 = F0, F1 = F1, Gm1 = Fm1, G0 = F0, G1 = F1)
n <- 100; qr <- rep(0.5, n); qd <- rep(0.5, n); ed <- rep(1, n); kc <- cdf(rep(0.5, n), rep(0.5, n), rep(0.5, n))
qr[1:30] <- 0.001; qd[c(1:20, 31:35)] <- 0.001; ed[31] <- -1
set_k <- function(i, Fm1, F0, F1) { kc[i, ] <<- cdf(Fm1, F0, F1) }
set_k(1, 0.001, 0.2, 0.6)   # contains 0 and 1, not -1: Sc.Partial
set_k(2, 0.001, 0.002, 0.6) # contains 1 only: Sc.Dominant
set_k(3, 0.3, 0.9, 0.999)   # contains -1 and 0: Se.Partial
set_k(4, 0.3, 0.999, 0.9999) # contains -1 only: Se.Dominant
set_k(5, 0, 0, 0.001)       # entirely beyond 1: Overdominant (Delta_dom > 0)
set_k(6, 0, 0.001, 0.999)   # entirely between 0 and 1: Sc.Partial
set_k(7, 0.4, 0.5, 0.6)     # contains -1, 0 and 1: by the sign of Delta_dom
set_k(8, 0.001, 0.5, 0.999) # contains 0 only, within (-1, 1): Additive
cls <- classify_dom(qr, qd, ed, kc)
stopifnot(abs(attr(cls, "fcr_level") - 0.05 * 20 / 100) < 1e-15,
          identical(cls[1:8], c("Sc.Partial", "Sc.Dominant", "Se.Partial", "Se.Dominant", "Overdominant", "Sc.Partial", "Overdominant", "Additive")),
          all(cls[21:30] == "Additive"), cls[31] == "Underdominant", all(cls[32:35] == "Overdominant"), all(cls[36:100] == "Conserved"))
cls2 <- classify_dom(replace(qr, 1, NA), qd, ed, kc); stopifnot(is.na(cls2[1]))

## ---- co-expression (Methods: Co-expression) ----
## the contrasts of a pair's correlation follow the definitions on the cis and trans cells; the shared shift is the decile median of the observed dom, subtracted
## from dom and from both inheritance contrasts and held fixed in the replicates; the accumulator and the k tables match a direct computation; the bootstrap keeps
## each hybrid cell's partition and pairs its alleles; the null draws return four spectra
set.seed(5); pg <- 12; mk <- function(n, r) { f <- rnorm(n); z <- matrix(rnorm(pg * n), pg); z[1:4, ] <- z[1:4, ] + r * rep(f, each = 4); dimnames(z) <- list(paste0("g", 1:pg), NULL); z }
nH <- 200; cisz <- sort(sample(nH, 90)); trz <- setdiff(seq_len(nH), cisz); Hc <- mk(nH, 1); He <- mk(nH, 1)
RZ <- list(MIX.SC = mk(150, 0.8), MIX.SE = mk(170, 0.3), HYC.SC = Hc[, cisz], HYC.SE = He[, cisz], HYT.SC = Hc[, trz], HYT.SE = He[, trz])
upz <- which(upper.tri(diag(pg))); binz <- rep(1:3, length.out = length(upz)); attr(RZ, "dom_bins") <- binz
cd <- coexpr_decompose(RZ, shift = NULL); c0 <- coexpr_decompose(RZ, bins = NULL); attr(RZ, "dom_shift") <- cd$shift
R_ <- function(x) shrink_cor(t(x))
stopifnot(max(abs(cd$cis - (R_(RZ$HYC.SC) - R_(RZ$HYC.SE)))) < 1e-12, max(abs(cd$trans - (cd$total - (R_(RZ$HYT.SC) - R_(RZ$HYT.SE))))) < 1e-12,
          max(abs((cd$inh_c - (cd$dom - cd$trans / 2))[upz])) < 1e-12, max(abs((cd$inh_e - (cd$dom + cd$trans / 2))[upz])) < 1e-12,
          all(abs(tapply(cd$dom[upz], binz, median)) < 1e-12), max(abs((c0$dom - cd$dom)[upz] - cd$shift[binz])) < 1e-12, isSymmetric(cd$dom))
## the fixed shift: a decomposition with the observed shift equals the observed one, and a replicate subtracts the same shift
stopifnot(max(abs(coexpr_decompose(RZ)$dom - cd$dom)) < 1e-12)
dz <- make_coexpr_draws(150, 170, nH, 6, 1, cis = cisz)
stopifnot(all(vapply(dz, function(d) length(d$HC) + length(d$HT) == nH && all(d$HC <= length(cisz)) && all(d$HT <= length(trz)), logical(1))))
bz <- lapply(dz, coexpr_bootstrap_one, resid = RZ)
d1 <- dz[[1]]; r1 <- coexpr_decompose(list(MIX.SC = RZ$MIX.SC[, d1$SC], MIX.SE = RZ$MIX.SE[, d1$SE], HYC.SC = RZ$HYC.SC[, d1$HC], HYC.SE = RZ$HYC.SE[, d1$HC],
                                           HYT.SC = RZ$HYT.SC[, d1$HT], HYT.SE = RZ$HYT.SE[, d1$HT]), bins = NULL)
stopifnot(max(abs(bz[[1]]$dom - (r1$dom[upz] - cd$shift[binz]))) < 1e-12)
acc <- coexpr_acc_fold(coexpr_acc_fold(coexpr_acc_new(length(upz)), bz[1:3]), bz[4:6])
basez <- data.frame(gene_i = 1:length(upz))
pt <- coexpr_part_table("dom", cd, acc, basez); kt <- coexpr_k_table("k", cd, acc, basez)
kd <- sapply(bz, function(b) .ratio_k(b$dom, b$total)); ud <- is.na(kd)
stopifnot(acc$n == 6, max(abs(pt$se - apply(sapply(bz, `[[`, "dom"), 1, sd))) < 1e-10, all(kt$n == 6),
          max(abs(kt$F0 - rowMeans(ud | kd <= 0))) < 1e-12, max(abs(kt$G1 - rowMeans(!ud & kd < 1))) < 1e-12, max(abs(kt$est - .ratio_k(cd$dom[upz], cd$total[upz]))) < 1e-12)
nz <- coexpr_perm_one(list(perm = coexpr_perm_draw(1, length(cisz), length(trz), 320), dom = coexpr_dom_draw(1, 150, 170, length(trz))), 3, 150, RZ)
stopifnot(identical(names(nz), c("total", "cis", "trans", "dom")), all(lengths(nz) == 3))

## ---- the gate NAs a contrast's estimates, intervals and p-values ----
G2 <- GATES; G2["g1", "log_cis"] <- FALSE; G2["g1", "kbal_trans"] <- FALSE; G2["g1", "bsize_cis"] <- FALSE
b2 <- boot_contrasts_one("g1", E, FITS, M, draws, ET, gates = G2)
stopifnot(is.na(b2$disp_cis_est), is.na(b2$bsize_cis_se), is.na(b2$kbal_trans_est), !is.na(b2$mean_cis_est), !is.na(b2$disp_trans_est),
          is.na(b2$disp_cis_p), is.na(b2$rtheta_cis_p), is.na(b2$kbal_trans_p), is.na(b2$rkbal_trans_p), !is.na(b2$rdisp_cis_p), !is.na(b2$disp_trans_p), !is.na(b2$bsize_trans_p))
stopifnot(zero_fail_frac(c(0, 0.1, NA, 0.2)) == 0.5, zero_fail_frac(c(0.1, NaN)) == 0.5, zero_fail_frac(c(0.1, 0.2)) == 0, is.na(zero_fail_frac(numeric(0))))

cat("contrast function checks passed\n")
