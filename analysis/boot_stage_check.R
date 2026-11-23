###############################################################
### boot_stage_check.R
### Local unit checks of the dispersion boundary, the zero-tolerant sign test, the contrast-level shared shift and the staged bootstrap (Methods: Estimating phi,
### Failed and zero refits, Hybrid noise shift, Bootstrap tests of significance, Dominance), run before any cluster job:
###   boundary rule   on simulated Poisson and NB counts, phi from the zero rule plus the interior check is the maximum of the adjusted likelihood on a direct grid,
###                   and a search that returns no finite maximum is restarted and still finds it
###   sign forms      wherever every fit is positive, the product-form side of zero of every log2 contrast (mean, phi, beta, kappa; theta by reversal) equals
###                   the sign of the centered log2 contrast, for all eight modes and shifted or not; with zeros, a zero side gives the other side's sign
###                   and two zero sides an undefined replicate
###   identities      on one cell set, centered Delta_inh,c = Delta_dom - Delta_trans / 2 and Delta_H,Pc = Delta_dom - Delta_total / 2 to machine precision
###   Xi fixed        boot_contrasts_one(shift = ) applies the same Xi in every replicate: each shifted replicate of disp_mid is the unshifted one minus Xi, so
###                   the shifted p-value is the sign test of those values; genes beyond the calibrator range take the end value of the curve
###   raw factors     raw Delta_dom built with the per-quantity factors (2^Xi for phi and beta, 2^-Xi for theta, 2^-2Xi for kappa) has the sign of log2
###                   Delta_dom when the two parental values are equal and the two hybrid allele values are equal
###   fit status      the fitted, boundary and failed refits of every dataset add up to the number of replicates
###   staging         draw plans are reproducible and fresh across seeds and equal make_draws(); the mean family refits means alone and matches the full run;
###                   stage_genes() reads n0 back from p; merge_boot_stage() replaces one family's columns; the k interval at a coverage level
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(21)

## ---- boundary rule ----
for (i in 1:120) {
  n <- sample(c(60, 400), 1); e <- rlnorm(n, 0, 0.4); mu <- exp(runif(1, log(0.05), log(15))); phi <- sample(c(0, 0, 0.02, 0.3), 1)
  y <- if (phi == 0) rpois(n, mu * e) else rnbinom(n, mu = mu * e, size = 1 / phi)
  if (sum(y) == 0) next
  f <- fit_dispersion(y, e, 0, "NB"); ll <- .cr_loglik(y, sum(y) / sum(e) * e)
  grid <- c(0, exp(seq(log(3e-7), log(55), length.out = 3000))); v <- vapply(grid, ll, numeric(1))
  stopifnot(f[["converged"]] == 1, ll(f[["disp"]]) >= max(v) - 1e-6 * max(1, abs(max(v))))
  if (f[["disp"]] == 0) stopifnot(f[["interior"]] == 0, max(v) <= ll(0) + CR.INTERIOR.TOL * max(1, abs(ll(0))) + 1e-9)
}
## a search that returns no finite maximum is restarted from CR.RESTART.LOGPHI and the moment estimate, and reaches the same maximum
y <- rnbinom(300, mu = 3, size = 1 / 0.4); e <- rep(1, 300); f0 <- fit_dispersion(y, e, 0, "NB")
optimize <- function(...) stop("no finite maximum")
f1 <- fit_dispersion(y, e, 0, "NB")
rm(optimize)
stopifnot(f1[["restarted"]] == 1, f1[["converged"]] == 1, abs(log(f1[["disp"]] / f0[["disp"]])) < 1e-3)

## ---- sign forms ----
grp <- .CONTRAST.GROUPS
for (r in 1:200) {
  gm <- setNames(as.list(exp(rnorm(length(grp), 0, 1.5))), grp); gp <- setNames(as.list(exp(rnorm(length(grp), -2, 1.5))), grp)
  xi <- rnorm(1, 0, 0.7); gps <- gp; gps[startsWith(grp, "HY")] <- lapply(gps[startsWith(grp, "HY")], function(v) v * 2^(-xi))
  for (pp in list(gp, gps)) {
    v <- .all_contrasts(gm, pp); sg <- .all_signs(gm, pp, v)
    stopifnot(all(sg == sign(v)))
  }
  ## the shifted dispersion form: Delta_dom positive exactly when phi_Hc phi_He > 2^(2 Xi) phi_Pc phi_Pe (direct products), kappa with 2^(4 Xi) and reversed
  v <- .all_contrasts(gm, gps); sg <- .all_signs(gm, gps, v)
  stopifnot(sg[["disp_mid"]] == sign(gp$HYT.SC.N * gp$HYT.SE.N - 2^(2 * xi) * gp$MIX.SC * gp$MIX.SE),
            sg[["disp_inh_c"]] == sign(gp$HYT.SC.N - 2^xi * gp$MIX.SC),
            sg[["disp_par_se"]] == sign(gp$HYT.SC.N * gp$HYT.SE.N - 2^(2 * xi) * gp$MIX.SE^2),
            sg[["disp_trans"]] == sign(gp$MIX.SC * gp$HYT.SE.N - gp$MIX.SE * gp$HYT.SC.N),
            sg[["kbal_mid"]] == sign(2^(4 * xi) * gm$MIX.SC * gp$MIX.SC^2 * gm$MIX.SE * gp$MIX.SE^2 - gm$HYT.SC.K * gp$HYT.SC.K^2 * gm$HYT.SE.K * gp$HYT.SE.K^2),
            sg[["bsize_mid"]] == sign(gm$HYT.SC.Z * gp$HYT.SC.Z * gm$HYT.SE.Z * gp$HYT.SE.Z - 2^(2 * xi) * gm$MIX.SC * gp$MIX.SC * gm$MIX.SE * gp$MIX.SE))
}
## zeros: one zero side gives the other side's sign, two zero sides an undefined replicate; the log2 value agrees where it is defined
gp0 <- setNames(as.list(rep(0.2, length(grp))), grp); gm0 <- setNames(as.list(rep(2, length(grp))), grp)
z <- gp0; z$MIX.SE <- 0
sg <- .all_signs(gm0, z, .all_contrasts(gm0, z))
stopifnot(sg[["disp_total"]] == 1, sg[["disp_trans"]] == 1, sg[["disp_mid"]] == 1, sg[["kbal_total"]] == -1, sg[["bsize_total"]] == 1)
z$HYT.SC.N <- 0                                                     # trans: phi_Pc phi_He > phi_Pe phi_Hc is 0.04 against 0: positive
sg <- .all_signs(gm0, z, .all_contrasts(gm0, z)); stopifnot(sg[["disp_trans"]] == 1, is.na(sg[["disp_mid"]]))
z2 <- gp0; z2$MIX.SC <- 0; z2$MIX.SE <- 0
sg <- .all_signs(gm0, z2, .all_contrasts(gm0, z2)); stopifnot(is.na(sg[["disp_total"]]), is.na(sg[["kbal_total"]]), is.na(sg[["rtheta_total"]]), sg[["rdisp_total"]] == 0)
m0 <- gm0; m0$MIX.SE <- 0                                            # a zero mean is handled like a zero dispersion
sg <- .all_signs(m0, gp0, .all_contrasts(m0, gp0)); stopifnot(sg[["mean_total"]] == 1, sg[["kbal_total"]] == -1)
## the sign test with the zero-tolerant sides: undefined and tied replicates on the null side
stopifnot(sign_boot_p(c(1, 1, NA, 0, -1), 1) == min(1, 2 * (1 + 3) / 6))

## ---- identities and the raw factors ----
for (r in 1:50) {
  gm <- setNames(as.list(exp(rnorm(length(grp)))), grp); gp <- setNames(as.list(exp(rnorm(length(grp), -2))), grp)
  xi <- rnorm(1, 0, 0.7); gps <- gp; gps[startsWith(grp, "HY")] <- lapply(gps[startsWith(grp, "HY")], function(v) v * 2^(-xi))
  a <- .all_contrasts(gm, gps)
  for (q in c("mean", "disp", "bsize", "kbal"))
    stopifnot(abs(a[[paste0(q, "_inh_c")]] - (a[[paste0(q, "_mid")]] - a[[paste0(q, "_trans")]] / 2)) < 1e-12,
              abs(a[[paste0(q, "_inh_e")]] - (a[[paste0(q, "_mid")]] + a[[paste0(q, "_trans")]] / 2)) < 1e-12,
              abs(a[[paste0(q, "_par_sc")]] - (a[[paste0(q, "_mid")]] - a[[paste0(q, "_total")]] / 2)) < 1e-12,
              abs(a[[paste0(q, "_par_se")]] - (a[[paste0(q, "_mid")]] + a[[paste0(q, "_total")]] / 2)) < 1e-12)
  ## raw factors: equal parental values and equal hybrid allele values, so the raw and log2 Delta_dom share a sign for phi, theta, beta and kappa
  e <- gm; f <- gp
  for (k in grp) { e[[k]] <- if (startsWith(k, "HY")) gm$HYT.SC else gm$MIX.SC; f[[k]] <- if (startsWith(k, "HY")) gp$HYT.SC else gp$MIX.SC }
  fs <- f; fs[startsWith(grp, "HY")] <- lapply(fs[startsWith(grp, "HY")], function(v) v * 2^(-xi))
  a <- .all_contrasts(e, fs)
  stopifnot(sign(a[["rdisp_mid"]]) == sign(a[["disp_mid"]]), sign(a[["rtheta_mid"]]) == -sign(a[["disp_mid"]]),
            sign(a[["rbsize_mid"]]) == sign(a[["bsize_mid"]]), sign(a[["rkbal_mid"]]) == sign(a[["kbal_mid"]]))
}

## ---- a small synthetic design for boot_contrasts_one() ----
G <- 3; nH <- 300; nP <- 200; gn <- paste0("g", seq_len(G)); mu <- c(2, 5, 0.3); phi <- c(0.3, 0.5, 0.05)
lib <- rlnorm(nH, log(2e4), 0.35); si <- rbeta(nH, 1e3, 1e3); libh <- round(lib)
nb <- function(m, ph) matrix(rnbinom(length(m), mu = m, size = 1 / rep(ph, ncol(m))), nrow(m), dimnames = list(gn, NULL))
Hsc <- nb(outer(mu, lib * si / 2e4), phi); Hse <- nb(outer(mu, lib * (1 - si) / 2e4), phi * 1.5)
eP <- rlnorm(nP, 0, 0.3); Msc <- nb(outer(mu, eP), phi); Mse <- nb(outer(mu * 0.7, eP), phi)
Lc <- rbinom(nH, libh, si); eta <- split_eta(Lc / libh, libh); oc <- Lc / 2e4; oe <- (libh - Lc) / 2e4
ix <- list(MEAN = split_indices_by_depth(nH, 0.4), DISP = split_indices_by_depth(nH, 0.35), KBAL = split_indices_by_depth(nH, 0.3), BSIZE = split_indices_by_depth(nH, 0.33))
M <- list(MIX.SC = Msc, MIX.SE = Mse, HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = oc, HYB.SE = oe, HYB.COMB = libh / 2e4)
ET <- list(MIX.SC = 0, MIX.SE = 0, HYB.SC = eta$Hc, HYB.SE = eta$He, HYB.COMB = 0)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K", BSIZE = ".Z")); cells <- ix[[sp]][[side]]
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, cells]; E[[nm]] <- (if (al == "SC") oc else oe)[cells]; ET[[nm]] <- (if (al == "SC") eta$Hc else eta$He)[cells]
}
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], eta = ET[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
splits <- list(HYC = ix$MEAN$c, HYC.N = ix$DISP$c, HYC.K = ix$KBAL$c, HYC.Z = ix$BSIZE$c); NC <- c(MIX.SC = nP, MIX.SE = nP, HYB.COMB = nH)

## ---- draw plans ----
pl <- make_draw_plan(NC, 30, 4, splits); pl2 <- make_draw_plan(NC, 30, 5, splits)
stopifnot(identical(plan_draw(pl, 7), plan_draw(pl, 7)), identical(make_draws(NC, 30, 4, splits), lapply(1:30, plan_draw, plan = pl)),
          !identical(plan_draw(pl, 7)$HYB, plan_draw(pl2, 7)$HYB), !identical(plan_draw(pl, 7)$HYB, plan_draw(pl, 8)$HYB))
set.seed(99); r1 <- runif(1); set.seed(99); invisible(plan_draw(pl, 3)); r2 <- runif(1); stopifnot(r1 == r2)     # a draw leaves the caller's stream alone

## ---- Xi fixed, fit status, families ----
XI <- 0.6; SH <- setNames(rep(XI, G), gn)
b0 <- boot_contrasts_one("g1", E, FITS, M, pl, ET, shift = NULL)
ba <- boot_contrasts_one("g1", E, FITS, M, pl, ET, shift = SH)
## the shifted dispersion replicates of disp_mid: the unshifted replicates minus Xi, recomputed here from the same draws
dmid <- vapply(1:30, function(b) { d <- plan_draw(pl, b)
  f <- lapply(setNames(grp, grp), function(g) { i <- d[[.DRAW.KEY[[g]]]]; .fit_one(M[[g]]["g1", i], E[[g]][i], if (length(ET[[g]]) > 1) ET[[g]][i] else ET[[g]]) })
  ph <- lapply(f, function(x) if (is.finite(x[["disp"]])) unname(x[["disp"]]) else 0); mu_ <- lapply(f, function(x) unname(x[["mu"]]))
  .all_contrasts(mu_, ph)[["disp_mid"]] }, numeric(1))
stopifnot(abs(ba$adisp_mid_p - sign_boot_p(sign(dmid - XI), sign(ba$adisp_mid_est))) < 1e-12, abs(b0$disp_mid_p - sign_boot_p(sign(dmid), sign(b0$disp_mid_est))) < 1e-12,
          isTRUE(all.equal(ba$adisp_mid_lo, .boot_summary(dmid - XI)[["lo"]])), abs(ba$adisp_mid_est - (b0$disp_mid_est - XI)) < 1e-12)
for (g in grp) stopifnot(ba[[paste0("fs_", g, "_fit")]] + ba[[paste0("fs_", g, "_bnd")]] + ba[[paste0("fs_", g, "_fail")]] == 30)
## the mean family refits means alone, with the same estimates and replicates as the full run; the dispersion family keeps every non-mean column
bm <- boot_contrasts_one("g1", E, FITS, M, pl, ET, shift = SH, family = "mean"); bd <- boot_contrasts_one("g1", E, FITS, M, pl, ET, shift = SH, family = "disp")
stopifnot(all(grepl("^mean_|^gene$|^n_boot$", names(bm))), isTRUE(all.equal(unlist(bm[-1]), unlist(ba[names(bm)][-1]))),
          isTRUE(all.equal(unlist(bd[-1]), unlist(ba[names(bd)][-1]))), !any(grepl("^mean_|^w|^cor_", names(bd))), all(boot_family_cols(names(ba), "disp") %in% names(bd)))
## the k interval at a coverage level: with a = 0.05 it is the 95% interval
bk <- boot_contrasts_one("g1", E, FITS, M, pl, ET, shift = SH, family = "disp", k_level = c(k_disp = 0.05, kt_disp = 0.2))
stopifnot(all(grepl("^kt?_|^gene$|^n_boot$", names(bk))), identical(bk$k_disp_flo, bk$k_disp_lo), identical(bk$k_disp_fhi, bk$k_disp_hi),
          is.na(bk$kt_disp_flo) || bk$kt_disp_flo >= bk$kt_disp_lo, !("k_mean_flo" %in% names(bk)))
## the dominance degree's interval rule at a level: an undefined replicate is -Inf below and +Inf above
rc <- .ratio_cdf(1, 2, c(1:38, NaN, NaN), rep(2, 40), a = 0.2); rc2 <- .ratio_cdf(1, 2, c(1:38, NaN, NaN), rep(2, 40), a = 0.05)
stopifnot(rc[["flo"]] == 2, rc[["fhi"]] == 36, rc2[["flo"]] == -Inf, rc2[["fhi"]] == Inf)

## ---- stage_genes() and merge_boot_stage() ----
bc <- data.frame(gene = c("a", "b", "c"), n_boot_mean = 1000, n_boot_disp = c(1000, 1000, 10000), mean_total_est = 1, mean_total_se = 0.1, mean_total_lo = 0, mean_total_hi = 2,
                 mean_total_p = c(2 * 40 / 1001, 2 * 41 / 1001, 1), disp_total_est = 1, disp_total_se = 1, disp_total_lo = 0, disp_total_hi = 1,
                 disp_total_p = c(1, 2 * 10 / 1001, 2 * 30 / 10001), stringsAsFactors = FALSE)
stopifnot(identical(stage_genes(bc, "mean", 40), "a"), identical(stage_genes(bc, "disp", 40, 1000), "b"), identical(stage_genes(bc, "disp", 50, 10000), "c"))
lt <- data.frame(gene = "a", n_boot = 10000, mean_total_est = 1, mean_total_se = 0.05, mean_total_lo = 0.5, mean_total_hi = 1.5, mean_total_p = 0.001)
mg <- merge_boot_stage(bc, lt, "mean")
stopifnot(mg$mean_total_p[1] == 0.001, mg$n_boot_mean[1] == 10000, mg$n_boot_mean[2] == 1000, mg$disp_total_p[1] == 1, identical(mg[2:3, ], bc[2:3, ]))
lt$mean_total_est <- 2; stopifnot(inherits(try(merge_boot_stage(bc, lt, "mean"), silent = TRUE), "try-error"))

cat("boundary, sign-form, shift and staging checks passed\n")
