###############################################################
### checks_sim.R
### Simulation checks of the methods document (docs/methods.md), one block per check ID of the validation handoff: A1, A2, A4 and C7 (simulation arms),
### B2 (simulation arm), B3, B4, C4, D5, D6, E1, E3 (simulation arm), F1(a), F2, F3 and F4, plus the packaging of the null calibration array (D1 to D4,
### D7, D10) and a reduced local run of it. Every check simulates from the shared generator (sim_hybrid_dataset(), R/functions_checks.R) on the observed
### design when the Section 2 checkpoint exists, and on stand-in offsets with the observed coefficients of variation otherwise (recorded in each result).
### Each result goes to results/checks/<ID>.rds with a row in results/checks/summary.csv; figures go to figures/extra/.
###
### Random numbers: every grid row or dataset draws from its own L'Ecuyer-CMRG stream (check_streams()), so a result does not depend on the number of
### cores. CHECK.SIZE scales the replicate counts (1 = the sizes below). RUN.CHECKS selects the checks to run.
### Runs locally; no cluster output is needed. The calibration array itself runs on the cluster (null_calibration_grid.R, null_calibration.sub).
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/checks_sim.R"
CORES        <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
CHECK.SIZE   <- 1
SEED.CHECKS  <- 2026
RUN.CHECKS   <- c("A2", "A1", "A7", "A9", "A10", "B2", "B3", "B4", "E1", "F1", "F2", "F3", "F4", "C4", "A4", "E3", "D5", "D6", "CAL", "CALMERGE")
set.seed(SEED.SECTION + 20)

## The observed design: offsets, hybrid depth factors and libraries, c-bar and sigma2_c from the Section 2 checkpoint when it exists.
DESIGN.ENV <- new.env()
if (file.exists(ckpt_path(2))) load(ckpt_path(2), envir = DESIGN.ENV)
POOLS <- check_design_pools(DESIGN.ENV)
cat(sprintf("Design pools: %s (Pc %d, Pe %d, hybrid %d cells; c-bar %.4f, sigma2_c %.3g)\n", POOLS$source, length(POOLS$O_Pc), length(POOLS$O_Pe),
            length(POOLS$S), POOLS$cbar, POOLS$var_c))
rm(DESIGN.ENV)
dir.create(CHECK.FIG.DIR, recursive = TRUE, showWarnings = FALSE)

##############################################################################
## A2. Unit tests of the slope rule and of the optimizer (Methods: Estimating phi)
# (a) The adjusted slope at phi = 0 has the closed form 0.5 [sum (y - m)^2 - sum m + sum m^2 / sum m]; a forward finite difference of the adjusted
# log-likelihood at phi = 1e-8 must agree with it in sign. The difference is taken on the stable form of the likelihood; the lgamma form that the fit
# uses is differenced too, to report its rounding error at that phi. (b) The fitted phi must be the maximum of the adjusted likelihood: a 2000-point grid
# over log phi on the search range must not find a better interior point.
if ("A2" %in% RUN.CHECKS) {
  n <- ceiling(1e4 * CHECK.SIZE); S <- check_streams(n, SEED.CHECKS + 2)
  A2 <- do.call(rbind, mclapply(seq_len(n), a2_one, streams = S, mc.cores = CORES))
  ok <- !A2$empty
  check_save("A2", A2, "run",
             headline = c(sprintf("datasets %d (empty %d)", sum(ok), sum(!ok)),
                          sprintf("(a) sign disagreements closed form vs stable finite difference: %d; median rel. error %.2g, max %.2g",
                                  sum(A2$sign_disagree[ok]), median(A2$rel_err[ok]), max(A2$rel_err[ok])),
                          sprintf("lgamma form (.cr_loglik) finite difference at 1e-8: sign disagreements %d, median rel. error %.2g",
                                  sum(sign(A2$slope_fd_lgamma[ok]) != sign(A2$slope_cf[ok])), median(abs(A2$slope_fd_lgamma[ok] - A2$slope_cf[ok]) / abs(A2$slope_cf[ok]))),
                          sprintf("(b) grid finds a better maximum: %d of %d, interior %d (largest log-likelihood gain %.3g; the others at an end of the search range)",
                                  sum(A2$grid_better[ok]), sum(ok), sum(A2$grid_better[ok] & A2$grid_interior[ok]), max(A2$ll_grid[ok] - A2$ll_fit[ok]))),
             mc_error = "exact counts", script = CHECK.SCRIPT)
}

##############################################################################
## A1. Cox-Reid against the unadjusted profile ML and the moment estimators (Methods: Estimating phi; Appendix A)
# The unadjusted profile ML at the ratio mean carries a bias of order 1 / (N Mbar phi) relative to phi; the Cox-Reid term removes its first-order part.
# Relative bias times N Mbar phi is flat in N for the unadjusted fit and near zero after the adjustment.
if ("A1" %in% RUN.CHECKS) {
  G <- expand.grid(N = c(800, 1600, 3200), Mbar = c(0.05, 0.2, 1, 5, 50), phi = c(0, 0.01, 0.1, 0.5, 2))
  S <- check_streams(nrow(G), SEED.CHECKS + 1)
  A1 <- do.call(rbind, mclapply(seq_len(nrow(G)), a1_row, grid = G, nrep = ceiling(500 * CHECK.SIZE), offsets = POOLS$O_Pc, streams = S, mc.cores = CORES))
  A1$bias_over_sd <- A1$bias / A1$sd; A1$bias_z <- A1$bias / A1$bias_se
  SB <- summary_by(A1[A1$phi > 0, ], "estimator", "bias_over_sd", z = "bias_z", keys = c("N", "Mbar", "phi"))
  hl <- sprintf("%s: median bias/SD %.3f, max |bias|/SD %.3f at %s, cells with |bias| > 2 MC SE %d of %d", SB$level, SB$median, SB$max_abs, SB$at, SB$n_z2, SB$n)
  SP <- summary_by(A1[A1$phi > 0 & A1$Mbar * A1$phi >= 0.05 & A1$estimator %in% c("cox_reid", "profile_ml"), ], "estimator", "rel_bias_x_NMphi")
  slope_pat <- sprintf("%s: relative bias x N Mbar phi, median %.2f (range %.2f to %.2f)", SP$level, SP$median, SP$min, SP$max)
  check_save("A1", A1, "run", headline = c(hl, slope_pat, sprintf("offsets %s", POOLS$source)),
             mc_error = sprintf("bias_se and sd_se per cell (%d genes per cell); median bias_se/SD %.3f", ceiling(500 * CHECK.SIZE), median(A1$bias_se / A1$sd, na.rm = TRUE)),
             script = CHECK.SCRIPT)
  fig_pdf("extra/check_A1_relative_bias.pdf", 9, 4)
  par(mfrow = c(1, 2))
  for (k in c("profile_ml", "cox_reid")) {
    x <- A1[A1$estimator == k & A1$phi > 0, ]
    plot(log10(x$N * x$Mbar * x$phi), x$bias / x$phi, pch = 19, col = match(x$N, c(800, 1600, 3200)), xlab = "log10(N Mbar phi)", ylab = "relative bias of phi",
         main = k, ylim = c(-0.5, 0.5)); abline(h = 0, lty = 2)
    legend("topright", legend = paste("N =", c(800, 1600, 3200)), col = 1:3, pch = 19, bty = "n")
  }
  dev.off()
}

##############################################################################
## A7 (simulation arm). Bias in phi from including the gene's own counts in the library that sets its offset (Methods: Sequencing depth offset)
# The data arm (checks_data.R) reports the largest share of a cell's library that one gene takes, to read against this curve.
if ("A7" %in% RUN.CHECKS) {
  G <- expand.grid(share = c(1e-4, 1e-3, 5e-3, 0.01, 0.02, 0.05), phi = c(0.1, 0.5)); S <- check_streams(nrow(G), SEED.CHECKS + 17)
  A7 <- do.call(rbind, mclapply(seq_len(nrow(G)), a7_self_row, grid = G, nrep = ceiling(300 * CHECK.SIZE), lib = POOLS$L, n_cells = 1000, streams = S, mc.cores = CORES))
  check_save("A7_sim", A7, "run", headline = sprintf("library share %g, phi %g: relative bias of phi from self-inclusion %.4f (MC SE %.4f)", A7$share, A7$phi, A7$rel_bias_incl_vs_excl, A7$mc_se),
             mc_error = "MC SE per row", script = CHECK.SCRIPT)
}

##############################################################################
## A9 (Cox-Reid revisit). Does truncation at zero explain the upward bias of the ML estimators in sparse cells?
# If the untruncated moment estimator is unbiased where the truncated estimators are biased, truncation explains the bias.
if ("A9" %in% RUN.CHECKS) {
  G <- expand.grid(N = c(800, 1600, 3200), Mbar = c(0.05, 0.1, 0.2, 0.5), phi = c(0.05, 0.1, 0.2, 0.5), offsets = c("equal", "observed"), stringsAsFactors = FALSE)
  S <- check_streams(nrow(G), SEED.CHECKS + 19)
  A9 <- do.call(rbind, mclapply(seq_len(nrow(G)), a9cr_row, grid = G, nrep = ceiling(500 * CHECK.SIZE), offsets = POOLS$O_Pc, streams = S, mc.cores = CORES))
  A9$bias_z <- A9$bias / A9$bias_se
  SB <- summary_by(A9, "estimator", "bias", z = "bias_z", keys = c("N", "Mbar", "phi", "offsets"))
  check_save("A9", A9, "run", headline = c(sprintf("%s: bias median %.4f, range %.4f to %.4f (largest at %s), |bias| > 2 MC SE in %d of %d cells", SB$level, SB$median, SB$min, SB$max, SB$at, SB$n_z2, SB$n),
                                           sprintf("P(phi-hat = 0), Cox-Reid: %.3f to %.3f", min(A9$p_zero[A9$estimator == "cox_reid"]), max(A9$p_zero[A9$estimator == "cox_reid"]))),
             mc_error = "bias_se per cell", script = CHECK.SCRIPT)
}

##############################################################################
## A10. The covariance of mean and dispersion errors on shared cells (Methods: Separate partitions for each quantity)
# Formula, derived for the moment estimator with equal offsets: Cov[mu_A, phi_B] = (N_AB / (N_A N_B)) mu (E[(a - 1)^3] - 2 phi^2). The check also measures it
# for the Appendix A, unadjusted ML and Cox-Reid estimators and for the observed offset spread, where no formula has been derived. N = 400 per set.
if ("A10" %in% RUN.CHECKS) {
  G <- expand.grid(family = c("gamma", "lognormal", "bimodal"), overlap = c(0, 0.5, 1), offsets = c("equal", "observed"), mu = 2, phi = 0.5, stringsAsFactors = FALSE)
  S <- check_streams(nrow(G), SEED.CHECKS + 20)
  A10 <- do.call(rbind, mclapply(seq_len(nrow(G)), a10_row, grid = G, nrep = ceiling(20000 * CHECK.SIZE), N = 400, offsets = POOLS$O_Pc, streams = S, mc.cores = CORES))
  A10$z_vs_formula <- (A10$N_cov - A10$N_cov_formula) / A10$N_cov_mc_se
  x <- A10[A10$estimator == "moment" & A10$offsets == "equal", ]
  check_save("A10", A10, "run", headline = c(sprintf("moment, equal offsets, %s overlap %.1f: N cov %.3f (MC SE %.3f) against formula %.3f", x$family, x$overlap, x$N_cov, x$N_cov_mc_se, x$N_cov_formula),
                                             sprintf("%s %s offsets %s overlap %.1f: N cov %.3f (MC SE %.3f), cor %.3f", A10$estimator, A10$offsets, A10$family, A10$overlap, A10$N_cov, A10$N_cov_mc_se, A10$cor)[A10$estimator != "moment" | A10$offsets != "equal"]),
             mc_error = sprintf("MC SE per row; cor SE %.4f", 1 / sqrt(ceiling(20000 * CHECK.SIZE) - 3)), script = CHECK.SCRIPT)
}

##############################################################################
## E1. Boundary likelihood-ratio test, upper bound and classes (Methods: Informative estimates)
# (a) At phi = 0 the statistic is an equal mixture of a point mass at zero and a chi-square with one degree of freedom, so the test at 2.71 has size 5%.
# (b) The upper likelihood bound is where the adjusted likelihood ratio from phi-hat reaches the same critical value (2.71, a one-sided 95% bound; the
# document says "upper 95% likelihood bound" without naming the critical value, see mismatches.md); its coverage is P(phi <= bound). (c) The class shares
# when the true phi Mbar is near the near-zero threshold 0.1.
if ("E1" %in% RUN.CHECKS) {
  Ga <- expand.grid(N = c(200, 800, 1600), Mbar = c(0.05, 0.2, 1, 5), cv = c(0.3, 0.7, 1), phi = c(0, 0.01, 0.1, 0.5))
  Gc <- do.call(rbind, lapply(c(0.2, 1, 5), function(m) data.frame(N = 800, Mbar = m, cv = 0.7, phi = c(0.02, 0.05, 0.1, 0.2, 0.5) / m)))
  G <- rbind(Ga, Gc); S <- check_streams(nrow(G), SEED.CHECKS + 3)
  E1 <- do.call(rbind, mclapply(seq_len(nrow(G)), e1_row, grid = G, nrep = ceiling(1000 * CHECK.SIZE), streams = S, mc.cores = CORES))
  E1$part <- c(rep("ab", nrow(Ga)), rep("c", nrow(Gc)))
  z <- E1[E1$part == "ab" & E1$phi == 0, ]; cvg <- E1[E1$part == "ab" & E1$phi > 0, ]; cc <- E1[E1$part == "c", ]
  check_save("E1", E1, "run",
             headline = c(sprintf("(a) size at phi = 0: median %.4f, range %.4f to %.4f over %d cells (N, Mbar, CV)", median(z$reject), min(z$reject), max(z$reject), nrow(z)),
                          sprintf("(b) coverage of the upper bound: median %.4f, range %.4f to %.4f", median(cvg$cover_upper), min(cvg$cover_upper), max(cvg$cover_upper)),
                          sprintf("(c) phi Mbar = 0.1: class2 share %s, class3 share %s at Mbar %s", paste(sprintf("%.3f", cc$class2[abs(cc$phiMbar - 0.1) < 1e-9]), collapse = "/"),
                                  paste(sprintf("%.3f", cc$class3[abs(cc$phiMbar - 0.1) < 1e-9]), collapse = "/"), paste(cc$Mbar[abs(cc$phiMbar - 0.1) < 1e-9], collapse = "/")),
                          sprintf("phi Mbar = 0.02 called class1 (resolved above zero): %s", paste(sprintf("%.3f", cc$class1[abs(cc$phiMbar - 0.02) < 1e-9]), collapse = "/"))),
             mc_error = sprintf("binomial SE of a 5%% rate with %d genes: %.4f", ceiling(1000 * CHECK.SIZE), binom_se(0.05, ceiling(1000 * CHECK.SIZE))), script = CHECK.SCRIPT)
}

##############################################################################
## B2 (simulation arm). The law of total variance estimator of sigma2_c (Methods: Transcriptome share corrections)
if ("B2" %in% RUN.CHECKS) {
  G <- expand.grid(var_c = c(0, 5e-6, 1.9e-5, 1e-4), N = c(800, 1600)); S <- check_streams(nrow(G), SEED.CHECKS + 4)
  B2 <- do.call(rbind, mclapply(seq_len(nrow(G)), b2_row, grid = G, nrep = ceiling(2000 * CHECK.SIZE), lib = POOLS$L, cbar = POOLS$cbar, streams = S, mc.cores = CORES))
  x <- B2[B2$quantity == "sigma2_c (untruncated)", ]
  check_save("B2_sim", B2, "run",
             headline = c(sprintf("sigma2_c = %g, N = %d: bias %.3g (MC SE %.2g), SD %.3g", x$var_c, x$N, x$bias, x$bias_se, x$sd), sprintf("libraries %s", POOLS$source)),
             mc_error = "bias_se per row", script = CHECK.SCRIPT)
}

##############################################################################
## B3. Closed loop of the transcriptome share corrections (Methods: Transcriptome share corrections)
if ("B3" %in% RUN.CHECKS) {
  G <- expand.grid(phi0 = c(0, 0.005, 0.05, 0.5, 5), Mbar = c(0.1, 1, 10), var_c = c(POOLS$var_c, 1e-4)); S <- check_streams(nrow(G), SEED.CHECKS + 5)
  B3 <- do.call(rbind, mclapply(seq_len(nrow(G)), b3_row, grid = G, nrep = ceiling(200 * CHECK.SIZE), pools = POOLS, streams = S, mc.cores = CORES))
  B3$bias_over_sd <- B3$bias / B3$sd
  x <- B3[B3$phi0 > 0, ]
  check_save("B3", B3, "run",
             headline = c(sprintf("phi0 > 0: median |bias|/SD %.3f, max %.3f (phi0=%g Mbar=%g var_c=%g)", median(abs(x$bias_over_sd), na.rm = TRUE), max(abs(x$bias_over_sd), na.rm = TRUE),
                                  x$phi0[which.max(abs(x$bias_over_sd))], x$Mbar[which.max(abs(x$bias_over_sd))], x$var_c[which.max(abs(x$bias_over_sd))]),
                          sprintf("share truncated by the correction: max %.3f", max(B3$truncated_by_correction)),
                          sprintf("relative bias of the corrected means: max |Hc| %.4f, max |He| %.4f", max(abs(B3$mean_rel_bias_Hc)), max(abs(B3$mean_rel_bias_He)))),
             mc_error = sprintf("bias_se per row; median mean_rel_bias_Hc_se %.4f", median(B3$mean_rel_bias_Hc_se)), script = CHECK.SCRIPT)
}

##############################################################################
## B4. NB-ML against the moment estimators when the factor t a is not gamma (Notes for revision: ML fit of non-gamma expression factors)
# pi is the allele's mean share (t = c / pi); sigma2_c 1.9e-5 is the observed value and 1e-3 an inflated one at which t is not negligible.
if ("B4" %in% RUN.CHECKS) {
  G <- expand.grid(family = c("gamma", "lognormal", "bimodal"), Mbar = c(0.5, 2, 10), phi = c(0.05, 0.3, 1), pi = c(0.25, 0.5, 0.75), var_c = c(POOLS$var_c, 1e-3),
                   stringsAsFactors = FALSE)
  S <- check_streams(nrow(G), SEED.CHECKS + 6)
  B4 <- do.call(rbind, mclapply(seq_len(nrow(G)), b4_row, grid = G, nrep = ceiling(300 * CHECK.SIZE), N = 1600, offsets = POOLS$S / 2, streams = S, mc.cores = CORES))
  B4$rel_bias <- B4$bias / B4$truth_var_ta; B4$bias_over_sd <- B4$bias / B4$sd
  B4$family_estimator <- paste(B4$family, B4$estimator)
  SB <- summary_by(B4, "family_estimator", "rel_bias", keys = c("Mbar", "phi", "pi", "var_c")); SD_ <- summary_by(B4, "family_estimator", "bias_over_sd")
  hl <- sprintf("%s: median relative bias %.3f (range %.3f to %.3f, largest at %s), median bias/SD %.2f", SB$level, SB$median, SB$min, SB$max, SB$at, SD_$median)
  check_save("B4", B4, "run", headline = hl, mc_error = sprintf("bias_se per row; median bias_se/truth %.4f", median(B4$bias_se / B4$truth_var_ta)), script = CHECK.SCRIPT)
}

##############################################################################
## B4 addendum. Cox-Reid against the Appendix A estimator at 20,000 cells with equal offsets, and at the observed hybrid depths and cell count
if ("B4" %in% RUN.CHECKS) {
  G <- expand.grid(family = c("gamma", "lognormal", "bimodal"), Mbar = c(1, 5), phi = c(0.3, 1), pi = 0.5, var_c = POOLS$var_c, stringsAsFactors = FALSE)
  S <- check_streams(2 * nrow(G), SEED.CHECKS + 21)
  B4A <- rbind(cbind(design = "20000 cells, equal offsets", do.call(rbind, mclapply(seq_len(nrow(G)), b4_row, grid = G, nrep = ceiling(100 * CHECK.SIZE), N = 20000, offsets = 1,
                                                                                      streams = S[seq_len(nrow(G))], mc.cores = CORES))),
               cbind(design = "observed hybrid depths", do.call(rbind, mclapply(seq_len(nrow(G)), b4_row, grid = G, nrep = ceiling(300 * CHECK.SIZE), N = length(POOLS$S), offsets = POOLS$S / 2,
                                                                                  streams = S[nrow(G) + seq_len(nrow(G))], mc.cores = CORES))))
  B4A$rel_bias <- B4A$bias / B4A$truth_var_ta; B4A$rel_bias_se <- B4A$bias_se / B4A$truth_var_ta
  x <- B4A[B4A$estimator %in% c("cox_reid", "appendixA") & B4A$family != "gamma", ]
  check_save("B4_addendum", B4A, "run", headline = sprintf("%s, %s, mean %g phi %g, %s: relative bias %.3f (MC SE %.3f)", x$design, x$family, x$Mbar, x$phi, x$estimator, x$rel_bias, x$rel_bias_se),
             mc_error = "MC SE per row", script = CHECK.SCRIPT)
}

##############################################################################
## F1(a). Appendix A estimator with estimated weights (Appendix A, final paragraphs)
if ("F1" %in% RUN.CHECKS) {
  G <- expand.grid(N = c(800, 1600), M = c(0.5, 2, 5), phi = c(0.3, 1), cv = c(0.5, 1), counts = c("NB", "Poisson-lognormal"), stringsAsFactors = FALSE)
  S <- check_streams(nrow(G), SEED.CHECKS + 7)
  F1 <- do.call(rbind, mclapply(seq_len(nrow(G)), f1_row, grid = G, nrep = ceiling(500 * CHECK.SIZE), streams = S, mc.cores = CORES))
  x <- F1[F1$estimator == "appendixA", ]
  check_save("F1_sim", F1, "run",
             headline = c(sprintf("Appendix A as documented: median |bias| %.4f, max |bias| %.4f (relative %.3f), cells with |bias| > 2 MC SE %d of %d", median(abs(x$bias)), max(abs(x$bias)),
                                  max(abs(x$bias / x$phi)), sum(abs(x$bias / x$bias_se) > 2), nrow(x)),
                          sprintf("convergence %.4f to %.4f, boundary %.4f to %.4f", min(x$converged), max(x$converged), min(x$boundary), max(x$boundary)),
                          sprintf("code's weighted estimator: max |bias| %.4f", max(abs(F1$bias[F1$estimator == "weighted_code"])))),
             mc_error = sprintf("bias_se per row; median %.4f", median(F1$bias_se)), script = CHECK.SCRIPT)
}

##############################################################################
## F2. Appendix B asymptotics with the current corrections (Appendix B: Variance estimates; summed counts)
# Poisson counts, shifted exponential deviations, true eta and mean rates. Each setting runs F2.REPS independent replicates of 8e6 cells, so the error of
# every limit is read against its Monte Carlo SE (SD over replicates / sqrt(R)) instead of a fixed tolerance. Two settings inflate sigma2_c (0.005, 0.01)
# so the share terms are large enough to resolve; the third uses the observed sigma2_c, where the share terms are far below any resolvable error. The
# sampling SD of each estimator at the observed number of hybrid cells (N_H) says how large an unresolved error could be relative to the data's own noise.
if ("F2" %in% RUN.CHECKS) {
  SET <- data.frame(cbar = c(0.5024, 0.45, 0.5024), var_c = c(0.005, 0.01, POOLS$var_c), s2d = c(0.1, 0.05, 0.15), s2q = c(0.2, 0.3, 0.1), mu_c = c(5, 2, 1), mu_e = c(3, 4, 1))
  F2.REPS <- max(2, ceiling(20 * CHECK.SIZE)); F2.N <- 8e6; N.H <- length(POOLS$S)
  TASKS <- expand.grid(rep = seq_len(F2.REPS), setting = seq_len(nrow(SET))); S <- check_streams(nrow(TASKS), SEED.CHECKS + 8)
  F2.REP <- do.call(rbind, mclapply(seq_len(nrow(TASKS)), f2_rep, tasks = TASKS, settings = SET, n = F2.N, streams = S, mc.cores = CORES))
  F2 <- do.call(rbind, lapply(seq_len(nrow(SET)), f2_summary, settings = SET, reps = F2.REP, n = F2.N, n_obs = N.H))
  check_save("F2", list(summary = F2, replicates = F2.REP), "run",
             headline = c(sprintf("setting %d (sigma2_c %.2g) %s: share term %.2g, error %.2g (MC SE %.1g, z %.1f), sampling SD at N_H = %d: %.2g",
                                  F2$setting, F2$var_c, F2$quantity, F2$share_term, F2$error, F2$mc_se, F2$z, N.H, F2$sd_at_n_obs),
                          sprintf("|z| > 2 in %d of %d", sum(abs(F2$z) > 2), nrow(F2))),
             mc_error = sprintf("%d replicates of %.0e cells per setting; MC SE per row", F2.REPS, F2.N), script = CHECK.SCRIPT)
}

##############################################################################
## F3 and F4. Appendix B at the observed N, and with unequal intrinsic variances between alleles
if (any(c("F3", "F4") %in% RUN.CHECKS)) {
  G <- rbind(expand.grid(omega = c(0.1, 0.25, 0.5, 0.75, 0.9), phi = c(0.1, 0.5), Mbar = c(0.5, 2, 10), q_ratio = 1),
             expand.grid(omega = c(0.25, 0.75), phi = 0.3, Mbar = 2, q_ratio = c(1, 2, 4)))
  S <- check_streams(nrow(G), SEED.CHECKS + 9)
  F34 <- do.call(rbind, mclapply(seq_len(nrow(G)), f34_row, grid = G, nrep = ceiling(300 * CHECK.SIZE), pools = POOLS, streams = S, mc.cores = CORES))
  F34$bias_over_sd <- F34$bias / F34$sd; F34$bias_z <- F34$bias / F34$bias_se
  F3 <- F34[F34$q_ratio == 1 & !(F34$phi == 0.3), ]; F4 <- F34[F34$phi == 0.3, ]
  check_save("F3", F3, "run",
             headline = with(list(SB = summary_by(F3, "quantity", "bias_over_sd", z = "bias_z", keys = c("omega", "phi", "Mbar")), SZ = summary_by(F3, "quantity", "p_zero", keys = c("omega", "phi", "Mbar"))),
                             sprintf("%s: median bias/SD %.3f, max |bias|/SD %.3f at %s, |bias| > 2 MC SE in %d of %d; share truncated at zero up to %.3f at %s",
                                     SB$level, SB$median, SB$max_abs, SB$at, SB$n_z2, SB$n, SZ$max, SZ$at)),
             mc_error = sprintf("bias_se per row (%d genes)", ceiling(300 * CHECK.SIZE)), script = CHECK.SCRIPT)
  x <- F4[F4$quantity == "delta", ]
  check_save("F4", F4, "run", headline = sprintf("delta error at q_ratio %g, omega %g: bias %.4f (MC SE %.4f), truth %.4f", x$q_ratio, x$omega, x$bias, x$bias_se, x$truth),
             mc_error = "bias_se per row", script = CHECK.SCRIPT)
}

##############################################################################
## C4. Internal pilot invariance (Methods: Estimating partition fractions; Wittes and Brittain 1990)
# Under the null the cis and trans tests should keep their size whether f is estimated from bootstrap variances of the same cells or fixed in advance.
if ("C4" %in% RUN.CHECKS) {
  nd <- ceiling(40 * CHECK.SIZE); S <- check_streams(nd, SEED.CHECKS + 10)
  C4 <- do.call(rbind, mclapply(seq_len(nd), c4_dataset, G = 100, mu_grid = c(0.3, 1, 3, 10), phi = 0.3, B_pilot = 50, n_perm = 500, f_fixed = 0.4, pools = POOLS,
                                streams = S, mc.cores = CORES))
  tab <- do.call(rbind, lapply(c("p_cis_est", "p_trans_est", "p_cis_fixed", "p_trans_fixed"), function(k) data.frame(arm = k, size_table(C4[[k]]))))
  check_save("C4", list(pvalues = C4, size = tab), "run",
             headline = c(sprintf("%s: rejection at 0.05 %.4f (SE %.4f), KS p %.3g", tab$arm, tab$reject_0.05, tab$se_0.05, tab$ks_p),
                          sprintf("estimated f: median %.3f (range %.3f to %.3f)", median(unique(C4$f_est)), min(C4$f_est), max(C4$f_est))),
             mc_error = "binomial SE per arm", script = CHECK.SCRIPT)
}

##############################################################################
## A4 and C7 (simulation arms). Mean and dispersion errors are uncorrelated to first order (Methods: Separate partitions for each quantity)
if ("A4" %in% RUN.CHECKS) {
  G <- expand.grid(Mbar = c(0.2, 1, 5), phi = c(0.05, 0.3, 1), omega = c(0.25, 0.75)); S <- check_streams(nrow(G), SEED.CHECKS + 11)
  A4 <- do.call(rbind, mclapply(seq_len(nrow(G)), a4_row, grid = G, nrep = ceiling(1000 * CHECK.SIZE), pools = POOLS, streams = S, mc.cores = CORES))
  A4$z <- A4$cor / A4$se
  SA <- summary_by(A4[A4$pair != "share_phi_zero" & (A4$scale == "raw" | A4$pair == "cis_dmu_dphi"), ], "pair", "cor", z = "z", keys = c("Mbar", "phi", "omega"))
  hl <- sprintf("%s: median cor %.3f, range %.3f to %.3f (largest |cor| at %s), |z| > 2 in %d of %d", SA$level, SA$median, SA$min, SA$max, SA$at, SA$n_z2, SA$n)
  check_save("A4_sim", A4[!grepl("^C7", A4$pair), ], "run", headline = hl[!grepl("^C7", hl)], mc_error = sprintf("Fisher SE 1/sqrt(n-3), n = %d", ceiling(1000 * CHECK.SIZE)), script = CHECK.SCRIPT)
  # The first-order covariance of mean and dispersion errors on shared cells is mu (E[(a - 1)^3] - 2 phi^2) / N: zero for a gamma factor, not otherwise.
  G <- expand.grid(family = c("gamma", "lognormal", "bimodal", "product"), mu = c(1, 5), phi = c(0.3, 1), stringsAsFactors = FALSE); S <- check_streams(nrow(G), SEED.CHECKS + 18)
  A4M <- do.call(rbind, mclapply(seq_len(nrow(G)), a4_moment_row, grid = G, nrep = ceiling(4000 * CHECK.SIZE), N = 2000, streams = S, mc.cores = CORES))
  check_save("A4_theory", A4M, "run", headline = sprintf("%s mu %g phi %g: E[(a-1)^3] %.2f (gamma %.2f), Cov formula %.2e, simulated %.2e (MC SE %.1e); cor moment %.3f, Cox-Reid %.3f",
                                                         A4M$family, A4M$mu, A4M$phi, A4M$third_moment, A4M$third_moment_gamma, A4M$cov_formula, A4M$cov_simulated, A4M$cov_mc_se, A4M$cor_moment, A4M$cor_cox_reid),
             mc_error = sprintf("cor SE %.3f", 1 / sqrt(ceiling(4000 * CHECK.SIZE) - 3)), script = CHECK.SCRIPT)
  check_save("C7_sim", A4[grepl("^C7", A4$pair), ], "run", headline = hl[grepl("^C7", hl)], mc_error = sprintf("Fisher SE 1/sqrt(n-3), n = %d", ceiling(1000 * CHECK.SIZE)), script = CHECK.SCRIPT)
}

##############################################################################
## E3 (simulation arm). Coverage of the 95% percentile intervals with the -Inf mapping (Methods: bootstrap intervals)
if ("E3" %in% RUN.CHECKS) {
  G <- expand.grid(Mbar = c(0.5, 2), phi = c(0.1, 0.5), ratio = c(1, 2), mu_ratio = 1.5); S <- check_streams(nrow(G), SEED.CHECKS + 12)
  E3 <- do.call(rbind, mclapply(seq_len(nrow(G)), e3_row, grid = G, nrep = ceiling(200 * CHECK.SIZE), B = 200, offsets = POOLS$O_Pc, n_cells = 800, streams = S, mc.cores = CORES))
  check_save("E3_sim", E3, "run",
             headline = with(list(SC = summary_by(E3, "quantity", "coverage", keys = c("Mbar", "phi", "ratio")), SZ = summary_by(E3, "quantity", "mean_zero_share")),
                             sprintf("%s: coverage median %.3f, range %.3f to %.3f; mean zero share up to %.4f", SC$level, SC$median, SC$min, SC$max, SZ$max)),
             mc_error = sprintf("binomial SE of 95%% coverage with %d intervals: %.4f", ceiling(200 * CHECK.SIZE), binom_se(0.95, ceiling(200 * CHECK.SIZE))), script = CHECK.SCRIPT)
}

##############################################################################
## D5. Raw dispersion null: projections (a) and the reuse of the log2 draws for total and cis (b) (Methods: Null for the raw dispersion contrasts)
if ("D5" %in% RUN.CHECKS) {
  n <- ceiling(1e4 * CHECK.SIZE); S <- check_streams(n, SEED.CHECKS + 13)
  D5a <- as.data.frame(t(vapply(seq_len(n), d5_case, numeric(5), streams = S)))
  G <- expand.grid(Mbar = c(1, 5), phi = c(0.05, 0.3, 1)); S <- check_streams(nrow(G), SEED.CHECKS + 14)
  D5b <- do.call(rbind, mclapply(seq_len(nrow(G)), d5b_row, grid = G, n_genes = ceiling(150 * CHECK.SIZE), n_sim = 200, offsets = POOLS$O_Pc, n_cells = 800, streams = S, mc.cores = CORES))
  check_save("D5", list(a = D5a, b = D5b), "run",
             headline = c(sprintf("(a) %d cases: raw |sum lambda phi0| > 1e-12 in %d, negative null values %d; log2 projection errors > 1e-12: %d; kinetic: %d", n,
                                  sum(D5a$raw_target_err > 1e-12 & D5a$raw_reachable == 1), sum(D5a$raw_negative), sum(D5a$log2_err > 1e-12, na.rm = TRUE), sum(D5a$kin_err > 1e-12, na.rm = TRUE)),
                          sprintf("(b) %s Mbar %g phi %g: size %.3f (SE %.3f)", D5b$arm, D5b$Mbar, D5b$phi, D5b$reject_0.05, D5b$se_0.05)),
             mc_error = "binomial SE per arm", script = CHECK.SCRIPT)
}

##############################################################################
## D6. Midparent mean null: the null allele means average to the midparent and keep the cis contrast (Methods: Midparent)
# The construction of null_contrasts_one(): fitted allele means mu*_r (MU of HYB.SC and HYB.SE), corrected by .correct_means(), and each allele's expected
# counts M-hat scaled by mu_MP / mu_H. The expected counts are passed through the same corrections, so the null mean contrast must be zero and the cis
# contrast unchanged, for random fitted means and c-bar.
if ("D6" %in% RUN.CHECKS) {
  n <- ceiling(1e4 * CHECK.SIZE); S <- check_streams(n, SEED.CHECKS + 15)
  D6 <- t(vapply(seq_len(n), d6_case, numeric(2), streams = S))
  check_save("D6", D6, "code", headline = sprintf("%d random cases: |mid null contrast| > 1e-12 in %d, |change of cis| > 1e-12 in %d (size: D1 mean_mid in the calibration array)",
                                                  n, sum(D6[, "mid"] > 1e-12), sum(D6[, "cis"] > 1e-12)), mc_error = "exact counts", script = CHECK.SCRIPT)
}

##############################################################################
## Null calibration array (D1 to D4, D7, D10): inputs for cluster/scripts/null_calibration_grid.R and a reduced local run
# The cluster job runs every task of calib_task_grid() with the sizes below; the local run takes a few tasks at reduced size to confirm the pipeline
# runs end to end and to give first numbers, marked "reduced".
if ("CAL" %in% RUN.CHECKS) {
  CAL.GRID <- calib_task_grid(); CAL.G <- 100; CAL.N.SIM <- 1000; CAL.N.PERM <- 1000; CAL.SEED <- SEED.CHECKS + 16
  CAL.STREAMS <- check_streams(nrow(CAL.GRID), CAL.SEED)
  save(CAL.GRID, CAL.G, CAL.N.SIM, CAL.N.PERM, CAL.SEED, CAL.STREAMS, POOLS, file = file.path(INPUT.DIR, "null_calibration_inputs.rda"))
  cat(sprintf("Saved null_calibration_inputs.rda: %d tasks x %d genes, N.SIM %d, N.PERM %d\n", nrow(CAL.GRID), CAL.G, CAL.N.SIM, CAL.N.PERM))
  LOCAL.TASKS <- which(CAL.GRID$mean == 1 & CAL.GRID$phi == 0.3 & CAL.GRID$omega == 0.75 & CAL.GRID$family == "gamma" &
                         (CAL.GRID$design %in% c("equal", "means") | (CAL.GRID$design == "level" & CAL.GRID$par == 1.5) | (CAL.GRID$design == "dispratio" & CAL.GRID$par == 4)))
  CAL.LOCAL <- do.call(rbind, lapply(LOCAL.TASKS, calib_dataset, grid = CAL.GRID, G = ceiling(40 * CHECK.SIZE), n_sim = 200, n_perm = 400, pools = POOLS, streams = CAL.STREAMS, cores = CORES))
  CAL.LOCAL.SIZE <- calib_size_summary(CAL.LOCAL)
  check_save("D1_local", list(pvalues = CAL.LOCAL, size = CAL.LOCAL.SIZE), "reduced",
             headline = c(sprintf("%s %s: reject at 0.05 %.3f (SE %.3f, n %d)", CAL.LOCAL.SIZE$design, CAL.LOCAL.SIZE$test, CAL.LOCAL.SIZE$reject_0.05, CAL.LOCAL.SIZE$se_0.05, CAL.LOCAL.SIZE$n)),
             mc_error = "binomial SE per test; reduced size (N.SIM 200, N.PERM 400)", script = CHECK.SCRIPT)
}

##############################################################################
## Null calibration array results (D1 to D4, D7, D10), read from the cluster outputs
# Run after the array finishes and its null_calibration_output_<k>of<K>.rda files are copied to OUTPUT.DIR. The task keys must cover the grid of the inputs.
if ("CALMERGE" %in% RUN.CHECKS) {
  CAL.FILES <- list.files(OUTPUT.DIR, "^null_calibration_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
  if (!length(CAL.FILES)) {
    for (id in c("D1", "D2", "D3", "D4", "D7", "D10")) check_save(id, NULL, "not run", headline = "no null_calibration_output files in OUTPUT.DIR (cluster job pending)", script = CHECK.SCRIPT)
  } else {
    load(file.path(INPUT.DIR, "null_calibration_inputs.rda"))
    CAL <- list()
    for (f in CAL.FILES) { e <- new.env(); load(f, envir = e); CAL[[f]] <- e$CAL.PART }
    CAL <- do.call(rbind, CAL)
    stopifnot("calibration outputs do not cover every task of the inputs" = setequal(unique(CAL$task), CAL.GRID$task))
    gam <- CAL$family == "gamma"
    D1 <- calib_size_summary(CAL[gam & CAL$design %in% c("equal", "means"), ], strata = c("mean", "phi"))
    D1.MAIN <- D1[!grepl("^w", D1$test), ]
    check_save("D1", D1, "run", headline = c(size_headline(D1.MAIN[D1.MAIN$design == "equal", ], "equal design, NB, by mean and phi"),
                                             size_headline(D1.MAIN[D1.MAIN$design == "means", ], "means design (dispersion nulls), NB, by mean and phi"),
                                             sprintf("KS p < 0.01 in %d of %d tests", sum(D1.MAIN$ks_p < 0.01), nrow(D1.MAIN))),
               mc_error = "binomial SE per test and stratum", script = CHECK.SCRIPT)
    D2 <- calib_size_summary(CAL[CAL$design == "dispratio", ])
    check_save("D2", D2, "run", headline = vapply(sort(unique(D2$par)), function(r) size_headline(D2[D2$par == r, ], paste("phi ratio", r)), ""), mc_error = "binomial SE per test", script = CHECK.SCRIPT)
    D3 <- calib_size_summary(CAL[CAL$design == "level", ], strata = "omega")
    check_save("D3", D3, "run", headline = sprintf("k %g, omega %g, %s: size %.3f (SE %.3f)", D3$par, D3$omega, ifelse(grepl("paired", D3$test), "arm B (pairs kept)", "arm A (pipeline)"), D3$reject_0.05, D3$se_0.05),
               mc_error = "binomial SE per row", script = CHECK.SCRIPT)
    D4 <- list(family = calib_size_summary(CAL[CAL$design == "equal" & !gam, ]), s2d_known = calib_size_summary(CAL[CAL$design == "s2d_known", ]),
               s2d_error = aggregate(cbind(err = sigma_d_hat - sigma_d_true) ~ mean + phi + omega, CAL[CAL$design == "equal" & gam, ], function(x) c(median = median(x), mad = mad(x))),
               private_truncated = aggregate(private_truncated ~ mean + phi + omega, CAL[CAL$design == "equal" & gam, ], mean))
    check_save("D4", D4, "run", headline = c(size_headline(D4$family[D4$family$family == "lognormal" & !grepl("^w", D4$family$test), ], "(a) lognormal truth"),
                                             size_headline(D4$family[D4$family$family == "bimodal" & !grepl("^w", D4$family$test), ], "(a) bimodal truth"),
                                             size_headline(D4$s2d_known[!grepl("^w", D4$s2d_known$test), ], "(b) true sigma2_d in the null"),
                                             sprintf("(c) private variance truncated: share of genes from %.3f to %.3f across strata", min(D4$private_truncated$private_truncated), max(D4$private_truncated$private_truncated))),
               mc_error = "binomial SE per test", script = CHECK.SCRIPT)
    D7 <- D1[grepl("^(bsize|kbal|kbalsup)_", D1$test) & D1$design == "equal", ]
    check_save("D7", D7, "run", headline = c(size_headline(D7[grepl("^kbal_", D7$test), ], "kappa, f_kappa partition (shared cells)"), size_headline(D7[grepl("^kbalsup_", D7$test), ], "kappa, supplementary f_mu and f_phi"),
                                             size_headline(D7[grepl("^bsize_", D7$test), ], "burst size")), mc_error = "binomial SE per test", script = CHECK.SCRIPT)
    D10 <- D1[grepl("^w", D1$test), ]
    check_save("D10", D10, "run", headline = c(size_headline(D10[D10$design == "equal", ], "weighted estimator, equal design"), size_headline(D10[D10$design == "means", ], "weighted estimator, means design")),
               mc_error = "binomial SE per test", script = CHECK.SCRIPT)
  }
}
