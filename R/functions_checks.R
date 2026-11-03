###############################################################
### functions_checks.R
### Functions of the validation checks of the methods document (docs/methods.md): the shared simulation generator, the result writer and the
### per-row units of work that the check scripts apply across their grids. Sourced after R/functions.R by analysis/checks/checks_sim.R,
### analysis/checks/checks_data.R, analysis/checks/checks_text.R and cluster/scripts/null_calibration_grid.R.
### Notation follows the document: parents Pc and Pe, hybrid alleles Hc and He, dispersion phi, burst frequency theta = 1 / phi, burst size zeta,
### kinetic balance kappa, transcriptome share c with mean c-bar and variance sigma2_c, multiplier variances eta.
###############################################################

## ============================================================
## Results, Monte Carlo error and random streams
## ============================================================

## CHECK.DIR holds one <ID>.rds per check and summary.csv; CHECK.FIG.DIR the check figures.
CHECK.DIR     <- file.path("results", "checks")
CHECK.FIG.DIR <- file.path("figures", "extra")

## check_save: writes one check's result to <dir>/<id>.rds and its one-row summary to <dir>/summary.csv, replacing an earlier row of the same ID, so a
## rerun of one check updates its row and keeps the others. headline and mc_error are short text: the numbers the criterion reads and their Monte Carlo
## error. status is "run", "reduced" (run at reduced size), "not run" (inputs missing, with the reason in headline) or "code" (assertions).
check_save <- function(id, result, status, headline, mc_error = "", script, dir = CHECK.DIR) {
  stopifnot(is.character(id), length(id) == 1, is.character(status), length(status) == 1)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(dir, paste0(id, ".rds")))
  row <- data.frame(ID = id, status = status, headline = paste(headline, collapse = "; "), mc_error = paste(mc_error, collapse = "; "),
                    script = script, date = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE)
  f <- file.path(dir, "summary.csv")
  old <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else row[0, ]
  out <- rbind(old[old$ID != id, , drop = FALSE], row)
  write.csv(out[order(out$ID), ], f, row.names = FALSE)
  cat(sprintf("[%s] %s: %s\n", id, status, row$headline))
  invisible(row)
}

## binom_se: the binomial standard error of a rejection rate p estimated from n tests.
binom_se <- function(p, n) sqrt(p * (1 - p) / n)

## mc_stats: bias, SD and RMSE of a set of Monte Carlo estimates against the truth, each with its Monte Carlo standard error, and the share of estimates
## at zero with its binomial SE. The SE of the SD uses the normal-theory sd / sqrt(2 (n - 1)); the SE of the RMSE is the delta-method SE of the root of
## the mean squared error. Non-finite estimates are dropped and counted.
mc_stats <- function(est, truth) {
  n_all <- length(est); est <- est[is.finite(est)]; n <- length(est)
  if (n < 3) return(c(n = n, n_dropped = n_all - n, bias = NA, bias_se = NA, sd = NA, sd_se = NA, rmse = NA, rmse_se = NA, p_zero = NA, p_zero_se = NA))
  e <- est - truth; s <- sd(est); mse <- mean(e^2)
  c(n = n, n_dropped = n_all - n, bias = mean(e), bias_se = s / sqrt(n), sd = s, sd_se = s / sqrt(2 * (n - 1)),
    rmse = sqrt(mse), rmse_se = sd(e^2) / sqrt(n) / (2 * sqrt(mse) + 1e-300), p_zero = mean(est == 0), p_zero_se = binom_se(mean(est == 0), n))
}

## check_streams: n independent L'Ecuyer-CMRG random number streams from one seed, one per unit of work. A unit that installs its own stream
## (with_stream()) gives the same draws however the units are spread over cores or array tasks.
check_streams <- function(n, seed) {
  old_kind <- RNGkind()[1]; had <- exists(".Random.seed", envir = .GlobalEnv); if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  s <- vector("list", n); cur <- get(".Random.seed", envir = .GlobalEnv)
  for (i in seq_len(n)) { s[[i]] <- cur; cur <- parallel::nextRNGStream(cur) }
  RNGkind(old_kind); if (had) assign(".Random.seed", old, envir = .GlobalEnv)
  s
}

## with_stream: evaluates expr with the given L'Ecuyer-CMRG stream as the random state and restores the caller's state afterwards.
with_stream <- function(stream, expr) {
  had <- exists(".Random.seed", envir = .GlobalEnv); if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv))
  assign(".Random.seed", stream, envir = .GlobalEnv)
  expr
}

## ============================================================
## Shared simulation generator (handoff section 3)
## ============================================================
## The observed design: offsets of the parental cells, depth factors S_H and libraries L_H of the hybrid cells, c-bar and sigma2_c. check_design_pools()
## takes them from the Section 2 checkpoint objects when they are loaded and otherwise builds stand-in pools whose offset coefficients of variation are
## those reported for the data (Pc 0.72, Pe 0.87, hybrid 0.44; Methods, Notes for revision). The stand-in median hybrid library is a placeholder.
check_design_pools <- function(env = .GlobalEnv, n_stand = c(Pc = 1000, Pe = 1000, H = 2000), lib_h_median = 2e4) {
  have <- all(vapply(c("EXPO.MIX.SC", "EXPO.MIX.SE", "EXPO.HYB", "HYB.SC", "HYB.SE"), exists, logical(1), envir = env))
  if (have) {
    L <- if (exists("SHARE.CELLS", envir = env)) get("SHARE.CELLS", env)$lib else colSums(get("HYB.SC", env)) + colSums(get("HYB.SE", env))
    return(list(O_Pc = get("EXPO.MIX.SC", env), O_Pe = get("EXPO.MIX.SE", env), S = get("EXPO.HYB", env), L = L,
                cbar = if (exists("C.BAR", envir = env)) get("C.BAR", env) else 0.5024,
                var_c = if (exists("VAR.C", envir = env)) get("VAR.C", env) else 1.9e-5, source = "observed"))
  }
  ln <- function(n, cv) { s <- sqrt(log(1 + cv^2)); exp(rnorm(n, -s^2 / 2, s)) }
  S <- ln(n_stand[["H"]], 0.44) * 2
  list(O_Pc = ln(n_stand[["Pc"]], 0.72), O_Pe = ln(n_stand[["Pe"]], 0.87), S = S, L = round(S * lib_h_median / 2),
       cbar = 0.5024, var_c = 1.9e-5, source = "stand-in")
}

## sim_design: the cells of one simulated experiment, shared by every gene of it. Parental offsets and hybrid cells (S_H with its library L_H) are
## resampled from the pools; each hybrid cell gets a true share c_i (beta, mean cbar, variance var_c) and an observed share c-hat_i = Binomial(L_i, c_i) / L_i.
sim_design <- function(pools, n = c(Pc = length(pools$O_Pc), Pe = length(pools$O_Pe), H = length(pools$S)), cbar = pools$cbar, var_c = pools$var_c) {
  ih <- sample.int(length(pools$S), n[["H"]], replace = TRUE)
  L <- pmax(round(pools$L[ih]), 2); c_i <- sim_draw_share(n[["H"]], cbar, var_c)
  list(O_Pc = sample(pools$O_Pc, n[["Pc"]], replace = TRUE), O_Pe = sample(pools$O_Pe, n[["Pe"]], replace = TRUE),
       S = pools$S[ih], L = L, c = c_i, chat = rbinom(n[["H"]], L, c_i) / L, cbar = cbar, var_c = var_c)
}

## sim_hybrid_dataset: the counts of one gene on the cells of a design (sim_design()). mu and phi are named c(Pc, Pe, Hc, He): per-genome rates and
## expression-factor variances. Parents: Y ~ Poisson(O mu a), a with mean 1 and variance phi. Hybrid alleles: Hc ~ Poisson(S_H c_i mu_Hc d_i q_Hci),
## He ~ Poisson(S_H (1 - c_i) mu_He d_i q_Hei), each allele offset S_H / 2. The shared extrinsic factor d has variance s2d, (1 - omega) times the smaller allele
## dispersion unless given, and each private factor q has the variance (1 + phi_r) / (1 + s2d) - 1 that gives the allele the factor variance phi_r. The
## covariance of the two allele factors is then s2d and omega = 1 - s2d / phi is the intrinsic fraction of the document. family is gamma, lognormal or
## bimodal (sim_pos_factor()), with equal first two moments. q_ratio > 1 gives the Hc allele a private variance q_ratio times that of He at the same mean
## private variance (Appendix B robustness). Returns the counts, the offsets and the true extrinsic and private variances.
sim_hybrid_dataset <- function(design, mu, phi, omega = 0.5, family = "gamma", s2d = NULL, q_ratio = 1) {
  stopifnot(all(c("Pc", "Pe", "Hc", "He") %in% names(mu)), all(c("Pc", "Pe", "Hc", "He") %in% names(phi)), all(phi >= 0), omega >= 0, omega <= 1)
  nH <- length(design$S)
  if (is.null(s2d)) s2d <- (1 - omega) * min(phi[["Hc"]], phi[["He"]])
  vq <- c(Hc = max((1 + phi[["Hc"]]) / (1 + s2d) - 1, 0), He = max((1 + phi[["He"]]) / (1 + s2d) - 1, 0))
  if (q_ratio != 1) { m <- mean(vq); vq <- c(Hc = 2 * m * q_ratio / (1 + q_ratio), He = 2 * m / (1 + q_ratio)) }
  d <- sim_pos_factor(nH, family, 1, s2d)
  y <- list(Pc = rpois(length(design$O_Pc), design$O_Pc * mu[["Pc"]] * sim_pos_factor(length(design$O_Pc), family, 1, phi[["Pc"]])),
            Pe = rpois(length(design$O_Pe), design$O_Pe * mu[["Pe"]] * sim_pos_factor(length(design$O_Pe), family, 1, phi[["Pe"]])),
            Hc = rpois(nH, design$S * design$c * mu[["Hc"]] * d * sim_pos_factor(nH, family, 1, vq[["Hc"]])),
            He = rpois(nH, design$S * (1 - design$c) * mu[["He"]] * d * sim_pos_factor(nH, family, 1, vq[["He"]])))
  list(y = y, O = list(Pc = design$O_Pc, Pe = design$O_Pe, Hc = design$S / 2, He = design$S / 2), s2d = s2d, s2q = vq)
}

## ============================================================
## Estimators compared with the pipeline's
## ============================================================

## .cr_loglik_stable: the Cox-Reid adjusted profile log-likelihood written as in the Methods (Estimating phi), sum_j<y log(1 + j phi) + y log M - log y!
## - (1/phi + y) log(1 + phi M) - 0.5 log sum M / (1 + phi M), with log1p. Every term stays of the size of the counts as phi -> 0, so a finite difference at
## phi = 1e-8 is exact to rounding; the lgamma form of .cr_loglik() differences two numbers of size log Gamma(1 / phi) there.
.cr_loglik_stable <- function(y, m) {
  ymax <- max(y); lfy <- lgamma(y + 1)
  function(phi) {
    if (phi == 0) return(sum(y * log(m) - m - lfy) - 0.5 * log(sum(m)))
    cs <- c(0, cumsum(log1p(phi * seq_len(max(ymax, 1)) - phi)))      # cs[k + 1] = sum_{j < k} log(1 + j phi)
    sum(cs[y + 1] + y * log(m) - lfy - (1 / phi + y) * log1p(phi * m)) - 0.5 * log(sum(m / (1 + phi * m)))
  }
}

## profile_ml_phi: the unadjusted profile maximum likelihood estimate of phi at the ratio mean, the same search as the pipeline without the Cox-Reid term.
## Its slope at phi = 0 is 0.5 [sum (y - m)^2 - sum y], so phi = 0 when that is not positive.
profile_ml_phi <- function(y, O) {
  if (sum(y) == 0) return(NA_real_)
  m <- sum(y) / sum(O) * O
  if (sum((y - m)^2) - sum(y) <= 0) return(0)
  k_ll <- function(lp) { k <- exp(-lp); sum(lgamma(y + k) - lgamma(k) + k * log(k / (k + m)) + y * log(m / (k + m))) }
  exp(optimize(k_ll, -rev(LOG.THETA.RANGE), maximum = TRUE)$maximum)
}

## appendixA_phi: the weighted method of moments of Appendix A as the document states it. With p_i = O_i / sum O, A_i = M_i (1 - p_i),
## B_i = M_i^2 (1 - 2 p_i) + p_i^2 sum M^2 and weights w_i = 1 / (1 + phi M_i)^2, phi-check = sum w ((y - M)^2 - A) / sum w B and
## phi-hat = (1 + 1 / sum y) phi-check / (1 - h phi-check), h = sum p^2, restricted to [0, upper], iterated from phi = 0 until the change is below tol.
## Returns phi (truncated), raw (the last untruncated value), converged, iter and boundary.
appendixA_phi <- function(y, O, maxit = 100, tol = 1e-8, upper = exp(-LOG.THETA.RANGE[[1]])) {
  out <- c(phi = NA_real_, raw = NA_real_, converged = NA_real_, iter = NA_real_, boundary = NA_real_)
  if (sum(y) == 0 || length(y) < 2) return(out)
  p <- O / sum(O); M <- sum(y) * p; A <- M * (1 - p); B <- M^2 * (1 - 2 * p) + p^2 * sum(M^2); h <- sum(p^2); r2 <- (y - M)^2
  phi <- 0; conv <- 0; raw <- NA_real_
  for (it in seq_len(maxit)) {
    w <- 1 / (1 + phi * M)^2
    chk <- sum(w * (r2 - A)) / sum(w * B)
    raw <- if (h * chk < 1) (1 + 1 / sum(y)) * chk / (1 - h * chk) else upper
    new <- min(max(raw, 0), upper)
    if (abs(new - phi) <= tol * (1 + phi)) { phi <- new; conv <- 1; break }
    phi <- new
  }
  c(phi = phi, raw = raw, converged = conv, iter = it, boundary = as.numeric(phi == 0))
}

## ============================================================
## Simulation check rows (checks_sim.R)
## ============================================================

## a1_row: grid row i of A1 (N, Mbar, phi). NREP genes are simulated as NB (gamma factor) counts with offsets resampled from the pool and fitted with the
## Cox-Reid estimator (pipeline), the unadjusted profile ML at the same mean, the code's weighted moment estimator and the Appendix A estimator.
## Returns the Monte Carlo statistics of each against phi, plus the relative bias times N Mbar phi, the scale on which the 1 / (N Mbar phi) bias is flat.
a1_row <- function(i, grid, nrep, offsets, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; est <- matrix(NA_real_, nrep, 4, dimnames = list(NULL, c("cox_reid", "profile_ml", "weighted_code", "appendixA")))
  for (r in seq_len(nrep)) {
    O <- sample(offsets, g$N, replace = TRUE); mu <- g$Mbar / mean(O)
    y <- rpois(g$N, O * mu * sim_pos_factor(g$N, "gamma", 1, g$phi))
    est[r, ] <- c(fit_dispersion(y, O)[["disp"]], profile_ml_phi(y, O), fit_dispersion(y, O, method = "weighted")[["disp"]], appendixA_phi(y, O)[["phi"]])
  }
  do.call(rbind, lapply(colnames(est), function(k) {
    s <- mc_stats(est[, k], g$phi)
    data.frame(g, estimator = k, t(s), rel_bias_x_NMphi = if (g$phi > 0) s[["bias"]] / g$phi * g$N * g$Mbar * g$phi else NA_real_, row.names = NULL)
  }))
})

## a2_one: one random dataset of A2. (a) The closed-form adjusted slope at phi = 0, 0.5 [sum (y - m)^2 - sum m + sum m^2 / sum m], against a forward finite
## difference at h = 1e-8 of the adjusted log-likelihood in its stable form (and of the lgamma form the optimizer uses, for its rounding error).
## (b) The fitted phi (fit_dispersion()) against a 2000-point grid over log phi on the search range: the grid finds a better interior maximum when its best
## adjusted log-likelihood exceeds that at the fitted phi by more than tol, and that point is interior when it is not an end of the grid.
a2_one <- function(i, streams, h = 1e-8, n_grid = 2000, tol = 1e-7) with_stream(streams[[i]], {
  N <- sample(c(30, 100, 300, 800, 1600), 1); Mbar <- exp(runif(1, log(0.02), log(50)))
  phi <- if (runif(1) < 0.3) 0 else exp(runif(1, log(1e-3), log(5))); cv <- runif(1, 0, 1.2)
  s <- sqrt(log(1 + cv^2)); O <- exp(rnorm(N, -s^2 / 2, s))
  y <- rpois(N, O * Mbar * sim_pos_factor(N, "gamma", 1, phi))
  if (sum(y) == 0) return(data.frame(i = i, N = N, Mbar = Mbar, phi = phi, cv = cv, empty = TRUE, slope_cf = NA, slope_fd = NA, slope_fd_lgamma = NA,
                                     sign_disagree = NA, rel_err = NA, phi_hat = NA, ll_fit = NA, ll_grid = NA, phi_grid = NA, grid_better = NA))
  m <- sum(y) / sum(O) * O
  ll <- .cr_loglik_stable(y, m); llg <- .cr_loglik(y, m)
  cf <- 0.5 * (sum((y - m)^2) - sum(m) + sum(m^2) / sum(m))
  fd <- (ll(h) - ll(0)) / h; fdg <- (llg(h) - llg(0)) / h
  fit <- fit_dispersion(y, O); ph <- fit[["disp_app"]]
  lgrid <- seq(-LOG.THETA.RANGE[[2]], -LOG.THETA.RANGE[[1]], length.out = n_grid)
  vg <- vapply(exp(lgrid), ll, numeric(1)); j <- which.max(vg); lf <- ll(ph)
  data.frame(i = i, N = N, Mbar = Mbar, phi = phi, cv = cv, empty = FALSE, slope_cf = cf, slope_fd = fd, slope_fd_lgamma = fdg,
             sign_disagree = sign(cf) != sign(fd) && abs(cf) > 1e-6 * (1 + sum(y)), rel_err = abs(fd - cf) / max(abs(cf), 1e-12),
             phi_hat = ph, ll_fit = lf, ll_grid = vg[j], phi_grid = exp(lgrid[j]), grid_better = vg[j] > lf + tol * max(1, abs(lf)),
             grid_interior = j > 1 && j < n_grid)
})

## e1_row: grid row i of E1 (N, Mbar, offset CV, phi). NREP genes simulated as NB with that phi; dispersion_class_stats() gives the boundary LRT, the upper
## likelihood bound and the class. Returns the rejection rate of the LRT (the size when phi = 0), the coverage of the upper bound (phi <= bound) and the
## class shares, each with its binomial SE.
e1_row <- function(i, grid, nrep, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; s <- sqrt(log(1 + g$cv^2))
  st <- t(vapply(seq_len(nrep), function(r) {
    O <- exp(rnorm(g$N, -s^2 / 2, s)); dispersion_class_stats(rpois(g$N, O * g$Mbar * sim_pos_factor(g$N, "gamma", 1, g$phi)), O)
  }, numeric(5)))
  ok <- is.finite(st[, "lrt"]); n <- sum(ok); rj <- mean(st[ok, "lrt"] > qchisq(1 - 2 * CLASS.ALPHA, 1)); cv_ <- mean(g$phi <= st[ok, "disp_upper"])
  cls <- vapply(1:3, function(k) mean(st[ok, "class"] == k), numeric(1))
  data.frame(g, n = n, reject = rj, reject_se = binom_se(rj, n), cover_upper = cv_, cover_se = binom_se(cv_, n),
             class1 = cls[1], class2 = cls[2], class3 = cls[3], class_se_max = max(binom_se(cls, n)), phiMbar = g$phi * g$Mbar, row.names = NULL)
})

## b2_row: grid row i of B2 (sigma2_c, N). NREP sets of N hybrid cells with libraries resampled from the pool, true shares drawn as beta with variance
## sigma2_c and observed shares as Binomial(L, c) / L, passed to share_terms() (the law of total variance estimator with L - 1). Returns the Monte Carlo
## statistics of sigma2_c-hat and eta_Hc-hat against the truth.
b2_row <- function(i, grid, nrep, lib, cbar, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  est <- t(vapply(seq_len(nrep), function(r) {
    L <- pmax(round(sample(lib, g$N, replace = TRUE)), 2); c_i <- sim_draw_share(g$N, cbar, g$var_c); st <- share_terms(rbinom(g$N, L, c_i) / L, L)
    c(var_c = st$var_obs - st$var_binom, var_c_trunc = st$var_c, eta_c = st$r_c)
  }, numeric(3)))
  rbind(data.frame(g, quantity = "sigma2_c (untruncated)", t(mc_stats(est[, "var_c"], g$var_c)), row.names = NULL),
        data.frame(g, quantity = "sigma2_c (truncated)", t(mc_stats(est[, "var_c_trunc"], g$var_c)), row.names = NULL),
        data.frame(g, quantity = "eta_Hc", t(mc_stats(est[, "eta_c"], g$var_c / cbar^2)), row.names = NULL))
})

## b3_row: grid row i of B3 (phi0, Mbar, sigma2_c). Closed loop of the share correction: both alleles of NREP genes are simulated with true dispersion
## phi0 on a fresh design whose observed shares give eta-hat (share_terms()); each allele is refitted with the pipeline estimator, corrected with
## share_r = eta-hat and truncated at zero; each fitted mean is divided by its abundance factor 2 c-bar-hat or 2 (1 - c-bar-hat). Returns the statistics of
## the corrected phi against phi0, the shares at zero and moved to zero by the correction, and the relative bias of the corrected means.
b3_row <- function(i, grid, nrep, pools, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  res <- t(vapply(seq_len(nrep), function(r) {
    ds <- sim_design(pools, var_c = g$var_c); st <- share_terms(ds$chat, ds$L); mu <- g$Mbar / mean(ds$S / 2)
    x <- sim_hybrid_dataset(ds, mu = c(Pc = mu, Pe = mu, Hc = mu, He = mu), phi = c(Pc = 0, Pe = 0, Hc = g$phi0, He = g$phi0), omega = 1)
    fc <- fit_dispersion(x$y$Hc, x$O$Hc, st$r_c); fe <- fit_dispersion(x$y$He, x$O$He, st$r_e)
    c(phi_c = fc[["disp"]], phi_e = fe[["disp"]], app_c = fc[["disp_app"]], app_e = fe[["disp_app"]],
      mu_c = fc[["mu"]] / (2 * st$cbar) / mu - 1, mu_e = fe[["mu"]] / (2 * (1 - st$cbar)) / mu - 1)
  }, numeric(6)))
  phi <- c(res[, "phi_c"], res[, "phi_e"]); app <- c(res[, "app_c"], res[, "app_e"])
  data.frame(g, t(mc_stats(phi, g$phi0)), truncated_by_correction = mean(phi == 0 & app > 0), truncated_se = binom_se(mean(phi == 0 & app > 0), length(phi)),
             mean_rel_bias_Hc = mean(res[, "mu_c"]), mean_rel_bias_Hc_se = sd(res[, "mu_c"]) / sqrt(nrep),
             mean_rel_bias_He = mean(res[, "mu_e"]), mean_rel_bias_He_se = sd(res[, "mu_e"]) / sqrt(nrep), row.names = NULL)
})

## b4_row: grid row i of B4 (family, Mbar, phi, allele weight pi, sigma2_c). The allele's count factor is t a, with t = c / pi the share multiplier of an
## allele whose mean share is pi (c beta with mean pi and variance sigma2_c) and a of the given family with variance phi. The truth is Var[t a] from 1e6
## latent draws. NREP genes are fitted without the share correction by the Cox-Reid NB fit, the code's weighted moment estimator and the Appendix A
## estimator, which all target Var[t a]. Returns each estimator's statistics against the truth.
b4_row <- function(i, grid, nrep, N, offsets, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  fac <- function(n) (sim_draw_share(n, g$pi, g$var_c) / g$pi) * sim_pos_factor(n, as.character(g$family), 1, g$phi)
  truth <- var(fac(1e6))
  est <- t(vapply(seq_len(nrep), function(r) {
    O <- sample(offsets, N, replace = TRUE); y <- rpois(N, O * g$Mbar / mean(O) * fac(N))
    c(cox_reid = fit_dispersion(y, O)[["disp"]], weighted_code = fit_dispersion(y, O, method = "weighted")[["disp"]], appendixA = appendixA_phi(y, O)[["phi"]])
  }, numeric(3)))
  do.call(rbind, lapply(colnames(est), function(k) data.frame(g, N = N, truth_var_ta = truth, estimator = k, t(mc_stats(est[, k], truth)), row.names = NULL)))
})

## a4_row: grid row i of the simulation arms of A4 and C7 (Mbar, phi, omega). NREP genes on one design: the Hc and He alleles of every hybrid cell share
## the extrinsic factor. Per replicate, mu and phi of Hc are fitted on the same cells (the kinetic balance partition), on half-overlapping halves (mean and
## noise partitions sharing cells) and on disjoint halves; the cis contrasts of the mean and the dispersion use the same cells; C7 pairs the mean of one
## allele with the dispersion of the other. Correlations across replicates on the raw scale (phi at zero included) and on the log2 scale (finite only),
## with the Fisher SE 1 / sqrt(n - 3).
a4_row <- function(i, grid, nrep, pools, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; ds <- sim_design(pools); nH <- length(ds$S); mu <- g$Mbar / mean(ds$S / 2) / (2 * ds$cbar)
  half <- seq_len(nH) <= nH / 2; quart <- seq_len(nH) <= nH / 4 | seq_len(nH) > 3 * nH / 4
  res <- t(vapply(seq_len(nrep), function(r) {
    x <- sim_hybrid_dataset(ds, mu = c(Pc = mu, Pe = mu, Hc = mu, He = mu), phi = c(Pc = g$phi, Pe = g$phi, Hc = g$phi, He = g$phi), omega = g$omega)
    fc <- fit_dispersion(x$y$Hc, x$O$Hc); fe <- fit_dispersion(x$y$He, x$O$He)
    c(mu_c = fc[["mu"]], phi_c = fc[["disp"]], mu_e = fe[["mu"]], phi_e = fe[["disp"]],
      mu_c_A = fit_dispersion(x$y$Hc[half], x$O$Hc[half])[["mu"]], phi_c_B = fit_dispersion(x$y$Hc[quart], x$O$Hc[quart])[["disp"]],
      phi_c_D = fit_dispersion(x$y$Hc[!half], x$O$Hc[!half])[["disp"]])
  }, numeric(7)))
  cr <- function(a, b, log2_scale) {
    if (log2_scale) { a <- log2(a); b <- log2(b) }
    ok <- is.finite(a) & is.finite(b); n <- sum(ok)
    c(cor = if (n > 3) cor(a[ok], b[ok]) else NA_real_, se = if (n > 3) 1 / sqrt(n - 3) else NA_real_, n = n)
  }
  pairs <- list(within_same_cells = c("mu_c", "phi_c"), within_half_overlap = c("mu_c_A", "phi_c_B"), within_disjoint = c("mu_c_A", "phi_c_D"),
                C7_mu_Hc_phi_He = c("mu_c", "phi_e"), C7_mu_He_phi_Hc = c("mu_e", "phi_c"))
  out <- do.call(rbind, lapply(names(pairs), function(p) do.call(rbind, lapply(c(FALSE, TRUE), function(l2)
    data.frame(g, pair = p, scale = if (l2) "log2" else "raw", t(cr(res[, pairs[[p]][1]], res[, pairs[[p]][2]], l2)), row.names = NULL)))))
  dcis <- cr(log2(res[, "mu_c"] / res[, "mu_e"]), log2(res[, "phi_c"] / res[, "phi_e"]), FALSE)
  rbind(out, data.frame(g, pair = "cis_dmu_dphi", scale = "log2", t(dcis), row.names = NULL), data.frame(g, pair = "share_phi_zero", scale = "raw",
        cor = mean(c(res[, "phi_c"], res[, "phi_e"]) == 0), se = binom_se(mean(c(res[, "phi_c"], res[, "phi_e"]) == 0), 2 * nrep), n = 2 * nrep, row.names = NULL))
})

## e3_row: grid row i of the simulation arm of E3 (Mbar, phi, phi ratio). NREP experiments of two parental groups with known log2 contrasts of phi, burst
## size (log2 mu + log2 phi) and kinetic balance (-log2 mu - 2 log2 phi); each gets B bootstrap replicates (cells resampled within group, indices drawn
## first), refitted, with a zero phi mapped to -Inf on the log2 scale; the 95% percentile interval and SE come from .boot_summary(). Returns coverage of the
## truth with its binomial SE, the median interval width and the mean share of replicates at zero.
e3_row <- function(i, grid, nrep, B, offsets, n_cells, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; truth <- c(disp = log2(g$ratio), bsize = log2(g$mu_ratio) + log2(g$ratio), kbal = -log2(g$mu_ratio) - 2 * log2(g$ratio))
  res <- t(vapply(seq_len(nrep), function(r) {
    O1 <- sample(offsets, n_cells, replace = TRUE); O2 <- sample(offsets, n_cells, replace = TRUE); mu1 <- g$Mbar / mean(O1)
    y1 <- rpois(n_cells, O1 * mu1 * g$mu_ratio * sim_pos_factor(n_cells, "gamma", 1, g$phi * g$ratio))
    y2 <- rpois(n_cells, O2 * mu1 * sim_pos_factor(n_cells, "gamma", 1, g$phi))
    idx <- lapply(seq_len(B), function(b) list(sample.int(n_cells, n_cells, TRUE), sample.int(n_cells, n_cells, TRUE)))
    D <- t(vapply(idx, function(ix) {
      f1 <- fit_dispersion(y1[ix[[1]]], O1[ix[[1]]]); f2 <- fit_dispersion(y2[ix[[2]]], O2[ix[[2]]])
      lm_ <- log2(f1[["mu"]] / f2[["mu"]]); ld <- log2(f1[["disp"]]) - log2(f2[["disp"]])
      c(disp = ld, bsize = lm_ + ld, kbal = -lm_ - 2 * log2(f1[["disp"]]) + 2 * log2(f2[["disp"]]))
    }, numeric(3)))
    unlist(lapply(names(truth), function(k) { s <- .boot_summary(D[, k]); setNames(c(s[["lo"]] <= truth[[k]] & truth[[k]] <= s[["hi"]], s[["hi"]] - s[["lo"]], s[["zero"]]),
                                                                                  paste0(k, c("_cover", "_width", "_zero"))) }))
  }, numeric(9)))
  do.call(rbind, lapply(names(truth), function(k) {
    cv_ <- mean(res[, paste0(k, "_cover")], na.rm = TRUE); n <- sum(is.finite(res[, paste0(k, "_cover")]))
    data.frame(g, quantity = k, coverage = cv_, coverage_se = binom_se(cv_, n), median_width = median(res[, paste0(k, "_width")], na.rm = TRUE),
               mean_zero_share = mean(res[, paste0(k, "_zero")], na.rm = TRUE), n = n, row.names = NULL)
  }))
})

## f1_row: grid row i of F1(a) (N, M, phi, offset CV, counts NB or Poisson-lognormal). NREP genes fitted with the Appendix A estimator as documented and with
## the code's weighted estimator. Returns their statistics against phi (the factor variance), with the convergence and boundary rates of Appendix A.
f1_row <- function(i, grid, nrep, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; s <- sqrt(log(1 + g$cv^2)); fam <- if (g$counts == "NB") "gamma" else "lognormal"
  res <- t(vapply(seq_len(nrep), function(r) {
    O <- exp(rnorm(g$N, -s^2 / 2, s)); y <- rpois(g$N, O * g$M * sim_pos_factor(g$N, fam, 1, g$phi)); a <- appendixA_phi(y, O)
    c(appendixA = a[["phi"]], conv = a[["converged"]], bnd = a[["boundary"]], weighted_code = fit_dispersion(y, O, method = "weighted")[["disp"]])
  }, numeric(4)))
  rbind(data.frame(g, estimator = "appendixA", t(mc_stats(res[, "appendixA"], g$phi)), converged = mean(res[, "conv"]), boundary = mean(res[, "bnd"]), row.names = NULL),
        data.frame(g, estimator = "weighted_code", t(mc_stats(res[, "weighted_code"], g$phi)), converged = NA, boundary = NA, row.names = NULL))
})

## f2_setting: F2, the large-N limits of Appendix B. n cells with depth S_H = 1, Poisson counts, additive expression factors a_r = 1 + d + q_r with shifted
## exponential deviations (d = sd_d (E - 1), E ~ Exp(1), and likewise q, so each has mean 0 and the stated variance), true shares c ~ beta(c-bar, sigma2_c)
## and the true mean rates and eta in every estimator. Compares the star estimators with their stated limits, the corrected estimators with sigma2_d and
## sigma2_q, and the variance of the summed factor a*_H with the formula for phi*_Sigma; the Monte Carlo SE of each star estimator is the SD of its
## per-cell term over sqrt(n). Cells are processed in chunks to bound memory.
f2_setting <- function(i, settings, n, chunk = 1e6, streams) with_stream(streams[[i]], {
  s <- settings[i, ]; cb <- s$cbar; v <- s$var_c
  eta_c <- v / cb^2; eta_e <- v / (1 - cb)^2; eta_ce <- v / (cb * (1 - cb)); eta_bar <- (eta_c + eta_e) / 2
  acc <- c(dd = 0, qq = 0, a1 = 0, a2 = 0, n = 0, dd2 = 0, qq2 = 0)
  for (k in seq_len(ceiling(n / chunk))) {
    m <- min(chunk, n - (k - 1) * chunk)
    c_i <- sim_draw_share(m, cb, v); d <- sqrt(s$s2d) * (rexp(m) - 1)
    ac <- 1 + d + sqrt(s$s2q) * (rexp(m) - 1); ae <- 1 + d + sqrt(s$s2q) * (rexp(m) - 1)
    yc <- rpois(m, c_i * s$mu_c * ac); ye <- rpois(m, (1 - c_i) * s$mu_e * ae)
    xc <- yc / (cb * s$mu_c); xe <- ye / ((1 - cb) * s$mu_e)                   # scaled counts with the true rates (O mu* = S c-bar mu)
    pi_ <- cb * s$mu_c / (cb * s$mu_c + (1 - cb) * s$mu_e)
    astar <- pi_ * (c_i / cb) * ac + (1 - pi_) * ((1 - c_i) / (1 - cb)) * ae
    td <- xc * xe; tq <- 0.5 * (xc - xe)^2 - 0.5 * (yc / (cb * s$mu_c)^2 + ye / ((1 - cb) * s$mu_e)^2)
    acc <- acc + c(sum(td), sum(tq), sum(astar), sum(astar^2), m, sum(td^2), sum(tq^2))
  }
  d_star <- acc[["dd"]] / acc[["n"]] - 1; q_star <- acc[["qq"]] / acc[["n"]]; var_astar <- acc[["a2"]] / acc[["n"]] - (acc[["a1"]] / acc[["n"]])^2
  phiH <- s$s2d + s$s2q
  lim_d <- s$s2d - eta_ce * (1 + s$s2d)
  lim_q <- s$s2q + 0.5 * (1 + phiH) * v * (1 / cb^2 + 1 / (1 - cb)^2) + eta_ce * (1 + s$s2d)
  sd2 <- (1 + d_star) / (1 - eta_ce) - 1; sq2 <- (q_star - (eta_bar + eta_ce) * (1 + sd2)) / (1 + eta_bar)
  pi_ <- cb * s$mu_c / (cb * s$mu_c + (1 - cb) * s$mu_e)
  phiSig <- s$s2d + (pi_^2 + (1 - pi_)^2) * s$s2q
  form <- phiSig + (1 + s$s2d) * (pi_^2 * eta_c + (1 - pi_)^2 * eta_e - 2 * pi_ * (1 - pi_) * eta_ce) + s$s2q * (pi_^2 * eta_c + (1 - pi_)^2 * eta_e)
  se_d <- sqrt((acc[["dd2"]] / acc[["n"]] - (acc[["dd"]] / acc[["n"]])^2) / acc[["n"]]); se_q <- sqrt((acc[["qq2"]] / acc[["n"]] - q_star^2) / acc[["n"]])
  data.frame(s, n = n, d_star = d_star, d_star_limit = lim_d, d_star_err = d_star - lim_d, d_star_mc_se = se_d, q_star_mc_se = se_q, q_star = q_star, q_star_limit = lim_q, q_star_err = q_star - lim_q,
             sigma2_d_hat = sd2, sigma2_d_err = sd2 - s$s2d, sigma2_q_hat = sq2, sigma2_q_err = sq2 - s$s2q,
             var_astar = var_astar, var_astar_formula = form, var_astar_err = var_astar - form, row.names = NULL)
})

## f34_row: grid row i of F3 and F4 (omega, phi, Mbar, intrinsic variance ratio q_ratio). NREP genes on one design with the observed N, fitted means and the
## pipeline's ploidy_components() (share terms from the simulated observed shares). The truth: sigma2_d = s2d (the allele covariance), the mean
## private variance s2q in the additive form of Appendix B (Var[a] - Cov = (1 + s2d) times the multiplicative private variance), and
## delta = log2(phi_H / phi_Sigma) with phi_H = s2d + mean s2q and phi_Sigma = s2d + pi^2 s2q_c + (1 - pi)^2 s2q_e (equal to the document's delta when
## q_ratio = 1). Returns bias, SD and the truncation shares of each estimate.
f34_row <- function(i, grid, nrep, pools, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; ds <- sim_design(pools); st <- share_terms(ds$chat, ds$L); mu <- g$Mbar / mean(ds$S / 2) / (2 * ds$cbar)
  sims <- lapply(seq_len(nrep), function(r) sim_hybrid_dataset(ds, mu = c(Pc = mu, Pe = mu, Hc = mu, He = mu), phi = c(Pc = g$phi, Pe = g$phi, Hc = g$phi, He = g$phi),
                                                                omega = g$omega, q_ratio = g$q_ratio))
  yc <- do.call(rbind, lapply(sims, function(x) x$y$Hc)); ye <- do.call(rbind, lapply(sims, function(x) x$y$He)); rownames(yc) <- rownames(ye) <- paste0("r", seq_len(nrep))
  pc <- ploidy_components(yc, ye, ds$S / 2, st)
  s2d <- sims[[1]]$s2d; s2q <- (1 + s2d) * sims[[1]]$s2q; pi_ <- ds$cbar / (ds$cbar + (1 - ds$cbar))
  phiH <- s2d + mean(s2q); phiS <- s2d + pi_^2 * s2q[["Hc"]] + (1 - pi_)^2 * s2q[["He"]]
  truth <- c(sigma_d = s2d, sigma_q = mean(s2q), delta = log2(phiH / phiS), delta_raw = phiH - phiS)
  do.call(rbind, lapply(names(truth), function(k) data.frame(g, quantity = k, truth = truth[[k]], t(mc_stats(pc[[k]], truth[[k]])), row.names = NULL)))
})

## c4_dataset: one experiment of C4 under the null (no cis and no trans effect: every per-genome rate equal and every dispersion equal). G genes on one
## design; the cis fraction f is estimated as in analysis.R 2.1 from B_pilot bootstrap replicates of all cells (Var of the hybrid log2 allele ratio and of
## the parental log2 ratio, pooled by medians across genes) or fixed in advance (f_fixed). Each gene's log2 mean cis and trans contrasts are tested with
## the pipeline's permutation nulls (cis: rescaled allele counts swapped within cells; trans: each parent pooled with its hybrid allele, offsets carrying
## the abundance factor and k-hat), n_perm permutations, add-one p-values. Returns the p-values of both arms and the estimated f.
c4_dataset <- function(i, G, mu_grid, phi, B_pilot, n_perm, f_fixed, pools, streams) with_stream(streams[[i]], {
  ds <- sim_design(pools); nH <- length(ds$S); cb <- ds$cbar
  genes <- lapply(seq_len(G), function(j) { m <- sample(mu_grid, 1); sim_hybrid_dataset(ds, mu = c(Pc = m, Pe = m, Hc = m, He = m), phi = c(Pc = phi, Pe = phi, Hc = phi, He = phi)) })
  l2 <- function(y, o) log2(sum(y) / sum(o))
  ## pilot: bootstrap variances of the hybrid allele log2 ratio (all hybrid cells) and of the parental log2 ratio
  ix <- lapply(seq_len(B_pilot), function(b) list(h = sample.int(nH, nH, TRUE), c = sample.int(length(ds$O_Pc), length(ds$O_Pc), TRUE), e = sample.int(length(ds$O_Pe), length(ds$O_Pe), TRUE)))
  AB <- t(vapply(genes, function(x) {
    hb <- vapply(ix, function(k) l2(x$y$Hc[k$h], x$O$Hc[k$h]) - l2(x$y$He[k$h], x$O$He[k$h]), numeric(1))
    pb <- vapply(ix, function(k) l2(x$y$Pc[k$c], x$O$Pc[k$c]) - l2(x$y$Pe[k$e], x$O$Pe[k$e]), numeric(1))
    c(A = nH * var(hb[is.finite(hb)]), B = var(pb[is.finite(pb)]))
  }, numeric(2)))
  r <- median(AB[, "B"], na.rm = TRUE) * nH / median(AB[, "A"], na.rm = TRUE); f_est <- fstar_from_r(r)
  pv <- function(f) {
    sp <- split_indices_by_depth(nH, f); cc <- sp$c; tt <- sp$t
    vapply(genes, function(x) {
      yc <- x$y$Hc[cc] / (2 * cb); ye <- x$y$He[cc] / (2 * (1 - cb)); oc <- x$O$Hc[cc]
      obs_cis <- l2(yc, oc) - l2(ye, oc)
      sw <- matrix(runif(length(cc) * n_perm) < 0.5, length(cc))
      null_cis <- log2(colSums(ifelse(sw, ye, yc)) / sum(oc)) - log2(colSums(ifelse(sw, yc, ye)) / sum(oc))
      mPc <- sum(x$y$Pc) / sum(x$O$Pc); mPe <- sum(x$y$Pe) / sum(x$O$Pe)
      mHc <- sum(x$y$Hc[tt]) / sum(x$O$Hc[tt]) / (2 * cb); mHe <- sum(x$y$He[tt]) / sum(x$O$He[tt]) / (2 * (1 - cb))
      k <- sqrt(mHc * mHe / (mPc * mPe))
      otc <- k * cb * 2 * x$O$Hc[tt]; ote <- k * (1 - cb) * 2 * x$O$He[tt]
      obs_tr <- (log2(mPc) - log2(mPe)) - (log2(mHc) - log2(mHe))
      pc_y <- c(x$y$Pc, x$y$Hc[tt]); pc_o <- c(x$O$Pc, otc); pe_y <- c(x$y$Pe, x$y$He[tt]); pe_o <- c(x$O$Pe, ote)
      nPc <- length(x$y$Pc); nPe <- length(x$y$Pe)
      null_tr <- vapply(seq_len(n_perm), function(b) { a <- sample.int(length(pc_y), nPc); e <- sample.int(length(pe_y), nPe)
        (l2(pc_y[a], pc_o[a]) - l2(pe_y[e], pe_o[e])) - (l2(pc_y[-a], pc_o[-a]) - l2(pe_y[-e], pe_o[-e])) }, numeric(1))
      c(cis = perm_pval(obs_cis, null_cis), trans = perm_pval(obs_tr, null_tr))
    }, numeric(2))
  }
  pe <- pv(f_est); pf <- pv(f_fixed)
  data.frame(dataset = i, gene = seq_len(G), f_est = f_est, f_fixed = f_fixed, p_cis_est = pe["cis", ], p_trans_est = pe["trans", ], p_cis_fixed = pf["cis", ], p_trans_fixed = pf["trans", ])
})

## d5_case: D5(a), one random case of the null projections: null_phi_raw() must zero sum(lambda phi0) whenever that is reachable with non-negative
## values and return no negative value; null_phi_log2() and null_phi_kinetic() must zero their contrasts exactly on the log2 scale.
d5_case <- function(i, streams) with_stream(streams[[i]], {
  md <- sample(names(.MODE.WEIGHTS), 1)
  lam <- .MODE.WEIGHTS[[md]]; ph <- exp(rnorm(length(lam), log(0.2), 1.5)); if (runif(1) < 0.2) ph[sample.int(length(ph), 1)] <- 0
  r <- null_phi_raw(ph, lam); pos_reach <- sum(lam > 0) > 0 && sum(lam < 0) > 0
  out <- c(raw_target_err = abs(sum(lam * r)), raw_negative = sum(r < 0), raw_reachable = pos_reach, log2_err = NA_real_, kin_err = NA_real_)
  if (all(ph > 0)) {
    out[["log2_err"]] <- abs(sum(lam * log2(null_phi_log2(ph, lam))))
    dmu <- rnorm(1); k <- null_phi_kinetic(ph, lam, dmu, -1, -2); out[["kin_err"]] <- abs(-dmu - 2 * sum(lam * log2(k)))
  }
  out
})

## d5b_row: D5(b), the raw total-contrast null of two parental groups whose true dispersions are equal (phi) and whose fitted values differ by noise. Arm
## "reuse" simulates the null at the geometric mean of the two fits (the log2 null, reused by the pipeline for total and cis), arm "separate" at their
## arithmetic mean (the raw projection). n_sim draws per gene, add-one two-sided p-values of the observed raw difference. Returns the size table of each arm.
d5b_row <- function(i, grid, n_genes, n_sim, offsets, n_cells, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  P <- t(vapply(seq_len(n_genes), function(j) {
    O1 <- sample(offsets, n_cells, TRUE); O2 <- sample(offsets, n_cells, TRUE); mu <- g$Mbar / mean(O1)
    y1 <- rpois(n_cells, O1 * mu * sim_pos_factor(n_cells, "gamma", 1, g$phi)); y2 <- rpois(n_cells, O2 * mu * sim_pos_factor(n_cells, "gamma", 1, g$phi))
    f1 <- fit_dispersion(y1, O1); f2 <- fit_dispersion(y2, O2); obs <- f1[["disp"]] - f2[["disp"]]
    if (!(f1[["disp"]] > 0 && f2[["disp"]] > 0)) return(c(p_reuse = NA_real_, p_separate = NA_real_))
    nul <- function(p0) vapply(seq_len(n_sim), function(b) fit_dispersion(rpois(n_cells, O1 * f1[["mu"]] * rgamma1(n_cells, p0)), O1)[["disp"]] -
                                                             fit_dispersion(rpois(n_cells, O2 * f2[["mu"]] * rgamma1(n_cells, p0)), O2)[["disp"]], numeric(1))
    c(p_reuse = perm_pval(obs, nul(sqrt(f1[["disp"]] * f2[["disp"]]))), p_separate = perm_pval(obs, nul((f1[["disp"]] + f2[["disp"]]) / 2)))
  }, numeric(2)))
  rbind(data.frame(g, arm = "reuse (geometric centre)", size_table(P[, "p_reuse"]), row.names = NULL),
        data.frame(g, arm = "separate (arithmetic centre)", size_table(P[, "p_separate"]), row.names = NULL))
})

## d6_case: D6, one random case of the midparent mean null construction of null_contrasts_one(): fitted allele means mu*_r (MU of HYB.SC and HYB.SE) and
## parental means, corrected by .correct_means(), and each allele's expected counts scaled by mu_MP / mu_H. Passed through the same corrections, the null
## means must give a zero midparent contrast and leave the cis contrast unchanged. Returns both absolute errors.
d6_case <- function(i, streams) with_stream(streams[[i]], {
  cb <- runif(1, 0.3, 0.7); gv <- list(MIX.SC = exp(rnorm(1)), MIX.SE = exp(rnorm(1)), HYB.SC = exp(rnorm(1)), HYB.SE = exp(rnorm(1)))
  cm <- .correct_means(gv, cb); sc_f <- 0.5 * (cm$MIX.SC + cm$MIX.SE) / (0.5 * (cm$HYB.SC + cm$HYB.SE))
  nul <- .correct_means(list(MIX.SC = gv$MIX.SC, MIX.SE = gv$MIX.SE, HYB.SC = gv$HYB.SC * sc_f, HYB.SE = gv$HYB.SE * sc_f), cb)
  c(mid = abs(unname(.contrast_value("mid", lapply(nul, function(x) c(MU = x)), "MU"))), cis = abs(log2(nul$HYB.SC / nul$HYB.SE) - log2(cm$HYB.SC / cm$HYB.SE)))
})

## summary_by: one row per level of `by` with the median, minimum and maximum of column col, the number of rows whose column z exceeds 2 in absolute value
## (a Monte Carlo z-score, NA when z is not given) and the grid values of the row with the largest absolute col. Used for the headline of each check.
summary_by <- function(df, by, col, z = NULL, keys = NULL) {
  do.call(rbind, lapply(split(df, df[[by]]), function(x) {
    v <- x[[col]]; j <- which.max(abs(v))
    data.frame(level = x[[by]][1], median = median(v, na.rm = TRUE), min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE), max_abs = abs(v[j]),
               n_z2 = if (is.null(z)) NA_integer_ else sum(abs(x[[z]]) > 2, na.rm = TRUE), n = sum(is.finite(v)),
               at = if (is.null(keys)) "" else paste(paste0(keys, "=", unlist(x[j, keys])), collapse = " "), stringsAsFactors = FALSE)
  }))
}

## size_table: rejection rates of a vector of null p-values at each level, with binomial SEs and the Kolmogorov-Smirnov p-value of uniformity.
size_table <- function(p, levels = c(0.05, 0.01, 0.001)) {
  p <- p[is.finite(p)]; n <- length(p)
  if (n == 0) return(data.frame(n = 0))
  r <- vapply(levels, function(a) mean(p <= a), numeric(1))
  data.frame(n = n, t(setNames(r, paste0("reject_", levels))), t(setNames(binom_se(r, n), paste0("se_", levels))),
             ks_p = suppressWarnings(ks.test(p, "punif")$p.value))
}

## ============================================================
## Null calibration datasets (D1 to D4, D7, D10; cluster/scripts/null_calibration_grid.R)
## ============================================================
## calib_task_grid: the tasks of the calibration array. The common grid is mean {0.2, 1, 5} x phi {0.05, 0.3, 1} x omega {0.25, 0.75}; each design sets
## which nulls are true:
##   equal      every group shares one per-genome mean and one phi (every null of D1, D7 and D10 is true);
##   means      group means differ, one phi (every dispersion null is true, D1 dispersion arms);
##   dispratio  equal means, phi_Pe = phi_He = r phi_Pc = r phi_Hc (D2: the mean permutation nulls with unequal dispersions);
##   level      hybrid rates k times the parental rates, parental ratio 2 (D3: the trans mean null, k-hat held fixed; also run at omega = 0.1);
##   lognormal, bimodal   the equal design with non-gamma factors (D4a);
##   s2d_known  the equal design with the true sigma2_d passed to the null instead of the per-gene estimate (D4b).
calib_task_grid <- function() {
  base <- expand.grid(mean = c(0.2, 1, 5), phi = c(0.05, 0.3, 1), omega = c(0.25, 0.75))
  add <- function(design, par = NA_real_, family = "gamma", omega = NULL) { b <- base; if (!is.null(omega)) { b$omega <- omega; b <- unique(b) }; cbind(b, design = design, par = par, family = family) }
  out <- rbind(add("equal"), add("means"), add("dispratio", 2), add("dispratio", 4), add("dispratio", 8),
               do.call(rbind, lapply(c(0.5, 0.75, 1, 1.5, 2), function(k) rbind(add("level", k), add("level", k, omega = 0.1)))),
               add("equal", family = "lognormal"), add("equal", family = "bimodal"), add("s2d_known"))
  out$family <- as.character(out$family); out$design <- as.character(out$design)
  out$task <- seq_len(nrow(out)); out
}

## calib_dataset: one task of the calibration array. G synthetic genes share one design (observed group sizes and offsets); each gene's per-genome mean is
## the grid mean times a lognormal spread (SD 0.3 on the log scale) and its dispersion the grid phi times the same spread, then the design's effects are
## applied. The datasets mirror analysis.R Section 2: share terms from the observed shares, ploidy_components() for sigma2_d and the summed-count
## correction, three hybrid partitions at f, the NB fits of every dataset with the share corrections, permutation draws, and null_contrasts_one() with both
## estimators. D3's pairing-preserving arm B of the trans mean null (perm_trans_pool_paired()) is evaluated on the same genes. D4c's private-variance
## truncation is recorded per gene: the null simulation of an allele sets its private variance to zero when 1 + phi(1 + eta) + eta < 1 + sigma2_d.
calib_dataset <- function(task, grid, G, n_sim, n_perm, pools, f = 0.4, streams, cores = 1) with_stream(streams[[task]], {
  tk <- grid[task, ]; ds <- sim_design(pools); nH <- length(ds$S)
  spread <- exp(rnorm(G, 0, 0.3)); spread2 <- exp(rnorm(G, 0, 0.3))
  mus <- lapply(seq_len(G), function(j) {
    m <- tk$mean * spread[j] / mean(ds$O_Pc)                                    # per-genome rate: tk$mean counts per parental cell on average
    switch(tk$design, means = m * c(Pc = 1, Pe = 0.5, Hc = 1.6, He = 0.8), level = m * c(Pc = 1, Pe = 0.5, Hc = tk$par, He = tk$par * 0.5), c(Pc = m, Pe = m, Hc = m, He = m))
  })
  phis <- lapply(seq_len(G), function(j) { p <- tk$phi * spread2[j]; if (tk$design == "dispratio") c(Pc = p, Pe = p * tk$par, Hc = p, He = p * tk$par) else c(Pc = p, Pe = p, Hc = p, He = p) })
  sims <- lapply(seq_len(G), function(j) sim_hybrid_dataset(ds, mus[[j]], phis[[j]], omega = tk$omega, family = tk$family))
  gn <- sprintf("g%04d", seq_len(G))
  mat <- function(k) { x <- do.call(rbind, lapply(sims, function(s) s$y[[k]])); rownames(x) <- gn; x }
  Hsc <- mat("Hc"); Hse <- mat("He")
  st <- share_terms(ds$chat, ds$L); SR <- share_r_groups(st$r_c, st$r_e); eH <- ds$S / 2
  pc <- ploidy_components(Hsc, Hse, eH, st)
  ix <- list(MEAN = split_indices_by_depth(nH, f), DISP = split_indices_by_depth(nH, f), KBAL = split_indices_by_depth(nH, f))
  M <- list(MIX.SC = mat("Pc"), MIX.SE = mat("Pe"), HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
  E <- list(MIX.SC = ds$O_Pc, MIX.SE = ds$O_Pe, HYB.SC = eH, HYB.SE = eH, HYB.COMB = ds$S)
  for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
    nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K"))
    M[[nm]] <- (if (al == "SC") Hsc else Hse)[, ix[[sp]][[side]], drop = FALSE]; E[[nm]] <- eH[ix[[sp]][[side]]]
  }
  FITS <- lapply(setNames(names(M), names(M)), function(k) {
    fr <- fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]], share_sub = if (k == "HYB.COMB") pc$share_sub else 0)
    data.frame(MU = fr$MU, DISP = fr$DISP, CLASS = fr$CLASS, row.names = rownames(fr))
  })
  s2d_true <- vapply(sims, `[[`, numeric(1), "s2d")
  share_obs <- list(cbar = st$cbar, r = SR, sigma_d = setNames(if (tk$design == "s2d_known") s2d_true else pc$sigma_d, gn))
  perms <- lapply(seq_len(n_perm), perm_label_draw, hyc = ix$MEAN$c, hyt = ix$MEAN$t, nHYB = nH, nSC = length(ds$O_Pc), nSE = length(ds$O_Pe))
  run <- function(g) null_contrasts_one(g, expos = E, fits = FITS, mats = M, perms = perms, share_obs = share_obs, n_sim = n_sim, seed = 1e6 * task, methods = c("NB", "weighted"))
  PR <- do.call(rbind, if (cores > 1) parallel::mclapply(gn, run, mc.cores = cores, mc.preschedule = FALSE) else lapply(gn, run))
  ## D3 arm B: the trans mean null with hybrid cells kept in pairs (the same hybrid cells move to the pseudo-parent side of both pools)
  pB <- vapply(gn, trans_null_paired_p, numeric(1), M = M, E = E, cb = st$cbar, n_perm = n_perm)
  ## D4c: genes whose null simulation of a cis-partition allele truncates the private variance at zero
  vq_trunc <- vapply(gn, function(g) any(1 + FITS$HYC.SC.N[g, "DISP"] * (1 + SR[["HYC.SC.N"]]) + SR[["HYC.SC.N"]] < 1 + share_obs$sigma_d[[g]],
                                          1 + FITS$HYC.SE.N[g, "DISP"] * (1 + SR[["HYC.SE.N"]]) + SR[["HYC.SE.N"]] < 1 + share_obs$sigma_d[[g]], na.rm = TRUE), logical(1))
  cbind(tk[rep(1, G), ], gene = gn, mu_gene = vapply(mus, `[[`, numeric(1), "Pc"), phi_gene = vapply(phis, `[[`, numeric(1), "Pc"),
        sigma_d_true = s2d_true, sigma_d_hat = pc$sigma_d, private_truncated = vq_trunc, mean_trans_p_paired = pB, PR, row.names = NULL)
})

## trans_null_paired_p: D3 arm B for one gene, the pairing-preserving trans mean null. The pools are those of the pipeline (each parent with the hybrid
## allele of its genome from the trans cells, offsets carrying k-hat and the abundance factor); perm_trans_pool_paired() sends the same hybrid cells to the
## pseudo-parent side of both pools, so the two alleles of a hybrid cell stay together and keep their shared extrinsic noise. Add-one two-sided p-value.
trans_null_paired_p <- function(g, M, E, cb, n_perm) {
  y_pc <- M$MIX.SC[g, ]; y_pe <- M$MIX.SE[g, ]; y_hc <- M$HYT.SC[g, ]; y_he <- M$HYT.SE[g, ]
  l2 <- function(y, o) { s <- sum(y); if (s > 0) log2(s / sum(o)) else NA_real_ }
  mPc <- sum(y_pc) / sum(E$MIX.SC); mPe <- sum(y_pe) / sum(E$MIX.SE)
  mHc <- sum(y_hc) / sum(E$HYT.SC) / (2 * cb); mHe <- sum(y_he) / sum(E$HYT.SE) / (2 * (1 - cb))
  if (!all(is.finite(c(mPc, mPe, mHc, mHe)) & c(mPc, mPe, mHc, mHe) > 0)) return(NA_real_)
  k <- sqrt(mHc * mHe / (mPc * mPe))
  pc_y <- c(y_pc, y_hc); pc_o <- c(E$MIX.SC, k * cb * 2 * E$HYT.SC); pe_y <- c(y_pe, y_he); pe_o <- c(E$MIX.SE, k * (1 - cb) * 2 * E$HYT.SE)
  nSC <- length(y_pc); nSE <- length(y_pe)
  obs <- (log2(mPc) - log2(mPe)) - (log2(mHc) - log2(mHe))
  nul <- vapply(seq_len(n_perm), function(b) {
    p <- perm_trans_pool_paired(nSC, nSE, length(y_hc)); a <- p$sc[seq_len(nSC)]; e <- p$se[seq_len(nSE)]
    (l2(pc_y[a], pc_o[a]) - l2(pe_y[e], pe_o[e])) - (l2(pc_y[-a], pc_o[-a]) - l2(pe_y[-e], pe_o[-e]))
  }, numeric(1))
  perm_pval(obs, nul)
}

## CALIB.TRUE.NULLS: the p-value columns whose null hypothesis each calibration design makes true. Every design built on equal means and dispersions
## (equal, lognormal, bimodal, s2d_known) makes every null true; means makes the dispersion nulls true; dispratio the mean nulls (the dispersions of the
## pooled groups differ, D2); level the trans mean null of both arms (D3).
.CAL.MODES <- "(total|cis|trans|mid|par_sc|par_se)"
CALIB.TRUE.NULLS <- list(
  all   = paste0("^(mean|disp|rdisp|bsize|kbal|kbalsup|wdisp|wrdisp|wbsize|wkbal|wkbalsup)_", .CAL.MODES, "_p$"),
  means = paste0("^(disp|rdisp|wdisp|wrdisp)_", .CAL.MODES, "_p$"),
  dispratio = paste0("^mean_", .CAL.MODES, "_p$"),
  level = "^mean_trans_p(_paired)?$")

## calib_size_summary: the size of every true null in calibration output PR (rows of calib_dataset()), by design, family, design parameter and test, and
## optionally by further strata (for example mean and phi): rejection at 0.05, 0.01 and 0.001 with binomial SEs and the KS p-value of uniformity.
calib_size_summary <- function(PR, strata = NULL) {
  key <- ifelse(PR$design %in% c("equal", "s2d_known"), "all", PR$design)
  by <- c("design", "family", "par", strata)
  do.call(rbind, lapply(split(seq_len(nrow(PR)), do.call(paste, c(PR[by], sep = "|"))), function(ix) {
    cols <- grep(CALIB.TRUE.NULLS[[key[ix[1]]]], names(PR), value = TRUE)
    do.call(rbind, lapply(cols, function(cn) data.frame(PR[ix[1], by, drop = FALSE], test = cn, size_table(PR[ix, cn]), row.names = NULL)))
  })) |> `rownames<-`(NULL)
}

## size_headline: the headline of a size table: tests whose rejection rate at 0.05 lies more than two binomial SEs from 0.05, the range of the rates and
## the worst test.
size_headline <- function(tab, label) {
  tab <- tab[tab$n > 0, ]
  if (!nrow(tab)) return(sprintf("%s: no tests", label))
  off <- abs(tab$reject_0.05 - 0.05) > 2 * pmax(tab$se_0.05, binom_se(0.05, tab$n)); j <- which.max(abs(tab$reject_0.05 - 0.05))
  sprintf("%s: %d tests, rejection at 0.05 from %.3f to %.3f, %d beyond 2 SE of 0.05 (worst %s %s at %.3f, n %d); at 0.01 up to %.4f, at 0.001 up to %.4f",
          label, nrow(tab), min(tab$reject_0.05), max(tab$reject_0.05), sum(off), tab$test[j], paste(tab[j, intersect(c("design", "family", "par", "mean", "phi", "omega"), names(tab))], collapse = "/"),
          tab$reject_0.05[j], tab$n[j], max(tab$reject_0.01), max(tab$reject_0.001))
}

## ============================================================
## Data check helpers (checks_data.R)
## ============================================================

## data_boot_gene: B bootstrap replicates of one gene's fits in every contrast dataset, from pre-drawn indices (make_draws()). Each replicate refits mu and
## phi with the pipeline estimator; with refit_share the share terms are re-estimated on the replicate's resampled hybrid cells (share_cells: c-hat and L
## in HYB order), otherwise the observed eta is applied. Returns a B x (2 x datasets) matrix of mu and phi.
data_boot_gene <- function(g, mats, expos, eta, draws, refit_share = TRUE, share_cells = NULL, groups = .CONTRAST.GROUPS) {
  out <- matrix(NA_real_, length(draws), 2 * length(groups), dimnames = list(NULL, c(paste0("mu.", groups), paste0("phi.", groups))))
  for (b in seq_along(draws)) {
    d <- draws[[b]]; r <- eta
    if (refit_share) { hb <- d$HYB; st <- share_terms(share_cells$chat[hb], share_cells$lib[hb]); r <- share_r_groups(st$r_c, st$r_e) }
    for (k in groups) {
      ix <- d[[.DRAW.KEY[[k]]]]; f <- fit_dispersion(mats[[k]][g, ix], expos[[k]][ix], r[[k]])
      out[b, paste0("mu.", k)] <- f[["mu"]]; out[b, paste0("phi.", k)] <- f[["disp"]]
    }
  }
  out
}

## boot_gene_stats: the summaries of one gene's bootstrap replicates (data_boot_gene()) that A4, A6, C7 and D8 read. Per dataset: cor(mu, phi) on the raw
## scale (a zero phi included) and on the log2 scale (finite only), cor(log2 zeta, log2 theta) and cor(log2 kappa, log2 mu) with zeta = mu phi,
## theta = 1 / phi and kappa = theta / zeta, and the SE ratios SE(log2 kappa) / SE(log2 theta) and SE(log phi) / SE(log mu) (SDs of finite log values).
## Per reported mode: cor(Delta mu, Delta phi) of the mean contrast (mean partition) with the dispersion contrast of its noise source and of its kinetic
## balance source, and the SE of the log2 dispersion contrast (finite replicates). The mean rates are corrected at cbar for the contrasts.
boot_gene_stats <- function(D, cbar, groups = .CONTRAST.GROUPS) {
  cf <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) > 3 && sd(a[ok]) > 0 && sd(b[ok]) > 0) cor(a[ok], b[ok]) else NA_real_ }
  sdf <- function(a) { a <- a[is.finite(a)]; if (length(a) > 3) sd(a) else NA_real_ }
  out <- list()
  for (k in groups) {
    mu <- D[, paste0("mu.", k)]; ph <- D[, paste0("phi.", k)]; lm_ <- log2(mu); lp <- ifelse(ph > 0, log2(ph), NA_real_)
    out[[paste0("cor_mu_phi_raw.", k)]] <- cf(mu, ph); out[[paste0("cor_mu_phi_log2.", k)]] <- cf(lm_, lp)
    out[[paste0("cor_zeta_theta.", k)]] <- cf(lm_ + lp, -lp); out[[paste0("cor_kappa_mu.", k)]] <- cf(-lm_ - 2 * lp, lm_)
    out[[paste0("se_ratio_kappa_theta.", k)]] <- sdf(-lm_ - 2 * lp) / sdf(-lp); out[[paste0("se_ratio_logphi_logmu.", k)]] <- sdf(lp) / sdf(lm_)
  }
  mc <- lapply(seq_len(nrow(D)), function(b) .mode_contrasts(as.list(D[b, paste0("mu.", groups)]) |> setNames(groups), as.list(D[b, paste0("phi.", groups)]) |> setNames(groups), cbar))
  arr <- function(md, col) vapply(mc, function(m) m[md, col], numeric(1))
  for (md in .OUT_MODES) {
    dm <- arr(md, "mean")
    out[[paste0("cor_dmu_dphi_noise.", md)]] <- cf(dm, arr(.DISP_SOURCE[[md]], "disp"))
    out[[paste0("cor_dmu_dphi_kbal.", md)]] <- cf(arr(.KBAL_SOURCE[[md]], "mean"), arr(.KBAL_SOURCE[[md]], "disp"))
    out[[paste0("se_disp.", md)]] <- sdf(arr(.DISP_SOURCE[[md]], "disp"))
  }
  out$cor_mu_Hc_phi_He_kbal <- cf(D[, "mu.HYC.SC.K"], D[, "phi.HYC.SE.K"]); out$cor_mu_He_phi_Hc_kbal <- cf(D[, "mu.HYC.SE.K"], D[, "phi.HYC.SC.K"])
  unlist(out)
}

## tercile_phi_gene: A7 data arm for one gene in one group: phi fitted within each depth tercile of the cells, with the depth offset and with a constant
## offset. With an adequate offset the fit with offsets is flat across terciles; the constant-offset fit absorbs the depth variation within a tercile.
tercile_phi_gene <- function(g, mat, expo, terc) {
  y <- as.numeric(mat[g, ])
  c(vapply(1:3, function(t) fit_dispersion(y[terc == t], expo[terc == t])[["disp"]], numeric(1)),
    vapply(1:3, function(t) fit_dispersion(y[terc == t], rep(1, sum(terc == t)))[["disp"]], numeric(1)))
}

## a7_self_row: A7 simulation arm, grid row i (library share s of the gene, phi). The gene's counts y are a share s of the cell's library on average;
## the offset built from the library including y (L = L_other + y, as the data do) is compared with the offset from the other genes only. Returns the
## relative bias of phi-hat with the self-inclusive offset against the exclusive one, with its Monte Carlo SE.
a7_self_row <- function(i, grid, nrep, lib, n_cells, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  d <- t(vapply(seq_len(nrep), function(r) {
    Lo <- sample(lib, n_cells, TRUE) * (1 - g$share); Oo <- Lo / median(Lo)
    y <- rpois(n_cells, Lo * g$share / (1 - g$share) * sim_pos_factor(n_cells, "gamma", 1, g$phi)); Oi <- (Lo + y) / median(Lo + y)
    c(excl = fit_dispersion(y, Oo)[["disp"]], incl = fit_dispersion(y, Oi)[["disp"]])
  }, c(excl = 0, incl = 0)))
  rb <- d[, "incl"] / d[, "excl"] - 1; rb <- rb[is.finite(rb)]
  data.frame(g, rel_bias_incl_vs_excl = mean(rb), mc_se = sd(rb) / sqrt(length(rb)), median_excl = median(d[, "excl"]), median_incl = median(d[, "incl"]), row.names = NULL)
})

## subsample_cis_var_gene: C2 for one gene: for each subsample of hybrid cells (a list of index vectors, one per fraction), the bootstrap variance of the cis
## log2 contrasts of the mean and of phi (both alleles refitted on the same resampled cells; boot holds B resample vectors per subsample). The hybrid part
## of the trans contrast is the same allele ratio on the trans cells, so its variance follows the same 1 / n law.
subsample_cis_var_gene <- function(g, yc, ye, expo, eta, subs, boot) {
  t(vapply(seq_along(subs), function(j) {
    v <- t(vapply(boot[[j]], function(ix) { ii <- subs[[j]][ix]; fc <- fit_dispersion(yc[g, ii], expo[ii], eta[["HYB.SC"]]); fe <- fit_dispersion(ye[g, ii], expo[ii], eta[["HYB.SE"]])
      c(mean = log2(fc[["mu"]] / fe[["mu"]]), disp = if (fc[["disp"]] > 0 && fe[["disp"]] > 0) log2(fc[["disp"]] / fe[["disp"]]) else NA_real_) }, c(mean = 0, disp = 0)))
    c(n = length(subs[[j]]), var_mean = var(v[is.finite(v[, "mean"]), "mean"]), var_disp = var(v[is.finite(v[, "disp"]), "disp"]))
  }, c(n = 0, var_mean = 0, var_disp = 0)))
}

## split_mean_calls: C6, one random partition of the hybrid cells at fraction f within depth deciles: the cis and trans log2 mean contrasts of every gene
## (corrected allele means of each subset) and their calls |estimate / se| > 1.96 with the observed bootstrap SEs of BOOT.CONTRASTS (a Wald approximation of
## the permutation tests, enough to read call stability). Returns a genes x 2 logical matrix of calls.
split_mean_calls <- function(s, genes, yc, ye, expo, depth, f, mu_pc, mu_pe, cbar, se_cis, se_trans, streams) with_stream(streams[[s]], {
  dec <- cut(rank(depth, ties.method = "first"), 10, labels = FALSE)
  cc <- unlist(lapply(split(seq_along(depth), dec), function(ix) ix[sample.int(length(ix), round(f * length(ix)))])); tt <- setdiff(seq_along(depth), cc)
  m <- function(y, ix, fac) rowSums(y[genes, ix, drop = FALSE]) / sum(expo[ix]) / fac
  cis <- log2(m(yc, cc, 2 * cbar) / m(ye, cc, 2 * (1 - cbar))); trans <- log2(mu_pc / mu_pe) - log2(m(yc, tt, 2 * cbar) / m(ye, tt, 2 * (1 - cbar)))
  cbind(cis = abs(cis / se_cis) > 1.96, trans = abs(trans / se_trans) > 1.96)
})

## col_quartiles: the 25%, 50% and 75% quantiles across genes (rows) of the selected columns of a statistics matrix, one row per column.
col_quartiles <- function(M, cols) t(apply(M[, cols, drop = FALSE], 2, quantile, c(0.25, 0.5, 0.75), na.rm = TRUE))

## var_ratio_gene: A5 for one gene in one group, the variance ratio of the ratio mean mu-hat to the NB maximum likelihood mean mu-tilde at the fitted
## M_i and phi, sum M (1 + phi M) sum M / (1 + phi M) / (sum M)^2, and the document's cap 1 + CV^2(O) (population CV of the offsets).
var_ratio_gene <- function(g, expo, fit) {
  M <- fit[g, "MU"] * expo; phi <- fit[g, "DISP"]
  if (!is.finite(phi) || !is.finite(fit[g, "MU"]) || fit[g, "MU"] <= 0) return(NA_real_)
  sum(M * (1 + phi * M)) * sum(M / (1 + phi * M)) / sum(M)^2
}

## ============================================================
## Text audit helpers (checks_text.R)
## ============================================================

## crossref_record: the Crossref metadata of one DOI (https://api.crossref.org/works/<doi>): title, author family names, container title, volume, first page
## or article number, and year (print, else online, else issued). retrieved is FALSE with the error message when the record cannot be read.
crossref_record <- function(doi) {
  x <- tryCatch(jsonlite::fromJSON(paste0("https://api.crossref.org/works/", utils::URLencode(doi, reserved = FALSE)))$message, error = function(e) conditionMessage(e))
  if (is.character(x)) return(data.frame(doi = doi, retrieved = FALSE, error = x, title = NA, authors = NA, journal = NA, volume = NA, page = NA, year = NA))
  yr <- c(x$`published-print`$`date-parts`[1], x$`published-online`$`date-parts`[1], x$issued$`date-parts`[1])
  data.frame(doi = doi, retrieved = TRUE, error = "", title = paste(x$title, collapse = " "), authors = paste(x$author$family, collapse = ", "),
             journal = paste(x$`container-title`, collapse = " "), volume = if (is.null(x$volume)) NA else x$volume,
             page = if (!is.null(x$page)) sub("-.*", "", x$page) else if (!is.null(x$`article-number`)) x$`article-number` else NA,
             year = if (length(yr) && length(yr[[1]])) yr[[1]][1] else NA)
}

## reference_match: compares one parsed citation (a row with authors, year, title, journal, volume, page) with its Crossref record: year, volume and first
## page equal, every cited family name among the record's authors, and the share of the cited title's words found in the record's title.
reference_match <- function(cit, rec) {
  words <- function(x) unique(strsplit(gsub("[^a-z0-9 ]", " ", tolower(x)), " +")[[1]])
  fam <- trimws(sub(" [A-Z]+$", "", trimws(strsplit(sub(" et al\\.?", "", cit$authors), ",")[[1]])))
  tw <- setdiff(words(cit$title), ""); rw <- words(rec$title)
  data.frame(doi = cit$doi, retrieved = rec$retrieved, year_ok = isTRUE(as.character(cit$year) == as.character(rec$year)),
             volume_ok = isTRUE(as.character(cit$volume) == as.character(rec$volume)) || is.na(cit$volume),
             page_ok = isTRUE(sub("^e", "", tolower(cit$page)) == sub("^e", "", tolower(rec$page))) || is.na(cit$page),
             authors_ok = all(vapply(fam, function(f) grepl(f, rec$authors, fixed = TRUE), logical(1))),
             title_word_share = if (length(tw)) mean(tw %in% rw) else NA_real_, record_title = rec$title, record_journal = rec$journal, record_year = rec$year,
             record_volume = rec$volume, record_page = rec$page, record_authors = rec$authors, error = rec$error)
}
