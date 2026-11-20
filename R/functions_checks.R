###############################################################
### functions_checks.R
### Functions of the validation checks of the methods document (docs/methods.md): the shared simulation generator, the result writer and the
### per-row units of work that the check scripts apply across their grids. Sourced after R/functions.R by analysis/checks/checks_sim.R,
### analysis/checks/checks_data.R, analysis/checks/checks_text.R and the check scripts of analysis/checks/ and analysis/doc_checks.R.
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

## rgamma1: n gamma draws with mean 1 and variance v (shape = rate = 1 / v); all ones when v is 0. The expression factor of the check simulations.
rgamma1 <- function(n, v) if (v > 0) rgamma(n, shape = 1 / v, rate = 1 / v) else rep(1, n)

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

## a4_moment_row: A4 theory check, grid row i (factor family, mu, phi). With equal offsets and the mean and the dispersion fitted on the same cells, the
## first-order covariance of the ratio mean and the moment estimate (S^2 - Ybar) / Ybar^2 is mu (E[(a - 1)^3] - 2 phi^2) / N, zero only when the factor has
## the gamma's third central moment. Compares that value (E[(a - 1)^3] from 4e6 latent draws) with the simulated covariance, and reports the correlations of
## the mean with the moment estimate and with the Cox-Reid estimate. Family "product" is the shared-times-private gamma factor of the generator (Var d = phi / 4).
a4_moment_row <- function(i, grid, nrep, N, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  fac <- function(n) if (g$family == "product") rgamma1(n, g$phi / 4) * rgamma1(n, (1 + g$phi) / (1 + g$phi / 4) - 1) else sim_pos_factor(n, as.character(g$family), 1, g$phi)
  a3 <- mean((fac(4e6) - 1)^3)
  X <- t(vapply(seq_len(nrep), function(r) { y <- rpois(N, g$mu * fac(N)); m <- mean(y); c(m, (var(y) - m) / m^2, fit_dispersion(y, rep(1, N))[["disp"]]) }, numeric(3)))
  cp <- (X[, 1] - mean(X[, 1])) * (X[, 2] - mean(X[, 2]))
  data.frame(g, N = N, third_moment = a3, third_moment_gamma = 2 * g$phi^2, cov_formula = g$mu * (a3 - 2 * g$phi^2) / N, cov_simulated = mean(cp) * nrep / (nrep - 1),
             cov_mc_se = sd(cp) / sqrt(nrep), cor_moment = cor(X[, 1], X[, 2]), cor_cox_reid = cor(X[, 1], X[, 3]), cor_se = 1 / sqrt(nrep - 3), row.names = NULL)
})

## a9cr_row: A9 (Cox-Reid revisit), grid row i (N, Mbar, phi, offsets "equal" or "observed"). In sparse cells both ML estimators are biased upward; the
## question is whether truncation at zero explains it. NREP genes are fitted with the unadjusted profile ML and the Cox-Reid estimate (both truncated at
## zero) and with the Appendix A moment estimator, kept untruncated (disp_raw) and truncated. Returns bias, SD and P(zero) of each against phi.
a9cr_row <- function(i, grid, nrep, offsets, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  est <- t(vapply(seq_len(nrep), function(r) {
    O <- if (g$offsets == "equal") rep(1, g$N) else sample(offsets, g$N, replace = TRUE); mu <- g$Mbar / mean(O)
    y <- rpois(g$N, O * mu * sim_pos_factor(g$N, "gamma", 1, g$phi)); w <- fit_dispersion(y, O, method = "weighted")
    c(profile_ml = profile_ml_phi(y, O), cox_reid = fit_dispersion(y, O)[["disp"]], moment_untruncated = w[["disp_raw"]], moment_truncated = w[["disp"]])
  }, c(profile_ml = 0, cox_reid = 0, moment_untruncated = 0, moment_truncated = 0)))
  do.call(rbind, lapply(colnames(est), function(k) data.frame(g, estimator = k, t(mc_stats(est[, k], g$phi)), row.names = NULL)))
})

## a10_row: A10, grid row i (factor family, overlap o, offsets "equal" or "observed", mu, phi). The mean is fitted on N cells A and the dispersion on N
## cells B that share o N of them. Returns N times the covariance of mu-hat_A with phi-hat_B for the plain moment estimator (S^2 - Ybar) / Ybar^2 on B,
## the Appendix A estimator, the unadjusted profile ML and the Cox-Reid estimate, each with its Monte Carlo SE and the correlation, next to the formula
## o mu (E[(a - 1)^3] - 2 phi^2), derived for the moment estimator with equal offsets.
a10_row <- function(i, grid, nrep, N, offsets, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; nsh <- round(g$overlap * N); ntot <- 2 * N - nsh
  a3 <- mean((sim_pos_factor(4e6, as.character(g$family), 1, g$phi) - 1)^3)
  X <- t(vapply(seq_len(nrep), function(r) {
    O <- if (g$offsets == "equal") rep(1, ntot) else sample(offsets, ntot, replace = TRUE); mu <- g$mu / mean(O)
    y <- rpois(ntot, O * mu * sim_pos_factor(ntot, as.character(g$family), 1, g$phi)); A <- seq_len(N); B <- (N - nsh) + seq_len(N)
    yb <- y[B]; ob <- O[B]; mb <- sum(yb) / sum(ob)
    c(mu_A = sum(y[A]) / sum(O[A]) * mean(O[A]) , moment = (var(yb / ob * mean(ob)) - mean(yb)) / mean(yb)^2, appendixA = fit_dispersion(yb, ob, method = "weighted")[["disp"]],
      profile_ml = profile_ml_phi(yb, ob), cox_reid = fit_dispersion(yb, ob)[["disp"]])
  }, c(mu_A = 0, moment = 0, appendixA = 0, profile_ml = 0, cox_reid = 0)))
  do.call(rbind, lapply(c("moment", "appendixA", "profile_ml", "cox_reid"), function(k) {
    cp <- (X[, "mu_A"] - mean(X[, "mu_A"])) * (X[, k] - mean(X[, k]))
    data.frame(g, N = N, estimator = k, N_cov = N * mean(cp), N_cov_mc_se = N * sd(cp) / sqrt(nrep), N_cov_formula = g$overlap * g$mu * (a3 - 2 * g$phi^2),
               cor = cor(X[, "mu_A"], X[, k]), cor_se = 1 / sqrt(nrep - 3), row.names = NULL)
  }))
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

## f2_rep: one replicate of F2, the large-N limits of Appendix B, for task i (a setting and a replicate number). n cells with depth S_H = 1, Poisson counts,
## additive expression factors a_r = 1 + d + q_r with shifted exponential deviations (d = sd_d (E - 1), E ~ Exp(1), and likewise q, each with mean 0 and
## the stated variance), true shares c ~ beta(c-bar, sigma2_c), and the true mean rates and eta in every estimator. Returns the star estimators, the
## corrected estimators, the variance of the summed factor a*_H, and the SDs of the per-cell terms of the two star estimators, from which the sampling SD
## of each estimator at any number of cells follows (SD / sqrt(N)). Cells are processed in chunks to bound memory.
f2_rep <- function(i, tasks, settings, n, chunk = 1e6, streams) with_stream(streams[[i]], {
  s <- settings[tasks$setting[i], ]; cb <- s$cbar; v <- s$var_c
  eta_c <- v / cb^2; eta_e <- v / (1 - cb)^2; eta_ce <- v / (cb * (1 - cb)); eta_bar <- (eta_c + eta_e) / 2
  pi_ <- cb * s$mu_c / (cb * s$mu_c + (1 - cb) * s$mu_e)
  acc <- c(dd = 0, qq = 0, a1 = 0, a2 = 0, n = 0, dd2 = 0, qq2 = 0)
  for (k in seq_len(ceiling(n / chunk))) {
    m <- min(chunk, n - (k - 1) * chunk)
    c_i <- sim_draw_share(m, cb, v); d <- sqrt(s$s2d) * (rexp(m) - 1)
    ac <- 1 + d + sqrt(s$s2q) * (rexp(m) - 1); ae <- 1 + d + sqrt(s$s2q) * (rexp(m) - 1)
    yc <- rpois(m, c_i * s$mu_c * ac); ye <- rpois(m, (1 - c_i) * s$mu_e * ae)
    xc <- yc / (cb * s$mu_c); xe <- ye / ((1 - cb) * s$mu_e)                   # scaled counts with the true rates (O mu* = S c-bar mu)
    td <- xc * xe; tq <- 0.5 * (xc - xe)^2 - 0.5 * (yc / (cb * s$mu_c)^2 + ye / ((1 - cb) * s$mu_e)^2)
    astar <- pi_ * (c_i / cb) * ac + (1 - pi_) * ((1 - c_i) / (1 - cb)) * ae
    acc <- acc + c(sum(td), sum(tq), sum(astar), sum(astar^2), m, sum(td^2), sum(tq^2))
  }
  d_star <- acc[["dd"]] / acc[["n"]] - 1; q_star <- acc[["qq"]] / acc[["n"]]
  sd2 <- (1 + d_star) / (1 - eta_ce) - 1; sq2 <- (q_star - (eta_bar + eta_ce) * (1 + sd2)) / (1 + eta_bar)
  c(setting = tasks$setting[i], d_star = d_star, q_star = q_star, sigma2_d = sd2, sigma2_q = sq2,
    var_astar = acc[["a2"]] / acc[["n"]] - (acc[["a1"]] / acc[["n"]])^2,
    sd_term_d = sqrt(acc[["dd2"]] / acc[["n"]] - (acc[["dd"]] / acc[["n"]])^2), sd_term_q = sqrt(acc[["qq2"]] / acc[["n"]] - q_star^2))
})

## f2_summary: F2 for one setting from its replicates (rows of f2_rep() output). For each quantity: its target (the stated limit of a star estimator, the
## true variance for a corrected one, the formula for Var[a*_H]), the share term the correction removes (target minus the share-free value), the mean
## estimate over replicates, its error with the Monte Carlo SE (SD over replicates / sqrt(R)) and z, and the sampling SD of the estimator at n_obs cells.
## The error is resolved only when it exceeds about two Monte Carlo SEs; an error below that is also below the sampling SD at n_obs by sqrt(R n / n_obs).
f2_summary <- function(k, settings, reps, n, n_obs) {
  s <- settings[k, ]; x <- reps[reps[, "setting"] == k, , drop = FALSE]; cb <- s$cbar; v <- s$var_c
  eta_c <- v / cb^2; eta_e <- v / (1 - cb)^2; eta_ce <- v / (cb * (1 - cb))
  pi_ <- cb * s$mu_c / (cb * s$mu_c + (1 - cb) * s$mu_e); phiH <- s$s2d + s$s2q
  lim_d <- s$s2d - eta_ce * (1 + s$s2d)
  lim_q <- s$s2q + 0.5 * (1 + phiH) * v * (1 / cb^2 + 1 / (1 - cb)^2) + eta_ce * (1 + s$s2d)
  phiSig <- s$s2d + (pi_^2 + (1 - pi_)^2) * s$s2q
  form <- phiSig + (1 + s$s2d) * (pi_^2 * eta_c + (1 - pi_)^2 * eta_e - 2 * pi_ * (1 - pi_) * eta_ce) + s$s2q * (pi_^2 * eta_c + (1 - pi_)^2 * eta_e)
  tg <- data.frame(quantity = c("d_star", "q_star", "sigma2_d", "sigma2_q", "var_astar"), target = c(lim_d, lim_q, s$s2d, s$s2q, form),
                   share_term = c(lim_d - s$s2d, lim_q - s$s2q, 0, 0, form - phiSig),
                   sd_term = c(median(x[, "sd_term_d"]), median(x[, "sd_term_q"]), median(x[, "sd_term_d"]), median(x[, "sd_term_q"]), NA))
  R <- nrow(x); est <- colMeans(x[, tg$quantity, drop = FALSE]); se <- apply(x[, tg$quantity, drop = FALSE], 2, sd) / sqrt(R)
  data.frame(setting = k, s[rep(1, nrow(tg)), ], tg, replicates = R, cells_per_replicate = n, estimate = est, error = est - tg$target, mc_se = se,
             z = (est - tg$target) / se, sd_at_n_obs = tg$sd_term / sqrt(n_obs), n_obs = n_obs, row.names = NULL)
}

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
## the parental log2 ratio, pooled by medians across genes) or fixed in advance (f_fixed). Each gene's log2 mean cis and trans contrasts (allele means
## divided by their abundance factors) are tested with the pipeline's bootstrap sign test (sign_boot_p()): B replicates resample each parent and the
## undivided hybrid cells, each copy keeping its cis or trans membership. Returns the p-values of both arms and the estimated f.
c4_dataset <- function(i, G, mu_grid, phi, B_pilot, B, f_fixed, pools, streams) with_stream(streams[[i]], {
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
  bx <- lapply(seq_len(B), function(b) list(h = sample.int(nH, nH, TRUE), c = sample.int(length(ds$O_Pc), length(ds$O_Pc), TRUE), e = sample.int(length(ds$O_Pe), length(ds$O_Pe), TRUE)))
  pv <- function(f) {
    sp <- split_indices_by_depth(nH, f); is_c <- seq_len(nH) %in% sp$c
    ## contr: the log2 mean cis and trans contrasts from the parental cells pc, pe and the hybrid cells h (with repeats), split by membership
    contr <- function(x, pc, pe, h) {
      hc <- h[is_c[h]]; ht <- h[!is_c[h]]
      cis <- (l2(x$y$Hc[hc], x$O$Hc[hc]) - log2(2 * cb)) - (l2(x$y$He[hc], x$O$He[hc]) - log2(2 * (1 - cb)))
      tr  <- (l2(x$y$Pc[pc], x$O$Pc[pc]) - l2(x$y$Pe[pe], x$O$Pe[pe])) - ((l2(x$y$Hc[ht], x$O$Hc[ht]) - log2(2 * cb)) - (l2(x$y$He[ht], x$O$He[ht]) - log2(2 * (1 - cb))))
      c(cis = cis, trans = tr)
    }
    vapply(genes, function(x) {
      obs <- contr(x, seq_along(x$y$Pc), seq_along(x$y$Pe), seq_len(nH))
      rep_ <- vapply(bx, function(k) contr(x, k$c, k$e, k$h), numeric(2))
      c(cis = if (is.finite(obs[["cis"]])) sign_boot_p(rep_["cis", ], sign(obs[["cis"]])) else NA_real_,
        trans = if (is.finite(obs[["trans"]])) sign_boot_p(rep_["trans", ], sign(obs[["trans"]])) else NA_real_)
    }, numeric(2))
  }
  pe <- pv(f_est); pf <- pv(f_fixed)
  data.frame(dataset = i, gene = seq_len(G), f_est = f_est, f_fixed = f_fixed, p_cis_est = pe["cis", ], p_trans_est = pe["trans", ], p_cis_fixed = pf["cis", ], p_trans_fixed = pf["trans", ])
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
  if (n == 0) return(data.frame(n = 0, t(setNames(rep(NA_real_, length(levels)), paste0("reject_", levels))), t(setNames(rep(NA_real_, length(levels)), paste0("se_", levels))), ks_p = NA_real_))
  r <- vapply(levels, function(a) mean(p <= a), numeric(1))
  data.frame(n = n, t(setNames(r, paste0("reject_", levels))), t(setNames(binom_se(r, n), paste0("se_", levels))),
             ks_p = suppressWarnings(ks.test(p, "punif")$p.value))
}

## estimator_pair_gene: one gene's fits for the estimator comparison and leak test (analysis/checks/estimator_leak.R). For each group (Pc, Pe, Hc, He) the
## ratio mean and the dispersion by the Cox-Reid NB fit and by the Appendix A moment estimator, both with the group's share correction, in two arms:
## "same", mean and dispersions on all cells of the group, and "disjoint", the mean on half A and both dispersions on half B (halves: one index vector
## per group), so the errors of the mean and of the dispersion share no cells. trim > 0 leaves out of every fit the cells with the largest scaled
## counts y / O (the top trim fraction of the cells in that fit, at least one), which shows how much a few cells carry each estimate. Returns a named
## vector arm.group.quantity (quantity mu, mbar, phi_cr, phi_mom).
estimator_pair_gene <- function(g, mats, expos, eta, halves, trim = 0) {
  out <- c()
  for (k in names(mats)) {
    y <- as.numeric(mats[[k]][g, ]); e <- expos[[k]]; A <- halves[[k]]; B <- setdiff(seq_along(e), A)
    for (arm in c("same", "disjoint")) {
      im <- if (arm == "same") seq_along(e) else A; ip <- if (arm == "same") seq_along(e) else B
      if (trim > 0) for (s in c("im", "ip")) {
        ix <- get(s); drop <- ix[order(y[ix] / e[ix], decreasing = TRUE)][seq_len(max(1, ceiling(trim * length(ix))))]
        assign(s, setdiff(ix, drop))
      }
      mu <- sum(y[im]) / sum(e[im])
      out[paste(arm, k, c("mu", "mbar", "phi_cr", "phi_mom"), sep = ".")] <-
        c(mu, mu * mean(e[ip]), fit_dispersion(y[ip], e[ip], eta[[k]], "NB")[["disp"]], fit_dispersion(y[ip], e[ip], eta[[k]], "weighted")[["disp"]])
    }
  }
  out
}

## leak_contrasts: the cis and trans contrasts of one arm of estimator_pair_gene() output (fits, genes x columns), as a long data frame: x = Delta log2 mu
## (allele means divided by their abundance factors 2 cbar and 2 (1 - cbar); trans = parental log2 ratio minus allele log2 ratio) and, for each estimator,
## the dispersion contrast on the log2 scale (scale "log2", a group with phi = 0 gives NA) or the raw difference of phi (scale "raw", phi = 0 included);
## diff = Cox-Reid minus moment. tercile splits the genes by their mean count per cell averaged over the four groups.
leak_contrasts <- function(fits, arm, cbar, scale = c("log2", "raw")) {
  scale <- match.arg(scale)
  F <- fits[, startsWith(colnames(fits), paste0(arm, ".")), drop = FALSE]; colnames(F) <- sub(paste0("^", arm, "\\."), "", colnames(F))
  lmu <- list(Pc = log2(F[, "Pc.mu"]), Pe = log2(F[, "Pe.mu"]), Hc = log2(F[, "Hc.mu"] / (2 * cbar)), He = log2(F[, "He.mu"] / (2 * (1 - cbar))))
  tr <- function(v) if (scale == "log2") ifelse(v > 0, log2(v), NA_real_) else v
  ph <- lapply(c(cr = "phi_cr", mom = "phi_mom"), function(q) lapply(c(Pc = "Pc", Pe = "Pe", Hc = "Hc", He = "He"), function(k) tr(F[, paste0(k, ".", q)])))
  terc <- cut(rank(rowMeans(F[, paste0(c("Pc", "Pe", "Hc", "He"), ".mbar")]), ties.method = "first"), 3, labels = c("low", "mid", "high"))
  do.call(rbind, lapply(c("cis", "trans"), function(md) {
    d <- function(l) if (md == "cis") l$Hc - l$He else (l$Pc - l$Pe) - (l$Hc - l$He)
    data.frame(gene = rownames(F), arm = arm, scale = scale, contrast = md, dlog2mu = d(lmu), dphi_cr = d(ph$cr), dphi_mom = d(ph$mom),
               diff = d(ph$cr) - d(ph$mom), tercile = terc, row.names = NULL)
  }))
}

## linear_excess_gene: L6 for one gene in one group, a test for a variance term linear in the mean that the NB lacks. Under Var(Y_i) = (1 + a) M_i + phi M_i^2
## with M_i = O_i mu and mu fitted as sum(y) / sum(O), E[(y_i - M_i)^2] = (1 + a) A_i + phi B_i with A_i = M_i (1 - p_i), B_i = M_i^2 (1 - 2 p_i) + p_i^2 sum M^2
## and p_i = O_i / sum O (Appendix A). So z_i = (y_i - M_i)^2 - A_i has expectation a A_i + phi B_i, and a and phi follow from a regression of z on A and B
## without intercept, identified by the spread of the offsets. Weights 1 / [M_i (1 + phi0 M_i)]^2 with phi0 the Cox-Reid fit; SEs by the sandwich (HC0).
## An NB gene has a = 0; a > 0 gives the NB fit phi-hat of about phi + a / M, a dispersion that falls with the mean. beta (named by gene) replaces the
## offsets by O^beta, the gene's own depth scaling (depth_slope_gene()). Returns a, its SE, phi_lin, mbar.
linear_excess_gene <- function(g, mat, expo, beta = NULL) {
  y <- as.numeric(mat[g, ]); if (!is.null(beta)) expo <- expo^beta[[g]]; mu <- sum(y) / sum(expo)
  if (!(mu > 0)) return(c(a = NA, a_se = NA, phi_lin = NA, mbar = NA))
  M <- mu * expo; p <- expo / sum(expo); A <- M * (1 - p); B <- M^2 * (1 - 2 * p) + p^2 * sum(M^2)
  phi0 <- max(fit_dispersion(y, expo)[["disp"]], 0); w <- 1 / (M * (1 + phi0 * M))^2; z <- (y - M)^2 - A
  X <- cbind(A, B); XtWX <- crossprod(X * w, X); bhat <- solve(XtWX, crossprod(X * w, z)); r <- as.vector(z - X %*% bhat)
  V <- solve(XtWX) %*% crossprod(X * (w * r)) %*% solve(XtWX)
  c(a = bhat[1], a_se = sqrt(V[1, 1]), phi_lin = bhat[2], mbar = mu * mean(expo))
}

## leak_regression: the leak test for one contrast. y = Delta log2 phi (Cox-Reid) - Delta log2 phi (moment) regressed across genes on x = Delta log2 mu.
## Returns the OLS slope with its SE, a 95% interval from B gene-resampling draws (indices drawn first), the Spearman correlation, the slopes of each
## estimator's own Delta log2 phi on x, and the number of genes.
leak_regression <- function(x, y, y_cr, y_mom, B = 1000, seed = 1) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; y_cr <- y_cr[ok]; y_mom <- y_mom[ok]; n <- length(x)
  if (n < 10) return(data.frame(genes = n, slope = NA, slope_se = NA, boot_lo = NA, boot_hi = NA, spearman = NA, slope_cr = NA, slope_mom = NA))
  f <- summary(lm(y ~ x))$coefficients
  idx <- with_local_seed(seed, lapply(seq_len(B), function(b) sample.int(n, n, TRUE)))
  bs <- vapply(idx, function(i) unname(coef(lm(y[i] ~ x[i]))[2]), numeric(1))
  data.frame(genes = n, slope = f[2, 1], slope_se = f[2, 2], boot_lo = unname(quantile(bs, 0.025)), boot_hi = unname(quantile(bs, 0.975)),
             spearman = cor(x, y, method = "spearman"), slope_cr = unname(coef(lm(y_cr ~ x))[2]), slope_mom = unname(coef(lm(y_mom ~ x))[2]))
}

## ============================================================
## Data check helpers (checks_data.R)
## ============================================================

## data_boot_gene: B bootstrap replicates of one gene's fits in every contrast dataset, from pre-drawn indices (make_draws()). Each replicate refits mu and
## phi with the pipeline estimator; with refit_share the share terms are re-estimated on the replicate's resampled hybrid cells (share_cells: c-hat and L
## in HYB order), otherwise the observed eta is applied; the summed counts HYB.COMB, when in groups, take the gene's share_sub (named by gene). Returns a
## B x (2 x datasets) matrix of mu and phi.
data_boot_gene <- function(g, mats, expos, eta, draws, refit_share = TRUE, share_cells = NULL, groups = .CONTRAST.GROUPS, share_sub = NULL) {
  out <- matrix(NA_real_, length(draws), 2 * length(groups), dimnames = list(NULL, c(paste0("mu.", groups), paste0("phi.", groups))))
  for (b in seq_along(draws)) {
    d <- draws[[b]]; r <- eta
    if (refit_share) { hb <- d$HYB; st <- share_terms(share_cells$chat[hb], share_cells$lib[hb]); r <- share_r_groups(st$r_c, st$r_e) }
    for (k in groups) {
      ix <- d[[.DRAW.KEY[[k]]]]; f <- fit_dispersion(mats[[k]][g, ix], expos[[k]][ix], r[[k]], "NB", if (k == "HYB.COMB" && !is.null(share_sub)) share_sub[[g]] else 0)
      out[b, paste0("mu.", k)] <- f[["mu"]]; out[b, paste0("phi.", k)] <- f[["disp"]]
    }
  }
  out
}

## boot_gene_stats: the summaries of one gene's bootstrap replicates (data_boot_gene()) that A4, A6, C7 and D8 read. Per dataset: cor(mu, phi) on the raw
## scale (a zero phi included) and on the log2 scale (finite only), cor(log2 zeta, log2 theta) and cor(log2 kappa, log2 mu) with zeta = mu phi,
## theta = 1 / phi and kappa = theta / zeta, and the SE ratios SE(log2 kappa) / SE(log2 theta) and SE(log phi) / SE(log mu) (SDs of finite log values).
## Per reported mode: cor(Delta mu, Delta phi) of the mean contrast (mean partition) with the dispersion contrast of its noise source and of its kinetic
## balance source and of its burst size source (each from one partition's cells), and the SE of the log2 dispersion contrast (finite replicates). The mean
## rates are corrected at cbar for the contrasts.
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
    out[[paste0("cor_dmu_dphi_bsize.", md)]] <- cf(arr(.BSIZE_SOURCE[[md]], "mean"), arr(.BSIZE_SOURCE[[md]], "disp"))
    out[[paste0("se_disp.", md)]] <- sdf(arr(.DISP_SOURCE[[md]], "disp"))
  }
  out$cor_mu_Hc_phi_He_kbal <- cf(D[, "mu.HYC.SC.K"], D[, "phi.HYC.SE.K"]); out$cor_mu_He_phi_Hc_kbal <- cf(D[, "mu.HYC.SE.K"], D[, "phi.HYC.SC.K"])
  out$cor_mu_Hc_phi_He_bsize <- cf(D[, "mu.HYC.SC.Z"], D[, "phi.HYC.SE.Z"]); out$cor_mu_He_phi_Hc_bsize <- cf(D[, "mu.HYC.SE.Z"], D[, "phi.HYC.SC.Z"])
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

## factor_moment_table: A11, the second and third moments of each gene's expression factor in one group, from factorial moments. For Y_i ~ Poisson(M_i a_i)
## with M_i = O_i mu-hat, E[Y(Y - 1)] = M^2 E[a^2] and E[Y(Y - 1)(Y - 2)] = M^3 E[a^3], so the ratios of sums m2 = sum y(y - 1) / sum M^2 and
## m3 = sum y(y - 1)(y - 2) / sum M^3 estimate E[a^2] and E[a^3] without any Poisson noise. phi_fm = m2 - 1 is the factor variance, k3 = m3 - 3 m2 + 2 its
## third central moment, and excess = k3 - 2 phi_fm^2 its departure from a gamma factor of the same variance, the quantity that sets the covariance of
## the mean and dispersion errors on shared cells (mu (E[(a - 1)^3] - 2 phi^2) / N). cols resamples cells (a bootstrap draw); genes are processed in blocks.
factor_moment_table <- function(mat, expo, genes = rownames(mat), cols = seq_along(expo), block = 500) {
  e <- expo[cols]; out <- vector("list", 0)
  for (rows in split(genes, ceiling(seq_along(genes) / block))) {
    y <- as.matrix(mat[rows, cols, drop = FALSE]); mu <- rowSums(y) / sum(e)
    m2 <- rowSums(y * (y - 1)) / (mu^2 * sum(e^2)); m3 <- rowSums(y * (y - 1) * (y - 2)) / (mu^3 * sum(e^3))
    phi <- m2 - 1; k3 <- m3 - 3 * m2 + 2
    out[[length(out) + 1]] <- data.frame(gene = rows, mbar = mu * mean(e), phi_fm = phi, k3 = k3, k3_gamma = 2 * phi^2, excess = k3 - 2 * phi^2, row.names = NULL)
  }
  do.call(rbind, out)
}

## factor_excess_test: A11 test of one gene in one group against a gamma factor. The gene's counts are simulated B times as NB on the same offsets with its
## fitted mean and factorial-moment dispersion (a gamma factor of the same variance), and the excess E[(a - 1)^3] - 2 phi^2 is re-estimated in each draw
## (factor_moment_table()), so the null carries the estimator's own finite-sample bias and skew. Returns the observed excess, the mean of the null excess,
## the observed and null ratios k3 / (2 phi^2), the standardized deviation z = (excess - null mean) / null SD, and the two-sided add-one p-value.
factor_excess_test <- function(g, mat, expo, B) {
  obs <- factor_moment_table(mat, expo, genes = g)
  if (!is.finite(obs$excess) || obs$phi_fm <= 0) return(c(excess = NA, null_excess = NA, ratio = NA, null_ratio = NA, z = NA, p = NA))
  M <- obs$mbar / mean(expo) * expo
  sim <- matrix(rnbinom(B * length(expo), mu = rep(M, each = B), size = 1 / obs$phi_fm), B); rownames(sim) <- paste0("s", seq_len(B))
  nul <- factor_moment_table(sim, expo)
  c(excess = obs$excess, null_excess = mean(nul$excess, na.rm = TRUE), ratio = obs$k3 / obs$k3_gamma, null_ratio = mean(nul$k3 / nul$k3_gamma, na.rm = TRUE),
    z = (obs$excess - mean(nul$excess, na.rm = TRUE)) / sd(nul$excess, na.rm = TRUE),
    p = { d <- nul$excess - mean(nul$excess, na.rm = TRUE); d <- d[!is.na(d)]; (1 + sum(abs(d) >= abs(obs$excess - mean(nul$excess, na.rm = TRUE)))) / (1 + length(d)) })
}

## top_cell_tally: L7, for every cell of one group the number of genes for which it is among the top frac of cells, ranked by the upper-tail probability
## P(Y >= y) of its count under the gene's fitted NB (mean mu O_i, dispersion phi; Poisson when phi = 0). The tail probability puts cells of every depth on
## one scale, so under the model each cell is top for about frac x genes whatever its library; a cell top for many more genes carries counts that are
## extreme across the transcriptome. mu and phi are named by gene. Genes are processed in blocks.
top_cell_tally <- function(mat, expo, genes, frac, mu, phi, block = 200) {
  k <- max(1, ceiling(frac * length(expo))); tally <- integer(length(expo))
  for (rows in split(genes, ceiling(seq_along(genes) / block))) {
    y <- as.matrix(mat[rows, , drop = FALSE]); M <- outer(mu[rows], expo); ph <- phi[rows]
    lp <- matrix(NA_real_, nrow(y), ncol(y))
    for (j in seq_along(rows)) lp[j, ] <- if (is.finite(ph[j]) && ph[j] > 0) pnbinom(y[j, ] - 1, size = 1 / ph[j], mu = M[j, ], lower.tail = FALSE, log.p = TRUE) else ppois(y[j, ] - 1, M[j, ], lower.tail = FALSE, log.p = TRUE)
    top <- apply(lp, 1, function(v) order(v)[seq_len(k)])
    tally <- tally + tabulate(as.vector(top), nbins = length(expo))
  }
  tally
}

## depth_slope_gene: L8, how one gene's counts scale with cell depth in one group: the slope beta of log E[y] = alpha + beta log O (quasi-Poisson fit). The
## library-size offset assumes beta = 1, every gene a fixed share of the library; beta < 1 or > 1 means the gene's share falls or rises with depth, so the
## offset leaves a depth-dependent misfit in the counts that the NB reads as extra variance. Returns beta, its SE and the mean count per cell.
depth_slope_gene <- function(g, mat, expo) {
  y <- as.numeric(mat[g, ]); lo <- log(expo)
  if (sum(y) == 0) return(c(beta = NA, beta_se = NA, mbar = 0))
  f <- tryCatch(suppressWarnings(glm.fit(cbind(1, lo), y, family = quasipoisson())), error = function(e) NULL)
  if (is.null(f) || !f$converged) return(c(beta = NA, beta_se = NA, mbar = mean(y)))
  disp <- sum(f$weights * f$residuals^2) / f$df.residual
  V <- tryCatch(solve(crossprod(cbind(1, lo) * sqrt(f$weights))) * disp, error = function(e) NULL)
  c(beta = unname(f$coefficients[2]), beta_se = if (is.null(V)) NA else sqrt(V[2, 2]), mbar = mean(y))
}

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

## ============================================================
## State-conditioned dispersion (state_conditioned.R)
## ============================================================

## state_cr_glm: the NB dispersion of counts y with log E[y] = off + X b, phi the maximizer over log(phi) of the Cox-Reid adjusted profile likelihood
## sum NB loglik(y; m(b_phi), phi) - 1/2 log det(X' W X), W = m / (1 + phi m), with b_phi from IRLS at each phi (warm-started from the previous phi). With X a
## single column of ones this is the pipeline's Cox-Reid fit, except that the mean is the NB maximum likelihood value instead of the ratio sum(y) /
## sum(exp(off)); state_conditioned.R fits both models with this function so that the covariates are the only difference. phi = 0 (the Poisson limit)
## when the adjusted likelihood at 0 is at least its interior maximum. Returns phi, the adjusted log likelihood at phi and the coefficients.
state_cr_glm <- function(y, off, X, lrange = -rev(LOG.THETA.RANGE), maxit = 50, tol = 1e-8) {
  stopifnot(is.numeric(y), all(y >= 0), length(off) == length(y), is.matrix(X), nrow(X) == length(y))
  b <- c(log(sum(y) / sum(exp(off))), rep(0, ncol(X) - 1)); lfy <- lgamma(y + 1)
  irls <- function(phi) {
    for (it in seq_len(maxit)) {
      lp <- off + drop(X %*% b); m <- exp(lp); w <- m / (1 + phi * m)
      bn <- tryCatch(drop(solve(crossprod(X, w * X), crossprod(X, w * (lp - off + (y - m) / m)))), error = function(e) NULL)
      if (is.null(bn) || any(!is.finite(bn))) return(FALSE)
      done <- max(abs(bn - b)) < tol; b <<- bn
      if (done) break
    }
    TRUE
  }
  adj_ll <- function(phi) {
    if (!irls(phi)) return(-Inf)
    m <- exp(off + drop(X %*% b)); w <- m / (1 + phi * m)
    ll <- if (phi == 0) sum(y * log(m) - m - lfy) else { k <- 1 / phi; sum(lgamma(y + k) - lgamma(k) - lfy + k * log(k / (k + m)) + y * log(m / (k + m))) }
    ld <- determinant(crossprod(X, w * X), logarithm = TRUE)$modulus
    ll - 0.5 * as.numeric(ld)
  }
  if (sum(y) == 0) return(list(phi = NA_real_, adj_ll = NA_real_, coef = rep(NA_real_, ncol(X))))
  ll0 <- adj_ll(0); b0 <- b
  opt <- optimize(function(lp) adj_ll(exp(lp)), lrange, maximum = TRUE, tol = 1e-4)
  if (!is.finite(opt$objective) || ll0 >= opt$objective) { b <- b0; return(list(phi = 0, adj_ll = ll0, coef = b0)) }
  phi <- exp(opt$maximum); val <- adj_ll(phi)                 # refit at the maximum, so b belongs to it
  list(phi = phi, adj_ll = val, coef = b)
}

## state_gene: one gene in every group, its dispersion without and with the cell-state covariates. totals[[k]] is a cells x modules matrix of each module's
## summed counts per cell (in the hybrid cells both alleles together, so the state is the cell's and the two alleles share it); members[[s]] holds the
## genes of module s. The gene's own counts are left out of every module it belongs to, so its own variance never defines its state. Each covariate is
## log2((module counts + 0.5) / offset), standardized within the group; a constant covariate is dropped. Both fits take the group's share correction,
## (phi - eta) / (1 + eta) truncated at zero, as the pipeline does. Returns, per group, mu (ratio estimate), mbar, phi_pipe (fit_dispersion(), the
## pipeline fit), phi0 (state_cr_glm(), intercept only), phiZ (with the covariates), lr (twice the gain in adjusted log likelihood) and one coefficient
## per module (change in log mean per SD of the state score).
state_gene <- function(g, mats, expos, eta, totals, members) {
  out <- c()
  for (k in names(mats)) {
    y <- as.numeric(mats[[k]][g, ]); e <- expos[[k]]
    own <- if (k %in% c("Hc", "He")) as.numeric(mats$Hc[g, ]) + as.numeric(mats$He[g, ]) else y
    Z <- vapply(names(members), function(s) log2(totals[[k]][, s] - (if (g %in% members[[s]]) own else 0) + 0.5) - log2(e), numeric(length(e)))
    Z <- scale(Z); keep <- apply(Z, 2, function(v) all(is.finite(v)))
    f0 <- state_cr_glm(y, log(e), matrix(1, length(y), 1))
    fz <- state_cr_glm(y, log(e), cbind(1, Z[, keep, drop = FALSE]))
    corr <- function(p) if (is.finite(p)) max(0, (p - eta[[k]]) / (1 + eta[[k]])) else NA_real_
    beta <- setNames(rep(NA_real_, length(members)), names(members)); beta[keep] <- fz$coef[-1]
    mu <- sum(y) / sum(e)
    out[paste(k, c("mu", "mbar", "phi_pipe", "phi0", "phiZ", "lr", paste0("beta_", names(members))), sep = ".")] <-
      c(mu, mu * mean(e), fit_dispersion(y, e, eta[[k]], "NB")[["disp"]], corr(f0$phi), corr(fz$phi), 2 * (fz$adj_ll - f0$adj_ll), beta)
  }
  out
}

## ============================================================
## Estimator mean squared error (estimator_mse.R)
## ============================================================

## mse_cell: one cell of the estimator comparison, row i of grid (group, family, mu, phi). R datasets are simulated on the group's observed offsets
## (expos[[group]], the Section 2 depth factors), each count Poisson with mean O mu a and the expression factor a of mean one and variance phi drawn from the
## family (sim_pos_factor(): gamma gives NB counts, lognormal and bimodal factors depart from the NB in their third moment). Each dataset is fitted with the
## pipeline's Cox-Reid estimator (fit_dispersion(method = "NB")) and the Appendix A weighted method of moments (method = "weighted"), both truncated at zero as
## in the pipeline. The replicates are dealt into nbatch equal batches; in each batch the bias (mean estimate minus phi), the variance and the mean squared
## error of the estimates are computed, and the result is their mean over batches with its SE (SD over batches / sqrt(nbatch)). Relative versions divide the
## bias by phi and the variance and MSE by phi^2, so genes of different phi share one scale. Returns one row per estimator with the share of estimates at
## zero and the number of finite estimates.
mse_cell <- function(i, grid, expos, R, nbatch, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; e <- expos[[as.character(g$group)]]; n <- length(e)
  stopifnot(is.numeric(R), R >= 2 * nbatch, nbatch >= 2, g$mu > 0, g$phi > 0)
  est <- vapply(seq_len(R), function(r) { y <- rpois(n, g$mu * e * sim_pos_factor(n, as.character(g$family), 1, g$phi))
    c(NB = fit_dispersion(y, e, 0, "NB")[["disp"]], weighted = fit_dispersion(y, e, 0, "weighted")[["disp"]]) }, numeric(2))
  b <- rep(seq_len(nbatch), length.out = R)
  do.call(rbind, lapply(rownames(est), function(m) {
    x <- est[m, ]; ok <- is.finite(x)
    st <- vapply(seq_len(nbatch), function(k) { z <- x[b == k & ok]; c(bias = mean(z) - g$phi, var = var(z), mse = mean((z - g$phi)^2)) }, numeric(3))
    av <- rowMeans(st); se <- apply(st, 1, sd) / sqrt(nbatch); sc <- c(bias = g$phi, var = g$phi^2, mse = g$phi^2)
    data.frame(g, estimator = m, n_cells = n, mbar = g$mu * mean(e), n_ok = sum(ok), zero = mean(x[ok] == 0),
               bias = av[["bias"]], bias_se = se[["bias"]], var = av[["var"]], var_se = se[["var"]], mse = av[["mse"]], mse_se = se[["mse"]],
               rel_bias = av[["bias"]] / sc[["bias"]], rel_bias_se = se[["bias"]] / sc[["bias"]], rel_var = av[["var"]] / sc[["var"]], rel_var_se = se[["var"]] / sc[["var"]],
               rel_mse = av[["mse"]] / sc[["mse"]], rel_mse_se = se[["mse"]] / sc[["mse"]], row.names = NULL)
  }))
})

## mse_contrast_cell: one cell of the contrast-level comparison (M2 of estimator_mse.R), row i of grid (family, mu, phi, design, effect). Each of R
## datasets simulates the two parents on their observed offsets (expos$Pc, expos$Pe) and the two alleles of every hybrid cell on the allele offset (expos$H),
## all with mean count O mu. The alleles of a cell share an extrinsic factor of variance ext_share x the smaller allele dispersion and carry private factors
## that bring each allele to its dispersion ((1 + phi) = (1 + s2d)(1 + v_q)); every factor has mean one and the family's shape (sim_pos_factor()). The design
## sets the true log2 dispersion ratios: "null" all equal; "cis" alleles and parents both differ by 2^effect (cis = total = effect, trans = 0); "trans" only
## the parents differ (trans = total = effect, cis = 0). The hybrid cells are split as in the pipeline into cis cells and trans cells (hyb_split, positions in
## expos$H). Each dataset is fitted with both estimators, truncated at zero, and the three log2 contrasts are formed: total = log2 phi_Pc - log2 phi_Pe,
## cis on the cis cells, trans = total - the allele ratio on the trans cells. A contrast is defined when every group it uses has phi > 0 (a zero gives
## +-Inf); the bias, variance and MSE against the true contrast are computed over the defined datasets, averaged over nbatch batches with their SE, and the
## share of defined datasets is returned with them. Returns one row per estimator and contrast.
mse_contrast_cell <- function(i, grid, expos, hyb_split, R, nbatch, ext_share, streams) with_stream(streams[[i]], {
  g <- grid[i, ]; r <- 2^g$effect; fam <- as.character(g$family)
  stopifnot(R >= 2 * nbatch, nbatch >= 2, g$mu > 0, g$phi > 0, ext_share >= 0, ext_share <= 1)
  phi <- switch(as.character(g$design), null = c(Pc = g$phi, Pe = g$phi, Hc = g$phi, He = g$phi), cis = c(Pc = g$phi, Pe = g$phi / r, Hc = g$phi, He = g$phi / r),
                trans = c(Pc = g$phi, Pe = g$phi / r, Hc = g$phi, He = g$phi))
  truth <- c(total = log2(phi[["Pc"]] / phi[["Pe"]]), cis = log2(phi[["Hc"]] / phi[["He"]]),
             trans = log2(phi[["Pc"]] / phi[["Pe"]]) - log2(phi[["Hc"]] / phi[["He"]]))
  eP1 <- expos$Pc; eP2 <- expos$Pe; eH <- expos$H; nH <- length(eH); ic <- hyb_split$c; it <- hyb_split$t
  s2d <- ext_share * min(phi[["Hc"]], phi[["He"]]); vq <- function(p) max((1 + p) / (1 + s2d) - 1, 0)
  ests <- c("NB", "weighted"); cons <- c("total", "cis", "trans")
  est <- vapply(seq_len(R), function(rep) {
    yPc <- rpois(length(eP1), g$mu * eP1 * sim_pos_factor(length(eP1), fam, 1, phi[["Pc"]]))
    yPe <- rpois(length(eP2), g$mu * eP2 * sim_pos_factor(length(eP2), fam, 1, phi[["Pe"]]))
    d <- sim_pos_factor(nH, fam, 1, s2d)
    yHc <- rpois(nH, g$mu * eH * d * sim_pos_factor(nH, fam, 1, vq(phi[["Hc"]])))
    yHe <- rpois(nH, g$mu * eH * d * sim_pos_factor(nH, fam, 1, vq(phi[["He"]])))
    unlist(lapply(ests, function(m) {
      l <- log2(c(Pc = fit_dispersion(yPc, eP1, 0, m)[["disp"]], Pe = fit_dispersion(yPe, eP2, 0, m)[["disp"]],
                  Hcc = fit_dispersion(yHc[ic], eH[ic], 0, m)[["disp"]], Hec = fit_dispersion(yHe[ic], eH[ic], 0, m)[["disp"]],
                  Hct = fit_dispersion(yHc[it], eH[it], 0, m)[["disp"]], Het = fit_dispersion(yHe[it], eH[it], 0, m)[["disp"]]))
      c(l[["Pc"]] - l[["Pe"]], l[["Hcc"]] - l[["Hec"]], (l[["Pc"]] - l[["Pe"]]) - (l[["Hct"]] - l[["Het"]]))
    }))
  }, numeric(length(ests) * length(cons)))
  rownames(est) <- paste(rep(ests, each = length(cons)), cons, sep = ".")
  b <- rep(seq_len(nbatch), length.out = R)
  do.call(rbind, lapply(rownames(est), function(k) {
    x <- est[k, ]; ok <- is.finite(x); tr <- truth[[sub("^[^.]+\\.", "", k)]]
    st <- vapply(seq_len(nbatch), function(j) { z <- x[b == j & ok]; c(bias = mean(z) - tr, var = var(z), mse = mean((z - tr)^2)) }, numeric(3))
    av <- rowMeans(st, na.rm = TRUE); se <- apply(st, 1, sd, na.rm = TRUE) / sqrt(rowSums(is.finite(st)))   # batches with too few defined datasets drop out
    data.frame(g, estimator = sub("\\..*$", "", k), contrast = sub("^[^.]+\\.", "", k), truth = tr, defined = mean(ok),
               bias = av[["bias"]], bias_se = se[["bias"]], var = av[["var"]], var_se = se[["var"]], mse = av[["mse"]], mse_se = se[["mse"]], row.names = NULL)
  }))
})

## mse_rare_cell: one cell of M3 (estimator_mse.R), the type I error of the pipeline's dispersion test when a rare cell state occurs at the same frequency p in
## every group (row i of grid: mu, phi, p, J). Each of R datasets draws, per group, which cells are in the rare state (binomial with probability p per cell: the
## parents separately, and once per hybrid cell, so both alleles of a hybrid cell share it) and multiplies their expression factor by J; the bulk factor is
## gamma with dispersion phi (the hybrid alleles share an extrinsic part of variance ext_share x phi). All groups have the same mu, phi, p and J, so every
## true contrast is zero and any rejection is a false positive. The hybrid cells are split into cis and trans cells (hyb_split). For each estimator the
## observed log2 contrasts are total (Pc - Pe), cis (alleles of the cis cells) and trans (total minus the allele ratio of the trans cells). The test is the
## pipeline's bootstrap sign test (sign_boot_p()): B resamples of the parents and of the undivided hybrid cells (each copy keeping its cis or trans membership),
## every group refitted, a failed refit counting as zero and a zero as -Inf on the log2 scale. A contrast is tested, as in analysis.R, when every group it uses
## has an observed phi > 0. Returns, per estimator and contrast: the share of datasets tested, the rejection rate at 0.05 and 0.01 among them with binomial SEs,
## the ratio of the SD of the observed contrast over the tested datasets to their median bootstrap SD (above 1 means the bootstrap is too narrow), and the
## rejection rate among tested datasets whose rare-cell counts differ between the parents by at least two (rare_split).
mse_rare_cell <- function(i, grid, expos, hyb_split, R, B, ext_share, estimators, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  stopifnot(R >= 10, B >= 20, g$mu > 0, g$phi > 0, g$p >= 0, g$p < 1, g$J >= 1, all(estimators %in% c("NB", "weighted")))
  e1 <- expos$Pc; e2 <- expos$Pe; eH <- expos$H; n1 <- length(e1); n2 <- length(e2); nH <- length(eH); ic <- hyb_split$c; it <- hyb_split$t
  s2d <- ext_share * g$phi; vq <- max((1 + g$phi) / (1 + s2d) - 1, 0)
  cons <- c("total", "cis", "trans"); l2 <- function(x) ifelse(x == 0, -Inf, log2(x))
  fit6 <- function(y, ix, m) c(Pc = fit_dispersion(y$Pc[ix$Pc], e1[ix$Pc], 0, m)[["disp"]], Pe = fit_dispersion(y$Pe[ix$Pe], e2[ix$Pe], 0, m)[["disp"]],
                               Hcc = fit_dispersion(y$Hc[ix$c], eH[ix$c], 0, m)[["disp"]], Hec = fit_dispersion(y$He[ix$c], eH[ix$c], 0, m)[["disp"]],
                               Hct = fit_dispersion(y$Hc[ix$t], eH[ix$t], 0, m)[["disp"]], Hec2 = fit_dispersion(y$He[ix$t], eH[ix$t], 0, m)[["disp"]])
  contr <- function(f) { l <- l2(f); x <- c(l[["Pc"]] - l[["Pe"]], l[["Hcc"]] - l[["Hec"]], (l[["Pc"]] - l[["Pe"]]) - (l[["Hct"]] - l[["Hec2"]])); x[is.nan(x)] <- NA; x }
  used <- list(total = c("Pc", "Pe"), cis = c("Hcc", "Hec"), trans = c("Pc", "Pe", "Hct", "Hec2"))
  res <- array(NA_real_, c(length(estimators), length(cons), R, 4), dimnames = list(estimators, cons, NULL, c("p", "obs", "boot_sd", "tested")))
  split2 <- logical(R)
  for (r in seq_len(R)) {
    k1 <- rbinom(n1, 1, g$p); k2 <- rbinom(n2, 1, g$p); kH <- rbinom(nH, 1, g$p); split2[r] <- abs(sum(k1) - sum(k2)) >= 2
    d <- sim_pos_factor(nH, "gamma", 1, s2d) * ifelse(kH == 1, g$J, 1)
    y <- list(Pc = rpois(n1, g$mu * e1 * sim_pos_factor(n1, "gamma", 1, g$phi) * ifelse(k1 == 1, g$J, 1)),
              Pe = rpois(n2, g$mu * e2 * sim_pos_factor(n2, "gamma", 1, g$phi) * ifelse(k2 == 1, g$J, 1)),
              Hc = rpois(nH, g$mu * eH * d * sim_pos_factor(nH, "gamma", 1, vq)), He = rpois(nH, g$mu * eH * d * sim_pos_factor(nH, "gamma", 1, vq)))
    draws <- lapply(seq_len(B), function(b) { h <- sample.int(nH, nH, replace = TRUE)
      list(Pc = sample.int(n1, n1, replace = TRUE), Pe = sample.int(n2, n2, replace = TRUE), c = h[h %in% ic], t = h[h %in% it]) })
    for (m in estimators) {
      fo <- fit6(y, list(Pc = seq_len(n1), Pe = seq_len(n2), c = ic, t = it), m); obs <- contr(fo)
      fb <- vapply(draws, fit6, numeric(6), y = y, m = m); fb[!is.finite(fb)] <- 0; cb <- apply(fb, 2, contr)      # a failed refit counts as zero
      for (j in seq_along(cons)) {
        ok <- isTRUE(all(fo[used[[cons[j]]]] > 0))
        res[m, j, r, ] <- c(if (ok) sign_boot_p(cb[j, ], sign(obs[j])) else NA_real_, obs[j], sd(cb[j, is.finite(cb[j, ])]), ok)
      }
    }
  }
  do.call(rbind, lapply(estimators, function(m) do.call(rbind, lapply(cons, function(cn) {
    x <- res[m, cn, , ]; t <- x[, "tested"] == 1; pv <- x[t, "p"]; nt <- sum(t); r05 <- mean(pv < 0.05); r01 <- mean(pv < 0.01)
    ps <- x[t & split2, "p"]
    data.frame(g, estimator = m, contrast = cn, datasets = R, tested = mean(t), type1_05 = r05, type1_05_se = binom_se(r05, nt), type1_01 = r01, type1_01_se = binom_se(r01, nt),
               sd_ratio = sd(x[t, "obs"]) / median(x[t, "boot_sd"]), rare_split_share = mean(split2[t]), type1_05_rare_split = if (length(ps)) mean(ps < 0.05) else NA_real_,
               type1_05_rare_split_se = if (length(ps)) binom_se(mean(ps < 0.05), length(ps)) else NA_real_, row.names = NULL)
  }))))
})

## ============================================================
## State groups (analysis/checks/state_groups.R)
## ============================================================

## state_assign: one label per cell on each of three state axes, from the cell's module z-scores (one column per module, named by level). Metabolism:
## derepressed when the derepression score exceeds thr, otherwise respiratory when the respiration score exceeds thr, otherwise fermentative, the
## default of glucose-grown cells; derepression is checked first because derepressed cells also raise respiration. Cell cycle: the phase with the
## highest score, since every cell is at some phase. Stress: stressed when the stress score exceeds thr.
state_assign <- function(z, thr, cycle_levels) {
  stopifnot(is.matrix(z), is.numeric(thr), length(thr) == 1, is.finite(thr), all(c("Derepressed", "Respiratory", "Stressed", cycle_levels) %in% colnames(z)))
  data.frame(metabolism = factor(ifelse(z[, "Derepressed"] > thr, "Derepressed", ifelse(z[, "Respiratory"] > thr, "Respiratory", "Fermentative")),
                                 levels = c("Fermentative", "Respiratory", "Derepressed")),
             cycle      = factor(cycle_levels[max.col(z[, cycle_levels, drop = FALSE], ties.method = "first")], levels = cycle_levels),
             stress     = factor(ifelse(z[, "Stressed"] > thr, "Stressed", "Unstressed"), levels = c("Unstressed", "Stressed")),
             row.names  = rownames(z))
}

## ============================================================
## Midparent design shift (analysis/checks/mid_shift_check.R)
## ============================================================

## median_se: the large-sample standard error of the median of x, 1.2533 MAD / sqrt(n) (exact for normal data, the MAD being scaled to the SD).
median_se <- function(x) 1.2533 * mad(x) / sqrt(length(x))

## mid_shift_gene: one gene's log2 midparent contrasts d (the per-genome hybrid against the mean of the two parents; .mode_contrasts(), mode "mid") of the
## mean and of the dispersion in every cell-bootstrap replicate. draws holds the pre-drawn index sets (Pc, Pe, H: parental cells and undivided hybrid cells),
## drawn once for all genes so the median over genes can be taken replicate by replicate. Each replicate refits the share terms on its resampled hybrid cells
## and the four groups with the pipeline's NB estimator; a failed refit counts as a zero dispersion. Returns a B x 2 matrix (mean, disp); a zero dispersion
## gives +-Inf in disp, an undefined value NA.
mid_shift_gene <- function(g, mats, expos, draws, share_cells) {
  t(vapply(draws, function(d) {
    st <- share_terms(share_cells$chat[d$H], share_cells$lib[d$H]); sr <- share_r_groups(st$r_c, st$r_e)
    ix <- list(MIX.SC = d$Pc, MIX.SE = d$Pe, HYB.SC = d$H, HYB.SE = d$H)
    f <- lapply(setNames(names(ix), names(ix)), function(k) .fit_one(mats[[k]][g, ix[[k]]], expos[[k]][ix[[k]]], sr[[k]]))
    mu <- lapply(f, function(x) unname(x[["mu"]])); ph <- lapply(f, function(x) { v <- unname(x[["disp"]]); if (is.finite(v)) v else 0 })
    .mode_contrasts(mu, ph, st$cbar, modes = "mid")["mid", c("mean", "disp")]
  }, numeric(2)))
}

## ============================================================
## Document checks run by analysis/doc_checks.R
## ============================================================

## dominance_bias_cell: one cell of DOM1, the bias of the dominance contrasts (mid, par_sc, par_se) and of the full-cell total, cis and trans contrasts when
## the hybrid truly differs from its parents. Row i of grid gives the parental mean count per cell (mean), the dispersion (phi), the expression-factor
## family and omega. Each of G genes draws its parental per-genome rates (the Sc rate at mean / mean offset times a lognormal spread of SD 0.3, the Se rate
## 2^N(0, 0.7) times that) and dispersions (phi times the spread, Se 2^N(0, 0.5) times it), and its hybrid alleles depart from their parents by 2^N(0, 0.5)
## in rate and in dispersion, so every dominance contrast has a known nonzero truth. The cells are one sim_design() of the observed design; the observed
## shares give the share terms; every group is fitted with the pipeline's NB estimator and share corrections on all of its cells, and the contrasts are built
## by .mode_contrasts() from the fits and from the true values (allele rates times their abundance factors, so the correction recovers them). Returns, per
## quantity (mean, disp, bsize, kbal) and contrast, the mean of estimate - truth over genes with its SE, the SD over genes and bias / SD.
dominance_bias_cell <- function(i, grid, pools, G, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  stopifnot(G >= 10, g$mean > 0, g$phi > 0)
  ds <- sim_design(pools); st <- share_terms(ds$chat, ds$L); SR <- share_r_groups(st$r_c, st$r_e); cb <- st$cbar
  modes <- c("total", "cis", "trans", "mid", "par_sc", "par_se")
  est <- tru <- array(NA_real_, c(G, length(modes), 4), dimnames = list(NULL, modes, c("mean", "disp", "bsize", "kbal")))
  for (j in seq_len(G)) {
    s <- exp(rnorm(1, 0, 0.3)); mPc <- g$mean * s / mean(ds$O_Pc); mPe <- mPc * 2^rnorm(1, 0, 0.7)
    pPc <- g$phi * exp(rnorm(1, 0, 0.3)); pPe <- pPc * 2^rnorm(1, 0, 0.5)
    mu <- c(Pc = mPc, Pe = mPe, Hc = mPc * 2^rnorm(1, 0, 0.5), He = mPe * 2^rnorm(1, 0, 0.5))
    phi <- c(Pc = pPc, Pe = pPe, Hc = pPc * 2^rnorm(1, 0, 0.5), He = pPe * 2^rnorm(1, 0, 0.5))
    sim <- sim_hybrid_dataset(ds, mu, phi, omega = g$omega, family = g$family)
    f <- list(MIX.SC = .fit_one(sim$y$Pc, sim$O$Pc), MIX.SE = .fit_one(sim$y$Pe, sim$O$Pe),
              HYB.SC = .fit_one(sim$y$Hc, sim$O$Hc, SR[["HYB.SC"]]), HYB.SE = .fit_one(sim$y$He, sim$O$He, SR[["HYB.SE"]]))
    gm <- lapply(f, function(x) unname(x[["mu"]])); gp <- lapply(f, function(x) unname(x[["disp"]]))
    gm[c("HYC.SC", "HYT.SC", "HYC.SE", "HYT.SE")] <- gm[c("HYB.SC", "HYB.SC", "HYB.SE", "HYB.SE")]
    gp[c("HYC.SC", "HYT.SC", "HYC.SE", "HYT.SE")] <- gp[c("HYB.SC", "HYB.SC", "HYB.SE", "HYB.SE")]
    tm <- list(MIX.SC = mu[["Pc"]], MIX.SE = mu[["Pe"]], HYB.SC = 2 * cb * mu[["Hc"]], HYB.SE = 2 * (1 - cb) * mu[["He"]])
    tp <- list(MIX.SC = phi[["Pc"]], MIX.SE = phi[["Pe"]], HYB.SC = phi[["Hc"]], HYB.SE = phi[["He"]])
    tm[c("HYC.SC", "HYT.SC", "HYC.SE", "HYT.SE")] <- tm[c("HYB.SC", "HYB.SC", "HYB.SE", "HYB.SE")]
    tp[c("HYC.SC", "HYT.SC", "HYC.SE", "HYT.SE")] <- tp[c("HYB.SC", "HYB.SC", "HYB.SE", "HYB.SE")]
    oe <- .mode_contrasts(gm, gp, cb, modes); ot <- .mode_contrasts(tm, tp, cb, modes)
    est[j, , ] <- cbind(oe[, "mean"], oe[, "disp"], oe[, "mean"] + oe[, "disp"], -oe[, "mean"] - 2 * oe[, "disp"])
    tru[j, , ] <- cbind(ot[, "mean"], ot[, "disp"], ot[, "mean"] + ot[, "disp"], -ot[, "mean"] - 2 * ot[, "disp"])
  }
  out <- list()
  for (q in dimnames(est)[[3]]) for (md in modes) {
    d <- est[, md, q] - tru[, md, q]; fin <- is.finite(d)
    out[[paste(q, md)]] <- data.frame(g, quantity = q, contrast = md, genes = sum(fin), infinite = sum(!fin), bias = if (sum(fin) > 2) mean(d[fin]) else NA_real_,
                                      bias_se = if (sum(fin) > 2) sd(d[fin]) / sqrt(sum(fin)) else NA_real_, sd = if (sum(fin) > 2) sd(d[fin]) else NA_real_, row.names = NULL)
  }
  out <- do.call(rbind, out); out$bias_over_sd <- out$bias / out$sd; out
})

## rare_cond_cell: one cell of M4, whether the state-conditioned layer (analysis.R Section 2.5) catches the false calls that a rare state makes (M3). Row
## i of grid gives mu, phi, the rare-state frequency p, its fold J and the scenario: "shared" (the state at frequency p in every group, as in M3) or
## "parents" (the state in the parents only, the two-culture case of the data). Each dataset draws the state per cell (once per hybrid cell), multiplies the
## gene's expression factor by J in state cells, and draws the cell's program-module count (the 20 markers) at 40 counts per unit offset, J_mod times higher in
## state cells, with a gamma noise of variance 0.3; the score is the log2 module count per unit offset centred at the tube's median (each parent separately),
## as state_scores() builds it. The gene's true contrasts are all zero. Every group is fitted without and with the score (fit_dispersion(z = )), on all cells
## and on B resamples (parents and undivided hybrid cells, cis or trans membership kept). The unconditioned log2 dispersion contrasts (total, cis, trans) are
## tested with the bootstrap sign test (a failed refit counting as zero, no test where an observed fit is zero), as the pipeline does; a significant call is flagged when the conditioned 95% percentile
## interval contains zero and the interval of the paired change excludes zero. Returns, per contrast: the share tested, the type I error at 0.05, the share of
## false calls flagged, the type I error after flagged calls are set aside, the coverage of zero by the conditioned interval, and the share of all datasets
## flagged (the false-flag rate when p = 0).
rare_cond_cell <- function(i, grid, expos, hyb_split, R, B, ext_share = 0.5, J_mod = 20, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  stopifnot(R >= 10, B >= 20, g$mu > 0, g$phi > 0, g$p >= 0, g$p < 1, g$J >= 1, g$scenario %in% c("shared", "parents"))
  e1 <- expos$Pc; e2 <- expos$Pe; eH <- expos$H; n1 <- length(e1); n2 <- length(e2); nH <- length(eH); ic <- hyb_split$c; it <- hyb_split$t
  s2d <- ext_share * g$phi; vq <- max((1 + g$phi) / (1 + s2d) - 1, 0)
  cons <- c("total", "cis", "trans"); l2 <- function(x) ifelse(x == 0, -Inf, log2(x))
  score <- function(k, e) { L <- log2(rpois(length(e), 40 * e * ifelse(k == 1, J_mod, 1) * sim_pos_factor(length(e), "gamma", 1, 0.3)) + 0.5) - log2(e); L - median(L) }
  fitd <- function(y, e, z) if (is.null(z)) fit_dispersion(y, e, 0, "NB")[["disp"]] else fit_dispersion(y, e, 0, "NB", z = z)[["disp"]]
  fit6 <- function(y, z, ix, cond) { zz <- function(v, j) if (cond) v[j] else NULL
    c(Pc = fitd(y$Pc[ix$Pc], e1[ix$Pc], zz(z$Pc, ix$Pc)), Pe = fitd(y$Pe[ix$Pe], e2[ix$Pe], zz(z$Pe, ix$Pe)),
      Hcc = fitd(y$Hc[ix$c], eH[ix$c], zz(z$H, ix$c)), Hec = fitd(y$He[ix$c], eH[ix$c], zz(z$H, ix$c)),
      Hct = fitd(y$Hc[ix$t], eH[ix$t], zz(z$H, ix$t)), Het = fitd(y$He[ix$t], eH[ix$t], zz(z$H, ix$t))) }
  contr <- function(f) { l <- l2(f); x <- c(l[["Pc"]] - l[["Pe"]], l[["Hcc"]] - l[["Hec"]], (l[["Pc"]] - l[["Pe"]]) - (l[["Hct"]] - l[["Het"]])); x[is.nan(x)] <- NA; x }
  used <- list(total = c("Pc", "Pe"), cis = c("Hcc", "Hec"), trans = c("Pc", "Pe", "Hct", "Het"))
  res <- array(NA_real_, c(length(cons), R, 4), dimnames = list(cons, NULL, c("p", "tested", "flag", "cover")))
  for (r in seq_len(R)) {
    k1 <- rbinom(n1, 1, g$p); k2 <- rbinom(n2, 1, g$p); kH <- if (g$scenario == "shared") rbinom(nH, 1, g$p) else integer(nH)
    d <- sim_pos_factor(nH, "gamma", 1, s2d) * ifelse(kH == 1, g$J, 1)
    y <- list(Pc = rpois(n1, g$mu * e1 * sim_pos_factor(n1, "gamma", 1, g$phi) * ifelse(k1 == 1, g$J, 1)),
              Pe = rpois(n2, g$mu * e2 * sim_pos_factor(n2, "gamma", 1, g$phi) * ifelse(k2 == 1, g$J, 1)),
              Hc = rpois(nH, g$mu * eH * d * sim_pos_factor(nH, "gamma", 1, vq)), He = rpois(nH, g$mu * eH * d * sim_pos_factor(nH, "gamma", 1, vq)))
    z <- list(Pc = score(k1, e1), Pe = score(k2, e2), H = score(kH, 2 * eH))
    draws <- lapply(seq_len(B), function(b) { h <- sample.int(nH, nH, replace = TRUE)
      list(Pc = sample.int(n1, n1, replace = TRUE), Pe = sample.int(n2, n2, replace = TRUE), c = h[h %in% ic], t = h[h %in% it]) })
    all_ix <- list(Pc = seq_len(n1), Pe = seq_len(n2), c = ic, t = it)
    fu <- fit6(y, z, all_ix, FALSE); ou <- contr(fu); oc <- contr(fit6(y, z, all_ix, TRUE))
    bu <- vapply(draws, fit6, numeric(6), y = y, z = z, cond = FALSE); bc <- vapply(draws, fit6, numeric(6), y = y, z = z, cond = TRUE)
    bu[!is.finite(bu)] <- 0; bc[!is.finite(bc)] <- 0                                     # a failed refit counts as zero
    cu <- apply(bu, 2, contr); cc <- apply(bc, 2, contr)
    for (j in seq_along(cons)) {
      ok <- isTRUE(all(fu[used[[cons[j]]]] > 0))
      pv <- if (ok) sign_boot_p(cu[j, ], sign(ou[j])) else NA_real_
      ci <- quantile(cc[j, ], c(0.025, 0.975), type = 1, na.rm = TRUE, names = FALSE)
      di <- quantile(cu[j, ] - cc[j, ], c(0.025, 0.975), type = 1, na.rm = TRUE, names = FALSE)
      cover <- isTRUE(ci[1] <= 0 && ci[2] >= 0)
      res[j, r, ] <- c(pv, ok, isTRUE(ok && pv < 0.05 && cover && (di[1] > 0 || di[2] < 0)), cover)
    }
  }
  do.call(rbind, lapply(cons, function(cn) {
    x <- res[cn, , ]; t <- x[, "tested"] == 1; rej <- t & x[, "p"] < 0.05; nt <- sum(t)
    data.frame(g, contrast = cn, datasets = R, tested = mean(t), type1_05 = mean(x[t, "p"] < 0.05), type1_05_se = binom_se(mean(x[t, "p"] < 0.05), nt),
               flagged_of_false_calls = if (sum(rej)) mean(x[rej, "flag"] == 1) else NA_real_,
               type1_05_unflagged = mean(rej[t] & x[t, "flag"] == 0), cond_cover0 = mean(x[, "cover"] == 1), flag_rate = mean(x[, "flag"] == 1), row.names = NULL)
  }))
})

## ============================================================
## P1: the bootstrap sign test for the dispersion-family contrasts (analysis/checks/sign_boot_test.R, cluster/scripts/sign_boot_grid.R)
## ============================================================

## SIGN.BOOT.STATS: the 24 statistics one dataset yields, contrast x quantity x scale. The quantities are the dispersion phi, the burst frequency
## theta = 1 / phi, the burst size zeta = mu phi and the kinetic balance kappa = theta / zeta = 1 / (mu phi^2); the contrasts are total (Pc - Pe), cis
## (the alleles of the cis cells) and trans (total minus the allele contrast of the trans cells); the scales are log2 and raw.
SIGN.BOOT.STATS <- as.vector(outer(outer(c("total", "cis", "trans"), c("phi", "theta", "zeta", "kappa"), paste, sep = "."), c("log2", "raw"), paste, sep = "."))

## sign_boot_stats: the 24 statistics (SIGN.BOOT.STATS) from the fitted (mu, phi) of the six groups (rows Pc, Pe, Hcc, Hec, Hct, Het; columns mu, phi).
## A zero phi gives log2 phi = -Inf and theta = Inf; +Inf - Inf gives NaN, the undefined value the test counts on the null side.
sign_boot_stats <- function(f) {
  q <- cbind(phi = f[, "phi"], theta = 1 / f[, "phi"], zeta = f[, "mu"] * f[, "phi"], kappa = 1 / (f[, "mu"] * f[, "phi"]^2))
  out <- numeric(0)
  for (sc in c("log2", "raw")) {
    v <- if (sc == "log2") suppressWarnings(log2(q)) else q
    tot <- v["Pc", ] - v["Pe", ]; cis <- v["Hcc", ] - v["Hec", ]; tra <- tot - (v["Hct", ] - v["Het", ])
    out <- c(out, setNames(c(tot, cis, tra), paste(rep(c("total", "cis", "trans"), each = 4), colnames(q), sc, sep = ".")))
  }
  out[SIGN.BOOT.STATS]
}

## sign_boot_cell: one cell of P1, the type I error and power of the bootstrap sign test. Row i of grid gives the cells per group N (each parent N
## cells, the hybrid N cells split half cis and half trans by depth, as split_indices_by_depth() does), the mean count per parental cell (mean), the
## dispersion phi, an effect ratio and its target ("none": every group has phi, every true contrast is zero; "total": Pc has ratio x phi; "cis": the Sc
## allele has ratio x phi, in the cis and the trans cells) and R, the datasets of the cell. Every group has the same per-genome rate, so with ratio 1 every
## quantity's true contrast is zero. Offsets are drawn from the pools (the observed depths); the hybrid alleles share an extrinsic factor of variance
## ext_share x the smaller allele phi and each allele carries the private variance that gives it its phi (phi = 0: Poisson counts). Each dataset fits the
## six groups with the pipeline's NB estimator, then B bootstrap replicates resample each parent and the undivided hybrid cells (a copy keeps its cis or
## trans membership) and refit every group, a failed refit counting as phi = 0. Following the proposal: a dataset with a failed observed fit has no
## contrast; the log2 estimate exists only with every phi > 0; the test reads the observed sign of the log2 contrast (+-Inf allowed), or of the raw one
## when the log2 value is undefined. Returns per statistic: the share of datasets with a contrast and with a log2 estimate, the rejection rates at 0.05 and
## 0.01 (p <= alpha) under the proposal and under the tie-splitting variant with binomial SEs, the coverage of the true contrast by the 95% percentile
## interval (undefined replicates dropped; NA where the truth is undefined, log2 at phi = 0), and the mean shares of tied and undefined replicates.
sign_boot_cell <- function(i, grid, pools, B, ext_share = 0.5, streams) with_stream(streams[[i]], {
  g <- grid[i, ]
  stopifnot(g$N >= 50, g$mean > 0, g$phi >= 0, g$ratio >= 1, g$target %in% c("none", "total", "cis"), g$R >= 10, B >= 20)
  N <- g$N; mu <- g$mean / median(c(pools$O_Pc, pools$O_Pe))
  phi <- c(Pc = g$phi, Pe = g$phi, Hc = g$phi, He = g$phi)
  if (g$target == "total") phi[["Pc"]] <- g$phi * g$ratio
  if (g$target == "cis")   phi[["Hc"]] <- g$phi * g$ratio
  s2d <- ext_share * min(phi[["Hc"]], phi[["He"]]); vq <- pmax((1 + phi[c("Hc", "He")]) / (1 + s2d) - 1, 0)
  ## true statistics: per-genome rate mu in every group (the hybrid allele offset is half the cell's depth, so the allele rate is the per-genome rate)
  tf <- cbind(mu = mu, phi = c(Pc = phi[["Pc"]], Pe = phi[["Pe"]], Hcc = phi[["Hc"]], Hec = phi[["He"]], Hct = phi[["Hc"]], Het = phi[["He"]]))
  truth <- sign_boot_stats(tf)
  fit6 <- function(y, e, ix) { f <- rbind(Pc = fit_dispersion(y$Pc[ix$Pc], e$Pc[ix$Pc])[c("mu", "disp")], Pe = fit_dispersion(y$Pe[ix$Pe], e$Pe[ix$Pe])[c("mu", "disp")],
                                          Hcc = fit_dispersion(y$Hc[ix$c], e$H[ix$c])[c("mu", "disp")], Hec = fit_dispersion(y$He[ix$c], e$H[ix$c])[c("mu", "disp")],
                                          Hct = fit_dispersion(y$Hc[ix$t], e$H[ix$t])[c("mu", "disp")], Het = fit_dispersion(y$He[ix$t], e$H[ix$t])[c("mu", "disp")])
                               colnames(f) <- c("mu", "phi"); f }
  ns <- length(SIGN.BOOT.STATS)
  res <- array(NA_real_, c(ns, g$R, 7), dimnames = list(SIGN.BOOT.STATS, NULL, c("tested", "log2_est", "p_null", "p_split", "cover", "ties", "undef")))
  for (r in seq_len(g$R)) {
    e <- list(Pc = sample(pools$O_Pc, N, replace = TRUE), Pe = sample(pools$O_Pe, N, replace = TRUE), H = sort(sample(pools$O_H, N, replace = TRUE), decreasing = TRUE))
    sp <- split_indices_by_depth(N, 0.5)
    d <- sim_pos_factor(N, "gamma", 1, s2d)
    y <- list(Pc = rpois(N, mu * e$Pc * sim_pos_factor(N, "gamma", 1, phi[["Pc"]])), Pe = rpois(N, mu * e$Pe * sim_pos_factor(N, "gamma", 1, phi[["Pe"]])),
              Hc = rpois(N, mu * e$H * d * sim_pos_factor(N, "gamma", 1, vq[["Hc"]])), He = rpois(N, mu * e$H * d * sim_pos_factor(N, "gamma", 1, vq[["He"]])))
    fo <- fit6(y, e, list(Pc = seq_len(N), Pe = seq_len(N), c = sp$c, t = sp$t))
    tested <- all(is.finite(fo[, "phi"]))
    obs <- sign_boot_stats(fo)
    xb <- vapply(seq_len(B), function(b) { h <- sample.int(N, N, replace = TRUE)
      f <- fit6(y, e, list(Pc = sample.int(N, N, replace = TRUE), Pe = sample.int(N, N, replace = TRUE), c = h[h %in% sp$c], t = h[h %in% sp$t]))
      f[!is.finite(f[, "phi"]), "phi"] <- 0; sign_boot_stats(f) }, numeric(ns))
    for (k in seq_len(ns)) {
      st <- SIGN.BOOT.STATS[k]; x <- xb[k, ]
      raw_twin <- sub("\\.log2$", ".raw", st)
      s <- if (is.nan(obs[[k]]) || is.na(obs[[k]])) sign(obs[[raw_twin]]) else sign(obs[[k]])
      ci <- if (any(!is.na(x))) quantile(x, c(0.025, 0.975), type = 1, na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
      res[k, r, ] <- c(tested, tested && is.finite(obs[[k]]),
                       if (tested) sign_boot_p(x, s, "null") else NA_real_, if (tested) sign_boot_p(x, s, "split") else NA_real_,
                       if (tested && is.finite(truth[[k]]) && all(!is.na(ci))) as.numeric(ci[1] <= truth[[k]] && truth[[k]] <= ci[2]) else NA_real_,
                       mean(x == 0, na.rm = TRUE), mean(is.na(x)))
    }
  }
  do.call(rbind, lapply(SIGN.BOOT.STATS, function(st) {
    x <- res[st, , ]; t <- x[, "tested"] == 1; nt <- sum(t); parts <- strsplit(st, ".", fixed = TRUE)[[1]]
    rj <- function(col, a) if (nt) mean(x[t, col] <= a) else NA_real_
    data.frame(g, contrast = parts[1], quantity = parts[2], scale = parts[3], truth = truth[[st]], tested = mean(t), log2_estimate = mean(x[, "log2_est"] == 1),
               reject_05 = rj("p_null", 0.05), reject_05_se = binom_se(rj("p_null", 0.05), nt), reject_01 = rj("p_null", 0.01), reject_01_se = binom_se(rj("p_null", 0.01), nt),
               reject_05_split = rj("p_split", 0.05), reject_01_split = rj("p_split", 0.01), coverage = mean(x[, "cover"], na.rm = TRUE),
               tie_share = mean(x[, "ties"], na.rm = TRUE), undefined_share = mean(x[, "undef"], na.rm = TRUE), row.names = NULL)
  }))
})
