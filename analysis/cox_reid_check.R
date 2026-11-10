###############################################################
### cox_reid_check.R
### Local unit checks of the Cox-Reid adjusted dispersion fit (.fit_one(), fit_dispersion()) and of dispersion_bound_stats(): the adjusted slope at zero
### matches its closed form and the boundary rule, the fit recovers phi and flags a maximum at the end of the search range, the share corrections lower
### the fitted phi by the expected amount, simulated Poisson genes pass the boundary likelihood-ratio test at the nominal rate, and the one-sided 95%
### likelihood bounds cover the true phi at about 95% each. Also checks the matrix fitter on one core and on a PSOCK cluster.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(1)

## the adjusted slope at phi = 0+ equals one half of sum((y - m)^2) - sum(m) + sum(m^2) / sum(m), and reduces to the (N - 1) m rule at equal offsets
e <- rlnorm(200, 0, 0.5); y <- rnbinom(200, mu = 2 * e, size = 3); m <- sum(y) / sum(e) * e
ll <- .cr_loglik(y, m); h <- 1e-4
slope_num <- (4 * ll(h) - ll(2 * h) - 3 * ll(0)) / (2 * h)
slope_cf  <- 0.5 * (sum((y - m)^2) - sum(m) + sum(m^2) / sum(m))
stopifnot(abs(slope_num / slope_cf - 1) < 1e-3)
y2 <- rpois(100, 3); m2 <- rep(mean(y2), 100)
stopifnot(abs((sum((y2 - m2)^2) - sum(m2) + sum(m2^2) / sum(m2)) - (sum((y2 - m2)^2) - 99 * mean(y2))) < 1e-8)
## a non-positive slope puts the fit on the boundary, a positive one does not
stopifnot(fit_dispersion(y2, rep(1, 100))[["boundary"]] == as.numeric(sum((y2 - m2)^2) <= 99 * mean(y2)))

## the fit recovers phi, its maximum matches a grid search, and a maximum at the end of the range is flagged
n <- 300; e <- rlnorm(n, 0, 0.5)
est <- replicate(300, fit_dispersion(rnbinom(n, mu = 2 * e, size = 1 / 0.3), e)[["disp"]])
stopifnot(abs(mean(est) / 0.3 - 1) < 0.05)
y <- rnbinom(n, mu = 3 * e, size = 1 / 0.2); llg <- .cr_loglik(y, sum(y) / sum(e) * e); grid <- exp(seq(-8, 3, length.out = 4000))
stopifnot(abs(fit_dispersion(y, e)[["disp_app"]] / grid[which.max(vapply(grid, llg, numeric(1)))] - 1) < 0.02)
stopifnot(fit_dispersion(rnbinom(n, mu = 5 * e, size = 1e-3), e)[["hit_range"]] == 1, fit_dispersion(y, e)[["hit_range"]] == 0)

## the share corrections: (phi - r) / (1 + r), then the additive term, then truncation at zero
a <- fit_dispersion(y, e)[["disp"]]
stopifnot(abs(fit_dispersion(y, e, 0.01)[["disp"]] - (a - 0.01) / 1.01) < 1e-12, abs(fit_dispersion(y, e, 0, "NB", 0.01)[["disp"]] - (a - 0.01)) < 1e-12)
stopifnot(fit_dispersion(y, e, 0, "NB", 10)[["disp"]] == 0, fit_dispersion(y, e, 0, "NB", 10)[["boundary"]] == 1)

## the critical value of the likelihood bounds is the 0.90 quantile of chi-square(1), 2.7055, which makes each bound one-sided at 95%
stopifnot(round(BOUND.CRIT, 4) == 2.7055)

## Poisson genes: the likelihood-ratio statistic is zero about half the time and exceeds 2.71 about 5% of the time
lrt <- replicate(1000, dispersion_bound_stats(rpois(n, 5 * e), e)[["lrt"]])
stopifnot(abs(mean(lrt > qchisq(0.9, 1)) - 0.05) < 0.02, abs(mean(lrt == 0) - 0.5) < 0.06)
## NB genes: each one-sided bound covers the true phi about 95% of the time, the lower bound is zero whenever the LRT does not exceed the critical value, and
## the share correction is applied to both bounds
bd <- t(replicate(400, dispersion_bound_stats(rnbinom(n, mu = 3 * e, size = 1 / 0.3), e)))
stopifnot(abs(mean(bd[, "disp_upper"] >= 0.3) - 0.95) < 0.04, abs(mean(bd[, "disp_lower"] <= 0.3) - 0.95) < 0.04, all(bd[, "disp_lower"] <= bd[, "disp_upper"]))
st <- dispersion_bound_stats(rpois(n, 1 * e), e); stopifnot(st[["lrt"]] >= 0, st[["disp_upper"]] >= 0, st[["lrt"]] > BOUND.CRIT || st[["disp_lower"]] == 0)
y3 <- rnbinom(n, mu = 3 * e, size = 1 / 0.3); b0 <- dispersion_bound_stats(y3, e); b1 <- dispersion_bound_stats(y3, e, share_r = 0.01)
stopifnot(abs(b1[["disp_upper"]] - (b0[["disp_upper"]] - 0.01) / 1.01) < 1e-6, abs(b1[["disp_lower"]] - max((b0[["disp_lower"]] - 0.01) / 1.01, 0)) < 1e-6)

## the matrix fitter agrees on one core and on a PSOCK cluster, with a per-gene additive correction
mat <- t(sapply(1:12, function(i) rnbinom(n, mu = (0.5 + i / 3) * e, size = 1 / (0.05 * i))))
rownames(mat) <- paste0("g", 1:12)
f1 <- fit_counts_offset(mat, e, share_r = 1e-4, share_sub = seq(0, 0.011, length.out = 12))
cl <- parallel::makeCluster(2)
f2 <- fit_counts_offset(mat, e, cl = cl, share_r = 1e-4, share_sub = seq(0, 0.011, length.out = 12)); parallel::stopCluster(cl)
stopifnot(isTRUE(all.equal(f1, f2)), all(c("LRT", "DISP_LOWER", "DISP_UPPER", "MBAR", "HIT_RANGE", "DISP_APP") %in% names(f1)), !any(grepl("CLASS", names(f1))))
cat("Cox-Reid fit checks passed\n")
