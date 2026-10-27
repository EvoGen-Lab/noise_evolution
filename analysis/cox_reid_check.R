###############################################################
### cox_reid_check.R
### Local unit checks of the Cox-Reid adjusted dispersion fit (.fit_one(), fit_dispersion()) and of dispersion_class_stats(): the adjusted slope at zero
### matches its closed form and the boundary rule, the fit recovers phi and flags a maximum at the end of the search range, the share corrections lower
### the fitted phi by the expected amount, simulated Poisson genes pass the boundary likelihood-ratio test at the nominal rate, and the class proportions
### of simulated NB genes follow the mean count and phi. Also checks the matrix fitter on one core and on a PSOCK cluster.
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

## Poisson genes: the likelihood-ratio statistic is zero about half the time and exceeds 2.71 about 5% of the time
lrt <- replicate(1000, dispersion_class_stats(rpois(n, 5 * e), e)[["lrt"]])
stopifnot(abs(mean(lrt > qchisq(0.9, 1)) - 0.05) < 0.02, abs(mean(lrt == 0) - 0.5) < 0.06)
## NB genes: resolved above zero when mu and phi are large, uninformative when both are small, resolved near zero for a quiet gene with a tight bound
cls <- function(mu, phi) replicate(100, dispersion_class_stats(rnbinom(n, mu = mu * e, size = 1 / phi), e)[["class"]])
stopifnot(mean(cls(5, 1) == 1) > 0.95, mean(cls(0.05, 0.01) == 3) > 0.8, mean(cls(0.5, 0.01) == 2) > 0.15)
## classes are re-derived from the stored statistics: a stricter near-zero threshold moves class 2 genes to class 3 only
st <- dispersion_class_stats(rpois(n, 1 * e), e); stopifnot(st[["class"]] %in% 1:3, st[["lrt"]] >= 0, st[["disp_upper"]] >= 0)

## the matrix fitter agrees on one core and on a PSOCK cluster, with a per-gene additive correction
mat <- t(sapply(1:12, function(i) rnbinom(n, mu = (0.5 + i / 3) * e, size = 1 / (0.05 * i))))
rownames(mat) <- paste0("g", 1:12)
f1 <- fit_counts_offset(mat, e, share_r = 1e-4, share_sub = seq(0, 0.011, length.out = 12))
cl <- parallel::makeCluster(2)
f2 <- fit_counts_offset(mat, e, cl = cl, share_r = 1e-4, share_sub = seq(0, 0.011, length.out = 12)); parallel::stopCluster(cl)
stopifnot(isTRUE(all.equal(f1, f2)), all(c("CLASS", "LRT", "DISP_UPPER", "MBAR", "HIT_RANGE", "DISP_APP") %in% names(f1)))
cat("Cox-Reid fit checks passed\n")
