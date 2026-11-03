###############################################################
### fit_dispersion_check.R
### Local unit check of fit_dispersion(): the weighted method of moments of Appendix A recovers phi without systematic bias under
### uneven exposures, its allele-share correction recovers the gene's phi, a converged estimate is a fixed point of the Appendix A update
### and sparse counts return a recorded convergence status, and the NB wrapper .fit_one() keeps its c(mu, disp, theta) output.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(11)

## bias: mean of the corrected phi over replicates, for a few (n, mu, phi) with lognormal exposures
GRID <- data.frame(n = c(100, 300, 40), mu = c(1, 2, 10), phi = c(0.5, 0.2, 0.3))
for (i in seq_len(nrow(GRID))) {
  est <- replicate(300, { e <- rlnorm(GRID$n[i], 0, 0.4); y <- rnbinom(GRID$n[i], mu = GRID$mu[i] * e, size = 1 / GRID$phi[i])
                          fit_dispersion(y, e, 0, "weighted")[["disp"]] })
  stopifnot(!anyNA(est), abs(mean(est) / GRID$phi[i] - 1) < 0.1)
}

## share correction: counts carry a mean-1 share factor with variance r on top of the gene's phi
r <- 0.05; phi <- 0.3; n <- 300
est <- replicate(300, { y <- rpois(n, 3 * rgamma(n, 1 / r, 1 / r) * rgamma(n, 1 / phi, 1 / phi))
                        fit_dispersion(y, rep(1, n), r, "weighted")[["disp"]] })
stopifnot(abs(mean(est) / phi - 1) < 0.1)

## a converged estimate is a fixed point of the Appendix A update: weights 1 / (1 + phi M)^2 at the estimate, phi-check, then
## (1 + 1 / sum(y)) phi-check / (1 - h phi-check)
e <- rlnorm(400, 0, 0.6); y <- rnbinom(400, mu = 2 * e, size = 1 / 0.4); f <- fit_dispersion(y, e, 0, "weighted")
p <- e / sum(e); M <- sum(y) * p; A <- M * (1 - p); B <- M^2 * (1 - 2 * p) + p^2 * sum(M^2); h <- sum(p^2); w <- 1 / (1 + f[["disp"]] * M)^2
chk <- sum(w * ((y - M)^2 - A)) / sum(w * B)
stopifnot(f[["converged"]] == 1, abs((1 + 1 / sum(y)) * chk / (1 - h * chk) - f[["disp"]]) < 1e-6 * (1 + f[["disp"]]))

## sparse counts: a finite estimate with its convergence status recorded
y <- c(rep(0, 43), rep(1, 5), 2, 3); e <- rep(1, 50)
f <- fit_dispersion(y, e, 0, "weighted")
stopifnot(f[["converged"]] %in% c(0, 1), is.finite(f[["disp"]]))

## no excess dispersion: the boundary is reported as a valid fit with the signed estimate kept
f <- fit_dispersion(rpois(200, 5), rep(1, 200), 0, "weighted")
stopifnot(f[["boundary"]] == 1 || f[["disp"]] > 0, is.finite(f[["disp_raw"]]))

## the NB wrapper keeps its three named outputs
stopifnot(identical(names(.fit_one(rnbinom(50, mu = 3, size = 2), rep(1, 50))), c("mu", "disp", "theta")))
cat("fit_dispersion checks passed\n")
