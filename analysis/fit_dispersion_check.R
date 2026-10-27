###############################################################
### fit_dispersion_check.R
### Local unit check of fit_dispersion(): the weighted method of moments recovers phi without systematic bias under uneven
### exposures, a mean-1 share factor adds r (1 + phi) to the fitted phi, it finds a solution where the plain weight iteration
### cycles (sparse counts), and the NB wrapper .fit_one() keeps its c(mu, disp, theta) output.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(11)

## bias: mean of the fitted phi over replicates, for a few (n, mu, phi) with lognormal exposures
GRID <- data.frame(n = c(100, 300, 40), mu = c(1, 2, 10), phi = c(0.5, 0.2, 0.3))
for (i in seq_len(nrow(GRID))) {
  est <- replicate(300, { e <- rlnorm(GRID$n[i], 0, 0.4); y <- rnbinom(GRID$n[i], mu = GRID$mu[i] * e, size = 1 / GRID$phi[i])
                          fit_dispersion(y, e, "weighted")[["disp"]] })
  stopifnot(!anyNA(est), abs(mean(est) / GRID$phi[i] - 1) < 0.1)
}

## share factor: counts carry a mean-1 factor with variance r on top of the gene's phi, so the fitted dispersion is the apparent value phi + (1 + phi) r
r <- 0.05; phi <- 0.3; n <- 300
est <- replicate(300, { y <- rpois(n, 3 * rgamma(n, 1 / r, 1 / r) * rgamma(n, 1 / phi, 1 / phi))
                        fit_dispersion(y, rep(1, n), "weighted")[["disp"]] })
stopifnot(abs(mean(est) / (phi + (1 + phi) * r) - 1) < 0.1)

## sparse counts where the fixed-point iteration cycles between 0 and 1.21: the root solve returns a positive phi
y <- c(rep(0, 43), rep(1, 5), 2, 3); e <- rep(1, 50)
f <- fit_dispersion(y, e, "weighted")
stopifnot(f[["converged"]] == 1, is.finite(f[["disp"]]))

## no excess dispersion: the boundary is reported as a valid fit with the signed estimate kept
f <- fit_dispersion(rpois(200, 5), rep(1, 200), "weighted")
stopifnot(f[["boundary"]] == 1 || f[["disp"]] > 0, is.finite(f[["disp_raw"]]))

## the NB wrapper keeps its three named outputs
stopifnot(identical(names(.fit_one(rnbinom(50, mu = 3, size = 2), rep(1, 50))), c("mu", "disp", "theta")))
cat("fit_dispersion checks passed\n")
