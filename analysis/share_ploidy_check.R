###############################################################
### share_ploidy_check.R
### Local simulation check of share_terms() and ploidy_components() against the model of Methods Appendices C and D. Hybrid cells have a latent transcriptome
### share c_i, a library size and, for each gene, per-genome rates, a shared extrinsic deviation d (variance sigma_d) and private intrinsic deviations q_c, q_e
### (variance sigma_q), all shifted exponentials with mean zero; counts are Poisson. With the true share variance the corrected estimators recover sigma_d,
### sigma_q and the ploidy shifts, the star estimators are biased upward as Appendix C derives, and the share_sub term equals the variance that share variation
### adds to the summed expression factor (Appendix D). Then the observed-cell estimator share_terms() recovers the share variance.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(7)
N <- 400000
cbar <- 0.5024; var_c <- 4e-4                                         # inflated share variance so the corrections are visible
k <- cbar * (1 - cbar) / var_c - 1; ci <- rbeta(N, cbar * k, (1 - cbar) * k)
S <- rlnorm(N, 0, 0.3) * 1; O <- S / 2                                 # allele offset S_H / 2
shexp <- function(v) sqrt(v) * (rexp(N) - 1)                           # shifted exponential, mean 0, variance v
genes <- data.frame(mu_c = c(3, 6, 12), mu_e = c(3, 2, 20), sd = c(0.10, 0.3, 0.05), sq = c(0.30, 0.1, 0.4))
Yc <- Ye <- matrix(0, nrow(genes), N); truth <- list()
for (g in seq_len(nrow(genes))) {
  d <- shexp(genes$sd[g]); ac <- 1 + d + shexp(genes$sq[g]); ae <- 1 + d + shexp(genes$sq[g])
  Yc[g, ] <- rpois(N, S * ci * genes$mu_c[g] * ac); Ye[g, ] <- rpois(N, S * (1 - ci) * genes$mu_e[g] * ae)
  pi <- cbar * genes$mu_c[g] / (cbar * genes$mu_c[g] + (1 - cbar) * genes$mu_e[g]); wsq <- pi^2 + (1 - pi)^2
  w <- pi / cbar - (1 - pi) / (1 - cbar)
  truth[[g]] <- c(pi = pi, delta = log2((genes$sd[g] + genes$sq[g]) / (genes$sd[g] + wsq * genes$sq[g])), delta_raw = (1 - wsq) * genes$sq[g],
                  share_sub = var_c * ((1 + genes$sd[g]) * w^2 + genes$sq[g] * (pi^2 / cbar^2 + (1 - pi)^2 / (1 - cbar)^2)),
                  var_atilde = var(pi * ci / cbar * ac + (1 - pi) * (1 - ci) / (1 - cbar) * ae) - var(pi * ac + (1 - pi) * ae))   # same draws with and without share variation
}
share_true <- list(cbar = cbar, var_c = var_c, r_c = var_c / cbar^2, r_e = var_c / (1 - cbar)^2, xi = var_c / (cbar * (1 - cbar)))
pc <- ploidy_components(Yc, Ye, rep(1, N) * O[1] * 0 + O, share_true)
tr <- do.call(rbind, truth)
stopifnot(all(abs(pc$sigma_d - genes$sd) < 0.01), all(abs(pc$sigma_q - genes$sq) < 0.01), all(abs(pc$pi - tr[, "pi"]) < 0.002),
          all(abs(pc$delta - tr[, "delta"]) < 0.02), all(abs(pc$delta_raw - tr[, "delta_raw"]) < 0.01))
## the summed-hybrid term equals the variance that share variation adds to the summed expression factor (compared on common random numbers)
stopifnot(all(abs(pc$share_sub - tr[, "var_atilde"]) < 0.0003), all(abs(pc$share_sub - tr[, "share_sub"]) < 1e-4))
## ignoring the share variation (zero share terms) leaves the star estimators, which are biased upward in sigma_q as Appendix C derives
pc0 <- ploidy_components(Yc, Ye, O, list(cbar = cbar, var_c = 0, r_c = 0, r_e = 0, xi = 0))
stopifnot(all(pc0$sigma_q > pc$sigma_q), all(pc0$share_sub == 0))
## share_terms() recovers the share variance from the observed shares: binomial noise is removed with L - 1
lib <- round(S * 20000); chat <- rbinom(N, lib, ci) / lib
st <- share_terms(chat, lib)
stopifnot(abs(st$var_c / var_c - 1) < 0.05, abs(st$cbar - cbar) < 1e-3, abs(st$offset - log2(cbar / (1 - cbar))) < 0.01, st$var_obs > st$var_binom)
## share_r_groups() names every dataset
r <- share_r_groups(0.001, 0.002); stopifnot(identical(names(r), names(.DRAW.KEY)), r[["HYC.SC.N"]] == 0.001, r[["HYT.SE"]] == 0.002, r[["HYB.COMB"]] == 0, r[["MIX.SC"]] == 0)
cat("share and ploidy checks passed\n")
