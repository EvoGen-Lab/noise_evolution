###############################################################
### share_ploidy_check.R
### Local simulation check of split_eta() and ploidy_components() against the model of Methods Appendix D with the library split of Library split correction.
### Hybrid cells have a latent transcriptome share s_i, a hybrid library L_H (small, so the split error is visible) and an observed Sc library
### L_Hc ~ Binomial(L_H, s_i); each allele is offset by its own library over D. Each gene has per-genome rates, a shared extrinsic deviation v (variance
### sigma_v) and private intrinsic deviations q_c, q_e (variance sigma_q), all shifted exponentials with mean zero; counts are Poisson. With the split variances
### of split_eta() the corrected estimators recover sigma_v, sigma_q, pi and the ploidy shifts Lambda and Lambda_raw; without them (eta = 0, the star
### estimators) sigma_q is biased upward, as Appendix D derives.
### Needs   : R/functions.R. Writes nothing; stops with an error when a check fails.
###############################################################
suppressWarnings(suppressMessages(source("R/functions.R")))
set.seed(7)
N <- 400000; D <- 300
sbar <- 0.5024; var_s <- 4e-4
k <- sbar * (1 - sbar) / var_s - 1; si <- rbeta(N, sbar * k, (1 - sbar) * k)
L <- pmax(rpois(N, 300 * rlnorm(N, 0, 0.3)), 20)                      # small hybrid libraries: eta about 1 / L, large enough to resolve
Lc <- pmin(pmax(rbinom(N, L, si), 1), L - 1); shat <- Lc / L
eta <- split_eta(shat, L); oc <- Lc / D; oe <- (L - Lc) / D            # allele offsets L_Hc / D and L_He / D
stopifnot(all(abs(eta$Hc - (1 - shat) / (shat * L)) < 1e-15), all(abs(eta$x - (eta$Hc + eta$He - 1 / L)) < 1e-15))
shexp <- function(v) sqrt(v) * (rexp(N) - 1)                           # shifted exponential, mean 0, variance v
genes <- data.frame(mu_c = c(3, 6, 12), mu_e = c(3, 2, 20), sv = c(0.10, 0.3, 0.05), sq = c(0.30, 0.1, 0.4))
Yc <- Ye <- matrix(0, nrow(genes), N); truth <- list()
for (g in seq_len(nrow(genes))) {
  v <- shexp(genes$sv[g]); ac <- 1 + v + shexp(genes$sq[g]); ae <- 1 + v + shexp(genes$sq[g])
  Yc[g, ] <- rpois(N, L / D * si * genes$mu_c[g] * ac); Ye[g, ] <- rpois(N, L / D * (1 - si) * genes$mu_e[g] * ae)
  pi <- sum(L * si) * genes$mu_c[g] / (sum(L * si) * genes$mu_c[g] + sum(L * (1 - si)) * genes$mu_e[g]); wsq <- pi^2 + (1 - pi)^2
  truth[[g]] <- c(pi = pi, Lambda = log2((genes$sv[g] + genes$sq[g]) / (genes$sv[g] + wsq * genes$sq[g])), Lambda_raw = (1 - wsq) * genes$sq[g])
}
tr <- do.call(rbind, truth)
pc <- ploidy_components(Yc, Ye, oc, oe, eta)
stopifnot(all(abs(pc$sigma_v - genes$sv) < 0.01), all(abs(pc$sigma_q - genes$sq) < 0.01), all(abs(pc$pi - tr[, "pi"]) < 0.002),
          all(abs(pc$Lambda - tr[, "Lambda"]) < 0.02), all(abs(pc$Lambda_raw - tr[, "Lambda_raw"]) < 0.01), all(pc$good))
## without the split variances the star estimators keep the split error, which raises sigma_q (Appendix D)
pc0 <- ploidy_components(Yc, Ye, oc, oe, list(Hc = numeric(N), He = numeric(N), x = numeric(N)))
stopifnot(all(pc0$sigma_q - pc$sigma_q > 0.005), all(abs(pc0$sigma_q - genes$sq) > abs(pc$sigma_q - genes$sq)))
cat("split and ploidy checks passed\n")
