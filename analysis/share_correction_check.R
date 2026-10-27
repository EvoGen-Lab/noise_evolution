###############################################################
### share_correction_check.R
### Local simulation check of the share term in a hybrid allele's fitted dispersion. Hybrid allele counts are drawn with a known dispersion phi
### and a known cell-to-cell variance of the Sc share, offset by half the cell's hybrid depth as in analysis.R Section 2, and fitted with
### .fit_one(). A mean-1 share factor with variance r = var(c) / c-bar^2 (Hc) or var(c) / (1 - c-bar)^2 (He) raises the apparent dispersion to
### phi + (1 + phi) r, and the fits keep that value; the check compares the fitted median with the prediction.
###
### Run from the project root (R/functions.R is sourced when .fit_one() is not loaded). Uses the hybrid library sizes LIB.HYB
### when they are in the session (after Section 1.2) and a log-normal stand-in otherwise. The scenarios are the share variance
### of the data (VAR.C when it exists, else the value below) and two inflated ones, because at the observed variance the
### term is far below the sampling error of any fit and the check would show nothing.
### Writes <TABLE.DIR>/share_correction_check.csv when TABLE.DIR exists.
###############################################################

if (!exists(".fit_one")) source("R/functions.R")
set.seed(1)
NREP <- 100                                                  # genes simulated per cell of the grid
MU.GRID  <- c(0.3, 1, 3, 10, 30)                             # mean count per unit offset
PHI.GRID <- c(0.001, 0.01, 0.05, 0.2, 0.5, 1)                # true dispersion
C.BAR.SIM <- if (exists("C.BAR")) C.BAR else 0.5024
VAR.C.OBS <- if (exists("VAR.C")) VAR.C else 1.8627e-5       # variance of the true share in the data
SCENARIO  <- c(observed = VAR.C.OBS, inflated_sd_0.02 = 0.02^2, inflated_sd_0.05 = 0.05^2)
LIB <- if (exists("LIB.HYB")) LIB.HYB else rlnorm(1400, log(3e4), 0.7)
D <- median(LIB); S.H <- LIB / D                             # hybrid depth relative to the reference; each allele is offset by S.H / 2
n <- length(LIB)

RES <- list()
for (sc in names(SCENARIO)) {
  v <- SCENARIO[[sc]]
  r.c <- v / C.BAR.SIM^2; r.e <- v / (1 - C.BAR.SIM)^2
  ## Beta share with mean c-bar and variance v (a point mass when v = 0)
  k <- C.BAR.SIM * (1 - C.BAR.SIM) / v - 1
  for (mu in MU.GRID) for (phi in PHI.GRID) {
    est <- matrix(NA_real_, NREP, 2, dimnames = list(NULL, c("Hc", "He")))
    for (g in seq_len(NREP)) {
      cc <- if (v > 0) rbeta(n, C.BAR.SIM * k, (1 - C.BAR.SIM) * k) else rep(C.BAR.SIM, n)   # latent share of each cell
      a_c <- rgamma(n, 1 / phi, scale = phi); a_e <- rgamma(n, 1 / phi, scale = phi)          # gene factors, mean 1, variance phi
      y_c <- rpois(n, S.H / 2 * mu * (cc / C.BAR.SIM) * a_c)                                  # Hc: share factor c / c-bar
      y_e <- rpois(n, S.H / 2 * mu * ((1 - cc) / (1 - C.BAR.SIM)) * a_e)                      # He: share factor (1 - c) / (1 - c-bar)
      est[g, ] <- c(.fit_one(y_c, S.H / 2)[["disp"]], .fit_one(y_e, S.H / 2)[["disp"]])
    }
    pos <- function(x) x[is.finite(x) & x > 0]               # fits with a positive dispersion
    fit <- c(est[, "Hc"], est[, "He"]); r.mean <- (r.c + r.e) / 2
    RES[[length(RES) + 1]] <- data.frame(
      scenario = sc, r_c = r.c, mu = mu, phi_true = phi,
      med_fitted = median(pos(fit)), predicted = phi + (1 + phi) * r.mean,
      log2_bias_fitted = median(log2(pos(fit) / phi)), log2_bias_predicted = log2((phi + (1 + phi) * r.mean) / phi),
      zero_fitted = mean(!(fit > 0 & is.finite(fit))), row.names = NULL)
  }
}
SHARE.CHECK <- do.call(rbind, RES)
cat("Hybrid-allele dispersion fitted without a share term (medians over genes; log2_bias is log2(estimate / truth); predicted is phi + (1 + phi) r):\n")
print(SHARE.CHECK, digits = 3, row.names = FALSE)
cat("\nMedian absolute log2 bias of the fits over the cells with phi_true >= 0.05 and mu >= 1 (where a fit has information), with the predicted bias:\n")
print(aggregate(cbind(fitted = abs(log2_bias_fitted), predicted = abs(log2_bias_predicted)) ~ scenario, data = SHARE.CHECK[SHARE.CHECK$phi_true >= 0.05 & SHARE.CHECK$mu >= 1, ], FUN = median), digits = 3, row.names = FALSE)
if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) write.csv(SHARE.CHECK, file.path(TABLE.DIR, "share_correction_check.csv"), row.names = FALSE)
