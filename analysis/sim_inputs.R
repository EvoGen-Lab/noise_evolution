###############################################################
### sim_inputs.R
### Local script (seconds, no cluster). Estimates the inputs of the estimator simulation from the
### data and packages them, with the parameter grid, for cluster/scripts/estimator_sim.R.
###
### Needs   : MIX.SC, MIX.SE, HYB.SC, HYB.SE (genes x cells UMI matrices; the Section 1 checkpoint),
###           and R/functions.R sourced (.fit_one).
### Output  : <INPUT.DIR>/estimator_sim_inputs.rda
###
### Estimated from the data
###   library sizes, DEPTH.REF, cell numbers  -> the simulated depths and cell counts
###   C.BAR, VAR.C (sigma^2_c), VAR.CC (sigma^2_C|c)  -> the transcriptome share of the simulated hybrid
###   per-gene mean and dispersion quantiles -> the gene range the grid covers, and how many genes
###                                            sit near each grid value (the weights the decision uses)
###############################################################

if (!exists("INPUT.DIR")) INPUT.DIR <- "."
if (!exists(".fit_one")) source("R/functions.R")
if (!exists("HYB.SC")) load(ckpt_path(1))

## ---- settings -------------------------------------------------------------------------
N.BOOT.SHARE <- 1000                          # bootstrap resamples of cells for the share-variance intervals
SEED.SHARE   <- 1
NREP   <- 200                                 # simulated datasets per grid row
NBOOT  <- 100                                 # bootstrap resamples per simulated dataset
NTRUE  <- 5e7                                 # latent cells behind each truth (streamed in chunks; the Monte Carlo error of a log2 variance is about 0.01 for the heaviest-tailed family)
SEED.BASE <- 1
MU.PROBS  <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95)  # quantiles of the observed gene means that set the mean grid
PHI.PROBS <- c(0.05, 0.25, 0.5, 0.75, 0.95)     # quantiles of the observed dispersions that set the dispersion grid
OMEGA.GRID  <- c(0, 0.25, 0.5, 0.75, 1)       # intrinsic fraction
CBAR.GRID   <- c(0.3, 0.5, 0.7)               # average Sc transcriptome share
FRAC.GRID   <- c(1, 0.5, 0.25, 0.1)           # cell number as a fraction of the observed number
FAMILIES    <- c("gamma", "lognormal", "bimodal")
DIV.LOG2    <- c(dmu = 1, dphi = 1)           # Se-side shift in the cis and trans scenarios (log2)

## ---- depth: library sizes and cell numbers ---------------------------------------------
## The simulated depths are drawn from these observed vectors, so the long right tail of the
## library-size distribution enters the simulation as observed.
LIB.HYB   <- colSums(HYB.SC) + colSums(HYB.SE)          # one hybrid library holds both genomes
LIB.PAR   <- c(colSums(MIX.SC), colSums(MIX.SE))
DEPTH.REF <- median(c(LIB.HYB, LIB.PAR))
N.OBS     <- c(HYB = ncol(HYB.SC), MIX.SC = ncol(MIX.SC), MIX.SE = ncol(MIX.SE))
cat(sprintf("Library size: median %.0f, SD %.0f, max %.0f (hybrid median %.0f, parent median %.0f)\n",
            DEPTH.REF, sd(c(LIB.HYB, LIB.PAR)), max(c(LIB.HYB, LIB.PAR)), median(LIB.HYB), median(LIB.PAR)))
cat("Cells:", paste(names(N.OBS), N.OBS, collapse = ", "), "\n")

## ---- transcriptome share: c-bar, sigma^2_c, sigma^2_C|c ------------------------------------
## One calculated share per cell mixes true cell-to-cell variation with measurement error. Two
## log-odds shares from disjoint gene halves measure the share twice with independent errors, so
## their covariance across cells is the latent variance (on the log-odds scale) and the rest of the
## full-gene variance is measurement error. The logistic slope c-bar (1 - c-bar) maps both back to
## the share scale.
L.C   <- colSums(HYB.SC); L.E <- colSums(HYB.SE)
C.BAR <- mean(L.C / (L.C + L.E))

RANK  <- order(rowSums(HYB.SC) + rowSums(HYB.SE), decreasing = TRUE)   # alternate by expression rank: halves of similar depth
HALF1 <- RANK[seq(1, length(RANK), by = 2)]
HALF2 <- RANK[seq(2, length(RANK), by = 2)]
LO1    <- log(colSums(HYB.SC[HALF1, , drop = FALSE]) / colSums(HYB.SE[HALF1, , drop = FALSE]))
LO2    <- log(colSums(HYB.SC[HALF2, , drop = FALSE]) / colSums(HYB.SE[HALF2, , drop = FALSE]))
LO.ALL <- log(L.C / L.E)
OK     <- is.finite(LO1) & is.finite(LO2) & is.finite(LO.ALL)

VAR.LAMBDA <- cov(LO1[OK], LO2[OK])                  # sigma^2_ell, latent log-odds variance
VAR.ZETA   <- var(LO.ALL[OK]) - VAR.LAMBDA           # measurement variance of the log-odds
SLOPE2 <- C.BAR^2 * (1 - C.BAR)^2
VAR.C  <- SLOPE2 * VAR.LAMBDA                        # sigma^2_c
VAR.CC <- SLOPE2 * VAR.ZETA                          # sigma^2_C|c

## Intervals: resample cells, indices drawn together.
set.seed(SEED.SHARE)
IDX  <- which(OK)
BOOT.IDX <- matrix(sample(IDX, length(IDX) * N.BOOT.SHARE, replace = TRUE), length(IDX), N.BOOT.SHARE)
BOOT <- matrix(NA_real_, 2, N.BOOT.SHARE, dimnames = list(c("VAR.C", "VAR.CC"), NULL))
for (b in seq_len(N.BOOT.SHARE)) {
  lam <- cov(LO1[BOOT.IDX[, b]], LO2[BOOT.IDX[, b]])
  BOOT[, b] <- SLOPE2 * c(lam, var(LO.ALL[BOOT.IDX[, b]]) - lam)
}
FLOOR.CC <- C.BAR * (1 - C.BAR) * mean(1 / LIB.HYB)  # binomial sampling floor of sigma^2_C|c at these library sizes
SHARE.TABLE <- data.frame(estimate = c(VAR.C, VAR.CC, FLOOR.CC),
                          lo = c(apply(BOOT, 1, quantile, 0.025), NA), hi = c(apply(BOOT, 1, quantile, 0.975), NA),
                          row.names = c("sigma2_c", "sigma2_C|c", "binomial floor of sigma2_C|c"))
cat(sprintf("c-bar = %.4f\n", C.BAR)); print(signif(SHARE.TABLE, 3))
if (VAR.C <= 0 || VAR.CC <= 0) warning("a share variance estimate is not positive; the gene halves may not be independent, so check before using it")
if (VAR.CC < FLOOR.CC) warning("sigma^2_C|c is below its binomial floor; the simulation raises it to the floor")

## ---- gene range: means and dispersions of the parental fits ----------------------------
## mu is the mean count at median depth (sum y / sum S) and phi = 1 / theta from the pipeline's
## single estimator. The grid values are quantiles of the genes that have a finite phi, and each
## gene is assigned to its nearest grid value on the log scale to give the grid's weights.
FITS <- rbind(
  t(apply(MIX.SC, 1, .fit_one, expo = colSums(MIX.SC) / DEPTH.REF)),
  t(apply(MIX.SE, 1, .fit_one, expo = colSums(MIX.SE) / DEPTH.REF)))
GOOD   <- is.finite(FITS[, "mu"]) & FITS[, "mu"] > 0
MU.HAT <- FITS[GOOD, "mu"]
PHI.HAT <- 1 / FITS[GOOD & is.finite(FITS[, "disp"]) & FITS[, "disp"] > 0, "disp"]
MU.GRID  <- signif(quantile(MU.HAT, MU.PROBS, names = FALSE), 2)
PHI.GRID <- signif(quantile(PHI.HAT, PHI.PROBS, names = FALSE), 2)
MU.W  <- tabulate(max.col(-abs(outer(log(MU.HAT), log(MU.GRID), "-"))), length(MU.GRID)) / length(MU.HAT)
PHI.W <- tabulate(max.col(-abs(outer(log(PHI.HAT), log(PHI.GRID), "-"))), length(PHI.GRID)) / length(PHI.HAT)
cat(sprintf("Genes with a finite dispersion: %d of %d fits\n", length(PHI.HAT), length(MU.HAT)))
cat("Mean grid (per cell at median depth):", MU.GRID, "\n  weights:", round(MU.W, 3), "\n")
cat("Dispersion grid (phi):", PHI.GRID, "\n  weights:", round(PHI.W, 3), "\n")

## ---- grid ---------------------------------------------------------------------------------
## Main grid: no divergence between the species (every group shares mu and phi). Scenario grid:
## cis (hybrid alleles and parents differ) and trans (parents differ, hybrid alleles equal) on a
## reduced set of values. phi and mu of the Se side are shifted by DIV.LOG2. omega = 0 is left out of
## the scenario grid: a part with no intrinsic mean cannot carry the extra private variance of the Se side.
MAIN <- expand.grid(mu = MU.GRID, phi = PHI.GRID, omega = OMEGA.GRID, cbar = CBAR.GRID, frac = FRAC.GRID,
                    family = FAMILIES, scenario = "null", stringsAsFactors = FALSE)
MAIN$dmu <- 0; MAIN$dphi <- 0
SCEN <- expand.grid(mu = MU.GRID[unique(round(seq(1, length(MU.GRID), length.out = 3)))],      # three values spread over each grid
                    phi = PHI.GRID[unique(round(seq(1, length(PHI.GRID), length.out = 3)))], omega = c(0.25, 0.5, 1), cbar = c(0.5, 0.7), frac = 1,
                    family = FAMILIES, scenario = c("cis", "trans"), stringsAsFactors = FALSE)
SCEN$dmu <- DIV.LOG2[["dmu"]]; SCEN$dphi <- DIV.LOG2[["dphi"]]
GRID <- rbind(MAIN, SCEN)
## Weight of a row in the decision: the share of observed genes near its mean and dispersion values.
GRID$w <- MU.W[match(GRID$mu, MU.GRID)] * PHI.W[match(GRID$phi, PHI.GRID)]
## Rows with the same family, phi, omega, cbar and scenario share one latent truth (it does not depend on mu or the cell number).
TKEY <- paste(GRID$family, GRID$phi, GRID$omega, GRID$cbar, GRID$scenario, GRID$dmu, GRID$dphi)
GRID$tkey <- match(TKEY, unique(TKEY))
GRID$i <- seq_len(nrow(GRID))
cat(sprintf("Grid: %d rows (%d main, %d scenario), %d distinct truths\n", nrow(GRID), nrow(MAIN), nrow(SCEN), max(GRID$tkey)))

save(GRID, LIB.HYB, LIB.PAR, DEPTH.REF, N.OBS, C.BAR, VAR.C, VAR.CC, SHARE.TABLE, MU.GRID, PHI.GRID, MU.W, PHI.W,
     NREP, NBOOT, NTRUE, SEED.BASE, file = file.path(INPUT.DIR, "estimator_sim_inputs.rda"))
