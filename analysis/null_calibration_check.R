###############################################################
### null_calibration_check.R
### Local calibration check of the null distributions (Methods: Null distributions; Notes for revision). Synthetic genes are simulated with every allele
### dispersion built from a factor shared by the two alleles of a cell (variance SIGMA.D) and a private one, transcriptome share variation, and three designs:
###   equal   every group shares one per-genome mean and one dispersion: every mean, dispersion and kinetic null is true;
###   level   the hybrid alleles are the parental means times a common level factor (K.LEVEL), with different parental means: the trans mean null is true,
###           which tests the level-adjusted trans permutation;
###   means   every group shares one dispersion but the means differ between groups: every dispersion null is true on both scales.
### Each gene is tested with null_contrasts_one(); under a true null the p-values are uniform, so the share below 0.05 and below 0.20 should be near 0.05 and 0.20.
### Writes <TABLE.DIR>/null_calibration_check.csv when TABLE.DIR exists. Needs R/functions.R; the defaults take tens of minutes on a few cores.
###############################################################
suppressWarnings(suppressMessages({ if (!exists("null_contrasts_one")) source("R/functions.R"); library(parallel) }))
set.seed(5)
NG <- 60; N.SIM <- 150; N.PERM <- 400; NH <- 900; NP <- 600; SIGMA.D <- 0.05; K.LEVEL <- 1.4; CORES <- max(1, detectCores() - 1)
lib <- rlnorm(NH, log(2e4), 0.35); S <- lib / 2e4; ci <- rbeta(NH, 0.5024 * 4e3, 0.4976 * 4e3); eP <- rlnorm(NP, 0, 0.3)
## counts of one gene: per-genome means (Pc, Pe, Hc, He) and one dispersion phi for every group; the alleles share the cell's extrinsic factor
gene_counts <- function(mu, phi) {
  d <- rgamma1(NH, SIGMA.D); vq <- max((1 + phi) / (1 + SIGMA.D) - 1, 0)
  list(MIX.SC = rpois(NP, eP * mu[1] * rgamma1(NP, phi)), MIX.SE = rpois(NP, eP * mu[2] * rgamma1(NP, phi)),
       HYB.SC = rpois(NH, S * ci * mu[3] * d * rgamma1(NH, vq)), HYB.SE = rpois(NH, S * (1 - ci) * mu[4] * d * rgamma1(NH, vq)))
}
design <- list()
for (i in seq_len(NG)) { m <- exp(rnorm(1, log(4), 0.5)); design[[length(design) + 1]] <- list(set = "equal", mu = rep(m, 4), phi = exp(rnorm(1, log(0.3), 0.3))) }
for (i in seq_len(NG)) { m <- exp(rnorm(1, log(4), 0.5)); r <- 2; design[[length(design) + 1]] <- list(set = "level", mu = c(m, m / r, K.LEVEL * m, K.LEVEL * m / r), phi = exp(rnorm(1, log(0.3), 0.3))) }
for (i in seq_len(NG)) { m <- exp(rnorm(1, log(4), 0.5)); design[[length(design) + 1]] <- list(set = "means", mu = m * c(1, 0.5, 1.6, 0.8), phi = exp(rnorm(1, log(0.3), 0.3))) }
gn <- sprintf("g%03d", seq_along(design))
cnt <- lapply(design, function(d) gene_counts(d$mu, d$phi))
mat <- function(k) { x <- t(sapply(cnt, `[[`, k)); rownames(x) <- gn; x }
Hsc <- mat("HYB.SC"); Hse <- mat("HYB.SE")
libh <- round(lib); chat <- rbinom(NH, libh, ci) / libh
st <- share_terms(chat, libh); SR <- share_r_groups(st$r_c, st$r_e); eH <- S / 2
pc <- ploidy_components(Hsc, Hse, eH, st)
ix <- list(MEAN = split_indices_by_depth(NH, 0.4), DISP = split_indices_by_depth(NH, 0.4), KBAL = split_indices_by_depth(NH, 0.4))
M <- list(MIX.SC = mat("MIX.SC"), MIX.SE = mat("MIX.SE"), HYB.SC = Hsc, HYB.SE = Hse, HYB.COMB = Hsc + Hse)
E <- list(MIX.SC = eP, MIX.SE = eP, HYB.SC = eH, HYB.SE = eH, HYB.COMB = S)
for (sp in names(ix)) for (side in c("c", "t")) for (al in c("SC", "SE")) {
  nm <- paste0("HY", toupper(side), ".", al, switch(sp, MEAN = "", DISP = ".N", KBAL = ".K"))
  M[[nm]] <- (if (al == "SC") Hsc else Hse)[, ix[[sp]][[side]]]; E[[nm]] <- eH[ix[[sp]][[side]]]
}
FITS <- lapply(names(M), function(k) { f <- fit_counts_offset(M[[k]], E[[k]], share_r = SR[[k]]); data.frame(MU = f$MU, DISP = f$DISP, row.names = rownames(f)) })
names(FITS) <- names(M)
share_obs <- list(cbar = st$cbar, r = SR, sigma_d = setNames(pc$sigma_d, gn))
perms <- with_local_seed(2, lapply(seq_len(N.PERM), perm_label_draw, hyc = ix$MEAN$c, hyt = ix$MEAN$t, nHYB = NH, nSC = NP, nSE = NP))
PR.CAL <- do.call(rbind, mclapply(gn, null_contrasts_one, mc.cores = CORES, expos = E, fits = FITS, mats = M, perms = perms, share_obs = share_obs, n_sim = N.SIM, seed = 3))
PR.CAL$set <- vapply(design, `[[`, "", "set")
## the tests whose null is true in each design
TRUE.NULLS <- list(equal = grep("^(mean|disp|rdisp|bsize|kbal|kbalsup)_(total|cis|trans|mid|par_sc|par_se)_p$", names(PR.CAL), value = TRUE),
                   level = "mean_trans_p",
                   means = grep("^(disp|rdisp)_(total|cis|trans|mid|par_sc|par_se)_p$", names(PR.CAL), value = TRUE))
NULL.CALIBRATION <- do.call(rbind, lapply(names(TRUE.NULLS), function(st_) do.call(rbind, lapply(TRUE.NULLS[[st_]], function(cn) {
  p <- PR.CAL[PR.CAL$set == st_, cn]; data.frame(design = st_, test = cn, genes = sum(!is.na(p)), below_05 = mean(p < 0.05, na.rm = TRUE), below_20 = mean(p < 0.20, na.rm = TRUE)) }))))
print(NULL.CALIBRATION, digits = 3, row.names = FALSE)
if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) write.csv(NULL.CALIBRATION, file.path(TABLE.DIR, "null_calibration_check.csv"), row.names = FALSE)
