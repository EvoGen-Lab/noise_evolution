###############################################################
### estimator_sim_summary.R
### Local. Reads the estimator_sim_output_<k>of<K>.rda files from cluster/scripts/estimator_sim.R,
### verifies that every grid row arrived, reports the study's outputs, and applies the decision rule.
###
### Needs   : R/functions.R (sourced here if wmean() is not loaded), <INPUT.DIR>/estimator_sim_inputs.rda (GRID) and <OUTPUT.DIR>/estimator_sim_output_*.rda
### Writes  : tables under TABLE.DIR (when it exists) and a printed decision.
###
### Decision rule (log2 phi scale, null-scenario rows, targets Hc, He and the summed hybrid H):
###   A cell is FLAGGED for an estimator when |bias| exceeds BIAS.FRAC x its bootstrap SE and the bias is
###   distinguishable from zero (|bias| > 2 x its replicate and truth Monte Carlo error).
###   The flagged share is the share of observed genes (the grid weights w) that fall in flagged cells.
###   KEEP the MLE        : its flagged share is below KEEP.BELOW over the non-gamma families.
###   ADOPT the moment est.: the MLE share is at least ADOPT.ABOVE, the moment estimator's share is under half of
###                          it, and its RMSE is at most RMSE.RATIO times the MLE's in the MLE's flagged cells.
###   otherwise KEEP the MLE and add the correction terms, and report the flagged regions.
###############################################################

if (!exists("INPUT.DIR"))  INPUT.DIR  <- "."
if (!exists("OUTPUT.DIR")) OUTPUT.DIR <- "."
if (!exists("wmean")) source("R/functions.R")   # wmean() lives in functions.R
BIAS.FRAC   <- 0.25      # bias must stay under this fraction of the bootstrap SE
KEEP.BELOW  <- 0.05      # flagged share of genes under which the MLE is kept as is
ADOPT.ABOVE <- 0.20      # flagged share of genes at which a switch is considered
RMSE.RATIO  <- 1.25      # largest acceptable RMSE of the moment estimator relative to the MLE

## ---- load and verify ------------------------------------------------------------------------
load(file.path(INPUT.DIR, "estimator_sim_inputs.rda"))                # GRID and the settings
PART.FILES <- list.files(OUTPUT.DIR, pattern = "^estimator_sim_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(PART.FILES) > 0)
EST.LIST <- AGREE.LIST <- vector("list", length(PART.FILES)); ROW.DONE <- logical(nrow(GRID))
for (k in seq_along(PART.FILES)) {
  e <- new.env(); load(PART.FILES[k], envir = e)
  EST.LIST[[k]] <- e$SIM.EST; AGREE.LIST[[k]] <- e$SIM.AGREE; ROW.DONE[e$ROWS.PART] <- TRUE
}
stopifnot(all(ROW.DONE))                                              # every array task must have finished
SIM.EST   <- do.call(rbind, EST.LIST); SIM.AGREE <- do.call(rbind, AGREE.LIST)
stopifnot(!anyDuplicated(SIM.AGREE$i), nrow(SIM.EST) == nrow(GRID) * 20)   # 10 targets x 2 estimators per row

KEYS <- c("i", "mu", "phi", "omega", "cbar", "frac", "family", "scenario", "w")
EST   <- merge(SIM.EST, GRID[, KEYS], by = "i")
AGREE <- merge(SIM.AGREE, GRID[, KEYS], by = "i")
EST$mc_se <- AGREE$truth_mc_se[match(EST$i, AGREE$i)]

## ---- flags -----------------------------------------------------------------------------------
## Bias against the bootstrap SE, and whether the bias is distinguishable from zero given the Monte Carlo
## error of the replicate mean and of the truth.
EST$ratio   <- abs(EST$bias) / EST$se_boot
EST$bias_se <- sqrt(EST$sd_emp^2 / EST$n_ok + EST$mc_se^2)
EST$usable  <- is.finite(EST$bias) & is.finite(EST$se_boot) & EST$n_ok >= 0.5 * NREP
EST$flag    <- EST$usable & EST$ratio > BIAS.FRAC & abs(EST$bias) > 2 * EST$bias_se
NULLROWS <- EST[EST$scenario == "null", ]

## ---- 1. bias, RMSE, failure, bootstrap SE and coverage by target, estimator and family ---------
SPLIT <- split(NULLROWS, list(NULLROWS$target, NULLROWS$est, NULLROWS$family), drop = TRUE)
TAB.ALL <- data.frame(
  target = vapply(SPLIT, function(d) d$target[1], ""), est = vapply(SPLIT, function(d) d$est[1], ""), family = vapply(SPLIT, function(d) d$family[1], ""),
  fail = vapply(SPLIT, function(d) wmean(d$fail, d$w), 0), bias = vapply(SPLIT, function(d) wmean(d$bias, d$w), 0),
  abs_bias = vapply(SPLIT, function(d) wmean(abs(d$bias), d$w), 0), rmse = vapply(SPLIT, function(d) wmean(d$rmse, d$w), 0),
  se_boot = vapply(SPLIT, function(d) wmean(d$se_boot, d$w), 0), sd_emp = vapply(SPLIT, function(d) wmean(d$sd_emp, d$w), 0),
  cover = vapply(SPLIT, function(d) wmean(d$cover, d$w), 0), flagged = vapply(SPLIT, function(d) wmean(d$flag[d$usable], d$w[d$usable]), 0),
  row.names = NULL)
cat("\n== Weighted bias, RMSE, failure rate, bootstrap SE, coverage (null scenario, log2 phi) ==\n")
print(TAB.ALL[order(TAB.ALL$target, TAB.ALL$est, TAB.ALL$family), ], digits = 3, row.names = FALSE)

## ---- 1b. the same by grid cell: mean x dispersion x cell fraction x family -------------------------
## Averaging over the whole grid mixes very low-count genes, where every fit fails or is noisy, with
## well-measured ones. This table keeps the mean, dispersion and cell-fraction axes (averaged over omega and cbar) and
## is small enough to share.
CELLS <- NULLROWS[NULLROWS$target %in% c("Hc", "He", "H", "Pc"), ]
CELLS$flag_num <- as.numeric(CELLS$flag)
TAB.CELL <- aggregate(cbind(fail, bias, rmse, sd_emp, se_boot, cover, flag_num) ~ target + est + family + mu + phi + frac, data = CELLS,
                      FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
cat("\n== MLE and moment estimator by mean (averaged over the other axes), target Hc and H ==\n")
for (tg in c("Hc", "H")) {
  d <- TAB.CELL[TAB.CELL$target == tg, ]
  by_mu <- aggregate(cbind(fail, bias, rmse, se_boot, flag_num) ~ est + family + mu, data = d, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
  cat("target", tg, "\n"); print(by_mu[order(by_mu$est, by_mu$family, by_mu$mu), ], digits = 2, row.names = FALSE)
}

## ---- 2. decision -------------------------------------------------------------------------------
DEC <- NULLROWS[NULLROWS$target %in% c("Hc", "He", "H") & NULLROWS$family != "gamma" & NULLROWS$usable, ]
SHARE <- c(mle = wmean(DEC$flag[DEC$est == "mle"], DEC$w[DEC$est == "mle"]),
           mom = wmean(DEC$flag[DEC$est == "mom"], DEC$w[DEC$est == "mom"]))
MLE <- DEC[DEC$est == "mle", ]; MOM <- DEC[DEC$est == "mom", ]
BOTH <- merge(MLE[MLE$flag, c("i", "target", "rmse", "w")], MOM[, c("i", "target", "rmse")], by = c("i", "target"), suffixes = c("_mle", "_mom"))
RMSE.REL <- if (nrow(BOTH)) wmean(BOTH$rmse_mom / BOTH$rmse_mle, BOTH$w) else NA_real_
cat(sprintf("\n== Decision ==\nFlagged share of genes (Hc, He, H; lognormal and bimodal): MLE %.3f, moment %.3f\n", SHARE[["mle"]], SHARE[["mom"]]))
cat(sprintf("Moment RMSE relative to MLE in the MLE's flagged cells (weighted): %.2f\n", RMSE.REL))
DECISION <- if (SHARE[["mle"]] < KEEP.BELOW) "KEEP the MLE as is" else if (SHARE[["mle"]] >= ADOPT.ABOVE && SHARE[["mom"]] < SHARE[["mle"]] / 2 &&
              is.finite(RMSE.REL) && RMSE.REL <= RMSE.RATIO) "ADOPT the moment estimator" else "KEEP the MLE and add the correction terms; report the flagged regions"
cat("Decision:", DECISION, "\n")
cat("\nFlagged share of genes by family and target (MLE | moment):\n")
BYFAM <- NULLROWS[NULLROWS$target %in% c("Hc", "He", "H", "Pc") & NULLROWS$usable, ]
BYFAM$tag <- paste(BYFAM$family, BYFAM$target)
print(round(rbind(mle = vapply(split(BYFAM[BYFAM$est == "mle", ], BYFAM$tag[BYFAM$est == "mle"]), function(d) wmean(d$flag, d$w), 0),
                  mom = vapply(split(BYFAM[BYFAM$est == "mom", ], BYFAM$tag[BYFAM$est == "mom"]), function(d) wmean(d$flag, d$w), 0)), 3))

## Flagged regions of the MLE: share of the (omega, cbar, frac, family, target) cells flagged at each mean and dispersion.
REG <- MLE[MLE$family != "gamma", ]
REGION <- tapply(REG$flag, list(mu = REG$mu, phi = REG$phi), mean)
cat("\nShare of MLE cells flagged by mean (rows) and dispersion (columns), non-gamma families:\n"); print(round(REGION, 2))

## ---- 3. agreement of the latent quantities with the derived corrections -------------------------
NULLA <- AGREE[AGREE$scenario == "null", ]
AGREE.TAB <- data.frame(
  quantity = c("phi_Hc: latent Var(a h) vs phi + sigma2_h (1 + phi), log2", "sigma2_h: latent Var(h) / sigma2_C|c over cbar^2",
               "phi_H: latent Var(a_H) vs Appendix A, log2", "delta: simulated minus formula"),
  median = c(median(log2(NULLA$phi_hc / NULLA$pred_hc)), median(NULLA$sigma2_h_hc / NULLA$sigma2_h_hc_formula),
             median(log2(NULLA$phi_H / NULLA$phi_H_formula)), median(NULLA$delta_sim - NULLA$delta_formula, na.rm = TRUE)),
  max_abs = c(max(abs(log2(NULLA$phi_hc / NULLA$pred_hc))), max(abs(NULLA$sigma2_h_hc / NULLA$sigma2_h_hc_formula - 1)),
              max(abs(log2(NULLA$phi_H / NULLA$phi_H_formula))), max(abs(NULLA$delta_sim - NULLA$delta_formula), na.rm = TRUE)),
  mc_se = median(NULLA$truth_mc_se))
cat("\n== Agreement of the latent truth with the derived formulas (null scenario) ==\n"); print(AGREE.TAB, digits = 3, row.names = FALSE)
cat(sprintf("Achieved sigma2_C|c in the simulation / target: median %.3f\n", median(AGREE$var_cc_sim) / VAR.CC))

## ---- 4. contrasts: fitted against truth and against the predicted correction ----------------------
CON <- EST[EST$target %in% c("cis", "trans", "dpar_Pc", "dpar_Pe", "dom") & EST$usable, ]
CSPLIT <- split(CON, list(CON$target, CON$est, CON$scenario), drop = TRUE)
TAB.CON <- data.frame(
  target = vapply(CSPLIT, function(d) d$target[1], ""), est = vapply(CSPLIT, function(d) d$est[1], ""), scenario = vapply(CSPLIT, function(d) d$scenario[1], ""),
  truth = vapply(CSPLIT, function(d) wmean(d$truth, d$w), 0), pred = vapply(CSPLIT, function(d) wmean(d$pred, d$w), 0),
  bias = vapply(CSPLIT, function(d) wmean(d$bias, d$w), 0), bias_vs_pred = vapply(CSPLIT, function(d) wmean(d$bias_pred, d$w), 0),
  se_boot = vapply(CSPLIT, function(d) wmean(d$se_boot, d$w), 0), cover = vapply(CSPLIT, function(d) wmean(d$cover, d$w), 0), row.names = NULL)
cat("\n== Contrasts (log2 phi): fitted - truth, fitted - predicted, bootstrap SE, coverage ==\n")
print(TAB.CON[order(TAB.CON$target, TAB.CON$scenario, TAB.CON$est), ], digits = 3, row.names = FALSE)

if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) {
  write.csv(TAB.ALL, file.path(TABLE.DIR, "estimator_sim_by_family.csv"), row.names = FALSE)
  write.csv(TAB.CELL, file.path(TABLE.DIR, "estimator_sim_by_cell.csv"), row.names = FALSE)
  write.csv(TAB.CON, file.path(TABLE.DIR, "estimator_sim_contrasts.csv"), row.names = FALSE)
  write.csv(AGREE.TAB, file.path(TABLE.DIR, "estimator_sim_agreement.csv"), row.names = FALSE)
  write.csv(EST, file.path(TABLE.DIR, "estimator_sim_all_rows.csv"), row.names = FALSE)
}
