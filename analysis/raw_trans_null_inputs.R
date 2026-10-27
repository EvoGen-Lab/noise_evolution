###############################################################
### raw_trans_null_inputs.R
### Local. Builds the scenario grid of the trans-null validation and writes raw_trans_null_inputs.rda for
### cluster/scripts/raw_trans_null_grid.R. The study asks whether the pipeline's trans nulls (pooled relabeling with independent or
### paired allele shuffles) and a null-constrained resampling (bootstrap recentered on the observed value) hold their type I error for
### the trans contrast of phi on the raw (difference) scale and on the log2 scale, and what each costs in power.
###
### A scenario fixes the true phi of the four groups, so the raw trans null (phi_pc - phi_pe = phi_hc - phi_he) and the log2 trans
### null (equal ratios) are each known to hold or fail. Scenarios cover: all groups equal; a cis-type difference that parents and
### hybrid share (both nulls hold); a hybrid phi scale 1.5x above the parents (the log2 null holds, the raw null fails); equal
### parental and hybrid differences with unequal ratios (the raw null holds, the log2 null fails); a trans-only effect (both fail);
### hybrid mean reads 4x the parents' (pooling cells with different means); a low-mean, low-phi gene near the Poisson boundary; and
### unequal exposure spread between parents and hybrid cells.
###
### Writes  : <INPUT.DIR>/raw_trans_null_inputs.rda (RAW.TRANS.INPUTS). Submit cluster/submit/raw_trans_null.sub afterwards.
###############################################################

library(here); source(here("R", "setup.R"))

SCEN <- data.frame(
  scenario = c("all_equal", "shared_cis_difference", "hybrid_scale_1.5x", "equal_difference_unequal_ratio", "trans_only",
               "hybrid_mean_4x", "low_mean_low_phi", "exposure_cv_differs"),
  phi_pc = c(0.4, 0.4, 0.4, 0.5, 0.4, 0.4, 0.05, 0.4), phi_pe = c(0.4, 0.2, 0.2, 0.3, 0.2, 0.4, 0.05, 0.4),
  phi_hc = c(0.4, 0.4, 0.6, 0.8, 0.3, 0.4, 0.05, 0.4), phi_he = c(0.4, 0.2, 0.3, 0.6, 0.3, 0.4, 0.05, 0.4),
  mu_p   = c(2, 2, 2, 2, 2, 1, 0.3, 2),               mu_h   = c(2, 2, 2, 2, 2, 4, 0.3, 2),
  cv_p   = 0.8,                                       cv_h   = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.3),
  stringsAsFactors = FALSE)
## Each scenario at three cell numbers, so the effect of group size on calibration is visible.
GRID.RT <- do.call(rbind, lapply(c(100, 200, 400), function(n) cbind(SCEN, n_cells = n)))
GRID.RT$i <- seq_len(nrow(GRID.RT))

RAW.TRANS.INPUTS <- list(GRID = GRID.RT, SHARED = 0.5, NJ = 1000, NI = 1000, NBOOT = 1000, ALPHA = c(0.05, 0.01), SEED.BASE = SEED.SECTION + 1500)
save(RAW.TRANS.INPUTS, file = file.path(INPUT.DIR, "raw_trans_null_inputs.rda"))
cat(sprintf("raw trans null inputs written: %d rows (%d scenarios x 3 cell numbers), NJ = %d, NI = %d, NBOOT = %d\n",
            nrow(GRID.RT), nrow(SCEN), RAW.TRANS.INPUTS$NJ, RAW.TRANS.INPUTS$NI, RAW.TRANS.INPUTS$NBOOT))
cat("Submit raw_trans_null.sub (job array) now; run raw_trans_null_summary.R once every raw_trans_null_output_<k>of<K>.rda exists.\n")
