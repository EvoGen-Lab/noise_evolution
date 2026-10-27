###############################################################
### cis_direct_summary.R
### Local. Reads the cis_direct_output_<k>of<K>.rda files from cluster/scripts/cis_direct_sim.R, verifies that every
### selected grid row arrived, and compares two estimators of the cis contrast of dispersion: option 2 (an NB
### fit per allele, then the log2 ratio) and the joint profile likelihood.
###
### Needs   : R/functions.R (sourced here if wmean() is not loaded), <INPUT.DIR>/estimator_sim_inputs.rda (GRID,
###           its weights w) and <OUTPUT.DIR>/cis_direct_output_*.rda
### Reports : (1) how much of the cis contrast option 2 loses, and what the profile interval
###           recovers there; (2) coverage and width where option 2 works.
###############################################################

if (!exists("INPUT.DIR"))  INPUT.DIR  <- "."
if (!exists("OUTPUT.DIR")) OUTPUT.DIR <- "."
if (!exists("wmean")) source("R/functions.R")   # wmean() lives in functions.R

## ---- load and verify -------------------------------------------------------------------------
load(file.path(INPUT.DIR, "estimator_sim_inputs.rda"))                # GRID
PART.FILES <- list.files(OUTPUT.DIR, pattern = "^cis_direct_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(PART.FILES) > 0)
PARTS <- vector("list", length(PART.FILES)); ROWS.SEEN <- integer(0)
for (k in seq_along(PART.FILES)) {
  e <- new.env(); load(PART.FILES[k], envir = e)
  PARTS[[k]] <- e$CIS.DIRECT; ROWS.SEEN <- c(ROWS.SEEN, e$ROWS.PART)
}
SEL <- which(GRID$scenario %in% c("null", "cis") & GRID$cbar == 0.5 & GRID$frac %in% c(1, 0.25) & GRID$omega %in% c(0, 0.5, 1))
stopifnot(setequal(ROWS.SEEN, SEL), !anyDuplicated(ROWS.SEEN))        # every selected row, once
CIS <- do.call(rbind, PARTS)
stopifnot(setequal(CIS$i, SEL))
D <- merge(CIS, GRID[, c("i", "mu", "phi", "omega", "frac", "family", "scenario", "w")], by = "i")
D$miss <- D$e2_fail                                                    # share of datasets where option 2 gives no contrast

## ---- 1. recovery where option 2 fails ----------------------------------------------------------
## Weighted over the grid by the share of genes at each mean and dispersion (wmean with D$w; shares among failures are sum(w x share x missing) / sum(w x missing)). rec_* are shares among the datasets where
## option 2 fails; an interval "excludes zero" when it supports a cis difference, "one-sided" when only one end is bounded.
cat("== 1. What option 2 loses and what the direct estimators recover (gene-weighted) ==\n")
SPLIT <- split(D, list(D$scenario, D$family), drop = TRUE)
TAB.REC <- data.frame(
  scenario = vapply(SPLIT, function(d) d$scenario[1], ""), family = vapply(SPLIT, function(d) d$family[1], ""),
  e2_missing = vapply(SPLIT, function(d) wmean(d$miss, d$w), 0),
  rec_two_sided = vapply(SPLIT, function(d) sum(d$w * d$rec_two_sided * d$miss, na.rm = TRUE) / sum(d$w * d$miss), 0),
  rec_one_sided = vapply(SPLIT, function(d) sum(d$w * d$rec_one_sided * d$miss, na.rm = TRUE) / sum(d$w * d$miss), 0),
  rec_excl0 = vapply(SPLIT, function(d) sum(d$w * d$rec_excl0 * d$miss, na.rm = TRUE) / sum(d$w * d$miss), 0),
  rec_cover = vapply(SPLIT, function(d) sum(d$w * d$rec_cover * d$miss, na.rm = TRUE) / sum(d$w * d$miss), 0), row.names = NULL)
print(TAB.REC, digits = 3, row.names = FALSE)
cat(sprintf("\nAll rows: option 2 gives no contrast for %.3f of datasets; among those the profile interval is one-sided in %.3f,\ntwo-sided in %.3f, excludes zero in %.3f (covers the truth in %.3f).\n",
            wmean(D$miss, D$w), sum(D$w * D$rec_one_sided * D$miss, na.rm = TRUE) / sum(D$w * D$miss), sum(D$w * D$rec_two_sided * D$miss, na.rm = TRUE) / sum(D$w * D$miss), sum(D$w * D$rec_excl0 * D$miss, na.rm = TRUE) / sum(D$w * D$miss),
            sum(D$w * D$rec_cover * D$miss, na.rm = TRUE) / sum(D$w * D$miss)))
cat("Where it is recovered, by mean (rows) and dispersion (columns): share of option-2-failed datasets whose profile interval excludes 0\n")
print(round(tapply(D$rec_excl0, list(mu = D$mu, phi = D$phi), mean, na.rm = TRUE), 2))
cat("Share of datasets where option 2 fails, same layout:\n")
print(round(tapply(D$miss, list(mu = D$mu, phi = D$phi), mean, na.rm = TRUE), 2))

## ---- 2. coverage, bias and width where option 2 works -----------------------------------------
## The profile columns marked e2ok use the same datasets as option 2, so they compare like with like.
cat("\n== 2. Where option 2 works: bias, SE, coverage and interval width (gene-weighted, log2 units) ==\n")
CMP <- data.frame(
  quantity = c("option 2: bias", "option 2: empirical SD", "option 2: bootstrap SE", "option 2: coverage of the percentile interval", "option 2: median interval width",
               "profile: bias of the maximiser", "profile: coverage of the LR interval (same datasets)", "profile: median interval width (same datasets, two-sided)"),
  value = c(wmean(D$e2_bias, D$w), wmean(D$e2_sd, D$w), wmean(D$e2_se_boot, D$w), wmean(D$e2_cover, D$w), wmean(D$e2_width, D$w), wmean(D$p_bias_e2ok, D$w), wmean(D$p_cover_e2ok, D$w), wmean(D$p_width_e2ok, D$w)))
print(CMP, digits = 3, row.names = FALSE)
cat(sprintf("Profile interval width / option 2 width: %.2f\n", wmean(D$p_width_e2ok, D$w) / wmean(D$e2_width, D$w)))
cat("Coverage over all datasets, including those option 2 drops: profile", round(wmean(D$p_cover, D$w), 3), "\n")

if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) {
  write.csv(TAB.REC, file.path(TABLE.DIR, "cis_direct_recovery.csv"), row.names = FALSE)
  write.csv(D, file.path(TABLE.DIR, "cis_direct_all_rows.csv"), row.names = FALSE)
}
