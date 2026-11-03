###############################################################
### pilot_dispersion_quantiles.R
### Local. Pulls the dispersion estimates (phi, log10 scale) of the four pilot datasets Pc (MIX.SC), Pe (MIX.SE), Hc (HYB.SC)
### and He (HYB.SE) and reports their quantiles. The pilot job writes only bootstrap SEs, so the point estimates come from
### refitting the pilot inputs (PILOT.MATS, PILOT.EXPOS, PILOT.GENES) with the pipeline's own estimator .fit_one().
###
### Run from the project root, either after Section 2 (the PILOT.* objects are then in the session) or with
###   INPUT.DIR <- "folder holding gene_pilot_inputs.rda"
### set first. Needs R/functions.R (sourced when .fit_one() is not loaded).
### Writes <TABLE.DIR>/pilot_dispersion_quantiles.csv when TABLE.DIR exists.
###############################################################

if (!exists(".fit_one")) source("R/functions.R")
if (!exists("PILOT.MATS")) {
  if (!exists("INPUT.DIR")) INPUT.DIR <- "."
  load(file.path(INPUT.DIR, "gene_pilot_inputs.rda"))                   # PILOT.MATS, PILOT.EXPOS, PILOT.GENES, N.BOOT, SEED.BOOT
}
SETS <- c(Pc = "MIX.SC", Pe = "MIX.SE", Hc = "HYB.SC", He = "HYB.SE")   # label -> pilot dataset
QUANTILE.PROBS <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)

## .fit_one() returns the dispersion phi in the current code (disp = phi, theta = 1 / phi). The earlier code returned the NB
## size theta in disp, so phi is 1 / disp there. fit_dispersion() exists only in the current code.
PHI.IS.DISP <- exists("fit_dispersion")
cat(if (PHI.IS.DISP) "functions.R returns phi in disp\n" else "functions.R is the earlier version: disp is theta, phi = 1 / disp\n")

## One phi per pilot gene and dataset. A gene with no counts, a failed fit or phi = 0 (Poisson limit) has no log10 phi.
## Each dataset needs an offset per cell. Pilot inputs written by the earlier Section 2 hold one shared hybrid offset under
## the name HYB instead of one per allele (HYB.SC, HYB.SE); that offset is used for both alleles here and the script says so.
cat("PILOT.EXPOS holds:", paste(names(PILOT.EXPOS), collapse = ", "), "\n")
PHI <- lapply(SETS, function(ds) {
  m <- PILOT.MATS[[ds]][PILOT.GENES, , drop = FALSE]; e <- PILOT.EXPOS[[ds]]
  if (is.null(e) && !is.null(PILOT.EXPOS[["HYB"]]) && ds %in% c("HYB.SC", "HYB.SE")) {
    cat(sprintf("%s: no per-allele offset in PILOT.EXPOS, using the shared hybrid offset HYB\n", ds)); e <- PILOT.EXPOS[["HYB"]]
  }
  stopifnot("no offset for this dataset in PILOT.EXPOS" = !is.null(e), "offset length differs from the number of cells" = length(e) == ncol(m))
  d <- vapply(seq_len(nrow(m)), function(i) .fit_one(m[i, ], e)[["disp"]], numeric(1))
  setNames(if (PHI.IS.DISP) d else 1 / d, rownames(m))
})

## Quantiles of log10 phi over the genes with a positive finite phi, one row per dataset.
LOG10.PHI <- lapply(PHI, function(p) log10(p[is.finite(p) & p > 0]))
PILOT.DISP.QUANTILES <- do.call(rbind, lapply(names(SETS), function(k) {
  x <- LOG10.PHI[[k]]
  data.frame(set = k, dataset = SETS[[k]], genes = length(PHI[[k]]), n_phi = length(x), n_zero_or_failed = length(PHI[[k]]) - length(x),
             t(setNames(quantile(x, QUANTILE.PROBS, names = FALSE), paste0("q", QUANTILE.PROBS * 100))), row.names = NULL)
}))
cat("\nQuantiles of log10(phi), pilot genes with a positive finite phi (n_zero_or_failed genes have none):\n")
print(PILOT.DISP.QUANTILES, digits = 3, row.names = FALSE)
cat("\nSame quantiles on the phi scale (10^q):\n")
print(cbind(PILOT.DISP.QUANTILES[, c("set", "n_phi")], round(10^PILOT.DISP.QUANTILES[, grep("^q", names(PILOT.DISP.QUANTILES))], 4)), row.names = FALSE)

if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) write.csv(PILOT.DISP.QUANTILES, file.path(TABLE.DIR, "pilot_dispersion_quantiles.csv"), row.names = FALSE)
