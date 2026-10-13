###############################################################
### cis_option_decision.R
### Local. Turns the finished cis_direct_sim.R output into a recommendation between option 2 (an NB fit per allele,
### then the log2 ratio, with a percentile bootstrap) and option 4 (direct estimation of the cis contrast by the joint
### profile likelihood with a likelihood-ratio interval). It first runs cis_direct_summary.R, which loads and verifies
### the cluster output and builds the per-row table D, and then scores five criteria whose thresholds are set below.
###
### Steps
###   1. Check on the cluster that every array task wrote its file: cis_direct_output_<k>of<K>.rda for k = 1..K.
###   2. Copy the files and estimator_sim_inputs.rda to the local OUTPUT.DIR and INPUT.DIR (R/setup.R defines both).
###   3. Run this script from the project root:
###        library(here); source(here("R", "setup.R")); source(here("analysis", "cis_option_decision.R"))
###      Without setup.R, set INPUT.DIR and OUTPUT.DIR first and source R/functions.R.
###   4. Read the verdict block at the end. Every criterion is printed with its value, its threshold and pass or fail,
###      and the table is written to <TABLE.DIR>/cis_option_decision.csv.
###
### Scope: the simulation covers the cis contrast only. A direct trans estimator was not simulated, so a verdict for
### option 4 applies to cis and trans stays with option 2.
###############################################################

if (!exists("INPUT.DIR"))  INPUT.DIR  <- "."
if (!exists("OUTPUT.DIR")) OUTPUT.DIR <- "."
if (!exists("wmean")) source("R/functions.R")
source(if (exists("here")) here::here("analysis", "cis_direct_summary.R") else "analysis/cis_direct_summary.R")   # builds D, TAB.REC

## ---- thresholds (a criterion passes when its value meets the threshold) ----------------------------------------
LOSS.MIN      <- 0.05   # option 2 gives no cis contrast in at least this share of datasets (gene-weighted): direct estimation adds information
RECOVER.MIN   <- 0.50   # of those datasets, at least this share gets a two-sided (bounded) profile interval
COVER.TARGET  <- 0.95   # nominal coverage of the interval
COVER.TOL     <- 0.02   # the profile interval is calibrated when its coverage lies within this of the target
WIDTH.MAX     <- 1.00   # profile width / option 2 width on the datasets option 2 can fit (at most this)
BIAS.MAX      <- 0.10   # profile bias of the maximiser at most this many log2 units, and not worse than option 2 by more than BIAS.WORSE
BIAS.WORSE    <- 0.05

## ---- criteria -----------------------------------------------------------------------------------------------
LOSS      <- wmean(D$miss, D$w)
RECOVER   <- sum(D$w * D$rec_two_sided * D$miss, na.rm = TRUE) / sum(D$w * D$miss)   # a bounded interval: the contrast is located, not only signed
RECOVER.X0 <- sum(D$w * D$rec_excl0 * D$miss, na.rm = TRUE) / sum(D$w * D$miss)         # the interval excludes zero: a cis call is supported
COVER.P   <- wmean(D$p_cover, D$w)            # all datasets, including those option 2 drops
COVER.E2  <- wmean(D$e2_cover, D$w)           # datasets option 2 can fit
COVER.PE  <- wmean(D$p_cover_e2ok, D$w)       # profile on the same datasets
WIDTH     <- wmean(D$p_width_e2ok, D$w) / wmean(D$e2_width, D$w)
BIAS.P    <- abs(wmean(D$p_bias_e2ok, D$w)); BIAS.E2 <- abs(wmean(D$e2_bias, D$w))

DECISION <- data.frame(
  criterion = c("C1 option 2 loses the cis contrast in", "C2 profile gives a bounded interval where option 2 fails",
                "C3 profile coverage, all datasets", "C3b profile coverage, datasets option 2 can fit", "   option 2 coverage (reference)",
                "C4 profile width / option 2 width", "C5 profile bias (log2)", "   option 2 bias (log2)"),
  value = c(LOSS, RECOVER, COVER.P, COVER.PE, COVER.E2, WIDTH, BIAS.P, BIAS.E2),
  threshold = c(LOSS.MIN, RECOVER.MIN, COVER.TARGET, COVER.TARGET, COVER.TARGET, WIDTH.MAX, BIAS.MAX, NA),
  rule = c(">=", ">=", paste("within", COVER.TOL), paste("within", COVER.TOL), "reference", "<=", "<=", "reference"),
  stringsAsFactors = FALSE)
DECISION$pass <- c(LOSS >= LOSS.MIN, RECOVER >= RECOVER.MIN, abs(COVER.P - COVER.TARGET) <= COVER.TOL, abs(COVER.PE - COVER.TARGET) <= COVER.TOL,
                   NA, WIDTH <= WIDTH.MAX, BIAS.P <= BIAS.MAX && BIAS.P - BIAS.E2 <= BIAS.WORSE, NA)
cat("\n== Option 2 versus option 4 (direct profile likelihood) for the cis contrast ==\n")
print(DECISION, digits = 3, row.names = FALSE)
cat(sprintf("Where option 2 fails, the profile interval excludes zero in %.3f of datasets (a cis call is supported there).\n", RECOVER.X0))

## ---- where option 2 is the weak one -------------------------------------------------------------------------
## Regimes (mean x dispersion) where option 2 is missing a large share of contrasts or its interval under-covers,
## listed with the profile interval's coverage in the same cells.
REGIME <- aggregate(cbind(miss = D$miss, e2_cover = D$e2_cover, p_cover = D$p_cover) ~ mu + phi, data = D, FUN = mean, na.action = na.pass)
REGIME <- REGIME[order(-REGIME$miss), ]
cat("\nCells (mean, dispersion) ranked by the share of datasets where option 2 gives no cis contrast (unweighted means over rows):\n")
print(head(REGIME, 10), digits = 2, row.names = FALSE)

## ---- verdict ---------------------------------------------------------------------------------------------------
CALIBRATED <- isTRUE(DECISION$pass[3]) && isTRUE(DECISION$pass[4])
GAIN       <- isTRUE(DECISION$pass[1]) && isTRUE(DECISION$pass[2])
SHARPER    <- isTRUE(DECISION$pass[6]) && isTRUE(DECISION$pass[7])
VERDICT <- if (CALIBRATED && (GAIN || SHARPER)) {
  "ADOPT option 4 for the cis contrast: the profile interval is calibrated and either recovers contrasts option 2 loses (C1, C2) or is at least as sharp and unbiased on the datasets option 2 can fit (C4, C5)."
} else if (!CALIBRATED && GAIN) {
  "MIXED: option 2 loses many cis contrasts but the profile interval is not calibrated here; keep option 2 and use the profile estimate only to flag genes option 2 cannot call, until its calibration is understood."
} else {
  "KEEP option 2: the direct estimator does not recover enough lost contrasts or is not better where option 2 works."
}
cat("\nVERDICT:", VERDICT, "\n")
cat("Trans: no direct estimator was simulated, so trans stays with option 2 whatever the cis verdict.\n")

if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) {
  write.csv(DECISION, file.path(TABLE.DIR, "cis_option_decision.csv"), row.names = FALSE)
  write.csv(REGIME,   file.path(TABLE.DIR, "cis_option_regimes.csv"),  row.names = FALSE)
}
