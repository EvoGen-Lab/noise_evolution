###############################################################
### cis_option_decision.R
### Local. Turns the finished cis_direct_sim.R output into a recommendation between option 2 (an NB fit per allele,
### then the log2 ratio, with a percentile bootstrap) and option 4 (direct estimation of the cis contrast by the joint
### profile likelihood with a likelihood-ratio interval). It first runs cis_direct_summary.R, which loads and verifies
### the cluster output and builds the per-row table D, and then scores five criteria whose thresholds are set below.
###
### Steps
###   1. Check on the cluster that every array task wrote its file: cis_direct_output_<k>of<K>.rda for k = 1..K.
###   2. Copy the files and estimator_sim_inputs.rda into one local folder (or two: inputs and outputs).
###   3. Open R in the project root and run
###        INPUT.DIR  <- "path/to/folder with estimator_sim_inputs.rda"
###        OUTPUT.DIR <- "path/to/folder with the cis_direct_output files"
###        source("analysis/cis_option_decision.R")
###      The script needs only R/functions.R (it sources it when wmean() is not loaded). It does not need R/setup.R;
###      when setup.R has been run, its INPUT.DIR, OUTPUT.DIR and TABLE.DIR are used as they are.
###   4. Read the verdict block at the end. Every criterion is printed with its value, its threshold and pass or fail,
###      and the table is written to <TABLE.DIR>/cis_option_decision.csv.
###
### Scope: the simulation covers the cis contrast only. A direct trans estimator was not simulated, so a verdict for
### option 4 applies to cis and trans stays with option 2.
###############################################################

if (!exists("INPUT.DIR"))  INPUT.DIR  <- "."
if (!exists("OUTPUT.DIR")) OUTPUT.DIR <- "."
if (!exists("wmean")) source("R/functions.R")
source("analysis/cis_direct_summary.R")   # builds D, TAB.REC

## ---- thresholds (a criterion passes when its value meets the threshold) ----------------------------------------
LOSS.MIN      <- 0.05   # option 2 gives no cis contrast in at least this share of datasets (gene-weighted): direct estimation adds information
GAIN.MIN      <- 0.10   # where option 2 fails, the share of cis-scenario datasets whose profile interval excludes zero exceeds the null-scenario share by at least this much
COVER.TARGET  <- 0.95   # nominal coverage of the interval
COVER.TOL     <- 0.02   # calibrated: coverage no lower than the target minus this (over-coverage is conservative and passes)
WIDTH.MAX     <- 1.00   # profile width / option 2 width on the datasets option 2 can fit (at most this)
BIAS.MAX      <- 0.10   # profile bias of the maximiser at most this many log2 units, and not worse than option 2 by more than BIAS.WORSE
BIAS.WORSE    <- 0.05

## ---- criteria -----------------------------------------------------------------------------------------------
LOSS      <- wmean(D$miss, D$w)
## Where option 2 fails the profile interval is one-sided (one allele has phi-hat = 0, so one end is unbounded), so the
## measure of recovery is how often it excludes zero in cis scenarios (power) beyond how often it does in null scenarios (type I).
EX0.CIS   <- with(D[D$scenario == "cis", ],  sum(w * rec_excl0 * miss, na.rm = TRUE) / sum(w * miss))
EX0.NULL  <- with(D[D$scenario == "null", ], sum(w * rec_excl0 * miss, na.rm = TRUE) / sum(w * miss))
RECOVER   <- EX0.CIS - EX0.NULL
COVER.P   <- wmean(D$p_cover, D$w)            # all datasets, including those option 2 drops
COVER.E2  <- wmean(D$e2_cover, D$w)           # datasets option 2 can fit
COVER.PE  <- wmean(D$p_cover_e2ok, D$w)       # profile on the same datasets
WIDTH     <- wmean(D$p_width_e2ok, D$w) / wmean(D$e2_width, D$w)
BIAS.P    <- abs(wmean(D$p_bias_e2ok, D$w)); BIAS.E2 <- abs(wmean(D$e2_bias, D$w))

DECISION <- data.frame(
  criterion = c("C1 option 2 loses the cis contrast in", "C2 excludes-zero rate beyond null where option 2 fails",
                "C3 profile coverage, all datasets", "C3b profile coverage, datasets option 2 can fit", "   option 2 coverage (reference)",
                "C4 profile width / option 2 width", "C5 profile bias (log2)", "   option 2 bias (log2)"),
  value = c(LOSS, RECOVER, COVER.P, COVER.PE, COVER.E2, WIDTH, BIAS.P, BIAS.E2),
  threshold = c(LOSS.MIN, GAIN.MIN, COVER.TARGET, COVER.TARGET, COVER.TARGET, WIDTH.MAX, BIAS.MAX, NA),
  rule = c(">=", ">=", paste("above -", COVER.TOL), paste("above -", COVER.TOL), "reference", "<=", "<=", "reference"),
  stringsAsFactors = FALSE)
DECISION$pass <- c(LOSS >= LOSS.MIN, RECOVER >= GAIN.MIN, COVER.P >= COVER.TARGET - COVER.TOL, COVER.PE >= COVER.TARGET - COVER.TOL,
                   NA, WIDTH <= WIDTH.MAX, BIAS.P <= BIAS.MAX && BIAS.P - BIAS.E2 <= BIAS.WORSE, NA)
cat("\n== Option 2 versus option 4 (direct profile likelihood) for the cis contrast ==\n")
print(DECISION, digits = 3, row.names = FALSE)
cat(sprintf("Where option 2 fails, the profile interval excludes zero in %.3f of cis-scenario and %.3f of null-scenario datasets (type I).\n", EX0.CIS, EX0.NULL))

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
