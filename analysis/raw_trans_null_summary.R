###############################################################
### raw_trans_null_summary.R
### Local. Reads the raw_trans_null_output_<k>of<K>.rda files from cluster/scripts/raw_trans_null_grid.R, verifies that every scenario row
### arrived, and scores each null (pool_indep, pool_paired, boot_shift) on the raw and the log2 trans contrast of phi.
###
### Needs   : <INPUT.DIR>/raw_trans_null_inputs.rda (RAW.TRANS.INPUTS) and <OUTPUT.DIR>/raw_trans_null_output_*.rda
### Writes  : tables under TABLE.DIR (when it exists) and a printed verdict per null and scale.
###
### Scoring. Where a scale's null is true (truth = 0) the rejection rate at ALPHA[1] is the type I error. A null HOLDS on that scale when
### no true-null scenario rejects more than ALPHA[1] + 2 binomial SE of the NJ datasets, and is CONSERVATIVE where it rejects less than
### ALPHA[1] - 2 SE (valid, but power is lost). Where the null is false (truth != 0) the rejection rate is the power; the verdict reports its mean over
### those scenarios so a null that holds can be compared on power.
###############################################################

if (!exists("INPUT.DIR"))  INPUT.DIR  <- "."
if (!exists("OUTPUT.DIR")) OUTPUT.DIR <- "."

load(file.path(INPUT.DIR, "raw_trans_null_inputs.rda"))                 # RAW.TRANS.INPUTS
GRID.RT <- RAW.TRANS.INPUTS$GRID; NJ <- RAW.TRANS.INPUTS$NJ; ALPHA <- RAW.TRANS.INPUTS$ALPHA
PART.FILES <- list.files(OUTPUT.DIR, pattern = "^raw_trans_null_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(PART.FILES) > 0)
ROWS <- integer(0); RES <- list()
for (f in PART.FILES) { load(f); ROWS <- c(ROWS, ROWS.PART); RES <- c(RES, RESULTS.PART) }
stopifnot("every grid row must arrive exactly once" = identical(sort(ROWS), GRID.RT$i))
R <- do.call(rbind, RES)
R <- cbind(scenario = GRID.RT$scenario[R$i], n_cells = GRID.RT$n_cells[R$i], R)
R$null_true <- abs(R$truth) < 1e-12
SE <- sqrt(ALPHA[1] * (1 - ALPHA[1]) / NJ)
R$verdict <- ifelse(!R$null_true, "power", ifelse(R$rate_a1 > ALPHA[1] + 2 * SE, "ANTI-CONSERVATIVE", ifelse(R$rate_a1 < ALPHA[1] - 2 * SE, "conservative", "holds")))
R$null_to_obs_sd <- R$null_sd / R$sd_obs

cat("== Rejection rate at alpha =", ALPHA[1], "by scenario, null and scale (truth = true trans contrast on that scale) ==\n")
print(R[order(R$scale, R$scenario, R$n_cells, R$method), c("scale", "scenario", "n_cells", "method", "truth", "n_valid", "rate_a1", "rate_a2", "null_to_obs_sd", "verdict")],
      digits = 3, row.names = FALSE)

cat("\n== Verdict per null and scale: true-null scenarios (type I error) and false-null scenarios (power) ==\n")
VERDICT <- do.call(rbind, lapply(split(R, list(R$method, R$scale), drop = TRUE), function(d) {
  tn <- d[d$null_true, ]; al <- d[!d$null_true, ]
  data.frame(method = d$method[1], scale = d$scale[1],
             max_type1 = max(tn$rate_a1), n_anticonservative = sum(tn$verdict == "ANTI-CONSERVATIVE"), n_conservative = sum(tn$verdict == "conservative"),
             worst_scenario = tn$scenario[which.max(tn$rate_a1)], mean_power = mean(al$rate_a1), row.names = NULL)
}))
print(VERDICT[order(VERDICT$scale, VERDICT$method), ], digits = 3, row.names = FALSE)
cat("\nA null is adopted for a scale when n_anticonservative is 0; among those, the larger mean_power wins. The pooled nulls also test that\n",
    "the groups share one distribution, so a rejection in a scenario where the scale's own null holds shows what they add beyond it.\n", sep = "")

if (exists("TABLE.DIR") && dir.exists(TABLE.DIR)) {
  write.csv(R, file.path(TABLE.DIR, "raw_trans_null_all_rows.csv"), row.names = FALSE)
  write.csv(VERDICT, file.path(TABLE.DIR, "raw_trans_null_verdict.csv"), row.names = FALSE)
}
