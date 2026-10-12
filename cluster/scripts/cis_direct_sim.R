###############################################################
### cis_direct_sim.R
### Compares three estimators of the cis contrast of dispersion on the same simulated paired hybrid cells:
### option 2 (an NB maximum-likelihood fit per allele, then log2 phi_Hc - log2 phi_He, with a percentile
### bootstrap), a direct joint-likelihood estimate with a likelihood-ratio interval (cis_profile), and a
### direct paired moment estimate of phi_Hc - phi_He (cis_paired_moment). The question is how much of the
### cis contrast that option 2 loses to failed fits the direct estimators recover, and whether pairing the
### alleles tightens the contrast. A separate study; it reads no existing result.
###
### Inputs  : estimator_sim_inputs.rda, written by analysis/sim_inputs.R (the same file estimator_sim.R reads)
### Output  : cis_direct_output_<k>of<K>.rda for array task k of K, holding CIS.DIRECT (one row per grid row)
###           and ROWS.PART (the GRID rows this task computed). analysis/cis_direct_summary.R reassembles them.
###
### The rows are the null and cis scenarios at c-bar 0.5, intrinsic fractions 0, 0.5 and 1, and cell
### fractions 1 and 0.25. As in estimator_sim.R, phase 1 draws the latent truth of each truth key the task owns
### (keys are dealt out round-robin, so no truth is drawn twice) and phase 2 simulates the rows of those keys.
### Seeds come from SEED.BASE and the key or row index, so results do not depend on how keys are split across tasks.
###############################################################

library('parallel')

source("functions.R")

load(Sys.getenv("SIM_INPUTS", unset = "estimator_sim_inputs.rda"))
## GRID, LIB.HYB, LIB.PAR, DEPTH.REF, N.OBS, VAR.C, VAR.CC, NREP, NBOOT, NTRUE, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE <- sprintf("cis_direct_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the whole array only
## computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

cl <- makeCluster(NUM.CORES, type = "FORK")

## ---- the keys and rows of this task ----
SEL       <- which(GRID$scenario %in% c("null", "cis") & GRID$cbar == 0.5 & GRID$frac %in% c(1, 0.25) & GRID$omega %in% c(0, 0.5, 1))
ALL.KEYS  <- sort(unique(GRID$tkey[SEL]))
MY.KEYS   <- ALL.KEYS[(seq_along(ALL.KEYS) - 1L) %% N.ARRAY == ARRAY.ID - 1L]
MY.ROWS   <- SEL[GRID$tkey[SEL] %in% MY.KEYS]
N.POINTS  <- length(MY.ROWS)
N.UPDATES <- 40

cat(sprintf("cis direct simulation start: array task %d of %d, %d truth keys, %d rows, NREP=%d, NBOOT=%d, NTRUE=%.0e, %d cores\n",
            ARRAY.ID, N.ARRAY, length(MY.KEYS), N.POINTS, NREP, NBOOT, NTRUE, NUM.CORES))
flush.console()

## ---- phase 1: latent truth of each key ----
## estimator_sim_truth: the population variances of the latent factors of one truth key (here Var(a h) of the two
## alleles is what the cis contrast targets), drawn from the latent model and not from the Appendix A formulas.
t0 <- Sys.time()
TRUTH <- parLapply(cl, MY.KEYS, estimator_sim_truth, GRID = GRID, LIB.HYB = LIB.HYB, LIB.PAR = LIB.PAR, DEPTH.REF = DEPTH.REF,
                   NTRUE = NTRUE, SEED.BASE = SEED.BASE, VAR.C = VAR.C, VAR.CC = VAR.CC)
names(TRUTH) <- MY.KEYS
cat(sprintf("[%s] truth done: %d keys in %.1f min\n", format(Sys.time(), "%H:%M:%S"), length(MY.KEYS),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
flush.console()

## ---- phase 2: the rows, chunked so the master can report progress ----
CHUNK.SIZE <- max(NUM.CORES, ceiling(N.POINTS / N.UPDATES))
chunks     <- split(MY.ROWS, ceiling(seq_along(MY.ROWS) / CHUNK.SIZE))
results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## cis_direct_row: one GRID row. It simulates NREP hybrid datasets, applies the three cis estimators to the same
  ## cells (option 2 with a cell bootstrap, the joint profile likelihood, the paired moment difference) and scores
  ## each against the latent cis truth: failure and bias, SE, coverage and interval width, and for the datasets
  ## where option 2 fails, how often the profile interval is one- or two-sided and excludes zero.
  results[[k]] <- parLapply(cl, chunks[[k]], cis_direct_row, GRID = GRID, TRUTH = TRUTH, LIB.HYB = LIB.HYB, DEPTH.REF = DEPTH.REF,
                            N.OBS = N.OBS, NREP = NREP, NBOOT = NBOOT, SEED.BASE = SEED.BASE, VAR.C = VAR.C)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (N.POINTS - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS, 100 * done / N.POINTS, el, eta))
  flush.console()
}

CIS.DIRECT <- do.call(rbind, do.call(c, results))
ROWS.PART  <- MY.ROWS
save(CIS.DIRECT, ROWS.PART, file = OUT.FILE)

cat(sprintf("cis direct simulation done: %d rows in %.1f min\n",
            length(ROWS.PART), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
