###############################################################
### estimator_sim.R
### Simulation study of the dispersion estimators: does the NB maximum-likelihood fit
### (.fit_one()) recover the dispersion of a count-scale factor that is not gamma, and how does
### a Poisson-corrected moment estimator compare on the same simulated cells? Runs one
### estimator_sim_row() per GRID row on a FORK cluster and saves the per-row summaries. Mirrors
### the structure of power_grid.R. A separate study: no existing result depends on it.
###
### Inputs  : estimator_sim_inputs.rda, written by analysis/sim_inputs.R
###           (GRID, LIB.HYB, LIB.PAR, DEPTH.REF, N.OBS, VAR.C, EPS.POOL, VAR.EPS, NREP, NBOOT, NTRUE, SEED.BASE, ...)
### Output  : estimator_sim_output_<k>of<K>.rda for array task k of K, holding SIM.EST (one row per
###           grid row x target x estimator), SIM.AGREE (latent quantities beside the formula
###           predictions, one row per grid row) and ROWS.PART (the GRID rows this task computed).
###           analysis/estimator_sim_summary.R reassembles them.
###
### Two phases per task. Phase 1 draws the latent truth of each truth key the task owns (rows that
### share family, dispersion, intrinsic fraction, mean share and scenario share one truth). Keys are dealt
### out round-robin to the array tasks, and a task runs every grid row of its keys, so no truth is drawn
### twice. Phase 2 simulates the replicates of each row against its truth. Each truth seeds itself from
### SEED.BASE and its key, each row from SEED.BASE and its index, so results are identical however the
### keys are split across tasks.
### NOTE: the simulation models the earlier Section 2 offsets (scran haplotype size factors with an estimation error), not the
### library-size offsets and allele-share correction the pipeline now uses; it has not been updated.
###############################################################

library('parallel')

source("functions.R")

load(Sys.getenv("SIM_INPUTS", unset = "estimator_sim_inputs.rda"))
## GRID, LIB.HYB, LIB.PAR, DEPTH.REF, N.OBS, VAR.C, EPS.POOL, VAR.EPS, NREP, NBOOT, NTRUE, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE <- sprintf("estimator_sim_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the whole array only
## computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

cl <- makeCluster(NUM.CORES, type = "FORK")

## ---- the keys and rows of this task ----
N.UPDATES <- 40
MY.KEYS   <- which((seq_len(max(GRID$tkey)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
MY.ROWS   <- which(GRID$tkey %in% MY.KEYS)
N.POINTS  <- length(MY.ROWS)

cat(sprintf("estimator simulation start: array task %d of %d, %d truth keys, %d rows, NREP=%d, NBOOT=%d, NTRUE=%.0e, %d cores\n",
            ARRAY.ID, N.ARRAY, length(MY.KEYS), N.POINTS, NREP, NBOOT, NTRUE, NUM.CORES))
cat(sprintf("inputs: sigma2_c=%.3g; offset error variances Pc=%.3g Pe=%.3g Hc=%.3g He=%.3g; depth median %.0f; cells HYB=%d, MIX.SC=%d, MIX.SE=%d\n",
            VAR.C, VAR.EPS[["PC"]], VAR.EPS[["PE"]], VAR.EPS[["HC"]], VAR.EPS[["HE"]], DEPTH.REF, N.OBS[["HYB"]], N.OBS[["MIX.SC"]], N.OBS[["MIX.SE"]]))
flush.console()

## ---- phase 1: latent truth of each key ----
## estimator_sim_truth: the population variances of the latent factors of one truth key (Var(a h) for the
## allele fits, Var(a_H) for the summed hybrid, Var(a) for the parents), accumulated over streamed chunks of
## drawn cells, with their Monte Carlo error. Truth comes from the draws, not from the Appendix A formulas.
t0 <- Sys.time()
TRUTH <- parLapply(cl, MY.KEYS, estimator_sim_truth, GRID = GRID, LIB.HYB = LIB.HYB, LIB.PAR = LIB.PAR, DEPTH.REF = DEPTH.REF,
                   NTRUE = NTRUE, SEED.BASE = SEED.BASE, VAR.C = VAR.C, EPS.POOL = EPS.POOL)
names(TRUTH) <- MY.KEYS
cat(sprintf("[%s] truth done: %d keys in %.1f min; worst Monte Carlo SE of a log2 variance: %.4f\n",
            format(Sys.time(), "%H:%M:%S"), length(MY.KEYS), as.numeric(difftime(Sys.time(), t0, units = "mins")),
            max(vapply(TRUTH, `[[`, numeric(1), "mc_se"))))
flush.console()

## ---- phase 2: the rows, chunked so the master can report progress ----
CHUNK.SIZE <- max(NUM.CORES, ceiling(N.POINTS / N.UPDATES))
chunks     <- split(MY.ROWS, ceiling(seq_along(MY.ROWS) / CHUNK.SIZE))
results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## estimator_sim_row: one GRID row (mean x dispersion x intrinsic fraction x mean share x cell fraction x
  ## factor family x scenario). It simulates NREP datasets from the latent model, fits the MLE and the moment
  ## estimator on Hc, He, the summed hybrid and both parents, bootstraps the cells for standard errors and
  ## coverage, and compares everything with the latent truth and the Appendix A predictions.
  results[[k]] <- parLapply(cl, chunks[[k]], estimator_sim_row, GRID = GRID, TRUTH = TRUTH, LIB.HYB = LIB.HYB, LIB.PAR = LIB.PAR,
                            DEPTH.REF = DEPTH.REF, N.OBS = N.OBS, NREP = NREP, NBOOT = NBOOT, SEED.BASE = SEED.BASE,
                            VAR.C = VAR.C, EPS.POOL = EPS.POOL)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (N.POINTS - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS,
              100 * done / N.POINTS, el, eta))
  flush.console()
}

ROW.RESULTS <- do.call(c, results)
SIM.EST     <- do.call(rbind, lapply(ROW.RESULTS, `[[`, "est"))
SIM.AGREE   <- do.call(rbind, lapply(ROW.RESULTS, `[[`, "agree"))
ROWS.PART   <- MY.ROWS
save(SIM.EST, SIM.AGREE, ROWS.PART, file = OUT.FILE)   # the row tables carry everything the truths contribute

cat(sprintf("estimator simulation done: %d rows in %.1f min\n",
            length(ROWS.PART), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
