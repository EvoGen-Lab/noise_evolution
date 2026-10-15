###############################################################
### raw_trans_null_grid.R
### Reads the trans-null validation inputs and runs raw_trans_null_row() for each scenario row on a FORK cluster: the type I error and
### power of three nulls (pooled relabeling with independent or paired allele shuffles, and a recentered bootstrap) for the trans
### contrast of phi on the raw and the log2 scale. A separate study; it reads no pipeline result. Mirrors power_modes_grid.R.
###
### Inputs  : raw_trans_null_inputs.rda, written by analysis/raw_trans_null_inputs.R
###           (RAW.TRANS.INPUTS: GRID, SHARED, NJ, NI, NBOOT, ALPHA, SEED.BASE)
### Output  : raw_trans_null_output_<k>of<K>.rda for array task k of K, holding ROWS.PART (the GRID rows this task computed) and
###           RESULTS.PART (one raw_trans_null_row() data.frame per row). analysis/raw_trans_null_summary.R reassembles them.
### Each row is seeded by its own index inside raw_trans_null_row(), so results do not depend on how rows are split across tasks.
###############################################################

library('parallel')

source("functions.R")
source("functions_power.R")
write_pkg_versions("raw_trans_null_grid")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("raw_trans_null_inputs.rda")
list2env(RAW.TRANS.INPUTS, environment())   # GRID, SHARED, NJ, NI, NBOOT, ALPHA, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
ARRAY.ID  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE  <- sprintf("raw_trans_null_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the whole array only computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

cl <- makeCluster(NUM.CORES, type = "FORK")

N.UPDATES  <- 20
MY.ROWS    <- which((seq_len(nrow(GRID)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
N.POINTS   <- length(MY.ROWS)
CHUNK.SIZE <- max(NUM.CORES, ceiling(N.POINTS / N.UPDATES))
chunks     <- split(MY.ROWS, ceiling(seq_along(MY.ROWS) / CHUNK.SIZE))

cat(sprintf("raw trans null grid start: array task %d of %d, %d rows, NJ=%d, NI=%d, NBOOT=%d, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, N.POINTS, NJ, NI, NBOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## raw_trans_null_row: one scenario row. It simulates NJ datasets with known phi in the four groups, tests the raw and the log2 trans
  ## contrast of each against three nulls (pooled independent, pooled paired, recentered bootstrap) and returns the rejection rates at
  ## ALPHA together with the observed contrast's mean and SD and the nulls' SD.
  results[[k]] <- parLapply(cl, chunks[[k]], raw_trans_null_row, GRID = GRID, SHARED = SHARED, NJ = NJ, NI = NI, NBOOT = NBOOT, ALPHA = ALPHA, SEED.BASE = SEED.BASE)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS, 100 * done / N.POINTS, el, (N.POINTS - done) * (el / done)))
  flush.console()
}

RESULTS.PART <- unlist(results, recursive = FALSE)
ROWS.PART    <- MY.ROWS
save(RESULTS.PART, ROWS.PART, file = OUT.FILE)
cat(sprintf("raw trans null grid done: %d rows in %.1f min\n", length(RESULTS.PART), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
