###############################################################
### power_modes_grid.R
### Reads the cis/trans power-grid inputs, runs one simulated power
### estimate per grid row (mean reads x size x phi x f x design x depth
### ratio) on a FORK cluster, saves the output. Each row covers every
### SIZE.RATIO value, the cis contrast and the trans contrast under both
### null versions on paired-allele hybrid data (power_modes_row()).
### Mirrors power_grid.R.
###
### Inputs  : power_modes_inputs.rda (MODES.INPUTS, a list holding GRID,
###           MEAN.READS, SIZE, PHI, FRAC, DESIGN, DEPTH.RATIO, SIZE.RATIO,
###           EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE)
### Output  : power_modes_output_<k>of<K>.rda for array task k of K, holding
###           ROWS.PART (the GRID rows this task computed) and RESULTS.PART
###           (one power_modes_row() result per row). Rows are dealt out
###           round-robin, and analysis/power_analysis.R Section 10.10
###           reassembles them in GRID order.
###
### Progress is printed by the master after each chunk of rows,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')

source("functions.R")
source("functions_power.R")

load("power_modes_inputs.rda")
list2env(MODES.INPUTS, environment())
## GRID, MEAN.READS, SIZE, PHI, FRAC, DESIGN, DEPTH.RATIO, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

## Job-array position. Each row is seeded by its own index inside power_modes_row(),
## so results are identical however the rows are split across tasks.
ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE <- sprintf("power_modes_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the whole array only
## computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

cl <- makeCluster(NUM.CORES, type = "FORK")

## ---- chunk the grid rows so the master can report progress ----
N.UPDATES  <- 40
MY.ROWS    <- which((seq_len(nrow(GRID)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
N.POINTS   <- length(MY.ROWS)
CHUNK.SIZE <- max(NUM.CORES, ceiling(N.POINTS / N.UPDATES))
chunks     <- split(MY.ROWS, ceiling(seq_along(MY.ROWS) / CHUNK.SIZE))

cat(sprintf("power modes grid start: array task %d of %d, %d rows (each covering %d SIZE.RATIO values, 3 modes), NJ=%d, NI=%d, PI1=%.2f, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, N.POINTS, length(SIZE.RATIO), NJ, NI, PI1, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## power_modes_row: power of the cis and trans contrasts for one row of GRID across every SIZE.RATIO
  ## value; the unit of work that power_modes_grid.R distributes with parLapply. Paired-allele hybrid data
  ## (shared extrinsic factor with fraction phi), the hybrid cells split into cis cells (f) and trans cells
  ## (1 - f), parents unpaired; cis is nulled by the within-cell allele swap and trans by pooled relabelings in
  ## two versions (independent allele shuffles and a pairing-preserving swap), all through the helpers
  ## the gene-level permutation uses. Returns power by mode, ratio and axis (mu, disp, cv2) at a BH FDR,
  ## with ratio 1 as the type I error, and the SD of the observed contrast under the null.
  results[[k]] <- parLapply(cl, chunks[[k]], power_modes_row, ALPHA = ALPHA, CELL.RATIO = CELL.RATIO, DEPTH.RATIO = DEPTH.RATIO, DESIGN = DESIGN, EXPOSURE.CV = EXPOSURE.CV, FRAC = FRAC, GRID = GRID, MEAN.READS = MEAN.READS, N.MIX = N.MIX, NI = NI, NJ = NJ, PHI = PHI, PI1 = PI1, SEED.BASE = SEED.BASE, SIZE = SIZE, SIZE.RATIO = SIZE.RATIO)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (N.POINTS - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS,
              100 * done / N.POINTS, el, eta))
  flush.console()
}

RESULTS.PART <- unlist(results, recursive = FALSE)
ROWS.PART    <- MY.ROWS
save(RESULTS.PART, ROWS.PART, file = OUT.FILE)

cat(sprintf("power modes grid done: %d rows in %.1f min\n",
            length(RESULTS.PART), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
