###############################################################
### power_grid.R
### Reads the power-analysis grid inputs, runs one simulated
### power estimate per (MEAN.READS, N.CELLS, SIZE) row on a FORK
### cluster (each row covers every SIZE.RATIO value internally,
### sharing draws across that sweep), saves the output. Mirrors
### the structure of gene_boot.R.
###
### Inputs  : power_inputs.rda
### Output  : power_output_<k>of<K>.rda for array task k of K, holding
###           ROWS.PART (the GRID rows this task computed) and POWER.PART
###           (those rows x length(SIZE.RATIO) columns). Rows are dealt out
###           round-robin, so every task gets a similar mix of cheap and
###           costly rows, and analysis.R Section 10.2 reassembles GRID order.
###
### Progress is printed by the master after each chunk of rows,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')

source("functions.R")

load("power_inputs.rda")
## GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

## Job-array position. Each row is seeded by its own index inside power_grid_row(),
## so results are identical however the rows are split across tasks.
ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE <- sprintf("power_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
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

cat(sprintf("power grid start: array task %d of %d, %d rows (each covering %d SIZE.RATIO values), NJ=%d, NI=%d, PI1=%.2f, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, N.POINTS, length(SIZE.RATIO), NJ, NI, PI1, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- parLapply(cl, chunks[[k]], power_grid_row,
                            GRID = GRID, MEAN.READS = MEAN.READS, N.CELLS = N.CELLS,
                            SIZE = SIZE, SIZE.RATIO = SIZE.RATIO,
                            EXPOSURE.CV = EXPOSURE.CV, CELL.RATIO = CELL.RATIO, ALPHA = ALPHA,
                            NJ = NJ, NI = NI, SEED.BASE = SEED.BASE, PI1 = PI1, N.MIX = N.MIX)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (N.POINTS - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS,
              100 * done / N.POINTS, el, eta))
  flush.console()
}

POWER.PART <- do.call(rbind, do.call(c, results))
ROWS.PART  <- MY.ROWS
save(POWER.PART, ROWS.PART, file = OUT.FILE)

cat(sprintf("power grid done: %d rows x %d SIZE.RATIO values in %.1f min\n",
            nrow(POWER.PART), ncol(POWER.PART), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
