###############################################################
### power_grid.R
### Reads the power-analysis grid inputs, runs one simulated
### power estimate per (MEAN.READS, N.CELLS, SIZE) row on a FORK
### cluster (each row covers every SIZE.RATIO value internally,
### sharing draws across that sweep), saves the output. Mirrors
### the structure of gene_boot.R.
###
### Inputs  : power_inputs.rda
### Output  : power_output.rda (POWER.MAT, nrow(GRID) rows x
###           length(SIZE.RATIO) columns, in GRID row order)
###
### Progress is printed by the master after each chunk of rows,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')

source("Functions.R")

load("power_inputs.rda")
## GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, SEED.BASE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")

## ---- chunk the grid rows so the master can report progress ----
N.UPDATES  <- 40
N.POINTS   <- nrow(GRID)
CHUNK.SIZE <- max(NUM.CORES, ceiling(N.POINTS / N.UPDATES))
chunks     <- split(seq_len(N.POINTS), ceiling(seq_len(N.POINTS) / CHUNK.SIZE))

cat(sprintf("power grid start: %d rows (each covering %d SIZE.RATIO values), NJ=%d, NI=%d, %d cores, %d chunks\n",
            N.POINTS, length(SIZE.RATIO), NJ, NI, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- parLapply(cl, chunks[[k]], power_grid_row,
                            GRID = GRID, MEAN.READS = MEAN.READS, N.CELLS = N.CELLS,
                            SIZE = SIZE, SIZE.RATIO = SIZE.RATIO,
                            EXPOSURE.CV = EXPOSURE.CV, CELL.RATIO = CELL.RATIO, ALPHA = ALPHA,
                            NJ = NJ, NI = NI, SEED.BASE = SEED.BASE)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (N.POINTS - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d rows (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, N.POINTS,
              100 * done / N.POINTS, el, eta))
  flush.console()
}

POWER.MAT <- do.call(rbind, do.call(c, results))
save(POWER.MAT, file = "power_output.rda")

cat(sprintf("power grid done: %d rows x %d SIZE.RATIO values in %.1f min\n",
            nrow(POWER.MAT), ncol(POWER.MAT), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
