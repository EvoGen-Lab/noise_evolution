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

## power_grid_row: power of one grid row i (mean reads x cell count x burst-frequency size) over the SIZE.RATIO sweep. The reference group's
## counts and fit are drawn once per replicate and reused across the sweep, so the sweep isolates the effect of SIZE.RATIO.
power_grid_row <- function(i, ALPHA, CELL.RATIO, EXPOSURE.CV, GRID, MEAN.READS, N.CELLS, N.MIX, NI, NJ, PI1, SEED.BASE, SIZE, SIZE.RATIO) {
  ## size_log2_ratio: log2 ratio of the NB size (disp, the burst-frequency axis) between two
  ## .fit_one() results, NA if either side is non-positive or non-finite.
  size_log2_ratio <- function(fit_a, fit_b) {
    a <- fit_a[["disp"]]; b <- fit_b[["disp"]]
    if (is.finite(a) && a > 0 && is.finite(b) && b > 0) log2(a) - log2(b) else NA_real_
  }

  row <- GRID[i, ]
  MEAN.READS.X <- MEAN.READS[row$m]
  N.SC.X       <- N.CELLS[row$n]
  N.SE.X       <- round(N.SC.X * CELL.RATIO)
  SIZE.1.X     <- SIZE[row$p]

  set.seed(SEED.BASE + 1e6 + row$n)  #shared draws, keyed only by cell count
  sdlog <- sqrt(log(1 + EXPOSURE.CV^2))
  EXPO.X.LIST <- lapply(seq_len(NJ), function(j) rlnorm(N.SC.X, meanlog = -0.5*sdlog^2, sdlog = sdlog))
  EXPO.Y.LIST <- lapply(seq_len(NJ), function(j) rlnorm(N.SE.X, meanlog = -0.5*sdlog^2, sdlog = sdlog))
  PERM.LIST   <- lapply(seq_len(NI), function(k) sample.int(N.SC.X + N.SE.X))

  set.seed(SEED.BASE + i)  #row-specific draws (X depends on m, p, n)
  X.LIST  <- vector("list", NJ)
  FA.LIST <- vector("list", NJ)
  for (j in seq_len(NJ)) {
    X.LIST[[j]]  <- rnbinom(n = N.SC.X, size = SIZE.1.X, mu = MEAN.READS.X*EXPO.X.LIST[[j]])
    FA.LIST[[j]] <- .fit_one(X.LIST[[j]], EXPO.X.LIST[[j]])  #shared across SIZE.RATIO below
  }

  P.ALL <- matrix(NA_real_, NJ, length(SIZE.RATIO))
  for (qi in seq_along(SIZE.RATIO)) {
    SIZE.2.X <- SIZE.1.X / SIZE.RATIO[qi]
    for (j in seq_len(NJ)) {
      EXPO.Y <- EXPO.Y.LIST[[j]]
      Y  <- rnbinom(n = N.SE.X, size = SIZE.2.X, mu = MEAN.READS.X*EXPO.Y)
      fb <- .fit_one(Y, EXPO.Y)   #same estimator as the null below
      OBS <- size_log2_ratio(FA.LIST[[j]], fb)

      XY <- c(X.LIST[[j]], Y); EXPO.XY <- c(EXPO.X.LIST[[j]], EXPO.Y)
      NULL.DIST <- numeric(NI)
      for (k in seq_len(NI)) {
        sp <- .fit_split(XY, EXPO.XY, PERM.LIST[[k]], N.SC.X)  #same estimator as the observed contrast
        NULL.DIST[k] <- size_log2_ratio(sp$a, sp$b)
      }
      P.ALL[j, qi] <- perm_pval(OBS, NULL.DIST)
    }
  }

  i0     <- which(SIZE.RATIO == 1)
  P.NULL <- P.ALL[, i0]
  n_alt  <- min(NJ, round(NJ * PI1 / (1 - PI1)))
  vapply(seq_along(SIZE.RATIO), function(qi) {
    if (qi == i0) return(mean(P.NULL < ALPHA, na.rm = TRUE))
    mean(replicate(N.MIX, {
      q <- p.adjust(c(P.NULL, P.ALL[sample.int(NJ, n_alt), qi]), method = "BH")
      mean(q[NJ + seq_len(n_alt)] < ALPHA, na.rm = TRUE)
    }), na.rm = TRUE)
  }, numeric(1))
}

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
  ## power_grid_row: power for one (MEAN.READS, N.CELLS, SIZE) row of GRID
  ## across every SIZE.RATIO value at once; the unit of work that
  ## power_grid.R distributes with parLapply, returning one power value per
  ## SIZE.RATIO.
  ##
  ## Design: two independent groups with equal mean reads, a reference (Sc)
  ## group of N.CELLS.X cells with NB size SIZE and a second group of
  ## round(N.CELLS.X * CELL.RATIO) cells with size SIZE / SIZE.RATIO.
  ## Exposures and permutation index sets are seeded from row$n (the cell
  ## count) alone, so rows sharing a cell count draw identical sequences. The
  ## reference group's counts and MLE fit are drawn once per replicate and
  ## reused across the SIZE.RATIO sweep, so the sweep isolates the effect of
  ## SIZE.RATIO.
  ##
  ## Calling power at a Benjamini-Hochberg FDR: every dataset gets a
  ## permutation p-value (NI shuffles of the pooled count/exposure pairs).
  ## For each SIZE.RATIO above 1, the p-values of a random subset of
  ## true-difference datasets (fraction PI1 of the pooled set) are pooled
  ## with the SIZE.RATIO = 1 datasets, which supply the null. BH is applied
  ## to the pooled set and power is the fraction of true-difference datasets
  ## called at q < ALPHA. Averaging over N.MIX random subsets uses every
  ## simulated dataset without further fitting. The SIZE.RATIO = 1 entry is
  ## the raw false-positive rate at p < ALPHA, a calibration check of the
  ## permutation test; SIZE.RATIO therefore includes 1 exactly once.
  results[[k]] <- parLapply(cl, chunks[[k]], power_grid_row, ALPHA = ALPHA, CELL.RATIO = CELL.RATIO, EXPOSURE.CV = EXPOSURE.CV, GRID = GRID, MEAN.READS = MEAN.READS, N.CELLS = N.CELLS, N.MIX = N.MIX, NI = NI, NJ = NJ, PI1 = PI1, SEED.BASE = SEED.BASE, SIZE = SIZE, SIZE.RATIO = SIZE.RATIO)
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
