###############################################################
### gene_pilot.R
### Reads the pilot inputs, runs the per-gene internal-pilot
### bootstrap (SE of log mu and log disp, undivided hybrid and
### parents) on a FORK cluster, saves the output. Mirrors the
### structure of gene_boot.R exactly; this is the same
### kind of per-gene, per-replicate bootstrap, just over four
### datasets instead of thirteen, and used only to set the Nc/Nt
### split fraction rather than to test any gene for significance.
###
### Inputs  : gene_pilot_inputs.rda  (PILOT.MATS, PILOT.EXPOS,
###           PILOT.GENES, N.BOOT, SEED.BOOT)
### Output  : gene_pilot_output.rda  (PILOT.SE, one row per gene)
###############################################################

library('parallel')
library('MASS')

source("functions.R")
load("gene_pilot_inputs.rda")  # PILOT.MATS, PILOT.EXPOS, PILOT.GENES, N.BOOT, SEED.BOOT

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(PILOT.GENES) / N.UPDATES))
chunks     <- split(PILOT.GENES, ceiling(seq_along(PILOT.GENES) / CHUNK.SIZE))

cat(sprintf("pilot start: %d genes, B=%d, %d cores, %d chunks\n",
            length(PILOT.GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- parLapply(cl, chunks[[k]], pilot_split_se_one,
                            mats = PILOT.MATS, expos = PILOT.EXPOS, B = N.BOOT, seed = SEED.BOOT)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(PILOT.GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(PILOT.GENES),
              100 * done / length(PILOT.GENES), el, eta))
  flush.console()
}

PILOT.SE <- do.call(rbind, unlist(results, recursive = FALSE))
save(PILOT.SE, file = "gene_pilot_output.rda")

cat(sprintf("pilot done: %d genes in %.1f min\n",
            nrow(PILOT.SE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
