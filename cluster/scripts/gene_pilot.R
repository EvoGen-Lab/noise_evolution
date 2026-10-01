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
### Output  : gene_pilot_output.rda  (PILOT.SE, one row per gene: bootstrap SEs of
###           log mu and log disp for the four datasets, plus the bootstrap covariance of
###           the two hybrid alleles' estimates, which share their resampled cells)
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
  ## One gene's bootstrap SE of log(mu) and log(size) for the four datasets needed to compute
  ## A and B, plus the bootstrap covariance between the two hybrid alleles. HYB.SC and HYB.SE
  ## are refit on the SAME resampled cell indices in each replicate because the two alleles
  ## are measured in the same cells; their log estimates therefore co-vary across replicates
  ## (HYB_logmu_cov, HYB_logdisp_cov), and the variance of the allele contrast is
  ## var(SC) + var(SE) - 2 cov(SC, SE), which estimate_f_star() uses for A.
  ##
  ## The RNG is seeded from a hash of the gene's own name and the resamples are
  ## drawn inside this function. Each gene's draws therefore depend only on its
  ## name, so results are identical for any gene chunking or core count of the
  ## FORK cluster that distributes the genes.
  results[[k]] <- parLapply(cl, chunks[[k]], pilot_split_se_one, B = N.BOOT, expos = PILOT.EXPOS, mats = PILOT.MATS, seed = SEED.BOOT)
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
