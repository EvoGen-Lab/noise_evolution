###############################################################
### gene_boot.R
### Reads the bootstrap inputs, runs the paired per-gene
### bootstrap on a FORK cluster, saves the output. Mirrors the
### structure of Permutation_Norm0_V2.R.
###
### Inputs  : gene_boot_inputs.rda (DRAWS.REPS: one set of pre-drawn
###           resamples per bootstrap replicate)
### Output  : gene_boot_output_<k>of<K>.rda for array task k of K, holding
###           BOOT.CONTRASTS (one row per gene) for replicate k. Task k
###           uses DRAWS.REPS[[k]]; replicate 1 is the reported bootstrap
###           and the others check seed adequacy.
###
### Progress is printed by the master after each chunk of genes,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')
library('MASS')          # glm.nb ships with base R; set lib= if your cluster needs it

source("functions.R")
write_pkg_versions("gene_boot")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("gene_boot_inputs.rda")     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS.REPS, N.BOOT

## Job-array position: task k bootstraps replicate k. A finished task leaves its
## output file, so resubmitting the array only computes the missing replicates.
ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
if (N.ARRAY != length(DRAWS.REPS))
  stop(sprintf("array has %d tasks but the inputs hold %d replicates", N.ARRAY, length(DRAWS.REPS)))
DRAWS    <- DRAWS.REPS[[ARRAY.ID]]
OUT.FILE <- sprintf("gene_boot_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40                                              # roughly this many progress lines
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("bootstrap start: replicate %d of %d, %d genes, B=%d, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, length(GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's paired bootstrap: refits all 13 groups on each pre-drawn resample, forms every mode's
  ## mean, bfreq (NB size) and CV2 contrast per replicate, and returns the point estimate, bootstrap SE,
  ## boundary fraction and mean-bfreq draw correlation for the eight reportable modes. Replicates where
  ## any needed fit is non-finite give NA for that contrast and are left out of its SD.
  results[[k]] <- parLapply(cl, chunks[[k]], boot_contrasts_one, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS, draws = DRAWS)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(GENES),
              100 * done / length(GENES), el, eta))
  flush.console()
}

BOOT.CONTRASTS <- do.call(rbind, unlist(results, recursive = FALSE))
save(BOOT.CONTRASTS, file = OUT.FILE)

cat(sprintf("bootstrap done: %d genes in %.1f min [replicate %d]\n",
            nrow(BOOT.CONTRASTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            ARRAY.ID))

stopCluster(cl)
