###############################################################
### gene_perm.R
### Reads the same inputs as the bootstrap, runs the paired
### per-gene permutation null on a FORK cluster, saves the output.
### Mirrors gene_boot.R.
###
### Inputs  : gene_perm_inputs.rda  (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS GENES, ...)
### Output  : gene_perm_output.rda
###           PERM.RESULTS, one row per gene: mean_<mode>_obs / _p,
###           bfreq_<mode>_obs / _p (the NB dispersion parameter,
###           reported under its biological name), bsize_<mode>_obs / _p
###           (mean minus bfreq, its own null), kbal_<mode>_obs / _p
###           (bfreq minus bsize, its own null), and cv2_<mode>_obs / _p
###           (1/mean + 1/bfreq, its own null).
###           For dpar_sc and dpar_se, the noise quantities also get
###           _p_ploidy and _p_ind columns, testing the contrast against
###           the ploidy-expected value using this gene's own shift and
###           using the independent-allele shift of 1, plus the
###           ploidy_shift column recording the shift used.
###
### Each contrast is nulled by its own label shuffle. total, cis, and
### trans are zero-centred under the null they test. dpar tests
### exchangeability of hybrid and parent cells, which is the null the
### ploidy-shifted observed value is read against (ploidy_shift() and
### ploidy_adjust_dpar() in Functions.R).
###############################################################

library('parallel')
library('MASS')

source("Functions.R")
load("gene_perm_inputs.rda")     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, SEED.PERM, PLOIDY.SHIFT

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
cl <- makeCluster(4L, type = "PSOCK")
clusterExport(cl, c("CONTRAST.MATS", "CONTRAST.EXPOS", "CONTRAST.FITS", "PERMS", "PLOIDY.SHIFT", "permute_contrasts_one"))
clusterEvalQ(cl, { library(MASS); source("Functions.R") })

N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("permutation start: %d genes, NPERM=%d, %d cores, %d chunks\n",
            length(GENES), N.PERM, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
 results[[k]] <- parLapply(cl, chunks[[k]], permute_contrasts_one,
                            mats = CONTRAST.MATS, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, perms = PERMS, ploidy_shift = PLOIDY.SHIFT)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(GENES),
              100 * done / length(GENES), el, eta))
  flush.console()
}

PERM.RESULTS <- do.call(rbind, unlist(results, recursive = FALSE))
save(PERM.RESULTS, file = "gene_perm_output.rda")

cat(sprintf("permutation done: %d genes in %.1f min\n",
            nrow(PERM.RESULTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
