###############################################################
### gene_boot.R
### Reads the bootstrap inputs, runs the paired per-gene bootstrap on a FORK cluster, saves the output. Two stages (Methods: Shared shift), set by the
### environment variable GENE_BOOT_STAGE (default 1):
###
### Stage 1 (job array, gene_boot.sub)
###   Inputs  : gene_boot_inputs.rda (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, N.BOOT, DRAWS.REPS: one set of pre-drawn resamples per bootstrap
###             replicate; CONTRAST.ETAS: the split variance of every cell of every dataset, 0 for the parents, which every replicate reads for its resampled
###             cells; GATES: the contrasts each gene's observed fits define)
###   Output  : gene_boot_output_<k>of<K>.rda for array task k of K, holding BOOT.CONTRASTS (one row per gene) for replicate k. Task k uses DRAWS.REPS[[k]];
###             replicate 1 is the reported bootstrap and the others check seed adequacy. Each row holds every contrast's estimate, bootstrap SE and
###             percentile interval, the sign-test p-value of every NB contrast (sign_boot_p(); a failed refit counts as a zero dispersion) and the moment
###             estimates with their flags. analysis.R fits the shared shift once from replicate 1's observed mid log2 dispersion contrasts.
### Stage 2 (one task, gene_boot_shift.sub; after analysis.R has written the shift)
###   Inputs  : gene_boot_inputs.rda and gene_shift_inputs.rda (DOM.SHIFT: the shared shift Xi at each gene's log2 mu_H from the observed curve, named by
###             gene, held fixed in every replicate)
###   Output  : gene_boot_shift_output.rda, holding BOOT.CONTRASTS for replicate 1 with every hybrid-allele dispersion divided by 2^Xi: the stage-1 columns
###             (unchanged), the shifted contrasts a<quantity>_<mode>_* and the degrees of dominance k_<quantity>_* and kt_<quantity>_*.
###
### Progress is printed by the master after each chunk of genes, so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')
library('MASS')          # glm.nb ships with base R; set lib= if your cluster needs it

source("functions.R")
STAGE <- as.integer(Sys.getenv("GENE_BOOT_STAGE", unset = "1"))
stopifnot(STAGE %in% 1:2)
write_pkg_versions(if (STAGE == 1) "gene_boot" else "gene_boot_shift")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("gene_boot_inputs.rda")     # GENES is the inference set (every gene for which at least one contrast is defined). CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS.REPS, N.BOOT, CONTRAST.ETAS, GATES
SHIFT <- NULL
if (STAGE == 2) {
  load("gene_shift_inputs.rda")  # DOM.SHIFT
  stopifnot(is.numeric(DOM.SHIFT), all(names(DOM.SHIFT) %in% GENES), all(is.finite(DOM.SHIFT)))
  SHIFT <- DOM.SHIFT
}

## Job-array position: stage-1 task k bootstraps replicate k. A finished task leaves its output file, so resubmitting the array only computes the missing
## replicates. Stage 2 runs replicate 1 alone.
ARRAY.ID <- if (STAGE == 1) as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1")) else 1L
N.ARRAY  <- if (STAGE == 1) as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1")) else length(DRAWS.REPS)
if (N.ARRAY != length(DRAWS.REPS))
  stop(sprintf("array has %d tasks but the inputs hold %d replicates", N.ARRAY, length(DRAWS.REPS)))
DRAWS    <- DRAWS.REPS[[ARRAY.ID]]
OUT.FILE <- if (STAGE == 1) sprintf("gene_boot_output_%dof%d.rda", ARRAY.ID, N.ARRAY) else "gene_boot_shift_output.rda"
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40                                              # roughly this many progress lines
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("bootstrap start: stage %d, replicate %d of %d, %d genes, B=%d, %d cores, %d chunks\n",
            STAGE, ARRAY.ID, N.ARRAY, length(GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's paired bootstrap: refits all 18 contrast groups on each pre-drawn resample with the split variances of its cells, forms every contrast per
  ## replicate, and returns the estimate, bootstrap SE, percentile interval and sign-test p-value of every contrast (stage 2: also after the shared shift).
  results[[k]] <- parLapply(cl, chunks[[k]], boot_contrasts_one, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS, draws = DRAWS, etas = CONTRAST.ETAS,
                            gates = GATES, methods = c("NB", "weighted"), shift = SHIFT)
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
