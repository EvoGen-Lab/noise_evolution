###############################################################
### gene_boot.R
### Reads the bootstrap inputs, runs the paired per-gene
### bootstrap on a FORK cluster, saves the output. Mirrors the
### structure of Permutation_Norm0_V2.R.
###
### Inputs  : gene_boot1_inputs.rda, or gene_boot<tag>_inputs.rda for a
###           given tag arg (e.g. "2" for a second-seed adequacy check)
### Output  : gene_boot1_output.rda, or gene_boot<tag>_output.rda to match
###           (BOOT.CONTRASTS, one row per gene)
###
### Progress is printed by the master after each chunk of genes,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')
library('MASS')          # glm.nb ships with base R; set lib= if your cluster needs it

source("Functions.R")

## Tag arg picks which input/output pair to use, so the same script
## serves both the primary run and any additional-seed adequacy check
## without duplicating the file. Defaults to "1", the primary run.
##   Rscript gene_boot.R    -> gene_boot1_inputs.rda / gene_boot1_output.rda
##   Rscript gene_boot.R 2  -> gene_boot2_inputs.rda / gene_boot2_output.rda
ARGS <- commandArgs(trailingOnly = TRUE)
TAG  <- if (length(ARGS) >= 1) ARGS[1] else "1"
IN.FILE  <- sprintf("gene_boot%s_inputs.rda", TAG)
OUT.FILE <- sprintf("gene_boot%s_output.rda", TAG)

load(IN.FILE)     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS, N.BOOT, SEED.BOOT

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40                                              # roughly this many progress lines
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("bootstrap start: %d genes, B=%d, %d cores, %d chunks\n",
            length(GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- parLapply(cl, chunks[[k]], boot_contrasts_one,
                            mats = CONTRAST.MATS, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, draws = DRAWS)
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

cat(sprintf("bootstrap done: %d genes in %.1f min [tag %s]\n",
            nrow(BOOT.CONTRASTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            TAG))

stopCluster(cl)
