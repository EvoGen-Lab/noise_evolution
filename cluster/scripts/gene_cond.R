###############################################################
### gene_cond.R
### State-conditioned sensitivity (analysis.R Section 2.5): every gene's contrasts refitted with and without the cell-level score of the derepression and
### stress program in the NB mean, on the bootstrap draws of the reported replicate, on a FORK cluster.
###
### Inputs  : gene_cond_inputs.rda (CONTRAST.MATS, CONTRAST.EXPOS, GENES, DRAWS: the pre-drawn resamples of bootstrap replicate 1; SHARE.CELLS and SHARE.OBS:
###           the transcriptome share inputs that every replicate refits; GATES: the contrasts each gene's fits define; STATE.COND: per tube each cell's
###           summed program count and depth factor, the program genes and the pseudocount; COND.KEY: the fingerprint of the inputs; N.COND.ARRAY)
### Output  : gene_cond_output_<k>of<K>.rda for array task k of K, holding COND.CONTRASTS (one row per gene of block k, cond_contrasts_one()) and COND.KEY.
###           The genes are split into K consecutive blocks; a finished task leaves its file, so resubmitting the array computes only the missing blocks.
###
### Progress is printed by the master after each chunk of genes, so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')

source("functions.R")
write_pkg_versions("gene_cond")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("gene_cond_inputs.rda")      # CONTRAST.MATS, CONTRAST.EXPOS, GENES, DRAWS, SHARE.OBS, SHARE.CELLS, GATES, STATE.COND, COND.KEY, N.COND.ARRAY

## Job-array position: task k runs gene block k.
ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
if (N.ARRAY != N.COND.ARRAY)
  stop(sprintf("array has %d tasks but the inputs expect %d gene blocks", N.ARRAY, N.COND.ARRAY))
OUT.FILE <- sprintf("gene_cond_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }
BLOCK <- GENES[cut(seq_along(GENES), N.ARRAY, labels = FALSE) == ARRAY.ID]

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
cl <- makeCluster(NUM.CORES, type = "FORK")

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(BLOCK) / N.UPDATES))
chunks     <- split(BLOCK, ceiling(seq_along(BLOCK) / CHUNK.SIZE))

cat(sprintf("state-conditioned bootstrap start: block %d of %d, %d genes, B=%d, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, length(BLOCK), length(DRAWS), NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene: every group refitted on every resample without and with the program score in the mean, and every contrast built from both sets of fits.
  results[[k]] <- parLapply(cl, chunks[[k]], cond_contrasts_one, expos = CONTRAST.EXPOS, mats = CONTRAST.MATS, draws = DRAWS, share_cells = SHARE.CELLS,
                            share_obs = SHARE.OBS, state = STATE.COND, gates = GATES)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(BLOCK) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(BLOCK), 100 * done / length(BLOCK), el, eta))
  flush.console()
}

COND.CONTRASTS <- do.call(rbind, unlist(results, recursive = FALSE))
save(COND.CONTRASTS, COND.KEY, file = OUT.FILE)
cat(sprintf("state-conditioned bootstrap done: %d genes in %.1f min [block %d]\n", nrow(COND.CONTRASTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins")), ARRAY.ID))

stopCluster(cl)
