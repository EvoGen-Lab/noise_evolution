###############################################################
### coexpr_boot.R
### Reads the co-expression bootstrap inputs, runs the B cell
### resamples on a FORK cluster, assembles and saves CB. Mirrors
### the structure of gene_boot.R, but chunks over draws rather
### than over genes, since the gene set here is fixed (CO.GENES)
### and the expensive step is one shrink_cor per dataset per draw.
###
### Inputs  : coexpr_boot_inputs.rda (RESID, COEXPR.POINT, DRAWS.COEXPR.REPS:
###           one set of pre-drawn resamples per bootstrap replicate)
### Output  : coexpr_boot_output_<k>of<K>.rda for array task k of K, holding
###           CB (CB$total, CB$cis, CB$trans, CB$dom, CB$inh_c, CB$inh_e: estimate, SE, z and p per
###           pair; CB$k and CB$kt: the degrees of dominance k = 2 dom / total and k_trans = 2 dom / trans with
###           their bootstrap distribution at -1, 0 and 1; CB$lambda) for replicate k. Task k uses DRAWS.COEXPR.REPS[[k]];
###           replicate 1 is the reported bootstrap and the others check seed
###           adequacy.
###
### Progress is printed by the master after each chunk of draws,
### so it lands in the job log even though the workers are forked.
###
### Memory note: each chunk's draws are folded into a running sum and
### sum-of-squares and the k counts immediately (coexpr_acc_fold()), not accumulated into
### one B x pairs matrix at the end. Holding all B draws simultaneously
### scales as O(B x pairs), and pairs scales as p^2, so a larger gene set
### can exceed available memory well before the CPU work finishes, which
### is exactly what happened at the previous, smaller gene count's scale
### once B and pairs both grew.
###############################################################

library('parallel')

source("functions.R")
write_pkg_versions("coexpr_boot")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("coexpr_boot_inputs.rda")     # RESID, COEXPR.POINT, DRAWS.COEXPR.REPS

## Job-array position: task k bootstraps replicate k, the same convention gene_boot.R
## uses. A finished task leaves its output file, so resubmitting the array only
## computes the missing replicates.
ARRAY.ID    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
if (N.ARRAY != length(DRAWS.COEXPR.REPS))
  stop(sprintf("array has %d tasks but the inputs hold %d replicates", N.ARRAY, length(DRAWS.COEXPR.REPS)))
DRAWS.COEXPR <- DRAWS.COEXPR.REPS[[ARRAY.ID]]
OUTPUT.FILE <- sprintf("coexpr_boot_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
if (file.exists(OUTPUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUTPUT.FILE)); quit(save = "no") }

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")

WORKER.PIDS <- unlist(clusterEvalQ(cl, Sys.getpid()))
cat(sprintf("requested %d cores, got %d worker PIDs: %s\n",
            NUM.CORES, length(WORKER.PIDS), paste(WORKER.PIDS, collapse = ", ")))
flush.console()

## ---- chunk the draws so the master can report progress ----
N.UPDATES  <- 40                                              # roughly this many progress lines
B          <- length(DRAWS.COEXPR)
CHUNK.SIZE <- max(NUM.CORES, ceiling(B / N.UPDATES))
chunks     <- split(seq_len(B), ceiling(seq_len(B) / CHUNK.SIZE))

cat(sprintf("coexpression bootstrap start: replicate %d of %d, %d genes, B=%d, %d cores, %d chunks\n",
            ARRAY.ID, N.ARRAY, nrow(RESID$MIX.SC), B, NUM.CORES, length(chunks)))
flush.console()

## Streaming accumulator (coexpr_acc_new(), coexpr_acc_fold()): per contrast the running sum and sum of squares, and the counts that give the bootstrap
## distribution of k and k_trans, so peak memory stays at O(chunk size x pairs) and the raw draws are discarded.
ACC <- coexpr_acc_new(choose(nrow(RESID$MIX.SC), 2))

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One bootstrap draw per task. A draw (from make_coexpr_draws()) resamples the columns of each RESID dataset; HYB.SC and HYB.SE take the same H draw
  ## because they are the same cells, and the shared shift of dom is recomputed in every draw on the fixed deciles. Each task returns every contrast at the
  ## upper-triangle pair positions, folded into the accumulator.
  chunk.results <- parLapply(cl, DRAWS.COEXPR[chunks[[k]]], coexpr_bootstrap_one, resid = RESID)
  ACC <- coexpr_acc_fold(ACC, chunk.results)
  rm(chunk.results); gc(FALSE)   # drop this chunk's raw draws before the next one

  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (B - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d draws (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, B, 100 * done / B, el, eta))
  flush.console()
}

## Finishes the accumulator into the CB structure: for each contrast a data.frame of gene_i/gene_j/est/se/z/p (est the point estimate, se the bootstrap SD
## from the running sums, z = est / se and p the two-sided normal tail), the k and k_trans tables, the decile shifts and lambda.
CB <- local({
  gn <- rownames(RESID$MIX.SC); p <- length(gn)
  ij <- arrayInd(which(upper.tri(matrix(0, p, p))), c(p, p))
  base <- data.frame(gene_i = gn[ij[, 1]], gene_j = gn[ij[, 2]])
  parts <- lapply(setNames(.COEXPR.PARTS, .COEXPR.PARTS), coexpr_part_table, point = COEXPR.POINT, acc = ACC, base = base)
  c(parts, list(k = coexpr_k_table("k", COEXPR.POINT, ACC, base), kt = coexpr_k_table("kt", COEXPR.POINT, ACC, base),
                shift = COEXPR.POINT$shift, lambda = COEXPR.POINT$lambda))
})
save(CB, file = OUTPUT.FILE)

cat(sprintf("coexpression bootstrap done: %d pairs in %.1f min [replicate %d]\n",
            nrow(CB$total),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            ARRAY.ID))

stopCluster(cl)
