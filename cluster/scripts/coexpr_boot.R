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
###           CB (CB$total, CB$cis, CB$trans, CB$dpar_sc, CB$dpar_se,
###           CB$lambda) for replicate k. Task k uses DRAWS.COEXPR.REPS[[k]];
###           replicate 1 is the reported bootstrap and the others check seed
###           adequacy.
###
### Progress is printed by the master after each chunk of draws,
### so it lands in the job log even though the workers are forked.
###
### Memory note: each chunk's draws are folded into a running sum and
### sum-of-squares immediately (coexpr_acc_update), not accumulated into
### one B x pairs matrix at the end. Holding all B draws simultaneously
### scales as O(B x pairs), and pairs scales as p^2, so a larger gene set
### can exceed available memory well before the CPU work finishes, which
### is exactly what happened at the previous, smaller gene count's scale
### once B and pairs both grew.
###############################################################

library('parallel')

source("functions.R")
write_pkg_versions("coexpr_boot")   # R and package versions of this job, compared locally by check_pkg_versions()

## The five decompositions the co-expression bootstrap tracks, in the order each draw returns them.
.COEXPR_PARTS <- c("total", "cis", "trans", "dpar_sc", "dpar_se")

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

## Streaming accumulator for the co-expression bootstrap SE. Holding every draw (a B x pairs matrix
## per contrast) costs O(B x pairs) memory and pairs grows as p^2. A standard deviation needs only
## sum(x) and sum(x^2) per pair, so folding in each chunk's draws keeps peak memory at
## O(chunk_size x pairs) and the raw draws can be discarded. p: number of genes (nrow of a RESID
## matrix), used to build the upper-triangle index once. Tracks the parts in .COEXPR_PARTS.
ACC <- local({
  p <- nrow(RESID$MIX.SC)
  up <- which(upper.tri(matrix(0, p, p)))
  z  <- numeric(length(up))
  sums <- unlist(lapply(.COEXPR_PARTS, function(nm) setNames(list(z, z), paste0(c("sum_", "sumsq_"), nm))), recursive = FALSE)
  c(list(up = up, n = 0L), sums)
})

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One bootstrap draw per task. A draw (from make_coexpr_draws()) resamples the columns of each RESID
  ## dataset. HYB.SC, HYB.SE and HYB.COMB take the same H draw because they are the same cells, and the
  ## per-gene ploidy factors are fixed, so every draw rescales Rhyb identically. Each task returns
  ## total/cis/trans/dpar_sc/dpar_se at the upper-triangle pair positions.
  chunk.results <- parLapply(cl, DRAWS.COEXPR[chunks[[k]]], coexpr_bootstrap_one, resid = RESID)
  ## Fold this chunk's draws into the running sums and draw count; the raw draws are then discarded.
  ACC <- local({
    acc <- ACC
    acc$n <- acc$n + length(chunk.results)
    for (nm in .COEXPR_PARTS) {
      m <- do.call(rbind, lapply(chunk.results, `[[`, nm))
      acc[[paste0("sum_", nm)]]   <- acc[[paste0("sum_", nm)]]   + colSums(m)
      acc[[paste0("sumsq_", nm)]] <- acc[[paste0("sumsq_", nm)]] + colSums(m^2)
    }
    acc
  })
  rm(chunk.results); gc(FALSE)   # drop this chunk's raw draws before the next one

  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (B - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d draws (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, B, 100 * done / B, el, eta))
  flush.console()
}

## Finishes the accumulator into the CB structure: for each part a data.frame of
## gene_i/gene_j/est/se/z/p, plus lambda. est is the point estimate pt, se the bootstrap SD from the
## running sums, z = est/se and p the two-sided normal tail.
CB <- local({
  acc <- ACC
  gn <- rownames(RESID$MIX.SC); p <- length(gn)
  ij <- arrayInd(acc$up, c(p, p))
  base <- data.frame(gene_i = gn[ij[, 1]], gene_j = gn[ij[, 2]])
  parts <- lapply(setNames(.COEXPR_PARTS, .COEXPR_PARTS), coexpr_part_table, point = COEXPR.POINT, acc = acc, base = base)
  c(parts, list(lambda = COEXPR.POINT$lambda))
})
save(CB, file = OUTPUT.FILE)

cat(sprintf("coexpression bootstrap done: %d pairs in %.1f min [replicate %d]\n",
            nrow(CB$total),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            ARRAY.ID))

stopCluster(cl)
