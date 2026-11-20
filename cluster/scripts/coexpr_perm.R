###############################################################
### coexpr_perm.R
### Reads the co-expression null inputs, runs the joint null draws on a FORK cluster, and saves the top-N.KEEP squared-eigenvalue spectrum of total, cis,
### trans and dom for every draw, for rank-matched testing (real rank k against the null's own rank k) of every candidate axis (Methods: Co-expression, Axes of
### divergence). Each draw's result is a handful of numbers (4 x N.KEEP), collected directly and stacked into matrices at the end.
###
### The total null pools the parental cells and reassigns them to two groups of the original sizes, the cis null swaps the residuals of the two alleles in a
### random subset of hybrid cells, and the trans null is their difference. The dominance null works within each genome: the Sc parent's cells are pooled with
### the Sc residuals of the hybrid cells and the Se parent's with the Se residuals, each pool is reassigned to parental and hybrid groups of the original sizes,
### and dom is recomputed with the same shared-shift adjustment as the observed matrix.
###
### Inputs  : coexpr_perm_inputs.rda  (RESID with its pair deciles attr "dom_bins", N.SC, N.SE, N.KEEP, DRAWS.PERM.COEXPR, DRAWS.DOM.COEXPR)
### Output  : coexpr_perm_output.rda  (NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS, NULL.DOM.RANKS: B x N.KEEP matrices, column k = the draw's
###                                    k-th largest-magnitude squared eigenvalue)
###############################################################

library('parallel')

source("functions.R")
write_pkg_versions("coexpr_perm")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

## Optional args: Rscript coexpr_perm.R [input.rda] [output.rda]
ARGS        <- commandArgs(trailingOnly = TRUE)
INPUT.FILE  <- if (length(ARGS) >= 1) ARGS[1] else "coexpr_perm_inputs.rda"
OUTPUT.FILE <- if (length(ARGS) >= 2) ARGS[2] else "coexpr_perm_output.rda"

load(INPUT.FILE)     # RESID, N.SC, N.SE, N.KEEP, DRAWS.PERM.COEXPR, DRAWS.DOM.COEXPR
stopifnot(length(DRAWS.PERM.COEXPR) == length(DRAWS.DOM.COEXPR), !is.null(attr(RESID, "dom_bins")))

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")

WORKER.PIDS <- unlist(clusterEvalQ(cl, Sys.getpid()))
cat(sprintf("requested %d cores, got %d worker PIDs: %s\n",
            NUM.CORES, length(WORKER.PIDS), paste(WORKER.PIDS, collapse = ", ")))
flush.console()

## ---- chunk the draws so the master can report progress ----
N.UPDATES  <- 40
B          <- length(DRAWS.PERM.COEXPR)
DRAWS      <- Map(list, perm = DRAWS.PERM.COEXPR, dom = DRAWS.DOM.COEXPR)
CHUNK.SIZE <- max(NUM.CORES, ceiling(B / N.UPDATES))
chunks     <- split(seq_len(B), ceiling(seq_len(B) / CHUNK.SIZE))

cat(sprintf("coexpression null start: %d genes, B=%d, %d cores, %d chunks\n",
            nrow(RESID$MIX.SC), B, NUM.CORES, length(chunks)))
flush.console()

RESULTS <- vector("list", B)

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One draw of the four null spectra (coexpr_perm_one()), for rank-matched testing of every candidate axis: observed rank k against the null's own rank k.
  chunk.results <- parLapply(cl, DRAWS[chunks[[k]]], coexpr_perm_one, n_keep = N.KEEP, nSC = N.SC, resid = RESID)
  RESULTS[chunks[[k]]] <- chunk.results
  rm(chunk.results); gc(FALSE)

  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (B - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d draws (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, B, 100 * done / B, el, eta))
  flush.console()
}

NULL.TOTAL.RANKS   <- do.call(rbind, lapply(RESULTS, `[[`, "total"))
NULL.CIS.RANKS     <- do.call(rbind, lapply(RESULTS, `[[`, "cis"))
NULL.TRANS.RANKS   <- do.call(rbind, lapply(RESULTS, `[[`, "trans"))
NULL.DOM.RANKS     <- do.call(rbind, lapply(RESULTS, `[[`, "dom"))

save(NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS, NULL.DOM.RANKS, file = OUTPUT.FILE)

cat(sprintf("coexpression null done: %d draws x %d ranks in %.1f min\n",
            B, N.KEEP,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
