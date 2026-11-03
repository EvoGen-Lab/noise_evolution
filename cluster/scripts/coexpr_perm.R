###############################################################
### coexpr_perm.R
### Reads the co-expression null inputs, runs the joint
### total/cis/par_sc/par_se draws on a FORK cluster, and saves the
### top-N.KEEP squared-eigenvalue spectrum of total, cis, trans, par_sc,
### and par_se for every draw. Used for rank-matched testing (real rank
### k against the null's own rank k) of every candidate axis in all five
### decompositions, replacing the earlier single shared-threshold
### approach. Mirrors the structure of coexpr_boot.R, but each
### draw's result is a handful of numbers (5 x N.KEEP) rather than a
### full pairwise matrix, so no incremental accumulator is needed;
### results are collected directly and stacked into matrices at the end.
###
### The total and cis nulls are permutation nulls. The par_sc and par_se
### nulls are bootstrap nulls with H0 imposed (Beran and Srivastava 1985):
### each group resamples its own cells after recoloring to a common
### correlation structure on the per-genome (ploidy-rescaled) scale, and the
### observed contrast's recipe runs unchanged on every resample. The
### recolored matrices are built once per contrast before the workers fork.
###
### Inputs  : coexpr_perm_inputs.rda  (RESID, N.SC, N.SE, N.KEEP,
###                                    DRAWS.PERM.COEXPR, DRAWS.NULL.COEXPR,
###                                    PLOIDY.F)
### Output  : coexpr_perm_output.rda  (NULL.TOTAL.RANKS, NULL.CIS.RANKS,
###                                    NULL.TRANS.RANKS, NULL.PAR.SC.RANKS,
###                                    NULL.PAR.SE.RANKS: B x N.KEEP
###                                    matrices, column k = the draw's
###                                    k-th largest-magnitude squared
###                                    eigenvalue)
###############################################################

library('parallel')

source("functions.R")
write_pkg_versions("coexpr_perm")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

## Optional args: Rscript coexpr_perm.R [input.rda] [output.rda]
ARGS        <- commandArgs(trailingOnly = TRUE)
INPUT.FILE  <- if (length(ARGS) >= 1) ARGS[1] else "coexpr_perm_inputs.rda"
OUTPUT.FILE <- if (length(ARGS) >= 2) ARGS[2] else "coexpr_perm_output.rda"

load(INPUT.FILE)     # RESID, N.SC, N.SE, N.KEEP, DRAWS.PERM.COEXPR, DRAWS.NULL.COEXPR, PLOIDY.F
stopifnot(length(DRAWS.PERM.COEXPR) == length(DRAWS.NULL.COEXPR))

## Recolored cells of the two par nulls, one setup per contrast. The identity check prints the contrast
## the recipe gives on the full recolored matrices (near zero) next to the differential-shrinkage
## remainder it should equal. Forked workers inherit PAR.SETUP without copying it.
PAR.SETUP <- list(sc = coexpr_null_par_setup(RESID$MIX.SC, RESID$HYB.COMB, PLOIDY.F),
                   se = coexpr_null_par_setup(RESID$MIX.SE, RESID$HYB.COMB, PLOIDY.F))
for (nm in names(PAR.SETUP)) {
  cat(sprintf("par_%s identity check:\n", nm)); print(round(PAR.SETUP[[nm]]$check, 4))
}
flush.console()

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")

WORKER.PIDS <- unlist(clusterEvalQ(cl, Sys.getpid()))
cat(sprintf("requested %d cores, got %d worker PIDs: %s\n",
            NUM.CORES, length(WORKER.PIDS), paste(WORKER.PIDS, collapse = ", ")))
flush.console()

## ---- chunk the draws so the master can report progress ----
N.UPDATES  <- 40
B          <- length(DRAWS.PERM.COEXPR)
DRAWS      <- Map(list, perm = DRAWS.PERM.COEXPR, null = DRAWS.NULL.COEXPR)
CHUNK.SIZE <- max(NUM.CORES, ceiling(B / N.UPDATES))
chunks     <- split(seq_len(B), ceiling(seq_len(B) / CHUNK.SIZE))

cat(sprintf("coexpression null start: %d genes, B=%d, %d cores, %d chunks\n",
            nrow(RESID$MIX.SC), B, NUM.CORES, length(chunks)))
flush.console()

RESULTS <- vector("list", B)

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One draw of the five null spectra, for rank-matched testing of every candidate axis (observed rank k
  ## is compared with the null's own rank k). draw: list(perm, null) from DRAWS.PERM.COEXPR and
  ## DRAWS.NULL.COEXPR. resid: the RESID list (including HYB.COMB). nSC: Sc-parent cell count, used to
  ## split the pooled, reshuffled parent cells back into groups of the original sizes (total). par_setup:
  ## the recolored cells of the par nulls. n_keep: ranks retained per decomposition (15, the top-15
  ## candidate window). Returns the top n_keep squared eigenvalues by magnitude, descending, for all five
  ## decompositions; this is the unit of work a cluster worker does.
  chunk.results <- parLapply(cl, DRAWS[chunks[[k]]], coexpr_perm_one, n_keep = N.KEEP, nSC = N.SC, resid = RESID, par_setup = PAR.SETUP)
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
NULL.PAR.SC.RANKS <- do.call(rbind, lapply(RESULTS, `[[`, "par_sc"))
NULL.PAR.SE.RANKS <- do.call(rbind, lapply(RESULTS, `[[`, "par_se"))

save(NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS,
     NULL.PAR.SC.RANKS, NULL.PAR.SE.RANKS, file = OUTPUT.FILE)

cat(sprintf("coexpression null done: %d draws x %d ranks in %.1f min\n",
            B, N.KEEP,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
