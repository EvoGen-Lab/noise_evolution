###############################################################
### coexpr_perm.R
### Reads the co-expression permutation-null inputs, runs the joint
### total/cis/dpar_sc/dpar_se draws on a FORK cluster, and saves the
### top-N.KEEP squared-eigenvalue spectrum of total, cis, trans, dpar_sc,
### and dpar_se for every draw. Used for rank-matched testing (real rank
### k against the null's own rank k) of every candidate axis in all five
### decompositions, replacing the earlier single shared-threshold
### approach. Mirrors the structure of coexpr_boot.R, but each
### draw's result is a handful of numbers (5 x N.KEEP) rather than a
### full pairwise matrix, so no incremental accumulator is needed;
### results are collected directly and stacked into matrices at the end.
###
### Inputs  : coexpr_perm_inputs.rda  (RESID, N.SC, N.SE, N.KEEP,
###                                    DRAWS.PERM.COEXPR)
### Output  : coexpr_perm_output.rda  (NULL.TOTAL.RANKS, NULL.CIS.RANKS,
###                                    NULL.TRANS.RANKS, NULL.DPAR.SC.RANKS,
###                                    NULL.DPAR.SE.RANKS: B x N.KEEP
###                                    matrices, column k = the draw's
###                                    k-th largest-magnitude squared
###                                    eigenvalue)
###############################################################

library('parallel')

source("functions.R")

## coexpr_perm_one: one permutation draw of the five null spectra, for rank-matched testing of every candidate axis. draw: one element of
## make_coexpr_perm_draws(). resid: the RESID list (including HYB.COMB). nSC, nSE: parent cell counts, used to split each pooled, reshuffled
## pool back into groups of the original sizes. n_keep: ranks retained per decomposition. Returns the top n_keep squared eigenvalues by
## magnitude, descending, for all five decompositions.
coexpr_perm_one <- function(draw, n_keep, nSC, nSE, resid) {
  pooled   <- cbind(resid$MIX.SC, resid$MIX.SE)
  perm.sc  <- pooled[, draw$idx[seq_len(nSC)]]
  perm.se  <- pooled[, draw$idx[-seq_len(nSC)]]
  perm.hsc <- resid$HYB.SC; perm.hse <- resid$HYB.SE
  perm.hsc[, draw$swap] <- resid$HYB.SE[, draw$swap]
  perm.hse[, draw$swap] <- resid$HYB.SC[, draw$swap]

  d <- coexpr_decompose(list(MIX.SC = perm.sc, MIX.SE = perm.se, HYB.SC = perm.hsc, HYB.SE = perm.hse))
  topk <- function(m) {
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    (ev[order(abs(ev), decreasing = TRUE)][seq_len(n_keep)])^2
  }

  ## The dpar nulls are zero-centred exchangeability nulls on the raw pooled cells (no ploidy
  ## rescale); the observed dpar matrices (COEXPR.POINT) carry the rescale from coexpr_decompose().
  ## dpar_sc's null pools Sc-parent and allele-summed hybrid cells and reshuffles them into
  ## pseudo-Sc / pseudo-hybrid groups of the original sizes, giving the leading eigenvalues of
  ## "condition A minus condition B" when condition carries no information. dpar_se's null does the
  ## same with Se-parent and hybrid cells.
  pooled.sc <- cbind(resid$MIX.SC, resid$HYB.COMB)
  perm.a.sc <- pooled.sc[, draw$idx_dpar_sc[seq_len(nSC)]]
  perm.b.sc <- pooled.sc[, draw$idx_dpar_sc[-seq_len(nSC)]]
  dpar_sc_null <- topk(shrink_cor(t(perm.b.sc)) - shrink_cor(t(perm.a.sc)))

  pooled.se <- cbind(resid$MIX.SE, resid$HYB.COMB)
  perm.a.se <- pooled.se[, draw$idx_dpar_se[seq_len(nSE)]]
  perm.b.se <- pooled.se[, draw$idx_dpar_se[-seq_len(nSE)]]
  dpar_se_null <- topk(shrink_cor(t(perm.b.se)) - shrink_cor(t(perm.a.se)))

  list(total = topk(d$total), cis = topk(d$cis), trans = topk(d$trans), dpar_sc = dpar_sc_null, dpar_se = dpar_se_null)
}

## Optional args: Rscript coexpr_perm.R [input.rda] [output.rda]
ARGS        <- commandArgs(trailingOnly = TRUE)
INPUT.FILE  <- if (length(ARGS) >= 1) ARGS[1] else "coexpr_perm_inputs.rda"
OUTPUT.FILE <- if (length(ARGS) >= 2) ARGS[2] else "coexpr_perm_output.rda"

load(INPUT.FILE)     # RESID, N.SC, N.SE, N.KEEP, DRAWS.PERM.COEXPR

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")

WORKER.PIDS <- unlist(clusterEvalQ(cl, Sys.getpid()))
cat(sprintf("requested %d cores, got %d worker PIDs: %s\n",
            NUM.CORES, length(WORKER.PIDS), paste(WORKER.PIDS, collapse = ", ")))
flush.console()

## ---- chunk the draws so the master can report progress ----
N.UPDATES  <- 40
B          <- length(DRAWS.PERM.COEXPR)
CHUNK.SIZE <- max(NUM.CORES, ceiling(B / N.UPDATES))
chunks     <- split(seq_len(B), ceiling(seq_len(B) / CHUNK.SIZE))

cat(sprintf("coexpression permutation null start: %d genes, B=%d, %d cores, %d chunks\n",
            nrow(RESID$MIX.SC), B, NUM.CORES, length(chunks)))
flush.console()

RESULTS <- vector("list", B)

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One permutation draw of the five null spectra, for rank-matched testing of every candidate axis
  ## (observed rank k is compared with the null's own rank k). draw: one element of
  ## make_coexpr_perm_draws(). resid: the RESID list (including HYB.COMB). nSC, nSE: parent cell counts,
  ## used to split each pooled, reshuffled pool back into groups of the original sizes. n_keep: ranks
  ## retained per decomposition (15, the top-15 candidate window). Returns the top n_keep squared
  ## eigenvalues by magnitude, descending, for all five decompositions; this is the unit of work a
  ## cluster worker does.
  chunk.results <- parLapply(cl, DRAWS.PERM.COEXPR[chunks[[k]]], coexpr_perm_one, n_keep = N.KEEP, nSC = N.SC, nSE = N.SE, resid = RESID)
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
NULL.DPAR.SC.RANKS <- do.call(rbind, lapply(RESULTS, `[[`, "dpar_sc"))
NULL.DPAR.SE.RANKS <- do.call(rbind, lapply(RESULTS, `[[`, "dpar_se"))

save(NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS,
     NULL.DPAR.SC.RANKS, NULL.DPAR.SE.RANKS, file = OUTPUT.FILE)

cat(sprintf("coexpression permutation null done: %d draws x %d ranks in %.1f min\n",
            B, N.KEEP,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
