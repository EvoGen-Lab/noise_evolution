###############################################################
### coexpr_boot.R
### Reads the co-expression bootstrap inputs, runs the B cell
### resamples on a FORK cluster, assembles and saves CB. Mirrors
### the structure of gene_boot.R, but chunks over draws rather
### than over genes, since the gene set here is fixed (CO.GENES)
### and the expensive step is one shrink_cor per dataset per draw.
###
### Inputs  : coexpr_boot1_inputs.rda, or coexpr_boot<tag>_inputs.rda for
###           a given tag arg (e.g. "2" for a second-seed adequacy check)
###           (RESID, COEXPR.POINT, DRAWS.COEXPR)
### Output  : coexpr_boot1_output.rda, or coexpr_boot<tag>_output.rda to
###           match (CB: CB$total, CB$cis, CB$trans, CB$lambda)
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

## One bootstrap draw of the co-expression decomposition. draw: one element of make_coexpr_draws().
## resid: the RESID list (genes x cells per dataset, including HYB.COMB). HYB.SC, HYB.SE and HYB.COMB
## are resampled with the same H draw because they are the same cells; the per-gene ploidy factors
## are fixed, so every draw rescales Rhyb identically. Returns total/cis/trans/dpar_sc/dpar_se at the
## upper-triangle pair positions; this is the unit of work a cluster worker does.
coexpr_bootstrap_one <- function(draw, resid) {
  d <- coexpr_decompose(list(
    MIX.SC = resid$MIX.SC[, draw$SC], MIX.SE = resid$MIX.SE[, draw$SE],
    HYB.SC = resid$HYB.SC[, draw$H],  HYB.SE = resid$HYB.SE[, draw$H],
    HYB.COMB = resid$HYB.COMB[, draw$H]),
    ploidy_f = attr(resid, "ploidy_f"))   # fixed per-gene factors, so each draw rescales Rhyb the same way
  p  <- nrow(resid$MIX.SC)
  up <- which(upper.tri(matrix(0, p, p)))
  list(total = d$total[up], cis = d$cis[up], trans = d$trans[up], dpar_sc = d$dpar_sc[up], dpar_se = d$dpar_se[up])
}

## The five decompositions the co-expression bootstrap tracks, in the order coexpr_bootstrap_one()
## returns them.
.COEXPR_PARTS <- c("total", "cis", "trans", "dpar_sc", "dpar_se")

## Streaming accumulator for the co-expression bootstrap SE. Holding every draw (a B x pairs matrix
## per contrast) costs O(B x pairs) memory and pairs grows as p^2. A standard deviation needs only
## sum(x) and sum(x^2) per pair, so folding in each chunk's draws keeps peak memory at
## O(chunk_size x pairs) and the raw draws can be discarded. p: number of genes (nrow of a RESID
## matrix), used to build the upper-triangle index once. Tracks the parts in .COEXPR_PARTS.
coexpr_acc_init <- function(p) {
  up <- which(upper.tri(matrix(0, p, p)))
  z  <- numeric(length(up))
  sums <- unlist(lapply(.COEXPR_PARTS, function(nm) setNames(list(z, z), paste0(c("sum_", "sumsq_"), nm))), recursive = FALSE)
  c(list(up = up, n = 0L), sums)
}

## Folds one chunk of coexpr_bootstrap_one() results (a list of per-draw lists) into the running
## sums and draw count; the caller can discard chunk_results afterwards.
coexpr_acc_update <- function(acc, chunk_results) {
  acc$n <- acc$n + length(chunk_results)
  for (nm in .COEXPR_PARTS) {
    m <- do.call(rbind, lapply(chunk_results, `[[`, nm))
    acc[[paste0("sum_", nm)]]   <- acc[[paste0("sum_", nm)]]   + colSums(m)
    acc[[paste0("sumsq_", nm)]] <- acc[[paste0("sumsq_", nm)]] + colSums(m^2)
  }
  acc
}

## Finishes the accumulator into the CB structure: for each part a data.frame of
## gene_i/gene_j/est/se/z/p, plus lambda. est is the point estimate pt, se the bootstrap SD from the
## running sums, z = est/se and p the two-sided normal tail.
coexpr_acc_finalize <- function(resid, pt, acc) {
  p  <- nrow(resid$MIX.SC); gn <- rownames(resid$MIX.SC)
  ij <- arrayInd(acc$up, c(p, p))
  base <- data.frame(gene_i = gn[ij[, 1]], gene_j = gn[ij[, 2]])

  se_of <- function(s, ss) sqrt(pmax(0, (ss - s^2 / acc$n) / (acc$n - 1)))
  mk <- function(est, s, ss) { se <- se_of(s, ss)
    data.frame(base, est = est, se = se, z = est / se, p = 2 * pnorm(-abs(est / se))) }

  parts <- lapply(setNames(.COEXPR_PARTS, .COEXPR_PARTS), function(nm)
    mk(pt[[nm]][acc$up], acc[[paste0("sum_", nm)]], acc[[paste0("sumsq_", nm)]]))
  c(parts, list(lambda = pt$lambda))
}

## Tag arg picks which input/output pair to use, so the same script
## serves both the primary run and any additional-seed adequacy check
## without duplicating the file. Defaults to "1", the primary run, the
## same convention gene_boot.R uses.
##   Rscript coexpr_boot.R    -> coexpr_boot1_inputs.rda / coexpr_boot1_output.rda
##   Rscript coexpr_boot.R 2  -> coexpr_boot2_inputs.rda / coexpr_boot2_output.rda
ARGS        <- commandArgs(trailingOnly = TRUE)
TAG         <- if (length(ARGS) >= 1) ARGS[1] else "1"
INPUT.FILE  <- sprintf("coexpr_boot%s_inputs.rda", TAG)
OUTPUT.FILE <- sprintf("coexpr_boot%s_output.rda", TAG)

load(INPUT.FILE)     # RESID, COEXPR.POINT, DRAWS.COEXPR

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

cat(sprintf("coexpression bootstrap start: %d genes, B=%d, %d cores, %d chunks\n",
            nrow(RESID$MIX.SC), B, NUM.CORES, length(chunks)))
flush.console()

ACC <- coexpr_acc_init(nrow(RESID$MIX.SC))

t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  chunk.results <- parLapply(cl, DRAWS.COEXPR[chunks[[k]]], coexpr_bootstrap_one,
                             resid = RESID)
  ACC <- coexpr_acc_update(ACC, chunk.results)
  rm(chunk.results); gc(FALSE)   # drop this chunk's raw draws before the next one

  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (B - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d draws (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, B, 100 * done / B, el, eta))
  flush.console()
}

CB <- coexpr_acc_finalize(RESID, COEXPR.POINT, ACC)
save(CB, file = OUTPUT.FILE)

cat(sprintf("coexpression bootstrap done: %d pairs in %.1f min [tag %s]\n",
            nrow(CB$total),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            TAG))

stopCluster(cl)
