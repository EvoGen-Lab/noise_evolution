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
### ploidy_adjust_dpar() in functions.R).
###############################################################

library('parallel')
library('MASS')

source("functions.R")
load("gene_perm_inputs.rda")     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, SEED.PERM, PLOIDY.SHIFT (when present)

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
## Ploidy shift. When permute_contrasts_one() accepts a ploidy_shift argument, the
## dpar contrasts are read against the ploidy-expected value and PLOIDY.SHIFT must
## come with the inputs. The two are checked together so the shift can never be
## dropped silently. Forked workers inherit PLOIDY.SHIFT from this process.
HAS.PLOIDY.ARG <- "ploidy_shift" %in% names(formals(permute_contrasts_one))
HAS.PLOIDY.VAL <- exists("PLOIDY.SHIFT")
if (HAS.PLOIDY.ARG != HAS.PLOIDY.VAL)
  stop(sprintf("ploidy shift mismatch: permute_contrasts_one() takes ploidy_shift = %s, PLOIDY.SHIFT in inputs = %s",
               HAS.PLOIDY.ARG, HAS.PLOIDY.VAL))
EXTRA.ARGS <- if (HAS.PLOIDY.ARG) list(ploidy_shift = PLOIDY.SHIFT) else list()
cat(sprintf("ploidy shift active: %s\n", HAS.PLOIDY.ARG))
if (!HAS.PLOIDY.ARG) cat("note: this functions.R has no ploidy_shift, so PERM.RESULTS will hold no _p_ploidy or _p_ind columns\n")
flush.console()

## Forked workers inherit PERMS, the fits, and the loaded functions from this
## process without copying them, so one shared PERMS object serves every core.
## mc.preschedule = FALSE hands genes to workers one at a time, which keeps all
## cores busy when genes differ in cost.

N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("permutation start: %d genes, NPERM=%d, %d cores, %d chunks\n",
            length(GENES), N.PERM, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- do.call(mclapply, c(list(X = chunks[[k]], FUN = permute_contrasts_one,
                           mats = CONTRAST.MATS, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, perms = PERMS,
                           mc.cores = NUM.CORES, mc.preschedule = FALSE), EXTRA.ARGS))
  bad <- vapply(results[[k]], function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop(sprintf("%d gene(s) failed in chunk %d, first: %s", sum(bad), k, chunks[[k]][which(bad)[1]]))
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

