###############################################################
### gene_perm.R
### Reads the same inputs as the bootstrap, runs the paired
### per-gene permutation null on a FORK cluster, saves the output.
### Mirrors gene_boot.R.
###
### Inputs  : gene_perm_inputs.rda  (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, NULL.DRAWS,
###           SHARE.OBS: the observed transcriptome share quantities (mean share, share corrections applied to the cis nulls, per-gene ploidy shifts kept by
###           every null), GATES: the contrasts each gene qualifies for, ...)
### Output  : gene_perm_output.rda
###           PERM.RESULTS, one row per gene: mean_<mode>_obs / _p,
###           disp_<mode>_obs / _p (the NB dispersion phi), bsize_<mode>_obs / _p
###           (mean plus disp, its own null), kbal_<mode>_obs / _p
###           (minus bsize minus disp, its own null), and cv2_<mode>_obs / _p
###           (1/mean + phi, its own null).
###           For par_sc and par_se, the noise quantities also get
###           _p_ploidy and _p_ind columns, testing the contrast against
###           the ploidy-expected value using this gene's own shift and
###           using the independent-allele shift of 1, plus the
###           ploidy_shift column recording the shift used.
###
### Each contrast is nulled by its own label shuffle. total, cis, and
### trans are zero-centred under the null they test. par tests
### exchangeability of hybrid and parent cells, which is the null the
### ploidy-shifted observed value is read against (ploidy_shift() and
### the ploidy adjustment of the par contrasts() in functions.R).
###############################################################

library('parallel')
library('MASS')

source("functions.R")
write_pkg_versions("gene_perm")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("gene_perm_inputs.rda")     # GENES is the inference set (every gene for which at least one comparison is eligible), not only the matched set. CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, SEED.PERM, SHARE.OBS, GATES, NULL.DRAWS

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
## The mid and par contrasts are read against the ploidy-expected value, so SHARE.OBS (with the per-gene shifts) must come with the
## inputs. Forked workers inherit it from this process.
stopifnot("SHARE.OBS is missing from gene_perm_inputs.rda" = exists("SHARE.OBS"), "shift" %in% names(SHARE.OBS))

## Forked workers inherit PERMS, the fits, and the loaded functions from this process without copying
## them, so one shared PERMS object serves every core. mc.preschedule = FALSE hands genes to workers
## one at a time, which keeps all cores busy when genes differ in cost.

N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("permutation start: %d genes, NPERM=%d, %d cores, %d chunks\n",
            length(GENES), N.PERM, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's permutation null: refits the relabeled groups for every mode and permutation, compares
  ## the observed contrasts (mean, disp, CV2, bsize, kbal) with the null by a two-sided permutation p, and
  ## adds ploidy-adjusted p-values for the par contrasts. bsize and kbal nulls recombine the mean and
  ## disp null draws of the same permutation; for cis and trans those draws come from relabelings of the
  ## same cells (perm_label_draw()).
  results[[k]] <- mclapply(chunks[[k]], permute_contrasts_one, mc.cores = NUM.CORES, mc.preschedule = FALSE, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS, perms = PERMS, share_obs = SHARE.OBS, seed = SEED.PERM, gates = GATES, resamples = NULL.DRAWS, methods = c("NB", "weighted"))
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

