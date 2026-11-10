###############################################################
### gene_perm.R
### The null generator: reads the inputs, builds every gene's null distributions on a FORK cluster
### (null_contrasts_one() in functions.R), and saves the p-values.
###
### Inputs  : gene_perm_inputs.rda  (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, N.SIM, SEED.PERM,
###           NULL.DRAWS: the pre-drawn resamples of the shifted bootstrap nulls (make_draws()), SHARE.CELLS: each hybrid cell's observed
###           Sc share and library size, from which every null draw refits the share terms, SHARE.OBS: the observed transcriptome share
###           quantities (mean share, share corrections of the datasets, each gene's extrinsic variance sigma_d for the mid mean simulation),
###           GATES: the contrasts each gene's fits define)
### Output  : gene_perm_output.rda
###           PERM.RESULTS, one row per gene: <quantity>_<mode>_obs / _p for mean, disp, rdisp, bsize, kbal and their weighted
###           versions, the share of infinite null draws (_infrac) and k_trans, the level factor of the trans mean null. cv2 and the raw
###           burst quantities carry the observed value only. analysis.R keeps a p-value only when every group of the contrast is resolved.
###
### The log2 mean contrasts are tested by permutation (total: pooled parents; cis: swapped rescaled allele counts; trans: each parent
### pooled with the same-genome allele of the trans cells, with the level factor k in the allele offsets; par: pooled rescaled hybrid
### sums and one parent) and by simulation for mid. The dispersion, burst size and kinetic balance contrasts are tested against shifted
### bootstrap nulls (N.SIM draws): the cells are resampled, every group is refitted with its share terms refit on the resampled hybrid cells,
### and each refitted dispersion is moved by the ratio of its null value to its fitted value.
###############################################################

library('parallel')
library('MASS')

source("functions.R")
write_pkg_versions("gene_perm")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("gene_perm_inputs.rda")     # GENES is the inference set (every gene for which at least one contrast is defined). CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, N.SIM, SEED.PERM, NULL.DRAWS, SHARE.CELLS, SHARE.OBS, GATES

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
## The mid mean simulation needs each gene's extrinsic variance (SHARE.OBS$sigma_d) and the shifted bootstrap nulls the resamples and the per-cell share data.
## Forked workers inherit them from this process.
stopifnot("SHARE.OBS is missing from gene_perm_inputs.rda" = exists("SHARE.OBS"), "sigma_d" %in% names(SHARE.OBS),
          "NULL.DRAWS is missing from gene_perm_inputs.rda" = exists("NULL.DRAWS"), "SHARE.CELLS is missing from gene_perm_inputs.rda" = exists("SHARE.CELLS"))

## Forked workers inherit PERMS, NULL.DRAWS, the fits, and the loaded functions from this process without copying
## them, so one shared copy serves every core. mc.preschedule = FALSE hands genes to workers
## one at a time, which keeps all cores busy when genes differ in cost.

N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("null start: %d genes, NPERM=%d, NSIM=%d, %d cores, %d chunks\n",
            length(GENES), N.PERM, N.SIM, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's null distributions: permutations for the log2 mean contrasts, a simulation for the mid mean, shifted bootstrap nulls for every
  ## dispersion and kinetic contrast, each draw refitted with the observed estimators and share corrections, and two-sided p-values.
  results[[k]] <- mclapply(chunks[[k]], null_contrasts_one, mc.cores = NUM.CORES, mc.preschedule = FALSE, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS,
                           perms = PERMS, draws = NULL.DRAWS, share_cells = SHARE.CELLS, share_obs = SHARE.OBS, n_sim = N.SIM, seed = SEED.PERM, gates = GATES, methods = c("NB", "weighted"))
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

cat(sprintf("null generation done: %d genes in %.1f min\n",
            nrow(PERM.RESULTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

