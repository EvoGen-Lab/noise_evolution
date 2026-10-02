###############################################################
### cluster_stability.R
### Section 7.3 cluster-stability bootstrap and Section 7.4 marker
### enrichment, run on forked workers. Mirrors gene_boot.R.
###
### Inputs  : cluster_stability_inputs.rda
###           (CSTAB.INPUTS, CSTAB.MARKER.OBJS, DS.LABELS, KEGG.DATA)
### Output  : cluster_stability_output.rda
###           (CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY)
###
### Every (dataset, resolution, replicate) is one job, so all cores
### stay busy across datasets of very different size. Resamples were
### drawn locally, which makes the result independent of core count
### and job order.
###############################################################

suppressPackageStartupMessages({
  library(parallel)
  library(Matrix)
  library(Seurat)
  library(mclust)          # adjustedRandIndex()
  library(clusterProfiler)
  library(org.Sc.sgd.db)
})

source("functions.R")
write_pkg_versions("cluster_stability")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("cluster_stability_inputs.rda")

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
TASKS <- CSTAB.INPUTS$tasks
B     <- CSTAB.INPUTS$key$B

## ---- Stage 1: bootstrap replicates ----
## b varies fastest, so replicates of one task stay contiguous and in order
JOBS <- expand.grid(b = seq_len(B), k = seq_len(nrow(TASKS)))

## Chunked so the master prints progress to the job log
chunks <- split(seq_len(nrow(JOBS)), ceiling(seq_len(nrow(JOBS)) / max(NUM.CORES, ceiling(nrow(JOBS) / 40))))
cat(sprintf("stability bootstrap start: %d tasks x B=%d = %d fits, %d cores\n", nrow(TASKS), B, nrow(JOBS), NUM.CORES)); flush.console()
ari <- numeric(0); t0 <- Sys.time()
for (ch in chunks) {
  r   <- mclapply(ch, bootstrap_ari_job, mc.cores = NUM.CORES, mc.preschedule = FALSE, inputs = CSTAB.INPUTS, jobs = JOBS, tasks = TASKS)
  ari <- c(ari, vapply(r, function(x) if (is.numeric(x) && length(x) == 1) x else NA_real_, numeric(1)))   # a lost worker stays an NA slot
  el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] %d / %d fits  elapsed %.1f min  eta %.1f min\n", format(Sys.time(), "%H:%M:%S"),
              length(ari), nrow(JOBS), el, (nrow(JOBS) - length(ari)) * el / length(ari))); flush.console()
}
CSTAB.ARI <- split(ari, factor(TASKS$task[JOBS$k], levels = TASKS$task))
cat(sprintf("replicates completed: %d / %d\n", sum(is.finite(ari)), length(ari)))

## ---- Stage 2: marker enrichment at each dataset's final resolution ----

## Marker genes and enrichment for every cluster of each dataset, keyed by cluster ID (see dataset_marker_enrichment()).
CSTAB.MARKERS <- lapply(setNames(names(CSTAB.MARKER.OBJS), names(CSTAB.MARKER.OBJS)), dataset_marker_enrichment, ari = CSTAB.ARI, inputs = CSTAB.INPUTS, marker_objs = CSTAB.MARKER.OBJS, labels = DS.LABELS, kegg_data = KEGG.DATA, cores = NUM.CORES)

CSTAB.KEY <- CSTAB.INPUTS$key
save(CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY, file = "cluster_stability_output.rda")
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
