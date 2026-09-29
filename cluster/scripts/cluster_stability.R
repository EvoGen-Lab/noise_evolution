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

source("Functions.R")
load("cluster_stability_inputs.rda")

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
TASKS <- CSTAB.INPUTS$tasks
B     <- CSTAB.INPUTS$key$B

## ---- Stage 1: bootstrap replicates ----
## b varies fastest, so replicates of one task stay contiguous and in order
JOBS <- expand.grid(b = seq_len(B), k = seq_len(nrow(TASKS)))

run_job <- function(j) {
  tk <- TASKS[JOBS$k[j], ]
  d  <- CSTAB.INPUTS$data[[tk$dataset]]
  tryCatch(boot_ari_one(d$counts, CSTAB.INPUTS$idx[[tk$dataset]][, JOBS$b[j]], CSTAB.INPUTS$ref[[tk$task]],
                        d$nfeatures, d$dims_n, tk$res, d$metric),
           error = function(e) { message(sprintf("%s replicate %d: %s", tk$task, JOBS$b[j], conditionMessage(e))); NA_real_ })
}

## Chunked so the master prints progress to the job log
chunks <- split(seq_len(nrow(JOBS)), ceiling(seq_len(nrow(JOBS)) / max(NUM.CORES, ceiling(nrow(JOBS) / 40))))
cat(sprintf("stability bootstrap start: %d tasks x B=%d = %d fits, %d cores\n", nrow(TASKS), B, nrow(JOBS), NUM.CORES)); flush.console()
ari <- numeric(0); t0 <- Sys.time()
for (ch in chunks) {
  r   <- mclapply(ch, run_job, mc.cores = NUM.CORES, mc.preschedule = FALSE)
  ari <- c(ari, vapply(r, function(x) if (is.numeric(x) && length(x) == 1) x else NA_real_, numeric(1)))   # a lost worker stays an NA slot
  el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("[%s] %d / %d fits  elapsed %.1f min  eta %.1f min\n", format(Sys.time(), "%H:%M:%S"),
              length(ari), nrow(JOBS), el, (nrow(JOBS) - length(ari)) * el / length(ari))); flush.console()
}
CSTAB.ARI <- split(ari, factor(TASKS$task[JOBS$k], levels = TASKS$task))
cat(sprintf("replicates completed: %d / %d\n", sum(is.finite(ari)), length(ari)))

## ---- Stage 2: marker enrichment at each dataset's final resolution ----
## Clusters within a dataset are spread over the forked workers
par_apply <- function(X, FUN) mclapply(X, FUN, mc.cores = NUM.CORES, mc.preschedule = FALSE)
CSTAB.MARKERS <- lapply(setNames(names(CSTAB.MARKER.OBJS), names(CSTAB.MARKER.OBJS)), function(d) {
  obj <- assemble_cluster_stability(CSTAB.INPUTS, CSTAB.ARI, d, CSTAB.MARKER.OBJS[[d]])$final_obj
  cat(sprintf("marker enrichment: %s, %d clusters\n", DS.LABELS[[d]], length(levels(Idents(obj))))); flush.console()
  cluster_marker_enrichment(obj, DS.LABELS[[d]], kegg_data = KEGG.DATA, apply_fun = par_apply)
})

CSTAB.KEY <- CSTAB.INPUTS$key
save(CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY, file = "cluster_stability_output.rda")
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
