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

## boot_ari_one: one bootstrap replicate. Resamples the columns of counts
## by idx, applies the same preparation as to_seurat_counts() in analysis.R
## (underscore to dash in gene names, CsparseMatrix) once per replicate,
## reruns the Normalize/HVG/Scale/PCA/Neighbors/Clusters pipeline at the
## original settings, and returns the adjusted Rand index against
## ref_clusters on the same resampled cells. Self-contained, so it runs
## equally well serially or on a forked worker (cluster_stability.R).
boot_ari_one <- function(counts, idx, ref_clusters, nfeatures, dims_n, resolution, metric = "manhattan") {
  boot_counts <- counts[, idx]
  rownames(boot_counts) <- gsub("_", "-", rownames(counts), fixed = TRUE)
  colnames(boot_counts) <- make.unique(colnames(counts)[idx])
  boot_counts <- as(boot_counts, "CsparseMatrix")
  boot_obj <- CreateSeuratObject(counts = boot_counts)
  boot_obj <- NormalizeData(boot_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  boot_obj <- FindVariableFeatures(boot_obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  boot_obj <- ScaleData(boot_obj, features = rownames(boot_obj), verbose = FALSE)
  boot_obj <- RunPCA(boot_obj, features = VariableFeatures(boot_obj), verbose = FALSE)
  boot_obj <- FindNeighbors(boot_obj, reduction = "pca", dims = 1:dims_n, annoy.metric = metric, verbose = FALSE)
  boot_obj <- FindClusters(boot_obj, resolution = resolution, verbose = FALSE)
  adjustedRandIndex(as.integer(ref_clusters[idx]), as.integer(Idents(boot_obj)))
}
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
## For each cluster: marker genes against the rest of the SAME dataset's cells (FindMarkers with
## ident.2 left at its default) and GO/KEGG over-representation of the up and down markers
## (run_enrichment), for however many clusters the dataset has. The object carries its own final,
## validated clustering in Idents, so the question answered is whether the clustering validated for
## THIS dataset corresponds to distinguishable biology; no GSEA is computed. A dataset with fewer than
## two clusters gives NULL (with a message); otherwise the entry is a list of
## markers/up/down/up_enrich/down_enrich/background/cluster_ids, one element per cluster, keyed by cluster ID.
CSTAB.MARKERS <- lapply(setNames(names(CSTAB.MARKER.OBJS), names(CSTAB.MARKER.OBJS)), function(d) {
  obj <- assemble_cluster_stability(CSTAB.INPUTS, CSTAB.ARI, d, CSTAB.MARKER.OBJS[[d]])$final_obj
  cat(sprintf("marker enrichment: %s, %d clusters\n", DS.LABELS[[d]], length(levels(Idents(obj))))); flush.console()
  cluster_ids <- sort(unique(as.character(Idents(obj))))
  if (length(cluster_ids) < 2) {
    cat(sprintf("%s: only one cluster found; skipping the per-cluster marker/enrichment comparison.\n", DS.LABELS[[d]]))
    return(NULL)
  }

  markers_list <- par_apply(cluster_ids, function(cc) {
    m <- suppressWarnings(FindMarkers(obj, ident.1 = cc))
    m[order(m$avg_log2FC, decreasing = TRUE), ]
  })
  names(markers_list) <- cluster_ids

  up_list   <- lapply(markers_list, function(m) m[abs(m$avg_log2FC) > log2(1.25) & -log10(m$p_val_adj) > 20 & m$avg_log2FC > 0, ])
  down_list <- lapply(markers_list, function(m) m[abs(m$avg_log2FC) > log2(1.25) & -log10(m$p_val_adj) > 20 & m$avg_log2FC < 0, ])

  ## avg_log2FC is selected by name: FindMarkers column order differs across Seurat versions.
  gene_vec <- function(m) { v <- m[["avg_log2FC"]]; names(v) <- row.names(m); v }
  background_list <- lapply(markers_list, gene_vec)
  up_genes_list    <- lapply(up_list,   gene_vec)
  down_genes_list  <- lapply(down_list, gene_vec)

  enrich <- setNames(par_apply(cluster_ids, function(cc) list(
    up   = run_enrichment(names(up_genes_list[[cc]]),   names(background_list[[cc]]), kegg_data = KEGG.DATA),
    down = run_enrichment(names(down_genes_list[[cc]]), names(background_list[[cc]]), kegg_data = KEGG.DATA))), cluster_ids)
  up_enrich   <- lapply(enrich, `[[`, "up")
  down_enrich <- lapply(enrich, `[[`, "down")

  list(markers = markers_list, up = up_list, down = down_list,
       up_enrich = up_enrich, down_enrich = down_enrich,
       background = background_list, cluster_ids = cluster_ids)
})

CSTAB.KEY <- CSTAB.INPUTS$key
save(CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY, file = "cluster_stability_output.rda")
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
