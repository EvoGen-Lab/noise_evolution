###############################################################
### gene_split.R
### Reads the count matrices for all 13 datasets (5 split-
### independent, 4 at the f_mean split, 4 at the f_disp split),
### runs the per-gene NB point-estimate fit for each dataset on a
### FORK cluster, saves the output. These are the observed fits
### used for GENES filtering and for the point estimates in BURST.CONTRASTS/PR,
### not a bootstrap or permutation; each gene is fit once per
### dataset, just spread across cores instead of run one at a time
### locally.
###
### Inputs  : gene_split_inputs.rda  (SPLIT.FIT.MATS, SPLIT.FIT.EXPOS,
###           two named lists with matching names, one entry per
###           dataset)
### Output  : gene_split_output.rda  (FIT.MIX.SC, FIT.MIX.SE,
###           FIT.HYB.SC, FIT.HYB.SE, FIT.HYB.COMB, FIT.HYC.SC,
###           FIT.HYC.SE, FIT.HYT.SC, FIT.HYT.SE, FIT.HYC.SC.N,
###           FIT.HYC.SE.N, FIT.HYT.SC.N, FIT.HYT.SE.N)
###############################################################

library('parallel')
library('MASS')

source("Functions.R")
load("gene_split_inputs.rda")  # SPLIT.FIT.MATS, SPLIT.FIT.EXPOS

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

FIT.RESULTS <- vector("list", length(SPLIT.FIT.MATS))
names(FIT.RESULTS) <- names(SPLIT.FIT.MATS)

t0 <- Sys.time()
for (nm in names(SPLIT.FIT.MATS)) {
  cat(sprintf("[%s] fitting %-10s %d genes, %d cells\n",
              format(Sys.time(), "%H:%M:%S"), nm,
              nrow(SPLIT.FIT.MATS[[nm]]), ncol(SPLIT.FIT.MATS[[nm]])))
  flush.console()
  FIT.RESULTS[[nm]] <- fit_counts_offset_parallel(SPLIT.FIT.MATS[[nm]], SPLIT.FIT.EXPOS[[nm]], cl)
}

for (nm in names(FIT.RESULTS)) assign(paste0("FIT.", nm), FIT.RESULTS[[nm]])
save(list = paste0("FIT.", names(FIT.RESULTS)), file = "gene_split_output.rda")

cat(sprintf("split fits done: %d datasets in %.1f min\n",
            length(FIT.RESULTS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
