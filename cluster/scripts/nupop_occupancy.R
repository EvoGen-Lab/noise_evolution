###############################################################
### nupop_occupancy.R
### Reads the NuPoP inputs, predicts nucleosome occupancy for
### every chromosome that carries a promoter in both species, and
### saves each gene's promoter occupancy track. Mirrors the
### structure of gene_boot.R.
###
### Inputs  : nupop_inputs.rda (NUPOP.INPUTS.SC, NUPOP.INPUTS.SE)
### Output  : nupop_output.rda (NUPOP.OCC.SC, NUPOP.OCC.SE, one
###           occupancy vector per gene, with the promoter coords
###           and any reduced-flank or unscored regions attached)
###
### Each window runs in its own forked process, so the master keeps
### going and prints one progress line per round of windows.
### Requires the Bioconductor package NuPoP on the cluster:
###   BiocManager::install("NuPoP")
###############################################################

library('parallel')
library('NuPoP')

source("functions.R")
write_pkg_versions("nupop_occupancy")   # R and package versions of this job, compared locally by check_pkg_versions()

load("nupop_inputs.rda")     # NUPOP.INPUTS.SC, NUPOP.INPUTS.SE

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

t0 <- Sys.time()
cat(sprintf("NuPoP start: %d cores\n", NUM.CORES))

cat(sprintf("Sc: %d chromosomes, %.1f Mb\n", length(NUPOP.INPUTS.SC$chroms), sum(nchar(NUPOP.INPUTS.SC$chroms)) / 1e6))
NUPOP.OCC.SC <- nupop_occupancy_cluster(NUPOP.INPUTS.SC, cores = NUM.CORES)
report_regions(NUPOP.OCC.SC, "Sc")

cat(sprintf("Se: %d chromosomes, %.1f Mb\n", length(NUPOP.INPUTS.SE$chroms), sum(nchar(NUPOP.INPUTS.SE$chroms)) / 1e6))
NUPOP.OCC.SE <- nupop_occupancy_cluster(NUPOP.INPUTS.SE, cores = NUM.CORES)
report_regions(NUPOP.OCC.SE, "Se")

save(NUPOP.OCC.SC, NUPOP.OCC.SE, file = "nupop_output.rda")

cat(sprintf("NuPoP done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
