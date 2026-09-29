###############################################################
### go_enrich.R
### Every Section 8 enrichment, run on forked workers:
### 8.4 over-representation for GO.SETS, 8.6 rank-based gseGO,
### 8.7 over-representation for the intrinsic-fraction sets.
###
### Inputs  : go_enrich_inputs.rda (GO.INPUTS, KEGG.DATA)
### Output  : go_enrich_output.rda
###           (GO.ENRICH, INTR.GO, INTR.GO.CLEAN, GO.GSE, GO.KEY)
###
### GO simplification and gseGO permutations dominate run time, and
### each set is independent, so each set is one job. KEGG runs offline
### from KEGG.DATA downloaded on the local side.
###############################################################

suppressPackageStartupMessages({
  library(parallel)
  library(BiocParallel)
  library(clusterProfiler)
  library(org.Sc.sgd.db)
})
register(SerialParam())   # one process per job, parallelism comes from mclapply

source("Functions.R")
load("go_enrich_inputs.rda")
NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

## ---- Job list: rank-based runs first, since they take longest ----
ora_jobs <- function(sets, universe, group)
  lapply(setNames(names(sets), paste0(group, "::", names(sets))),
         function(s) list(kind = "ora", genes = sets[[s]], universe = universe))
gse_jobs <- unlist(lapply(names(GO.INPUTS$GSE.LISTS), function(q)
  lapply(setNames(c("BP", "MF", "CC"), sprintf("GSE::GO.GSE.%s.%s", q, c("BP", "MF", "CC"))),
         function(ont) list(kind = "gse", ranks = GO.INPUTS$GSE.LISTS[[q]], ont = ont))), recursive = FALSE)
JOBS <- c(gse_jobs,
          ora_jobs(GO.INPUTS$GO.SETS,         GO.INPUTS$GO.UNIVERSE,         "GO"),
          ora_jobs(GO.INPUTS$INTR.SETS,       GO.INPUTS$INTR.UNIVERSE,       "INTR"),
          ora_jobs(GO.INPUTS$INTR.SETS.CLEAN, GO.INPUTS$INTR.UNIVERSE.CLEAN, "INTR.CLEAN"))

cat(sprintf("enrichment start: %d jobs (%d rank-based), %d cores\n", length(JOBS), length(gse_jobs), NUM.CORES)); flush.console()
t0 <- Sys.time()
RES <- mclapply(seq_along(JOBS), function(k)
  tryCatch(go_enrich_one(JOBS[[k]], KEGG.DATA, GO.INPUTS$GO.QVAL, GO.INPUTS$SEED.GO + k),
           error = function(e) { message(sprintf("%s: %s", names(JOBS)[k], conditionMessage(e))); NULL }),
  mc.cores = NUM.CORES, mc.preschedule = FALSE)
names(RES) <- names(JOBS)
cat(sprintf("enrichment done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- Unpack into the objects Section 8 reads ----
pick <- function(group) { r <- RES[startsWith(names(RES), paste0(group, "::"))]; setNames(r, sub("^[^:]+::", "", names(r))) }
GO.ENRICH     <- pick("GO")[names(GO.INPUTS$GO.SETS)]
INTR.GO       <- pick("INTR")[names(GO.INPUTS$INTR.SETS)]
INTR.GO.CLEAN <- pick("INTR.CLEAN")[names(GO.INPUTS$INTR.SETS.CLEAN)]
GO.GSE        <- pick("GSE")
GO.KEY        <- GO.INPUTS

save(GO.ENRICH, INTR.GO, INTR.GO.CLEAN, GO.GSE, GO.KEY, file = "go_enrich_output.rda")
