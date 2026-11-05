###############################################################
### artefact_null_grid.R
### The simulation arm of the tests of why the data depart from the NB (analysis/checks/artefact_tests.R): T0, the gene-matched NB null of the
### trimming, depth-stratum and thinning statistics with the per-cell burden of T1 and the nulls of T4 and T11, over chunks of genes; and T8 and T12,
### simulated datasets under the null-divergence scenarios A and B and the artefact scenarios S1 to S5, each fitted and summarized like the data.
###
### Inputs  : artefact_inputs.rda (ART.MATS, ART.EXPOS, ART.ETA, ART.TERC, ART.FITS, ART.HALVES, ART.GENES, ART.CFG, ART.TASKS, ART.CHUNKS,
###           ART.SCEN.GRID, ART.STREAMS.GENE, ART.STREAMS.SCEN), written by artefact_tests.R
### Output  : artefact_output_<k>of<K>.rda for array task k of K, holding ART.PART: the results of the tasks of ART.TASKS this array task ran
###           (round-robin), keyed by task id. Every gene and every scenario draws from its own L'Ecuyer-CMRG stream, so the results do not
###           depend on how the tasks are split. artefact_tests.R reads the outputs and checks the task keys.
###############################################################

library('parallel')
library('MASS')

source("functions.R")
source("functions_checks.R")
source("functions_artefacts.R")
write_pkg_versions("artefact_null_grid")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("artefact_inputs.rda")

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
ARRAY.ID  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE  <- sprintf("artefact_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the array only computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

MY.TASKS <- which((seq_len(nrow(ART.TASKS)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
cat(sprintf("artefact grid: array task %d of %d, %d tasks, %d cores\n", ARRAY.ID, N.ARRAY, length(MY.TASKS), NUM.CORES)); flush.console()

## A T0 task splits its chunk of genes over the cores (one sub-chunk per core, each returning summed burden matrices); a scenario task spreads its genes
## over the cores inside art_scenario_task().
ART.PART <- list(); t0 <- Sys.time()
for (i in MY.TASKS) {
  tk <- ART.TASKS[i, ]
  ART.PART[[as.character(i)]] <- if (tk$kind == "T0") {
    sub <- split(ART.CHUNKS[[tk$index]], (seq_along(ART.CHUNKS[[tk$index]]) - 1) %% NUM.CORES)     # one sub-chunk per core, round-robin
    mclapply(sub, art_null_chunk, mats = ART.MATS, expos = ART.EXPOS, eta = ART.ETA, frac = ART.CFG$frac, terc = ART.TERC, R = ART.CFG$R, R4 = ART.CFG$R4,
             R11 = ART.CFG$R11, thin_p = ART.CFG$thin_p, streams = ART.STREAMS.GENE, mc.cores = NUM.CORES, mc.preschedule = FALSE)
  } else {
    art_scenario_task(tk$index, grid = ART.SCEN.GRID, fits = ART.FITS, expos = ART.EXPOS, eta = ART.ETA, halves = ART.HALVES, genes = ART.GENES, frac = ART.CFG$frac,
                      cbar = ART.CFG$cbar, mbar_breaks = ART.CFG$mbar_breaks, streams = ART.STREAMS.SCEN, cores = NUM.CORES)
  }
  cat(sprintf("[%s] task %d (%s %d) done, %.1f min elapsed\n", format(Sys.time(), "%H:%M:%S"), i, tk$kind, tk$index, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
save(ART.PART, file = OUT.FILE)
cat("artefact grid done\n")
