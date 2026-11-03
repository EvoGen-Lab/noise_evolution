###############################################################
### null_calibration_grid.R
### Null calibration array for the validation checks D1 to D4, D7 and D10 (analysis/checks/checks_sim.R). Each task of
### calib_task_grid() builds one synthetic experiment of CAL.G genes from the shared generator on the observed design, with the
### null of its design true, and runs the pipeline's own fits and null generator (null_contrasts_one(), both estimators) on it.
###
### Inputs  : null_calibration_inputs.rda (CAL.GRID, CAL.G, CAL.N.SIM, CAL.N.PERM, CAL.SEED, CAL.STREAMS, POOLS), written by
###           checks_sim.R block CAL
### Output  : null_calibration_output_<k>of<K>.rda for array task k of K, holding CAL.TASKS (the grid rows this task ran) and
###           CAL.PART (one row per gene: the task's design, the gene's true values and every p-value). Tasks are dealt out
###           round-robin; each task draws from its own L'Ecuyer-CMRG stream, so the results do not depend on how the tasks are split.
###           checks_sim.R block CALMERGE reads the outputs and checks the task keys.
###############################################################

library('parallel')
library('MASS')

source("functions.R")
source("functions_checks.R")
write_pkg_versions("null_calibration_grid")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("null_calibration_inputs.rda")

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
ARRAY.ID  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE  <- sprintf("null_calibration_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the array only computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

CAL.TASKS <- which((seq_len(nrow(CAL.GRID)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
cat(sprintf("null calibration: array task %d of %d, %d grid tasks x %d genes, N.SIM %d, N.PERM %d, %d cores\n",
            ARRAY.ID, N.ARRAY, length(CAL.TASKS), CAL.G, CAL.N.SIM, CAL.N.PERM, NUM.CORES))
flush.console()

## The genes of a task are spread over the cores inside calib_dataset() (mclapply, one gene at a time), so the progress line
## after each task lands in the job log.
CAL.PART <- list(); t0 <- Sys.time()
for (j in seq_along(CAL.TASKS)) {
  CAL.PART[[j]] <- calib_dataset(CAL.TASKS[j], grid = CAL.GRID, G = CAL.G, n_sim = CAL.N.SIM, n_perm = CAL.N.PERM, pools = POOLS, streams = CAL.STREAMS, cores = NUM.CORES)
  cat(sprintf("[%s] task %d (%d of %d) done, %.1f min elapsed\n", format(Sys.time(), "%H:%M:%S"), CAL.TASKS[j], j, length(CAL.TASKS),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
CAL.PART <- do.call(rbind, CAL.PART)
stopifnot(setequal(unique(CAL.PART$task), CAL.TASKS))
save(CAL.TASKS, CAL.PART, file = OUT.FILE)
cat("null calibration done\n")
