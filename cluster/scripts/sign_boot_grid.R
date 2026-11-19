###############################################################
### sign_boot_grid.R
### Cluster array of check P1 (analysis/checks/sign_boot_test.R): the bootstrap sign test of the dispersion-family contrasts on every cell of P1.GRID.
###
### Inputs  : sign_boot_inputs.rda (P1.GRID, P1.B, P1.POOLS, P1.STREAMS), written by sign_boot_test.R block INPUTS
### Output  : sign_boot_output_<k>of<K>.rda for array task k of K, holding P1.TASKS (the grid cells this task ran) and P1.PART (one row per cell and
###           statistic, sign_boot_cell()). Cells are dealt round-robin and spread over the cores; each cell draws from its own L'Ecuyer-CMRG stream,
###           so the results do not depend on how the cells are split. sign_boot_test.R block MERGE reads the outputs and checks that they cover the grid.
###############################################################

library('parallel')

source("functions.R")
source("functions_checks.R")
write_pkg_versions("sign_boot_grid")   # R and package versions of this job, compared with renv.lock by check_pkg_versions()

load("sign_boot_inputs.rda")

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
ARRAY.ID  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
OUT.FILE  <- sprintf("sign_boot_output_%dof%d.rda", ARRAY.ID, N.ARRAY)
## A finished task leaves its output file, so resubmitting the array only computes the tasks that are still missing.
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

P1.TASKS <- which((seq_len(nrow(P1.GRID)) - 1L) %% N.ARRAY == ARRAY.ID - 1L)
cat(sprintf("sign boot test: array task %d of %d, %d cells, B = %d, %d cores\n", ARRAY.ID, N.ARRAY, length(P1.TASKS), P1.B, NUM.CORES))
flush.console()

t0 <- Sys.time()
P1.PART <- do.call(rbind, mclapply(P1.TASKS, sign_boot_cell, grid = P1.GRID, pools = P1.POOLS, B = P1.B, streams = P1.STREAMS,
                                   mc.cores = NUM.CORES, mc.preschedule = FALSE))
stopifnot(setequal(unique(P1.PART$task), P1.TASKS))
save(P1.TASKS, P1.PART, file = OUT.FILE)
cat(sprintf("sign boot test done: %d cells in %.1f min\n", length(P1.TASKS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
