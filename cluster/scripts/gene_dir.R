###############################################################
### gene_dir.R
### The direction bootstrap of the classification by direction (analysis.R Section 2.4; spec "Classification by direction"), on a FORK cluster, one job
### array over gene blocks per stage, set by GENE_DIR_STAGE:
###   stage 1  every gene, 1,000 replicates (DIR.PLANS$stage1)
###   stage 2  10,000 fresh replicates for the genes whose family has a margin test with n0 < 40 at stage 1 (DIR.PLANS$stage2)
###   stage 3  100,000 fresh replicates for the genes whose family has a margin test with n0 < 50 at stage 2 (DIR.PLANS$stage3)
###   stage 4  the theta intervals of the divergent genes at level 1 - 2 alpha_FCR, from the 100,000 draws of the stage-3 plan (THETA.TAIL)
### Every replicate resamples each group, repartitions the resampled hybrid cells at random into cis, trans, A and B sets and refits (dir_draw(),
### dir_contrasts_one()); the mean family refits means alone.
### Inputs  : gene_dir_inputs.rda (DIR.MATS, DIR.EXPOS, DIR.ETAS, DIR.FITS: the parents and the two alleles on all cells; GENES; DIR.PLANS; DIR.DELTA; DOM.SHIFT,
###           NULL when the shared shift is not applied) and gene_dir_stage<s>_inputs.rda (DIR.JOBS: the genes of each family, mean and disp; DIR.PLAN: the
###           plan's name; DIR.KEY; N.DIR.ARRAY; THETA.TAIL, NULL but in stage 4)
### Output  : gene_dir_stage<s>_output_<k>of<K>.rda holding DIR.OUT (one row per gene, quantity and scheme) and DIR.KEY. A finished task leaves its file, so
###           resubmitting the array computes only the missing blocks.
###############################################################

library('parallel')
library('MASS')

source("functions.R")
STAGE <- as.integer(Sys.getenv("GENE_DIR_STAGE", unset = "1"))
stopifnot(STAGE %in% 1:4)
write_pkg_versions(sprintf("gene_dir_stage%d", STAGE))   # compared with renv.lock by check_pkg_versions()

load("gene_dir_inputs.rda")                     # DIR.MATS, DIR.EXPOS, DIR.ETAS, DIR.FITS, GENES, DIR.PLANS, DIR.DELTA, DOM.SHIFT
load(sprintf("gene_dir_stage%d_inputs.rda", STAGE))   # DIR.JOBS, DIR.PLAN, DIR.KEY, N.DIR.ARRAY, THETA.TAIL
stopifnot(is.list(DIR.JOBS), all(names(DIR.JOBS) %in% c("mean", "disp")), all(unlist(DIR.JOBS) %in% GENES), DIR.PLAN %in% names(DIR.PLANS))

ARRAY.ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
if (N.ARRAY != N.DIR.ARRAY) stop(sprintf("array has %d tasks but the inputs expect %d gene blocks", N.ARRAY, N.DIR.ARRAY))
OUT.FILE <- sprintf("gene_dir_stage%d_output_%dof%d.rda", STAGE, ARRAY.ID, N.ARRAY)
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

PLAN <- DIR.PLANS[[DIR.PLAN]]
DIR.OUT <- NULL
for (fam in names(DIR.JOBS)) {
  g <- DIR.JOBS[[fam]]; block <- if (length(g)) g[cut(seq_along(g), N.ARRAY, labels = FALSE) == ARRAY.ID] else character(0)   # a family with no genes runs nothing
  if (!length(block)) next
  CHUNK.SIZE <- max(NUM.CORES, ceiling(length(block) / 40))
  chunks <- split(block, ceiling(seq_along(block) / CHUNK.SIZE))
  cat(sprintf("direction bootstrap start: stage %d, block %d of %d, family %s, %d genes, B=%d, %d cores\n", STAGE, ARRAY.ID, N.ARRAY, fam, length(block), PLAN$B, NUM.CORES))
  flush.console()
  t0 <- Sys.time(); done <- 0
  for (k in seq_along(chunks)) {
    ## One gene: every replicate's random split and refits, the margin tests, theta and p of every quantity and scheme of the family
    res <- parLapply(cl, chunks[[k]], dir_contrasts_one, expos = DIR.EXPOS, fits = DIR.FITS, mats = DIR.MATS, etas = DIR.ETAS, plan = PLAN, delta = DIR.DELTA,
                     family = fam, shift = DOM.SHIFT, theta_tail = THETA.TAIL)
    DIR.OUT <- rbind(DIR.OUT, do.call(rbind, res))
    done <- done + length(chunks[[k]]); el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] %d / %d genes  elapsed %.1f min  eta %.1f min\n", format(Sys.time(), "%H:%M:%S"), done, length(block), el, (length(block) - done) * el / done))
    flush.console()
  }
}
save(DIR.OUT, DIR.KEY, file = OUT.FILE)
cat(sprintf("direction bootstrap done: stage %d, block %d of %d, %d rows\n", STAGE, ARRAY.ID, N.ARRAY, if (is.null(DIR.OUT)) 0L else nrow(DIR.OUT)))
stopCluster(cl)
