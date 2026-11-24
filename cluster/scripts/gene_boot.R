###############################################################
### gene_boot.R
### Reads the bootstrap inputs, runs the paired per-gene bootstrap on a FORK cluster, saves the output. The replicates are drawn in stages (Methods: Bootstrap
### tests of significance), set by the environment variable GENE_BOOT_STAGE (default 1). Every replicate's indices come from a draw plan (make_draw_plan()):
### draw b of a plan is drawn from its own seed when it is used, so the jobs ship no index files and every gene of a stage sees the same draws.
###
### Stage 1 (job array over seed replicates, gene_boot.sub)
###   Inputs  : gene_boot_inputs.rda (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, N.BOOT, BOOT.PLANS: the draw plans of every stage, $rep holding one
###             plan of N.BOOT draws per seed replicate; CONTRAST.ETAS: the split variance of every cell of every dataset, 0 for the parents, which every
###             replicate reads for its resampled cells; GATES: the contrasts each gene's observed fits define)
###   Output  : gene_boot_output_<k>of<K>.rda for array task k of K, holding BOOT.CONTRASTS (one row per gene) for seed replicate k, unshifted, 1,000 replicates.
###             Replicate 1 is the first stage; the others check seed adequacy. analysis.R picks the calibrators from the stage-1 BH q-values of the log2 total
###             and trans dispersion contrasts and fits the shared shift Xi once.
### Stage 2 (one task, gene_boot_shift.sub; after analysis.R has written gene_shift_inputs.rda)
###   Inputs  : gene_boot_inputs.rda and gene_shift_inputs.rda (DOM.SHIFT: Xi at each gene's log2 mu_H, named by gene, held fixed in every replicate)
###   Output  : gene_boot_shift_output.rda, holding BOOT.CONTRASTS for the stage-1 draws of replicate 1 with the shift: the stage-1 columns (unchanged), the
###             shifted contrasts a<quantity>_<mode>_* and the degrees of dominance k_<quantity>_* and kt_<quantity>_*.
### Stages 3 and 4 (job arrays over gene blocks, gene_boot_stage.sub; after analysis.R has written gene_stage<s>_inputs.rda)
###   Inputs  : gene_boot_inputs.rda, gene_shift_inputs.rda and gene_stage<s>_inputs.rda (STAGE.JOBS: the genes of each family, mean and disp, that continue;
###             STAGE.PLAN: the name of their plan in BOOT.PLANS, fresh draws of 10,000 (stage 3) or 100,000 (stage 4); STAGE.KEY; N.STAGE.ARRAY)
###   Output  : gene_boot_stage<s>_output_<k>of<K>.rda, holding BOOT.MEAN and BOOT.DISP (boot_contrasts_one(family = "mean" / "disp")) for gene block k of each
###             family, and STAGE.KEY.
### Stage 5 (job array over gene blocks, gene_boot_stage.sub; after the final q-values)
###   Inputs  : as stage 4, with STAGE.JOBS$disp the genes classed by k_P or k_T and K.LEVEL the false coverage-statement rate of each k family
###   Output  : gene_boot_stage5_output_<k>of<K>.rda, holding BOOT.K (the degrees of dominance with their intervals at the coverage level) and STAGE.KEY.
###
### Progress is printed by the master after each chunk of genes, so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')
library('MASS')          # glm.nb ships with base R; set lib= if your cluster needs it

source("functions.R")
STAGE <- as.integer(Sys.getenv("GENE_BOOT_STAGE", unset = "1"))
stopifnot(STAGE %in% 1:5)
write_pkg_versions(c("gene_boot", "gene_boot_shift", "gene_boot_stage3", "gene_boot_stage4", "gene_boot_stage5")[STAGE])   # compared with renv.lock by check_pkg_versions()

load("gene_boot_inputs.rda")     # GENES is the inference set (every gene for which at least one contrast is defined). CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, BOOT.PLANS, N.BOOT, CONTRAST.ETAS, GATES
SHIFT <- NULL
if (STAGE >= 2) {
  load("gene_shift_inputs.rda")  # DOM.SHIFT
  stopifnot(is.numeric(DOM.SHIFT), all(names(DOM.SHIFT) %in% GENES), all(is.finite(DOM.SHIFT)))
  SHIFT <- DOM.SHIFT
}
K.LEVEL <- NULL
if (STAGE >= 3) {
  load(sprintf("gene_stage%d_inputs.rda", STAGE))   # STAGE.JOBS, STAGE.PLAN, STAGE.KEY, N.STAGE.ARRAY (stage 5 also K.LEVEL)
  stopifnot(is.list(STAGE.JOBS), all(names(STAGE.JOBS) %in% c("mean", "disp")), all(unlist(STAGE.JOBS) %in% GENES), STAGE.PLAN %in% names(BOOT.PLANS))
}

## Job-array position. Stage 1: task k bootstraps seed replicate k. Stage 2: replicate 1 alone. Stages 3 to 5: task k runs gene block k of each family. A
## finished task leaves its output file, so resubmitting an array computes only the missing tasks.
ARRAY.ID <- if (STAGE == 2) 1L else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
N.ARRAY  <- if (STAGE == 2) 1L else as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = "1"))
N.EXPECT <- switch(as.character(STAGE), "1" = length(BOOT.PLANS$rep), "2" = 1L, N.STAGE.ARRAY)
if (N.ARRAY != N.EXPECT) stop(sprintf("array has %d tasks but the inputs expect %d", N.ARRAY, N.EXPECT))
OUT.FILE <- switch(as.character(STAGE), "1" = sprintf("gene_boot_output_%dof%d.rda", ARRAY.ID, N.ARRAY), "2" = "gene_boot_shift_output.rda",
                   sprintf("gene_boot_stage%d_output_%dof%d.rda", STAGE, ARRAY.ID, N.ARRAY))
if (file.exists(OUT.FILE)) { cat(sprintf("%s already exists, nothing to do\n", OUT.FILE)); quit(save = "no") }

## The work of this task: one entry per family, each a gene list, a draw plan and the arguments of boot_contrasts_one(). Stages 1 and 2 run every column
## ("all") with the moment estimator beside the NB fit; later stages refit one family at a time with the NB estimator alone (the moment columns and the
## mean-dispersion correlations stay those of the first stage).
JOBS <- if (STAGE <= 2) list(all = list(genes = GENES, plan = BOOT.PLANS$rep[[ARRAY.ID]], family = "all", methods = c("NB", "weighted"))) else
  lapply(setNames(names(STAGE.JOBS), names(STAGE.JOBS)), function(fam) { g <- STAGE.JOBS[[fam]]
    list(genes = if (length(g)) g[cut(seq_along(g), N.ARRAY, labels = FALSE) == ARRAY.ID] else character(0), plan = BOOT.PLANS[[STAGE.PLAN]], family = fam, methods = "NB") })

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

OUT <- list()
for (fam in names(JOBS)) {
  job <- JOBS[[fam]]
  ## ---- chunk the genes so the master can report progress ----
  N.UPDATES  <- 40                                              # roughly this many progress lines
  CHUNK.SIZE <- max(NUM.CORES, ceiling(length(job$genes) / N.UPDATES))
  chunks     <- split(job$genes, ceiling(seq_along(job$genes) / CHUNK.SIZE))
  cat(sprintf("bootstrap start: stage %d, task %d of %d, family %s, %d genes, B=%d, %d cores, %d chunks\n",
              STAGE, ARRAY.ID, N.ARRAY, fam, length(job$genes), job$plan$B, NUM.CORES, length(chunks)))
  flush.console()
  results <- vector("list", length(chunks))
  t0 <- Sys.time(); done <- 0
  for (k in seq_along(chunks)) {
    ## One gene's paired bootstrap: refits the contrast groups of its family on each draw of the plan with the split variances of its cells, forms every
    ## contrast per replicate, and returns the estimate, bootstrap SE, percentile interval and sign-test p-value of every contrast (with the shift from stage 2).
    results[[k]] <- parLapply(cl, chunks[[k]], boot_contrasts_one, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS, draws = job$plan,
                              etas = CONTRAST.ETAS, gates = GATES, methods = job$methods, shift = SHIFT, family = job$family, k_level = K.LEVEL)
    done <- done + length(chunks[[k]])
    el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta  <- if (el > 0) (length(job$genes) - done) * (el / done) else NA_real_
    cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
                format(Sys.time(), "%H:%M:%S"), done, length(job$genes), 100 * done / max(1, length(job$genes)), el, eta))
    flush.console()
  }
  OUT[[fam]] <- if (length(results)) do.call(rbind, unlist(results, recursive = FALSE)) else NULL
  cat(sprintf("family %s done: %d genes in %.1f min\n", fam, length(job$genes), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

if (STAGE <= 2) {
  BOOT.CONTRASTS <- OUT$all
  save(BOOT.CONTRASTS, file = OUT.FILE)
} else if (STAGE < 5) {
  BOOT.MEAN <- OUT$mean; BOOT.DISP <- OUT$disp
  save(BOOT.MEAN, BOOT.DISP, STAGE.KEY, file = OUT.FILE)
} else {
  BOOT.K <- OUT$disp
  save(BOOT.K, STAGE.KEY, file = OUT.FILE)
}
cat(sprintf("bootstrap done: stage %d, task %d of %d\n", STAGE, ARRAY.ID, N.ARRAY))

stopCluster(cl)
