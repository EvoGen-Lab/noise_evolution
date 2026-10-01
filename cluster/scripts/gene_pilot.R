###############################################################
### gene_pilot.R
### Reads the pilot inputs, runs the per-gene internal-pilot
### bootstrap (SE of log mu and log disp, undivided hybrid and
### parents) on a FORK cluster, saves the output. Mirrors the
### structure of gene_boot.R exactly; this is the same
### kind of per-gene, per-replicate bootstrap, just over four
### datasets instead of thirteen, and used only to set the Nc/Nt
### split fraction rather than to test any gene for significance.
###
### Inputs  : gene_pilot_inputs.rda  (PILOT.MATS, PILOT.EXPOS,
###           PILOT.GENES, N.BOOT, SEED.BOOT)
### Output  : gene_pilot_output.rda  (PILOT.SE, one row per gene: bootstrap SEs of
###           log mu and log disp for the four datasets, plus the bootstrap covariance of
###           the two hybrid alleles' estimates, which share their resampled cells)
###############################################################

library('parallel')
library('MASS')

source("functions.R")

load("gene_pilot_inputs.rda")  # PILOT.MATS, PILOT.EXPOS, PILOT.GENES, N.BOOT, SEED.BOOT

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(PILOT.GENES) / N.UPDATES))
chunks     <- split(PILOT.GENES, ceiling(seq_along(PILOT.GENES) / CHUNK.SIZE))

cat(sprintf("pilot start: %d genes, B=%d, %d cores, %d chunks\n",
            length(PILOT.GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's bootstrap SE of log(mu) and log(size) for the four datasets needed to compute
  ## A and B, plus the bootstrap covariance between the two hybrid alleles. HYB.SC and HYB.SE
  ## are refit on the SAME resampled cell indices in each replicate because the two alleles
  ## are measured in the same cells; their log estimates therefore co-vary across replicates
  ## (HYB_logmu_cov, HYB_logdisp_cov), and the variance of the allele contrast is
  ## var(SC) + var(SE) - 2 cov(SC, SE), which estimate_f_star() uses for A.
  ##
  ## The RNG is seeded from a hash of the gene's own name and the resamples are
  ## drawn inside this function. Each gene's draws therefore depend only on its
  ## name, so results are identical for any gene chunking or core count of the
  ## FORK cluster that distributes the genes.
  results[[k]] <- parLapply(cl, chunks[[k]], function(g) {
    set.seed(SEED.BOOT + sum(utf8ToInt(g)))
    n.p.sc <- length(PILOT.EXPOS$MIX.SC); n.p.se <- length(PILOT.EXPOS$MIX.SE); n.h <- length(PILOT.EXPOS$HYB)
    logmu <- matrix(NA_real_, N.BOOT, 4, dimnames = list(NULL, c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")))
    logsz <- logmu
    for (b in seq_len(N.BOOT)) {
      i.p.sc <- sample.int(n.p.sc, n.p.sc, replace = TRUE)
      i.p.se <- sample.int(n.p.se, n.p.se, replace = TRUE)
      i.h    <- sample.int(n.h,    n.h,    replace = TRUE)
      f.mc <- .fit_one(PILOT.MATS$MIX.SC[g, i.p.sc], PILOT.EXPOS$MIX.SC[i.p.sc])
      f.me <- .fit_one(PILOT.MATS$MIX.SE[g, i.p.se], PILOT.EXPOS$MIX.SE[i.p.se])
      f.hc <- .fit_one(PILOT.MATS$HYB.SC[g, i.h],    PILOT.EXPOS$HYB[i.h])
      f.he <- .fit_one(PILOT.MATS$HYB.SE[g, i.h],    PILOT.EXPOS$HYB[i.h])
      if (is.finite(f.mc["mu"])   && f.mc["mu"]   > 0) logmu[b, "MIX.SC"] <- log(f.mc["mu"])
      if (is.finite(f.me["mu"])   && f.me["mu"]   > 0) logmu[b, "MIX.SE"] <- log(f.me["mu"])
      if (is.finite(f.hc["mu"])   && f.hc["mu"]   > 0) logmu[b, "HYB.SC"] <- log(f.hc["mu"])
      if (is.finite(f.he["mu"])   && f.he["mu"]   > 0) logmu[b, "HYB.SE"] <- log(f.he["mu"])
      if (is.finite(f.mc["disp"]) && f.mc["disp"] > 0) logsz[b, "MIX.SC"] <- log(f.mc["disp"])
      if (is.finite(f.me["disp"]) && f.me["disp"] > 0) logsz[b, "MIX.SE"] <- log(f.me["disp"])
      if (is.finite(f.hc["disp"]) && f.hc["disp"] > 0) logsz[b, "HYB.SC"] <- log(f.hc["disp"])
      if (is.finite(f.he["disp"]) && f.he["disp"] > 0) logsz[b, "HYB.SE"] <- log(f.he["disp"])
    }
    allele_cov <- function(m) {
      ok <- is.finite(m[, "HYB.SC"]) & is.finite(m[, "HYB.SE"])
      if (sum(ok) > 2) cov(m[ok, "HYB.SC"], m[ok, "HYB.SE"]) else NA_real_
    }
    data.frame(
      gene = g,
      MIX.SC_logmu_se   = sd(logmu[, "MIX.SC"], na.rm = TRUE),
      MIX.SE_logmu_se   = sd(logmu[, "MIX.SE"], na.rm = TRUE),
      HYB.SC_logmu_se   = sd(logmu[, "HYB.SC"], na.rm = TRUE),
      HYB.SE_logmu_se   = sd(logmu[, "HYB.SE"], na.rm = TRUE),
      MIX.SC_logdisp_se = sd(logsz[, "MIX.SC"], na.rm = TRUE),
      MIX.SE_logdisp_se = sd(logsz[, "MIX.SE"], na.rm = TRUE),
      HYB.SC_logdisp_se = sd(logsz[, "HYB.SC"], na.rm = TRUE),
      HYB.SE_logdisp_se = sd(logsz[, "HYB.SE"], na.rm = TRUE),
      HYB_logmu_cov     = allele_cov(logmu),
      HYB_logdisp_cov   = allele_cov(logsz),
      row.names = NULL, check.names = FALSE)
  })
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(PILOT.GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(PILOT.GENES),
              100 * done / length(PILOT.GENES), el, eta))
  flush.console()
}

PILOT.SE <- do.call(rbind, unlist(results, recursive = FALSE))
save(PILOT.SE, file = "gene_pilot_output.rda")

cat(sprintf("pilot done: %d genes in %.1f min\n",
            nrow(PILOT.SE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

stopCluster(cl)
