###############################################################
### gene_boot.R
### Reads the bootstrap inputs, runs the paired per-gene
### bootstrap on a FORK cluster, saves the output. Mirrors the
### structure of Permutation_Norm0_V2.R.
###
### Inputs  : gene_boot1_inputs.rda, or gene_boot<tag>_inputs.rda for a
###           given tag arg (e.g. "2" for a second-seed adequacy check)
### Output  : gene_boot1_output.rda, or gene_boot<tag>_output.rda to match
###           (BOOT.CONTRASTS, one row per gene)
###
### Progress is printed by the master after each chunk of genes,
### so it lands in the job log even though the workers are forked.
###############################################################

library('parallel')
library('MASS')          # glm.nb ships with base R; set lib= if your cluster needs it

source("functions.R")

## boot_contrasts_one: one gene's bootstrap contrasts across all modes (the unit of work a cluster worker does). draws holds every resample
## index set, drawn before any fit so replicates are independent.
boot_contrasts_one <- function(g, expos, fits, mats, draws) {
  ## The set of draws-list keys a mode's contrast is built from, e.g. "trans"
  ## touches the MIX.SC/MIX.SE and HYT draws, "trans_n" the MIX.SC/MIX.SE and
  ## HYT.N draws. Two modes share draws exactly when this set matches.
  .mode_draw_keys <- function(mode) sort(unique(unname(.DRAW.KEY[.MODES[[mode]]])))

  B <- length(draws); modes <- names(.MODES)
  mcol <- paste0("m.", modes); scol <- paste0("s.", modes); ccol <- paste0("c.", modes)
  M <- matrix(NA_real_, B, 3 * length(modes), dimnames = list(NULL, c(mcol, scol, ccol)))
  for (b in seq_len(B)) {
    d <- draws[[b]]
    fg <- function(mat, ex, idx) .fit_one(mat[g, idx], ex[idx])
    f <- list(
      MIX.SC   = fg(mats$MIX.SC,   expos$MIX.SC, d$MIX.SC),
      MIX.SE   = fg(mats$MIX.SE,   expos$MIX.SE, d$MIX.SE),
      HYC.SC   = fg(mats$HYC.SC,   expos$HYC,    d$HYC),
      HYC.SE   = fg(mats$HYC.SE,   expos$HYC,    d$HYC),
      HYT.SC   = fg(mats$HYT.SC,   expos$HYT,    d$HYT),
      HYT.SE   = fg(mats$HYT.SE,   expos$HYT,    d$HYT),
      HYC.SC.N = fg(mats$HYC.SC.N, expos$HYC.N,  d$HYC.N),
      HYC.SE.N = fg(mats$HYC.SE.N, expos$HYC.N,  d$HYC.N),
      HYT.SC.N = fg(mats$HYT.SC.N, expos$HYT.N,  d$HYT.N),
      HYT.SE.N = fg(mats$HYT.SE.N, expos$HYT.N,  d$HYT.N),
      HYB.COMB = fg(mats$HYB.COMB, expos$HYB,    d$HYB),
      HYB.SC   = fg(mats$HYB.SC,   expos$HYB,    d$HYB),
      HYB.SE   = fg(mats$HYB.SE,   expos$HYB,    d$HYB))
    gv.mu <- lapply(f, function(z) unname(z[["mu"]]))
    gv.bf <- lapply(f, function(z) unname(z[["disp"]]))
    gv.cv <- lapply(f, function(z) .cv2_of(z[["mu"]], z[["disp"]]))
    for (k in seq_along(modes)) {
      grp <- .MODES[[modes[k]]]
      mu <- unlist(gv.mu[grp]); bf <- unlist(gv.bf[grp]); cv <- unlist(gv.cv[grp])
      if (all(is.finite(mu) & mu > 0)) M[b, mcol[k]] <- .contrast_value(modes[k], gv.mu, "MU")
      if (all(is.finite(bf) & bf > 0)) M[b, scol[k]] <- .contrast_value(modes[k], gv.bf, "BFREQ")
      if (all(is.finite(cv) & cv > 0)) M[b, ccol[k]] <- .contrast_value(modes[k], gv.cv, "CV2")
    }
  }
  se   <- apply(M, 2, sd, na.rm = TRUE)
  bdry <- setNames(vapply(modes, function(md) mean(is.na(M[, paste0("s.", md)])), numeric(1)), modes)
  point <- function(q, gv_fun) {
    gv <- setNames(lapply(names(fits), function(grp) gv_fun(fits[[grp]][g, ])), names(fits))
    vapply(modes, function(md) .contrast_value(md, gv, q), numeric(1))
  }
  pm <- point("MU",    function(row) row[["MU"]])
  ps <- point("BFREQ", function(row) row[["DISP"]])
  pc <- point("CV2",   function(row) .cv2_of(row[["MU"]], row[["DISP"]]))
  out <- list(gene = g)
  for (md in .OUT_MODES) {
    smd <- .BFREQ_SOURCE[[md]]
    out[[paste0("mean_", md, "_est")]]   <- unname(pm[md])
    out[[paste0("mean_", md, "_se")]]    <- unname(se[paste0("m.", md)])
    out[[paste0("bfreq_", md, "_est")]]  <- unname(ps[smd])
    out[[paste0("bfreq_", md, "_se")]]   <- unname(se[paste0("s.", smd)])
    out[[paste0("bfreq_", md, "_bdry")]] <- unname(bdry[smd])
    out[[paste0("cv2_", md, "_est")]]    <- unname(pc[smd])
    out[[paste0("cv2_", md, "_se")]]     <- unname(se[paste0("c.", smd)])
    ## cor_<md>: correlation of the mean and bfreq bootstrap draws, defined when both contrasts are
    ## built from the same resampled datasets (matching .mode_draw_keys). For cis and trans the bfreq
    ## value comes from the noise split (HYC.N/HYT.N) and the mean from the mean split (HYC/HYT),
    ## which are resampled with separate index vectors, so cor is NA and downstream code uses 0.
    out[[paste0("cor_", md)]] <- if (identical(.mode_draw_keys(md), .mode_draw_keys(smd))) {
      x <- M[, paste0("m.", md)]; y <- M[, paste0("s.", md)]; ok <- is.finite(x) & is.finite(y)
      if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
    } else NA_real_
  }
  .r2 <- function(c1, c2) {
	x <- M[, c1]; y <- M[, c2]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
  }
  out$cor_dpar_mean  <- .r2("m.dpar_sc", "m.dpar_se")
  out$cor_dpar_bfreq <- .r2("s.dpar_sc", "s.dpar_se")
  data.frame(out, row.names = NULL, check.names = FALSE)
}

## Tag arg picks which input/output pair to use, so the same script
## serves both the primary run and any additional-seed adequacy check
## without duplicating the file. Defaults to "1", the primary run.
##   Rscript gene_boot.R    -> gene_boot1_inputs.rda / gene_boot1_output.rda
##   Rscript gene_boot.R 2  -> gene_boot2_inputs.rda / gene_boot2_output.rda
ARGS <- commandArgs(trailingOnly = TRUE)
TAG  <- if (length(ARGS) >= 1) ARGS[1] else "1"
IN.FILE  <- sprintf("gene_boot%s_inputs.rda", TAG)
OUT.FILE <- sprintf("gene_boot%s_output.rda", TAG)

load(IN.FILE)     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS, N.BOOT, SEED.BOOT

NUM.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))

cl <- makeCluster(NUM.CORES, type = "FORK")
clusterEvalQ(cl, library(MASS))

## ---- chunk the genes so the master can report progress ----
N.UPDATES  <- 40                                              # roughly this many progress lines
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("bootstrap start: %d genes, B=%d, %d cores, %d chunks\n",
            length(GENES), N.BOOT, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  ## One gene's paired bootstrap: refits all 13 groups on each pre-drawn resample, forms every mode's
  ## mean, bfreq (NB size) and CV2 contrast per replicate, and returns the point estimate, bootstrap SE,
  ## boundary fraction and mean-bfreq draw correlation for the eight reportable modes. Replicates where
  ## any needed fit is non-finite give NA for that contrast and are left out of its SD.
  results[[k]] <- parLapply(cl, chunks[[k]], boot_contrasts_one, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, mats = CONTRAST.MATS, draws = DRAWS)
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(GENES),
              100 * done / length(GENES), el, eta))
  flush.console()
}

BOOT.CONTRASTS <- do.call(rbind, unlist(results, recursive = FALSE))
save(BOOT.CONTRASTS, file = OUT.FILE)

cat(sprintf("bootstrap done: %d genes in %.1f min [tag %s]\n",
            nrow(BOOT.CONTRASTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            TAG))

stopCluster(cl)
