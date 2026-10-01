###############################################################
### gene_perm.R
### Reads the same inputs as the bootstrap, runs the paired
### per-gene permutation null on a FORK cluster, saves the output.
### Mirrors gene_boot.R.
###
### Inputs  : gene_perm_inputs.rda  (CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS GENES, ...)
### Output  : gene_perm_output.rda
###           PERM.RESULTS, one row per gene: mean_<mode>_obs / _p,
###           bfreq_<mode>_obs / _p (the NB dispersion parameter,
###           reported under its biological name), bsize_<mode>_obs / _p
###           (mean minus bfreq, its own null), kbal_<mode>_obs / _p
###           (bfreq minus bsize, its own null), and cv2_<mode>_obs / _p
###           (1/mean + 1/bfreq, its own null).
###           For dpar_sc and dpar_se, the noise quantities also get
###           _p_ploidy and _p_ind columns, testing the contrast against
###           the ploidy-expected value using this gene's own shift and
###           using the independent-allele shift of 1, plus the
###           ploidy_shift column recording the shift used.
###
### Each contrast is nulled by its own label shuffle. total, cis, and
### trans are zero-centred under the null they test. dpar tests
### exchangeability of hybrid and parent cells, which is the null the
### ploidy-shifted observed value is read against (ploidy_shift() and
### ploidy_adjust_dpar() in functions.R).
###############################################################

library('parallel')
library('MASS')

source("functions.R")

## One gene's permutation null: refits the relabeled groups for every mode and permutation, compares
## the observed contrasts (mean, bfreq, CV2, bsize, kbal) with the null by a two-sided permutation p, and
## adds ploidy-adjusted p-values for dpar when ploidy_shift is supplied. bsize and kbal nulls recombine
## the mean and bfreq null draws of the same permutation.
permute_contrasts_one <- function(g, mats, expos, fits, perms, ploidy_shift = NULL) {
  l2 <- function(x) log2(x); B <- length(perms); modes <- names(.MODES)
  Nm <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Ns <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Nc <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  gv   <- setNames(lapply(names(fits), function(grp) fits[[grp]][g, ]), names(fits))
  gvmu <- lapply(gv, function(z) z[["MU"]]); gvbf <- lapply(gv, function(z) z[["DISP"]])
  gvcv <- lapply(gv, function(z) .cv2_of(z[["MU"]], z[["DISP"]]))
  obs.m <- vapply(modes, function(m) .contrast_value(m, gvmu, "MU"),    numeric(1))
  obs.s <- vapply(modes, function(m) .contrast_value(m, gvbf, "BFREQ"), numeric(1))
  obs.c <- vapply(modes, function(m) .contrast_value(m, gvcv, "CV2"),   numeric(1))
  cMIXsc <- mats$MIX.SC[g, ]; cMIXse <- mats$MIX.SE[g, ]
  cHYCsc <- mats$HYC.SC[g, ]; cHYCse <- mats$HYC.SE[g, ]
  cHYTsc <- mats$HYT.SC[g, ]; cHYTse <- mats$HYT.SE[g, ]
  cHYCsc.N <- mats$HYC.SC.N[g, ]; cHYCse.N <- mats$HYC.SE.N[g, ]
  cHYTsc.N <- mats$HYT.SC.N[g, ]; cHYTse.N <- mats$HYT.SE.N[g, ]
  cHYBsc <- mats$HYB.SC[g, ]; cHYBse <- mats$HYB.SE[g, ]; cHYBc <- mats$HYB.COMB[g, ]
  mp.m <- l2((gvmu[["MIX.SC"]] + gvmu[["MIX.SE"]]) / 2)
  mp.s <- 0.5 * (l2(gvbf[["MIX.SC"]]) + l2(gvbf[["MIX.SE"]]))
  mp.c <- 0.5 * (l2(gvcv[["MIX.SC"]]) + l2(gvcv[["MIX.SE"]]))
  nSC <- length(cMIXsc); nSE <- length(cMIXse); nHYB <- length(cHYBc)
  rat <- function(a, b, q) { x <- a[q]; y <- b[q]
    if (is.finite(x) && x > 0 && is.finite(y) && y > 0) l2(x) - l2(y) else NA_real_ }
  ## CV2 contrast of two .fit_one()/.fit_split() results
  ratcv2 <- function(a, b) { x <- .cv2_of(a[["mu"]], a[["disp"]]); y <- .cv2_of(b[["mu"]], b[["disp"]])
    if (is.finite(x) && x > 0 && is.finite(y) && y > 0) l2(x) - l2(y) else NA_real_ }
  for (b in seq_len(B)) {
    p <- perms[[b]]
    sp <- .fit_split(c(cMIXsc, cMIXse), c(expos$MIX.SC, expos$MIX.SE), p$total, nSC)
    Nm[b,"total"] <- rat(sp$a, sp$b, "mu"); Ns[b,"total"] <- rat(sp$a, sp$b, "disp"); Nc[b,"total"] <- ratcv2(sp$a, sp$b)

    sc <- ifelse(p$cis, cHYCse, cHYCsc); se <- ifelse(p$cis, cHYCsc, cHYCse)
    fa <- .fit_one(sc, expos$HYC); fb <- .fit_one(se, expos$HYC)
    Nm[b,"cis"] <- rat(fa, fb, "mu"); Ns[b,"cis"] <- rat(fa, fb, "disp"); Nc[b,"cis"] <- ratcv2(fa, fb)

    ## Noise-split (f_disp) version of the cis null, used only for
    ## the bfreq_cis/cv2_cis output columns, mirroring cis_n in .MODES
    sc.n <- ifelse(p$cis_n, cHYCse.N, cHYCsc.N); se.n <- ifelse(p$cis_n, cHYCsc.N, cHYCse.N)
    fa.n <- .fit_one(sc.n, expos$HYC.N); fb.n <- .fit_one(se.n, expos$HYC.N)
    Nm[b,"cis_n"] <- rat(fa.n, fb.n, "mu"); Ns[b,"cis_n"] <- rat(fa.n, fb.n, "disp"); Nc[b,"cis_n"] <- ratcv2(fa.n, fb.n)

    sSC <- .fit_split(c(cMIXsc, cHYTsc), c(expos$MIX.SC, expos$HYT), p$transSC, nSC)
    sSE <- .fit_split(c(cMIXse, cHYTse), c(expos$MIX.SE, expos$HYT), p$transSE, nSE)
    Nm[b,"trans"] <- rat(sSC$a, sSE$a, "mu")   - rat(sSC$b, sSE$b, "mu")
    Ns[b,"trans"] <- rat(sSC$a, sSE$a, "disp") - rat(sSC$b, sSE$b, "disp")
    Nc[b,"trans"] <- ratcv2(sSC$a, sSE$a)      - ratcv2(sSC$b, sSE$b)

    ## Noise-split version of the trans null, used only for the
    ## bfreq_trans/cv2_trans output columns, mirroring trans_n in .MODES
    sSC.n <- .fit_split(c(cMIXsc, cHYTsc.N), c(expos$MIX.SC, expos$HYT.N), p$transSC_n, nSC)
    sSE.n <- .fit_split(c(cMIXse, cHYTse.N), c(expos$MIX.SE, expos$HYT.N), p$transSE_n, nSE)
    Nm[b,"trans_n"] <- rat(sSC.n$a, sSE.n$a, "mu")   - rat(sSC.n$b, sSE.n$b, "mu")
    Ns[b,"trans_n"] <- rat(sSC.n$a, sSE.n$a, "disp") - rat(sSC.n$b, sSE.n$b, "disp")
    Nc[b,"trans_n"] <- ratcv2(sSC.n$a, sSE.n$a)      - ratcv2(sSC.n$b, sSE.n$b)

    synth <- cMIXsc[p$dom_i] + cMIXse[p$dom_j]
    es    <- expos$MIX.SC[p$dom_i] + expos$MIX.SE[p$dom_j]
    fs <- .fit_one(synth, es)
    if (is.finite(fs["mu"])   && fs["mu"]   > 0) Nm[b,"dom"] <- l2(fs["mu"])   - mp.m
    if (is.finite(fs["disp"]) && fs["disp"] > 0) Ns[b,"dom"] <- l2(fs["disp"]) - mp.s
    fs.cv <- .cv2_of(fs[["mu"]], fs[["disp"]])
    if (is.finite(fs.cv) && fs.cv > 0) Nc[b,"dom"] <- l2(fs.cv) - mp.c
    d1 <- .fit_split(c(cHYBc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$dparSC, nHYB)
    Nm[b,"dpar_sc"] <- rat(d1$a, d1$b, "mu"); Ns[b,"dpar_sc"] <- rat(d1$a, d1$b, "disp"); Nc[b,"dpar_sc"] <- ratcv2(d1$a, d1$b)
    d2 <- .fit_split(c(cHYBc, cMIXse), c(expos$HYB, expos$MIX.SE), p$dparSE, nHYB)
    Nm[b,"dpar_se"] <- rat(d2$a, d2$b, "mu"); Ns[b,"dpar_se"] <- rat(d2$a, d2$b, "disp"); Nc[b,"dpar_se"] <- ratcv2(d2$a, d2$b)
    i1 <- .fit_split(c(cHYBsc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$inhSC, nHYB)
    Nm[b,"inh_sc"] <- rat(i1$a, i1$b, "mu"); Ns[b,"inh_sc"] <- rat(i1$a, i1$b, "disp"); Nc[b,"inh_sc"] <- ratcv2(i1$a, i1$b)
    i2 <- .fit_split(c(cHYBse, cMIXse), c(expos$HYB, expos$MIX.SE), p$inhSE, nHYB)
    Nm[b,"inh_se"] <- rat(i2$a, i2$b, "mu"); Ns[b,"inh_se"] <- rat(i2$a, i2$b, "disp"); Nc[b,"inh_se"] <- ratcv2(i2$a, i2$b)
  }
  out <- list(gene = g)
  for (md in .OUT_MODES) {
    smd <- .BFREQ_SOURCE[[md]]
    out[[paste0("mean_", md, "_obs")]]  <- unname(obs.m[md])
    out[[paste0("mean_", md, "_p")]]    <- perm_pval(obs.m[md], Nm[, md])
    out[[paste0("bfreq_", md, "_obs")]] <- unname(obs.s[smd])
    out[[paste0("bfreq_", md, "_p")]]   <- perm_pval(obs.s[smd], Ns[, smd])
    out[[paste0("cv2_", md, "_obs")]]   <- unname(obs.c[smd])
    out[[paste0("cv2_", md, "_p")]]     <- perm_pval(obs.c[smd], Nc[, smd])

    ## bsize null: mean null minus bfreq null, permutation by permutation. Where md == smd the two
    ## draws share one relabeling, so their correlation carries through. For cis and trans they come
    ## from separate relabelings (mean split vs noise split) and are independent.
    bs_null <- Nm[, md] - Ns[, smd]
    obs.bs  <- unname(obs.m[md]) - unname(obs.s[smd])
    out[[paste0("bsize_", md, "_obs")]] <- obs.bs
    out[[paste0("bsize_", md, "_p")]]   <- perm_pval(obs.bs, bs_null)

    ## kbal = bfreq - bsize; its null is the bfreq null minus the bsize null from the same draws, so
    ## kbal gets a permutation p-value for every mode, usable by classify_reg()/classify_dom().
    kbal_null <- Ns[, smd] - bs_null
    obs.kbal  <- unname(obs.s[smd]) - obs.bs
    out[[paste0("kbal_", md, "_obs")]] <- obs.kbal
    out[[paste0("kbal_", md, "_p")]]   <- perm_pval(obs.kbal, kbal_null)

    ## Ploidy-adjusted p-values for the hybrid-versus-parent noise contrasts.
    ## HYB.COMB sums two alleles, which raises bfreq by the shift s (see
    ## ploidy_shift). Removing s from the observed contrast and reading it
    ## against the same zero-centred exchangeability null tests hybrid noise
    ## at the per-genome scale of the parents. _p_ploidy uses this gene's own
    ## s and _p_ind uses s = 1, the fully independent equal-allele bound.
    ## bfreq falls by s, bsize rises by s, kbal falls by 2s, and cv2 rebuilds
    ## with its latent term scaled by 2^s. The mean axis needs no adjustment.
    if (!is.null(ploidy_shift) && md %in% c("dpar_sc", "dpar_se")) {
      par.nm <- if (md == "dpar_sc") "MIX.SC" else "MIX.SE"
      mh <- gvmu[["HYB.COMB"]]; kh <- gvbf[["HYB.COMB"]]; cvp <- gvcv[[par.nm]]
      adj <- function(s) {
        cvh <- if (is.finite(s) && is.finite(mh) && mh > 0 && is.finite(kh) && kh > 0) 1 / mh + 2^s / kh else NA_real_
        cv2 <- if (is.finite(cvh) && cvh > 0 && is.finite(cvp) && cvp > 0) l2(cvh) - l2(cvp) else NA_real_
        list(bfreq = unname(obs.s[smd]) - s, bsize = obs.bs + s, kbal = obs.kbal - 2 * s, cv2 = cv2)
      }
      nulls <- list(bfreq = Ns[, smd], bsize = bs_null, kbal = kbal_null, cv2 = Nc[, smd])
      s.own <- unname(ploidy_shift[g])
      for (tag in c("ploidy", "ind")) {
        ob <- adj(if (tag == "ploidy") s.own else 1)
        for (q in names(nulls)) out[[paste0(q, "_", md, "_p_", tag)]] <- perm_pval(ob[[q]], nulls[[q]])
      }
    }
  }
  if (!is.null(ploidy_shift)) out$ploidy_shift <- unname(ploidy_shift[g])
  data.frame(out, row.names = NULL, check.names = FALSE)
}
load("gene_perm_inputs.rda")     # CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, SEED.PERM, PLOIDY.SHIFT (when present)

NUM.CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores() - 1L)))
## Ploidy shift. When permute_contrasts_one() accepts a ploidy_shift argument, the
## dpar contrasts are read against the ploidy-expected value and PLOIDY.SHIFT must
## come with the inputs. The two are checked together so the shift can never be
## dropped silently. Forked workers inherit PLOIDY.SHIFT from this process.
HAS.PLOIDY.ARG <- "ploidy_shift" %in% names(formals(permute_contrasts_one))
HAS.PLOIDY.VAL <- exists("PLOIDY.SHIFT")
if (HAS.PLOIDY.ARG != HAS.PLOIDY.VAL)
  stop(sprintf("ploidy shift mismatch: permute_contrasts_one() takes ploidy_shift = %s, PLOIDY.SHIFT in inputs = %s",
               HAS.PLOIDY.ARG, HAS.PLOIDY.VAL))
EXTRA.ARGS <- if (HAS.PLOIDY.ARG) list(ploidy_shift = PLOIDY.SHIFT) else list()
cat(sprintf("ploidy shift active: %s\n", HAS.PLOIDY.ARG))
if (!HAS.PLOIDY.ARG) cat("note: this functions.R has no ploidy_shift, so PERM.RESULTS will hold no _p_ploidy or _p_ind columns\n")
flush.console()

## Forked workers inherit PERMS, the fits, and the loaded functions from this
## process without copying them, so one shared PERMS object serves every core.
## mc.preschedule = FALSE hands genes to workers one at a time, which keeps all
## cores busy when genes differ in cost.

N.UPDATES  <- 40
CHUNK.SIZE <- max(NUM.CORES, ceiling(length(GENES) / N.UPDATES))
chunks     <- split(GENES, ceiling(seq_along(GENES) / CHUNK.SIZE))

cat(sprintf("permutation start: %d genes, NPERM=%d, %d cores, %d chunks\n",
            length(GENES), N.PERM, NUM.CORES, length(chunks)))
flush.console()

results <- vector("list", length(chunks))
t0 <- Sys.time(); done <- 0
for (k in seq_along(chunks)) {
  results[[k]] <- do.call(mclapply, c(list(X = chunks[[k]], FUN = permute_contrasts_one,
                           mats = CONTRAST.MATS, expos = CONTRAST.EXPOS, fits = CONTRAST.FITS, perms = PERMS,
                           mc.cores = NUM.CORES, mc.preschedule = FALSE), EXTRA.ARGS))
  bad <- vapply(results[[k]], function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop(sprintf("%d gene(s) failed in chunk %d, first: %s", sum(bad), k, chunks[[k]][which(bad)[1]]))
  done <- done + length(chunks[[k]])
  el   <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta  <- if (el > 0) (length(GENES) - done) * (el / done) else NA_real_
  cat(sprintf("[%s] %d / %d genes (%.1f%%)  elapsed %.1f min  eta %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), done, length(GENES),
              100 * done / length(GENES), el, eta))
  flush.console()
}

PERM.RESULTS <- do.call(rbind, unlist(results, recursive = FALSE))
save(PERM.RESULTS, file = "gene_perm_output.rda")

cat(sprintf("permutation done: %d genes in %.1f min\n",
            nrow(PERM.RESULTS),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

