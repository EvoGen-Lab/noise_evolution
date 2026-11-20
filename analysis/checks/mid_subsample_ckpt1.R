###############################################################
### mid_subsample_ckpt1.R
### MID4: does the hybrid's lower share of zero dispersion fits (6 to 8% against 18 to 20% in the parents thinned to allele depth, MID3) hold at equal
### cell numbers? The hybrid has more cells (about 2,800 against about 800 to 1,200 per parent), and with more cells a small true phi is less often
### estimated at exactly zero, so the gap could come from sample size. Each draw matches both depth and cell number:
###   parents   every UMI kept with probability p (p_c = median Sc allele library / median Pc library, p_e likewise for Se), as in MID3;
###   hybrid    a random subsample of the hybrid cells, of the Pc cell count for Hc and of the Pe cell count for He, without replacement; the hybrid
###             counts are not thinned, so each allele keeps its own depth.
### Offsets: parents their thinned library over the unthinned median parental library; hybrid alleles half the cell's hybrid depth. Fits: the pipeline's
### Cox-Reid NB, with R.C and R.E on the allele dispersions. Per draw and pair (Pc with Hc, Pe with He): the cell numbers, the median count per cell, the
### share of genes at phi = 0, the median log2 phi of the positive fits, and the same on the genes positive in both groups of the pair; and the log2 and
### raw d with thinned parents and subsampled hybrid alleles, overall and by decile of the full parental mean count. Also the hybrid at full cell number
### with the same thinned parents, so the effect of the subsampling alone shows. Genes as MID1 to MID3 (mean count >= MID4.MIN.COUNT per cell in every
### unthinned group). Section 1 checkpoint only. Writes results/checks/MID4_pairs.csv, MID4_d.csv, MID4.rds and a row in summary.csv.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel); library(Matrix)

CHECK.SCRIPT <- "analysis/checks/mid_subsample_ckpt1.R"
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.MID4 <- 6026
if (!exists("MID4.REPS"))      MID4.REPS      <- 3     # draws, each with its own thinning and its own hybrid subsample
if (!exists("MID4.MIN.COUNT")) MID4.MIN.COUNT <- 0.1

stopifnot("mid_subsample_ckpt1.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, C.BAR, R.C, R.E
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

M4.GENE.IDS <- rownames(MIX.SC)[Reduce(`&`, lapply(list(MIX.SC, MIX.SE, HYB.SC, HYB.SE), function(m) rowMeans(m) >= MID4.MIN.COUNT))]
M4.P <- c(Pc = median(colSums(HYB.SC)) / median(colSums(MIX.SC)), Pe = median(colSums(HYB.SE)) / median(colSums(MIX.SE)))
M4.N <- c(Pc = ncol(MIX.SC), Pe = ncol(MIX.SE), H = ncol(HYB.SC))
stopifnot(all(M4.P > 0 & M4.P < 1), "the hybrid must have at least as many cells as each parent to subsample" = M4.N[["H"]] >= max(M4.N[c("Pc", "Pe")]))
cat(sprintf("%d genes; cells Pc %d, Pe %d, hybrid %d; thinning Pc %.3f, Pe %.3f\n", length(M4.GENE.IDS), M4.N[["Pc"]], M4.N[["Pe"]], M4.N[["H"]], M4.P[["Pc"]], M4.P[["Pe"]]))

DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
M4.HYB.EXPO <- (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF / 2
M4.MEAN.COUNT <- 0.5 * (rowMeans(MIX.SC[M4.GENE.IDS, , drop = FALSE]) + rowMeans(MIX.SE[M4.GENE.IDS, , drop = FALSE]))

## Hybrid allele fits on all cells (the comparison without subsampling)
M4.HYB.FULL <- sapply(c(Hc = "Hc", He = "He"), function(k) { m <- if (k == "Hc") HYB.SC else HYB.SE; eta <- if (k == "Hc") R.C else R.E
  unlist(mclapply(M4.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ]), M4.HYB.EXPO, eta, "NB")[["disp"]], mc.cores = CORES)) })
rownames(M4.HYB.FULL) <- M4.GENE.IDS

## Each draw: thinned parents and hybrid alleles subsampled to each parent's cell count, every index and thinning drawn on the draw's own stream
M4.DRAWS <- lapply(seq_len(MID4.REPS), function(r) with_local_seed(SEED.MID4 + r, {
  par <- sapply(c(Pc = "Pc", Pe = "Pe"), function(k) {
    m <- if (k == "Pc") MIX.SC else MIX.SE
    if (inherits(m, "sparseMatrix")) { m <- as(m, "CsparseMatrix"); m@x <- as.numeric(rbinom(length(m@x), as.integer(m@x), M4.P[[k]])); m <- drop0(m) }
    else m[] <- rbinom(length(m), as.integer(m), M4.P[[k]])
    e <- colSums(m) / DEPTH.REF
    unlist(mclapply(M4.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ]), e, 0, "NB")[["disp"]], mc.cores = CORES))
  })
  sub <- list(Hc = sort(sample.int(M4.N[["H"]], M4.N[["Pc"]])), He = sort(sample.int(M4.N[["H"]], M4.N[["Pe"]])))
  hyb <- sapply(c(Hc = "Hc", He = "He"), function(k) { m <- if (k == "Hc") HYB.SC else HYB.SE; eta <- if (k == "Hc") R.C else R.E; ix <- sub[[k]]
    unlist(mclapply(M4.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ix]), M4.HYB.EXPO[ix], eta, "NB")[["disp"]], mc.cores = CORES)) })
  rownames(par) <- rownames(hyb) <- M4.GENE.IDS
  list(par = par, hyb = hyb, sub = sub)
}))

## Pairs: each thinned parent against its allele, subsampled to the parent's cell count and at full cell number
MID4.PAIRS <- do.call(rbind, lapply(seq_len(MID4.REPS), function(r) do.call(rbind, lapply(list(c("Pc", "Hc", "HYB.SC", "MIX.SC"), c("Pe", "He", "HYB.SE", "MIX.SE")), function(pr) {
  dr <- M4.DRAWS[[r]]; pv <- dr$par[, pr[1]]
  do.call(rbind, lapply(c("subsampled", "all cells"), function(v) {
    hv <- if (v == "subsampled") dr$hyb[, pr[2]] else M4.HYB.FULL[, pr[2]]
    nh <- if (v == "subsampled") length(dr$sub[[pr[2]]]) else M4.N[["H"]]
    both <- is.finite(pv) & is.finite(hv) & pv > 0 & hv > 0
    data.frame(draw = r, parent = pr[1], allele = pr[2], hybrid_cells = v, cells_parent = M4.N[[pr[1]]], cells_allele = nh,
               median_count_parent_thinned = median(M4.P[[pr[1]]] * rowMeans(get(pr[4])[M4.GENE.IDS, , drop = FALSE])),
               median_count_allele = median(rowMeans(get(pr[3])[M4.GENE.IDS, , drop = FALSE])),
               share_zero_parent = mean(pv == 0, na.rm = TRUE), share_zero_allele = mean(hv == 0, na.rm = TRUE),
               median_log2_phi_parent_pos = median(log2(pv[pv > 0]), na.rm = TRUE), median_log2_phi_allele_pos = median(log2(hv[hv > 0]), na.rm = TRUE),
               genes_both_positive = sum(both), median_log2_ratio_both = median(log2(hv[both] / pv[both])), row.names = NULL)
  }))
}))))
cat("\nThinned parent against its allele at the parent's cell count (subsampled) and at the full hybrid cell number:\n"); print(MID4.PAIRS, digits = 3, row.names = FALSE)

## d with thinned parents and subsampled hybrid alleles, by decile of the full parental mean count
MID4.D <- do.call(rbind, lapply(seq_len(MID4.REPS), function(r) {
  ph <- cbind(M4.DRAWS[[r]]$par, M4.DRAWS[[r]]$hyb); l <- suppressWarnings(log2(ph)); l[!is.finite(l)] <- NA
  d2 <- 0.5 * (l[, "Hc"] + l[, "He"]) - 0.5 * (l[, "Pc"] + l[, "Pe"]); dr <- 0.5 * (ph[, "Hc"] + ph[, "He"]) - 0.5 * (ph[, "Pc"] + ph[, "Pe"])
  dc <- cut(M4.MEAN.COUNT, unique(quantile(M4.MEAN.COUNT, 0:10 / 10)), include.lowest = TRUE)
  rbind(data.frame(draw = r, bin = "all", genes = sum(is.finite(d2)), median_d_log2 = median(d2, na.rm = TRUE), median_d_raw = median(dr, na.rm = TRUE), row.names = NULL),
        data.frame(draw = r, bin = levels(dc), genes = as.vector(tapply(is.finite(d2), dc, sum)), median_d_log2 = as.vector(tapply(d2, dc, median, na.rm = TRUE)),
                   median_d_raw = as.vector(tapply(dr, dc, median, na.rm = TRUE)), row.names = NULL))
}))
cat("\nMidparent d of the dispersion, thinned parents and subsampled hybrid alleles:\n"); print(MID4.D[MID4.D$bin == "all", ], digits = 3, row.names = FALSE)
print(MID4.D[MID4.D$draw == 1, ], digits = 3, row.names = FALSE)

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(MID4.PAIRS, file.path(CHECK.DIR, "MID4_pairs.csv"), row.names = FALSE)
write.csv(MID4.D, file.path(CHECK.DIR, "MID4_d.csv"), row.names = FALSE)
S4 <- MID4.PAIRS[MID4.PAIRS$hybrid_cells == "subsampled", ]
check_save("MID4", list(pairs = MID4.PAIRS, d = MID4.D, hyb_full = M4.HYB.FULL, draws = M4.DRAWS, p = M4.P, n = M4.N,
                        settings = list(reps = MID4.REPS, min_count = MID4.MIN.COUNT)), "run",
           headline = c(sprintf("draw %d, %s (%d cells) vs %s (%d cells): phi = 0 in %.3f vs %.3f; median log2 phi on genes positive in both: allele - parent %.3f (%d genes)",
                                S4$draw, S4$parent, S4$cells_parent, S4$allele, S4$cells_allele, S4$share_zero_parent, S4$share_zero_allele,
                                S4$median_log2_ratio_both, S4$genes_both_positive),
                        sprintf("draw %d: log2 d %.3f, raw d %.4f", MID4.D$draw[MID4.D$bin == "all"], MID4.D$median_d_log2[MID4.D$bin == "all"], MID4.D$median_d_raw[MID4.D$bin == "all"])),
           mc_error = sprintf("%d draws of thinning and subsampling", MID4.REPS), script = CHECK.SCRIPT)
