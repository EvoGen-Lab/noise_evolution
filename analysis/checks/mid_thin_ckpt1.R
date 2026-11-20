###############################################################
### mid_thin_ckpt1.R
### MID3: is the hybrid's higher dispersion (MID1: median log2 d of about 1.4, falling with expression) a depth effect? Each hybrid allele has about 0.6
### counts per cell against 2.6 to 5.3 in the parents (MID2). Binomial thinning of UMIs leaves the NB dispersion unchanged in theory (a thinned NB is NB
### with the same phi and a scaled mean; check A8), so thinning each parent to its allele's depth puts the parents at the hybrid's counts with their own
### true phi. If the shift is depth (estimator bias in sparse cells, and the selection of positive estimates that the log2 scale needs), the thinned
### parents' phi rises toward the hybrid's and d falls toward zero; if the hybrid is truly noisier, d stays.
### Thinning: every UMI of Pc is kept with probability p_c = median Sc allele library / median Pc library, and of Pe with p_e = median Se allele library /
### median Pe library (allele library: the cell's column sum of HYB.SC or HYB.SE), so each parent matches the depth of its own allele. MID3.REPS draws.
### Offsets are recomputed from the thinned libraries, as the pipeline would; the hybrid fits are the unthinned ones throughout.
### Per draw and overall: per group the share of genes at phi = 0 and the median log2 phi of the positive fits; the log2 d (geometric midpoints) and the
### raw d (arithmetic) with full and with thinned parents, on the genes where each is defined and on the genes where both are, overall and by decile of
### the full mean count. Fits are the pipeline's Cox-Reid NB with the library split correction of the allele dispersions (SPLIT.ETA). Genes need a mean count of
### at least MID3.MIN.COUNT per cell in every unthinned group, as MID1. Section 1 checkpoint only.
### Writes results/checks/MID3_thin.csv, MID3_thin_groups.csv, MID3.rds and a row in summary.csv.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel); library(Matrix)

CHECK.SCRIPT <- "analysis/checks/mid_thin_ckpt1.R"
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
SEED.MID3 <- 5026
if (!exists("MID3.REPS"))      MID3.REPS      <- 3     # thinning draws
if (!exists("MID3.MIN.COUNT")) MID3.MIN.COUNT <- 0.1   # minimum mean count per cell in every unthinned group

stopifnot("mid_thin_ckpt1.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

## Genes and thinning fractions
M3.GENE.IDS <- rownames(MIX.SC)[Reduce(`&`, lapply(list(MIX.SC, MIX.SE, HYB.SC, HYB.SE), function(m) rowMeans(m) >= MID3.MIN.COUNT))]
M3.P <- c(Pc = median(colSums(HYB.SC)) / median(colSums(MIX.SC)), Pe = median(colSums(HYB.SE)) / median(colSums(MIX.SE)))
stopifnot("an allele library is not shallower than its parent's; thinning needs p < 1" = all(M3.P > 0 & M3.P < 1))
cat(sprintf("%d genes; thinning fractions: Pc %.3f, Pe %.3f\n", length(M3.GENE.IDS), M3.P[["Pc"]], M3.P[["Pe"]]))

## Fits of one group: Cox-Reid phi (split corrected) of every gene, against offsets from the group's own libraries over the reference depth.
## Offsets as in Section 2: parents library / DEPTH.REF, hybrid alleles the library of their own genome / DEPTH.REF. DEPTH.REF stays the unthinned median parental
## library, so a thinned parent's offsets scale with its thinning and its fitted mean falls by p while phi is unchanged.
DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
M3.FIT.HYB <- sapply(c(Hc = "Hc", He = "He"), function(k) { m <- if (k == "Hc") HYB.SC else HYB.SE; eta <- if (k == "Hc") SPLIT.ETA$Hc else SPLIT.ETA$He
  e <- colSums(m) / DEPTH.REF
  unlist(mclapply(M3.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ]), e, eta, "NB")[["disp"]], mc.cores = CORES)) })
rownames(M3.FIT.HYB) <- M3.GENE.IDS
M3.FIT.PAR.FULL <- sapply(c(Pc = "Pc", Pe = "Pe"), function(k) { m <- if (k == "Pc") MIX.SC else MIX.SE; e <- colSums(m) / DEPTH.REF
  unlist(mclapply(M3.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ]), e, 0, "NB")[["disp"]], mc.cores = CORES)) })
rownames(M3.FIT.PAR.FULL) <- M3.GENE.IDS

## Thinned parental fits, one matrix per draw. Each draw thins the whole matrix (every gene's UMIs, so the libraries and offsets thin with them) on its own stream.
M3.FIT.PAR.THIN <- lapply(seq_len(MID3.REPS), function(r) with_local_seed(SEED.MID3 + r, {
  out <- sapply(c(Pc = "Pc", Pe = "Pe"), function(k) {
    m <- if (k == "Pc") MIX.SC else MIX.SE
    if (inherits(m, "sparseMatrix")) { m <- as(m, "CsparseMatrix"); m@x <- as.numeric(rbinom(length(m@x), as.integer(m@x), M3.P[[k]])); m <- drop0(m) }
    else m[] <- rbinom(length(m), as.integer(m), M3.P[[k]])
    e <- colSums(m) / DEPTH.REF
    unlist(mclapply(M3.GENE.IDS, function(g) fit_dispersion(as.numeric(m[g, ]), e, 0, "NB")[["disp"]], mc.cores = CORES))
  })
  rownames(out) <- M3.GENE.IDS; out
}))

## Per-gene d, full and thinned parents, on both scales
M3.MEAN.COUNT <- 0.5 * (rowMeans(MIX.SC[M3.GENE.IDS, , drop = FALSE]) + rowMeans(MIX.SE[M3.GENE.IDS, , drop = FALSE]))
M3.D <- lapply(c(list(full = M3.FIT.PAR.FULL), setNames(M3.FIT.PAR.THIN, paste0("thin", seq_len(MID3.REPS)))), function(par) {
  ph <- cbind(par, M3.FIT.HYB); l <- suppressWarnings(log2(ph)); l[!is.finite(l)] <- NA
  data.frame(d_log2 = 0.5 * (l[, "Hc"] + l[, "He"]) - 0.5 * (l[, "Pc"] + l[, "Pe"]),
             d_raw = 0.5 * (ph[, "Hc"] + ph[, "He"]) - 0.5 * (ph[, "Pc"] + ph[, "Pe"]), row.names = M3.GENE.IDS)
})

## Per group: share of genes at phi = 0 and median log2 phi of the positive fits
MID3.GROUPS <- do.call(rbind, c(
  lapply(c("Hc", "He"), function(k) { v <- M3.FIT.HYB[, k]; data.frame(group = k, version = "full", share_zero = mean(v == 0, na.rm = TRUE),
                                                                            median_log2_phi_pos = median(log2(v[v > 0]), na.rm = TRUE), row.names = NULL) }),
  lapply(c("Pc", "Pe"), function(k) { v <- M3.FIT.PAR.FULL[, k]; data.frame(group = k, version = "full", share_zero = mean(v == 0, na.rm = TRUE),
                                                                                median_log2_phi_pos = median(log2(v[v > 0]), na.rm = TRUE), row.names = NULL) }),
  unlist(lapply(seq_len(MID3.REPS), function(r) lapply(c("Pc", "Pe"), function(k) { v <- M3.FIT.PAR.THIN[[r]][, k]
    data.frame(group = k, version = paste0("thin", r), share_zero = mean(v == 0, na.rm = TRUE), median_log2_phi_pos = median(log2(v[v > 0]), na.rm = TRUE), row.names = NULL) })),
    recursive = FALSE)))
cat("\nPer group: share of genes at phi = 0 and median log2 phi of the positive fits:\n"); print(MID3.GROUPS, digits = 3, row.names = FALSE)

## d with full and thinned parents: on the genes where each is defined, on the genes where both are, and by decile of the full mean count
MID3.D <- do.call(rbind, lapply(seq_len(MID3.REPS), function(r) {
  f <- M3.D$full; t <- M3.D[[paste0("thin", r)]]; both <- is.finite(f$d_log2) & is.finite(t$d_log2)
  dc <- cut(M3.MEAN.COUNT, unique(quantile(M3.MEAN.COUNT, 0:10 / 10)), include.lowest = TRUE)
  rbind(data.frame(draw = r, bin = "all", genes_full = sum(is.finite(f$d_log2)), genes_thin = sum(is.finite(t$d_log2)), genes_both = sum(both),
                   median_d_log2_full = median(f$d_log2, na.rm = TRUE), median_d_log2_thin = median(t$d_log2, na.rm = TRUE),
                   median_d_log2_full_both = median(f$d_log2[both]), median_d_log2_thin_both = median(t$d_log2[both]),
                   median_d_raw_full = median(f$d_raw, na.rm = TRUE), median_d_raw_thin = median(t$d_raw, na.rm = TRUE), row.names = NULL),
        do.call(rbind, lapply(levels(dc), function(b) { i <- dc %in% b; bi <- i & both
          data.frame(draw = r, bin = b, genes_full = sum(i & is.finite(f$d_log2)), genes_thin = sum(i & is.finite(t$d_log2)), genes_both = sum(bi),
                     median_d_log2_full = median(f$d_log2[i], na.rm = TRUE), median_d_log2_thin = median(t$d_log2[i], na.rm = TRUE),
                     median_d_log2_full_both = median(f$d_log2[bi]), median_d_log2_thin_both = median(t$d_log2[bi]),
                     median_d_raw_full = median(f$d_raw[i], na.rm = TRUE), median_d_raw_thin = median(t$d_raw[i], na.rm = TRUE), row.names = NULL) })))
}))
cat("\nMidparent d of the dispersion with full and with thinned parents (draw 1 by decile of the full mean count; every draw overall):\n")
print(MID3.D[MID3.D$draw == 1, ], digits = 3, row.names = FALSE)
print(MID3.D[MID3.D$bin == "all", ], digits = 3, row.names = FALSE)

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(MID3.D, file.path(CHECK.DIR, "MID3_thin.csv"), row.names = FALSE)
write.csv(MID3.GROUPS, file.path(CHECK.DIR, "MID3_thin_groups.csv"), row.names = FALSE)
MID3.ALL <- MID3.D[MID3.D$bin == "all", ]
check_save("MID3", list(d = MID3.D, groups = MID3.GROUPS, p = M3.P, hyb = M3.FIT.HYB, par_full = M3.FIT.PAR.FULL, par_thin = M3.FIT.PAR.THIN,
                        settings = list(reps = MID3.REPS, min_count = MID3.MIN.COUNT)), "run",
           headline = c(sprintf("thinning Pc %.3f, Pe %.3f", M3.P[["Pc"]], M3.P[["Pe"]]),
                        sprintf("draw %d: log2 d full %.3f, thinned %.3f (genes with both: %.3f -> %.3f, %d genes); raw d full %.4f, thinned %.4f",
                                MID3.ALL$draw, MID3.ALL$median_d_log2_full, MID3.ALL$median_d_log2_thin, MID3.ALL$median_d_log2_full_both,
                                MID3.ALL$median_d_log2_thin_both, MID3.ALL$genes_both, MID3.ALL$median_d_raw_full, MID3.ALL$median_d_raw_thin)),
           mc_error = sprintf("%d thinning draws", MID3.REPS), script = CHECK.SCRIPT)
