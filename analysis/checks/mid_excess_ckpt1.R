###############################################################
### mid_excess_ckpt1.R
### MID2: is the shared rise of the hybrid's dispersion over the midparent (MID1: median log2 d of about 1.4 for phi, falling with expression) a variance
### term linear in the mean, which the NB lacks, or a true rise of phi? From the Section 1 checkpoint alone.
### Model per gene and group: Var(Y_i) = (1 + a) M_i + phi M_i^2, M_i = O_i mu (linear_excess_gene(), the L6 regression of the leak tests). A term a M_i
### with a > 0 adds variance that grows like the mean, the signature of counting (allele assignment, UMI handling); the NB fit then reads it as a
### dispersion of about phi + a / M, which falls with expression. A culture or other biological difference instead moves phi itself, at every
### expression level. So:
###   1. a per group: the median a and the share of genes with a > 2 SE, for Pc, Pe, Hc, He. Counting in the hybrid shows as a larger a in Hc and He.
###   2. d from the quadratic part: d_lin = (log2 phi_lin,Hc + log2 phi_lin,He) / 2 - (log2 phi_lin,Pc + log2 phi_lin,Pe) / 2, on the genes with every
###      phi_lin > 0, against the Cox-Reid d (MID1's) on the same genes, overall and by decile of mean count. If d_lin falls to near zero and is flat, the
###      linear term explains the shift; if d_lin stays near the Cox-Reid d, phi itself is higher in the hybrid.
###   3. the raw-scale version, phi_lin midpoints as arithmetic means (no positivity needed), and the share of genes with a negative phi_lin per group.
### Offsets as in Section 2 (parents: library over the median parental library; hybrid alleles: the library of their own genome over the same reference); the
### allele dispersions carry the library split correction (split variances SPLIT.ETA), and phi_lin has the split bias delta of the gene's Cox-Reid fit subtracted. Genes need a mean count of at least MID2.MIN.COUNT per cell in every group.
### Writes results/checks/MID2_excess.csv (per group and per decile summaries), MID2_genes.csv, MID2.rds and a row in summary.csv.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(parallel)

CHECK.SCRIPT <- "analysis/checks/mid_excess_ckpt1.R"
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
if (!exists("MID2.MIN.COUNT")) MID2.MIN.COUNT <- 0.1      # minimum mean count per cell in every group, as MID1

stopifnot("mid_excess_ckpt1.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
M2.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
M2.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF, Pe = colSums(MIX.SE) / DEPTH.REF, Hc = colSums(HYB.SC) / DEPTH.REF, He = colSums(HYB.SE) / DEPTH.REF)
M2.ETA   <- list(Pc = 0, Pe = 0, Hc = SPLIT.ETA$Hc, He = SPLIT.ETA$He)
M2.GENE.IDS <- rownames(MIX.SC)[Reduce(`&`, lapply(M2.MATS, function(m) rowMeans(m) >= MID2.MIN.COUNT))]
cat(sprintf("%d genes with a mean count of at least %g per cell in every group\n", length(M2.GENE.IDS), MID2.MIN.COUNT))

## Per gene and group: the linear coefficient a with its SE, phi_lin and the Cox-Reid phi (both split corrected) and the mean count per cell
M2.FIT <- lapply(setNames(names(M2.MATS), names(M2.MATS)), function(k) {
  le <- do.call(rbind, mclapply(M2.GENE.IDS, linear_excess_gene, mat = M2.MATS[[k]], expo = M2.EXPOS[[k]], mc.cores = CORES))
  cr <- t(vapply(M2.GENE.IDS, function(g) fit_dispersion(M2.MATS[[k]][g, ], M2.EXPOS[[k]], M2.ETA[[k]], "NB")[c("disp", "split_bias")], numeric(2)))
  data.frame(a = le[, "a"], a_se = le[, "a_se"], phi_lin = le[, "phi_lin"] - ifelse(is.finite(cr[, "split_bias"]), cr[, "split_bias"], 0), phi_cr = cr[, "disp"], mbar = le[, "mbar"],
             row.names = M2.GENE.IDS)
})

## 1. The linear term per group
MID2.GROUPS <- do.call(rbind, lapply(names(M2.FIT), function(k) { f <- M2.FIT[[k]]; ok <- is.finite(f$a) & is.finite(f$a_se) & f$a_se > 0
  data.frame(group = k, genes = sum(ok), median_a = median(f$a[ok]), share_a_above_2se = mean(f$a[ok] / f$a_se[ok] > 2), share_a_below_m2se = mean(f$a[ok] / f$a_se[ok] < -2),
             median_phi_cr = median(f$phi_cr, na.rm = TRUE), median_phi_lin = median(f$phi_lin, na.rm = TRUE), share_phi_lin_negative = mean(f$phi_lin < 0, na.rm = TRUE),
             median_mbar = median(f$mbar, na.rm = TRUE), row.names = NULL) }))
cat("\nLinear variance term a and dispersions per group (Var = (1 + a) M + phi M^2):\n"); print(MID2.GROUPS, digits = 3, row.names = FALSE)

## 2 and 3. d from the Cox-Reid phi and from phi_lin, per gene
MID2.GENES <- local({
  P <- lapply(c("phi_cr", "phi_lin"), function(v) sapply(M2.FIT, function(f) f[[v]])); names(P) <- c("cr", "lin")
  l2 <- lapply(P, function(m) { x <- suppressWarnings(log2(m)); x[!is.finite(x)] <- NA; x })
  geo <- function(l) 0.5 * (l[, "Hc"] + l[, "He"]) - 0.5 * (l[, "Pc"] + l[, "Pe"])
  data.frame(gene = M2.GENE.IDS, mean_count = rowMeans(sapply(M2.MATS[c("Pc", "Pe")], function(m) rowMeans(m[M2.GENE.IDS, , drop = FALSE]))),
             d_cr = geo(l2$cr), d_lin = geo(l2$lin),
             d_cr_raw = 0.5 * (P$cr[, "Hc"] + P$cr[, "He"]) - 0.5 * (P$cr[, "Pc"] + P$cr[, "Pe"]),
             d_lin_raw = 0.5 * (P$lin[, "Hc"] + P$lin[, "He"]) - 0.5 * (P$lin[, "Pc"] + P$lin[, "Pe"]),
             a_hyb = 0.5 * (M2.FIT$Hc$a + M2.FIT$He$a), a_par = 0.5 * (M2.FIT$Pc$a + M2.FIT$Pe$a), stringsAsFactors = FALSE)
})
MID2.BOTH <- is.finite(MID2.GENES$d_cr) & is.finite(MID2.GENES$d_lin)     # every phi_cr and phi_lin positive
MID2.DECILE <- local({
  x <- MID2.GENES[MID2.BOTH, ]; dc <- cut(x$mean_count, unique(quantile(x$mean_count, 0:10 / 10)), include.lowest = TRUE)
  out <- do.call(rbind, lapply(levels(dc), function(b) { i <- dc %in% b
    data.frame(mean_count_bin = b, genes = sum(i), median_d_cr = median(x$d_cr[i]), median_d_lin = median(x$d_lin[i]),
               median_a_hyb = median(x$a_hyb[i], na.rm = TRUE), median_a_par = median(x$a_par[i], na.rm = TRUE), row.names = NULL) }))
  rbind(data.frame(mean_count_bin = "all", genes = nrow(x), median_d_cr = median(x$d_cr), median_d_lin = median(x$d_lin),
                   median_a_hyb = median(x$a_hyb, na.rm = TRUE), median_a_par = median(x$a_par, na.rm = TRUE)), out)
})
cat(sprintf("\nlog2 d of the dispersion, Cox-Reid against the quadratic part phi_lin, on the %d genes with every phi positive (of %d):\n", sum(MID2.BOTH), nrow(MID2.GENES)))
print(MID2.DECILE, digits = 3, row.names = FALSE)
MID2.RAW <- data.frame(genes = sum(is.finite(MID2.GENES$d_cr_raw) & is.finite(MID2.GENES$d_lin_raw)),
                       median_d_cr_raw = median(MID2.GENES$d_cr_raw, na.rm = TRUE), median_d_lin_raw = median(MID2.GENES$d_lin_raw, na.rm = TRUE))
cat("\nRaw scale (arithmetic midpoints, every gene with finite values):\n"); print(MID2.RAW, digits = 3, row.names = FALSE)

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(rbind(cbind(table = "groups", MID2.GROUPS[, c("group", "genes", "median_a", "share_a_above_2se")], mean_count_bin = NA, median_d_cr = NA, median_d_lin = NA),
                cbind(table = "deciles", group = NA, genes = MID2.DECILE$genes, median_a = MID2.DECILE$median_a_hyb, share_a_above_2se = NA,
                      mean_count_bin = MID2.DECILE$mean_count_bin, median_d_cr = MID2.DECILE$median_d_cr, median_d_lin = MID2.DECILE$median_d_lin)),
          file.path(CHECK.DIR, "MID2_excess.csv"), row.names = FALSE)
write.csv(MID2.GENES, file.path(CHECK.DIR, "MID2_genes.csv"), row.names = FALSE)
check_save("MID2", list(groups = MID2.GROUPS, deciles = MID2.DECILE, raw = MID2.RAW, genes = MID2.GENES, fits = M2.FIT, settings = list(min_count = MID2.MIN.COUNT)),
           "run", headline = c(sprintf("%s: median a %.3f, a > 2 SE in %.2f of genes, phi_lin < 0 in %.2f", MID2.GROUPS$group, MID2.GROUPS$median_a,
                                       MID2.GROUPS$share_a_above_2se, MID2.GROUPS$share_phi_lin_negative),
                               sprintf("log2 d, %d genes: Cox-Reid %.3f, quadratic part %.3f; raw d: Cox-Reid %.4f, quadratic part %.4f", MID2.DECILE$genes[1],
                                       MID2.DECILE$median_d_cr[1], MID2.DECILE$median_d_lin[1], MID2.RAW$median_d_cr_raw, MID2.RAW$median_d_lin_raw)),
           mc_error = "point estimates; a with its sandwich SE per gene", script = CHECK.SCRIPT)
