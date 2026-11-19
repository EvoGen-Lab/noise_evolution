###############################################################
### gene_mu_phi_table.R
### One row per gene: its ORF, its gene name and the point estimates of the mean rate mu and the dispersion phi in each group (Pc, Pe, Hc, He), from
### the Section 1 checkpoint alone. The fits are the pipeline's (fit_dispersion(), Cox-Reid NB) on all cells of each group, with the offsets of
### analysis.R Section 2: a parent's library over the median parental library, and half the cell's hybrid depth for each hybrid allele. The allele
### dispersions carry the transcriptome share corrections (R.C, R.E) and the allele means are divided by their abundance factors 2 c-bar and
### 2 (1 - c-bar), so every mu is a per-genome rate at the reference depth. phi = 0 is a fit at the Poisson boundary; NA is no counts or a failed fit.
### Writes results/tables/gene_mu_phi_checkpoint1.csv (TABLE.DIR).
###############################################################
library(here); source(here("R", "setup.R"))
library(parallel); library(org.Sc.sgd.db)

CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core

stopifnot("gene_mu_phi_table.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE (cells after QC, genes after the union filter), C.BAR, R.C, R.E

## Offsets and share corrections of Section 2
DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
TAB.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
TAB.EXPOS <- list(Pc = colSums(MIX.SC) / DEPTH.REF, Pe = colSums(MIX.SE) / DEPTH.REF,
                  Hc = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF / 2, He = (colSums(HYB.SC) + colSums(HYB.SE)) / DEPTH.REF / 2)
TAB.ETA   <- c(Pc = 0, Pe = 0, Hc = R.C, He = R.E)
TAB.ABUND <- c(Pc = 1, Pe = 1, Hc = 2 * C.BAR, He = 2 * (1 - C.BAR))     # abundance factors: a hybrid allele mean divided by it is a per-genome rate
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)))

GENE.FIT.TABLE <- data.frame(orf = rownames(MIX.SC), stringsAsFactors = FALSE)
GENE.FIT.TABLE$gene <- unname(suppressMessages(AnnotationDbi::mapIds(org.Sc.sgd.db, keys = GENE.FIT.TABLE$orf, keytype = "ORF", column = "GENENAME", multiVals = "first")))
GENE.FIT.TABLE$gene <- ifelse(is.na(GENE.FIT.TABLE$gene), GENE.FIT.TABLE$orf, GENE.FIT.TABLE$gene)   # ORFs without a standard name keep their ORF
for (k in names(TAB.MATS)) {
  f <- do.call(rbind, mclapply(seq_len(nrow(TAB.MATS[[k]])), function(i) fit_dispersion(TAB.MATS[[k]][i, ], TAB.EXPOS[[k]], TAB.ETA[[k]], "NB")[c("mu", "disp")],
                               mc.cores = CORES))
  GENE.FIT.TABLE[[paste0(k, "_mu")]]  <- f[, "mu"] / TAB.ABUND[[k]]
  GENE.FIT.TABLE[[paste0(k, "_phi")]] <- f[, "disp"]
  cat(sprintf("%s: %d genes fitted, %d at phi = 0, %d without an estimate\n", k, sum(is.finite(f[, "disp"])), sum(f[, "disp"] %in% 0), sum(!is.finite(f[, "disp"]))))
}

dir.create(TABLE.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(GENE.FIT.TABLE, file.path(TABLE.DIR, "gene_mu_phi_checkpoint1.csv"), row.names = FALSE)
print(head(GENE.FIT.TABLE), digits = 3, row.names = FALSE)
