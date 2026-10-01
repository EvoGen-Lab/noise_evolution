##############################################################################
## SINGLE-CELL DATA ANALYSIS                                                ##
##############################################################################
###   1. DATA PROCESSING
###     1.1 Data Import
###   2. NEGATIVE BINOMIAL FITS
###     2.1 Internal split fraction
###     2.2 Split-dependent NB fits
###     2.3 Per-gene contrasts, bootstrap SEs, and permutation null
###     2.4 Bootstrap adequacy check
###   3. REGULATORY AND DOMINANCE RESULTS
###     3.1 Regulatory classification
###     3.2 Significance histograms for total, cis, and trans effects
###     3.3 Dominance classification
###     3.4 Mean against burst frequency and burst size, parents and hybrids
###     3.5 Class relationship heatmaps
###     3.6 Gene-identity overlap
###   4. CO-EXPRESSION
###     4.1 Gene set, residuals, and point-estimate correlation matrices
###     4.2 Bootstrap for coexpression
###     4.3 Bootstrap reliability check
###     4.4 Bootstrap adequacy check
###     4.5 Pair-level classification
###     4.6 Co-expression permutation null
###     4.7 Rank check and candidate axes
###     4.8 Mixture model and GO enrichment for candidate axes
###     4.9 Permutation validation
###     4.10 Total Axis 1
###     4.11 Total Axis 2
###     4.12 Burst frequency/size consistency among co-expressed genes
###     4.13 Cis and trans candidate axes
###   5. INTRINSIC / EXTRINSIC NOISE
###     5.1 Reliability calibration
###     5.2 Depth confound check
###     5.3 Per-gene fractions
###     5.4 Regulatory and dominance class summaries
###     5.5 Burst kinetics relationship 
###   6. PROMOTER ARCHITECTURE
###     6.1 Genomes, annotations, and per-species promoter scores
###     6.2 Validation of individual features
###     6.3 Between-species divergence in architecture
###     6.4 Directional concordance of features
###     6.5 Candidate genes
###   7. BROAD CELLULAR DIFFERENCES IN EXPRESSION AND REGULATION
###     7.1 Build Seurat objects
###     7.2 Normalize, select features, and reduce dimensionality
###     7.3 Cluster each dataset
###     7.4 Cluster composition: identity, marker enrichment, and consistency
###     7.5 Within/between-cluster noise partitioning
###     7.6 Cell-cycle and metabolic module scoring
###     7.7 Covariate-noise diagnostic
###     7.8 Species composition comparison and confound bound
###   8. GO ENRICHMENT
###     8.1 Overall parental divergence
###     8.2 Regulatory classes, split by direction
###     8.3 Any-cis / any-trans, pooled across class
###     8.4 Enrichment for each set
###     8.5 Summary table
###     8.6 Gene set enrichment
###     8.7 Intrinsic / extrinsic noise
###   9. EXTERNAL NOISE VALIDATION
###     9.1 Load external protein datasets
###     9.2 NB-implied noise and per-source merge
###     9.3 Published single-cell RNA-seq sources
###     9.4 Correlations among data sets
###     9.5 Mean-adjusted noise
###     9.6 Diagnostic plots
###   10. POWER ANALYSIS
###     10.1 Grid definition and cluster submission
###     10.2 Reshape cluster output into the POWER array
###     10.3 Figure 1: power vs mean expression, panels over N.CELLS x SIZE.RATIO
###     10.4 Supplementary diagnostic grids
###     10.5 Figure 2: power vs burst frequency, panels over N.CELLS x SIZE.RATIO
###     10.6 Figure 3: power vs mean expression, panels over N.CELLS x SIZE
###     10.7 Heatmap A
###     10.8 Heatmap B
###
### Main Figure Locations
### Figure 1: Section 3.1        Figure 5: Section 3.4
### Figure 2: Section 3.1        Figure 6: Section 3.4
### Figure 3: Section 3.4        Figure 7: Section 3.4
### Figure 4: Section 3.4        Figure 8: Section 5.4
### Figure 11: Section 4.7
### Figure numbers 9 and 10 are held for the power-analysis figures
### drawn in Section 10.
###############################################################

# renv pins every package to the version recorded in renv.lock. Run
# renv::restore() once on a new machine
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::status()

# here() anchors every path to the project root so that functions and data load correctly
library(here)
stopifnot(file.exists(here("R", "functions.R")))
FUNCTIONS.FILE <- here("R", "functions.R")
source(FUNCTIONS.FILE)

# data/ holds everything read from disk. cluster_inputs and cluster_outputs
# carry *_inputs.rda and *_output.rda to and from the SLURM cluster.
SC.DIR       <- here("data", "single_cell")
GENOME.DIR   <- here("data", "genomes")
EXTERNAL.DIR <- here("data", "external")
INPUT.DIR    <- here("data", "cluster_inputs")
OUTPUT.DIR   <- here("data", "cluster_outputs")

# results/ holds everything this script writes. checkpoints hold the
# section{N}_checkpoint.rda files, console holds a copy of what each section
# prints to the R console (section{N}_console.txt), and tables holds .csv files.
CHECKPOINT.DIR <- here("results", "checkpoints")
CONSOLE.DIR    <- here("results", "console")
TABLE.DIR      <- here("results", "tables")

# figures/ holds main text figures, extended data/supplement figures, and everything
# else (diagnostic and exploratory plots)
FIGURE.DIR <- here("figures")
for (d in c(INPUT.DIR, OUTPUT.DIR, CHECKPOINT.DIR, CONSOLE.DIR, TABLE.DIR,
            file.path(FIGURE.DIR, c("main", "extended", "extra")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Packages
suppressWarnings(suppressPackageStartupMessages({
  library(MASS)
  library(parallel)
  library(future)
  library(ggplot2)
  library(Seurat)
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(enrichplot)
  library(mixtools)
  library(Biostrings)
  library(rtracklayer)
  library(cluster)
  library(mclust)
  library(readxl)
  library(data.table)
}))

# Classes
REG.CLASS <- c("Conserved", "Cis", "Trans", "Cis + Trans", "Compensatory")
DOM.CLASS <- c("Conserved", "Sc.Dominant", "Se.Dominant", "Overdominant", "Underdominant", "Additive")

## Color palettes
# Species anchors
SPECIES.COLOR <- c(Sc = "#B5533C", Se = "#3E6B7A")
# Regulatory classes
COLOR.LIST.1 <- c("#D4D4CF","#2B2F42","#640B14","#E0A526","#5F8F5A"); names(COLOR.LIST.1) <- REG.CLASS
# Dominance classes
COLOR.LIST.2 <- c("#D4D4CF", SPECIES.COLOR[["Sc"]], SPECIES.COLOR[["Se"]], "#7A2E4E","#D9A7BF","#8C8C86"); names(COLOR.LIST.2) <- DOM.CLASS
# Diverging heatmap ramp (warm = enriched, cool = depleted
COLOR.LIST.3 <- colorRampPalette(c("#70453F","#A2623A","#C8893A","#DDB264","#F7F7F7","#D6ABD6","#BB82C3","#9D57AA","#6F3C79"))(11)
# Warm ramp with monotone lightness
SEQ.ANCHORS <- c("#F6F1E4","#E8C77A","#C98A3A","#8A4A3A","#33384A")
COLOR.SEQ <- colorRampPalette(SEQ.ANCHORS)(100)
# Cell-cycle phases
COLOR.PHASE <- c(G1 = "#A9CBF2", S = "#3F63BF", G2M = "#132057")
# Power-curve lines follow the same anchors as COLOR.SEQ, minus the cream
# end that vanishes on white
POWER.COLOR <- SEQ.ANCHORS[-1]
# Plum ramp for cluster labels
COLOR.CLUSTER <- c("#D7ACD7", "#A159AF", "#512B59")
# Neutral colors for thresholds and reference lines
COLOR.ACCENT <- "#3A3A3A"
COLOR.GREY <- c(light = "#DADAD5", mid = "#B0B0AB", dark = "#6E6E6A")

console_start(1)
##############################################################################
## 1. DATA PROCESSING                                                       ##
##############################################################################
## 1.1 Data
# MIX is a mixture of S.eubayanus (SE) and S.cerevisiae (SC)
# HYB is a hybrid of S.eubayanus (SE) x S.cerevisiae (SC)

# Read in count data matrices
MIX.SC <- as.matrix(read.csv(file = file.path(SC.DIR, "MixSc_for_genes_percell_pergene_count.csv"), row.names = 1))
MIX.SE <- as.matrix(read.csv(file = file.path(SC.DIR, "MixSe_for_genes_percell_pergene_count.csv"), row.names = 1))
HYB.SC <- as.matrix(read.csv(file = file.path(SC.DIR, "HybSc_for_genes_percell_pergene_count.csv"), row.names = 1))
HYB.SE <- as.matrix(read.csv(file = file.path(SC.DIR, "HybSe_for_genes_percell_pergene_count.csv"), row.names = 1))

# Remove cells that have substantial reads from both species in the mixed
# parent sample. For cells with only a small fraction of reads from the
# minority species, remove from the minority species only
MIX.SHARED.CELLS <- intersect(colnames(MIX.SC), colnames(MIX.SE))
MIX.SC.SHARED.SUM <- colSums(MIX.SC[, MIX.SHARED.CELLS, drop = FALSE])
MIX.SE.SHARED.SUM <- colSums(MIX.SE[, MIX.SHARED.CELLS, drop = FALSE])
MIX.MINORITY.FRAC <- pmin(MIX.SC.SHARED.SUM, MIX.SE.SHARED.SUM) / (MIX.SC.SHARED.SUM + MIX.SE.SHARED.SUM)

# Sweep of minority-read-fraction thresholds
MIX.FRAC.SWEEP <- 10^seq(log10(0.001), log10(0.5), length.out = 200)
MINORITY.FRAC.THRESHOLD <- 0.005
MIX.FRAC.SWEEP.N.EXCEED <- sapply(MIX.FRAC.SWEEP, function(t) sum(MIX.MINORITY.FRAC > t))
fig_pdf("extra/S_mix_shared_cell_contamination_sweep.pdf", 6, 5)
plot(MIX.FRAC.SWEEP * 100, MIX.FRAC.SWEEP.N.EXCEED, log = "x", type = "l",
     xlab = "Minority-species read fraction threshold (%, log scale)",
     ylab = "# of shared cells exceeding threshold",
     main = "Contamination fraction among cells shared between MIX.SC and MIX.SE")
rect(MINORITY.FRAC.THRESHOLD * 100, par("usr")[3], 10^par("usr")[2], par("usr")[4],
     col = adjustcolor(COLOR.ACCENT, alpha.f = 0.1), border = NA)
abline(v = MINORITY.FRAC.THRESHOLD * 100, col = COLOR.ACCENT, lwd = 2, lty = 2)
legend("topright",
       legend = sprintf("cutoff = %.1f%% (%d cells removed as doublets)",
                         MINORITY.FRAC.THRESHOLD * 100, sum(MIX.MINORITY.FRAC > MINORITY.FRAC.THRESHOLD)),
       col = COLOR.ACCENT, lwd = 2, lty = 2, bty = "n")
dev.off()

# Cells with substantial reads from both species: drop from both matrices
MIX.DOUBLET.CELLS <- MIX.SHARED.CELLS[MIX.MINORITY.FRAC > MINORITY.FRAC.THRESHOLD]
# Cells with only minor cross-contamination: drop from the minority species only
MIX.MINOR.CONTAM.CELLS <- MIX.SHARED.CELLS[MIX.MINORITY.FRAC <= MINORITY.FRAC.THRESHOLD]
MIX.MINOR.CONTAM.IN.SC <- MIX.MINOR.CONTAM.CELLS[MIX.SC.SHARED.SUM[MIX.MINOR.CONTAM.CELLS] < MIX.SE.SHARED.SUM[MIX.MINOR.CONTAM.CELLS]]
MIX.MINOR.CONTAM.IN.SE <- MIX.MINOR.CONTAM.CELLS[MIX.SE.SHARED.SUM[MIX.MINOR.CONTAM.CELLS] < MIX.SC.SHARED.SUM[MIX.MINOR.CONTAM.CELLS]]

MIX.SC <- MIX.SC[, !(colnames(MIX.SC) %in% c(MIX.DOUBLET.CELLS, MIX.MINOR.CONTAM.IN.SC))]
MIX.SE <- MIX.SE[, !(colnames(MIX.SE) %in% c(MIX.DOUBLET.CELLS, MIX.MINOR.CONTAM.IN.SE))]

# Combine mitochondrial directions
HYB.SC['MT_F',] <- HYB.SC['MT_F',]+HYB.SC['MT_R',]; HYB.SC <- HYB.SC[!(row.names(HYB.SC) %in% c('MT_R')),]
HYB.SE['MT_F',] <- HYB.SE['MT_F',]+HYB.SE['MT_R',]; HYB.SE <- HYB.SE[!(row.names(HYB.SE) %in% c('MT_R')),]
MIX.SC['MT_F',] <- MIX.SC['MT_F',]+MIX.SC['MT_R',]; MIX.SC <- MIX.SC[!(row.names(MIX.SC) %in% c('MT_R')),]
MIX.SE['MT_F',] <- MIX.SE['MT_F',]+MIX.SE['MT_R',]; MIX.SE <- MIX.SE[!(row.names(MIX.SE) %in% c('MT_R')),]

# Mitochondrial ratios in hybrid
# MITO.IGNORE flags hybrid cells with zero Sc or Se mitochondrial reads.
# The two histograms plot the per-cell log2(Sc/Se) mitochondrial read ratio, 
# and total Sc read depth split by whether a cell was flagged in MITO.IGNORE (red).
HYB.MITO.CELLS <- intersect(colnames(HYB.SC), colnames(HYB.SE))
HYB.SC.MITO <- HYB.SC['MT_F', HYB.MITO.CELLS]
HYB.SE.MITO <- HYB.SE['MT_F', HYB.MITO.CELLS]
MITO.IGNORE <- which(HYB.SC.MITO == 0 | HYB.SE.MITO == 0)
MITO.RATIO  <- log2(HYB.SC.MITO[-MITO.IGNORE]/HYB.SE.MITO[-MITO.IGNORE])
fig_pdf("extra/S_mito_ratio_hybrid.pdf", 6, 5)
par(mfrow=c(1,1))
hist(MITO.RATIO,breaks=40,xlab="log2(Sc/Se) mitochondrial ratio in hybrid",main="")
abline(v=mean(MITO.RATIO),col=COLOR.ACCENT,lwd=2,lty=2)
dev.off()
cat("Mitochondrial log2(Sc/Se) ratio in hybrid, quantiles:\n")
print(round(quantile(MITO.RATIO, probs = c(0, 0.25, 0.5, 0.75, 1)), 3))

fig_pdf("extra/S_mito_ignore_depth.pdf", 6, 5)
hist(colSums(HYB.SC[,HYB.MITO.CELLS[-MITO.IGNORE]]),breaks=seq(0,80000,2500), xlab="Hybrid Sc Reads/cell",main="")
hist(colSums(HYB.SC[,HYB.MITO.CELLS[MITO.IGNORE]]),breaks=seq(0,80000,2500),col=COLOR.GREY[["dark"]],add=TRUE)
dev.off()

# Remove mitochondrial reads from remainder of analyses.
MIX.SC <- MIX.SC[rownames(MIX.SC) != "MT_F", , drop = FALSE]
MIX.SE <- MIX.SE[rownames(MIX.SE) != "MT_F", , drop = FALSE]
HYB.SC <- HYB.SC[rownames(HYB.SC) != "MT_F", , drop = FALSE]
HYB.SE <- HYB.SE[rownames(HYB.SE) != "MT_F", , drop = FALSE]

# Remove low-quality cells and low-information genes
# Thresholds is based on each dataset's depth and cell number
CELL.MAD.K     <- 3     # keep cells above median - 3 MAD on the log10 scale
GENE.LAMBDA0   <- 0.10  # mean reads per cell required at the shallowest dataset
GENE.CELL.FRAC <- 0.10  # detected-cell floor as a fraction of the smallest dataset

# Hybrid alleles come from the same cells, so they share one QC decision
HYB.CELLS <- intersect(colnames(HYB.SC), colnames(HYB.SE))
HYB.SC <- HYB.SC[, HYB.CELLS]; HYB.SE <- HYB.SE[, HYB.CELLS]

## Cell QC from precomputed library size (lib). A vector input lets hybrid alleles
## be pooled so each cell receives one keep/drop decision. Only the lower tail is
## trimmed (median - k MAD on the log10 scale), which retains large G2 cells that
## carry more RNA. min_reads is a floor for empty barcodes that rarely binds in
## typical data. Returns list(keep = logical per cell, lib_cut = cutoff used).
# Cells are filtered first, on all genes, so depth reflects the full library
QC.STATS <- list(
  MIX.SC = list(lib = colSums(MIX.SC), det = colSums(MIX.SC > 0)),
  MIX.SE = list(lib = colSums(MIX.SE), det = colSums(MIX.SE > 0)),
  HYB    = list(lib = colSums(HYB.SC) + colSums(HYB.SE),
                det = colSums(HYB.SC > 0) + colSums(HYB.SE > 0)))
# Genes-detected floor was dropped entirely (too aggressive, especially on
# the already-small MIX.SC dataset); library size alone now decides
QC.CELLS <- lapply(QC.STATS, qc_cell_cutoff, k = CELL.MAD.K)
QC.TITLES <- c(MIX.SC = "Sc Parent", MIX.SE = "Se Parent", HYB = "Hybrid (alleles pooled)")
fig_pdf("extra/S_cell_count_threshold_qc.pdf", 12, 4)
par(mfrow=c(1,3))
for (nm in names(QC.STATS)) {
  plot(QC.STATS[[nm]]$lib, QC.STATS[[nm]]$det, log = "xy", pch=19, cex=0.6, xlab="Count per Cell", ylab="Number of Genes", main=QC.TITLES[[nm]])
  abline(v=QC.CELLS[[nm]]$lib_cut, col=COLOR.ACCENT, lty=2)
}
dev.off()
MIX.SC <- MIX.SC[, QC.CELLS$MIX.SC$keep]; MIX.SE <- MIX.SE[, QC.CELLS$MIX.SE$keep]
HYB.SC <- HYB.SC[, QC.CELLS$HYB$keep];    HYB.SE <- HYB.SE[, QC.CELLS$HYB$keep]

# Genes are filtered second, using the retained cells and one relative abundance floor
GENE.QC <- qc_gene_keep(list(MIX.SC = MIX.SC, MIX.SE = MIX.SE, HYB.SC = HYB.SC, HYB.SE = HYB.SE),
                        lambda0 = GENE.LAMBDA0, cell_frac = GENE.CELL.FRAC)
fig_pdf("extra/S_gene_count_threshold_qc.pdf", 8, 8)
par(mfrow=c(2,2))
GENE.QC.MATS <- list(MIX.SC = MIX.SC, MIX.SE = MIX.SE, HYB.SC = HYB.SC, HYB.SE = HYB.SE)
GENE.QC.TITLES <- c(MIX.SC = "Sc Parent", MIX.SE = "Se Parent", HYB.SC = "Sc Hybrid", HYB.SE = "Se Hybrid")
for (nm in names(GENE.QC.MATS)) {
  m <- GENE.QC.MATS[[nm]]
  plot(rowSums(m) + 1, rowSums(m > 0), log = "x", pch=19, cex=0.6, ylim=c(0,ncol(m)), xlab="Count per Gene (+1)", ylab="Number of Cells with >=1 Count", main=GENE.QC.TITLES[[nm]])
  abline(v=GENE.QC$p_min * sum(m) + 1, col=COLOR.ACCENT, lty=2); abline(h=GENE.QC$n_min, col=COLOR.ACCENT, lty=2)
}
dev.off()

# Subset to the same set of genes across all samples, in the same order
GENE.SET <- GENE.QC$genes
MIX.SC <- MIX.SC[GENE.SET,]; MIX.SE <- MIX.SE[GENE.SET,]
HYB.SC <- HYB.SC[GENE.SET,]; HYB.SE <- HYB.SE[GENE.SET,]

# Comparability check. Retained depth should differ by dataset while the
# implied mean count floor scales with depth.
QC.SUMMARY <- data.frame(
  n_cells  = c(ncol(MIX.SC), ncol(MIX.SE), ncol(HYB.SC), ncol(HYB.SE)),
  reads_per_cell = c(sum(MIX.SC)/ncol(MIX.SC), sum(MIX.SE)/ncol(MIX.SE), sum(HYB.SC)/ncol(HYB.SC), sum(HYB.SE)/ncol(HYB.SE)),
  mean_count_floor = unname(GENE.QC$exp_mean),
  row.names = c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE"))
cat("Cell and gene QC summary\n"); print(round(QC.SUMMARY, 3))
cat(sprintf("Genes retained: %d, minimum detected cells: %d\n", length(GENE.SET), GENE.QC$n_min))

# Checkpoint
save(MIX.SC, MIX.SE, HYB.SC, HYB.SE, QC.SUMMARY, file = ckpt_path(1))

console_start(2)
##############################################################################
## 2. NEGATIVE BINOMIAL FITS                                                ##
##############################################################################
# Order hybrid cells by total reads
HYB.ORD <- order(colSums(HYB.SC) + colSums(HYB.SE), decreasing = TRUE)

# Reorder the two hybrid alleles by depth, then sum them into a single
# pooled hybrid matrix for analyses that don't need the alleles kept separate
HYB.SC   <- HYB.SC[, HYB.ORD]
HYB.SE   <- HYB.SE[, HYB.ORD]
HYB.COMB <- HYB.SC + HYB.SE

# Per-cell relative depth, used as an offset in the negative binomial fits
DEPTH.REF <- median(c(colSums(MIX.SC), colSums(MIX.SE), colSums(HYB.SC) + colSums(HYB.SE)))
EXPO.MIX.SC <- colSums(MIX.SC)/DEPTH.REF
EXPO.MIX.SE <- colSums(MIX.SE)/DEPTH.REF
EXPO.HYB    <- (colSums(HYB.SC) + colSums(HYB.SE))/DEPTH.REF

# Seeds and resample counts for bootstrap, permutation, and coexpression steps
N.BOOT <- 1000; SEED.BOOT <- 1            # bootstrap resamples
N.PERM <- 10000; SEED.PERM <- 1          # permutation shuffles
N.COEXPR <- 5000; SEED.COEXPR <- 1        # Coexpression

## 2.1 Internal split fraction
# Builds a rough gene set from raw count for an exploratory NB fit. Used to 
# choose the split of cells that should be used so that cis and trans estimates
# have similar power to detect significant differences
PILOT.MATS  <- list(MIX.SC = MIX.SC, MIX.SE = MIX.SE, HYB.SC = HYB.SC, HYB.SE = HYB.SE)
PILOT.EXPOS <- list(MIX.SC = EXPO.MIX.SC, MIX.SE = EXPO.MIX.SE, HYB = EXPO.HYB)
PILOT.GENES <- qc_gene_keep(PILOT.MATS, lambda0 = 0.001, cell_frac = 0.10)$genes

## ---- Cluster round trip: Rscript gene_pilot.R ----
## Reads gene_pilot_inputs.rda (saved below), writes gene_pilot_output.rda (PILOT.SE)
save(PILOT.MATS, PILOT.EXPOS, PILOT.GENES, N.BOOT, SEED.BOOT, file = file.path(INPUT.DIR, "gene_pilot_inputs.rda"))
load(file.path(OUTPUT.DIR, "gene_pilot_output.rda"))
## ---- end cluster round trip ----

## Pooled f* for the mean and noise axes. A and B are pooled across
## genes BEFORE the ratio and sqrt, not averaged per-gene afterward,
## because squaring a bootstrap SE roughly doubles its relative
## error, and a per-gene ratio of two such squared, independently
## noisy quantities amplifies sampling noise, worse for size than for
## mean, since dispersion estimates are inherently noisier to begin
## with. Median pooling is used because both A and B are right-
## skewed: low expression inflates the mean axis, low true
## overdispersion inflates the noise axis. Neither is a reason to
## drop a gene, since both just mean less information for that gene's
## Sc/Se contrast, the same way a small true effect size makes
## detection harder without indicating a problem with the gene.
##
## A is the per-cell variance of the hybrid allele contrast, n_h times
## var(SC) + var(SE) - 2 cov(SC, SE): the two alleles are measured in the same
## cells, so their estimates co-vary and the contrast variance is smaller than
## the sum of the two marginal variances. B is the parental term; the parents
## are different cells, so its two variances add.
## pilot: PILOT.SE table from gene_pilot.R; n_h: number of hybrid cells.
## Returns f_mean / f_disp (split fractions), r_mean / r_disp (B * Nh / A), the
## number of genes pooled on each axis, and the median allele correlation
## (cor_mean / cor_disp) that the covariance term removes.
# Estimate split fraction for mean and noise
SPLIT.FRAC <- local({
  pilot <- PILOT.SE
  n_h <- ncol(HYB.SC)
  ## Closed-form balance point. With A the hybrid-driven variance
  ## coefficient and B the fixed, non-tunable parent-driven term,
  ## setting SE_cis(f) = SE_trans(f) gives a quadratic in f whose root
  ## in (0, 0.5] is f* = [(2+r) - sqrt(r^2+4)] / (2r), r = B*Nh/A.
  ## r -> 0 (parent term negligible) recovers f* -> 0.5, the even
  ## split, correctly, since there is then no asymmetry to correct for.

  cov_cols <- c("HYB_logmu_cov", "HYB_logdisp_cov")
  if (!all(cov_cols %in% names(pilot)))
    stop("PILOT.SE has no hybrid-allele covariance columns (", paste(cov_cols, collapse = ", "),
         "); rerun gene_pilot.R so it records them.")
  A_mean <- n_h * (pilot$HYB.SC_logmu_se^2   + pilot$HYB.SE_logmu_se^2   - 2 * pilot$HYB_logmu_cov)
  B_mean <-        pilot$MIX.SC_logmu_se^2   + pilot$MIX.SE_logmu_se^2
  A_disp <- n_h * (pilot$HYB.SC_logdisp_se^2 + pilot$HYB.SE_logdisp_se^2 - 2 * pilot$HYB_logdisp_cov)
  B_disp <-        pilot$MIX.SC_logdisp_se^2 + pilot$MIX.SE_logdisp_se^2

  ok_mean <- is.finite(A_mean) & is.finite(B_mean) & A_mean > 0
  ok_disp <- is.finite(A_disp) & is.finite(B_disp) & A_disp > 0

  r_mean <- median(B_mean[ok_mean]) * n_h / median(A_mean[ok_mean])
  r_disp <- median(B_disp[ok_disp]) * n_h / median(A_disp[ok_disp])

  list(f_mean = fstar_from_r(r_mean), f_disp = fstar_from_r(r_disp),
       r_mean = r_mean, r_disp = r_disp,
       n_genes_mean = sum(ok_mean), n_genes_disp = sum(ok_disp),
       cor_mean = median((pilot$HYB_logmu_cov   / (pilot$HYB.SC_logmu_se   * pilot$HYB.SE_logmu_se))[ok_mean],   na.rm = TRUE),
       cor_disp = median((pilot$HYB_logdisp_cov / (pilot$HYB.SC_logdisp_se * pilot$HYB.SE_logdisp_se))[ok_disp], na.rm = TRUE))
})

## 2.2 Split-dependent NB fits
# Partitions the hybrid cells: 
# HYC/HYT at SPLIT.FRAC$f_mean, used for the mean (MU) axis.
# HYC.N/HYT.N at SPLIT.FRAC$f_disp, used for the noise (DISP) axis.
SPLIT.IDX.MEAN <- split_indices_by_depth(ncol(HYB.SC), SPLIT.FRAC$f_mean)
SPLIT.IDX.DISP <- split_indices_by_depth(ncol(HYB.SC), SPLIT.FRAC$f_disp)

HYC.SC   <- HYB.SC[, SPLIT.IDX.MEAN$c]; HYT.SC   <- HYB.SC[, SPLIT.IDX.MEAN$t]
HYC.SE   <- HYB.SE[, SPLIT.IDX.MEAN$c]; HYT.SE   <- HYB.SE[, SPLIT.IDX.MEAN$t]
EXPO.HYC    <- EXPO.HYB[SPLIT.IDX.MEAN$c];    EXPO.HYT    <- EXPO.HYB[SPLIT.IDX.MEAN$t]

HYC.SC.N <- HYB.SC[, SPLIT.IDX.DISP$c]; HYT.SC.N <- HYB.SC[, SPLIT.IDX.DISP$t]
HYC.SE.N <- HYB.SE[, SPLIT.IDX.DISP$c]; HYT.SE.N <- HYB.SE[, SPLIT.IDX.DISP$t]
EXPO.HYC.N  <- EXPO.HYB[SPLIT.IDX.DISP$c];    EXPO.HYT.N  <- EXPO.HYB[SPLIT.IDX.DISP$t]

# SPLIT.FIT.MATS/SPLIT.FIT.EXPOS hold 13 datasets:
# 5 split-independent (MIX.SC, MIX.SE, HYB.SC, HYB.SE, HYB.COMB)
# 4 at f_mean (HYC/HYT x SC/SE)
# 4 at f_disp (HYC.N/HYT.N x SC/SE)
SPLIT.FIT.MATS <- list(
  MIX.SC = MIX.SC, MIX.SE = MIX.SE, HYB.SC = HYB.SC, HYB.SE = HYB.SE, HYB.COMB = HYB.COMB,
  HYC.SC = HYC.SC, HYC.SE = HYC.SE, HYT.SC = HYT.SC, HYT.SE = HYT.SE,
  HYC.SC.N = HYC.SC.N, HYC.SE.N = HYC.SE.N, HYT.SC.N = HYT.SC.N, HYT.SE.N = HYT.SE.N)

SPLIT.FIT.EXPOS <- list(
  MIX.SC = EXPO.MIX.SC, MIX.SE = EXPO.MIX.SE, HYB.SC = EXPO.HYB, HYB.SE = EXPO.HYB, HYB.COMB = EXPO.HYB,
  HYC.SC = EXPO.HYC, HYC.SE = EXPO.HYC, HYT.SC = EXPO.HYT, HYT.SE = EXPO.HYT,
  HYC.SC.N = EXPO.HYC.N, HYC.SE.N = EXPO.HYC.N, HYT.SC.N = EXPO.HYT.N, HYT.SE.N = EXPO.HYT.N)

## ---- Local parallel cluster: NB fits for the 13 split datasets ----
NUM.CORES.LOCAL <- max(1, detectCores() - 1)
cl <- makeCluster(NUM.CORES.LOCAL, type = "PSOCK")

t0 <- Sys.time()
SPLIT.FIT.RESULTS <- vector("list", length(SPLIT.FIT.MATS))
names(SPLIT.FIT.RESULTS) <- names(SPLIT.FIT.MATS)
for (nm in names(SPLIT.FIT.MATS)) {
cat(sprintf("[%s] fitting %-10s %d genes, %d cells\n", format(Sys.time(), "%H:%M:%S"), nm, nrow(SPLIT.FIT.MATS[[nm]]), ncol(SPLIT.FIT.MATS[[nm]])))
  flush.console()
  SPLIT.FIT.RESULTS[[nm]] <- fit_counts_offset(SPLIT.FIT.MATS[[nm]], SPLIT.FIT.EXPOS[[nm]], cl = cl)
}
cat(sprintf("split fits done: %d datasets in %.1f min\n", length(SPLIT.FIT.RESULTS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
stopCluster(cl)
## ---- end local parallel cluster ----

# Reorders the 13 fits into split-independent, then f_mean and f_disp
SPLIT.FITS <- SPLIT.FIT.RESULTS[c("MIX.SC", "MIX.SE", "HYC.SC", "HYC.SE", "HYT.SC", "HYT.SE", "HYC.SC.N", "HYC.SE.N", "HYT.SC.N", "HYT.SE.N", "HYB.SC", "HYB.SE", "HYB.COMB")]

NCELLS <- c(MIX.SC = ncol(MIX.SC), MIX.SE = ncol(MIX.SE), HYC.SC = ncol(HYC.SC), HYC.SE = ncol(HYC.SE), HYT.SC = ncol(HYT.SC), HYT.SE = ncol(HYT.SE), HYC.SC.N = ncol(HYC.SC.N), HYC.SE.N = ncol(HYC.SE.N), HYT.SC.N = ncol(HYT.SC.N), HYT.SE.N = ncol(HYT.SE.N), HYB.SC = ncol(HYB.SC), HYB.SE = ncol(HYB.SE), HYB.COMB = ncol(HYB.COMB))

## 2.3 Per-gene contrasts, bootstrap SEs, and permutation null

## Genes whose fit passes in every dataset of `fits`: a finite DISP below THETA.CAP (a detectable,
## non-degenerate dispersion), a raw mean count above the dataset's depth-scaled floor, and
## detection in at least min_expr_frac of the smallest dataset's cells. depth is a named vector of
## mean reads per cell for each dataset. The mean floor min_mean applies at the shallowest dataset
## and scales up in proportion to depth; the detected-cell floor is one absolute number of cells,
## so every dataset faces the same information requirement. The fit tables must be gene-aligned.
## Returns the passing gene names.
# Gene filter: requires an NB fit with a finite, non-degenerate
# dispersion estimate, a mean count above a floor that scales with dataset
# depth, and detection in an absolute number of cells set from the smallest
# dataset, in every one of the 13 datasets
SPLIT.DEPTH <- vapply(SPLIT.FIT.MATS, function(m) sum(m) / ncol(m), numeric(1))
GENES <- local({
  fits <- SPLIT.FITS
  ncells <- NCELLS
  depth <- SPLIT.DEPTH
  min_mean <- 0.001
  min_expr_frac <- 0.1
  groups <- names(fits)
  stopifnot(all(groups %in% names(ncells)), all(groups %in% names(depth)))
  genes <- rownames(fits[[1]])
  for (g in groups)
    if (!identical(rownames(fits[[g]]), genes))
      stop("fit frames are not gene-aligned; reorder to a common gene set first")
  floor_mean <- min_mean * depth[groups] / min(depth[groups])
  n_min <- ceiling(min_expr_frac * min(ncells[groups]))
  pass <- Reduce(`&`, lapply(groups, gene_pass_group, fits = fits, floor_mean = floor_mean, n_min = n_min))
  genes[which(pass)]
})

# Create bootstrap and permutation input files  
CONTRAST.MATS <- list(
  MIX.SC = MIX.SC[GENES, ], MIX.SE = MIX.SE[GENES, ],
  HYC.SC = HYC.SC[GENES, ], HYC.SE = HYC.SE[GENES, ],
  HYT.SC = HYT.SC[GENES, ], HYT.SE = HYT.SE[GENES, ],
  HYC.SC.N = HYC.SC.N[GENES, ], HYC.SE.N = HYC.SE.N[GENES, ],
  HYT.SC.N = HYT.SC.N[GENES, ], HYT.SE.N = HYT.SE.N[GENES, ],
  HYB.SC = HYB.SC[GENES, ], HYB.SE = HYB.SE[GENES, ], HYB.COMB = HYB.COMB[GENES, ])

CONTRAST.EXPOS <- list(MIX.SC = EXPO.MIX.SC, MIX.SE = EXPO.MIX.SE,
             HYC = EXPO.HYC, HYT = EXPO.HYT,
             HYC.N = EXPO.HYC.N, HYT.N = EXPO.HYT.N,
             HYB = EXPO.HYB)

CONTRAST.FITS <- lapply(SPLIT.FITS, contrast_fit_frame, genes = GENES)

# Expected rise in bfreq when HYB.COMB sums the two hybrid alleles, per gene
# in log2 units (0 to 1). The permutation job reads the dpar noise contrasts
# against it, and the BURST.CONTRASTS block below removes it from the dpar estimates.
PLOIDY.SHIFT <- ploidy_shift(CONTRAST.MATS, CONTRAST.EXPOS)[GENES]

## Pre-draws every permutation label vector (one list entry per permutation, all modes) from one
## seeded stream, so each gene is tested against the same relabelings. total/trans/dpar/inh shuffle
## pooled cell labels, cis swaps alleles within each hybrid cell, dom pairs random parent cells.
DRAWS <- make_draws(NCELLS, N.BOOT, SEED.BOOT)
PERMS <- local({
  ncells <- NCELLS
  NPERM <- N.PERM
  seed <- SEED.PERM
  set.seed(seed)
  nSC <- ncells[["MIX.SC"]]; nSE <- ncells[["MIX.SE"]]
  nHYC <- ncells[["HYC.SC"]]; nHYT <- ncells[["HYT.SC"]]
  nHYC.N <- ncells[["HYC.SC.N"]]; nHYT.N <- ncells[["HYT.SC.N"]]
  nHYB <- ncells[["HYB.COMB"]]
  lapply(seq_len(NPERM), perm_label_draw, nHYB = nHYB, nHYC = nHYC, nHYC.N = nHYC.N, nHYT = nHYT, nHYT.N = nHYT.N, nSC = nSC, nSE = nSE)
})

# gene_boot.R and gene_perm.R for the cluster
save(CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS, N.BOOT, SEED.BOOT, file = file.path(INPUT.DIR, "gene_boot1_inputs.rda"))
save(CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, PERMS, N.PERM, SEED.PERM, PLOIDY.SHIFT, file = file.path(INPUT.DIR, "gene_perm_inputs.rda"))

# DRAWS reassigned to a second seed for checking adequcey of bootstrap
DRAWS <- make_draws(NCELLS, N.BOOT, SEED.BOOT + 1)
save(CONTRAST.MATS, CONTRAST.EXPOS, CONTRAST.FITS, GENES, DRAWS, N.BOOT, file = file.path(INPUT.DIR, "gene_boot2_inputs.rda"))

## ---- Cluster round trip: Rscript gene_boot.R / gene_boot.R 2 / gene_perm.R ----
## gene_boot.R    reads gene_boot1_inputs.rda, writes gene_boot1_output.rda (BOOT.CONTRASTS)
## gene_boot.R 2  reads gene_boot2_inputs.rda, writes gene_boot2_output.rda (BOOT.CONTRASTS, loaded in 2.4)
## gene_perm.R    reads gene_perm_inputs.rda,  writes gene_perm_output.rda  (PERM.RESULTS)
load(file.path(OUTPUT.DIR, "gene_boot1_output.rda"))           # BOOT.CONTRASTS
load(file.path(OUTPUT.DIR, "gene_perm_output.rda"))            # PERM.RESULTS
## ---- end cluster round trip (gene_boot2_output.rda is loaded in 2.4) ----

## One row per mode: n, raw and attenuation-corrected mean-bfreq correlation, and its bootstrap CI.
## Benjamini-Hochberg FDR across genes for every permutation p-value
## column. Each contrast and quantity (for example mean_cis or
## bfreq_dpar_sc) is its own family of genes, so a q-value answers
## "what fraction of genes called in THIS contrast are expected to be
## false". p.adjust() leaves NAs in place, so genes with undefined fits
## are not counted as tests. Downstream classifiers read the _q columns.
## Puts the dpar noise estimates on the per-genome scale of the haploid
## parents. df: BURST.CONTRASTS after add_burst_contrasts(). fits:
## CONTRAST.FITS. shift: ploidy_shift() output. Each adjusted column keeps its
## raw twin as <col>_raw, and the column ploidy_shift records the shift applied.
## Mean columns pass through untouched. SEs stay as bootstrapped.
# BURST.CONTRASTS holds, per gene, the bootstrap point estimate and SE for the mean
# and dispersion (disp) contrasts in each of the eight modes, plus two derived
# quantities per mode: bfreq (burst frequency, the disp contrast
# itself) and bsize (burst size, mean minus disp in log2 space)
# The dpar noise estimates then leave the HYB.COMB allele-summing shift behind
# (ploidy_adjust_dpar), with the raw values kept in the _est_raw columns.
BURST.CONTRASTS <- local({
  df <- add_burst_contrasts(BOOT.CONTRASTS)
  fits <- CONTRAST.FITS
  shift <- PLOIDY.SHIFT
  s   <- unname(shift[df$gene])
  hyb <- fits$HYB.COMB[df$gene, ]
  cv_h <- 1 / hyb$MU + 2^s / hyb$DISP                # latent term scaled by 2^s
  df$ploidy_shift <- s
  for (sp in c("sc", "se")) {
    par  <- fits[[if (sp == "sc") "MIX.SC" else "MIX.SE"]][df$gene, ]
    cv_p <- 1 / par$MU + 1 / par$DISP
    for (q in c("bfreq", "bsize", "kbal", "cv2")) df[[paste0(dpar_est_col(q, sp = sp), "_raw")]] <- df[[dpar_est_col(q, sp = sp)]]
    df[[dpar_est_col("bfreq", sp = sp)]] <- df[[dpar_est_col("bfreq", sp = sp)]] - s
    df[[dpar_est_col("bsize", sp = sp)]] <- df[[dpar_est_col("bsize", sp = sp)]] + s
    df[[dpar_est_col("kbal", sp = sp)]]  <- df[[dpar_est_col("kbal", sp = sp)]]  - 2 * s
    df[[dpar_est_col("cv2", sp = sp)]]   <- ifelse(is.finite(cv_h) & cv_h > 0 & is.finite(cv_p) & cv_p > 0, log2(cv_h) - log2(cv_p), NA_real_)
  }
  df
})
PR <- PERM.RESULTS[match(BURST.CONTRASTS$gene, PERM.RESULTS$gene), ]   
# Benjamini-Hochberg FDR across genes, one family per contrast and quantity.
# Adds a _q column beside every _p column.
PR <- local({
  df <- PR
  method <- "BH"
  pc <- grep("_p(_ploidy|_ind)?$", names(df), value = TRUE)
  for (nm in pc) df[[sub("_p(?=(_ploidy|_ind)?$)", "_q", nm, perl = TRUE)]] <- p.adjust(df[[nm]], method = method)
  df
})
local({
  modes <- c("cis", "trans", "total")
  B <- 2000
  ## Gene-resampling bootstrap CI for statistics of eiv_components(). The CI is NA when more than half
  ## of the draws are non-finite.
  do.call(rbind, lapply(modes, eiv_mode_ci_row, contrasts = BURST.CONTRASTS, B = B))
})

# Diagnostic panel, before classification exists: top row cis vs trans
# (mean, burst frequency, burst size), middle row mean vs burst
# frequency (total, cis, trans), bottom row burst kinetics (total, cis,
# trans)
fig_pdf("extra/contrasts_panels.pdf", 13, 13, useDingbats = TRUE)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (panel in list(c("cis_trans", "mean"), c("cis_trans", "bfreq"), c("cis_trans", "bsize"),
                   c("mean_bfreq", "total"), c("mean_bfreq", "cis"), c("mean_bfreq", "trans"),
                   c("burst_kinetics", "total"), c("burst_kinetics", "cis"), c("burst_kinetics", "trans")))
  plot_contrast_scatter(BURST.CONTRASTS, panel[1], panel[2])
dev.off()

## 2.4 Bootstrap adequacy check
# Compares per-gene bootstrap SEs between the two seeds saved in 2.3
BURST.CONTRASTS2 <- local({ load(file.path(OUTPUT.DIR, "gene_boot2_output.rda")); add_burst_contrasts(BOOT.CONTRASTS) })

fig_pdf("extra/S_gene_seed_compare.pdf", 6, 9)
par(mfrow = c(3, 2))
SEED.CHECK <- setNames(lapply(c("total", "cis", "trans"), gene_seed_check, bc1 = BURST.CONTRASTS, bc2 = BURST.CONTRASTS2), c("total", "cis", "trans"))
dev.off()

# Expect correlation 0.9+ and a tight SE ratio 
do.call(rbind, lapply(names(SEED.CHECK), seed_check_rows, checks = SEED.CHECK, count_label = "n_genes"))

# Checkpoint: the fitted NB models, per-gene contrasts, and every Section 2
# object that later sections read.
save(GENES, HYB.COMB, CONTRAST.FITS, CONTRAST.MATS, CONTRAST.EXPOS, SPLIT.FITS,
     BOOT.CONTRASTS, BURST.CONTRASTS, PERM.RESULTS, PR, PLOIDY.SHIFT,
     HYB.SC, HYB.SE, EXPO.MIX.SC, EXPO.MIX.SE, EXPO.HYB, N.COEXPR, SEED.COEXPR,
     file = ckpt_path(2))

console_start(3)
##############################################################################
## 3. REGULATORY AND DOMINANCE RESULTS                                      ##
##############################################################################
# Built on BURST.CONTRASTS (effects) and PERM.RESULTS (p-values)

## 3.1 Regulatory classification: mean, burst frequency, burst size,
## kinetic balance
# REG.VEC$kbal classifies kbal (= bfreq - bsize), the rotated quantity that
# replaces a direct bfreq-vs-bsize comparison everywhere one would
# otherwise be tempted below: bfreq and bsize are not independent
# (bsize = mean - bfreq by construction), so classifying and comparing
# them side by side mostly re-tests the size of bfreq's own variance
# across genes, not biology. kbal and mean are the well-posed pair.
QUANTITIES <- c("mean", "bfreq", "bsize", "kbal")
REG.VEC <- lapply(setNames(QUANTITIES, QUANTITIES), function(q) reg_class_vec(BURST.CONTRASTS, PR, q))

## Figure 1: cis vs trans scatter above class counts, mean/burst
# frequency/burst size. Scatter row on top, class-count barplots below
fig_pdf("main/01_cis_trans.pdf", 13, 10)
par(mfrow = c(2, 3), mar = c(5, 4.5, 2, 1))
plot_cis_trans_class(BURST.CONTRASTS, PR, "mean")
plot_cis_trans_class(BURST.CONTRASTS, PR, "bfreq")
plot_cis_trans_class(BURST.CONTRASTS, PR, "bsize")
for (q in c("mean", "bfreq", "bsize")) class_count_barplot(REG.VEC[[q]], REG.CLASS, COLOR.LIST.1, QUANTITY.LABEL[[q]])
dev.off()

# BH FDR level for the class-overlap heatmaps below
OVERLAP.FDR <- 0.01

## Figure 2: regulatory class overlap heatmaps, all pairings.
# Cell text = fold enrichment (obs/exp); * = BH-corrected q < OVERLAP.FDR
# (0.01), adjusted across the cells of each heatmap.
# Color ramp = log2(obs/exp), diverging through white at zero.
class_overlap_triptych("main/02_overlap_regulatory.pdf", REG.VEC, REG.CLASS, width = 21, height = 7, mar = c(8, 8, 1, 1))

## 3.2 Significance histograms for total, cis, and trans effects
fig_pdf("extra/S_reg_sig_hist.pdf", 12, 10)
par(mfrow = c(3, 3))
SIG.HIST.YMAX <- list(total = c(mean = 300, bfreq = 600, bsize = 300), cis = c(mean = 400, bfreq = 600, bsize = 400), trans = c(mean = 400, bfreq = 400, bsize = 400))
for (md in c("total", "cis", "trans"))
  for (q in c("mean", "bfreq", "bsize")) sig_hist_panel(BURST.CONTRASTS, PR, md, q, SIG.HIST.YMAX[[md]][[q]])
dev.off()

## 3.3 Dominance classification: mean, burst frequency, burst size,
## kinetic balance
# dpar contrasts compare the combined hybrid (both alleles, offset by the
# full diploid library) with one haploid parent. The offset places both on a
# share of library scale, so an additive gene sits at the midparent on the
# mean axis and the library doubling cancels. On the noise axis the sum of two
# alleles averages their intrinsic noise, so bfreq rises in the hybrid and cv2
# falls. The BURST.CONTRASTS block removes that shift from the dpar estimates
DOM.VEC <- lapply(setNames(QUANTITIES, QUANTITIES), function(q) dom_class_vec(BURST.CONTRASTS, PR, q))
CLASS.VEC <- list(REG = REG.VEC, DOM = DOM.VEC)

# Sensitivity to the size of the shift. The _q_ind calls use s = 1 for every
# gene, the upper bound for fully independent equal alleles. Classes that hold
# under both bases do not depend on how private the allele noise is.
DOM.BFREQ.IND <- dom_class_vec(BURST.CONTRASTS, PR, "bfreq", basis = "ind")
PLOIDY.CHECK <- data.frame(
  median_shift = median(PLOIDY.SHIFT, na.rm = TRUE),
  raw_sc = median(BURST.CONTRASTS$bfreq_dpar_sc_est_raw, na.rm = TRUE),
  adj_sc = median(BURST.CONTRASTS$bfreq_dpar_sc_est,     na.rm = TRUE),
  raw_se = median(BURST.CONTRASTS$bfreq_dpar_se_est_raw, na.rm = TRUE),
  adj_se = median(BURST.CONTRASTS$bfreq_dpar_se_est,     na.rm = TRUE))
PLOIDY.CHECK   
table(own = DOM.VEC$bfreq, ind = DOM.BFREQ.IND)

## 3.4 Mean against burst frequency and burst size, parents and hybrids
# Burst size is the composite of the two (bsize = mean - bfreq in log2
# space)
# The mean axis of dpar shares the parental scale, so the three views compare
# directly.
fig_pdf("extra/S_mean_vs_bfreq_bsize_all_modes.pdf", 13, 18)
par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2, 1))
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "total",   reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bfreq", main = "mean vs burst frequency: parents")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "total",   reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bsize", main = "mean vs burst size: parents")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "dpar_sc", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bfreq", main = "mean vs burst frequency: hybrid vs Sc")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "dpar_sc", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bsize", main = "mean vs burst size: hybrid vs Sc")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "dpar_se", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bfreq", main = "mean vs burst frequency: hybrid vs Se")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "dpar_se", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bsize", main = "mean vs burst size: hybrid vs Se")
dev.off()

## Figure 3, supplement: mean against burst frequency and burst size,
# per-class slopes for the parental (total) contrast. Points neutral
# gray. Solid lines mark regulatory classes, dashed lines mark dominance. 
fig_pdf(file.path("extended", "S_03_mean_vs_bfreq_bsize.pdf"), 13, 6)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "total", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bfreq")
plot_mean_bfreq_class(BURST.CONTRASTS, PR, "total", reg_class = REG.VEC$mean, dom_class = DOM.VEC$mean, y_quantity = "bsize")
dev.off()

## Figure 5: dominance scatter above class counts, mean/burst frequency/burst size. 
fig_pdf("main/05_dominance.pdf", 13, 10)
par(mfrow = c(2, 3), mar = c(5, 4.5, 2, 1))
plot_dom_class(BURST.CONTRASTS, PR, "mean",  frame = "parent")
plot_dom_class(BURST.CONTRASTS, PR, "bfreq", frame = "parent")
plot_dom_class(BURST.CONTRASTS, PR, "bsize", frame = "parent")
for (q in c("mean", "bfreq", "bsize")) class_count_barplot(DOM.VEC[[q]], DOM.CLASS, COLOR.LIST.2, QUANTITY.LABEL[[q]])
dev.off()

## Figure 6: dominance class overlap heatmaps, all three pairings.
class_overlap_triptych("main/06_overlap_dominance.pdf", DOM.VEC, DOM.CLASS, width = 22.5, height = 7.5, mar = c(9, 9, 1, 1))

## Figure 4: rotated burst kinetics. Left, net mean change (x) against kinetic balance (y = bfreq - bsize)
## with SE bars, read from the mean_ and kbal_ columns add_burst_contrasts() stores. The rotation
## separates the two burst kinetics: movement along x is a change in total mean, movement along y is a
## shift between frequency and size. A gene is significant when its kinetic balance differs from 0 by a
## two-sided z-test (nominal p < sig); significant points are drawn black over grey n.s. points. KS keeps
## gene, y, sy, p and direction (sig_pos, sig_neg or ns) for the companion barplot at right.
fig_pdf("main/04_burst_kinetics.pdf", 9, 5)
par(mfrow = c(1, 2), mar = c(5, 4.5, 2, 1))
KS <- local({
  mode <- "total"
  sig <- 0.05
  main <- NULL
  bar_col <- adjustcolor(COLOR.GREY[["dark"]], 0.25)
  sig_col <- "black"
  ns_col <- COLOR.GREY[["mid"]]
  x <- BURST.CONTRASTS[[paste0("mean_", mode, "_est")]]; sx <- BURST.CONTRASTS[[paste0("mean_", mode, "_se")]]
  y <- BURST.CONTRASTS[[paste0("kbal_", mode, "_est")]]; sy <- BURST.CONTRASTS[[paste0("kbal_", mode, "_se")]]
  ok <- is.finite(x) & is.finite(y) & is.finite(sx) & is.finite(sy)
  genes <- BURST.CONTRASTS$gene[ok]; x <- x[ok]; sx <- sx[ok]; y <- y[ok]; sy <- sy[ok]
  p   <- 2 * pnorm(-abs(y) / sy)
  dir <- ifelse(!is.finite(p) | p >= sig, "ns", ifelse(y > 0, "sig_pos", "sig_neg"))
  ks  <- data.frame(gene = genes, y = y, sy = sy, p = p, direction = dir, stringsAsFactors = FALSE)
  sig_idx <- dir != "ns"
  if (is.null(main)) main <- paste0("burst kinetics: ", mode)
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = .sym(x, sx), ylim = .sym(y, sy), xlab = "net mean change (log2)", ylab = "frequency - amplitude (kinetic balance)", main = main)
  abline(h = 0, v = 0, col = COLOR.GREY[["dark"]])
  segments(x - sx, y, x + sx, y, col = bar_col)        # bars first
  segments(x, y - sy, x, y + sy, col = bar_col)
  points(x[!sig_idx], y[!sig_idx], pch = 16, cex = 0.5, col = ns_col)  # ns below
  points(x[ sig_idx], y[ sig_idx], pch = 16, cex = 0.5, col = sig_col)  # sig on top
  legend("topleft", legend = c(paste0("sig (p<", sig, ")"), "n.s."), col = c(sig_col, ns_col), pch = 16, bty = "n", cex = 0.8)
  invisible(ks)
})   # returns the gene/y/sy/p/direction table
bp_counts <- table(factor(KS$direction, levels = c("sig_pos","sig_neg","ns")))
barplot(bp_counts,
        col    = c(sig_pos = "black", sig_neg = COLOR.GREY[["dark"]], ns = COLOR.GREY[["light"]]),
        names.arg = c("sig > 0", "sig < 0", "n.s."),
        ylab   = "# of genes", las = 1, border = NA)
dev.off()

## Figure 7 (supplement): ggplot2 violins of kinetic balance by regulatory class (left panel) and
## dominance class (right panel). Values range roughly -1 to 1 (log2 units).
KBAL <- BURST.CONTRASTS$kbal_total_est   # kinetic balance, total; add_burst_contrasts() already derives this
fig_pdf(file.path("extended", "S_07_violins.pdf"), 8, 4.5)
print(local({
  value <- KBAL
  reg_class <- REG.VEC$mean
  dom_class <- DOM.VEC$mean
  ylab <- "kinetic balance (burst frequency - amplitude)"
  stopifnot(requireNamespace("ggplot2", quietly=TRUE))
  df <- rbind(
    data.frame(value=value, class=factor(reg_class,levels=REG.CLASS), panel="Regulatory"),
    data.frame(value=value, class=factor(dom_class, levels=DOM.CLASS), panel="Dominance"))
  df <- df[is.finite(df$value) & !is.na(df$class), ]
  df$panel <- factor(df$panel, levels=c("Regulatory","Dominance"))
  cols <- c(COLOR.LIST.1, COLOR.LIST.2[setdiff(names(COLOR.LIST.2),names(COLOR.LIST.1))])
  ggplot2::ggplot(df, ggplot2::aes(class, value, fill=class)) +
    ggplot2::geom_hline(yintercept=0, colour=COLOR.GREY[["mid"]], linewidth=0.3) +
    ggplot2::geom_violin(scale="width", trim=TRUE, colour=COLOR.GREY[["dark"]], linewidth=0.3) +
    ggplot2::geom_boxplot(width=0.12, outlier.size=0.3, fill="white", colour=COLOR.GREY[["dark"]], linewidth=0.3) +
    ggplot2::facet_wrap(~panel, scales="free_x") +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::labs(x=NULL, y=ylab) +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35, hjust=1), legend.position="none", strip.background=ggplot2::element_blank())
}))
dev.off()

## 3.5 Class relationship heatmaps
# log2 observed/expected, BH-corrected significance. 
HEATMAP.PANELS <- data.frame(
  y_kind = c(rep("REG", 3), rep("DOM", 3), "REG", "REG", "REG"), y_q = c("mean", "mean", "mean", "mean", "mean", "mean", "mean", "bfreq", "kbal"),
  x_kind = c(rep("REG", 3), rep("DOM", 3), "DOM", "DOM", "DOM"), x_q = c("bfreq", "bsize", "kbal", "bfreq", "bsize", "kbal", "mean", "bfreq", "kbal"))
class_heatmap_grid("extra/S_class_heatmaps.pdf", HEATMAP.PANELS, CLASS.VEC, OVERLAP.FDR, width = 13, height = 13, mfrow = c(3, 3))

## 3.6 Gene-identity overlap
# class_identity_overlap() computes Cohen's kappa (chance-corrected
# agreement between two class vectors across the whole table), a
# permutation null on kappa, and a per-class Jaccard overlap on gene
# membership.
REG.OVERLAP.MEAN.BFREQ  <- class_identity_overlap(clean_reg(REG.VEC$mean),  clean_reg(REG.VEC$bfreq), levels = REG.CLASS, nperm = 2000)
REG.OVERLAP.MEAN.BSIZE  <- class_identity_overlap(clean_reg(REG.VEC$mean),  clean_reg(REG.VEC$bsize), levels = REG.CLASS, nperm = 2000)
REG.OVERLAP.BFREQ.BSIZE <- class_identity_overlap(clean_reg(REG.VEC$bfreq), clean_reg(REG.VEC$bsize), levels = REG.CLASS, nperm = 2000)
REG.OVERLAP.MEAN.KBAL   <- class_identity_overlap(clean_reg(REG.VEC$mean),  clean_reg(REG.VEC$kbal),  levels = REG.CLASS, nperm = 2000)

# Structural bfreq-bsize correlation: for one mode, compares the observed bfreq-vs-bsize correlation with the
# one forced by bsize = mean - bfreq. Algebra gives Cov(bfreq, bsize) = Cov(bfreq, mean) - Var(bfreq).
# Setting Cov(bfreq, mean) = 0 yields the correlation expected with no biological coupling,
# rho_null = -Var(bfreq) / sqrt(Var(bfreq) * Var(bsize)). Vm, Vs and Cms come from eiv_components(),
# so the attenuation-corrected observed correlation is compared with that null on the same scale.
# STRUCT.BFREQ.BSIZE is a flat list with the fixed fields listed in STRUCT.FIELDS below.
STRUCT.FIELDS <- c("n", "Vm", "Vf", "Vs_bsize", "rho_bfreq_bsize_observed", "rho_bfreq_bsize_null", "excess_over_null")
STRUCT.BFREQ.BSIZE      <- local({
  mode <- "total"
  e <- eiv_components(BURST.CONTRASTS, mode)
  Vm <- unname(e["Vm"]); Vf <- unname(e["Vs"]); Cmf <- unname(e["Cms"])   # Vs/Cms here are bfreq's, not bsize's
  Vs_bsize <- Vm + Vf - 2 * Cmf
  rho_obs  <- suppressWarnings((Cmf - Vf) / sqrt(Vf * Vs_bsize))
  rho_null <- suppressWarnings(-Vf / sqrt(Vf * Vs_bsize))
  list(n = unname(e["n"]), Vm = Vm, Vf = Vf, Vs_bsize = Vs_bsize,
       rho_bfreq_bsize_observed = rho_obs, rho_bfreq_bsize_null = rho_null,
       excess_over_null = rho_obs - rho_null)
})
stopifnot("bfreq_bsize_structural() fields differ from STRUCT.FIELDS; update the function and the call sites together" =
          identical(names(STRUCT.BFREQ.BSIZE), STRUCT.FIELDS))

REG.OVERLAP.SUMMARY <- rbind(
  summarize_class_overlap(REG.OVERLAP.MEAN.BFREQ,  label = "Mean vs burst frequency"),
  summarize_class_overlap(REG.OVERLAP.MEAN.BSIZE,  label = "Mean vs burst size"),
  summarize_class_overlap(REG.OVERLAP.BFREQ.BSIZE, label = "Burst frequency vs burst size (structural, see STRUCT.BFREQ.BSIZE)"),
  summarize_class_overlap(REG.OVERLAP.MEAN.KBAL,   label = "Mean vs frequency-size balance"))
print(REG.OVERLAP.SUMMARY)
cat(sprintf("Bfreq-vs-bsize correlation: observed (attenuation-corrected) rho = %.3f, expected from the mean/bfreq variance imbalance alone (Cov(mean,bfreq) = 0) rho = %.3f, excess over that null = %.3f\n",
  STRUCT.BFREQ.BSIZE[["rho_bfreq_bsize_observed"]], STRUCT.BFREQ.BSIZE[["rho_bfreq_bsize_null"]], STRUCT.BFREQ.BSIZE[["excess_over_null"]]))
write.csv(REG.OVERLAP.SUMMARY, file.path(TABLE.DIR, "gene_identity_overlap_mean_bfreq_bsize.csv"), row.names = FALSE)

fig_pdf("extra/S_reg_overlap_kappa_null.pdf", 20, 5)
par(mfrow = c(1, 4))
for (ov in list(list(ov = REG.OVERLAP.MEAN.BFREQ,  lab = "mean vs burst frequency"),
                list(ov = REG.OVERLAP.MEAN.BSIZE,  lab = "mean vs burst size"),
                list(ov = REG.OVERLAP.BFREQ.BSIZE, lab = "burst frequency vs burst size (structural)"),
                list(ov = REG.OVERLAP.MEAN.KBAL,   lab = "mean vs frequency-size balance"))) {
  hist(ov$ov$kappa_null, breaks = 40, col = COLOR.GREY[["light"]], border = "white", main = paste0("Permutation null: ", ov$lab), xlab = "kappa (reshuffled labels)", xlim = range(c(ov$ov$kappa_null, ov$ov$kappa)))
  abline(v = ov$ov$kappa, col = COLOR.ACCENT, lwd = 2, lty = 2)
  legend("topright", bty = "n", cex = 0.8, legend = sprintf("observed kappa = %.3f (p = %.4f)", ov$ov$kappa, ov$ov$kappa_p))
}
dev.off()

# Checkpoint
save(REG.VEC, DOM.VEC,
     REG.OVERLAP.MEAN.BFREQ, REG.OVERLAP.MEAN.BSIZE, REG.OVERLAP.BFREQ.BSIZE, REG.OVERLAP.MEAN.KBAL,
     STRUCT.BFREQ.BSIZE, REG.OVERLAP.SUMMARY,
     file = ckpt_path(3))

console_start(4)
##############################################################################
## 4. CO-EXPRESSION                                                         ##
##############################################################################
## 4.1 Gene set, residuals, and point-estimate correlation matrices
# RESID holds, for each gene and cell, the standardized residual from the
# fitted NB model: (observed count - fitted mean) / fitted SD. This is
# the part of a cell's count NOT explained by that gene's own average
# expression level. Correlating
# these residuals across genes measures co-fluctuation: whether cells
# that are unexpectedly high for gene A also tend to be unexpectedly
# high for gene B.
#
# Reliability (RHO.ALL) is, for each gene, the fraction of its total
# count variance that is real cell-to-cell expression variation rather
# than Poisson counting noise: rho = mu / (mu + k), from the fitted NB
# mean and dispersion. A gene with low rho is dominated by counting
# noise, so its residual mostly reflects that noise rather than true
# co-fluctuation, and correlating two such genes would be measuring
# noise against noise. CO.GENES keeps every gene reliable enough for
# that correlation to be meaningful: the top 250 by expression, plus any
# other gene with rho > 0.3.
RHO.ALL <- gene_reliability(CONTRAST.FITS[c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")], GENES)
CO.GENES <- union(GENES[order(CONTRAST.FITS$MIX.SC$MU, decreasing = TRUE)][1:250], names(RHO.ALL)[RHO.ALL > 0.3])
length(CO.GENES)

# Both hybrid alleles summed, allele identity ignored, the same
# HYB.COMB and EXPO.HYB combination already used for the single-gene
# "hybrid total" (dpar_sc/dpar_se)
RESID <- list(
  MIX.SC = nb_residuals(MIX.SC[CO.GENES, ], EXPO.MIX.SC, CONTRAST.FITS$MIX.SC[CO.GENES, ]),
  MIX.SE = nb_residuals(MIX.SE[CO.GENES, ], EXPO.MIX.SE, CONTRAST.FITS$MIX.SE[CO.GENES, ]),
  HYB.SC = nb_residuals(HYB.SC[CO.GENES, ], EXPO.HYB,    CONTRAST.FITS$HYB.SC[CO.GENES, ]),
  HYB.SE = nb_residuals(HYB.SE[CO.GENES, ], EXPO.HYB,    CONTRAST.FITS$HYB.SE[CO.GENES, ]),
  HYB.COMB = nb_residuals(HYB.COMB[CO.GENES, ], EXPO.HYB, CONTRAST.FITS$HYB.COMB[CO.GENES, ]))

## Per-gene factor f that puts a HYB.COMB residual correlation on the
## per-genome scale of the haploid parents. Adjusted correlation of genes i
## and j is R_ij * f_i * f_j.
## Derivation. The residual of gene i is (y - mu) / sd, so the correlation of
## two genes is mu_i * mu_j * cov(L_i, L_j) / (sd_i * sd_j), where L is the
## latent deviation. Private allele noise does not correlate across genes, so
## the shared term e carries the whole covariance. Summing alleles leaves that
## covariance unchanged and only lowers each gene's latent variance by
## g = 2^-s (see ploidy_shift), which raises theta by 1/g, so the
## haploid-equivalent theta is g * theta. The covariance stays fixed, so the
## haploid-equivalent correlation equals the hybrid's times sd_h / sd_c for
## each gene, with sd_h^2 = mu + mu^2 / theta and sd_c^2 = mu + mu^2 / (g *
## theta). That ratio squared is
## f^2 = g * (mu + theta) / (mu + g * theta), which is at most 1 and equals g
## when counts are large. Adjusted hybrid correlations therefore shrink toward
## zero, in step with the lower latent variance of a single genome.
## fit: CONTRAST.FITS$HYB.COMB rows for the genes (MU, DISP); shift: PLOIDY.SHIFT for the same genes.
## A gene with no finite shift or theta gets f = 1 (no adjustment) and is counted below.
# Summing the two hybrid alleles lowers each gene's latent variance and leaves
# the between-gene covariance alone, so HYB.COMB residual correlations run
# higher than a haploid genome's. PLOIDY.F gives each gene the factor f that
# rescales the hybrid correlation to R_ij * f_i * f_j.
PLOIDY.F <- local({
  fit <- CONTRAST.FITS$HYB.COMB[CO.GENES, ]
  shift <- PLOIDY.SHIFT[CO.GENES]
  g  <- 2^(-shift); mu <- fit$MU; th <- fit$DISP
  f2 <- g * (mu + th) / (mu + g * th)
  f2[!is.finite(th) & is.finite(g)] <- 1     # Poisson limit, no latent noise to rescale
  f  <- sqrt(f2)
  f[!is.finite(f)] <- 1
  setNames(f, rownames(fit))
})
attr(RESID, "ploidy_f") <- PLOIDY.F
cat(sprintf("ploidy factor: median f = %.3f, %d of %d genes unadjusted\n", median(PLOIDY.F), sum(!is.finite(PLOIDY.SHIFT[CO.GENES])), length(CO.GENES)))

# COEXPR.POINT is the point estimate of the pairwise residual correlation
# matrix, decomposed into total (parents), cis (allele-specific within
# the hybrid), trans (the remainder), and dpar_sc/dpar_se (the hybrid's
# own correlation, alleles summed, against each parent separately)
COEXPR.POINT <- coexpr_decompose(RESID)

# Asks if the total correlation matrix is well described by a single shared
# factor or is pair-specific with no such structure?
fig_pdf("extra/S_coexpr_rank_check.pdf", 5, 5)
RANK.CHECK <- coexpr_rank_check(COEXPR.POINT$total, k = 1)
dev.off()
RANK.CHECK$r2 

## 4.2 Bootstrap for coexpression
# Two independent seeds at the same B
DRAWS.COEXPR <- make_coexpr_draws(ncol(RESID$MIX.SC), ncol(RESID$MIX.SE), ncol(RESID$HYB.SC), N.COEXPR, SEED.COEXPR)
save(RESID, COEXPR.POINT, DRAWS.COEXPR, file = file.path(INPUT.DIR, "coexpr_boot1_inputs.rda"))

DRAWS.COEXPR <- make_coexpr_draws(ncol(RESID$MIX.SC), ncol(RESID$MIX.SE), ncol(RESID$HYB.SC), N.COEXPR, SEED.COEXPR + 1)
save(RESID, COEXPR.POINT, DRAWS.COEXPR, file = file.path(INPUT.DIR, "coexpr_boot2_inputs.rda"))

## ---- Cluster round trip: Rscript coexpr_boot.R / coexpr_boot.R 2 ----
## coexpr_boot.R    reads coexpr_boot1_inputs.rda, writes coexpr_boot1_output.rda (CB)
## coexpr_boot.R 2  reads coexpr_boot2_inputs.rda, writes coexpr_boot2_output.rda (CB, loaded in 4.4)
load(file.path(OUTPUT.DIR, "coexpr_boot1_output.rda"))   # CB ($total, $cis, $trans, $lambda)
## ---- end cluster round trip (coexpr_boot2_output.rda is loaded in 4.4) ----

## 4.3 Bootstrap reliability check

## Checks that the observed bootstrap SE of co-expression pairs tracks the NB-predicted attenuation
## factor sqrt(rho_i * rho_j), to guide a rho floor for gene-set inclusion. Draws the SE-vs-attenuation
## scatter with binned medians on the open device.
## pairs: one of CB$total/cis/trans (data.frame with gene_i, gene_j, se). rho: named vector from
## gene_reliability() covering all genes in pairs. floors: candidate rho cutoffs, each reported with
## its retained pair count and median SE so the floor can be read off where SE stabilizes.
# Tests whether pairwise bootstrap SE tracks the NB-predicted
# attenuation factor sqrt(rho_i * rho_j).
RHO.TOTAL <- gene_reliability(CONTRAST.FITS[c("MIX.SC", "MIX.SE")], GENES)
fig_pdf("extra/S_coexpr_reliability_check.pdf", 7, 6)
REL.CHECK <- local({
  pairs <- CB$total
  rho <- RHO.TOTAL[CO.GENES]
  n_bins <- 10
  floors <- seq(0.1, 0.8, by = 0.1)
  attn <- sqrt(rho[pairs$gene_i] * rho[pairs$gene_j])
  ok   <- is.finite(attn) & is.finite(pairs$se)
  attn <- attn[ok]; se <- pairs$se[ok]

  bins  <- cut(attn, breaks = quantile(attn, seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  bin_x <- tapply(attn, bins, median)
  bin_y <- tapply(se,   bins, median)
  bin_n <- tapply(se,   bins, length)

  floor_summary <- do.call(rbind, lapply(floors, se_floor_row, attn = attn, se = list(se = se)))

 plot(attn, se, pch = 19, cex = 0.4, col = COLOR.GREY[["mid"]], xlab = expression(sqrt(rho[i] * rho[j])), ylab = "bootstrap SE", main = "Co-expression SE vs predicted reliability")
  lines(bin_x, bin_y, type = "b", pch = 19, lwd = 2)

  list(spearman = cor(attn, se, method = "spearman"),
       bins = data.frame(x = bin_x, y = bin_y, n = bin_n),
       floor_summary = floor_summary)
})
dev.off()

REL.CHECK$spearman        
REL.CHECK$floor_summary

## 4.4 Bootstrap adequacy check

## Verifies that a loaded coexpr_boot*_output.rda CB object has the pair count CO.GENES implies,
## so a stale or mismatched cluster output is caught before it reaches CB.CLASS, the pair lists or the
## seed comparison. The same check is applied to every CB object loaded.
# Confirms the bootstrap SE has converged at N.COEXPR
EXPECTED.PAIRS <- choose(length(CO.GENES), 2)
CB2 <- local({ load(file.path(OUTPUT.DIR, "coexpr_boot2_output.rda")); CB })
local({
  cb <- CB2
  expected_pairs <- EXPECTED.PAIRS
  label <- "coexpr_boot2_output.rda"
  if (nrow(cb$total) != expected_pairs)
  stop(sprintf("%s has %d pairs but CO.GENES expects %d; rerun the matching cluster job", label, nrow(cb$total), expected_pairs))
})

## Two-seed adequacy check for the co-expression bootstrap, one panel per contrast. If SE has not
## converged at this B, SEs underestimated by chance in one run can make many pairs look spuriously
## significant. CB and CB2 are two CB objects at the same B and gene set that differ only in seed; rows
## share the same gene_i/gene_j order (true when both come from the same CO.GENES and upper.tri() call).
fig_pdf("extra/S_coexpr_seed_compare.pdf", 15, 3.2)
par(mfrow = c(1, 5))
SEED.CHECK <- setNames(lapply(c("total", "cis", "trans", "dpar_sc", "dpar_se"), coexpr_seed_check, cb1 = CB, cb2 = CB2), c("total", "cis", "trans", "dpar_sc", "dpar_se"))
dev.off()

# expect correlation 0.9+ and a tight SE ratio
do.call(rbind, lapply(names(SEED.CHECK), seed_check_rows, checks = SEED.CHECK, count_label = "n_pairs"))

## 4.5 Pair-level classification

## Pair-level regulatory classification. total, cis and trans p-values are BH-adjusted across all tested
## pairs before classifying, which controls the false discovery rate over the quarter-million
## simultaneous pair tests. CB: the coexpr_boot.R output.
## One row per pair, with the five-class call from cis/trans carried alongside so downstream steps
## read the same BH-based class.
# CB.CLASS gives the five-way regulatory class per pair
# CB.DOM.CLASS is the co-expression analog of the single-gene dominance classification
CB.CLASS <- local({
  sig <- 0.05
  padj_total <- p.adjust(CB$total$p, "BH")
  padj_cis   <- p.adjust(CB$cis$p,   "BH")
  padj_trans <- p.adjust(CB$trans$p, "BH")
  cls <- classify_reg(padj_cis, padj_trans, CB$cis$est, CB$trans$est, sig = sig)$class
  data.frame(gene_i = CB$total$gene_i, gene_j = CB$total$gene_j,
             total_est = CB$total$est, total_padj = padj_total,
             cis_est   = CB$cis$est,   cis_padj   = padj_cis,
             trans_est = CB$trans$est, trans_padj = padj_trans,
             class = cls, stringsAsFactors = FALSE)
})
## Pair-level dominance classification (CB.DOM.CLASS), the co-expression analog of classify_dom() at the
## single-gene level. Compares the hybrid's own pairwise correlation (both alleles summed, allele identity
## ignored; dpar_sc/dpar_se from the co-expression decomposition, on the per-genome scale set by PLOIDY.F)
## against each parent's pairwise correlation separately: does the hybrid pair's co-expression resemble
## Sc, resemble Se, both (Additive), neither in a consistent direction (Over/Underdominant), or match
## neither parent while the two parents also differ from each other. Same two-test, BH-adjusted-first
## pattern as CB.CLASS, fed classify_dom() instead of classify_reg().
CB.DOM.CLASS <- local({
  sig <- 0.05
  padj_dpar_sc <- p.adjust(CB$dpar_sc$p, "BH")
  padj_dpar_se <- p.adjust(CB$dpar_se$p, "BH")
  cls <- classify_dom(padj_dpar_sc, padj_dpar_se, CB$dpar_sc$est, CB$dpar_se$est, sig = sig)$class
  data.frame(gene_i = CB$dpar_sc$gene_i, gene_j = CB$dpar_sc$gene_j,
             dpar_sc_est = CB$dpar_sc$est, dpar_sc_padj = padj_dpar_sc,
             dpar_se_est = CB$dpar_se$est, dpar_se_padj = padj_dpar_se,
             class = cls, stringsAsFactors = FALSE)
})

fig_pdf("extra/S_coexpr_dom_class.pdf", 6, 11)
par(mfrow = c(2, 1), mar = c(5, 4.5, 2, 1))
plot_coexpr_scatter(CB.DOM.CLASS, "dominance")
class_count_barplot(CB.DOM.CLASS$class, DOM.CLASS, COLOR.LIST.2, "Co-expression pair dominance class", ylab = "# of pairs")
dev.off()

table(factor(CB.DOM.CLASS$class, levels = DOM.CLASS))

## 4.6 Co-expression permutation null

## Draws for the joint permutation null of the total, cis, trans, dpar_sc and dpar_se spectra
## (rank-matched axis testing); one draw yields all five null spectra. idx: a random re-split of the
## pooled parent cells into pseudo-Sc / pseudo-Se groups of the original sizes (total's null).
## swap: an independent per-hybrid-cell coin flip of which allele is labelled Sc vs Se (cis's null).
## idx_dpar_sc / idx_dpar_se: random re-splits of the pooled Sc-parent + hybrid cells (resp. Se-parent
## + hybrid cells) into pseudo-parent and pseudo-hybrid groups of the original sizes. The total and
## cis draws are independent, so trans null = total null - cis null carries the same independence as
## the real estimates; the dpar nulls have no identity linking them to the other spectra.
# Builds null distributions of candidate-axis magnitude for the total,
# cis, and trans matrices. Rank-matched testing needs the null resolved
# at every rank up to N.KEEP, the number of top candidate axes tested in Section 4.9
N.SC  <- ncol(RESID$MIX.SC)
N.SE  <- ncol(RESID$MIX.SE)
N.HYB <- ncol(RESID$HYB.SC)
N.KEEP        <- 15
N.PERM.COEXPR <- 10000
DRAWS.PERM.COEXPR <- local({
  nSC <- N.SC
  nSE <- N.SE
  nH <- N.HYB
  B <- N.PERM.COEXPR
  seed <- SEED.COEXPR
  set.seed(seed)
  n_tot <- nSC + nSE
  lapply(seq_len(B), coexpr_perm_draw, nH = nH, nSC = nSC, nSE = nSE, n_tot = n_tot)
})
save(RESID, N.SC, N.SE, N.KEEP, DRAWS.PERM.COEXPR, file = file.path(INPUT.DIR, "coexpr_perm_inputs.rda"))

## ---- Cluster round trip: Rscript coexpr_perm.R ----
## Reads coexpr_perm_inputs.rda (saved above), writes coexpr_perm_output.rda
## (NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS, NULL.DPAR.SC.RANKS,
## NULL.DPAR.SE.RANKS, loaded in 4.9)
## ---- end cluster round trip (coexpr_perm_output.rda is loaded in 4.9) ----

## 4.7 Rank check and candidate axes
# Runs coexpr_rank_check() on each of the five divergence matrices and
# identifies candidate axes: the top 40
# axes by |eigenvalue|, filtered to a 1% variance floor and the
# MIN.EFFECTIVE.GENES participation-ratio floor

## Figure 11: cis vs trans scatter above class counts
fig_pdf("main/11_coexpr_cis_trans.pdf", 6, 11)
par(mfrow = c(2, 1), mar = c(5, 4.5, 2, 1))
plot_coexpr_scatter(CB.CLASS, "cis_trans")
class_count_barplot(CB.CLASS$class, REG.CLASS, COLOR.LIST.1, "co-expression pair regulatory class", ylab = "# of pairs")
dev.off()

MIN.EFFECTIVE.GENES <- 10
COEXPR.AXES   <- c("total", "cis", "trans", "dpar_sc", "dpar_se")

RANK.CHECK.LIST <- list()
CANDIDATE.LIST  <- list()

## Candidate axes from the rank check: the top n_candidate (40) axes by |eigenvalue|, each with its
## share of total squared-eigenvalue variance (axis_var) and its effective number of driving genes
## (axis_pr, the participation ratio 1/sum(loading^4) of a unit-length loading vector). sig_axes keeps
## axes with axis_var >= 0.01 and axis_pr >= MIN.EFFECTIVE.GENES; a dropped axis is named on the console.
## The table lists the top N.KEEP axes.
for (ax in COEXPR.AXES) {
  fig_pdf(sprintf("extra/S_coexpr_rank_check_%s_own.pdf", ax), 5, 5)
  RANK.CHECK.LIST[[ax]] <- coexpr_rank_check(COEXPR.POINT[[ax]], k = 1)
  dev.off()
  cat(sprintf("%s: rank-1 R^2 = %.3f\n", ax, RANK.CHECK.LIST[[ax]]$r2))

  CANDIDATE.LIST[[ax]] <- local({
    rank_check <- RANK.CHECK.LIST[[ax]]
    n_candidate <- 40
    var_floor <- 0.01
    eff_genes_min <- MIN.EFFECTIVE.GENES
    label <- ax
    n_candidate <- min(n_candidate, length(rank_check$values) - 1)
    axis_order  <- order(abs(rank_check$values), decreasing = TRUE)[1:n_candidate]
    axis_var    <- setNames(rank_check$values[axis_order]^2 / sum(rank_check$values^2), axis_order)
    axis_pr     <- setNames(sapply(axis_order, function(k) 1 / sum(rank_check$vectors[, k]^4)), axis_order)

    sig_axes <- as.integer(names(axis_var)[axis_var >= var_floor])
    dropped  <- sig_axes[axis_pr[as.character(sig_axes)] < eff_genes_min]
    if (length(dropped) > 0)
    cat(sprintf("%s: dropping axis %d (eff_genes = %.1f, below floor of %d)\n", label, dropped, axis_pr[as.character(dropped)], eff_genes_min), sep = "")
    sig_axes <- setdiff(sig_axes, dropped)

    list(axis_var = axis_var, axis_pr = axis_pr, sig_axes = sig_axes,
         table = data.frame(axis = axis_order, var = axis_var, eff_genes = axis_pr)[1:N.KEEP, ])
  })
  cat(sprintf("\n-- %s: candidate axes --\n", ax))
  print(CANDIDATE.LIST[[ax]]$table)
cat(sprintf("%s: axes clearing the variance and participation-ratio floors: %s\n", ax, paste(CANDIDATE.LIST[[ax]]$sig_axes, collapse = ", ")))
}


## 4.8 Mixture model and GO enrichment for candidate axes 
# Fits a two-component mixture to each axis's loading vector, splits
# genes into two poles by posterior probability, and enriches each pole
# against CO.GENES

EXTRA.AXES.LIST <- list()

## Fits a two-component Gaussian mixture (mixtools::normalmixEM, random starts) to each candidate
## axis's gene loadings, assigns genes to the low or high pole by posterior probability (> 0.5), and
## runs BP GO enrichment on each pole (>= 5 genes) against co_genes. Writes one page per axis to
## pdf_path and returns a named list (axis<k>) with variance and participation-ratio stats, mixture
## parameters, the two gene sets and their enrichment results.
for (ax in COEXPR.AXES) {
  EXTRA.AXES.LIST[[ax]] <- local({
    rank_check <- RANK.CHECK.LIST[[ax]]
    mat <- COEXPR.POINT[[ax]]
    sig_axes <- CANDIDATE.LIST[[ax]]$sig_axes
    axis_var <- CANDIDATE.LIST[[ax]]$axis_var
    axis_pr <- CANDIDATE.LIST[[ax]]$axis_pr
    co_genes <- CO.GENES
    pdf_path <- sprintf(file.path(FIGURE.DIR, "extra/S_coexpr_%s_axis_mixtures.pdf"), ax)
    axis_label <- ax
    out <- list()
    pdf(pdf_path, width = 6, height = 5, useDingbats = FALSE)
    for (k in sig_axes) {
      load_k <- setNames(rank_check$vectors[, k], rownames(mat))
      mix_k  <- tryCatch(normalmixEM(load_k, k = 2), error = function(e) NULL)
      if (is.null(mix_k)) next

      plot(mix_k, loglik = FALSE, density = TRUE, xlab2 = sprintf("loading on %s %d", axis_label, k))

      lo <- which.min(mix_k$mu); hi <- which.max(mix_k$mu)
      post_lo  <- setNames(mix_k$posterior[, lo], names(load_k))
      post_hi  <- setNames(mix_k$posterior[, hi], names(load_k))
      genes_lo <- names(post_lo)[post_lo > 0.5]
      genes_hi <- names(post_hi)[post_hi > 0.5]

      enrich_lo <- axis_pole_enrichment(genes_lo, co_genes, "BP", min_genes = 5)
      enrich_hi <- axis_pole_enrichment(genes_hi, co_genes, "BP", min_genes = 5)

      out[[paste0("axis", k)]] <- list(
        var_explained = axis_var[as.character(k)], eff_genes = axis_pr[as.character(k)],
        mu = mix_k$mu[c(lo, hi)], sigma = mix_k$sigma[c(lo, hi)], lambda = mix_k$lambda[c(lo, hi)],
        genes_lo = genes_lo, genes_hi = genes_hi, enrich_lo = enrich_lo, enrich_hi = enrich_hi)
    }
    dev.off()
    out
  })

  cat(sprintf("\n-- %s: axis mixture components --\n", ax))
  print(round(sapply(EXTRA.AXES.LIST[[ax]], axis_mixture_summary), 3))
}

## 4.9 Permutation validation
# Validates the candidate axes of each matrix: the k-th largest
# candidate axis (by variance share) is tested against the null's own
# k-th largest squared eigenvalue, rank-matched 
load(file.path(OUTPUT.DIR, "coexpr_perm_output.rda"))   # NULL.TOTAL.RANKS, NULL.CIS.RANKS, NULL.TRANS.RANKS, NULL.DPAR.SC.RANKS, NULL.DPAR.SE.RANKS

NULL.RANKS.LIST <- list(total = NULL.TOTAL.RANKS, cis = NULL.CIS.RANKS, trans = NULL.TRANS.RANKS, dpar_sc = NULL.DPAR.SC.RANKS, dpar_se = NULL.DPAR.SE.RANKS)

VALIDATED.LIST <- list()

## Tests each of the top N.KEEP candidate axes against its own permutation null: the k-th largest
## candidate axis (by variance share) against the null's k-th largest squared eigenvalue (rank-matched,
## no single shared threshold). p-values use the add-one rule, BH q-values run across the candidate
## axes, and the mixture-derived extra axes (EXTRA.AXES.LIST) are restricted to axes with q < alpha.
## VALIDATED.LIST keeps that list with the p-value table and the validated axis numbers.
for (ax in COEXPR.AXES) {
  VALIDATED.LIST[[ax]] <- local({
    rank_check <- RANK.CHECK.LIST[[ax]]
    axis_var <- CANDIDATE.LIST[[ax]]$axis_var
    null_ranks <- NULL.RANKS.LIST[[ax]]
    extra_axes <- EXTRA.AXES.LIST[[ax]]
    n_top <- N.KEEP
    alpha <- 0.05
    total_ss       <- sum(rank_check$values^2)
    candidate_axes <- as.integer(names(axis_var)[1:n_top])
    candidate_raw  <- axis_var[as.character(candidate_axes)] * total_ss
    ## Add-one correction (as in class_identity_overlap()'s kappa_p),
    ## since a leading axis routinely beats every one of the permutation
    ## draws: without it, mean(null >= obs) reports an exact 0 that
    ## overstates precision no finite permutation count can support,
    ## rather than the true floor of 1 / (n_perm + 1)
    candidate_p    <- sapply(seq_along(candidate_raw), function(k) (1 + sum(null_ranks[, k] >= candidate_raw[k])) / (1 + nrow(null_ranks)))

    ## Benjamini-Hochberg FDR across the candidate axes of this matrix.
    ## Rank-matched nulls make the axis tests positively related, the
    ## setting in which BH keeps the false discovery rate at its nominal level.
    candidate_q    <- p.adjust(candidate_p, method = "BH")
    validated_axes <- candidate_axes[candidate_q < alpha]

    list(table = data.frame(axis = candidate_axes, raw = candidate_raw, p_value = candidate_p, q_value = candidate_q),
         validated_axes = validated_axes,
         extra_axes = extra_axes[intersect(paste0("axis", validated_axes), names(extra_axes))])
  })

  cat(sprintf("\n-- %s: axis validation --\n", ax))
  VT <- VALIDATED.LIST[[ax]]$table
  VT$p_value <- signif(VT$p_value, 4)   # default rounding hides the permutation floor (1 / (N.PERM.COEXPR + 1)) as 0.000
  VT$q_value <- signif(VT$q_value, 4)   # BH FDR across the candidate axes; validated axes are those with q below alpha
  print(VT)
  cat(sprintf("%s: validated axes = %s\n", ax, paste(VALIDATED.LIST[[ax]]$validated_axes, collapse = ", ")))
}

# Mixture and enrichment results for the validated total axes. 
EXTRA.AXES <- VALIDATED.LIST$total$extra_axes

# Enrichment results for every validated total axis, at enrichGO()'s
# default qvalueCutoff = 0.2
for (nm in names(EXTRA.AXES)) {
  cat(sprintf("\n-- total %s: enrichment, low-loading pole (top terms, qvalue < 0.20) --\n", nm))
  print_enrich_brief(EXTRA.AXES[[nm]]$enrich_lo, q = 0.2)
  cat(sprintf("\n-- total %s: enrichment, high-loading pole (top terms, qvalue < 0.20) --\n", nm))
  print_enrich_brief(EXTRA.AXES[[nm]]$enrich_hi, q = 0.2)
}

## 4.10 Total Axis 1
# Add MF, CC, and KEGG enrichment for Total axis 1
AXIS1 <- EXTRA.AXES$axis1
NEG.LOAD.GENES <- AXIS1$genes_lo
POS.LOAD.GENES <- AXIS1$genes_hi
length(NEG.LOAD.GENES); length(POS.LOAD.GENES); length(CO.GENES)

AXIS1.ENRICH <- lapply(list(NEG = list(genes = NEG.LOAD.GENES, bp = AXIS1$enrich_lo), POS = list(genes = POS.LOAD.GENES, bp = AXIS1$enrich_hi)), function(pole)
  c(list(BP = pole$bp), lapply(c(MF = "MF", CC = "CC", KEGG = "KEGG"), axis_pole_enrichment, genes = pole$genes, universe = CO.GENES)))

# Growth-rate check: log2 fold change in mean expression (Sc/Se) for
# each loading group
POS.MU.LOG2FC <- mu_log2fc(POS.LOAD.GENES, fits = CONTRAST.FITS)   # ribosome / translation
NEG.MU.LOG2FC <- mu_log2fc(NEG.LOAD.GENES, fits = CONTRAST.FITS)   # glycolysis / fermentation
summary(POS.MU.LOG2FC); summary(NEG.MU.LOG2FC)

fig_pdf("extra/S_coexpr_growth_check.pdf", 5, 5)
boxplot(list(POS = POS.MU.LOG2FC, NEG = NEG.MU.LOG2FC), ylab = "log2(Sc MU / Se MU)", main = "mean expression shift by loading group")
abline(h = 0, lty = 2, col = COLOR.GREY[["dark"]])
dev.off()

# Per-species mean expression level for the NEG group
NEG.MU.SC <- setNames(log2(CONTRAST.FITS$MIX.SC[NEG.LOAD.GENES, "MU"]), NEG.LOAD.GENES)
NEG.MU.SE <- setNames(log2(CONTRAST.FITS$MIX.SE[NEG.LOAD.GENES, "MU"]), NEG.LOAD.GENES)
summary(NEG.MU.SC); summary(NEG.MU.SE)

fig_pdf("extra/S_coexpr_neg_group_by_species.pdf", 5, 5)
boxplot(list(Sc = NEG.MU.SC, Se = NEG.MU.SE), ylab = "log2(MU)", main = "NEG group (glycolysis/fermentation): level by species")
dev.off()

NEG.LOG2FC.SORTED <- sort(NEG.MU.LOG2FC, decreasing = TRUE)
head(NEG.LOG2FC.SORTED, 15)   # most Sc-elevated / Se-depressed genes in the set

fig_pdf("extra/S_coexpr_neg_log2fc_hist.pdf", 6, 5)
hist(NEG.MU.LOG2FC, breaks = 30, xlab = "log2(Sc MU / Se MU)", main = "NEG group, per-gene expression shift")
abline(v = 0, lty = 2, col = COLOR.GREY[["dark"]])
dev.off()

# Enrichment for the genes at least 2-fold Sc-elevated within the NEG group
NEG.TAIL.GENES <- names(NEG.MU.LOG2FC)[NEG.MU.LOG2FC > 1]
length(NEG.TAIL.GENES)

NEG.TAIL.ENRICH <- lapply(c(BP = "BP", CC = "CC"), axis_pole_enrichment, genes = NEG.TAIL.GENES, universe = CO.GENES)

# Cis/trans decomposition of axis 1's eigenvalue
AXIS1.CT <- coexpr_axis_cis_trans(RANK.CHECK$loading1, COEXPR.POINT)
AXIS1.CT

## 4.11 Total Axis 2
# RESP.LOAD.GENES is axis 2's smaller pole (cellular respiration,
# oxidative phosphorylation, ion transport). AXIS2.BULK.GENES is the
# larger, translation-annotated pole.
RESP.LOAD.GENES  <- EXTRA.AXES$axis2$genes_hi
AXIS2.BULK.GENES <- EXTRA.AXES$axis2$genes_lo
length(RESP.LOAD.GENES); length(AXIS2.BULK.GENES)

# Per-species mean expression level for the respiration/OXPHOS group
RESP.MU.SC <- setNames(log2(CONTRAST.FITS$MIX.SC[RESP.LOAD.GENES, "MU"]), RESP.LOAD.GENES)
RESP.MU.SE <- setNames(log2(CONTRAST.FITS$MIX.SE[RESP.LOAD.GENES, "MU"]), RESP.LOAD.GENES)
summary(RESP.MU.SC); summary(RESP.MU.SE)
wilcox.test(RESP.MU.SC, RESP.MU.SE, paired = TRUE)

fig_pdf("extra/S_coexpr_resp_group_by_species.pdf", 5, 5)
boxplot(list(Sc = RESP.MU.SC, Se = RESP.MU.SE), ylab = "log2(MU)", main = "respiration/OXPHOS group: level by species")
dev.off()

# Cis/trans decomposition of axis 2's eigenvalue
AXIS2.CT <- coexpr_axis_cis_trans(RANK.CHECK$vectors[, 2], COEXPR.POINT)
AXIS2.CT

# For pairs within the ribosome group, within the
# glycolysis/fermentation group, and pairs crossing between them, how do
# they classify under the existing five-class scheme (CB.CLASS)?
CB.CLASS$grp_i <- loading_group(CB.CLASS$gene_i, pos_genes = POS.LOAD.GENES, neg_genes = NEG.LOAD.GENES)
CB.CLASS$grp_j <- loading_group(CB.CLASS$gene_j, pos_genes = POS.LOAD.GENES, neg_genes = NEG.LOAD.GENES)
CB.CLASS$pair_type <- with(CB.CLASS, ifelse(is.na(grp_i) | is.na(grp_j), NA, ifelse(grp_i == grp_j, paste0("within_", grp_i), "cross")))

table(CB.CLASS$pair_type, CB.CLASS$class)

## 4.12 Burst frequency/size consistency among co-expressed genes
# Does a gene's burst-frequency class agree with its burst-size class as
# often within the co-expressed subset as it does genome-wide, or does
# the co-expressed population behave differently? CO.MEAN.KBAL
# is the comparison that isolates a real shift, since kbal does not
# inherit the same structural link to bfreq that bsize does.
CO.IDX         <- match(CO.GENES, BURST.CONTRASTS$gene)
CO.VEC         <- lapply(REG.VEC, function(v) v[CO.IDX])

BFREQ.BSIZE.OVERLAP.CO <- class_identity_overlap(clean_reg(CO.VEC$bfreq), clean_reg(CO.VEC$bsize), levels = REG.CLASS, nperm = 2000)
MEAN.KBAL.OVERLAP.CO   <- class_identity_overlap(clean_reg(CO.VEC$mean),  clean_reg(CO.VEC$kbal),  levels = REG.CLASS, nperm = 2000)

cat(sprintf("Bfreq vs bsize class concordance (structural, see STRUCT.BFREQ.BSIZE): genome-wide = %.3f (n = %d, kappa = %.3f, p = %.4f), CO.GENES = %.3f (n = %d, kappa = %.3f, p = %.4f)\n",
  REG.OVERLAP.BFREQ.BSIZE$concordance, REG.OVERLAP.BFREQ.BSIZE$n_genes, REG.OVERLAP.BFREQ.BSIZE$kappa, REG.OVERLAP.BFREQ.BSIZE$kappa_p,
  BFREQ.BSIZE.OVERLAP.CO$concordance,  BFREQ.BSIZE.OVERLAP.CO$n_genes,  BFREQ.BSIZE.OVERLAP.CO$kappa,  BFREQ.BSIZE.OVERLAP.CO$kappa_p))
cat(sprintf("Mean vs frequency-size balance class concordance: genome-wide = %.3f (n = %d, kappa = %.3f, p = %.4f), CO.GENES = %.3f (n = %d, kappa = %.3f, p = %.4f)\n",
  REG.OVERLAP.MEAN.KBAL$concordance, REG.OVERLAP.MEAN.KBAL$n_genes, REG.OVERLAP.MEAN.KBAL$kappa, REG.OVERLAP.MEAN.KBAL$kappa_p,
  MEAN.KBAL.OVERLAP.CO$concordance,  MEAN.KBAL.OVERLAP.CO$n_genes,  MEAN.KBAL.OVERLAP.CO$kappa,  MEAN.KBAL.OVERLAP.CO$kappa_p))

class_heatmap_grid("extra/S_coexpr_burst_mechanism_heatmap.pdf",
                   data.frame(y_kind = "REG", y_q = c("bfreq", "mean"), x_kind = "REG", x_q = c("bsize", "kbal")),
                   list(REG = CO.VEC), OVERLAP.FDR, width = 10, height = 5, mfrow = c(1, 2), kind_words = c(REG = ""), suffix = " (CO.GENES)")

AXIS.POLE.GROUPS <- list(Axis1.Ribosome    = POS.LOAD.GENES,
                        Axis1.Glycolysis  = NEG.LOAD.GENES,
                        Axis2.Respiration = RESP.LOAD.GENES,
                        Axis2.Bulk        = AXIS2.BULK.GENES)

# Fraction of each pole classified any-trans on each burst-kinetics axis
t(sapply(AXIS.POLE.GROUPS, pole_trans_fractions, contrasts = BURST.CONTRASTS, bfreq_class = REG.VEC$bfreq, bsize_class = REG.VEC$bsize, kbal_class = REG.VEC$kbal))

## 4.13 Cis and trans candidate axes
# Runs the detailed enrichment on
# every cis and trans axis that cleared validation
for (ax_mode in c("cis", "trans")) {
  validated <- VALIDATED.LIST[[ax_mode]]$validated_axes
  extra <- VALIDATED.LIST[[ax_mode]]$extra_axes
  if (length(validated) == 0) {
    cat(sprintf("\n-- %s: no validated axes --\n", ax_mode))
    next
  }
  for (nm in validated) {
    key <- paste0("axis", nm)
    cat(sprintf("\n-- %s %s: enrichment, low-loading pole (top terms, qvalue < 0.20) --\n", ax_mode, nm))
    print_enrich_brief(extra[[key]]$enrich_lo, q = 0.2)
    cat(sprintf("\n-- %s %s: enrichment, high-loading pole (top terms, qvalue < 0.20) --\n", ax_mode, nm))
    print_enrich_brief(extra[[key]]$enrich_hi, q = 0.2)
  }
}

# Checkpoint
save(RESID, PLOIDY.F, COEXPR.POINT, CB, CB2, CB.CLASS, CB.DOM.CLASS,
     RANK.CHECK.LIST, CANDIDATE.LIST, EXTRA.AXES.LIST, VALIDATED.LIST,
     AXIS1.CT, AXIS2.CT, BFREQ.BSIZE.OVERLAP.CO, AXIS1.ENRICH, NEG.TAIL.ENRICH,
     file = ckpt_path(4))

console_start(5)
##############################################################################
## 5. INTRINSIC / EXTRINSIC NOISE (allele-level co-expression)             ##
##############################################################################
## 5.1 Reliability calibration
# Diagnostic only, run to determine reliability floor for the
# allele-pair correlation analysis. 
RESID.ALLELE <- list(
  HYB.SC = nb_residuals(HYB.SC[GENES, ], EXPO.HYB, CONTRAST.FITS$HYB.SC[GENES, ]),
  HYB.SE = nb_residuals(HYB.SE[GENES, ], EXPO.HYB, CONTRAST.FITS$HYB.SE[GENES, ]))

# Per-allele reliability
RHO.HYB.SC <- gene_reliability(CONTRAST.FITS["HYB.SC"], GENES)
RHO.HYB.SE <- gene_reliability(CONTRAST.FITS["HYB.SE"], GENES)

## sc, se: HYB.SC / HYB.SE residual matrices, genes x cells, same cell columns, built over all genes
## (the allele-pair correlation is per gene and needs no shared gene set). rho_sc, rho_se:
## gene_reliability() on the HYB.SC and HYB.SE fits separately, since the two alleles' reliabilities
## enter the attenuation sqrt(rho_Sc * rho_Se) individually. Returns the sampled genes, raw and
## disattenuated correlations with their SEs, the attenuation, and a floor summary.
fig_pdf("extra/S_intrinsic_reliability_check.pdf", 9, 4.5)
INTR.REL.CHECK <- local({
  sc <- RESID.ALLELE$HYB.SC
  se <- RESID.ALLELE$HYB.SE
  rho_sc <- RHO.HYB.SC
  rho_se <- RHO.HYB.SE
  n_per_bin <- 40
  n_bins <- 10
  B <- 300
  seed <- 1
  floors <- seq(0.05, 0.6, by = 0.05)
  ## A reliability-stratified sample of genes for the calibration check: up to n_per_bin genes from each
  ## quantile bin of the attenuation, so the low-reliability end (sparse among genes overall) is
  ## represented well enough to see whether the SE flattens there. Seeded for reproducibility.
  genes <- intersect(rownames(sc), rownames(se))
  attn  <- sqrt(rho_sc[genes] * rho_se[genes])
  samp  <- local({
    set.seed(seed)
    ok    <- is.finite(attn)
    bins  <- cut(attn[ok], breaks = quantile(attn[ok], seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
    genes <- names(attn)[ok]
    unname(unlist(tapply(genes, bins, function(g) sample(g, min(n_per_bin, length(g))))))
  })

  rho_obs <- row_cor(sc[samp, , drop = FALSE], se[samp, , drop = FALSE])
  ## Bootstrap SE of each row's allele-residual correlation: B hybrid-cell resamples per gene, with the
  ## SD of the resampled correlations returned per row. Looped per gene since it runs on the few hundred
  ## genes of the calibration sample.
  se_obs  <- local({
    sc <- sc[samp, , drop = FALSE]
    se <- se[samp, , drop = FALSE]
    set.seed(seed)
    n <- ncol(sc)
    vapply(seq_len(nrow(sc)), allele_cor_boot_se_row, numeric(1), B = B, n = n, sc = sc, se = se)
  })
  a       <- attn[samp]
  rho_true <- rho_obs / a
  ## Exact error propagation for a fixed-scale division: SE(x/a) = SE(x)/a
  se_true  <- se_obs / a

  bins  <- cut(a, breaks = quantile(a, seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  bin_x <- tapply(a, bins, median)

  floor_summary <- do.call(rbind, lapply(floors, se_floor_row, attn = a, se = list(se_obs = se_obs, se_true = se_true), count_label = "n_genes"))

  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  plot(a, se_obs, pch = 19, cex = 0.4, col = COLOR.GREY[["mid"]], xlab = expression(sqrt(rho[Sc] * rho[Se])), ylab = "bootstrap SE, raw correlation", main = "Raw allele correlation")
  lines(bin_x, tapply(se_obs, bins, median), type = "b", pch = 19, lwd = 2)

 plot(a, se_true, pch = 19, cex = 0.4, col = COLOR.GREY[["mid"]], xlab = expression(sqrt(rho[Sc] * rho[Se])), ylab = "bootstrap SE, disattenuated", main = "Disattenuated allele correlation")
  lines(bin_x, tapply(se_true, bins, median), type = "b", pch = 19, lwd = 2)

  list(genes = samp, rho_obs = rho_obs, se_obs = se_obs,
       rho_true = rho_true, se_true = se_true, attn = a,
       floor_summary = floor_summary)
})
dev.off()

INTR.REL.CHECK$floor_summary 

## 5.2 Depth confound check
# Both alleles of a gene share the same cell and therefore the same
# sequencing depth, so a positive allele-pair correlation could reflect
# depth-correlated structure the NB offset left in the residuals rather
# than shared biological noise. Check whether each allele's
# residuals correlate with total cell depth on their own, and whether the
# allele-pair correlation survives once depth is partialled out.
DEPTH.CELL <- EXPO.HYB
INTR.SAMP  <- INTR.REL.CHECK$genes

DEPTH.COR.SC <- sapply(INTR.SAMP, function(g) cor(RESID.ALLELE$HYB.SC[g, ], DEPTH.CELL))
DEPTH.COR.SE <- sapply(INTR.SAMP, function(g) cor(RESID.ALLELE$HYB.SE[g, ], DEPTH.CELL))

summary(DEPTH.COR.SC); summary(DEPTH.COR.SE)   # expect both centered near 0, no strong depth bias left in the residuals
mean(DEPTH.COR.SC > 0); mean(DEPTH.COR.SE > 0)   # expect close to 0.5 if depth bias is not systematic in one direction

## Correlation between two residual vectors with a third variable (per-cell
## depth) partialled out. Used to test whether an allele-pair correlation
## reflects shared biological noise or leftover depth structure the NB
## offset failed to remove; if depth is the driver, this collapses toward
## zero relative to the raw correlation, and if not, it tracks the raw
## correlation closely.
RHO.PARTIAL <- sapply(INTR.SAMP, partial_cor_depth, resid = RESID.ALLELE, depth = DEPTH.CELL)

fig_pdf("extra/S_intrinsic_depth_partial.pdf", 5, 5)
plot(INTR.REL.CHECK$rho_obs, RHO.PARTIAL, pch = 19, cex = 0.5, xlab = "raw allele correlation", ylab = "allele correlation, depth partialled out")
abline(0, 1, col = COLOR.ACCENT, lty = 2); abline(h = 0, lty = 2)
dev.off()

# How far the depth-partialled correlation moves from the raw correlation.
summary(RHO.PARTIAL - INTR.REL.CHECK$rho_obs)

## 5.3 Per-gene fractions
# Under a symmetric shared-factor model, where
# each allele's standardized residual noise splits into a component
# shared across both alleles and a component private to that allele, the
# shared fraction of variance equals the correlation between the two
# alleles. The disattenuated allele-pair correlation is therefore read
# directly as the extrinsic fraction, and its complement as the
# intrinsic fraction.
FLOOR <- 0.10
ATTN.FULL <- sqrt(RHO.HYB.SC[GENES] * RHO.HYB.SE[GENES])

keep        <- ATTN.FULL >= FLOOR
GENES.NOISE <- GENES[keep]

RHO.OBS.FULL  <- row_cor(RESID.ALLELE$HYB.SC[GENES.NOISE, , drop = FALSE], RESID.ALLELE$HYB.SE[GENES.NOISE, , drop = FALSE])
ATTN.KEEP     <- ATTN.FULL[keep]
RHO.TRUE.FULL <- RHO.OBS.FULL / ATTN.KEEP

# Disattenuation can push an estimate slightly outside [0, 1] from
# sampling noise alone, particularly near the floor
# Need to correct and flag these genes
EXTRINSIC.FRAC <- pmin(pmax(RHO.TRUE.FULL, 0), 1)
INTRINSIC.FRAC <- 1 - EXTRINSIC.FRAC

NOISE.DECOMP <- data.frame(gene = GENES.NOISE, rho_obs = RHO.OBS.FULL, attn = ATTN.KEEP, rho_true = RHO.TRUE.FULL, extrinsic_frac = EXTRINSIC.FRAC, intrinsic_frac = INTRINSIC.FRAC, clipped = RHO.TRUE.FULL < 0 | RHO.TRUE.FULL > 1)

# How many genes needed clipping? 
table(NOISE.DECOMP$clipped)

fig_pdf("extra/S_extrinsic_fraction_hist.pdf", 5, 5)
hist(NOISE.DECOMP$extrinsic_frac, breaks = 30, main = "Extrinsic noise fraction across genes", xlab = "estimated extrinsic fraction", col = COLOR.GREY[["mid"]], border = "white")
abline(v = median(NOISE.DECOMP$extrinsic_frac), col = COLOR.ACCENT, lwd = 2, lty = 2)
dev.off()

summary(NOISE.DECOMP$extrinsic_frac)

# Second copy of the decomposition excluding genes whose disattenuated
# correlation needed clipping to [0, 1]
NOISE.DECOMP.CLEAN <- NOISE.DECOMP[!NOISE.DECOMP$clipped, ]

## 5.4 Regulatory and dominance class summaries
# Mean and variance of the extrinsic fraction within each class of the
# existing five-class regulatory scheme and six-class dominance scheme,
# for mean, burst frequency, and burst size
class_idx       <- match(NOISE.DECOMP$gene,       BURST.CONTRASTS$gene)
class_idx_clean <- match(NOISE.DECOMP.CLEAN$gene, BURST.CONTRASTS$gene)

EXTFRAC.BY.CLASS <- class_stats(class_mean_var, NOISE.DECOMP$extrinsic_frac, CLASS.VEC, class_idx)

EXTFRAC.BY.CLASS.CLEAN <- class_stats(class_mean_var, NOISE.DECOMP.CLEAN$extrinsic_frac, CLASS.VEC, class_idx_clean)

EXTFRAC.BY.CLASS
EXTFRAC.BY.CLASS.CLEAN

# Ambiguous classifications can arise where a significant permutation
# p-value pairs with an unresolved bootstrap sign, on any of the three
# quantities. Exclude these from the ANOVA/Tukey tests and GO analysis.
AMBIG.GENES <- BURST.CONTRASTS$gene[Reduce(`|`, lapply(c(REG.VEC[c("mean", "bfreq", "bsize")], DOM.VEC[c("mean", "bfreq", "bsize")]), function(v) v == "Ambiguous"))]

NOISE.DECOMP.NOAMBIG       <- NOISE.DECOMP[!(NOISE.DECOMP$gene %in% AMBIG.GENES), ]
NOISE.DECOMP.CLEAN.NOAMBIG <- NOISE.DECOMP.CLEAN[!(NOISE.DECOMP.CLEAN$gene %in% AMBIG.GENES), ]

class_idx_noambig       <- match(NOISE.DECOMP.NOAMBIG$gene,       BURST.CONTRASTS$gene)
class_idx_clean_noambig <- match(NOISE.DECOMP.CLEAN.NOAMBIG$gene, BURST.CONTRASTS$gene)

# One-way ANOVA per classification, testing whether the class means
# reflect a real overall difference rather than sampling noise.
# Tukey pairwise comparisons follow automatically wherever the omnibus
# test clears sig 
ANOVA.BY.CLASS <- class_stats(class_anova, NOISE.DECOMP.NOAMBIG$extrinsic_frac, CLASS.VEC, class_idx_noambig)

ANOVA.BY.CLASS.CLEAN <- class_stats(class_anova, NOISE.DECOMP.CLEAN.NOAMBIG$extrinsic_frac, CLASS.VEC, class_idx_clean_noambig)

# Omnibus F, p, and eta-squared per classification
t(sapply(ANOVA.BY.CLASS,       function(a) c(F = a$f, df1 = a$df1, df2 = a$df2, p = a$p, eta_sq = a$eta_sq)))
t(sapply(ANOVA.BY.CLASS.CLEAN, function(a) c(F = a$f, df1 = a$df1, df2 = a$df2, p = a$p, eta_sq = a$eta_sq)))

# Pairwise Tukey tables
lapply(ANOVA.BY.CLASS,       `[[`, "tukey")
lapply(ANOVA.BY.CLASS.CLEAN, `[[`, "tukey")

## Figure 8: intrinsic / extrinsic fraction histogram, Poisson shot-noise corrected. The fraction is
## intr / (intr + max(extr, 0)) per gene, the share of allele-pair noise that is private to each allele.
## A gene with no measurable intrinsic noise has fraction 0; genes without positive allele means, or
## with no intrinsic and no extrinsic noise, are NA. INTR.FRAC is a numeric vector named by gene.
INTR.FRAC <- local({
  mats <- CONTRAST.MATS
  expos <- CONTRAST.EXPOS
  ie   <- intrinsic_extrinsic_components(mats, expos)
  frac <- ie$intr / (ie$intr + pmax(ie$extr, 0))
  frac[!ie$good | !is.finite(frac)] <- NA_real_
  setNames(frac, ie$gene)
})[BURST.CONTRASTS$gene]
INTR.FRAC[ATTN.FULL < FLOOR] <- NA_real_
fig_pdf("main/08_intrinsic_fraction.pdf", 7, 5)
par(mar = c(5, 4.5, 2, 1))
## Histogram of the intrinsic fraction (genes with a fraction in [0, 1]), with regulatory-class medians
## (triangles above) and dominance-class medians (inverted triangles below). The legend sits at topleft,
## clear of the tall bars on the right.
local({
  frac <- INTR.FRAC
  reg_class <- REG.VEC$mean
  dom_class <- DOM.VEC$mean
  brk <- 30
  main <- NULL
  ok <- is.finite(frac) & frac>=0 & frac<=1
  f  <- frac[ok]; rc <- reg_class[ok]; dc <- dom_class[ok]
  h  <- hist(f, breaks=seq(0,1,length.out=brk+1), plot=FALSE)
  top <- max(h$counts)
  plot(h, col=COLOR.GREY[["light"]], border="white", xlim=c(0, 1), ylim=c(-top*0.06, top*1.20), xlab="intrinsic / (intrinsic + extrinsic)", ylab="# of genes", main=if(is.null(main)) "" else main)
  rm <- class_median(rc, REG.CLASS, f = f); dm <- class_median(dc, DOM.CLASS, f = f)
  yR <- top*1.10; yD <- -top*0.04
  for (i in seq_along(rm))
    if (is.finite(rm[i])) segments(rm[i],0,rm[i],yR, col=COLOR.LIST.1[i], lty=3, lwd=0.7)
  points(rm, rep(yR,length(rm)), pch=17, col=COLOR.LIST.1, cex=1.1, xpd=NA)
  points(dm, rep(yD,length(dm)), pch=25, col=COLOR.LIST.2, bg=COLOR.LIST.2, cex=1.1, xpd=NA)
  legend("topleft", legend=c(REG.CLASS, NA, DOM.CLASS), col=c(COLOR.LIST.1, NA, COLOR.LIST.2), pch=c(rep(17, length(REG.CLASS)), NA, rep(25, length(DOM.CLASS))), pt.bg=c(rep(NA, length(REG.CLASS)), NA, COLOR.LIST.2), bty="n", cex=0.72)
})
dev.off()

## 5.5 Burst kinetics relationship
# Tests whether genes with a larger species difference in burst
# frequency or burst size also carry a different extrinsic fraction.
BFREQ.EXTFRAC.COR <- cor.test(abs(BURST.CONTRASTS$bfreq_total_est[class_idx]), NOISE.DECOMP$extrinsic_frac, method = "spearman", exact = FALSE)
BSIZE.EXTFRAC.COR <- cor.test(abs(BURST.CONTRASTS$bsize_total_est[class_idx]), NOISE.DECOMP$extrinsic_frac, method = "spearman", exact = FALSE)
KBAL.EXTFRAC.COR  <- cor.test(abs(BURST.CONTRASTS$kbal_total_est[class_idx]),  NOISE.DECOMP$extrinsic_frac, method = "spearman", exact = FALSE)
BFREQ.EXTFRAC.COR.CLEAN <- cor.test(abs(BURST.CONTRASTS$bfreq_total_est[class_idx_clean]), NOISE.DECOMP.CLEAN$extrinsic_frac, method = "spearman", exact = FALSE)
BSIZE.EXTFRAC.COR.CLEAN <- cor.test(abs(BURST.CONTRASTS$bsize_total_est[class_idx_clean]), NOISE.DECOMP.CLEAN$extrinsic_frac, method = "spearman", exact = FALSE)
KBAL.EXTFRAC.COR.CLEAN  <- cor.test(abs(BURST.CONTRASTS$kbal_total_est[class_idx_clean]),  NOISE.DECOMP.CLEAN$extrinsic_frac, method = "spearman", exact = FALSE)

BFREQ.EXTFRAC.COR;       BSIZE.EXTFRAC.COR;       KBAL.EXTFRAC.COR
BFREQ.EXTFRAC.COR.CLEAN; BSIZE.EXTFRAC.COR.CLEAN; KBAL.EXTFRAC.COR.CLEAN

# Checkpoint
save(NOISE.DECOMP, NOISE.DECOMP.CLEAN, NOISE.DECOMP.NOAMBIG, NOISE.DECOMP.CLEAN.NOAMBIG, AMBIG.GENES,
     EXTFRAC.BY.CLASS, EXTFRAC.BY.CLASS.CLEAN, ANOVA.BY.CLASS, ANOVA.BY.CLASS.CLEAN,
     BFREQ.EXTFRAC.COR, BSIZE.EXTFRAC.COR, KBAL.EXTFRAC.COR,
     BFREQ.EXTFRAC.COR.CLEAN, BSIZE.EXTFRAC.COR.CLEAN, KBAL.EXTFRAC.COR.CLEAN,
     file = ckpt_path(5))

console_start(6)
##############################################################################
## 6. PROMOTER ARCHITECTURE                                                 ##
##############################################################################
# Three sequence-based features: a position frequency matrix (Bucher
# 1990) score for the TATA box, the longest poly(dA:dT) tract, and
# predicted nucleosome occupancy from NuPoP (Xi et al. 2010) 

## 6.1 Genomes, annotations, and per-species promoter scores
# Promoters are defined as the sequence immediately upstream of each gene's start codon, 
# up to PROM.MAX.BP or the distance to the nearest neighboring gene,
# whichever is shorter.
PROM.MAX.BP <- 300
PROM.MIN.BP <- 50

GENOME.SC <- read_genome_fasta(file.path(GENOME.DIR, "S.cerevisiae_YPS1000_V3.2.fasta"))
GENOME.SE <- read_genome_fasta(file.path(GENOME.DIR, "S.eubayanus_CBS12357^T_V3.2.fasta"))
GENES.SC  <- read_gff_genes(file.path(GENOME.DIR, "S.cerevisiae_YPS1000_V3.2.gff"))
GENES.SE  <- read_gff_genes(file.path(GENOME.DIR, "S.eubayanus_CBS12357^T_V3.2.gff"))

# Mitochondrial genes are excluded
MITO.PATTERN <- "(^|_)(chr)?(mt|mito|mitochondrion|mitochondrial)$"
GENES.SC <- GENES.SC[!grepl(MITO.PATTERN, GENES.SC$seqid, ignore.case = TRUE), ]
GENES.SE <- GENES.SE[!grepl(MITO.PATTERN, GENES.SE$seqid, ignore.case = TRUE), ]

PROM.SC <- extract_promoters(GENOME.SC, GENES.SC, max_bp = PROM.MAX.BP, min_bp = PROM.MIN.BP)
PROM.SE <- extract_promoters(GENOME.SE, GENES.SE, max_bp = PROM.MAX.BP, min_bp = PROM.MIN.BP)
cat(sprintf("Sc: %d of %d promoters extended past the neighbor boundary to reach %d bp\n", sum(attr(PROM.SC, "extended")), length(PROM.SC), PROM.MIN.BP))
cat(sprintf("Se: %d of %d promoters extended past the neighbor boundary to reach %d bp\n", sum(attr(PROM.SE, "extended")), length(PROM.SE), PROM.MIN.BP))

SCORE.SC <- score_promoters(PROM.SC)
SCORE.SE <- score_promoters(PROM.SE)

# NuPoP occupancy
# NuPoP runs on the cluster. The inputs carry each species' chromosomes
# and promoter coordinates.
NUPOP.INPUTS.SC <- nupop_cluster_inputs(GENOME.SC, PROM.SC)
NUPOP.INPUTS.SE <- nupop_cluster_inputs(GENOME.SE, PROM.SE)
save(NUPOP.INPUTS.SC, NUPOP.INPUTS.SE, file = file.path(INPUT.DIR, "nupop_inputs.rda"))

## ---- Cluster round trip: Rscript nupop_occupancy.R ----
## Reads nupop_inputs.rda (saved above), writes nupop_output.rda (NUPOP.OCC.SC, NUPOP.OCC.SE)
load(file.path(OUTPUT.DIR, "nupop_output.rda"))
## ---- end cluster round trip ----
OCC.SC <- score_promoters_nupop(PROM.SC, NUPOP.OCC.SC)
OCC.SE <- score_promoters_nupop(PROM.SE, NUPOP.OCC.SE)
attr(OCC.SC, "failed_regions"); attr(OCC.SE, "failed_regions")
SCORE.SC <- merge(SCORE.SC, OCC.SC, by = "gene")
SCORE.SE <- merge(SCORE.SE, OCC.SE, by = "gene")

## 6.2 Validation: does each feature associate with burst frequency?
# Uses each species' own allele within the hybrid (SPLIT.FITS$HYC.SC.N/
# HYC.SE.N, Section 2.2). Each feature is tested against burst frequency (DISP) with
# mean expression (MU) as the covariate. Burst size is not tested separately: with MU as
# the covariate, log BSIZE = log MU - log DISP makes the BSIZE test the DISP test with the
# sign of the slope reversed.
TATA.SC      <- setNames(SCORE.SC$tata_score, SCORE.SC$gene)
TATA.SE      <- setNames(SCORE.SE$tata_score, SCORE.SE$gene)
POLYAT.SC    <- setNames(SCORE.SC$polyat_len, SCORE.SC$gene)
POLYAT.SE    <- setNames(SCORE.SE$polyat_len, SCORE.SE$gene)
OCC.SC.SCORE <- setNames(SCORE.SC$occ_score,  SCORE.SC$gene)
OCC.SE.SCORE <- setNames(SCORE.SE$occ_score,  SCORE.SE$gene)

NOISE.VALIDATE <- rbind(
  cbind(species = "Sc", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SC.N, TATA.SC,      "TATA score",                   response = "DISP")),
  cbind(species = "Sc", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SC.N, POLYAT.SC,    "poly(dA:dT) length",           response = "DISP")),
  cbind(species = "Sc", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SC.N, OCC.SC.SCORE, "nucleosome occupancy (NuPoP)", response = "DISP")),
  cbind(species = "Se", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SE.N, TATA.SE,      "TATA score",                   response = "DISP")),
  cbind(species = "Se", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SE.N, POLYAT.SE,    "poly(dA:dT) length",           response = "DISP")),
  cbind(species = "Se", parameter = "bfreq", architecture_noise_check(SPLIT.FITS$HYC.SE.N, OCC.SE.SCORE, "nucleosome occupancy (NuPoP)", response = "DISP")))
print(NOISE.VALIDATE)

## 6.3 Between-species divergence in architecture

## Merges Sc and Se promoter scores on gene identity (inner join: genes
## present in both tables, whether or not their scores are NA) and adds the
## between-species shift (Se minus Sc) in each feature. Columns carry _sc/_se
## suffixes, including tata_motif and polyat_tract, the best-matching
## sequences. tata_pos_delta and polyat_pos_delta give the shift in where
## each best match sits (bp upstream of the ATG): a large score or length
## shift with a small position shift points to a substitution at about the
## same site, while a large position shift with a small score or length
## shift is more consistent with an indel moving a similar element, or a
## jump to a different site. Both score tables must already carry the
## occ_score column (merge score_promoters_nupop()'s output first).
##
## occ_delta is the raw occupancy shift. occ_access_delta flips its sign
## so that a positive value means more accessible for all three features: a
## higher TATA score or longer poly(dA:dT) tract and a lower occupancy all
## mean a more open promoter. The direction tests, concordance bins and
## candidate tables (analysis.R Sections 6.2 to 6.5) use occ_access_delta
## for this reason.
# Per-gene promoter divergence between the species, merged on gene identity
ARCH <- local({
  sc_scores <- SCORE.SC
  se_scores <- SCORE.SE
  m <- merge(sc_scores, se_scores, by = "gene", suffixes = c("_sc", "_se"))
  m$tata_delta        <- m$tata_score_se - m$tata_score_sc
  m$tata_pos_delta    <- m$tata_pos_se   - m$tata_pos_sc
  m$polyat_delta      <- m$polyat_len_se - m$polyat_len_sc
  m$polyat_pos_delta  <- m$polyat_pos_se - m$polyat_pos_sc
  m$occ_delta         <- m$occ_score_se  - m$occ_score_sc
  m$occ_access_delta  <- -m$occ_delta
  m
})
ARCH <- ARCH[match(BURST.CONTRASTS$gene, ARCH$gene), ]

## 6.4 Directional concordance: does the promoter shift point the right way? (burst frequency, burst size, frequency-size balance)
# Tests whether the allele that gained more accessibility is the noisier
# allele in the hybrid, across all three features at once
# (promoter_direction_test()), separately for burst frequency, burst
# size, and kinetic balance.
PROM.DIRECTION.TEST.BFREQ <- promoter_direction_test(BURST.CONTRASTS, PR, ARCH, REG.VEC$bfreq, quantity = "bfreq")
PROM.DIRECTION.TEST.BSIZE <- promoter_direction_test(BURST.CONTRASTS, PR, ARCH, REG.VEC$bsize, quantity = "bsize")
PROM.DIRECTION.TEST.KBAL  <- promoter_direction_test(BURST.CONTRASTS, PR, ARCH, REG.VEC$kbal,  quantity = "kbal")
print(PROM.DIRECTION.TEST.BFREQ)
print(PROM.DIRECTION.TEST.BSIZE)
print(PROM.DIRECTION.TEST.KBAL)

# Does concordance rise toward the largest promoter shifts? Bins the
# same gene set by |delta| magnitude and checks whether concordance
# rises from the smallest bin toward the largest, for burst frequency,
# burst size, and kinetic balance
TATA.CONCORD.BINS.BFREQ   <- concordance_by_magnitude(ARCH$tata_delta,        BURST.CONTRASTS$bfreq_cis_est, REG.VEC$bfreq)
POLYAT.CONCORD.BINS.BFREQ <- concordance_by_magnitude(ARCH$polyat_delta,      BURST.CONTRASTS$bfreq_cis_est, REG.VEC$bfreq)
OCC.CONCORD.BINS.BFREQ    <- concordance_by_magnitude(ARCH$occ_access_delta,  BURST.CONTRASTS$bfreq_cis_est, REG.VEC$bfreq)
TATA.CONCORD.BINS.BSIZE   <- concordance_by_magnitude(ARCH$tata_delta,        BURST.CONTRASTS$bsize_cis_est, REG.VEC$bsize)
POLYAT.CONCORD.BINS.BSIZE <- concordance_by_magnitude(ARCH$polyat_delta,      BURST.CONTRASTS$bsize_cis_est, REG.VEC$bsize)
OCC.CONCORD.BINS.BSIZE    <- concordance_by_magnitude(ARCH$occ_access_delta,  BURST.CONTRASTS$bsize_cis_est, REG.VEC$bsize)
TATA.CONCORD.BINS.KBAL    <- concordance_by_magnitude(ARCH$tata_delta,        BURST.CONTRASTS$kbal_cis_est,  REG.VEC$kbal)
POLYAT.CONCORD.BINS.KBAL  <- concordance_by_magnitude(ARCH$polyat_delta,      BURST.CONTRASTS$kbal_cis_est,  REG.VEC$kbal)
OCC.CONCORD.BINS.KBAL     <- concordance_by_magnitude(ARCH$occ_access_delta,  BURST.CONTRASTS$kbal_cis_est,  REG.VEC$kbal)
cat("\n-- TATA score: concordance by magnitude, burst frequency --\n")
print(TATA.CONCORD.BINS.BFREQ)
cat("\n-- poly(dA:dT) length: concordance by magnitude, burst frequency --\n")
print(POLYAT.CONCORD.BINS.BFREQ)
cat("\n-- nucleosome occupancy (NuPoP): concordance by magnitude, burst frequency --\n")
print(OCC.CONCORD.BINS.BFREQ)
cat("\n-- TATA score: concordance by magnitude, burst size --\n")
print(TATA.CONCORD.BINS.BSIZE)
cat("\n-- poly(dA:dT) length: concordance by magnitude, burst size --\n")
print(POLYAT.CONCORD.BINS.BSIZE)
cat("\n-- nucleosome occupancy (NuPoP): concordance by magnitude, burst size --\n")
print(OCC.CONCORD.BINS.BSIZE)
cat("\n-- TATA score: concordance by magnitude, frequency-size balance --\n")
print(TATA.CONCORD.BINS.KBAL)
cat("\n-- poly(dA:dT) length: concordance by magnitude, frequency-size balance --\n")
print(POLYAT.CONCORD.BINS.KBAL)
cat("\n-- nucleosome occupancy (NuPoP): concordance by magnitude, frequency-size balance --\n")
print(OCC.CONCORD.BINS.KBAL)

fig_pdf("extra/S_promoter_concordance_by_magnitude.pdf", 13, 13)
par(mfrow = c(3, 3))
plot_concordance_by_magnitude(TATA.CONCORD.BINS.BFREQ,   main = "TATA (burst frequency)")
plot_concordance_by_magnitude(POLYAT.CONCORD.BINS.BFREQ, main = "poly(dA:dT) (burst frequency)")
plot_concordance_by_magnitude(OCC.CONCORD.BINS.BFREQ,    main = "nucleosome occupancy (burst frequency)")
plot_concordance_by_magnitude(TATA.CONCORD.BINS.BSIZE,   main = "TATA (burst size)")
plot_concordance_by_magnitude(POLYAT.CONCORD.BINS.BSIZE, main = "poly(dA:dT) (burst size)")
plot_concordance_by_magnitude(OCC.CONCORD.BINS.BSIZE,    main = "nucleosome occupancy (burst size)")
plot_concordance_by_magnitude(TATA.CONCORD.BINS.KBAL,    main = "TATA (frequency-size balance)")
plot_concordance_by_magnitude(POLYAT.CONCORD.BINS.KBAL,  main = "poly(dA:dT) (frequency-size balance)")
plot_concordance_by_magnitude(OCC.CONCORD.BINS.KBAL,     main = "nucleosome occupancy (frequency-size balance)")
dev.off()

## 6.5 Candidate genes (burst frequency, burst size, frequency-size balance)
# Intersects a significant cis component of divergence (Cis, Cis +
# Trans, or Compensatory) with a large shift in a promoter feature,
# ranked by |<quantity>_cis_est|. 
PROM.NOISE.CANDIDATES.BFREQ <- promoter_noise_candidates(BURST.CONTRASTS, PR, ARCH, REG.VEC$bfreq, quantity = "bfreq")
PROM.NOISE.CANDIDATES.BSIZE <- promoter_noise_candidates(BURST.CONTRASTS, PR, ARCH, REG.VEC$bsize, quantity = "bsize")
PROM.NOISE.CANDIDATES.KBAL  <- promoter_noise_candidates(BURST.CONTRASTS, PR, ARCH, REG.VEC$kbal,  quantity = "kbal")
write.csv(PROM.NOISE.CANDIDATES.BFREQ$tata,   file.path(TABLE.DIR, "tata_bfreq_candidates.csv"),   row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.BFREQ$polyat, file.path(TABLE.DIR, "polyat_bfreq_candidates.csv"), row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.BFREQ$occ,    file.path(TABLE.DIR, "occ_bfreq_candidates.csv"),    row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.BSIZE$tata,   file.path(TABLE.DIR, "tata_bsize_candidates.csv"),   row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.BSIZE$polyat, file.path(TABLE.DIR, "polyat_bsize_candidates.csv"), row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.BSIZE$occ,    file.path(TABLE.DIR, "occ_bsize_candidates.csv"),    row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.KBAL$tata,    file.path(TABLE.DIR, "tata_kbal_candidates.csv"),    row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.KBAL$polyat,  file.path(TABLE.DIR, "polyat_kbal_candidates.csv"),  row.names = FALSE)
write.csv(PROM.NOISE.CANDIDATES.KBAL$occ,     file.path(TABLE.DIR, "occ_kbal_candidates.csv"),     row.names = FALSE)

# Checkpoint
save(SCORE.SC, SCORE.SE, ARCH, NOISE.VALIDATE,
     PROM.DIRECTION.TEST.BFREQ, PROM.DIRECTION.TEST.BSIZE, PROM.DIRECTION.TEST.KBAL,
     TATA.CONCORD.BINS.BFREQ, POLYAT.CONCORD.BINS.BFREQ, OCC.CONCORD.BINS.BFREQ,
     TATA.CONCORD.BINS.BSIZE, POLYAT.CONCORD.BINS.BSIZE, OCC.CONCORD.BINS.BSIZE,
     TATA.CONCORD.BINS.KBAL, POLYAT.CONCORD.BINS.KBAL, OCC.CONCORD.BINS.KBAL,
     PROM.NOISE.CANDIDATES.BFREQ, PROM.NOISE.CANDIDATES.BSIZE, PROM.NOISE.CANDIDATES.KBAL,
     file = ckpt_path(6))

console_start(7)
##############################################################################
## 7. BROAD CELLULAR DIFFERENCES IN EXPRESSION AND REGULATION               ##
##############################################################################
## 7.1 Build Seurat objects
# Combines the parental and hybrid count matrices into merged matrices
# for joint clustering (PARENT, HYBRID, MERGE), with a per-sample
# prefix on cell names so no two samples share a column name once
# combined.
# The seven datasets (single samples, then the combined sets) share one pipeline below; YSC, HVG, PCS,
# RES.SWEEP and BOOT hold their Seurat object, variable-feature elbow, retained PCs, resolution sweep and
# bootstrap-validated clustering, each as a list named by dataset.
DS.NAMES   <- c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE", "PARENT", "HYBRID", "MERGE")
DS.LABELS  <- c(MIX.SC = "Sc parent", MIX.SE = "Se parent", HYB.SC = "Hybrid, Sc allele", HYB.SE = "Hybrid, Se allele",
                PARENT = "Both parents combined", HYBRID = "Hybrid combined", MERGE = "All four merged")
YSC <- list()
MIX.SC.M <- MIX.SC; colnames(MIX.SC.M) <- paste0("MIXSC_", colnames(MIX.SC.M))
MIX.SE.M <- MIX.SE; colnames(MIX.SE.M) <- paste0("MIXSE_", colnames(MIX.SE.M))
HYB.SC.M <- HYB.SC; colnames(HYB.SC.M) <- paste0("HYBSC_", colnames(HYB.SC.M))
HYB.SE.M <- HYB.SE; colnames(HYB.SE.M) <- paste0("HYBSE_", colnames(HYB.SE.M))
HYBRID <- cbind(HYB.SC.M,HYB.SE.M)
PARENT <- cbind(MIX.SC.M,MIX.SE.M)
MERGE <- cbind(PARENT,HYBRID)

# Convert a count matrix to a Seurat-ready sparse matrix with dash-
# delimited feature names.
YSC$MIX.SC <- CreateSeuratObject(counts = to_seurat_counts(MIX.SC))
YSC$MIX.SE <- CreateSeuratObject(counts = to_seurat_counts(MIX.SE))
YSC$HYB.SC <- CreateSeuratObject(counts = to_seurat_counts(HYB.SC))
YSC$HYB.SE <- CreateSeuratObject(counts = to_seurat_counts(HYB.SE))

YSC$PARENT <- CreateSeuratObject(counts = to_seurat_counts(PARENT)); YSC.cells <- Cells(YSC$PARENT)
Idents(object = YSC$PARENT, cells = YSC.cells[1:ncol(MIX.SC)]) <- "Sc Parent"
Idents(object = YSC$PARENT, cells = YSC.cells[(1+ncol(MIX.SC)):(ncol(MIX.SC)+ncol(MIX.SE))]) <- "Se Parent"

YSC$HYBRID <- CreateSeuratObject(counts = to_seurat_counts(HYBRID)); YSC.cells <- Cells(YSC$HYBRID)
Idents(object = YSC$HYBRID, cells = YSC.cells[1:ncol(HYB.SC)]) <- "Sc Hybrid"
Idents(object = YSC$HYBRID, cells = YSC.cells[(1+ncol(HYB.SC)):(ncol(HYB.SC)+ncol(HYB.SE))]) <- "Se Hybrid"

YSC$MERGE <- CreateSeuratObject(counts = to_seurat_counts(MERGE)); YSC.cells <- Cells(YSC$MERGE)
Idents(object = YSC$MERGE, cells = YSC.cells[1:ncol(MIX.SC)]) <- "Sc Parent"
Idents(object = YSC$MERGE, cells = YSC.cells[(1+ncol(MIX.SC)):(ncol(MIX.SC)+ncol(MIX.SE))]) <- "Se Parent"
Idents(object = YSC$MERGE, cells = YSC.cells[(1+ncol(MIX.SC)+ncol(MIX.SE)):(ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC))]) <- "Sc Hybrid"
Idents(object = YSC$MERGE, cells = YSC.cells[(1+ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC)):(ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC)+ncol(HYB.SE))]) <- "Se Hybrid"

## 7.2 Normalize, select features, and reduce dimensionality
# Log-normalize each dataset, select variable features by a
# data-driven elbow on the ranked standardized-variance curve, scale,
# and run PCA. MIN.CLUSTER.CELLS is the minimum number of cells any
# cluster may contain in the resolution sweep and bootstrap steps that follow.
MIN.CLUSTER.CELLS <- 50

PREP <- lapply(setNames(DS.NAMES, DS.NAMES), function(d) prepare_dataset(YSC[[d]], d, rownames(YSC$MERGE), FIGURE.DIR))
YSC <- lapply(PREP, `[[`, "obj")
HVG <- lapply(PREP, `[[`, "hvg")

# Number of PCs to retain per dataset (elbow_pcs()): the elbow of the ranked percent-variance curve.
ELBOW.INPUTS <- lapply(YSC, elbow_pcs)
PCS <- sapply(ELBOW.INPUTS, `[[`, "pcs")

fig_pdf("extra/S_pca_elbow.pdf", 6, 5)
for (nm in names(ELBOW.INPUTS)) {
  e <- ELBOW.INPUTS[[nm]]
  plot_df <- data.frame(pct = e$pct, cumu = e$cumu, rank = seq_along(e$pct))
  print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > e$pcs)) +
    geom_text() +
    scale_color_manual(values = c("TRUE" = COLOR.GREY[["mid"]], "FALSE" = COLOR.ACCENT), name = "beyond chosen PCs") +
    geom_vline(xintercept = 90, color = COLOR.GREY[["mid"]]) +
    geom_hline(yintercept = min(e$pct[e$pct > 5]), color = COLOR.GREY[["mid"]]) +
    ggtitle(nm) +
    theme_bw())
}
dev.off()

## 7.3 Cluster each dataset
# Clusters cells in PCA space restricted to the retained dimensions, at
# the coarsest resolution within tolerance of the observed maximum
# silhouette (sweep_cluster_resolution()), checks bootstrap stability
# against the sweep's highest-silhouette plateau
# (cluster_stability.R on the cluster), and keeps whichever candidate
# has the higher bootstrap mean ARI as the final clustering.
RES.SWEEP <- lapply(setNames(DS.NAMES, DS.NAMES), cluster_dataset, ysc = YSC, pcs = PCS, labels = DS.LABELS,
                    min_cells = MIN.CLUSTER.CELLS, figure_dir = FIGURE.DIR)
YSC <- lapply(RES.SWEEP, `[[`, "obj")

## ---- Cluster execution of the stability bootstrap (Section 7.3) ----
## CSTAB.INPUTS packages everything cluster_stability.R needs. For each dataset it records the
## candidate resolutions (chosen, plus the coarsest point of the highest-silhouette plateau), the
## reference partition at each candidate (computed here on the fitted neighbor graph, so the local and
## cluster sides share one labelling), sparse counts, and the pre-drawn resample matrix. key
## fingerprints the inputs so the returning output can be matched to them.
# Stability bootstrap on the cluster. The sweeps above fix the
# candidate resolutions and reference partitions. cluster_stability.R
# refits every replicate in parallel, picks each dataset's final
# resolution by bootstrap mean ARI, and runs the Section 7.4 marker
# enrichment on the four single datasets at that final resolution.
N.CLUSTER.BOOT <- 500; SEED.CLUSTER.BOOT <- 1
CSTAB.INPUTS <- local({
  sweeps <- RES.SWEEP
  counts <- setNames(mget(DS.NAMES), DS.NAMES)
  nfeatures <- sapply(HVG, `[[`, "n_features")
  dims_n <- PCS
  B <- N.CLUSTER.BOOT
  seed <- SEED.CLUSTER.BOOT
  metric <- "manhattan"
  ds <- names(sweeps)
  tasks <- do.call(rbind, lapply(ds, stability_task_rows, sweeps = sweeps))
  tasks$task <- sprintf("%s@%.2f", tasks$dataset, tasks$res)
  ref <- setNames(lapply(seq_len(nrow(tasks)), reference_partition, sweeps = sweeps, tasks = tasks), tasks$task)
  data <- setNames(lapply(ds, stability_dataset_inputs, counts = counts, dims_n = dims_n, metric = metric, nfeatures = nfeatures, sweeps = sweeps), ds)
  ## ---- Cluster-stability bootstrap ----
  ## A stable partition survives resampling of the cells. Each replicate
  ## resamples cells with replacement from the raw count matrix (a fitted
  ## Seurat object cannot represent a cell drawn twice), builds a fresh
  ## Seurat object with uniquified barcodes, and reruns Normalize /
  ## FindVariableFeatures / Scale / PCA / Neighbors / Clusters at the SAME
  ## nfeatures, dims, resolution and metric as the original fit, so the
  ## comparison isolates sampling variation. The resampled clustering is
  ## compared with the original labels of the same resampled cells (in draw
  ## order) by adjusted Rand index. A high mean ARI means the partition is
  ## reproducible rather than a boundary Louvain draws through continuous
  ## variation; a low or widely spread ARI marks the practical resolution
  ## limit for the dataset's cell count.
  ##
  ## Every resample is drawn up front, one matrix per dataset, from one seeded stream (boot_resample_matrix()).
  ## Fixing the draws before any Seurat call keeps each replicate an independent resample, and the same
  ## matrix serves every candidate resolution of a dataset, so the resolution comparison is paired.
  idx <- setNames(lapply(ds, boot_resample_matrix, B = B, counts = counts, seed = seed), ds)
  key <- list(tasks = tasks, cells = lapply(counts, colnames), B = B, seed = seed)
  list(tasks = tasks, ref = ref, data = data, idx = idx, key = key)
})
MARKER.DS <- c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")
CSTAB.MARKER.OBJS <- lapply(YSC[MARKER.DS], diet_for_markers)
KEGG.DATA <- kegg_local("sce")

## ---- Cluster round trip: Rscript cluster_stability.R ----
## Reads cluster_stability_inputs.rda (saved below), writes
## cluster_stability_output.rda (CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY)
save(CSTAB.INPUTS, CSTAB.MARKER.OBJS, DS.LABELS, KEGG.DATA, file = file.path(INPUT.DIR, "cluster_stability_inputs.rda"))

## Loads the cluster output and checks its key against the inputs. For each dataset the loop below then
## prints the bootstrap comparison table (assembled by assemble_cluster_stability()) and, if the
## bootstrap-validated final resolution differs from the one the resolution sweep chose, a note naming
## the override and the bootstrap mean ARI at each, so the switch is traceable in the log rather than
## silent.
load_cluster_output(file.path(OUTPUT.DIR, "cluster_stability_output.rda"), "cluster_stability.R")   # CSTAB.ARI, CSTAB.MARKERS, CSTAB.KEY
check_cluster_key(CSTAB.KEY, CSTAB.INPUTS$key, "cluster_stability_output.rda", "cluster_stability.R")
## ---- end cluster round trip ----
BOOT <- list()
for (d in DS.NAMES) {
  boot <- assemble_cluster_stability(CSTAB.INPUTS, CSTAB.ARI, d, RES.SWEEP[[d]]$obj)
  local({
    label <- DS.LABELS[[d]]
    cat(sprintf("%s: bootstrap comparison across candidate resolutions\n", label)); print(boot$table)
    chosen_row <- boot$table[boot$table$role == "chosen", ]
    if (boot$final_res != chosen_row$res) {
      final_ari <- boot$table$boot_mean_ari[boot$table$res == boot$final_res]
      cat(sprintf("%s: switching to resolution %.2f (bootstrap mean ARI %.3f vs %.3f at the originally chosen %.2f)\n",
                  label, boot$final_res, final_ari, chosen_row$boot_mean_ari, chosen_row$res))
    }
  })
  BOOT[[d]] <- boot
  YSC[[d]] <- boot$final_obj
}

## Robustness check on the annoy.metric choice. Clusters the same retained PCs under each of two metrics
## at the same resolution and records both partitions plus the adjusted Rand index between them. A high ARI means the metric does not change which cells
## group together and either is defensible; a low ARI makes the metric a
## decision to state and justify in Methods.
METRIC.CHECK.MIX.SC <- local({
  obj <- YSC$MIX.SC
  dims <- 1:PCS[["MIX.SC"]]
  resolution <- BOOT$MIX.SC$final_res
  metrics <- c("manhattan", "euclidean")
  cl <- lapply(metrics, metric_clusters, dims = dims, obj = obj, resolution = resolution)
  names(cl) <- metrics
  list(clusters = cl, ari = adjustedRandIndex(as.integer(cl[[1]]), as.integer(cl[[2]])))
})
cat(sprintf("Sc parent: Manhattan vs Euclidean ARI = %.3f\n", METRIC.CHECK.MIX.SC$ari))

fig_pdf("extra/S_umap_clustering_checks.pdf", 6, 5)
UMAP.TITLES <- c(MIX.SC = "Sc parent (mono-culture)", MIX.SE = "Se parent (mono-culture)", HYB.SC = "Hybrid, Sc allele counts",
                 HYB.SE = "Hybrid, Se allele counts", PARENT = "Both parents combined (Sc + Se)",
                 HYBRID = "Hybrid, both allele views combined", MERGE = "All four samples merged, clusters")
for (d in DS.NAMES) YSC[[d]] <- umap_dataset(YSC[[d]], PCS[[d]], UMAP.TITLES[[d]])
dev.off()

## 7.4 Cluster composition: identity, marker enrichment, and consistency
# Identifies which cells belong to which cluster within each dataset,
# tests what distinguishes each cluster from the rest of its own
# dataset by marker gene GO/KEGG enrichment, and checks how consistent
# clustering is between the hybrid allele views and their corresponding
# parents.
YSC.IDENTS <- as.numeric(Idents(YSC$MERGE)) - 1
SC.PAR.ID  <- 1:(ncol(MIX.SC))
SE.PAR.ID  <- (1+ncol(MIX.SC)):(ncol(MIX.SC)+ncol(MIX.SE))
SC.HYB.ID  <- (1+ncol(MIX.SC)+ncol(MIX.SE)):(ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC))
SE.HYB.ID  <- (1+ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC)):(ncol(MIX.SC)+ncol(MIX.SE)+ncol(HYB.SC)+ncol(HYB.SE))

SC.PAR.CLUSTERS <- cluster_sizes(SC.PAR.ID, idents = YSC.IDENTS, min_cells = MIN.CLUSTER.CELLS)
SE.PAR.CLUSTERS <- cluster_sizes(SE.PAR.ID, idents = YSC.IDENTS, min_cells = MIN.CLUSTER.CELLS)
SC.HYB.CLUSTERS <- cluster_sizes(SC.HYB.ID, idents = YSC.IDENTS, min_cells = MIN.CLUSTER.CELLS)
SE.HYB.CLUSTERS <- cluster_sizes(SE.HYB.ID, idents = YSC.IDENTS, min_cells = MIN.CLUSTER.CELLS)
cat(sprintf("Clusters with at least %d cells, by dataset:\n", MIN.CLUSTER.CELLS))
cat("Sc parent:", SC.PAR.CLUSTERS, "\nSe parent:", SE.PAR.CLUSTERS, "\nHybrid Sc allele:", SC.HYB.CLUSTERS, "\nHybrid Se allele:", SE.HYB.CLUSTERS, "\n")

RES.MIX.SC <- CSTAB.MARKERS$MIX.SC
plot_cluster_marker_enrichment(RES.MIX.SC, "Sc parent", file.path(FIGURE.DIR, "extra/S_mix_sc_cluster_marker_enrichment.pdf"))

RES.MIX.SE <- CSTAB.MARKERS$MIX.SE
plot_cluster_marker_enrichment(RES.MIX.SE, "Se parent", file.path(FIGURE.DIR, "extra/S_mix_se_cluster_marker_enrichment.pdf"))

RES.HYB.SC <- CSTAB.MARKERS$HYB.SC
plot_cluster_marker_enrichment(RES.HYB.SC, "Hybrid, Sc allele", file.path(FIGURE.DIR, "extra/S_hyb_sc_cluster_marker_enrichment.pdf"))

RES.HYB.SE <- CSTAB.MARKERS$HYB.SE
plot_cluster_marker_enrichment(RES.HYB.SE, "Hybrid, Se allele", file.path(FIGURE.DIR, "extra/S_hyb_se_cluster_marker_enrichment.pdf"))

# Pseudobulk correlation of each hybrid-allele cluster against each
# corresponding parent cluster.
fig_pdf("extra/S_sc_hybrid_cluster_consistency_scatter.pdf", 4*length(SC.HYB.CLUSTERS), 4*length(SC.PAR.CLUSTERS))
par(mfrow=c(length(SC.PAR.CLUSTERS), length(SC.HYB.CLUSTERS)))
for (pc in SC.PAR.CLUSTERS) for (hc in SC.HYB.CLUSTERS) {
  x <- rowSums(MIX.SC[,YSC.IDENTS[SC.PAR.ID] == pc]); y <- rowSums(HYB.SC[,YSC.IDENTS[SC.HYB.ID] == hc])
  plot(x, y, pch=19,cex=0.8,xlab=paste("Sc Cluster", pc),ylab=paste("Sc Hyb Cluster", hc))
  cat(sprintf("Sc cluster %d vs Sc hybrid cluster %d: cor = %.3f\n", pc, hc, cor(x, y)))
}
dev.off()

fig_pdf("extra/S_se_hybrid_cluster_consistency_scatter.pdf", 4*length(SE.HYB.CLUSTERS), 4*length(SE.PAR.CLUSTERS))
par(mfrow=c(length(SE.PAR.CLUSTERS), length(SE.HYB.CLUSTERS)))
for (pc in SE.PAR.CLUSTERS) for (hc in SE.HYB.CLUSTERS) {
  x <- rowSums(MIX.SE[,YSC.IDENTS[SE.PAR.ID] == pc]); y <- rowSums(HYB.SE[,YSC.IDENTS[SE.HYB.ID] == hc])
  plot(x, y, pch=19,cex=0.8,xlab=paste("Se Cluster", pc),ylab=paste("Se Hyb Cluster", hc))
  cat(sprintf("Se cluster %d vs Se hybrid cluster %d: cor = %.3f\n", pc, hc, cor(x, y)))
}
dev.off()

# Pseudobulk correlation of each Sc-parent cluster against all of Se,
# identifying which Sc-parent cluster is the species comparison point.
fig_pdf("extra/S_sc_cluster_vs_se_scatter.pdf", 6, 4*length(SC.PAR.CLUSTERS))
par(mfrow=c(length(SC.PAR.CLUSTERS),1))
SC.PAR.VS.SE.COR <- setNames(numeric(length(SC.PAR.CLUSTERS)), SC.PAR.CLUSTERS)
for (pc in SC.PAR.CLUSTERS) {
  x <- rowSums(MIX.SC[,YSC.IDENTS[SC.PAR.ID] == pc]); y <- rowSums(MIX.SE)
  plot(x, y, pch=19,cex=0.8,xlab=paste("Sc Cluster", pc),ylab="Se")
  SC.PAR.VS.SE.COR[as.character(pc)] <- cor(x, y)
  cat(sprintf("Sc cluster %d vs Se (all): cor = %.3f\n", pc, SC.PAR.VS.SE.COR[as.character(pc)]))
}
dev.off()

if (length(SE.PAR.CLUSTERS) == 1 && length(SC.PAR.VS.SE.COR) > 0) {
  SC.SPECIES.CLUSTER <- as.integer(names(SC.PAR.VS.SE.COR)[which.max(SC.PAR.VS.SE.COR)])
  SE.SPECIES.CLUSTER <- SE.PAR.CLUSTERS[1]
} else {
  cat(sprintf("Found %d Se-parent cluster(s) (expected 1); species-cluster identification skipped.\n", length(SE.PAR.CLUSTERS)))
  SC.SPECIES.CLUSTER <- NA_integer_
  SE.SPECIES.CLUSTER <- NA_integer_
}

# Consistent pairing between the Sc-allele and Se-allele hybrid
# clusters: for each Sc-hybrid cluster, the Se-hybrid cluster the same
# cells most often carry, kept as a pair only if that agreement runs
# both ways.
HYB.CROSSTAB <- table(Sc = YSC.IDENTS[SC.HYB.ID], Se = YSC.IDENTS[SE.HYB.ID])
sc_best_se <- apply(HYB.CROSSTAB, 1, function(row) as.integer(names(which.max(row))))
se_best_sc <- apply(HYB.CROSSTAB, 2, function(col) as.integer(names(which.max(col))))
CONSISTENT.PAIRS <- do.call(rbind, lapply(names(sc_best_se), consistent_pair_row, sc_best_se = sc_best_se, se_best_sc = se_best_sc))
N.CONSISTENT.PAIRS <- nrow(CONSISTENT.PAIRS)
cat(sprintf("%d mutually-consistent hybrid cluster pair(s) found. Sc-by-Se hybrid cluster crosstab:\n", N.CONSISTENT.PAIRS))
print(HYB.CROSSTAB)
print(CONSISTENT.PAIRS)

HYBRID.CONSISTENT.CLUSTER.CELLS <- which(mapply(function(sc, se) any(CONSISTENT.PAIRS$sc == sc & CONSISTENT.PAIRS$se == se), YSC.IDENTS[SC.HYB.ID], YSC.IDENTS[SE.HYB.ID]))
cat(sprintf("%d / %d hybrid cells (%.1f%%) fall in a consistent cluster pair\n", length(HYBRID.CONSISTENT.CLUSTER.CELLS), length(SC.HYB.ID), 100*length(HYBRID.CONSISTENT.CLUSTER.CELLS)/length(SC.HYB.ID)))

## 7.5 Within/between-cluster noise partitioning
# Partitions each gene's variance into within- and between-cluster
# components using each dataset's own final clustering, relates the
# within/between ratio to burst kinetics, compares the ratio between
# species, and relates it to the regulatory and dominance classification.
CL.MIX.SC <- Idents(YSC$MIX.SC)[colnames(CONTRAST.MATS$MIX.SC)]
WB.MIX.SC <- within_between_decomp(CONTRAST.MATS$MIX.SC, CONTRAST.EXPOS$MIX.SC, CL.MIX.SC)
cat("Sc parent: cells per cluster\n"); print(WB.MIX.SC$cluster_n)

fig_pdf("extra/S_within_between_hist_MIX.SC.pdf", 6, 5)
WB.HIST.MIX.SC <- plot_within_between_hist(WB.MIX.SC, main = "Sc parent: within/between-cluster variance ratio")
dev.off()
cat(sprintf("Sc parent: within/between ratio, n = %d genes, %.1f%% with within > between (ratio > 1)\n",
            WB.HIST.MIX.SC$n, 100 * WB.HIST.MIX.SC$frac_above_1))

BF.MIX.SC         <- CONTRAST.FITS$MIX.SC
MU.LOG2.MIX.SC    <- setNames(log2(BF.MIX.SC$MU), rownames(BF.MIX.SC))
BFREQ.LOG2.MIX.SC <- setNames(log2(BF.MIX.SC$DISP), rownames(BF.MIX.SC))
BSIZE.LOG2.MIX.SC <- setNames(log2(BF.MIX.SC$MU / BF.MIX.SC$DISP), rownames(BF.MIX.SC))

fig_pdf("extra/S_within_between_vs_burst_MIX.SC.pdf", 15, 5)
par(mfrow = c(1, 3), mar = c(5, 4.5, 2, 1))
WB.VS.MU.MIX.SC    <- plot_within_between_vs_quantity(WB.MIX.SC, MU.LOG2.MIX.SC,    xlab = "log2(mean)")
WB.VS.BFREQ.MIX.SC <- plot_within_between_vs_quantity(WB.MIX.SC, BFREQ.LOG2.MIX.SC, xlab = "log2(burst frequency)")
WB.VS.BSIZE.MIX.SC <- plot_within_between_vs_quantity(WB.MIX.SC, BSIZE.LOG2.MIX.SC, xlab = "log2(burst size)")
dev.off()

cat(sprintf("Sc parent: within/between vs mean,       n = %d, Spearman rho = %.3f\n", WB.VS.MU.MIX.SC$n,    WB.VS.MU.MIX.SC$rho))
cat(sprintf("Sc parent: within/between vs burst freq, n = %d, Spearman rho = %.3f\n", WB.VS.BFREQ.MIX.SC$n, WB.VS.BFREQ.MIX.SC$rho))
cat(sprintf("Sc parent: within/between vs burst size, n = %d, Spearman rho = %.3f\n", WB.VS.BSIZE.MIX.SC$n, WB.VS.BSIZE.MIX.SC$rho))

CL.MIX.SE <- Idents(YSC$MIX.SE)[colnames(CONTRAST.MATS$MIX.SE)]
WB.MIX.SE <- within_between_decomp(CONTRAST.MATS$MIX.SE, CONTRAST.EXPOS$MIX.SE, CL.MIX.SE)
cat("Se parent: cells per cluster\n"); print(WB.MIX.SE$cluster_n)

fig_pdf("extra/S_within_between_hist_MIX.SE.pdf", 6, 5)
WB.HIST.MIX.SE <- plot_within_between_hist(WB.MIX.SE, main = "Se parent: within/between-cluster variance ratio")
dev.off()
cat(sprintf("Se parent: within/between ratio, n = %d genes, %.1f%% with within > between (ratio > 1)\n",
            WB.HIST.MIX.SE$n, 100 * WB.HIST.MIX.SE$frac_above_1))

BF.MIX.SE         <- CONTRAST.FITS$MIX.SE
MU.LOG2.MIX.SE    <- setNames(log2(BF.MIX.SE$MU), rownames(BF.MIX.SE))
BFREQ.LOG2.MIX.SE <- setNames(log2(BF.MIX.SE$DISP), rownames(BF.MIX.SE))
BSIZE.LOG2.MIX.SE <- setNames(log2(BF.MIX.SE$MU / BF.MIX.SE$DISP), rownames(BF.MIX.SE))

fig_pdf("extra/S_within_between_vs_burst_MIX.SE.pdf", 15, 5)
par(mfrow = c(1, 3), mar = c(5, 4.5, 2, 1))
WB.VS.MU.MIX.SE    <- plot_within_between_vs_quantity(WB.MIX.SE, MU.LOG2.MIX.SE,    xlab = "log2(mean)")
WB.VS.BFREQ.MIX.SE <- plot_within_between_vs_quantity(WB.MIX.SE, BFREQ.LOG2.MIX.SE, xlab = "log2(burst frequency)")
WB.VS.BSIZE.MIX.SE <- plot_within_between_vs_quantity(WB.MIX.SE, BSIZE.LOG2.MIX.SE, xlab = "log2(burst size)")
dev.off()

cat(sprintf("Se parent: within/between vs mean,       n = %d, Spearman rho = %.3f\n", WB.VS.MU.MIX.SE$n,    WB.VS.MU.MIX.SE$rho))
cat(sprintf("Se parent: within/between vs burst freq, n = %d, Spearman rho = %.3f\n", WB.VS.BFREQ.MIX.SE$n, WB.VS.BFREQ.MIX.SE$rho))
cat(sprintf("Se parent: within/between vs burst size, n = %d, Spearman rho = %.3f\n", WB.VS.BSIZE.MIX.SE$n, WB.VS.BSIZE.MIX.SE$rho))

## Cross-species scatter of within_between_decomp()'s ratio (log10), one
## point per ortholog gene pair, for two datasets whose gene sets are
## already ortholog-matched by name (e.g. MIX.SC vs MIX.SE, both indexed
## by the shared GENES ortholog-pair set built earlier in the pipeline,
## so no additional ortholog mapping is needed here). Spearman rho
## reported, same rank-based convention used throughout this section.
fig_pdf("extra/S_within_between_cross_species.pdf", 6, 6)
WB.CROSS.SPECIES <- local({
  wb_a <- WB.MIX.SC
  wb_b <- WB.MIX.SE
  lab_a <- "Sc"
  lab_b <- "Se"
  main <- "Within/between ratio: Sc vs Se parent"
  m  <- merge(wb_a$table[, c("gene", "ratio_within_between")], wb_b$table[, c("gene", "ratio_within_between")], by = "gene", suffixes = c("_a", "_b"))
  ok <- is.finite(m$ratio_within_between_a) & m$ratio_within_between_a > 0 & is.finite(m$ratio_within_between_b) & m$ratio_within_between_b > 0
  x   <- log10(m$ratio_within_between_a[ok]); y <- log10(m$ratio_within_between_b[ok])
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.null(main)) main <- sprintf("n = %d genes, Spearman rho = %.3f", sum(ok), rho)
  plot(x, y, pch = 16, cex = 0.4, col = adjustcolor(COLOR.GREY[["dark"]], 0.38),
       xlab = sprintf("log10(within / between), %s", lab_a), ylab = sprintf("log10(within / between), %s", lab_b), main = main)
  abline(0, 1, lty = 3, col = COLOR.GREY[["mid"]])
  abline(lm(y ~ x), col = COLOR.ACCENT, lty = 2)
  invisible(list(n = sum(ok), rho = rho))
})
dev.off()
cat(sprintf("Within/between ratio, Sc vs Se parent, n = %d genes, Spearman rho = %.3f\n",
            WB.CROSS.SPECIES$n, WB.CROSS.SPECIES$rho))

CLASS.LEVELS <- list(REG = REG.CLASS, DOM = DOM.CLASS)
CLASS.COLORS <- list(REG = COLOR.LIST.1, DOM = COLOR.LIST.2)
for (kind in c("REG", "DOM"))
  for (q in c("mean", "bfreq", "bsize"))
    report_within_between_by_class(WB.MIX.SC, WB.MIX.SE, CLASS.VEC[[kind]][[q]], BURST.CONTRASTS$gene,
      CLASS.LEVELS[[kind]], CLASS.COLORS[[kind]], "Sc parent", "Se parent",
      sprintf("%s class (%s)", c(REG = "regulatory", DOM = "dominance")[[kind]], QUANTITY.LABEL[[q]]),
      file.path(FIGURE.DIR, sprintf("extra/S_within_between_by_%s_%s_class.pdf", tolower(kind), q)))

## 7.6 Cell-cycle and metabolic module scoring
# Scores every cell for cell-cycle phase (three curated regulons: the
# G1/S "CLN2 cluster", the mitotic "CLB2 cluster", and the M/G1
# boundary wave) and for three metabolic gene sets (glycolysis,
# oxidative phosphorylation, ribosome biogenesis), then validates each dataset's clustering
# against these scores directly.
S.GENES   <- c("YMR199W", "YPL256C", "YPR120C", "YGR109C", "YBR088C", "YKL113C", "YAR007C")  # CLN1, CLN2, CLB5, CLB6, POL30, RAD27, RFA1
G2M.GENES <- c("YGR108W", "YPR119W", "YDL155W", "YMR001C", "YGL116W")                        # CLB1, CLB2, CLB3, CDC5, CDC20
MG1.GENES <- c("YLR079W", "YDR146C", "YLR131C")                                              # SIC1, SWI5, ACE2

YSC$MIX.SC <- score_cell_cycle_by_cluster(YSC$MIX.SC, "Sc parent", file.path(FIGURE.DIR, "extra/S_cell_cycle_scoring_MIX.SC.pdf"))
YSC$MIX.SE <- score_cell_cycle_by_cluster(YSC$MIX.SE, "Se parent", file.path(FIGURE.DIR, "extra/S_cell_cycle_scoring_MIX.SE.pdf"))
YSC$HYB.SC <- score_cell_cycle_by_cluster(YSC$HYB.SC, "Hybrid, Sc allele", file.path(FIGURE.DIR, "extra/S_cell_cycle_scoring_HYB.SC.pdf"))
YSC$HYB.SE <- score_cell_cycle_by_cluster(YSC$HYB.SE, "Hybrid, Se allele", file.path(FIGURE.DIR, "extra/S_cell_cycle_scoring_HYB.SE.pdf"))

METABOLIC.GENE.SETS <- list(
  Glycolysis = go_gene_set("GO:0006096"),  # glycolytic process
  OXPHOS     = go_gene_set("GO:0006119"),  # oxidative phosphorylation
  RiBi       = go_gene_set("GO:0042254"))  # ribosome biogenesis

YSC$MIX.SC <- score_modules_by_cluster(YSC$MIX.SC, METABOLIC.GENE.SETS, "Sc parent", file.path(FIGURE.DIR, "extra/S_metabolic_scoring_MIX.SC.pdf"))
YSC$MIX.SE <- score_modules_by_cluster(YSC$MIX.SE, METABOLIC.GENE.SETS, "Se parent", file.path(FIGURE.DIR, "extra/S_metabolic_scoring_MIX.SE.pdf"))
YSC$HYB.SC <- score_modules_by_cluster(YSC$HYB.SC, METABOLIC.GENE.SETS, "Hybrid, Sc allele", file.path(FIGURE.DIR, "extra/S_metabolic_scoring_HYB.SC.pdf"))
YSC$HYB.SE <- score_modules_by_cluster(YSC$HYB.SE, METABOLIC.GENE.SETS, "Hybrid, Se allele", file.path(FIGURE.DIR, "extra/S_metabolic_scoring_HYB.SE.pdf"))

## 7.7 Covariate-noise diagnostic
# Tests whether per-cell NB Pearson residual noise, from the already-
# fit exposure-offset NB model, depends on a cell's cell-cycle position
# (single continuous axis, first PC of S.Score/G2M.Score) or metabolic
# state (discrete, joint k-means on the three metabolic module scores).
# Run on the four datasets that feed CONTRAST.MATS directly.
CC.AXIS.MIX.SC <- cell_cycle_continuum(YSC$MIX.SC)
CC.AXIS.MIX.SE <- cell_cycle_continuum(YSC$MIX.SE)
CC.AXIS.HYB.SC <- cell_cycle_continuum(YSC$HYB.SC)
CC.AXIS.HYB.SE <- cell_cycle_continuum(YSC$HYB.SE)

MET.STATE.MIX.SC <- metabolic_state_cluster(YSC$MIX.SC)
MET.STATE.MIX.SE <- metabolic_state_cluster(YSC$MIX.SE)
MET.STATE.HYB.SC <- metabolic_state_cluster(YSC$HYB.SC)
MET.STATE.HYB.SE <- metabolic_state_cluster(YSC$HYB.SE)

COV.DIAG.MIX.SC <- covariate_noise_diagnostic(CONTRAST.MATS$MIX.SC, CONTRAST.EXPOS$MIX.SC, CONTRAST.FITS$MIX.SC, CC.AXIS.MIX.SC, MET.STATE.MIX.SC, "Sc parent")
COV.DIAG.MIX.SE <- covariate_noise_diagnostic(CONTRAST.MATS$MIX.SE, CONTRAST.EXPOS$MIX.SE, CONTRAST.FITS$MIX.SE, CC.AXIS.MIX.SE, MET.STATE.MIX.SE, "Se parent")
COV.DIAG.HYB.SC <- covariate_noise_diagnostic(CONTRAST.MATS$HYB.SC, CONTRAST.EXPOS$HYB,    CONTRAST.FITS$HYB.SC, CC.AXIS.HYB.SC, MET.STATE.HYB.SC, "Hybrid, Sc allele")
COV.DIAG.HYB.SE <- covariate_noise_diagnostic(CONTRAST.MATS$HYB.SE, CONTRAST.EXPOS$HYB,    CONTRAST.FITS$HYB.SE, CC.AXIS.HYB.SE, MET.STATE.HYB.SE, "Hybrid, Se allele")

fig_pdf("extra/S_covariate_noise_diagnostic.pdf", 9, 4.5)
plot_covariate_noise_diagnostic(COV.DIAG.MIX.SC, "Sc parent")
plot_covariate_noise_diagnostic(COV.DIAG.MIX.SE, "Se parent")
plot_covariate_noise_diagnostic(COV.DIAG.HYB.SC, "Hybrid, Sc allele")
plot_covariate_noise_diagnostic(COV.DIAG.HYB.SE, "Hybrid, Se allele")
dev.off()

## 7.8 Species composition comparison and confound bound
# Compares Sc and Se on the cell-cycle axis and the three metabolic
# module scores (Cohen's d, rank-sum test), and multiplies each
# composition difference by the matching 7.7 noise-association effect
# size to bound the fraction of residual variance a composition
# difference of that size could explain.
MET.SCORE.MIX.SC <- as.matrix(YSC$MIX.SC[[c("Glycolysis1", "OXPHOS1", "RiBi1")]])
MET.SCORE.MIX.SE <- as.matrix(YSC$MIX.SE[[c("Glycolysis1", "OXPHOS1", "RiBi1")]])
MET.SCORE.HYB.SC <- as.matrix(YSC$HYB.SC[[c("Glycolysis1", "OXPHOS1", "RiBi1")]])
MET.SCORE.HYB.SE <- as.matrix(YSC$HYB.SE[[c("Glycolysis1", "OXPHOS1", "RiBi1")]])

CC.SHARED.PARENT <- cell_cycle_continuum_shared(YSC$MIX.SC, YSC$MIX.SE)
CC.SHARED.HYBRID <- cell_cycle_continuum_shared(YSC$HYB.SC, YSC$HYB.SE)

COMP.BOUND.PARENT <- species_composition_report(CC.SHARED.PARENT$x1, CC.SHARED.PARENT$x2, MET.SCORE.MIX.SC, MET.SCORE.MIX.SE, COV.DIAG.MIX.SC, COV.DIAG.MIX.SE, "Sc vs Se parent")
COMP.BOUND.HYBRID <- species_composition_report(CC.SHARED.HYBRID$x1, CC.SHARED.HYBRID$x2, MET.SCORE.HYB.SC, MET.SCORE.HYB.SE, COV.DIAG.HYB.SC, COV.DIAG.HYB.SE, "Hybrid, Sc vs Se allele")

fig_pdf("extra/S_species_composition_bound.pdf", 8, 5)
par(mfrow = c(1, 2), mar = c(8, 4.5, 3, 1))
barplot(setNames(COMP.BOUND.PARENT$bound, COMP.BOUND.PARENT$axis), las = 2, col = COLOR.GREY[["mid"]], border = NA,
        ylab = "bound on expected noise shift (fraction of residual SD)", main = "Sc vs Se parent")
barplot(setNames(COMP.BOUND.HYBRID$bound, COMP.BOUND.HYBRID$axis), las = 2, col = COLOR.GREY[["mid"]], border = NA,
        ylab = "bound on expected noise shift (fraction of residual SD)", main = "Hybrid, Sc vs Se allele")
dev.off()

# Checkpoint
save(CC.AXIS.MIX.SC, CC.AXIS.MIX.SE, CC.AXIS.HYB.SC, CC.AXIS.HYB.SE,
     MET.STATE.MIX.SC, MET.STATE.MIX.SE, MET.STATE.HYB.SC, MET.STATE.HYB.SE,
     COV.DIAG.MIX.SC, COV.DIAG.MIX.SE, COV.DIAG.HYB.SC, COV.DIAG.HYB.SE,
     CC.SHARED.PARENT, CC.SHARED.HYBRID, COMP.BOUND.PARENT, COMP.BOUND.HYBRID,
     file = ckpt_path(7))

console_start(8)
##############################################################################
## 8. GO ENRICHMENT                                                         ##
##############################################################################
# Universe is BURST.CONTRASTS$gene throughout. GO.SIG is the
# permutation-p / class-call threshold used to build every gene set
# below.
GO.SIG  <- 0.05
GO.QVAL <- 0.2

## 8.1 Overall parental divergence: mean, burst frequency, burst size
# Baseline sets per quantity: any gene with a significant parent-vs-
# parent difference, split by direction
MEAN.LFC       <- setNames(BURST.CONTRASTS$mean_total_est, BURST.CONTRASTS$gene)
MEAN.ORD       <- sort(MEAN.LFC, decreasing = TRUE)
SIG.MEAN.TOTAL <- setNames(PR$mean_total_q, BURST.CONTRASTS$gene)[names(MEAN.ORD)]
MEAN.SC.UP     <- names(MEAN.ORD)[SIG.MEAN.TOTAL < GO.SIG & MEAN.ORD > 0]
MEAN.SE.UP     <- names(MEAN.ORD)[SIG.MEAN.TOTAL < GO.SIG & MEAN.ORD < 0]

BFREQ.LFC       <- setNames(BURST.CONTRASTS$bfreq_total_est, BURST.CONTRASTS$gene)
BFREQ.ORD       <- sort(BFREQ.LFC, decreasing = TRUE)
SIG.BFREQ.TOTAL <- setNames(PR$bfreq_total_q, BURST.CONTRASTS$gene)[names(BFREQ.ORD)]
BFREQ.SC.UP     <- names(BFREQ.ORD)[SIG.BFREQ.TOTAL < GO.SIG & BFREQ.ORD > 0]
BFREQ.SE.UP     <- names(BFREQ.ORD)[SIG.BFREQ.TOTAL < GO.SIG & BFREQ.ORD < 0]

BSIZE.LFC       <- setNames(BURST.CONTRASTS$bsize_total_est, BURST.CONTRASTS$gene)
BSIZE.ORD       <- sort(BSIZE.LFC, decreasing = TRUE)
SIG.BSIZE.TOTAL <- setNames(PR$bsize_total_q, BURST.CONTRASTS$gene)[names(BSIZE.ORD)]
BSIZE.SC.UP     <- names(BSIZE.ORD)[SIG.BSIZE.TOTAL < GO.SIG & BSIZE.ORD > 0]
BSIZE.SE.UP     <- names(BSIZE.ORD)[SIG.BSIZE.TOTAL < GO.SIG & BSIZE.ORD < 0]

GO.SETS <- list(Mean_Total_Sc  = MEAN.SC.UP,  Mean_Total_Se  = MEAN.SE.UP,
                Bfreq_Total_Sc = BFREQ.SC.UP, Bfreq_Total_Se = BFREQ.SE.UP,
                Bsize_Total_Sc = BSIZE.SC.UP, Bsize_Total_Se = BSIZE.SE.UP)

## 8.2 Regulatory classes, each split by direction: mean, burst
## frequency, burst size, frequency-size balance
# Twenty class sets: five classes (Conserved, Cis, Trans, Cis + Trans,
# Compensatory) x two directions (Sc-higher, Se-higher), for each of
# four quantities.
REG.SETS.MEAN  <- build_reg_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, quantity = "mean",  sig = GO.SIG)
REG.SETS.BFREQ <- build_reg_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, quantity = "bfreq", sig = GO.SIG)
REG.SETS.BSIZE <- build_reg_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, quantity = "bsize", sig = GO.SIG)
REG.SETS.KBAL  <- build_reg_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, quantity = "kbal",  sig = GO.SIG)
names(REG.SETS.MEAN)  <- paste0("Mean_",  names(REG.SETS.MEAN))
names(REG.SETS.BFREQ) <- paste0("Bfreq_", names(REG.SETS.BFREQ))
names(REG.SETS.BSIZE) <- paste0("Bsize_", names(REG.SETS.BSIZE))
names(REG.SETS.KBAL)  <- paste0("Kbal_",  names(REG.SETS.KBAL))
GO.SETS <- c(GO.SETS, REG.SETS.MEAN, REG.SETS.BFREQ, REG.SETS.BSIZE, REG.SETS.KBAL)

## 8.3 Any-cis / any-trans, pooled across class: mean, burst frequency,
## burst size, frequency-size balance
ANY.SETS.MEAN  <- c(build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "cis",   quantity = "mean",  sig = GO.SIG),
                    build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "trans", quantity = "mean",  sig = GO.SIG))
ANY.SETS.BFREQ <- c(build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "cis",   quantity = "bfreq", sig = GO.SIG),
                    build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "trans", quantity = "bfreq", sig = GO.SIG))
ANY.SETS.BSIZE <- c(build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "cis",   quantity = "bsize", sig = GO.SIG),
                    build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "trans", quantity = "bsize", sig = GO.SIG))
ANY.SETS.KBAL  <- c(build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "cis",   quantity = "kbal",  sig = GO.SIG),
                    build_component_go_sets(BURST.CONTRASTS, PR, universe = BURST.CONTRASTS$gene, component = "trans", quantity = "kbal",  sig = GO.SIG))
names(ANY.SETS.MEAN)  <- paste0("Mean_",  names(ANY.SETS.MEAN))
names(ANY.SETS.BFREQ) <- paste0("Bfreq_", names(ANY.SETS.BFREQ))
names(ANY.SETS.BSIZE) <- paste0("Bsize_", names(ANY.SETS.BSIZE))
names(ANY.SETS.KBAL)  <- paste0("Kbal_",  names(ANY.SETS.KBAL))
GO.SETS <- c(GO.SETS, ANY.SETS.MEAN, ANY.SETS.BFREQ, ANY.SETS.BSIZE, ANY.SETS.KBAL)

## 8.4 Run enrichment for every set
INTR.SETS       <- frac_group_sets(NOISE.DECOMP.NOAMBIG$gene,       NOISE.DECOMP.NOAMBIG$intrinsic_frac)
INTR.SETS.CLEAN <- frac_group_sets(NOISE.DECOMP.CLEAN.NOAMBIG$gene, NOISE.DECOMP.CLEAN.NOAMBIG$intrinsic_frac)
GO.INPUTS <- list(GO.SETS = GO.SETS, GO.UNIVERSE = BURST.CONTRASTS$gene,
                  INTR.SETS = INTR.SETS, INTR.UNIVERSE = NOISE.DECOMP.NOAMBIG$gene,
                  INTR.SETS.CLEAN = INTR.SETS.CLEAN, INTR.UNIVERSE.CLEAN = NOISE.DECOMP.CLEAN.NOAMBIG$gene,
                  GSE.LISTS = list(MEAN = MEAN.ORD, BFREQ = BFREQ.ORD, BSIZE = BSIZE.ORD),
                  GO.QVAL = GO.QVAL, SEED.GO = 1)
if (!exists("KEGG.DATA")) KEGG.DATA <- kegg_local("sce")

## ---- Cluster round trip: Rscript go_enrich.R ----
## Reads go_enrich_inputs.rda (saved below), writes go_enrich_output.rda
## (GO.ENRICH, INTR.GO, INTR.GO.CLEAN, GO.GSE, GO.KEY)
save(GO.INPUTS, KEGG.DATA, file = file.path(INPUT.DIR, "go_enrich_inputs.rda"))

load_cluster_output(file.path(OUTPUT.DIR, "go_enrich_output.rda"), "go_enrich.R")   # GO.ENRICH, INTR.GO, INTR.GO.CLEAN, GO.GSE, GO.KEY
check_cluster_key(GO.KEY, GO.INPUTS, "go_enrich_output.rda", "go_enrich.R")
## ---- end cluster round trip ----
cat(sprintf("GO/KEGG enrichment loaded for %d gene sets against a universe of %d genes\n", length(GO.ENRICH), length(BURST.CONTRASTS$gene)))


## 8.5 Summary table: significant terms per set, per ontology

## Summary table: one row per set, gene set size plus significant-term counts for BP/CC/MF/KEGG.
## GO.SETS and GO.ENRICH are the matched gene-set list and run_enrichment output (same names, same order).
GO.SUMMARY <- local({
  sets <- GO.SETS
  enrich <- GO.ENRICH
  q <- GO.QVAL
  data.frame(
    set     = names(sets),
    n_genes = vapply(sets, length, integer(1)),
    GO_BP   = vapply(lapply(enrich, `[[`, "BP"), n_sig_terms, integer(1), q = q),
    GO_CC   = vapply(lapply(enrich, `[[`, "CC"), n_sig_terms, integer(1), q = q),
    GO_MF   = vapply(lapply(enrich, `[[`, "MF"), n_sig_terms, integer(1), q = q),
    KEGG    = vapply(lapply(enrich, `[[`, "KEGG"), n_sig_terms, integer(1), q = q),
    row.names = NULL
  )
})
print(GO.SUMMARY)
write.csv(GO.SUMMARY, file.path(TABLE.DIR, "go_kegg_summary.csv"), row.names = FALSE)

## 8.6 Gene set enrichment (rank-based, no significance cutoff)
# GO.GSE holds GO.GSE.MEAN.BP through GO.GSE.BSIZE.CC from go_enrich.R
list2env(GO.GSE, envir = environment())

## 8.7 Intrinsic / extrinsic noise: GO enrichment by intrinsic fraction
# Low, average, and high intrinsic-fraction sets,
# enriched by go_enrich.R

sapply(INTR.GO,       function(s) sapply(s, n_sig_terms, q = GO.QVAL))
sapply(INTR.GO.CLEAN, function(s) sapply(s, n_sig_terms, q = GO.QVAL))

# Checkpoint
save(GO.SETS, GO.ENRICH, GO.SUMMARY,
     INTR.SETS, INTR.SETS.CLEAN, INTR.GO, INTR.GO.CLEAN,
     GO.GSE.MEAN.BP, GO.GSE.MEAN.MF, GO.GSE.MEAN.CC,
     GO.GSE.BFREQ.BP, GO.GSE.BFREQ.MF, GO.GSE.BFREQ.CC,
     GO.GSE.BSIZE.BP, GO.GSE.BSIZE.MF, GO.GSE.BSIZE.CC,
     file = ckpt_path(8))

console_start(9)
##############################################################################
## 9. EXTERNAL NOISE VALIDATION                                             ##
##############################################################################
# Checks the MIX.SC negative-binomial fit against seven independent
# published measurements of yeast expression noise: three protein-level
# flow cytometry studies (Newman 2006, Keren 2015, Stewart-Ornstein
# 2012) and four single-cell RNA-seq datasets (Gasch 2017,
# Nadal-Ribelles 2019, Jackson 2020, Jariani 2020). CONTRAST.FITS and
# GENES come from Section 2 above. The scRNA-seq sources are read,
# QC-filtered, and NB-fit here with the same offset model MIX.SC uses,
# so all eight datasets carry the same mean, CV^2, burst frequency, and
# burst size columns.

## 9.1 Newman 2006, Keren 2015, and Stewart-Ornstein 2012 datasets
# MIX.SC NB fit (mean and dispersion)
NB.SC <- data.frame(
  ORF  = GENES,
  MU   = as.numeric(CONTRAST.FITS$MIX.SC[GENES, "MU"]),
  DISP = as.numeric(CONTRAST.FITS$MIX.SC[GENES, "DISP"]),
  stringsAsFactors = FALSE
)

# Newman et al. 2006 (genome-wide GFP, YEPD)
# Supplementary Table 1. CV is reported as a percentage, so it is converted to a fraction
# and squared to land on the same CV^2 scale used everywhere else here.
NEWMAN.RAW <- read_excel(file.path(EXTERNAL.DIR, "Newman.2006.xls"),
                          sheet = "Sheet1", col_names = FALSE, skip = 4)

NEWMAN <- data.frame(
  ORF  = NEWMAN.RAW[[1]],
  Mean = as.numeric(NEWMAN.RAW[[4]]),
  CV2  = (as.numeric(NEWMAN.RAW[[24]]) / 100)^2,
  DM   = as.numeric(NEWMAN.RAW[[25]]),
  stringsAsFactors = FALSE
)
NEWMAN <- NEWMAN[!is.na(NEWMAN$ORF) & !is.na(NEWMAN$Mean) & !is.na(NEWMAN$CV2), ]

# Keren et al. 2015 (YFP promoter library, glucose plus amino acids)
# Processed supplemental file (post-gating, post-filtering, replicates
# united). Glucose plus amino acids is the closest match to standard rich-medium
# growth
KEREN.RAW <- read.delim(file.path(EXTERNAL.DIR, "Keren.2015.txt"), stringsAsFactors = FALSE)

KEREN <- data.frame(
  ORF  = KEREN.RAW$ORF,
  Mean = KEREN.RAW$MU_yfp_glu_plus_AA,
  CV2  = KEREN.RAW$CVsq_yfp_glu_plus_AA,
  Fano = KEREN.RAW$Fano_yfp_glu_plus_AA,
  stringsAsFactors = FALSE
)
KEREN <- KEREN[!is.na(KEREN$ORF) & !is.na(KEREN$Mean) & !is.na(KEREN$CV2), ]

# Stewart-Ornstein et al. 2012 (YFP fusions, total noise)
# Table S1. 
STEWART.RAW <- read_excel(file.path(EXTERNAL.DIR, "Stewart-Ornstein.2012.xls"), sheet = "TableS1")

STEWART <- data.frame(
  ORF  = gsub("'", "", STEWART.RAW$ORF),
  Mean = STEWART.RAW[["Mean Expression (AU)"]],
  CV2  = STEWART.RAW$Cvtot^2,
  stringsAsFactors = FALSE
)
STEWART <- STEWART[!is.na(STEWART$ORF) & !is.na(STEWART$Mean) & !is.na(STEWART$CV2), ]

## 9.2 NB-implied noise and per-source merge
# NB-implied CV^2 and Fano factor for MIX.SC
NB.SC$CV2   <- 1 / NB.SC$MU + 1 / NB.SC$DISP
NB.SC$Fano  <- 1 + NB.SC$MU / NB.SC$DISP
NB.SC$BFREQ <- NB.SC$DISP           # burst frequency, this pipeline's own DISP column
NB.SC$BSIZE <- NB.SC$MU / NB.SC$DISP  # burst size, equivalently Fano - 1

# Implied burst frequency and burst size for each external source, via
# add_burst_terms(). 
NEWMAN  <- add_burst_terms(NEWMAN)
KEREN   <- add_burst_terms(KEREN)
STEWART <- add_burst_terms(STEWART)

# Per-source join on ORF
EXT.MERGE <- list(
  Newman          = merge(NB.SC, NEWMAN,  by = "ORF"),
  Keren           = merge(NB.SC, KEREN,   by = "ORF"),
  StewartOrnstein = merge(NB.SC, STEWART, by = "ORF")
)

## 9.3 Published single-cell RNA-seq sources
# QC thresholds are set per source
GASCH.MIN.CELL.COUNT   <- 10000; GASCH.MIN.CELLS.EXPR   <- 5
NADAL.MIN.CELL.COUNT   <- 500;   NADAL.MIN.CELLS.EXPR   <- 5
JACKSON.MIN.CELL.COUNT <- 200;   JACKSON.MIN.CELLS.EXPR <- 20
JARIANI.MIN.CELL.COUNT <- 500;   JARIANI.MIN.CELLS.EXPR <- 10

# One PSOCK cluster for NB fits
N.CORES  <- 4
NOISE.CL <- parallel::makeCluster(N.CORES)

# Gasch et al. 2017 (Fluidigm C1, unstressed BY4741, ~80 cells)
# The Unstressed columns are selected from the header. Values
# are spike-in normalized rather than raw counts, so they are rounded to
# integers before fitting, and the ERCC spike-in rows are dropped. 
GASCH.FILE   <- file.path(EXTERNAL.DIR, "Gasch.2017.GSE102475_GASCH_NaCl-scRNAseq_NormData.txt")
GASCH.HEADER <- get_data_header(GASCH.FILE)
GASCH.KEEP   <- grepl("_Unstressed_", GASCH.HEADER)
GASCH.RAW    <- fread_matrix(GASCH.FILE, keep = GASCH.KEEP)
GASCH.RAW    <- GASCH.RAW[!grepl("^ERCC-", rownames(GASCH.RAW)), , drop = FALSE]
GASCH.MAT    <- to_numeric_matrix(GASCH.RAW, "Gasch")
GASCH.MAT    <- round(GASCH.MAT)
GASCH <- fit_scrna_source(GASCH.MAT, rownames(GASCH.MAT), "Gasch", GASCH.MIN.CELL.COUNT, GASCH.MIN.CELLS.EXPR, NOISE.CL)

# Nadal-Ribelles et al. 2019 (yscRNA-seq, unstressed BY4741, ~127 cells)
# Rows are TSS-level, so several rows can share one ORF and
# collapse_to_orf() sums them. 
NADAL.FILE   <- file.path(EXTERNAL.DIR, "Nadal-Ribelles.2019.41564_2018_346_MOESM3_ESM.txt")
NADAL.HEADER <- read_header_line(NADAL.FILE, skip = 1)
NADAL.RAW    <- fread(file = NADAL.FILE, sep = "\t", select = setdiff(NADAL.HEADER, "comGeneName"), skip = 1, header = TRUE)
NADAL.IDS    <- NADAL.RAW$geneName
NADAL.MAT    <- to_numeric_matrix(as.data.frame(NADAL.RAW[, setdiff(colnames(NADAL.RAW), "geneName"), with = FALSE]), "Nadal-Ribelles")
NADAL <- fit_scrna_source(NADAL.MAT, NADAL.IDS, "Nadal-Ribelles", NADAL.MIN.CELL.COUNT, NADAL.MIN.CELLS.EXPR, NOISE.CL)

# Jackson et al. 2020 (10x droplet, pooled genotypes, one condition, ~11,000 cells)
# Cells pool all 12
# genotypes (11 TF deletions plus wild type), so between-genotype
# differences in mean expression enter alongside cell-to-cell noise and
# these estimates are an upper bound.
JACKSON.COND   <- "YPD"
JACKSON.FILE   <- file.path(EXTERNAL.DIR, "Jackson.2020.GSE125162_ALL-fastqTomat0-Counts.tsv")
JACKSON.HEADER <- get_data_header(JACKSON.FILE)
JACKSON.HEADER.CLEAN <- trimws(gsub('^"|"$', "", JACKSON.HEADER))
JACKSON.ORF.KEEP <- !is.na(suppressMessages(map_to_orf(JACKSON.HEADER.CLEAN)))
if (!any(JACKSON.ORF.KEEP))
  message("No Jackson gene columns resolved to an ORF; first few raw column names: ",
          paste(head(JACKSON.HEADER, 5), collapse = ", "))

JACKSON.MAT <- fread_matrix(JACKSON.FILE, keep = JACKSON.ORF.KEEP)
message("Jackson conditions present: ",
        paste(unique(sub("^[0-9]+_", "", rownames(JACKSON.MAT))), collapse = ", "))
JACKSON.MAT <- JACKSON.MAT[grepl(paste0("_", JACKSON.COND, "$"), rownames(JACKSON.MAT)), , drop = FALSE]
JACKSON.MAT <- t(to_numeric_matrix(JACKSON.MAT, "Jackson"))
JACKSON <- fit_scrna_source(JACKSON.MAT, rownames(JACKSON.MAT), "Jackson", JACKSON.MIN.CELL.COUNT, JACKSON.MIN.CELLS.EXPR, NOISE.CL)

# Jariani et al. 2020 (10x droplet, continuous glucose growth, ~1,000 cells)
# Genes-by-cells, blank first header field, standard 10x barcode column
# names. Glucose 12h file
JARIANI.FILE <- file.path(EXTERNAL.DIR, "Jariani.2020.GSM4297055_processed_counts_glu_12h.txt")
JARIANI.MAT  <- fread_matrix(JARIANI.FILE)
JARIANI.MAT  <- to_numeric_matrix(JARIANI.MAT, "Jariani")
JARIANI <- fit_scrna_source(JARIANI.MAT, rownames(JARIANI.MAT), "Jariani", JARIANI.MIN.CELL.COUNT, JARIANI.MIN.CELLS.EXPR, NOISE.CL)

parallel::stopCluster(NOISE.CL)

NEW.SOURCES <- list(Gasch = GASCH, NadalRibelles = NADAL, Jackson = JACKSON, Jariani = JARIANI)
NEW.MERGE   <- lapply(NEW.SOURCES, merge, x = NB.SC, by = "ORF")

## 9.4 Correlations against MIX.SC and among external sources
# Mean-vs-mean and noise-vs-noise, each source against MIX.SC
EXT.CORR <- do.call(rbind, lapply(names(EXT.MERGE), noise_corr_rows, tables = EXT.MERGE, kind = "vs_mix"))
EXT.CORR <- EXT.CORR[, c("source", "statistic", "n", "rho", "p")]

# Newman's DM is a mean-corrected residual rather than a raw CV^2, so it
# is reported separately rather than folded into EXT.CORR above.
NEWMAN.DM.CORR <- cor_row(1 / EXT.MERGE$Newman$DISP, EXT.MERGE$Newman$DM)
NEWMAN.DM.CORR$source    <- "Newman"
NEWMAN.DM.CORR$statistic <- "DM vs 1/DISP"

print(EXT.CORR, row.names = FALSE)
print(NEWMAN.DM.CORR[, c("source", "statistic", "n", "rho", "p")], row.names = FALSE)

# Each single-cell RNA-seq source against MIX.SC
NEW.CORR <- do.call(rbind, lapply(names(NEW.MERGE), noise_corr_rows, tables = NEW.MERGE, kind = "vs_mix"))
NEW.CORR <- NEW.CORR[, c("source", "statistic", "n", "rho", "p")]
print(NEW.CORR, row.names = FALSE)

# Pairwise agreement among all seven published sources
EXT.RAW <- list(Newman = NEWMAN[, c("ORF", "Mean", "CV2", "BFREQ", "BSIZE")],
                Keren  = KEREN[,  c("ORF", "Mean", "CV2", "BFREQ", "BSIZE")],
                StewartOrnstein = STEWART[, c("ORF", "Mean", "CV2", "BFREQ", "BSIZE")])
NEW.RAW <- lapply(NEW.SOURCES, raw_source_table)
ALL.RAW <- c(EXT.RAW, NEW.RAW)

ALL.PAIRS <- combn(names(ALL.RAW), 2, simplify = FALSE)

ALL.CORR.PAIRWISE <- do.call(rbind, lapply(ALL.PAIRS, noise_corr_rows, tables = ALL.RAW, kind = "pairwise"))
ALL.CORR.PAIRWISE <- ALL.CORR.PAIRWISE[, c("comparison", "statistic", "n", "rho", "p")]
print(ALL.CORR.PAIRWISE, row.names = FALSE)

## 9.5 Mean-adjusted noise
# mean_adjusted_noise() expresses each gene's noise relative to genes of
# similar abundance in the same dataset, so the comparison below is of
# gene-specific noise with each dataset's own abundance trend removed.
NB.SC$CV2_ADJ <- mean_adjusted_noise(NB.SC$MU, NB.SC$CV2)
ALL.RAW <- lapply(ALL.RAW, add_cv2_adj)

# How much of each source's raw CV^2 variation follows from abundance
# alone.
MEAN.CV2.LINK <- sapply(c(list(MIX.SC = transform(NB.SC, Mean = MU)), ALL.RAW),
  function(d) cor(log(d$Mean), log(d$CV2), method = "spearman", use = "complete.obs"))
print(round(MEAN.CV2.LINK, 2))

# Agreement with MIX.SC in gene-specific noise, independent of agreement
# in expression level.
ADJ.CORR <- do.call(rbind, lapply(names(ALL.RAW), adj_corr_row, all_raw = ALL.RAW, nb_sc = NB.SC))
print(ADJ.CORR, row.names = FALSE)

## 9.6 Diagnostic plots
fig_pdf(file.path("extra", "external_noise_validation.pdf"), 8, 4, useDingbats = TRUE)
for (src in names(EXT.MERGE)) {
  d <- EXT.MERGE[[src]]
  par(mfrow = c(1, 2))
  plot(d$MU, d$Mean, log = "xy", pch = 19, cex = 0.5, col = "black",
       xlab = "MIX.SC mean (MU)", ylab = paste(src, "mean"),
       main = paste(src, "mean vs MIX.SC mean"))
  plot(d$CV2.x, d$CV2.y, log = "xy", pch = 19, cex = 0.5, col = "black",
       xlab = "MIX.SC implied CV^2", ylab = paste(src, "CV^2"),
       main = paste(src, "noise vs MIX.SC noise"))
}
dev.off()

# One page per single-cell RNA-seq source, showing the mean panel, the
# raw CV^2 panel, and the mean-adjusted panel together
fig_pdf(file.path("extra", "external_noise_validation_scrnaseq.pdf"), 12, 4, useDingbats = TRUE)
for (src in names(NEW.SOURCES)) {
  d <- merge(NB.SC, ALL.RAW[[src]], by = "ORF")
  par(mfrow = c(1, 3))

  plot(d$MU, d$Mean, log = "xy", pch = 19, cex = 0.5, col = "black",
       xlab = "MIX.SC mean (MU)", ylab = paste(src, "mean"),
       main = paste(src, "mean"))

  plot(d$CV2.x, d$CV2.y, log = "xy", pch = 19, cex = 0.5, col = "black",
       xlab = "MIX.SC CV^2", ylab = paste(src, "CV^2"),
       main = paste(src, "raw noise"))

  ok <- is.finite(d$CV2_ADJ.x) & is.finite(d$CV2_ADJ.y)
  plot(d$CV2_ADJ.x[ok], d$CV2_ADJ.y[ok], pch = 19, cex = 0.5, col = "black",
       xlab = "MIX.SC mean-adjusted log CV^2",
       ylab = paste(src, "mean-adjusted log CV^2"),
       main = paste(src, "mean-adjusted noise"))
  abline(h = 0, v = 0, col = COLOR.GREY[["mid"]], lty = 2)
}
dev.off()

save(NB.SC, EXT.MERGE, EXT.CORR, NEWMAN.DM.CORR, NEW.SOURCES, NEW.MERGE, NEW.CORR,
     ALL.RAW, ALL.CORR.PAIRWISE, MEAN.CV2.LINK, ADJ.CORR,
     file = ckpt_path(9))

cat(sprintf("figures written to figures/main, figures/extended, figures/extra  (%s)\n", date()))


console_start(10)
##############################################################################
## 10. POWER ANALYSIS                                                       ##
##############################################################################
# Run through 10.1 to build power_inputs.rda, then submit
# power.sub (a SLURM job array) to the cluster. Once every power_output_<k>of<K>.rda
# exists, run 10.2 onward to reshape them back into POWER and reproduce every plot.

set.seed(1)

# Significance level for power analysis. Power is judged at this Benjamini-Hochberg
# FDR (q < ALPHA), matching the gene-level analysis.
ALPHA <- 0.05
# Average number of reads per cell, before per-cell capture-depth scaling.
MEAN.READS <- c(0.25,0.5,1,2,4,8,16,32,64,128)
# Burst frequency (NB SIZE / dispersion parameter, theta) for the
# reference group.
SIZE <- c(0.5,1,2,4,8,16,32,64,128,256)
# Ratio of SIZE between the two groups. 1 confirms the false positive
# rate is controlled.
SIZE.RATIO <- c(1,1.1,1.25,1.5,1.75,2,3,4,6,8,12)

# Reference-group (Sc) cell count, evenly spaced.
N.CELLS <- c(400, 800, 1600, 2400)
# Se:Sc cell-count ratio. Set to 1, equal group sizes.
CELL.RATIO <- 1

# Per-cell capture-depth coefficient of variation.
EXPOSURE.CV <- 0.8

NJ <- 1000 # Simulated datasets per row
NI <- N.PERM # Permutations per dataset, matched to the gene-level analysis
PI1 <- if (exists("PR")) min(0.5, max(0.02, mean(PR$bfreq_total_q < ALPHA, na.rm = TRUE))) else 0.1
cat(sprintf("PI1 (fraction of genes called by BH at q < %.2f, parental burst frequency) = %.3f\n", ALPHA, PI1))
N.MIX <- 50 # Random mixtures averaged per grid cell
SEED.BASE <- 1

## 10.1 Grid definition and cluster submission
# Each row covers one (MEAN.READS, N.CELLS, SIZE) combination and returns
# power for every SIZE.RATIO value at once
M <- 1:length(MEAN.READS); N <- 1:length(N.CELLS)
P <- 1:length(SIZE);       Q <- 1:length(SIZE.RATIO)

GRID <- expand.grid(m = M, n = N, p = P)

## ---- Cluster round trip: SLURM job array power.sub (Rscript power_grid.R, one task per array index) ----
## Reads power_inputs.rda (saved below), writes one power_output_<k>of<K>.rda
## per array task (ROWS.PART, POWER.PART), reshaped into POWER in 10.2
save(GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE,
     file = file.path(INPUT.DIR, "power_inputs.rda"))

cat(sprintf("power grid inputs written: %d rows x %d SIZE.RATIO values\n", nrow(GRID), length(SIZE.RATIO)))
cat("Submit power.sub (job array) now; resume below once every power_output_<k>of<K>.rda exists.\n")
## ---- end cluster round trip (power_output_<k>of<K>.rda files are loaded in 10.2) ----

## 10.2 Reshape cluster output into the POWER array
# Each array task writes the GRID rows it computed (ROWS.PART) and their power values
# (POWER.PART). Rows are placed back in GRID order, so the result does not depend on
# how many array tasks were used.
PART.FILES <- list.files(OUTPUT.DIR, pattern = "^power_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
stopifnot(length(PART.FILES) > 0)
POWER.MAT <- matrix(NA_real_, nrow(GRID), length(SIZE.RATIO))
ROW.DONE  <- logical(nrow(GRID))
for (f in PART.FILES) {
  e <- new.env(); load(f, envir = e)
  POWER.MAT[e$ROWS.PART, ] <- e$POWER.PART
  ROW.DONE[e$ROWS.PART]    <- TRUE
}
stopifnot(all(ROW.DONE))   # every array task must have finished before reshaping

POWER <- array(NA_real_, dim = c(length(M), length(N), length(P), length(Q)))
for (i in seq_len(nrow(GRID))) {
  POWER[GRID$m[i], GRID$n[i], GRID$p[i], ] <- POWER.MAT[i, ]
}
save(POWER, file = file.path(OUTPUT.DIR, "power.rda"))

## Value of each POWER dimension: mean reads (M), cell count (N), burst-frequency SIZE (P), SIZE.RATIO (Q)
POWER.GRID <- list(M = MEAN.READS, N = N.CELLS, P = SIZE, Q = SIZE.RATIO)
LOG2.READS.AXIS <- list(at = c(-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), labels = -2:10, mgp = c(3, 0.4, 0))

## 10.3 Figure 1: power vs mean expression, panels over N.CELLS x SIZE.RATIO
power_line_figure("Power.Analysis.grid_mean_by_ncells.pdf", POWER, POWER.GRID, "N", "M", x = "Q", series = "P",
  main_fmt = "Mean=%2$s Sc N.Cells=%1$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  xlab = "SIZE Ratio", ticks = TRUE)
power_line_figure("Power.Analysis.1.pdf", POWER, POWER.GRID, "N", "Q", x = "M", series = "P",
  main_fmt = "SIZE Ratio = %2$s\nSc # Cells = %1$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  styled = c(LOG2.READS.AXIS, list(cex_axis = 0.6, xlab = expression('log'[2]*'(Average Reads per Cell)'), xlab_mgp = c(1.2, 1, 0), xlab_cex = 0.6)))

## 10.4 Supplementary diagnostic grids: every other pairing of the four
# grid axes (mean, burst frequency, SIZE.RATIO, N.CELLS)
power_line_figure("Power.Analysis.grid_ncells_by_mean.pdf", POWER, POWER.GRID, "Q", "M", x = "N", series = "P",
  main_fmt = "Ratio=%1$s Mean=%2$s", legend_labels = paste0("SIZE = ", SIZE), legend_title = "Burst frequency",
  log_x = FALSE, xlab = "Sc Cell Count")
power_line_figure("Power.Analysis.grid_logmean_by_size_ratio.pdf", POWER, POWER.GRID, "P", "Q", x = "M", series = "N",
  main_fmt = "SIZE=%1$s Ratio=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = "log2(Mean Reads)")
power_line_figure("Power.Analysis.grid_size_by_ratio.pdf", POWER, POWER.GRID, "M", "Q", x = "P", series = "N",
  main_fmt = "MEAN=%1$s Ratio=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = expression('log'[2]*'(SIZE)'))
power_line_figure("Power.Analysis.grid_sizeratio_by_mean.pdf", POWER, POWER.GRID, "P", "M", x = "Q", series = "N",
  main_fmt = "SIZE=%1$s Mean=%2$s", legend_labels = paste0("Sc Cells = ", N.CELLS), legend_title = "Sc Cell Count",
  xlab = "SIZE.RATIO", ticks = TRUE)
power_line_figure("Power.Analysis.grid_sizeratio_by_ncells_full.pdf", POWER, POWER.GRID, "P", "N", x = "Q", series = "M",
  main_fmt = "SIZE=%1$s Sc Cell Number=%2$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  xlab = "SIZE.RATIO", ticks = TRUE)
power_line_figure("Power.Analysis.grid_ncells_by_size_ratio_full.pdf", POWER, POWER.GRID, "P", "Q", x = "N", series = "M",
  main_fmt = "SIZE=%1$s Ratio=%2$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  log_x = FALSE, xlab = "Sc Cell Count")

## 10.5 Figure 2: power vs burst frequency, panels over N.CELLS x SIZE.RATIO
power_line_figure("Power.Analysis.2.pdf", POWER, POWER.GRID, "N", "Q", x = "P", series = "M", mar = c(4.5, 3, 2, 1),
  main_fmt = "SIZE Ratio = %2$s\nSc # Cells = %1$s", legend_labels = paste0("Mean = ", MEAN.READS), legend_title = "Reads/cell",
  styled = list(at = log2(SIZE), labels = SIZE, cex_axis = 0.55, mgp = c(3, 0.7, 0), las = 2, tcl = -0.3,
                xlab = "Burst Frequency (SIZE)", xlab_mgp = c(3.3, 1, 0), xlab_cex = 0.6))
power_line_figure("Power.Analysis.grid_size_by_ncells_mean.pdf", POWER, POWER.GRID, "N", "M", x = "P", series = "Q",
  main_fmt = "Sc Cell Number=%1$s Mean=%2$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  xlab = expression('log'[2]*'(SIZE)'))
power_line_figure("Power.Analysis.grid_ncells_by_size_mean.pdf", POWER, POWER.GRID, "P", "M", x = "N", series = "Q",
  main_fmt = "SIZE=%1$s Mean=%2$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  log_x = FALSE, xlab = "Sc Cell Count")

## 10.6 Figure 3: power vs mean expression, panels over N.CELLS x SIZE
power_line_figure("Power.Analysis.3.pdf", POWER, POWER.GRID, "N", "P", x = "M", series = "Q",
  main_fmt = "SIZE = %2$s\nSc # Cells = %1$s", legend_labels = paste0("Ratio = ", SIZE.RATIO), legend_title = "SIZE ratio",
  styled = c(LOG2.READS.AXIS, list(cex_axis = 0.55, xlab = expression('log'[2]*'(Average Reads per Cell)'), xlab_mgp = c(1.4, 1, 0), xlab_cex = 0.8)))

## 10.7 Heatmap A: power over mean expression x cell count, faceted by burst frequency
HM.Q <- which.min(abs(SIZE.RATIO - 2))
HM.COLS <- COLOR.SEQ

open_grid_pdf(file.path(FIGURE.DIR, "extended", "Power.Analysis.heatmap_mean_by_ncells.pdf"), ceiling(length(P)/2), 2, panel_w = 3.4, panel_h = 3.2)
for (p in P) {
  mat <- POWER[,,p,HM.Q]  # length(M) x length(N)
  image(log2(MEAN.READS), N.CELLS, mat, zlim = c(0,1), col = HM.COLS,
        xaxt = "n", yaxt = "n", xlab = "Mean Reads per Cell", ylab = "Sc Cell Count",
        main = paste0("SIZE = ", SIZE[p]), cex.main = 0.85)
  axis(1, at = log2(MEAN.READS), labels = MEAN.READS, las = 2, cex.axis = 0.55)
  axis(2, at = N.CELLS, labels = N.CELLS, cex.axis = 0.65)
  if (any(mat >= 0.8, na.rm = TRUE) && any(mat < 0.8, na.rm = TRUE)) {
    contour(log2(MEAN.READS), N.CELLS, mat, levels = 0.8, add = TRUE, col = "white", lwd = 1.5, labcex = 0.6)
  }
}
dev.off()

fig_pdf(file.path("extended", "Power.Analysis.heatmap_colorbar.pdf"), 4, 1.3, useDingbats = TRUE)
par(mar = c(2.2,1,1,1))
image(seq(0,1,length.out=100), 1, matrix(seq(0,1,length.out=100), ncol=1), col = HM.COLS,
      axes = FALSE, xlab = paste0("Power to detect SIZE.RATIO = ", round(SIZE.RATIO[HM.Q],2)), ylab = "")
axis(1, at = c(0,0.2,0.4,0.6,0.8,1.0), labels = c("0%","20%","40%","60%","80%","100%"))
dev.off()

## 10.8 Heatmap B: minimum detectable SIZE.RATIO, faceted by burst frequency

## Minimum detectable SIZE.RATIO: the smallest SIZE.RATIO (ratios in ascending order) at which power
## reaches 80%, one number per (MEAN.READS, N.CELLS, SIZE) in place of the full power-vs-ratio curve.
## Lower means more sensitive. NA when even the largest tested ratio misses the target, which image()
## renders as blank.
MDR <- array(NA_real_, dim = c(length(M), length(N), length(P)))
for (m in M) for (n in N) for (p in P) MDR[m,n,p] <- local({
  power_vec <- POWER[m, n, p, ]
  ratios <- SIZE.RATIO
  target <- 0.8
  hit <- which(power_vec >= target)
  if (length(hit) == 0) return(NA_real_)
  ratios[min(hit)]
})

MDR.COLS <- rev(COLOR.SEQ)
MDR.ZLIM <- log2(range(SIZE.RATIO))

open_grid_pdf(file.path(FIGURE.DIR, "extended", "Power.Analysis.min_detectable_ratio.pdf"), ceiling(length(P)/2), 2, panel_w = 3.4, panel_h = 3.2)
for (p in P) {
  mat <- log2(MDR[,,p])  # NA stays NA
  image(log2(MEAN.READS), N.CELLS, mat, zlim = MDR.ZLIM, col = MDR.COLS,
        xaxt = "n", yaxt = "n", xlab = "Mean Reads per Cell", ylab = "Sc Cell Count",
        main = paste0("SIZE = ", SIZE[p], "\n(blank = not detectable by ratio=12)"), cex.main = 0.7)
  axis(1, at = log2(MEAN.READS), labels = MEAN.READS, las = 2, cex.axis = 0.55)
  axis(2, at = N.CELLS, labels = N.CELLS, cex.axis = 0.65)
}
dev.off()

fig_pdf(file.path("extended", "Power.Analysis.min_detectable_ratio_colorbar.pdf"), 4, 1.3, useDingbats = TRUE)
par(mar = c(2.2,1,1,1))
zseq <- seq(MDR.ZLIM[1], MDR.ZLIM[2], length.out = 100)
image(zseq, 1, matrix(zseq, ncol=1), col = MDR.COLS,
      axes = FALSE, xlab = "Minimum detectable SIZE.RATIO (80% power)", ylab = "")
RATIO.TICKS <- SIZE.RATIO
axis(1, at = log2(RATIO.TICKS), labels = RATIO.TICKS, cex.axis = 0.7)
dev.off()

# Checkpoint
save(GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, PI1, N.MIX, SEED.BASE,
     POWER, MDR, file = ckpt_path(10))

console_stop()