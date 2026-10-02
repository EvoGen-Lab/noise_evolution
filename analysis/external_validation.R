###############################################################
### external_validation.R
### Section 9 of the pipeline: comparison of the MIX.SC noise
### estimates with published yeast noise data (three protein-level
### studies, four single-cell RNA-seq studies).
###
### Run after analysis.R has written the Section 2 checkpoint
### (section2_checkpoint.rda, which supplies GENES and CONTRAST.FITS).
### Functions live in R/functions_external.R (this analysis) and
### R/functions.R (shared NB fit, session helpers).
### Writes the Section 9 checkpoint and the validation figures.
###############################################################

# here() anchors every path to the project root, so R/setup.R loads from any working
# directory. setup.R checks renv, defines the project paths, class names and color
# palettes, and sources R/functions.R.
library(here)
source(here("R", "setup.R"))
source(here("R", "functions_external.R"))

suppressWarnings(suppressPackageStartupMessages({
  library(MASS)
  library(parallel)
  library(readxl)
  library(data.table)
}))

# GENES (the genes with contrasts) and CONTRAST.FITS (per-dataset NB mean and
# dispersion) come from the Section 2 checkpoint
SEC2 <- new.env()
load(ckpt_path(2), envir = SEC2)
stopifnot(all(c("GENES", "CONTRAST.FITS") %in% ls(SEC2)))
GENES         <- SEC2$GENES
CONTRAST.FITS <- SEC2$CONTRAST.FITS
rm(SEC2)

console_start(9)
set.seed(SEED.SECTION + 9)
##############################################################################
## 9. EXTERNAL NOISE VALIDATION                                             ##
##############################################################################
# Checks the MIX.SC negative-binomial fit against seven independent
# published measurements of yeast expression noise: three protein-level
# flow cytometry studies (Newman 2006, Keren 2015, Stewart-Ornstein
# 2012) and four single-cell RNA-seq datasets (Gasch 2017,
# Nadal-Ribelles 2019, Jackson 2020, Jariani 2020). CONTRAST.FITS and
# GENES come from the Section 2 checkpoint loaded above. The scRNA-seq sources are read,
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
# add_burst_terms(). The three protein sources report fluorescence in instrument units, so they use
# scale = "relative": BFREQ is NA and BSIZE (the implied Fano factor) carries rank information only,
# which leaves CV^2 and the mean-adjusted noise (9.5) as their unit-free comparisons.
NEWMAN  <- add_burst_terms(NEWMAN,  scale = "relative")
KEREN   <- add_burst_terms(KEREN,   scale = "relative")
STEWART <- add_burst_terms(STEWART, scale = "relative")

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
# The cluster is stopped whether or not a source fails to load or fit
tryCatch({
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
}, finally = parallel::stopCluster(NOISE.CL))

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

cat(sprintf("figures written to figures/extra  (%s)\n", date()))
