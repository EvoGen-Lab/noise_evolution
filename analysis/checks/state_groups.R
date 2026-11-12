###############################################################
### state_groups.R
### ST1: cell numbers in each of 24 candidate cell states, and where those states sit in the transcriptome. Runs from the Section 1 checkpoint alone.
### Three state axes, each level marked by about 10 genes chosen before any noise result (gene symbols mapped to ORFs through org.Sc.sgd.db):
###   metabolism  Fermentative (glycolysis and ethanol production), Respiratory (TCA cycle and respiratory chain), Derepressed (glyoxylate cycle,
###               fatty-acid beta-oxidation and gluconeogenesis, the non-fermentative program of the trans-leak cells)
###   cycle       G1S (SBF/MBF wave), S (histone genes), G2M (CLB2 cluster), MG1 (SIC1 cluster with the Ace2 daughter genes)
###   stress      Stressed (induced arm of the environmental stress response) or Unstressed
### Each module's score is the log2 of its summed counts per cell over the cell's depth factor (both alleles summed in the hybrid cells, so the state is
### the cell's). Scores are centered at each group's median (STATE.CENTER "group") or at the median of all cells ("pooled") and scaled by the pooled
### median absolute deviation, so a cell's robust z says how far it sits above the typical cell of its group. Assignment (state_assign()): a cell is
### Derepressed, else Respiratory, when that z exceeds STATE.Z, else Fermentative; its phase is the highest cycle z; it is Stressed when the stress z
### exceeds STATE.Z. Group centering compares each cell with its own culture, so a state that a whole culture shifts into is not counted as a rare
### state; pooled centering keeps such shifts but also carries the species differences in the marker genes' expression.
### Outputs: the 24-state counts and percentages per group (Pc, Pe, hybrid H) with the pooled frequency and the count it predicts in the smallest
### group, the counts of the rare levels at other thresholds (STATE.Z.SWEEP), each module's genes in the data and its mean count per cell, and the
### per-cell labels. PCA per group on the log-normalized counts of the most variable genes, the marker genes left out (PCA.EXCLUDE.MARKERS), so any
### separation of the states in the PCA is evidence from genes that did not define them; the share of each PC's variance that each axis explains
### (R squared). Points: color metabolism, shape cycle phase, closed symbols stressed and open symbols unstressed.
### Results go to results/checks/ST1.rds, results/checks/ST1_state_groups.csv and a row in summary.csv; the figure to
### figures/extra/check_ST1_state_pca.pdf.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R"))
library(Matrix); library(org.Sc.sgd.db); library(irlba)

CHECK.SCRIPT        <- "analysis/checks/state_groups.R"
STATE.Z             <- 2.5                 # robust z above which a cell is in a rare level (Derepressed, Respiratory, Stressed)
STATE.Z.SWEEP       <- c(2, 2.5, 3, 4)     # thresholds for the sensitivity table of the rare levels
STATE.CENTER        <- "group"             # "group" or "pooled" (see the header)
STATE.PSEUDO        <- 0.5                 # pseudocount of the module sums, so a cell with no module counts has a finite score
PCA.GENES           <- 2000                # most variable genes per group in the PCA
PCA.PCS             <- 10
PCA.EXCLUDE.MARKERS <- TRUE
set.seed(SEED.SECTION + 26)

STATE.SYMBOLS <- list(
  metabolism = list(Fermentative = c("HXK2", "PFK1", "FBA1", "TPI1", "TDH3", "PGK1", "ENO2", "CDC19", "PDC1", "ADH1"),
                    Respiratory  = c("CIT1", "ACO1", "KGD1", "SDH1", "CYT1", "QCR2", "COX4", "COX5A", "COX6", "ATP2"),
                    Derepressed  = c("ICL1", "MLS1", "CIT2", "MDH2", "ACS1", "FOX2", "POX1", "POT1", "PCK1", "FBP1")),
  cycle      = list(G1S = c("CLN1", "CLN2", "PCL1", "CLB5", "CLB6", "RNR1", "CDC21", "POL30", "RAD27", "RFA1"),
                    S   = c("HTA1", "HTA2", "HTB1", "HTB2", "HHT1", "HHT2", "HHF1", "HHF2", "HHO1"),
                    G2M = c("CLB1", "CLB2", "CDC5", "CDC20", "SWI5", "ACE2", "BUD4", "CHS2", "MOB1", "KIP2"),
                    MG1 = c("SIC1", "PCL9", "ASH1", "CDC6", "EGT2", "CTS1", "DSE1", "DSE2", "SCW11", "AMN1")),
  stress     = list(Stressed = c("HSP12", "HSP26", "HSP42", "HSP104", "CTT1", "TPS1", "TPS2", "GSY2", "DDR2", "GRE2")))
MET.LEVELS   <- c("Fermentative", "Respiratory", "Derepressed")
CYCLE.LEVELS <- names(STATE.SYMBOLS$cycle)
LEVEL.SYMBOLS <- unlist(unname(STATE.SYMBOLS), recursive = FALSE)   # one symbol vector per level, named by level

stopifnot("state_groups.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE (cells after QC, genes after the union filter, ortholog rows named by Sc ORF)
stopifnot(identical(rownames(HYB.SC), rownames(HYB.SE)), identical(colnames(HYB.SC), colnames(HYB.SE)))

## Gene symbols to ORFs. Every symbol is listed with its ORF and whether the ortholog row is in the data, so missing markers (no 1:1 ortholog, or
## removed by the gene filter) are visible.
SYMBOL.MAP <- data.frame(level = rep(names(LEVEL.SYMBOLS), lengths(LEVEL.SYMBOLS)), symbol = unlist(LEVEL.SYMBOLS), row.names = NULL)
SYMBOL.MAP$orf <- suppressMessages(AnnotationDbi::mapIds(org.Sc.sgd.db, keys = SYMBOL.MAP$symbol, keytype = "GENENAME", column = "ORF", multiVals = "first"))
SYMBOL.MAP$in_data <- SYMBOL.MAP$orf %in% rownames(MIX.SC)
MEMBERS <- lapply(setNames(names(LEVEL.SYMBOLS), names(LEVEL.SYMBOLS)), function(l) unique(SYMBOL.MAP$orf[SYMBOL.MAP$level == l & SYMBOL.MAP$in_data]))
print(SYMBOL.MAP[!SYMBOL.MAP$in_data, ], row.names = FALSE)
cat("Marker genes in the data per level:", paste(sprintf("%s %d", names(MEMBERS), lengths(MEMBERS)), collapse = ", "), "\n")
stopifnot("every state level needs at least 3 marker genes in the data" = all(lengths(MEMBERS) >= 3))

## Module sums and scores per group. The hybrid cell carries both alleles, summed; its depth factor is its full hybrid depth over the parental median.
DS        <- list(Pc = MIX.SC, Pe = MIX.SE, H = HYB.SC + HYB.SE)
DEPTH     <- lapply(DS, colSums)
DEPTH.REF <- median(c(DEPTH$Pc, DEPTH$Pe))
MODTOT    <- lapply(DS, function(x) vapply(MEMBERS, function(m) colSums(x[m, , drop = FALSE]), numeric(ncol(x))))
for (k in names(DS)) rownames(MODTOT[[k]]) <- colnames(DS[[k]])
SCORE     <- Map(function(t, d) log2((t + STATE.PSEUDO) / (d / DEPTH.REF)), MODTOT, DEPTH)

## Robust z: centered per group (or pooled), scaled by the pooled median absolute deviation of the centered scores.
CENTER <- if (STATE.CENTER == "group") lapply(SCORE, function(s) apply(s, 2, median)) else rep(list(apply(do.call(rbind, SCORE), 2, median)), length(SCORE))
CENTERED <- Map(function(s, c) sweep(s, 2, c), SCORE, CENTER)
SCALE <- apply(do.call(rbind, CENTERED), 2, mad)
stopifnot("a module score has zero spread; it cannot separate cells" = all(SCALE > 0))
Z <- lapply(CENTERED, function(s) sweep(s, 2, SCALE, "/"))

## How well each module is measured: its mean summed count per cell (the counting noise of the score is about 1 / that count) and the share of cells
## with no count at all.
MODULE.N <- do.call(rbind, lapply(names(DS), function(k) data.frame(group = k, level = names(MEMBERS), genes = lengths(MEMBERS), mean_count = colMeans(MODTOT[[k]]),
                                                                      pct_zero = 100 * colMeans(MODTOT[[k]] == 0), row.names = NULL)))
print(MODULE.N, digits = 3, row.names = FALSE)
cat("Spearman correlation of the module z-scores, all cells:\n")
print(round(cor(do.call(rbind, Z), method = "spearman"), 2))

## Labels at STATE.Z and the 24-state table: counts and percentages per group, the pooled frequency, and the count that frequency predicts in the
## smallest group (the expected count that the removal rule reads).
LAB  <- lapply(Z, state_assign, thr = STATE.Z, cycle_levels = CYCLE.LEVELS)
GRID <- expand.grid(stress = c("Unstressed", "Stressed"), cycle = CYCLE.LEVELS, metabolism = MET.LEVELS, stringsAsFactors = FALSE)[, 3:1]
COUNTS <- GRID
for (k in names(LAB)) {
  a <- LAB[[k]]
  n <- vapply(seq_len(nrow(GRID)), function(i) sum(a$metabolism == GRID$metabolism[i] & a$cycle == GRID$cycle[i] & a$stress == GRID$stress[i]), numeric(1))
  COUNTS[[paste0(k, "_n")]] <- n
  COUNTS[[paste0(k, "_pct")]] <- 100 * n / nrow(a)
}
N.CELLS <- vapply(LAB, nrow, numeric(1))
COUNTS$pooled_pct   <- 100 * rowSums(COUNTS[paste0(names(LAB), "_n")]) / sum(N.CELLS)
COUNTS$min_expected <- COUNTS$pooled_pct / 100 * min(N.CELLS)
cat(sprintf("\nCells per group: %s. Threshold z > %.1f, centering %s.\n", paste(sprintf("%s %d", names(N.CELLS), N.CELLS), collapse = ", "), STATE.Z, STATE.CENTER))
print(COUNTS, digits = 3, row.names = FALSE)
for (ax in c("metabolism", "cycle", "stress")) {
  cat(sprintf("\n%s, cells per group (%% of group):\n", ax))
  tab <- sapply(LAB, function(a) table(a[[ax]]))
  print(cbind(tab, round(100 * sweep(tab, 2, N.CELLS, "/"), 2)))
}

## The rare levels at other thresholds.
SWEEP <- list()
for (thr in STATE.Z.SWEEP) for (k in names(Z)) {
  a <- state_assign(Z[[k]], thr = thr, cycle_levels = CYCLE.LEVELS)
  SWEEP[[paste(thr, k)]] <- data.frame(z = thr, group = k, cells = nrow(a), Respiratory = sum(a$metabolism == "Respiratory"),
                                       Derepressed = sum(a$metabolism == "Derepressed"), Stressed = sum(a$stress == "Stressed"),
                                       states_below_5 = sum(table(interaction(a$metabolism, a$cycle, a$stress)) < 5), row.names = NULL)
}
SWEEP <- do.call(rbind, SWEEP)
cat("\nRare levels by threshold:\n"); print(SWEEP, row.names = FALSE)

##############################################################################
## PCA per group on log-normalized counts (each cell scaled to the parental median depth, log1p), the PCA.GENES most variable genes, each gene scaled to
## unit variance, so the total variance equals the number of genes and each PC's share is its variance over that number. With the marker genes left
## out, the states separate in the PCA only through other genes that co-vary with them.
PCA <- list(); PCA.R2 <- list()
MARKERS <- unique(unlist(MEMBERS))
for (k in names(DS)) {
  x <- as(DS[[k]], "CsparseMatrix")
  if (PCA.EXCLUDE.MARKERS) x <- x[setdiff(rownames(x), MARKERS), , drop = FALSE]
  nx <- x %*% Diagonal(x = DEPTH.REF / DEPTH[[k]])
  nx@x <- log1p(nx@x)
  m <- rowMeans(nx); v <- rowMeans(nx^2) - m^2
  top <- head(order(v, decreasing = TRUE), PCA.GENES)
  top <- top[v[top] > 0]
  pc <- with_local_seed(SEED.SECTION + 26, prcomp_irlba(t(nx[top, ]), n = PCA.PCS, center = TRUE, scale. = TRUE))
  PCA[[k]] <- list(x = pc$x, pct = 100 * pc$sdev^2 / length(top), genes = rownames(nx)[top])
  a <- LAB[[k]]; a$state24 <- interaction(a$metabolism, a$cycle, a$stress, drop = TRUE)
  PCA.R2[[k]] <- data.frame(group = k, pc = seq_len(PCA.PCS), pct_var = PCA[[k]]$pct,
                            vapply(c("metabolism", "cycle", "stress", "state24"), function(ax) apply(pc$x, 2, function(p) if (nlevels(droplevels(a[[ax]])) > 1) summary(lm(p ~ a[[ax]]))$r.squared else NA_real_),
                                   numeric(PCA.PCS)), row.names = NULL)
}
PCA.R2 <- do.call(rbind, PCA.R2)
cat("\nShare of each PC's variance explained by each axis (R squared); pct_var is the PC's share of the total:\n")
print(PCA.R2, digits = 3, row.names = FALSE)
## Variance of the top PCs explained by each axis, as a percentage of the total over the PCA genes (sum over PCs of pct_var x R squared).
PCA.AXIS <- do.call(rbind, lapply(split(PCA.R2, PCA.R2$group), function(d) data.frame(group = d$group[1], pcs = PCA.PCS, pct_var_pcs = sum(d$pct_var),
                                                                                       t(colSums(d$pct_var * d[c("metabolism", "cycle", "stress", "state24")], na.rm = TRUE)), row.names = NULL)))
cat("\nPercent of total variance (over the PCA genes) in the top PCs that each axis explains:\n"); print(PCA.AXIS, digits = 3, row.names = FALSE)

MET.COL      <- c(Fermentative = "#0072B2", Respiratory = "#009E73", Derepressed = "#D55E00")
PCH.OPEN     <- setNames(c(1, 0, 2, 5), CYCLE.LEVELS)       # unstressed
PCH.CLOSED   <- setNames(c(16, 15, 17, 18), CYCLE.LEVELS)   # stressed
GROUP.LABELS <- c(Pc = "Sc parent", Pe = "Se parent", H = "Hybrid (both alleles)")
fig_pdf("extra/check_ST1_state_pca.pdf", 10, 15)
layout(matrix(c(1:6, 7, 7), ncol = 2, byrow = TRUE), heights = c(1, 1, 1, 0.3))
par(mar = c(4.5, 4.5, 3, 1))
for (k in names(PCA)) {
  a <- LAB[[k]]; p <- PCA[[k]]
  o <- order(a$stress == "Stressed", a$metabolism != "Fermentative")   # common cells first, rare cells drawn on top
  pch <- ifelse(a$stress == "Stressed", PCH.CLOSED[as.character(a$cycle)], PCH.OPEN[as.character(a$cycle)])
  col <- ifelse(a$metabolism == "Fermentative" & a$stress == "Unstressed", setNames(adjustcolor(MET.COL, 0.35), names(MET.COL))[as.character(a$metabolism)], setNames(adjustcolor(MET.COL, 0.9), names(MET.COL))[as.character(a$metabolism)])
  for (pr in list(c(1, 2), c(3, 4))) {
    plot(p$x[o, pr[1]], p$x[o, pr[2]], pch = pch[o], col = col[o], cex = 0.5,
         xlab = sprintf("PC%d (%.1f%%)", pr[1], p$pct[pr[1]]), ylab = sprintf("PC%d (%.1f%%)", pr[2], p$pct[pr[2]]),
         main = sprintf("%s, %d cells%s", GROUP.LABELS[[k]], nrow(a), if (PCA.EXCLUDE.MARKERS) ", markers excluded" else ""))
  }
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", ncol = 4, bty = "n", cex = 1,
       legend = c(names(MET.COL), paste(CYCLE.LEVELS, "unstressed"), paste(CYCLE.LEVELS, "stressed")),
       col = c(MET.COL, rep(COLOR.ACCENT, 2 * length(CYCLE.LEVELS))), pch = c(rep(15, 3), PCH.OPEN, PCH.CLOSED))
dev.off()

dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(COUNTS, file.path(CHECK.DIR, "ST1_state_groups.csv"), row.names = FALSE)
check_save("ST1", list(counts = COUNTS, sweep = SWEEP, module_n = MODULE.N, symbols = SYMBOL.MAP, members = MEMBERS, labels = LAB, z = Z,
                       pca_r2 = PCA.R2, pca_axis = PCA.AXIS, pca = PCA,
                       settings = list(z = STATE.Z, center = STATE.CENTER, pseudo = STATE.PSEUDO, pca_genes = PCA.GENES, exclude_markers = PCA.EXCLUDE.MARKERS)), "run",
           headline = c(sprintf("z > %.1f, %s centering; cells %s", STATE.Z, STATE.CENTER, paste(sprintf("%s %d", names(N.CELLS), N.CELLS), collapse = ", ")),
                        sprintf("%s: Respiratory %d, Derepressed %d, Stressed %d; states under 5 cells %d of 24",
                                SWEEP$group[SWEEP$z == STATE.Z], SWEEP$Respiratory[SWEEP$z == STATE.Z], SWEEP$Derepressed[SWEEP$z == STATE.Z],
                                SWEEP$Stressed[SWEEP$z == STATE.Z], SWEEP$states_below_5[SWEEP$z == STATE.Z]),
                        sprintf("%s: top %d PCs hold %.1f%% of variance; metabolism %.2f, cycle %.2f, stress %.2f, 24 states %.2f",
                                PCA.AXIS$group, PCA.AXIS$pcs, PCA.AXIS$pct_var_pcs, PCA.AXIS$metabolism, PCA.AXIS$cycle, PCA.AXIS$stress, PCA.AXIS$state24)),
           mc_error = "none (deterministic, apart from the truncated PCA)", script = CHECK.SCRIPT)
