###############################################################
### quicklook_angles.R
### Quick-look angle classification of every gene from point estimates (handoff "Claude Code handoff: quick-look angle classification", claude.ai artifact
### FUpo4Xgrw8WjhvmQzCfQhp). A first look at the data, not the final analysis: model-based standard errors and delta-method intervals stand in for the
### bootstrap, all hybrid cells serve every contrast (no cis/trans partition), and the Hc and He estimates are taken as uncorrelated.
###   Quantities  mean mu, dispersion phi, burst frequency alpha = 1 / phi, burst size beta = mu phi, kinetic balance kappa = alpha / beta, all on log2.
###   Schemes     cis/trans (theta_R), inheritance (theta_P) and trans inheritance (theta_T): each gene is a point (h, v) in an oriented plane, classed by the
###               angle theta = atan2(v, h) of that point and its delta-method interval against the landmark lines.
### Needs   : the Section 1 checkpoint (MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA). The shared shift Xi (DOM.SHIFT) is read from the Section 2 checkpoint when
###           that exists and QL.USE.XI is TRUE; otherwise the script runs without it and says so.
### Writes  : analysis/output/quicklook_angles.csv (gene x quantity x scheme), quicklook_class_counts.csv, quicklook_angles.pdf and quicklook_inputs.txt.
### Changes nothing in the pipeline; its functions are in R/functions_quicklook.R. Runs locally in minutes.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_quicklook.R"))
library(parallel); library(ggplot2); library(patchwork)
set.seed(SEED.SECTION + 50)

QL.OUT.DIR <- here("analysis", "output"); dir.create(QL.OUT.DIR, recursive = TRUE, showWarnings = FALSE)
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
if (!exists("QL.SIG"))       QL.SIG       <- 0.05    # BH level of the Conserved gate
if (!exists("QL.MIN.COUNT")) QL.MIN.COUNT <- 0       # mean count per cell a gene needs in every group (0 keeps every gene)
if (!exists("QL.USE.XI"))    QL.USE.XI    <- TRUE    # divide the hybrid allele dispersions by 2^Xi when the Section 2 checkpoint holds DOM.SHIFT

## ---- Inputs: the observed counts of the four groups and the split variances of the hybrid cells ----
stopifnot("quicklook_angles.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE, SPLIT.ETA
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)),
          identical(dim(HYB.SC), dim(HYB.SE)), length(SPLIT.ETA$Hc) == ncol(HYB.SC), length(SPLIT.ETA$He) == ncol(HYB.SE))

## Offsets as in Section 2 (Methods: Library split correction): a parent's offset is its library over the median parental library D, and each hybrid allele is
## offset by the library of its own genome, L_Hc / D and L_He / D, so every mean is a per-genome rate on the parents' scale. The split variance of each hybrid
## cell (SPLIT.ETA) sets the split correction of the allele dispersions.
QL.D     <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
QL.MATS  <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
QL.EXPOS <- lapply(QL.MATS, function(m) colSums(m) / QL.D)
QL.ETA   <- list(Pc = 0, Pe = 0, Hc = SPLIT.ETA$Hc, He = SPLIT.ETA$He)
stopifnot(all(vapply(QL.EXPOS, function(e) all(e > 0), logical(1))))
QL.GENES <- rownames(MIX.SC)[Reduce(`&`, lapply(QL.MATS, function(m) rowMeans(m) >= QL.MIN.COUNT))]

## ---- Fits: the pipeline's estimator in each group, with the variances of log2 mu and log2 phi from the fit alone ----
## mu is the ratio estimate and phi the Cox-Reid fit with the boundary rule, the interior check and, for the alleles, the split correction (fit_dispersion()).
## All hybrid cells enter each allele fit.
t0 <- Sys.time()
QL.FIT <- lapply(setNames(names(QL.MATS), names(QL.MATS)), function(k)
  do.call(rbind, mclapply(QL.GENES, ql_fit_gene, mat = QL.MATS[[k]], expo = QL.EXPOS[[k]], eta = QL.ETA[[k]], mc.cores = CORES)))
for (k in names(QL.FIT)) rownames(QL.FIT[[k]]) <- QL.GENES
cat(sprintf("Fitted %d genes in 4 groups in %.1f min\n", length(QL.GENES), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
QL.STATUS.TAB <- sapply(QL.FIT, function(f) table(factor(QL.STATUS[f[, "status"]], levels = QL.STATUS)))
cat("Fit status by group:\n"); print(QL.STATUS.TAB)

## ---- The shared shift Xi (optional): the hybrid's genome-wide rise in dispersion, removed from the hybrid alleles before the contrasts ----
## Dividing each hybrid allele phi by 2^Xi moves the log2 allele values by -Xi for phi, +Xi for alpha, -Xi for beta and +2 Xi for kappa; the mean carries no
## shift. It cancels in cis and trans and moves Delta_dom. It rescales phi and its standard error alike, so the variances stay as they are.
QL.XI <- setNames(rep(0, length(QL.GENES)), QL.GENES); QL.XI.NOTE <- "Xi not applied: no Section 2 checkpoint with DOM.SHIFT (or QL.USE.XI is FALSE)"
if (QL.USE.XI && file.exists(ckpt_path(2))) local({
  e <- new.env(); load(ckpt_path(2), envir = e)
  if (exists("DOM.SHIFT", envir = e, inherits = FALSE)) {
    hit <- intersect(QL.GENES, names(e$DOM.SHIFT))
    QL.XI[hit] <<- e$DOM.SHIFT[hit]
    QL.XI.NOTE <<- sprintf("Xi applied from the Section 2 checkpoint (DOM.SHIFT) to %d of %d genes; the rest run without it", length(hit), length(QL.GENES))
  }
})
cat(QL.XI.NOTE, "\n")

## ---- Per-group log2 values and their variances for the five quantities ----
## log2 alpha = -log2 phi, log2 beta = log2 mu + log2 phi, log2 kappa = -log2 mu - 2 log2 phi. The NB mean and dispersion estimates are uncorrelated on the same
## cells, so the variances add: V(alpha) = V(phi), V(beta) = V(mu) + V(phi), V(kappa) = V(mu) + 4 V(phi).
QL.GRP <- c("Pc", "Pe", "Hc", "He")
QL.COL <- lapply(setNames(colnames(QL.FIT$Pc), colnames(QL.FIT$Pc)), function(col) sapply(QL.GRP, function(k) QL.FIT[[k]][, col]))
LM <- suppressWarnings(log2(QL.COL$mu)); LP <- suppressWarnings(log2(QL.COL$phi)); VM <- QL.COL$v_mu; VP <- QL.COL$v_phi; ST <- QL.COL$status
XI.H <- cbind(Pc = 0, Pe = 0, Hc = QL.XI, He = QL.XI)
QL.Q <- list(mean = list(x = LM, v = VM, xi = 0), phi = list(x = LP, v = VP, xi = -1), alpha = list(x = -LP, v = VP, xi = 1),
             beta = list(x = LM + LP, v = VM + VP, xi = -1), kappa = list(x = -LM - 2 * LP, v = VM + 4 * VP, xi = 2))
## A quantity is undefined for a gene with no counts in a group (every quantity) or with a failed or zero dispersion in a group (every quantity but the mean):
## its log2 contrasts are infinite or missing, so it is flagged and takes no angle.
QL.BAD <- list(
  mean = ifelse(rowSums(ST == 4) > 0, "No counts", NA_character_),
  disp = ifelse(rowSums(ST == 4) > 0, "No counts", ifelse(rowSums(ST == 3) > 0, "Failed fit", ifelse(rowSums(ST == 2) > 0, "Zero fit", NA_character_))))

## ---- Contrasts, Wald tests and the three planes for every quantity ----
QL.TABLE <- do.call(rbind, lapply(names(QL.Q), function(q) {
  qq <- QL.Q[[q]]; x <- qq$x + qq$xi * XI.H; bad <- if (q == "mean") QL.BAD$mean else QL.BAD$disp
  x[!is.na(bad), ] <- NA; colnames(x) <- QL.GRP; colnames(qq$v) <- QL.GRP; rownames(x) <- QL.GENES
  ct <- ql_contrasts(x, qq$v)
  cbind(quantity = q, do.call(rbind, lapply(names(QL.SCHEMES), ql_scheme_table, ct = ct, bad = bad, sig = QL.SIG)), xi = unname(QL.XI[rownames(ct)]) * (qq$xi != 0))
}))
write.csv(QL.TABLE, file.path(QL.OUT.DIR, "quicklook_angles.csv"), row.names = FALSE)

## ---- Class counts per quantity and scheme, the cis + trans classes shown split (by which effect is larger) and merged ----
QL.COUNTS <- do.call(rbind, lapply(c("class", "class_merged"), function(col) {
  tb <- as.data.frame(table(quantity = QL.TABLE$quantity, scheme = QL.TABLE$scheme, class = QL.TABLE[[col]], useNA = "ifany"), stringsAsFactors = FALSE)
  tb <- tb[tb$Freq > 0, ]; tb$version <- if (col == "class") "split" else "merged"; tb
}))
QL.COUNTS <- QL.COUNTS[!(QL.COUNTS$version == "merged" & QL.COUNTS$scheme != "cis_trans"), ]
QL.COUNTS <- QL.COUNTS[order(match(QL.COUNTS$quantity, names(QL.Q)), QL.COUNTS$scheme, QL.COUNTS$version, -QL.COUNTS$Freq), c("quantity", "scheme", "version", "class", "Freq")]
names(QL.COUNTS)[5] <- "genes"
cat("\nGenes per class:\n"); print(QL.COUNTS, row.names = FALSE)
write.csv(QL.COUNTS, file.path(QL.OUT.DIR, "quicklook_class_counts.csv"), row.names = FALSE)

## ---- Figure: theta against magnitude, one column per scheme and one row per quantity, landmarks dashed, Conserved genes in grey ----
QL.PLOT <- QL.TABLE[is.finite(QL.TABLE$theta) & QL.TABLE$fit_flag == "", ]
QL.PLOT$quantity <- factor(QL.PLOT$quantity, levels = names(QL.Q))
QL.TITLE <- c(cis_trans = "cis / trans (theta_R)", inheritance = "Inheritance (theta_P)", trans_inh = "trans inheritance (theta_T)")
QL.PANELS <- lapply(names(QL.SCHEMES), function(sc) {
  d <- QL.PLOT[QL.PLOT$scheme == sc, ]
  lv <- sort(setdiff(unique(d$class), c("Conserved", "No trans divergence")))
  cols <- c(setNames(scales::hue_pal()(length(lv)), lv), Conserved = "grey75", `No trans divergence` = "grey75")
  lm <- QL.SCHEMES[[sc]]$landmarks; lm <- if (QL.SCHEMES[[sc]]$circular) c(lm, -90) else lm
  ggplot(d, aes(theta, magnitude, colour = class)) +
    geom_vline(xintercept = lm, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.6) + facet_grid(quantity ~ ., scales = "free_y") +
    scale_x_continuous(limits = c(-90, 90), breaks = c(-90, -45, 0, 45, 90)) + scale_colour_manual(values = cols, name = NULL) +
    labs(title = QL.TITLE[[sc]], x = "theta (degrees)", y = "magnitude (log2)") + theme_bw(base_size = 8) +
    theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm")) + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.5, alpha = 1)))
})
ggsave(file.path(QL.OUT.DIR, "quicklook_angles.pdf"), wrap_plots(QL.PANELS, nrow = 1), width = 14, height = 13, useDingbats = FALSE)

## ---- What was used ----
QL.INPUTS <- c(sprintf("Section 1 checkpoint: %s", ckpt_path(1)),
               "Counts: MIX.SC (Pc), MIX.SE (Pe), HYB.SC (Hc), HYB.SE (He), all hybrid cells for every contrast",
               sprintf("Offsets: library / D, D = %.1f (median parental library); alleles on their own genome's library; split variances SPLIT.ETA", QL.D),
               "Fits: fit_dispersion(method = 'NB'): ratio mean, Cox-Reid phi with the boundary rule and interior check, split-corrected for Hc and He",
               sprintf("Genes: %d (mean count per cell >= %g in every group)", length(QL.GENES), QL.MIN.COUNT),
               QL.XI.NOTE,
               sprintf("Gate: both contrasts of a scheme with BH q >= %g (Wald tests, BH within contrast x quantity) -> Conserved", QL.SIG),
               "Intervals: theta +- 1.96 SE (delta method); Hc-He covariance set to zero")
writeLines(QL.INPUTS, file.path(QL.OUT.DIR, "quicklook_inputs.txt"))
cat("\n", paste(QL.INPUTS, collapse = "\n"), "\n", sep = "")
