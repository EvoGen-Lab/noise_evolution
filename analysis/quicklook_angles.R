###############################################################
### quicklook_angles.R
### Quick-look classification by direction of every gene from point estimates (spec "Classification by direction", claude.ai artifact Y9Vn84qUaxXJVf8ZebEJEr;
### it replaces the angle scheme of the first quick-look handoff). A first look at the data, not the final analysis: model-based standard errors and normal
### intervals stand in for the bootstrap, all hybrid cells serve every part (no split, so A and B covary through the shared hybrid cells), and the Hc and He
### estimates are taken as uncorrelated. The gate and the classes are the pipeline's (dir_classify_table()).
###   Quantities  mean mu, dispersion phi, burst frequency alpha = 1 / phi, burst size beta = mu phi, kinetic balance kappa = alpha / beta, all on log2.
###   Schemes     cis/trans, inheritance and trans inheritance: parts A and B from Pe toward Pc, diagonals S and D; relevance and equivalence tests against a
###               margin gate each gene; a divergent gene is classed by its direction theta = atan2(B, A) mod 180 and the interval of theta against bands of
###               +-DIR.W degrees around the class landmarks.
### Needs   : the Section 1 checkpoint (MIX.SC, MIX.SE, HYB.SC, HYB.SE; the split variances are computed from them). The shared shift Xi (DOM.SHIFT) is read from the Section 2 checkpoint when
###           that exists and QL.USE.XI is TRUE; otherwise the script runs without it and says so.
### Writes  : analysis/output/quicklook_angles.csv (gene x quantity x scheme), quicklook_class_counts.csv, quicklook_sensitivity.csv, quicklook_angles.pdf and
###           quicklook_inputs.txt.
### Changes nothing in the pipeline; its functions are in R/functions_quicklook.R. Runs locally in minutes.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_quicklook.R"))
library(parallel); library(ggplot2); library(patchwork)
set.seed(SEED.SECTION + 50)

QL.OUT.DIR <- here("analysis", "output"); dir.create(QL.OUT.DIR, recursive = TRUE, showWarnings = FALSE)
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)   # mclapply forks, so Windows runs on one core
if (!exists("QL.SIG"))       QL.SIG       <- 0.05    # BH level of the gate
if (!exists("QL.MIN.COUNT")) QL.MIN.COUNT <- 0       # mean count per cell a gene needs in every group (0 keeps every gene)
if (!exists("QL.USE.XI"))    QL.USE.XI    <- TRUE    # divide the hybrid allele dispersions by 2^Xi when the Section 2 checkpoint holds DOM.SHIFT

## ---- Inputs: the observed counts of the four groups and the split variances of the hybrid cells ----
stopifnot("quicklook_angles.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))      # MIX.SC, MIX.SE, HYB.SC, HYB.SE
stopifnot(identical(rownames(MIX.SC), rownames(HYB.SC)), identical(rownames(MIX.SE), rownames(MIX.SC)), identical(rownames(HYB.SE), rownames(HYB.SC)),
          identical(dim(HYB.SC), dim(HYB.SE)))
## The split variance of each hybrid cell follows from its counts alone, as in Section 1.2: the observed share s = L_Hc / L_H and the library L_H give
## eta_Hc = (1 - s) / (s L_H) and eta_He = s / ((1 - s) L_H) (split_eta()). Computing it here from the loaded matrices keeps it aligned with their columns and
## serves checkpoints saved before SPLIT.ETA was added to them.
QL.LIB.H  <- colSums(HYB.SC) + colSums(HYB.SE)
SPLIT.ETA <- split_eta(colSums(HYB.SC) / QL.LIB.H, QL.LIB.H)

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

## ---- The shared shift Xi (optional): the hybrid's genome-wide rise in noise, used for the inheritance schemes centred on it ----
## Centring Delta_dom is an open item, so both versions are reported: inheritance and trans_inh as measured, and inheritance_xi and trans_inh_xi with each hybrid
## log2 value moved by -Xi for phi and beta, +Xi for alpha and +2 Xi for kappa (the mean carries no shift). Xi cancels in cis and trans.
QL.XI <- setNames(rep(0, length(QL.GENES)), QL.GENES); QL.XI.NOTE <- "Xi not available (no Section 2 checkpoint with DOM.SHIFT, or QL.USE.XI is FALSE): the centred inheritance schemes are not reported"
if (QL.USE.XI && file.exists(ckpt_path(2))) local({
  e <- new.env(); load(ckpt_path(2), envir = e)
  if (exists("DOM.SHIFT", envir = e, inherits = FALSE)) {
    hit <- intersect(QL.GENES, names(e$DOM.SHIFT))
    QL.XI[hit] <<- e$DOM.SHIFT[hit]
    QL.XI.NOTE <<- sprintf("Xi from the Section 2 checkpoint (DOM.SHIFT) for %d of %d genes, used for the centred inheritance schemes (the rest centred on 0)", length(hit), length(QL.GENES))
  }
})
cat(QL.XI.NOTE, "\n")

## ---- Per-group log2 values and their variances for the five quantities ----
## log2 alpha = -log2 phi, log2 beta = log2 mu + log2 phi, log2 kappa = -log2 mu - 2 log2 phi. The NB mean and dispersion estimates are uncorrelated on the same
## cells, so the variances add: V(alpha) = V(phi), V(beta) = V(mu) + V(phi), V(kappa) = V(mu) + 4 V(phi).
QL.GRP <- c("Pc", "Pe", "Hc", "He")
QL.COL <- lapply(setNames(colnames(QL.FIT$Pc), colnames(QL.FIT$Pc)), function(col) sapply(QL.GRP, function(k) QL.FIT[[k]][, col]))
LM <- suppressWarnings(log2(QL.COL$mu)); LP <- suppressWarnings(log2(QL.COL$phi)); VM <- QL.COL$v_mu; VP <- QL.COL$v_phi; ST <- QL.COL$status
QL.Q <- list(mean = list(x = LM, v = VM), disp = list(x = LP, v = VP), alpha = list(x = -LP, v = VP), bsize = list(x = LM + LP, v = VM + VP), kbal = list(x = -LM - 2 * LP, v = VM + 4 * VP))
for (q in names(QL.Q)) { colnames(QL.Q[[q]]$x) <- QL.GRP; colnames(QL.Q[[q]]$v) <- QL.GRP; rownames(QL.Q[[q]]$x) <- QL.GENES }
## A gene with no counts in a group has no estimate for any quantity; a failed dispersion fit leaves every quantity but the mean undefined (failed fit); a zero
## dispersion makes the log2 parts infinite (zero fit: tested, but given no direction).
QL.FAILED <- list(mean = rowSums(ST == 4) > 0, disp = rowSums(ST == 4 | ST == 3) > 0)
QL.ZERO   <- list(mean = rep(FALSE, length(QL.GENES)), disp = rowSums(ST == 2) > 0 & !QL.FAILED$disp)

## ---- Margins: delta_phi = log2(1.15), and for every other quantity the same multiple of its median SE of S ----
QL.DELTA <- local({
  se <- vapply(names(QL.Q), function(q) median(sqrt(QL.Q[[q]]$v[, "Pc"] + QL.Q[[q]]$v[, "Pe"])[!QL.FAILED[[if (q == "mean") "mean" else "disp"]] & !QL.ZERO[[if (q == "mean") "mean" else "disp"]]], na.rm = TRUE), numeric(1))
  DIR.DELTA.PHI * se / se[["disp"]]
})
cat("Margins delta (log2), scaled by the median SE of S:\n"); print(round(QL.DELTA, 4))

## ---- Rows of the direction table, the gate, and the classes ----
## Pass 1 builds every row with 95% theta intervals and gates the genes; the false coverage-statement rate of each quantity and scheme, alpha_FCR = 0.05 N_div /
## N_test, then sets the level 1 - 2 alpha_FCR of the theta intervals of pass 2, which gives the classes. The centred inheritance schemes (_xi) need Xi.
QL.SCH <- c(.DIR.SCHEMES, if (any(QL.XI != 0)) c("inheritance_xi", "trans_inh_xi"))
QL.COMBOS <- expand.grid(quantity = names(QL.Q), scheme = QL.SCH, stringsAsFactors = FALSE)
QL.COMBOS <- QL.COMBOS[!(QL.COMBOS$quantity == "mean" & grepl("_xi$", QL.COMBOS$scheme)), ]
QL.TAB1 <- do.call(rbind, Map(ql_dir_combo, QL.COMBOS$quantity, QL.COMBOS$scheme, MoreArgs = list(Q = QL.Q, delta = QL.DELTA, zero = QL.ZERO, failed = QL.FAILED, xi = QL.XI)))
QL.PRE <- dir_classify_table(QL.TAB1, sig = QL.SIG)
QL.ALPHA <- tapply(QL.PRE$alpha_fcr, paste(QL.PRE$quantity, QL.PRE$scheme, sep = "|"), function(a) a[1])
QL.TAB <- do.call(rbind, Map(ql_dir_combo, QL.COMBOS$quantity, QL.COMBOS$scheme, MoreArgs = list(Q = QL.Q, delta = QL.DELTA, zero = QL.ZERO, failed = QL.FAILED, xi = QL.XI, alpha = QL.ALPHA)))
QL.TABLE <- dir_classify_table(QL.TAB, sig = QL.SIG)
write.csv(QL.TABLE, file.path(QL.OUT.DIR, "quicklook_angles.csv"), row.names = FALSE)

## ---- Class counts at the main margin and band, and the sensitivity of the counts to the margin multiplier and the band width ----
QL.COUNTS <- as.data.frame(table(quantity = QL.TABLE$quantity, scheme = QL.TABLE$scheme, class = QL.TABLE$class_label), stringsAsFactors = FALSE)
QL.COUNTS <- QL.COUNTS[QL.COUNTS$Freq > 0, ]; names(QL.COUNTS)[4] <- "genes"
QL.COUNTS <- QL.COUNTS[order(match(QL.COUNTS$quantity, names(QL.Q)), QL.COUNTS$scheme, -QL.COUNTS$genes), ]
cat("\nGenes per class (margin multiplier 1, w = ", DIR.W, " degrees):\n", sep = ""); print(QL.COUNTS, row.names = FALSE)
write.csv(QL.COUNTS, file.path(QL.OUT.DIR, "quicklook_class_counts.csv"), row.names = FALSE)
QL.SENS <- do.call(rbind, lapply(seq_along(DIR.DELTA.MULT), function(m) do.call(rbind, lapply(DIR.W.SENS, function(w) {
  ct <- dir_classify_table(QL.TAB, m = m, w = w, sig = QL.SIG); tb <- as.data.frame(table(quantity = ct$quantity, scheme = ct$scheme, class = ct$class), stringsAsFactors = FALSE)
  cbind(delta_mult = DIR.DELTA.MULT[m], w = w, tb[tb$Freq > 0, ]) }))))
write.csv(QL.SENS, file.path(QL.OUT.DIR, "quicklook_sensitivity.csv"), row.names = FALSE)

## ---- Figure: direction theta (0 to 180 degrees) against the size of (A, B), one column per scheme and one row per quantity; class bands shaded, the
## boundary landmark dashed, divergent genes coloured by class and the rest in grey ----
QL.PLOT <- QL.TABLE[is.finite(QL.TABLE$theta), ]
QL.PLOT$quantity <- factor(QL.PLOT$quantity, levels = names(QL.Q)); QL.PLOT$magnitude <- sqrt(QL.PLOT$A^2 + QL.PLOT$B^2)
QL.PLOT$colour <- ifelse(QL.PLOT$gate %in% "Divergent", QL.PLOT$class, "Not divergent")
QL.TITLE <- c(cis_trans = "cis / trans", inheritance = "Inheritance", trans_inh = "trans inheritance", inheritance_xi = "Inheritance, centred on Xi", trans_inh_xi = "trans inheritance, centred on Xi")
QL.PANELS <- lapply(intersect(names(QL.TITLE), unique(QL.PLOT$scheme)), function(sc) {
  d <- QL.PLOT[QL.PLOT$scheme == sc, ]; seg <- dir_segments(sub("_xi$", "", sc), DIR.W); bands <- seg[seg$band, ]
  lv <- sort(setdiff(unique(d$colour), "Not divergent"))
  cols <- c(setNames(scales::hue_pal()(length(lv)), lv), `Not divergent` = "grey80")
  bd <- rbind(bands, transform(bands[bands$from < 0, ], from = from + 180, to = to + 180))
  ggplot(d, aes(theta, magnitude, colour = colour)) +
    geom_rect(data = bd, aes(xmin = pmax(from, 0), xmax = pmin(to, 180), ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey92") +
    geom_vline(xintercept = if (sub("_xi$", "", sc) == "cis_trans") 45 else 135, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.6) + facet_grid(quantity ~ ., scales = "free_y") +
    scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, 45)) + scale_colour_manual(values = cols, name = NULL) +
    labs(title = QL.TITLE[[sc]], x = "direction theta (degrees)", y = "size of (A, B), log2") + theme_bw(base_size = 8) +
    theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm")) + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.5, alpha = 1)))
})
QL.FIGURE <- wrap_plots(QL.PANELS, nrow = 1)
ggsave(file.path(QL.OUT.DIR, "quicklook_angles.pdf"), QL.FIGURE, width = 3.2 * length(QL.PANELS) + 1, height = 13, useDingbats = FALSE)
cat(sprintf("Figure written to %s\n", normalizePath(file.path(QL.OUT.DIR, "quicklook_angles.pdf"), winslash = "/")))
## In an interactive session the figure is also drawn on the current device (the RStudio Plots pane); enlarge the pane or use Zoom.
if (interactive()) print(QL.FIGURE)

## ---- What was used ----
QL.INPUTS <- c(sprintf("Section 1 checkpoint: %s", ckpt_path(1)),
               "Counts: MIX.SC (Pc), MIX.SE (Pe), HYB.SC (Hc), HYB.SE (He), all hybrid cells for every part (no split)",
               sprintf("Offsets: library / D, D = %.1f (median parental library); alleles on their own genome's library; split variances from the hybrid counts (split_eta())", QL.D),
               "Fits: fit_dispersion(method = 'NB'): ratio mean, Cox-Reid phi with the boundary rule and interior check, split-corrected for Hc and He",
               sprintf("Genes: %d (mean count per cell >= %g in every group)", length(QL.GENES), QL.MIN.COUNT),
               QL.XI.NOTE,
               sprintf("Margins (log2): %s; sensitivity at multipliers %s", paste(sprintf("%s %.3f", names(QL.DELTA), QL.DELTA), collapse = ", "), paste(DIR.DELTA.MULT, collapse = ", ")),
               sprintf("Gate: relevance and equivalence tests (normal approximation), BH within quantity x scheme x contrast x test at %g", QL.SIG),
               sprintf("Classes: theta intervals at 1 - 2 alpha_FCR (delta method), bands +-%g degrees (sensitivity %s); Hc-He covariance set to zero", DIR.W, paste(DIR.W.SENS, collapse = ", ")))
writeLines(QL.INPUTS, file.path(QL.OUT.DIR, "quicklook_inputs.txt"))
cat("\n", paste(QL.INPUTS, collapse = "\n"), "\n", sep = "")
