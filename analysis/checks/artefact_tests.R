###############################################################
### artefact_tests.R
### Tests of why the data depart from the NB (handoff "tests of why the data depart from the NB", T0 to T12), run from the Section 1 checkpoint.
### Each test separates candidate causes of four findings of estimator_leak.R: Cox-Reid and moment estimates differ in the parents' highest-phi bins,
### trimming the top cells cuts phi in the hybrid alleles, a variance term linear in the mean is positive and grows with the mean, and Delta phi and
### Delta mu are negatively coupled across genes. The hypotheses are H1 (a few low-quality cells), H2 (an additive component outside the offset), H3
### (non-proportional scaling with depth), H4 (expression jackpots in real cell states), H5 (biological coupling) and H6 (outlier cells behind the
### estimator differences).
###
### Definitions, fixed before any result (checked against estimator_leak.R; see results/checks/T0_audit.md):
###   trimming    the top ceiling(0.005 n) cells of a fit by scaled count y / O are left out (5 of 834 Pc, 7 of 1240 Pe, 15 of 2804 hybrid cells);
###   own slope   the OLS slope of Delta log2 phi (Cox-Reid) on Delta log2 mu across genes; T8, T9 and T12 keep genes whose groups all have phi > 0 and
###               dispersion class 1 (the pipeline's gate for log2 contrasts);
###   leak slope  the slope of [Delta log2 phi (Cox-Reid) - Delta log2 phi (moment)] on Delta log2 mu, equal to the Cox-Reid own slope minus the moment
###               own slope on the same genes;
###   a           the coefficient of A_i = M_i (1 - p_i) in E[(y - M)^2 - A_i] = a A_i + phi B_i (linear_excess_gene()), a variance term linear in the mean.
### Rules fixed before any result:
###   T2  a cell is removed when its burden z = (observed - null mean) / null SD exceeds the (1 - T2.FWER) quantile, over the T0 replicates, of the largest
###       z among the group's cells in a replicate, with T2.FWER = 0.01: a familywise 1% rate per group that the null itself calibrates, because the
###       burden is discrete and skewed (a normal tail removed 4 to 7 cells per group from pure NB data); a hybrid cell flagged in either allele leaves both;
###   T7a every UMI is kept with probability 0.5; T7b parental UMIs are kept with the ratio of the median hybrid allele library to the median parental
###       library (UMIs, because reads are not in the project);
###   T10 gene sets from GO through org.Sc.sgd.db: cytosolic ribosome GO:0022626, ribosome biogenesis GO:0042254, nucleosome GO:0000786, cellular response
###       to heat GO:0034605; random sets matched on mean-count bin, 1000 permutations.
###   T11 the 20 cells with the largest scaled counts per gene and group, on at most T11.MAX.GENES genes per group.
### Simulations: T0 (gene-matched NB null, ART.R replicates), T8 (scenarios A and B) and T12 (S1 to S5) run in cluster/scripts/artefact_null_grid.R
### (cluster/submit/artefact_null.sub) from artefact_inputs.rda, written here; their outputs are read from OUTPUT.DIR. With ART.LOCAL = TRUE they run here
### on ART.LOCAL.GENES genes with ART.LOCAL.R replicates instead, and the affected results are marked "reduced".
### Results: results/checks/T0 ... T12 .rds, rows in summary.csv, results/checks/T0_audit.md and REPORT_artefacts.md; figures in figures/extra/.
###############################################################
library(here); source(here("R", "setup.R")); source(here("R", "functions_checks.R")); source(here("R", "functions_artefacts.R"))
library(parallel)

CHECK.SCRIPT    <- "analysis/checks/artefact_tests.R"
CORES           <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 1)
SEED.ART        <- 7026
ART.FRAC        <- 0.005    # trimming fraction
ART.R           <- 200      # T0 replicates per gene and group
ART.R4          <- 10       # T0 replicates that also get the T4 model fits
ART.R11         <- 5        # T0 replicates that also get the T11 influence
ART.THIN        <- 0.5      # T7a keep probability
ART.N.CHUNKS    <- 48       # T0 gene chunks of the array
ART.SCEN.REPS   <- 3        # replicates of scenarios A and B (T8)
ART.LOCAL       <- FALSE    # TRUE runs the simulations here on a gene subset
ART.LOCAL.GENES <- 300
ART.LOCAL.R     <- 50
T2.FWER         <- 0.01
T11.MAX.GENES   <- 1000
LEAK.MIN.MEAN   <- 0.02; LEAK.MIN.CELLS <- 50     # the gene filter of estimator_leak.R
MBAR.BREAKS     <- c(0, 0.1, 0.3, 1, 3, 10, Inf)
PHI.BREAKS      <- c(0, 0.05, 0.2, 0.5, 1, Inf)
set.seed(SEED.SECTION + 24)

stopifnot("artefact_tests.R needs the Section 1 checkpoint" = file.exists(ckpt_path(1)))
load(ckpt_path(1))
dir.create(CHECK.DIR, recursive = TRUE, showWarnings = FALSE); dir.create(CHECK.FIG.DIR, recursive = TRUE, showWarnings = FALSE)

## Groups, offsets and share corrections as in analysis.R Section 2 and estimator_leak.R; genes with data in all four groups.
D.REF    <- median(c(colSums(MIX.SC), colSums(MIX.SE)))
LIB      <- list(Pc = colSums(MIX.SC), Pe = colSums(MIX.SE), Hc = colSums(HYB.SC) + colSums(HYB.SE), He = colSums(HYB.SC) + colSums(HYB.SE))
ART.EXPOS <- list(Pc = LIB$Pc / D.REF, Pe = LIB$Pe / D.REF, Hc = LIB$Hc / D.REF / 2, He = LIB$He / D.REF / 2)
ART.ETA   <- c(Pc = 0, Pe = 0, Hc = R.C, He = R.E)
FULL      <- list(Pc = MIX.SC, Pe = MIX.SE, Hc = HYB.SC, He = HYB.SE)
ART.GENES <- rownames(MIX.SC)[Reduce(`&`, lapply(FULL, function(m) rowMeans(m) >= LEAK.MIN.MEAN & rowSums(m > 0) >= LEAK.MIN.CELLS))]
ART.MATS  <- lapply(FULL, function(m) m[ART.GENES, , drop = FALSE])
ART.TERC  <- lapply(ART.EXPOS, function(e) cut(rank(e, ties.method = "first"), 3, labels = FALSE))
ART.HALVES <- with_local_seed(SEED.ART + 1, local({
  h <- sample.int(ncol(HYB.SC), floor(ncol(HYB.SC) / 2))
  list(Pc = sample.int(ncol(MIX.SC), floor(ncol(MIX.SC) / 2)), Pe = sample.int(ncol(MIX.SE), floor(ncol(MIX.SE) / 2)), Hc = h, He = h)
}))
NCELLS <- vapply(ART.EXPOS, length, integer(1))
cat(sprintf("Artefact tests: %d genes; cells %s\n", length(ART.GENES), paste(names(NCELLS), NCELLS, collapse = ", ")))

##############################################################################
## Definitions audit and baseline audit (rule 7, T0 step 1)
writeLines(c("# T0 audit: definitions and the earlier NB baseline", "",
  "Definitions in estimator_leak.R (checked in the code):", "",
  sprintf("- Trimming (L3): `estimator_pair_gene(trim = 0.005)` leaves out the top `ceiling(0.005 n)` cells by scaled count y/O within each fit's cells, at least one. On these data that is %s cells (the handoff inferred 4, 6 and 14). Mean and dispersion are both refitted without them.",
          paste(sprintf("%d of %d %s", ceiling(ART.FRAC * NCELLS), NCELLS, names(NCELLS)), collapse = ", ")),
  "- Own slope: `leak_regression()` returns `slope_cr` and `slope_mom`, the OLS slopes of Delta log2 phi on Delta log2 mu for each estimator, on genes with phi > 0 in every group of the contrast. L2 applied no class gate.",
  "- Leak slope: the OLS slope of [Delta log2 phi (Cox-Reid) - Delta log2 phi (moment)] on Delta log2 mu. On the same genes it equals slope_cr - slope_mom exactly, because OLS is linear in the response.",
  "- a (L6): the coefficient of A_i = M_i (1 - p_i) in the weighted regression of z_i = (y_i - M_i)^2 - A_i on A_i and B_i (`linear_excess_gene()`), so Var(y) = (1 + a) M + phi M^2.", "",
  "The NB baseline that the console comparison used (L3 against synthetic data in this session):", "",
  "- It was not gene-matched. Genes had log-uniform mean rates from 0.05 to 20 and phi from 0.05 to 1.5, with 1000 cells per parent and 2000 hybrid cells, lognormal parental offsets (log SD 0.6 and 0.7) and hybrid depth factors (log SD 0.4), and gamma factors.",
  "- The comparison also binned by phi only, so genes of different mean count were compared. At low mean counts trimming removes most of the variance even under the NB, so that baseline understated the trimming effect expected for the sparse hybrid genes.",
  "- T0 replaces it with a null simulated for every gene and group from that gene's own mean, Cox-Reid phi and offsets."),
  file.path(CHECK.DIR, "T0_audit.md"))

##############################################################################
## Observed statistics of every gene and group (T0, T3, T7a observed arms) and the per-gene fits for the data signature (T8, T9, T12)
# T7a thins the whole matrices, so the libraries and the reference depth are recomputed from the thinned counts.
THIN <- with_local_seed(SEED.ART + 2, lapply(FULL, function(m) matrix(rbinom(length(m), m, ART.THIN), nrow(m), dimnames = dimnames(m))))
D.THIN <- median(c(colSums(THIN$Pc), colSums(THIN$Pe))); LIB.THIN.H <- colSums(THIN$Hc) + colSums(THIN$He)
EXPO.THIN <- list(Pc = colSums(THIN$Pc) / D.THIN, Pe = colSums(THIN$Pe) / D.THIN, Hc = LIB.THIN.H / D.THIN / 2, He = LIB.THIN.H / D.THIN / 2)
OBS <- list(); OBS.TRIM <- list()
for (k in ART.GROUPS) {
  st <- mclapply(ART.GENES, function(g) art_gene_stats(as.numeric(ART.MATS[[k]][g, ]), ART.EXPOS[[k]], ART.ETA[[k]], ART.FRAC, ART.TERC[[k]], as.numeric(THIN[[k]][g, ]), EXPO.THIN[[k]]), mc.cores = CORES)
  OBS[[k]] <- data.frame(gene = ART.GENES, group = k, do.call(rbind, st), row.names = NULL); OBS.TRIM[[k]] <- lapply(st, attr, "trim")
}
OBS.ALL <- do.call(rbind, OBS)
ART.CORE <- do.call(rbind, mclapply(ART.GENES, art_gene_core, mats = ART.MATS, expos = ART.EXPOS, eta = ART.ETA, halves = ART.HALVES, frac = ART.FRAC, mc.cores = CORES))
rownames(ART.CORE) <- ART.GENES
DATA.SIG <- art_signature(ART.CORE, C.BAR, MBAR.BREAKS, seed = SEED.ART + 3)
ART.FITS <- ART.CORE[, paste0(rep(ART.GROUPS, each = 2), c(".mu", ".phi"))]

##############################################################################
## Inputs of the simulation array (T0, T8, T12), and its results
ART.CFG <- list(frac = ART.FRAC, R = ART.R, R4 = ART.R4, R11 = ART.R11, thin_p = ART.THIN, cbar = C.BAR, mbar_breaks = MBAR.BREAKS)
ART.CHUNKS <- split(ART.GENES, (seq_along(ART.GENES) - 1) %% min(ART.N.CHUNKS, length(ART.GENES)))     # genes dealt round-robin into chunks
ART.SCEN.GRID <- rbind(data.frame(scenario = "A", par1 = NA, par2 = NA, rep = seq_len(ART.SCEN.REPS)), data.frame(scenario = "B", par1 = NA, par2 = NA, rep = seq_len(ART.SCEN.REPS)),
                       data.frame(scenario = "S1", par1 = c(0.05, 0.1, 0.2), par2 = NA, rep = 1), data.frame(scenario = "S2", par1 = c(0.8, 0.9, 1.1), par2 = NA, rep = 1),
                       data.frame(scenario = "S3", par1 = c(0.005, 0.005, 0.01, 0.01), par2 = c(5, 10, 5, 10), rep = 1), data.frame(scenario = "S4", par1 = c(0.05, 0.1, 0.2), par2 = NA, rep = 1),
                       data.frame(scenario = "S5", par1 = c(0.25, 0.5), par2 = NA, rep = 1))
ART.TASKS <- rbind(data.frame(kind = "T0", index = seq_along(ART.CHUNKS)), data.frame(kind = "SCEN", index = seq_len(nrow(ART.SCEN.GRID))))
ART.STREAMS.GENE <- setNames(check_streams(length(ART.GENES), SEED.ART + 4), ART.GENES)
ART.STREAMS.SCEN <- check_streams(nrow(ART.SCEN.GRID), SEED.ART + 5)
save(ART.MATS, ART.EXPOS, ART.ETA, ART.TERC, ART.FITS, ART.HALVES, ART.GENES, ART.CFG, ART.TASKS, ART.CHUNKS, ART.SCEN.GRID, ART.STREAMS.GENE, ART.STREAMS.SCEN,
     file = file.path(INPUT.DIR, "artefact_inputs.rda"))
ART.FILES <- list.files(OUTPUT.DIR, "^artefact_output_[0-9]+of[0-9]+\\.rda$", full.names = TRUE)
ART.STATUS <- "run"
if (length(ART.FILES)) {
  ART.RES <- list()
  for (f in ART.FILES) { e <- new.env(); load(f, envir = e); ART.RES <- c(ART.RES, e$ART.PART) }
  stopifnot("artefact outputs do not cover every task of the inputs" = setequal(as.integer(names(ART.RES)), seq_len(nrow(ART.TASKS))))
  T0.PARTS <- unlist(ART.RES[as.character(which(ART.TASKS$kind == "T0"))], recursive = FALSE)
  SCEN.RES <- ART.RES[as.character(which(ART.TASKS$kind == "SCEN"))]
} else if (ART.LOCAL) {
  ART.STATUS <- "reduced"
  LOCAL.GENES <- with_local_seed(SEED.ART + 6, sample(ART.GENES, min(ART.LOCAL.GENES, length(ART.GENES))))
  T0.PARTS <- mclapply(split(LOCAL.GENES, (seq_along(LOCAL.GENES) - 1) %% CORES), art_null_chunk, mats = ART.MATS, expos = ART.EXPOS, eta = ART.ETA, frac = ART.FRAC,
                       terc = ART.TERC, R = ART.LOCAL.R, R4 = 3, R11 = 2, thin_p = ART.THIN, streams = ART.STREAMS.GENE, mc.cores = CORES, mc.preschedule = FALSE)
  SCEN.RES <- lapply(seq_len(nrow(ART.SCEN.GRID)), art_scenario_task, grid = ART.SCEN.GRID, fits = ART.FITS, expos = ART.EXPOS, eta = ART.ETA, halves = ART.HALVES,
                     genes = LOCAL.GENES, frac = ART.FRAC, cbar = C.BAR, mbar_breaks = MBAR.BREAKS, streams = ART.STREAMS.SCEN, cores = CORES)
} else {
  T0.PARTS <- NULL; SCEN.RES <- NULL
  cat("No artefact_output files in OUTPUT.DIR: submit cluster/submit/artefact_null.sub with artefact_inputs.rda, or set ART.LOCAL = TRUE.\n")
}
if (!is.null(T0.PARTS)) {
  T0.STATS <- do.call(rbind, lapply(T0.PARTS, `[[`, "stats"))
  T0.BURDEN <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) Reduce(`+`, lapply(T0.PARTS, function(p) p$burden[[k]])))
  T0.OBS.BURDEN <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) Reduce(`+`, lapply(T0.PARTS, function(p) p$obs_burden[[k]])))
  T0.T4 <- do.call(rbind, lapply(T0.PARTS, `[[`, "t4")); T0.T11 <- do.call(rbind, lapply(T0.PARTS, `[[`, "t11"))
}
NOT.RUN <- "the simulation array has no output yet (artefact_null.sub) and ART.LOCAL is FALSE"

##############################################################################
## T0. Trimming against the gene-matched NB null
# For each gene the observed trimming statistics are placed in the gene's own simulated distribution. Under the NB, 5% of genes fall below the 5th
# percentile of log2(trimmed phi / phi) and 5% above the 95th percentile of the excess share carried by the trimmed cells.
if (is.null(T0.PARTS)) check_save("T0", NULL, "not run", headline = NOT.RUN, script = CHECK.SCRIPT) else {
  T0.STATS$mean_bin <- cut(T0.STATS$mbar, MBAR.BREAKS, right = FALSE); T0.STATS$phi_bin <- cut(T0.STATS$phi, PHI.BREAKS, right = FALSE)
  T0.STATS$trim_below_q05 <- T0.STATS$log2_trim < T0.STATS$null_trim_q05; T0.STATS$exf_above_q95 <- T0.STATS$excess_frac > T0.STATS$null_exf_q95
  T0 <- list()
  for (k in ART.GROUPS) for (b in levels(T0.STATS$mean_bin)) for (p in levels(T0.STATS$phi_bin)) {
    i <- T0.STATS$group == k & T0.STATS$mean_bin == b & T0.STATS$phi_bin == p & !is.na(T0.STATS$phi_bin)
    if (sum(i) < 10) next
    a <- T0.STATS$trim_below_q05[i]; x <- T0.STATS$exf_above_q95[i]
    T0[[length(T0) + 1]] <- data.frame(group = k, mean_count = b, phi = p, genes = sum(i), share_trim_below_null_q05 = mean(a, na.rm = TRUE), se_trim = binom_se(mean(a, na.rm = TRUE), sum(!is.na(a))),
                                       share_excess_above_null_q95 = mean(x, na.rm = TRUE), se_excess = binom_se(mean(x, na.rm = TRUE), sum(!is.na(x))),
                                       median_log2_trim = median(T0.STATS$log2_trim[i], na.rm = TRUE), median_null_log2_trim = median(T0.STATS$null_trim_q50[i], na.rm = TRUE),
                                       median_excess_frac = median(T0.STATS$excess_frac[i], na.rm = TRUE), median_null_excess_frac = median(T0.STATS$null_exf_q50[i], na.rm = TRUE), row.names = NULL)
  }
  T0 <- do.call(rbind, T0)
  T0.G <- aggregate(cbind(share_trim_below_null_q05 = trim_below_q05, share_excess_above_null_q95 = exf_above_q95) ~ group, T0.STATS, mean)
  T0.G$genes <- as.vector(table(T0.STATS$group)[T0.G$group])
  check_save("T0", list(bins = T0, groups = T0.G, genes = T0.STATS), ART.STATUS,
             headline = c(sprintf("%s: trimmed phi below the gene's null 5th percentile in %.3f, excess share above its 95th percentile in %.3f of %d genes (NB: 0.05 each)",
                                  T0.G$group, T0.G$share_trim_below_null_q05, T0.G$share_excess_above_null_q95, T0.G$genes)),
             mc_error = sprintf("%d replicates per gene; binomial SE per bin in the rds", if (ART.STATUS == "reduced") ART.LOCAL.R else ART.R), script = CHECK.SCRIPT)
  fig_pdf("extra/check_T0_excess_share.pdf", 10, 3)
  par(mfrow = c(1, 4), mar = c(4, 4, 2, 1))
  for (k in ART.GROUPS) { i <- T0.STATS$group == k
    plot(T0.STATS$null_exf_q50[i], T0.STATS$excess_frac[i], pch = 19, cex = 0.3, col = adjustcolor(COLOR.GREY[["dark"]], 0.4), xlim = c(0, 1), ylim = c(0, 1),
         xlab = "null median excess share", ylab = "observed excess share", main = k); abline(0, 1, col = COLOR.ACCENT) }
  dev.off()
}

##############################################################################
## T1. Per-cell outlier burden
# A cell's burden is the number of genes for which it is trimmed. Its expectation under the NB depends on its offset and comes from the T0 replicates:
# excess = observed - null mean, z = excess / null SD, and p = (1 + replicates with burden >= observed) / (1 + replicates).
T1.CELLS <- NULL
if (is.null(T0.PARTS)) check_save("T1", NULL, "not run", headline = NOT.RUN, script = CHECK.SCRIPT) else {
  RP.GENES <- if (requireNamespace("org.Sc.sgd.db", quietly = TRUE)) intersect(go_gene_set("GO:0022626", org.Sc.sgd.db::org.Sc.sgd.db), rownames(MIX.SC)) else character(0)
  DBL <- if (requireNamespace("scDblFinder", quietly = TRUE) && requireNamespace("SingleCellExperiment", quietly = TRUE)) "available" else "scDblFinder not installed"
  CK7 <- if (file.exists(ckpt_path(7))) { e7 <- new.env(); load(ckpt_path(7), envir = e7); e7 } else NULL
  DS7 <- c(Pc = "MIX.SC", Pe = "MIX.SE", Hc = "HYB.SC", He = "HYB.SE")
  T1 <- list(); T1.COR <- list(); T1.TOP <- list(); ZMAX <- list(); Z3.NULL <- list()
  for (k in ART.GROUPS) {
    m <- FULL[[k]]; obs <- T0.OBS.BURDEN[[k]]; nb <- T0.BURDEN[[k]]
    mu0 <- colMeans(nb); sd0 <- apply(nb, 2, sd); z <- ifelse(sd0 > 0, (obs - mu0) / sd0, NA_real_)     # a cell never trimmed in the null has no z
    znull <- sweep(sweep(nb, 2, mu0), 2, ifelse(sd0 > 0, sd0, NA_real_), "/")                                # the same z for every null replicate
    ZMAX[[k]] <- apply(znull, 1, max, na.rm = TRUE); Z3.NULL[[k]] <- mean(rowSums(znull > 3, na.rm = TRUE))
    p <- (1 + rowSums(t(nb) >= obs)) / (1 + nrow(nb))
    cells <- colnames(m); if (is.null(cells)) cells <- as.character(seq_len(ncol(m)))
    covs <- data.frame(library = LIB[[k]], detected = colSums(m > 0),
                       mito_share = if (exists("MITO.QC")) MITO.QC[[if (k %in% c("Hc", "He")) "HYB" else c(Pc = "MIX.SC", Pe = "MIX.SE")[[k]]]]$frac[cells] else NA_real_,
                       ribo_share = if (length(RP.GENES)) colSums(m[RP.GENES, , drop = FALSE]) / colSums(m) else NA_real_,
                       chat = if (k %in% c("Hc", "He") && length(C.HAT) == ncol(m)) C.HAT else NA_real_,
                       doublet_score = if (exists("MULTIPLET.SCORE")) MULTIPLET.SCORE[[if (k %in% c("Hc", "He")) "HYB" else c(Pc = "MIX.SC", Pe = "MIX.SE")[[k]]]]$score[cells] else NA_real_,
                       scdbl_score = if (DBL == "available") tryCatch(scDblFinder::scDblFinder(SingleCellExperiment::SingleCellExperiment(list(counts = m)))$scDblFinder.score, error = function(err) NA_real_) else NA_real_,
                       cc_axis = if (!is.null(CK7) && !is.null(CK7$CC.AXIS[[DS7[[k]]]])) CK7$CC.AXIS[[DS7[[k]]]][cells] else NA_real_,
                       met_state = if (!is.null(CK7) && !is.null(CK7$MET.STATE[[DS7[[k]]]])) as.integer(CK7$MET.STATE[[DS7[[k]]]][cells]) else NA_integer_)
    T1[[k]] <- data.frame(group = k, cell = cells, burden = obs, null_mean = mu0, null_sd = sd0, excess = obs - mu0, z = z, p = p, covs, row.names = NULL)
    for (cv in setdiff(names(covs), "met_state"))
      T1.COR[[paste(k, cv)]] <- data.frame(group = k, covariate = cv, t(art_spearman_ci(covs[[cv]], obs - mu0, seed = SEED.ART + 7)), row.names = NULL)
    T1.COR[[paste(k, "met_state")]] <- data.frame(group = k, covariate = "met_state (Kruskal-Wallis p)", cells = sum(is.finite(covs$met_state)),
                                                  rho = if (sum(is.finite(covs$met_state)) > 10) kruskal.test(obs - mu0, factor(covs$met_state))$p.value else NA, lo = NA, hi = NA)
    T1.TOP[[k]] <- head(T1[[k]][order(-T1[[k]]$z), ], 50)
  }
  T1.CELLS <- do.call(rbind, T1); T1.COR <- do.call(rbind, T1.COR); rownames(T1.COR) <- NULL
  T1.SUM <- list()
  for (k in ART.GROUPS) { x <- T1.CELLS[T1.CELLS$group == k, ]
    T1.SUM[[k]] <- data.frame(group = k, cells = nrow(x), expected_mean_burden = mean(x$null_mean), z_above_3 = sum(x$z > 3, na.rm = TRUE), expected_z_above_3 = round(Z3.NULL[[k]], 1), threshold_T2 = unname(quantile(ZMAX[[k]], 1 - T2.FWER)),
                              top1pct_share_of_excess = sum(sort(pmax(x$excess, 0), TRUE)[1:ceiling(0.01 * nrow(x))]) / sum(pmax(x$excess, 0))) }
  T1.SUM <- do.call(rbind, T1.SUM)
  x <- T1.COR[is.finite(T1.COR$rho) & !grepl("Kruskal", T1.COR$covariate), ]
  check_save("T1", list(cells = T1.CELLS, correlations = T1.COR, top50 = T1.TOP, summary = T1.SUM, blocked = c(scDblFinder = DBL, ribosomal = if (length(RP.GENES)) "GO:0022626" else "org.Sc.sgd.db not installed",
                                                                                                          cell_state = if (is.null(CK7)) "Section 7 checkpoint missing" else "used")), ART.STATUS,
             headline = c(sprintf("%s: cells with burden z > 3: %d (mean over NB replicates %.1f); top 1%% of cells hold %.3f of the positive excess", T1.SUM$group, T1.SUM$z_above_3, T1.SUM$expected_z_above_3,
                                  T1.SUM$top1pct_share_of_excess),
                          sprintf("%s excess burden vs %s: Spearman %.3f (95%% %.3f to %.3f)", x$group, x$covariate, x$rho, x$lo, x$hi),
                          sprintf("not available: %s", paste(c(if (DBL != "available") "scDblFinder", if (!length(RP.GENES)) "ribosomal share (org.Sc.sgd.db)", if (is.null(CK7)) "cell cycle and metabolic state (Section 7 checkpoint)",
                                                                if (!exists("MITO.QC")) "mitochondrial share and doublet score (rerun Section 1 to save MITO.QC and MULTIPLET.SCORE)"), collapse = ", "))),
             mc_error = "burden null from the T0 replicates; Spearman intervals from 500 cell resamples", script = CHECK.SCRIPT)
  fig_pdf("extra/check_T1_burden_library.pdf", 10, 3)
  par(mfrow = c(1, 4), mar = c(4, 4, 2, 1))
  for (k in ART.GROUPS) { x <- T1.CELLS[T1.CELLS$group == k, ]
    plot(log10(x$library), x$excess, pch = 19, cex = 0.4, col = adjustcolor(COLOR.GREY[["dark"]], 0.5), xlab = "log10 library", ylab = "excess burden", main = k); abline(h = 0, col = COLOR.ACCENT) }
  dev.off()
}

##############################################################################
## T2. Cell-level QC sensitivity: the fits on the cells left by the T2 rule (header), against all cells
if (is.null(T1.CELLS)) check_save("T2", NULL, "not run", headline = NOT.RUN, script = CHECK.SCRIPT) else {
  RM <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) which(T1.CELLS$z[T1.CELLS$group == k] > T1.SUM$threshold_T2[T1.SUM$group == k]))
  RM$Hc <- RM$He <- sort(union(RM$Hc, RM$He))
  KEEP.MATS <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) if (length(RM[[k]])) ART.MATS[[k]][, -RM[[k]], drop = FALSE] else ART.MATS[[k]])
  KEEP.EXPOS <- lapply(setNames(ART.GROUPS, ART.GROUPS), function(k) if (length(RM[[k]])) ART.EXPOS[[k]][-RM[[k]]] else ART.EXPOS[[k]])
  KEEP.HALVES <- with_local_seed(SEED.ART + 8, local({ h <- sample.int(ncol(KEEP.MATS$Hc), floor(ncol(KEEP.MATS$Hc) / 2))
    list(Pc = sample.int(ncol(KEEP.MATS$Pc), floor(ncol(KEEP.MATS$Pc) / 2)), Pe = sample.int(ncol(KEEP.MATS$Pe), floor(ncol(KEEP.MATS$Pe) / 2)), Hc = h, He = h) }))
  T2 <- list()
  for (arm in c("all cells", "retained cells")) {
    mats <- if (arm == "all cells") ART.MATS else KEEP.MATS; expos <- if (arm == "all cells") ART.EXPOS else KEEP.EXPOS; halves <- if (arm == "all cells") ART.HALVES else KEEP.HALVES
    names(mats) <- names(expos) <- names(halves) <- ART.GROUPS
    fits <- do.call(rbind, mclapply(ART.GENES, estimator_pair_gene, mats = mats, expos = expos, eta = ART.ETA, halves = halves, mc.cores = CORES)); rownames(fits) <- ART.GENES
    trimf <- do.call(rbind, mclapply(ART.GENES, estimator_pair_gene, mats = mats, expos = expos, eta = ART.ETA, halves = halves, trim = ART.FRAC, mc.cores = CORES)); rownames(trimf) <- ART.GENES
    for (k in ART.GROUPS) {
      cr <- fits[, paste0("same.", k, ".phi_cr")]; mo <- fits[, paste0("same.", k, ".phi_mom")]; bp <- cut(cr, PHI.BREAKS, right = FALSE); ok <- cr > 0 & mo > 0
      tr <- log2(trimf[, paste0("same.", k, ".phi_cr")] / cr)
      le <- t(vapply(ART.GENES, linear_excess_gene, numeric(4), mat = mats[[k]], expo = expos[[k]]))
      for (l in levels(bp)) { i <- ok & bp == l
        if (sum(i, na.rm = TRUE) < 10) next
        T2[[paste(arm, k, l)]] <- data.frame(cells = arm, group = k, phi_cr = l, genes = sum(i, na.rm = TRUE), removed = if (arm == "all cells") 0 else length(RM[[k]]),
                                             median_log2_cr_over_mom = median(log2(cr / mo)[i], na.rm = TRUE), median_log2_trim_cr = median(tr[i][is.finite(tr[i])]),
                                             median_a = median(le[i, "a"], na.rm = TRUE), share_a_above_2se = mean(le[i, "a"] / le[i, "a_se"] > 2, na.rm = TRUE), row.names = NULL) }
    }
    for (a2 in c("same", "disjoint")) { d <- leak_contrasts(fits, a2, C.BAR)
      for (md in c("cis", "trans")) { dd <- d[d$contrast == md, ]
        T2[[paste(arm, a2, md)]] <- data.frame(cells = arm, group = paste(a2, md), phi_cr = "slopes", genes = sum(is.finite(dd$dlog2mu) & is.finite(dd$diff)), removed = NA,
                                               leak_regression(dd$dlog2mu, dd$diff, dd$dphi_cr, dd$dphi_mom, 500, SEED.ART + 9)[c("slope", "boot_lo", "boot_hi", "slope_cr", "slope_mom")], row.names = NULL) } }
  }
  T2.BINS <- do.call(rbind, T2[!grepl("cis$|trans$", names(T2))]); T2.SLOPES <- do.call(rbind, T2[grepl("cis$|trans$", names(T2))]); rownames(T2.BINS) <- rownames(T2.SLOPES) <- NULL
  check_save("T2", list(removed = RM, bins = T2.BINS, slopes = T2.SLOPES, rule = sprintf("z above the %g quantile of the per-replicate maximum z (T1 threshold_T2)", 1 - T2.FWER)), ART.STATUS,
             headline = c(sprintf("removed cells: %s", paste(sprintf("%s %d", ART.GROUPS, vapply(RM, length, integer(1))), collapse = ", ")),
                          sprintf("%s, %s: leak slope %.4f (%.4f to %.4f); own slopes Cox-Reid %.3f, moment %.3f", T2.SLOPES$cells, T2.SLOPES$group, T2.SLOPES$slope, T2.SLOPES$boot_lo, T2.SLOPES$boot_hi,
                                  T2.SLOPES$slope_cr, T2.SLOPES$slope_mom),
                          sprintf("%s %s phi %s: log2(CR / moment) %.3f, trimming %.3f, a %.3f (a > 2 SE %.3f)", T2.BINS$cells, T2.BINS$group, T2.BINS$phi_cr, T2.BINS$median_log2_cr_over_mom,
                                  T2.BINS$median_log2_trim_cr, T2.BINS$median_a, T2.BINS$share_a_above_2se)),
             mc_error = "gene bootstrap (500) for the slopes", script = CHECK.SCRIPT)
}

##############################################################################
## T3. Dispersion by depth tercile, against the T0 null; genes with an expected count of at least 1 per cell in every tercile
T3 <- list()
for (k in ART.GROUPS) {
  x <- OBS[[k]]; ok1 <- x$mstrat1 >= 1 & x$mstrat2 >= 1 & x$mstrat3 >= 1; mb <- cut(x$mbar, MBAR.BREAKS, right = FALSE)
  nl <- if (!is.null(T0.PARTS)) T0.STATS[T0.STATS$group == k, ][match(x$gene, T0.STATS$gene[T0.STATS$group == k]), ] else NULL
  for (l in c("all", levels(mb))) { i <- ok1 & (l == "all" | mb == l) & is.finite(x$strat1)
    if (sum(i) < 10) next
    T3[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), genes_all_strata_ge1 = sum(ok1),
                                    obs_low = median(x$strat1[i], na.rm = TRUE), obs_mid = median(x$strat2[i], na.rm = TRUE), obs_high = median(x$strat3[i], na.rm = TRUE),
                                    null_low = if (is.null(nl)) NA else median(nl$null_strat1[i], na.rm = TRUE), null_mid = if (is.null(nl)) NA else median(nl$null_strat2[i], na.rm = TRUE),
                                    null_high = if (is.null(nl)) NA else median(nl$null_strat3[i], na.rm = TRUE), row.names = NULL) }
}
T3 <- do.call(rbind, T3)
x <- T3[T3$mean_count == "all", ]
check_save("T3", T3, if (is.null(T0.PARTS)) "run (no null)" else ART.STATUS,
           headline = sprintf("%s (%d genes with >= 1 count per cell in every tercile): median log2(phi tercile / phi all) low %.3f, mid %.3f, high depth %.3f; NB null %.3f, %.3f, %.3f",
                              x$group, x$genes, x$obs_low, x$obs_mid, x$obs_high, x$null_low, x$null_mid, x$null_high),
           mc_error = "null medians from the T0 replicates", script = CHECK.SCRIPT)

##############################################################################
## T4. Non-proportional scaling fits: E[y] = O mu, mu O^beta and mu (O + rho), each NB with phi free (art_nb_models())
T4.FITS <- list()
for (k in ART.GROUPS) T4.FITS[[k]] <- data.frame(group = k, gene = ART.GENES, mbar = OBS[[k]]$mbar, mu = OBS[[k]]$mu,
                                                 do.call(rbind, mclapply(ART.GENES, function(g) art_nb_models(as.numeric(ART.MATS[[k]][g, ]), ART.EXPOS[[k]]), mc.cores = CORES)), row.names = NULL)
T4.FITS <- do.call(rbind, T4.FITS)
T4 <- list()
for (src in c("data", "NB null")) {
  d <- if (src == "data") T4.FITS else if (!is.null(T0.PARTS)) transform(T0.T4, group = ART.GROUPS[group]) else NULL
  if (is.null(d)) next
  d$mean_bin <- cut(d$mbar, MBAR.BREAKS, right = FALSE)
  for (k in ART.GROUPS) for (l in c("all", levels(d$mean_bin))) { i <- d$group == k & (l == "all" | d$mean_bin == l) & is.finite(d$beta)
    if (sum(i) < 10) next
    T4[[paste(src, k, l)]] <- data.frame(source = src, group = k, mean_count = l, fits = sum(i), median_beta = median(d$beta[i]), share_beta_below = mean((d$beta[i] - 1) / d$beta_se[i] < -2, na.rm = TRUE),
                                         share_beta_above = mean((d$beta[i] - 1) / d$beta_se[i] > 2, na.rm = TRUE), median_rho = median(d$rho[i]), share_rho_2se = mean(d$rho[i] / d$rho_se[i] > 2, na.rm = TRUE),
                                         median_lr_beta = median(2 * (d$ll_b[i] - d$ll_a[i])), share_lr_beta_3.84 = mean(2 * (d$ll_b[i] - d$ll_a[i]) > 3.84),
                                         median_lr_rho = median(2 * (d$ll_c[i] - d$ll_a[i])), share_lr_rho_3.84 = mean(2 * (d$ll_c[i] - d$ll_a[i]) > 3.84),
                                         median_log2_phi_b_over_a = median(log2(d$phi_b[i] / d$phi_a[i])), median_log2_phi_c_over_a = median(log2(d$phi_c[i] / d$phi_a[i])), row.names = NULL) }
}
T4 <- do.call(rbind, T4)
T4.TREND <- list()
for (k in ART.GROUPS) { d <- T4.FITS[T4.FITS$group == k & T4.FITS$mu > 0 & is.finite(T4.FITS$phi_a), ]
  T4.TREND[[k]] <- data.frame(group = k, trend_a = unname(coef(lm(log2(d$phi_a) ~ log2(d$mu)))[2]), trend_b = unname(coef(lm(log2(d$phi_b) ~ log2(d$mu)))[2]),
                              trend_c = unname(coef(lm(log2(d$phi_c) ~ log2(d$mu)))[2])) }
T4.TREND <- do.call(rbind, T4.TREND)
x <- T4[T4$mean_count == "all", ]
check_save("T4", list(summary = T4, trends = T4.TREND, fits = T4.FITS), if (is.null(T0.PARTS)) "run (no null)" else ART.STATUS,
           headline = c(sprintf("%s %s: beta median %.3f (below 1 by 2 SE %.3f, above %.3f), rho = lambda/mu median %.3g (> 2 SE %.3f); LR > 3.84 for beta %.3f, for rho %.3f; log2 phi change %.3f (beta), %.3f (rho)",
                                x$source, x$group, x$median_beta, x$share_beta_below, x$share_beta_above, x$median_rho, x$share_rho_2se, x$share_lr_beta_3.84, x$share_lr_rho_3.84,
                                x$median_log2_phi_b_over_a, x$median_log2_phi_c_over_a),
                        sprintf("%s: trend of log2 phi on log2 mu, model a %.3f, b %.3f, c %.3f", T4.TREND$group, T4.TREND$trend_a, T4.TREND$trend_b, T4.TREND$trend_c)),
           mc_error = "Hessian SEs per gene; null from the first T0 replicates", script = CHECK.SCRIPT)

##############################################################################
## T5. Ambient RNA profile: the inputs are per-cell count tables (data/single_cell/*_for_genes_percell_pergene_count.csv) without empty droplets
check_save("T5", NULL, "blocked", headline = "blocked, raw matrix missing: the project reads filtered per-cell count tables and has no empty-droplet barcodes", script = CHECK.SCRIPT)

##############################################################################
## T6. Size factor alternatives (sensitivity only): (i) library over its median, (ii) scran deconvolution, (iii) library without the 50 most expressed genes
TOP50.G <- lapply(FULL, function(m) order(rowSums(m), decreasing = TRUE)[1:50]); TOP50.G$Hc <- TOP50.G$He <- order(rowSums(HYB.SC) + rowSums(HYB.SE), decreasing = TRUE)[1:50]
LIB.X <- list(Pc = colSums(MIX.SC[-TOP50.G$Pc, ]), Pe = colSums(MIX.SE[-TOP50.G$Pe, ]), Hc = colSums(HYB.SC[-TOP50.G$Hc, ]) + colSums(HYB.SE[-TOP50.G$Hc, ]))
D.X <- median(c(LIB.X$Pc, LIB.X$Pe))
SF <- list(library = ART.EXPOS, no_top50 = list(Pc = LIB.X$Pc / D.X, Pe = LIB.X$Pe / D.X, Hc = LIB.X$Hc / D.X / 2, He = LIB.X$Hc / D.X / 2))
if (requireNamespace("scran", quietly = TRUE)) {
  sfp <- with_local_seed(SEED.ART + 10, scran::calculateSumFactors(cbind(MIX.SC, MIX.SE)))
  sfh <- with_local_seed(SEED.ART + 11, scran::calculateSumFactors(HYB.SC + HYB.SE))
  SF$scran <- list(Pc = sfp[seq_len(ncol(MIX.SC))], Pe = sfp[ncol(MIX.SC) + seq_len(ncol(MIX.SE))], Hc = sfh / 2 * median(LIB$Hc) / D.REF / median(sfh), He = sfh / 2 * median(LIB$Hc) / D.REF / median(sfh))
}
T6 <- list(); T6.COR <- list()
for (sf in names(SF)) for (k in ART.GROUPS) {
  ph <- vapply(ART.GENES, function(g) fit_dispersion(as.numeric(ART.MATS[[k]][g, ]), SF[[sf]][[k]], ART.ETA[[k]])[["disp"]], numeric(1))
  le <- t(vapply(ART.GENES, linear_excess_gene, numeric(4), mat = ART.MATS[[k]], expo = SF[[sf]][[k]]))
  ok <- ph > 0 & OBS[[k]]$phi > 0; mb <- cut(OBS[[k]]$mbar, MBAR.BREAKS, right = FALSE)
  T6[[paste(sf, k)]] <- data.frame(size_factor = sf, group = k, median_log2_phi_vs_library = median(log2(ph[ok] / OBS[[k]]$phi[ok])), median_a = median(le[, "a"], na.rm = TRUE),
                                   share_a_above_2se = mean(le[, "a"] / le[, "a_se"] > 2, na.rm = TRUE), trend_log2phi_on_log2mu = unname(coef(lm(log2(ph[ph > 0]) ~ log2(OBS[[k]]$mu[ph > 0])))[2]), row.names = NULL)
  T6.COR[[paste(sf, k)]] <- data.frame(size_factor = sf, group = k, spearman_with_library = cor(SF[[sf]][[k]], ART.EXPOS[[k]], method = "spearman"))
}
T6 <- do.call(rbind, T6); T6.COR <- do.call(rbind, T6.COR)
check_save("T6", list(summary = T6, size_factor_correlation = T6.COR), "run",
           headline = c(sprintf("%s %s: log2(phi / phi with library) %.3f, a median %.3f (a > 2 SE %.3f), trend %.3f", T6$size_factor, T6$group, T6$median_log2_phi_vs_library, T6$median_a,
                                T6$share_a_above_2se, T6$trend_log2phi_on_log2mu),
                        sprintf("%s %s: Spearman with library size factors %.3f", T6.COR$size_factor, T6.COR$group, T6.COR$spearman_with_library),
                        if (is.null(SF$scran)) "scran not installed: deconvolution size factors not computed" else "scran deconvolution included"),
           mc_error = "exact", script = CHECK.SCRIPT)

##############################################################################
## T7a. Constant-p thinning (p = 0.5), against the T0 null thinned the same way (the null keeps the offsets; the data recompute them)
T7A <- list()
for (k in ART.GROUPS) {
  x <- OBS[[k]]; mb <- cut(x$mbar, MBAR.BREAKS, right = FALSE)
  nl <- if (!is.null(T0.PARTS)) T0.STATS[T0.STATS$group == k, ][match(x$gene, T0.STATS$gene[T0.STATS$group == k]), ] else NULL
  for (l in levels(mb)) { i <- mb == l
    if (sum(i) < 10) next
    T7A[[paste(k, l)]] <- data.frame(group = k, mean_count = l, genes = sum(i), median_log2_phi_thin = median(x$thin_log2phi[i][is.finite(x$thin_log2phi[i])]), median_mu_ratio = median(x$thin_mu_ratio[i], na.rm = TRUE),
                                     null_log2_phi_thin = if (is.null(nl)) NA else median(nl$null_thin_log2phi[i], na.rm = TRUE), null_mu_ratio = if (is.null(nl)) NA else median(nl$null_thin_mu_ratio[i], na.rm = TRUE), row.names = NULL) }
}
T7A <- do.call(rbind, T7A)
check_save("T7a", T7A, if (is.null(T0.PARTS)) "run (no null)" else ART.STATUS,
           headline = sprintf("%s mean count %s: log2(phi thinned / phi) %.3f (NB null %.3f), mu ratio %.3f (null %.3f)", T7A$group, T7A$mean_count, T7A$median_log2_phi_thin, T7A$null_log2_phi_thin,
                              T7A$median_mu_ratio, T7A$null_mu_ratio), mc_error = "one thinning draw of the data", script = CHECK.SCRIPT)

##############################################################################
## T7b. Depth-matched trans: parental UMIs thinned to the hybrid allele depth, then the own slopes (gated as in T8)
F.PC <- min(1, median(colSums(HYB.SC)) / median(LIB$Pc)); F.PE <- min(1, median(colSums(HYB.SE)) / median(LIB$Pe))     # a parent shallower than the hybrid allele is kept whole
TH.MATS <- ART.MATS
TH.MATS$Pc <- with_local_seed(SEED.ART + 12, matrix(rbinom(length(ART.MATS$Pc), ART.MATS$Pc, F.PC), nrow(ART.MATS$Pc), dimnames = dimnames(ART.MATS$Pc)))
TH.MATS$Pe <- with_local_seed(SEED.ART + 13, matrix(rbinom(length(ART.MATS$Pe), ART.MATS$Pe, F.PE), nrow(ART.MATS$Pe), dimnames = dimnames(ART.MATS$Pe)))
LIB.TH <- list(Pc = with_local_seed(SEED.ART + 14, colSums(matrix(rbinom(length(MIX.SC), MIX.SC, F.PC), nrow(MIX.SC)))), Pe = with_local_seed(SEED.ART + 15, colSums(matrix(rbinom(length(MIX.SE), MIX.SE, F.PE), nrow(MIX.SE)))))
D.TH <- median(c(LIB.TH$Pc, LIB.TH$Pe))
TH.EXPOS <- ART.EXPOS; TH.EXPOS$Pc <- LIB.TH$Pc / D.TH; TH.EXPOS$Pe <- LIB.TH$Pe / D.TH
TH.CORE <- do.call(rbind, mclapply(ART.GENES, art_gene_core, mats = TH.MATS, expos = TH.EXPOS, eta = ART.ETA, halves = ART.HALVES, frac = ART.FRAC, mc.cores = CORES)); rownames(TH.CORE) <- ART.GENES
TH.SIG <- art_signature(TH.CORE, C.BAR, MBAR.BREAKS, seed = SEED.ART + 16)
T7B <- rbind(cbind(data = "observed depth", DATA.SIG$slopes), cbind(data = "parents thinned to hybrid depth", TH.SIG$slopes))
check_save("T7b", T7B, "run",
           headline = c(sprintf("UMI keep fractions: Pc %.3f, Pe %.3f (UMIs; reads are not in the project)", F.PC, F.PE),
                        sprintf("%s, %s %s: own slope %.3f (95%% %.3f to %.3f), %d genes", T7B$data, T7B$arm, T7B$contrast, T7B$slope, T7B$lo, T7B$hi, T7B$genes)),
           mc_error = "gene bootstrap (200) per slope; one thinning draw", script = CHECK.SCRIPT)

##############################################################################
## T8 and T12. Null-divergence and artefact scenarios against the data signature
if (is.null(SCEN.RES)) { check_save("T8", NULL, "not run", headline = NOT.RUN, script = CHECK.SCRIPT); check_save("T12", NULL, "not run", headline = NOT.RUN, script = CHECK.SCRIPT) } else {
  SC.SLOPES <- do.call(rbind, lapply(SCEN.RES, `[[`, "slopes")); SC.GROUPS <- do.call(rbind, lapply(SCEN.RES, `[[`, "groups")); SC.BINS <- do.call(rbind, lapply(SCEN.RES, `[[`, "bins"))
  T8.AB <- aggregate(cbind(slope, lo, hi, genes) ~ scenario + arm + contrast, SC.SLOPES[SC.SLOPES$scenario %in% c("A", "B"), ], mean)
  T8.AB$source <- paste("scenario", T8.AB$scenario, "(mean of", ART.SCEN.REPS, "datasets)")
  T8 <- rbind(cbind(source = "data", DATA.SIG$slopes[, c("arm", "contrast", "genes", "slope", "lo", "hi")]), T8.AB[, c("source", "arm", "contrast", "genes", "slope", "lo", "hi")])
  T8.MC <- aggregate(slope ~ scenario + arm + contrast, SC.SLOPES[SC.SLOPES$scenario %in% c("A", "B"), ], function(v) sd(v) / sqrt(length(v)))
  check_save("T8", list(slopes = T8, mc_se = T8.MC, replicates = SC.SLOPES[SC.SLOPES$scenario %in% c("A", "B"), ]), ART.STATUS,
             headline = sprintf("%s, %s %s: own slope %.3f (%.3f to %.3f), %d genes", T8$source, T8$arm, T8$contrast, T8$slope, T8$lo, T8$hi, T8$genes),
             mc_error = if (any(is.finite(T8.MC$slope))) sprintf("MC SE of the scenario slopes across datasets up to %.3f", max(T8.MC$slope, na.rm = TRUE)) else "one dataset per scenario: no MC SE across datasets", script = CHECK.SCRIPT)
  T12.G <- rbind(cbind(scenario = "data", par1 = NA, par2 = NA, DATA.SIG$groups), SC.GROUPS[, c("scenario", "par1", "par2", names(DATA.SIG$groups))])
  T12.S <- rbind(cbind(scenario = "data", par1 = NA, par2 = NA, DATA.SIG$slopes), SC.SLOPES[, c("scenario", "par1", "par2", names(DATA.SIG$slopes))])
  T12.B <- rbind(cbind(scenario = "data", par1 = NA, par2 = NA, DATA.SIG$bins), SC.BINS[, c("scenario", "par1", "par2", names(DATA.SIG$bins))])
  check_save("T12", list(groups = T12.G, slopes = T12.S, bins = T12.B), ART.STATUS,
             headline = c(sprintf("%s %s %s %s: trend of log2 phi on log2 mu %.3f, Spearman(a, mean) %.3f", T12.G$scenario, T12.G$par1, T12.G$par2, T12.G$group, T12.G$trend_log2phi_on_log2mu, T12.G$spearman_a_mean),
                          sprintf("%s %s %s, %s %s: own slope %.3f", T12.S$scenario, T12.S$par1, T12.S$par2, T12.S$arm, T12.S$contrast, T12.S$slope)),
             mc_error = "one dataset per artefact scenario; see T8 for the replicate spread of scenario B (S0)", script = CHECK.SCRIPT)
}

##############################################################################
## T9. Coupling by information: own slopes by mean-count tercile and on class-1 genes with mean count >= 3 in every group; frequency-driven share
# Delta log2 phi = -Delta log2 mu for a change in burst frequency at fixed burst size and 0 for a change in burst size at fixed frequency, so -slope reads
# roughly as the frequency-driven share. This is a rough reading: it assumes the two kinds of change are uncorrelated with each other across genes.
MB.GENE <- rowMeans(ART.CORE[, paste0(ART.GROUPS, ".mbar")]); TERC.G <- cut(rank(MB.GENE, ties.method = "first"), 3, labels = c("low", "mid", "high"))
HI <- Reduce(`&`, lapply(ART.GROUPS, function(k) ART.CORE[, paste0(k, ".mbar")] >= 3 & ART.CORE[, paste0(k, ".class")] %in% 1))
T9 <- rbind(cbind(genes_in = "all", DATA.SIG$slopes),
            do.call(rbind, lapply(levels(TERC.G), function(t) cbind(genes_in = paste("mean-count tercile", t), art_signature(ART.CORE[TERC.G == t, ], C.BAR, MBAR.BREAKS, seed = SEED.ART + 17)$slopes))),
            cbind(genes_in = "class 1, mean count >= 3 in every group", art_signature(ART.CORE[HI, , drop = FALSE], C.BAR, MBAR.BREAKS, seed = SEED.ART + 18)$slopes),
            cbind(genes_in = "parents thinned to hybrid depth (T7b)", TH.SIG$slopes))
T9$frequency_share_rough <- pmin(pmax(-T9$slope, 0), 1)
check_save("T9", T9, "run", headline = sprintf("%s, %s %s: own slope %.3f (%.3f to %.3f), %d genes; rough frequency-driven share %.2f", T9$genes_in, T9$arm, T9$contrast, T9$slope, T9$lo, T9$hi, T9$genes,
                                               T9$frequency_share_rough), mc_error = "gene bootstrap (200) per slope", script = CHECK.SCRIPT)

##############################################################################
## T10. Biology of the high-count cells: coherent gene sets against mean-matched random sets, and cell state of the high-burden cells
if (!requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
  check_save("T10", NULL, "blocked", headline = "org.Sc.sgd.db not installed: the GO gene sets cannot be built", script = CHECK.SCRIPT)
} else {
  SETS <- list(ribosome = "GO:0022626", ribi = "GO:0042254", nucleosome = "GO:0000786", heat_response = "GO:0034605")
  SETS <- lapply(SETS, function(go) intersect(go_gene_set(go, org.Sc.sgd.db::org.Sc.sgd.db), ART.GENES))
  T10 <- list()
  for (k in ART.GROUPS) {
    trims <- setNames(OBS.TRIM[[k]], ART.GENES); mb <- cut(OBS[[k]]$mbar, MBAR.BREAKS, right = FALSE); names(mb) <- ART.GENES
    for (s in names(SETS)) { gs <- SETS[[s]]
      if (length(gs) < 5) { T10[[paste(k, s)]] <- data.frame(group = k, set = s, genes = length(gs), jaccard = NA, null_mean = NA, p = NA); next }
      obs <- art_jaccard_mean(trims[gs])
      nul <- with_local_seed(SEED.ART + 19, vapply(seq_len(1000), function(b) art_jaccard_mean(trims[unlist(lapply(split(gs, mb[gs]), function(v) sample(ART.GENES[mb == mb[v[1]]], length(v))))]), numeric(1)))
      T10[[paste(k, s)]] <- data.frame(group = k, set = s, genes = length(gs), jaccard = obs, null_mean = mean(nul), p = (1 + sum(nul >= obs)) / 1001, row.names = NULL) }
  }
  T10 <- do.call(rbind, T10)
  T10.STATE <- if (!is.null(T1.CELLS) && any(is.finite(T1.CELLS$cc_axis))) aggregate(cbind(cc_axis) ~ group + (z > 3), T1.CELLS, median) else "cell state not available (Section 7 checkpoint)"
  check_save("T10", list(sets = T10, state = T10.STATE), "run",
             headline = c(sprintf("%s %s (%d genes): mean Jaccard of trimmed cells %.4f against %.4f for mean-matched random sets, p %.3f", T10$group, T10$set, T10$genes, T10$jaccard, T10$null_mean, T10$p),
                          if (is.character(T10.STATE)) T10.STATE else "cell-cycle axis of high-burden cells in the rds"),
             mc_error = "1000 mean-matched random sets", script = CHECK.SCRIPT)
}

##############################################################################
## T11. Estimator influence: leave-one-out influence of the top cells, data against the T0 null
T11.G <- with_local_seed(SEED.ART + 20, sample(ART.GENES, min(T11.MAX.GENES, length(ART.GENES))))
T11.D <- list()
for (k in ART.GROUPS) T11.D[[k]] <- data.frame(group = k, gene = T11.G, mbar = OBS[[k]]$mbar[match(T11.G, ART.GENES)],
                                               do.call(rbind, mclapply(T11.G, function(g) art_influence(as.numeric(ART.MATS[[k]][g, ]), ART.EXPOS[[k]], ART.ETA[[k]]), mc.cores = CORES)), row.names = NULL)
T11.D <- do.call(rbind, T11.D)
T11 <- list()
for (src in c("data", "NB null")) {
  d <- if (src == "data") T11.D else if (!is.null(T0.PARTS)) transform(T0.T11, group = ART.GROUPS[group]) else NULL
  if (is.null(d)) next
  d$mean_bin <- cut(d$mbar, MBAR.BREAKS, right = FALSE)
  for (k in ART.GROUPS) for (l in levels(d$mean_bin)) { i <- d$group == k & d$mean_bin == l
    if (sum(i) < 10) next
    T11[[paste(src, k, l)]] <- data.frame(source = src, group = k, mean_count = l, fits = sum(i), max_share_cr = median(d$max_share_cr[i], na.rm = TRUE), top5_share_cr = median(d$top5_share_cr[i], na.rm = TRUE),
                                          max_share_mom = median(d$max_share_mom[i], na.rm = TRUE), top5_share_mom = median(d$top5_share_mom[i], na.rm = TRUE), row.names = NULL) }
}
T11 <- do.call(rbind, T11)
check_save("T11", list(summary = T11, genes = T11.D), if (is.null(T0.PARTS)) "run (no null)" else ART.STATUS,
           headline = sprintf("%s %s mean count %s: top-5 cells carry %.3f of phi (Cox-Reid) and %.3f (moment); largest single cell %.3f and %.3f", T11$source, T11$group, T11$mean_count,
                              T11$top5_share_cr, T11$top5_share_mom, T11$max_share_cr, T11$max_share_mom),
           mc_error = "medians over genes; null from the first T0 replicates; the L1 table before and after T2 is in T2", script = CHECK.SCRIPT)

##############################################################################
## REPORT_artefacts.md: one table per hypothesis with the tests that bear on it and their headline numbers; interpretation is left to the reader.
SUMMARY <- read.csv(file.path(CHECK.DIR, "summary.csv"), stringsAsFactors = FALSE)
HYP <- list(H1 = c("T0", "T1", "T2", "T11"), H2 = c("T3", "T4", "T5", "T12"), H3 = c("T3", "T4", "T6", "T7a", "T12"), H4 = c("T1", "T10", "T12"), H5 = c("T7b", "T8", "T9"), H6 = c("T2", "T11"))
HYP.TEXT <- c(H1 = "A few low-depth or low-quality cells carry counts beyond the NB", H2 = "An additive component outside the offset", H3 = "Counts scale non-proportionally with the offset",
              H4 = "The high-count cells are real expression jackpots", H5 = "The negative Delta phi / Delta mu coupling is biology", H6 = "The estimators differ because of outlier cells")
REP <- c("# Tests of why the data depart from the NB", "", sprintf("Generated %s by analysis/checks/artefact_tests.R. Numbers only; each test's full result is in results/checks/<ID>.rds.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
for (h in names(HYP)) {
  x <- SUMMARY[SUMMARY$ID %in% HYP[[h]], ]
  REP <- c(REP, sprintf("## %s. %s", h, HYP.TEXT[[h]]), "", "| Test | Status | Headline numbers |", "|---|---|---|",
           sprintf("| %s | %s | %s |", x$ID, x$status, gsub("\\|", "/", gsub("; ", "<br>", x$headline))), "")
}
REP <- c(REP, "## Blocked or limited", "", "- T5: the raw droplet matrix is not in the project.", "- T7b: reads are not in the project, so parental UMIs were thinned instead.",
         "- T1 and T10: scDblFinder, org.Sc.sgd.db and the Section 7 checkpoint are used when present; the T1 row lists what was missing.",
         "- Definitions: trimming leaves out ceiling(0.005 n) cells (5, 7 and 15), not the 4, 6 and 14 inferred in the handoff; see T0_audit.md.")
writeLines(REP, file.path(CHECK.DIR, "REPORT_artefacts.md"))
cat("Wrote", file.path(CHECK.DIR, "REPORT_artefacts.md"), "\n")
