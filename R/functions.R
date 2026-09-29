###############################################################
### Functions.R
### Function library for the cis/trans mean and noise analysis,
### sourced by Analysis.R and by the SLURM cluster
### scripts (gene_pilot.R, gene_split.R, gene_boot.R,
### gene_perm.R, coexpr_boot.R, coexpr_perm.R, nupop_occupancy.R,
### cluster_stability.R, go_enrich.R)
###
### Outline (function name - purpose), grouped by section:
###
###   1. OFFSET NEGATIVE-BINOMIAL FIT
###     neg_binom_fit_offset() - Offset NB fit for one gene: exposure-weighted mean plus glm.nb dispersion, with a Poisson-vs-NB pre-check.
###     fit_counts_offset() - Matrix version of neg_binom_fit_offset(): fits every gene (row) against its exposure vector.
###     boot_disp_logse() - Bootstrap SE of log(size) over cells, for one gene.
###     .fit_one() - Offset NB fit for one gene, used throughout the bootstrap and permutation machinery.
###     .fit_split() - Splits a pooled sample by a pre-drawn permutation and fits each half with .fit_one().
###   1b. INTERNAL-PILOT SPLIT FRACTION (f*)
###     pilot_split_se_one() - One gene's bootstrap SE of log(mu) and log(size), for the four split-independent pilot datasets.
###     pilot_split_se() - Serial wrapper running pilot_split_se_one() across a gene set, for a quick local check.
###     estimate_f_star() - Pooled optimal split fraction (f*) for the mean and noise axes, from the pilot SE curve.
###     split_indices_by_depth() - Depth-matched split of depth-ordered hybrid cells at an arbitrary fraction.
###     fit_counts_offset_row() - Single-gene row version of fit_counts_offset(), for parLapply distribution across genes.
###     fit_counts_offset_parallel() - Runs fit_counts_offset_row() across all genes on an already-open cluster.
###     raw_gene_prefilter() - Fit-free preliminary gene filter using only raw count summaries, ahead of the real NB fit.
###     .fstar_from_r() - Closed-form optimal split fraction f* from the parent/hybrid noise ratio r.
###   2. GENE FILTER (marginal information only; never on a contrast)
###     gene_pass_group() - Per-group pass/fail test (mean count and expression-fraction thresholds), used by build_gene_sets().
###     build_gene_sets() - Builds the final per-mode and full gene sets from per-dataset fits.
###     refine_by_boundary() - Bootstrap-checks candidate genes near the Poisson/NB boundary and reclassifies borderline DISP = Inf calls.
###     chk() - Calibration helper: flags a fit as usable based on convergence and finite, non-degenerate dispersion.
###     chk_prec() - Precision companion to chk(): flags fits whose SE is too imprecise to trust.
###   3. PAIRED PER-GENE BOOTSTRAP OF CONTRASTS
###     make_draws() - Builds the resampling index draws used by the paired per-gene bootstrap.
###     boot_contrasts_one() - One gene's bootstrap contrasts across all modes; the unit of work parLapply distributes in gene_boot.R.
###     boot_contrasts() - Serial, local wrapper running boot_contrasts_one() across all genes.
###     add_burst_contrasts() - Derives burst-frequency, burst-size, and kinetic-balance contrasts from the mean and size contrasts already in a data frame.
###     .contrast_value() - Looks up the right fitted quantity (mu or size, from the right split) for one contrast mode.
###   4. DISATTENUATED MEAN-DISP COUPLING (population summary)
###     eiv_components() - Errors-in-variables mean-size coupling for one mode: attenuation-corrected correlation and slope.
###     eiv_boot_ci() - Bootstrap confidence interval for one of eiv_components()'s statistics.
###     eiv_table() - Table of disattenuated mean-size coupling across a set of modes.
###   5. PLOTS (gene-level, square symmetric panels, SE bars)
###     .sym() - Symmetric axis limits spanning a vector of values plus their SE.
###     plot_cis_trans() - Cis vs trans scatter for one quantity (mean/size/bfreq/bsize), coloured by regulatory class, with SE bars.
###     plot_mean_bfreq() - Mean vs noise (size) scatter for one mode, coloured by regulatory/dominance class.
###     plot_burst_kinetics() - Burst kinetics scatter for one mode: net mean change against the stored kinetic-balance contrast.
###     sig_hist() - Significance-shaded histogram of a contrast, coloured by direction and permutation significance.
###   6. PERMUTATION NULL  (per-gene significance for the contrasts)
###     make_perms() - Mirrors the bootstrap but shuffles labels to null each contrast.
###     permute_contrasts_one() - One gene's permutation-null contrasts across all modes, including kinetic balance; the unit of work parLapply distributes in gene_perm.R.
###     permute_contrasts() - Serial, local wrapper running permute_contrasts_one() across all genes.
###   7. CO-EXPRESSION  (residual co-fluctuation, cis/trans decomposed)
###     nb_residuals() - Pearson residuals from the offset NB fit, genes x cells.
###     shrink_cor() - Analytic shrinkage of a correlation matrix toward the identity matrix.
###     coexpr_decompose() - Decomposes the point-estimate correlation structure into total/cis/trans/dpar_sc/dpar_se matrices.
###     coexpr_raw_cor() - Raw (unshrunk) per-dataset correlation matrices, kept for diagnostic comparison against the shrunk versions.
###     coexpr_axis_cis_trans() - Exact cis/trans decomposition of one eigenvector's eigenvalue, via total = cis + trans.
###     make_coexpr_draws() - Builds the cell-resampling draws used by the co-expression bootstrap.
###     coexpr_bootstrap_one() - One bootstrap draw of the co-expression decomposition, for one resampled cell set.
###     make_coexpr_perm_draws() - Builds the draws used by the co-expression permutation null.
###     coexpr_perm_one() - One permutation draw's top-N squared-eigenvalue spectrum for total/cis/trans/dpar_sc/dpar_se.
###     assemble_coexpr_bootstrap() - Builds CB (per-pair estimate/SE/z/p) from the point estimate and the streaming bootstrap accumulator.
###     coexpr_acc_init() - Initializes the streaming sum/sum-of-squares accumulator used by the co-expression bootstrap.
###     coexpr_acc_update() - Folds one chunk of coexpr_bootstrap_one() draws into the running accumulator.
###     coexpr_acc_finalize() - Converts the finished accumulator into CB's per-pair estimate/SE structure.
###     coexpr_bootstrap() - Serial, local wrapper running the co-expression bootstrap without the cluster round trip.
###     gene_reliability() - Per-gene reliability score: the fraction of a gene's variance that is real signal versus sampling noise.
###     check_coexpr_reliability() - Tests whether observed co-expression bootstrap SE tracks the NB-predicted attenuation factor.
###     row_cor() - Row-wise Pearson correlation between two same-shape matrices.
###     stratified_rho_sample() - Draws a reliability-stratified sample of genes for the intrinsic/extrinsic calibration check.
###     allele_cor_boot_se() - Bootstrap SE of one sampled gene's own allele-residual correlation.
###     check_intrinsic_reliability() - Tests whether the intrinsic/extrinsic decomposition's disattenuation behaves as the NB model predicts.
###     partial_cor_depth() - Correlation between two residual vectors with per-cell depth partialled out.
###     class_mean_var() - Mean and variance of a continuous score within each level of a class vector.
###     frac_group_sets() - Splits genes into low/average/high groups by a continuous score, for downstream GO enrichment by group.
###     class_anova() - One-way ANOVA testing whether a continuous score differs across class levels.
###     coexpr_class_table() - Pair-level regulatory classification (five-way) built from cis and trans p-values, BH-adjusted.
###     coexpr_dom_class_table() - Pair-level dominance classification, the co-expression analog of classify_dom().
###     coexpr_gene_degree() - Per-gene degree in the co-expression divergence network, how many.
###     plot_coexpr_cis_trans() - Cis vs trans scatter for co-expression pairs, coloured by pair-level regulatory class.
###     plot_coexpr_dom_class() - Dpar_sc vs dpar_se scatter for co-expression pairs, coloured by pair-level dominance class.
###     plot_coexpr_pair() - Residual scatter for one gene pair in each parent dataset, with per-dataset correlation in the legend.
###     check_coexpr_pairs() - Confirms a loaded CB object's pair count matches CO.GENES, catching a stale cluster output file.
###     seed_compare_core() - Shared scatter/correlation/ratio core for the two two-seed adequacy checks below.
###     coexpr_seed_compare() - Two-seed adequacy check for the co-expression bootstrap: compares CB between two independent seeds.
###     gene_seed_compare() - Two-seed adequacy check for the per-gene bootstrap: compares SE between two independent seeds.
###     coexpr_rank_check() - Rank-k eigendecomposition of a divergence matrix, with reconstruction R^2 against the observed values.
###     coexpr_candidate_axes() - Identifies candidate axes from a rank check: variance share, participation ratio, and floor screening.
###     coexpr_axis_mixtures() - Two-component Gaussian mixture fit and GO enrichment for each candidate axis's gene loadings.
###     coexpr_axis_validate() - Permutation-null validation of candidate axes, rank-matched against their own null spectrum.
###   8. REGULATORY AND DOMINANCE CLASSIFICATION  (offset level)
###     classify_reg() - Five-way regulatory classification (Conserved/Cis/Trans/Cis+Trans/Compensatory) from cis and trans p-values.
###     classify_dom() - Six-way dominance classification (hybrid vs each parent) from dpar_sc and dpar_se p-values.
###     clean_reg() - Collapses "Cis x Trans" into "Compensatory" and drops "Ambiguous" from a class vector.
###     reg_class_vec() - Vector form of classify_reg(), aligned to BURST.CONTRASTS's rows, for a given quantity.
###     dom_class_vec() - Vector form of classify_dom(), aligned to BURST.CONTRASTS's rows, for a given quantity.
###     class_overlap_heatmap() - Shared contingency-table/statistics core for the two class-overlap heatmaps below.
###     class_heatmap() - Log2(observed/expected) association heatmap between two class vectors, with BH-adjusted significance (levels auto-derived from data).
###     .heatmap_legend() - Draws class_heatmap()'s color-to-value legend strip beside the heatmap.
###     class_identity_overlap() - Cohen's kappa plus per-class Jaccard overlap between two class vectors on the same genes.
###     .cohen_kappa() - Cohen's kappa for a class-by-class contingency table, used by class_identity_overlap().
###     summarize_class_overlap() - One printable summary row (n, concordance, kappa, permutation p) per class-overlap comparison.
###     plot_cis_trans_class() - Cis vs trans scatter for one quantity, coloured by class, with SE bars drawn behind points.
###     plot_mean_bfreq_class() - Mean vs noise scatter with per-class regulatory and dominance slopes overlaid.
###   8b. GO / KEGG enrichment for the regulatory classification
###     build_reg_go_sets() - Builds one GO gene set per (class, direction) pair for a given classification quantity.
###     build_component_go_sets() - Builds pooled any-cis or any-trans GO gene sets, direction-matched on that component's own sign.
###     run_enrichment() - Runs GO (BP/MF/CC, simplified) and KEGG enrichment for one gene set against a fixed universe.
###     n_sig_terms() - Counts significant terms (q < threshold) in one enrichResult, 0 for a NULL or empty result.
###     summarize_go_sets() - One summary row per gene set: size plus significant-term counts across BP/MF/CC/KEGG.
###     plot_geneset_direction_stack() - Stacked barplot showing each class's membership split by direction within a gene set.
###     cluster_marker_enrichment() - Marker/enrichment comparison of each cluster against the rest of a dataset's own cells, generalized to however many clusters it has.
###     plot_cluster_marker_enrichment() - Writes cluster_marker_enrichment()'s per-cluster up/down enrichment to one pdf via barplot_enrich_pair().
###     report_cluster_marker_enrichment() - Console-text mirror of plot_cluster_marker_enrichment(): significant terms and p-values per cluster and direction, via print_enrich_brief().
###     score_cell_cycle_by_cluster() - Cell-cycle phase scoring (CellCycleScoring()/AddModuleScore()) extended to a dataset's own validated clustering; violin plot, phase-composition-by-cluster barplot, and console table.
###     go_gene_set() - Pulls all genes (including descendant terms) annotated to a GO term from org.Sc.sgd.db, for building an independently-sourced module-score gene set.
###     score_modules_by_cluster() - Continuous module-score validation (AddModuleScore()) of a cluster's marker-based identity against named, independently-sourced gene sets; violin plot, per-cluster mean +/- SE barplot, and console table.
###     cell_cycle_continuum() - Single continuous cell-cycle axis, first PC of S.Score/G2M.Score.
###     cell_cycle_continuum_shared() - Two-object version of cell_cycle_continuum(), pooled PCA so between-object means are comparable.
###     metabolic_state_cluster() - Discrete metabolic state from joint k-means clustering on the Glycolysis/OXPHOS/RiBi module scores, k chosen by silhouette.
###     covariate_noise_diagnostic() - Per-gene test of whether NB Pearson residual noise depends on cell-cycle position (Spearman) or metabolic state (Kruskal-Wallis, eta-squared).
###     plot_covariate_noise_diagnostic() - Histogram pair of covariate_noise_diagnostic()'s per-gene effect sizes.
###     species_composition_bound() - Cohen's d and rank-sum test between two species/allele views on one continuous axis, combined with a noise-association effect size into a bound on expected noise shift.
###     species_composition_report() - Runs species_composition_bound() across the cell-cycle axis and the three metabolic module scores for one species/allele pair.
###   9. PUBLICATION FIGURES
###     plot_class_overlap() - Mean-class x size-class enrichment heatmap (log2 obs/exp, BH-adjusted significance, explicit shared levels).
###     kbal_sig() - Bootstrap z-test on kinetic balance (burst frequency minus size) for one mode.
###     plot_burst_kinetics_sig() - Barplot of kbal_sig()'s z-test on kinetic balance, shaded by significance.
###     plot_dom_class() - Dominance scatter in the parent frame (or A/D rotation), coloured by dominance class.
###     plot_violins() - ggplot2 violin plot of a burst quantity, split by regulatory and dominance class.
###     intrinsic_fraction() - Two-allele intrinsic/extrinsic noise decomposition, Poisson-shot-noise corrected.
###     plot_intrinsic_hist() - Histogram of intrinsic fraction with regulatory and dominance class medians marked.
###   10. PROMOTER ARCHITECTURE (TATA box, poly(dA:dT), nucleosome occupancy)
###     read_genome_fasta() - Reads a genome FASTA into a named DNAStringSet, one sequence per chromosome.
###     read_gff_genes() - Reads a GFF3 annotation into a per-gene coordinate table.
###     extract_promoters() - Extracts each gene's promoter sequence, bounded by its upstream neighbor.
###     tata_box_score() - Scores a promoter sequence's best TATA-box PWM match in the canonical location window.
###     poly_at_tract() - Finds the longest poly(dA:dT) run in a promoter sequence.
###     score_promoters() - Applies both the TATA PWM score and poly(dA:dT) tract length to every promoter in a set.
###     nupop_cluster_inputs() - Packages one species' chromosomes and promoter coordinates for the NuPoP cluster job.
###     nupop_predict_window() - Predicts NuPoP occupancy for one sequence window in its own temporary folder.
###     nupop_run_windows() - Runs a table of flanked windows, one forked process per window.
###     nupop_occupancy_cluster() - Tiles, scores and bisects every chromosome on the cluster and returns promoter occupancy tracks.
###     nupop_window_score() - Reduces each gene's full NuPoP occupancy track to one comparable window score.
###     score_promoters_nupop() - Scores each promoter from the cluster occupancy tracks after confirming they match the current promoters.
###     promoter_divergence() - Merges Sc and Se promoter scores on gene identity and computes the between-species difference.
###     promoter_direction_test() - Tests whether a promoter feature's between-species direction matches the direction of cis divergence.
###     concordance_by_magnitude() - Splits an any-cis gene set into magnitude bins and computes concordance within each bin.
###     plot_concordance_by_magnitude() - Bar plot of concordance_by_magnitude()'s output, one bar per magnitude bin.
###     promoter_noise_candidates() - Flags genes with a significant cis component and a concordant promoter-architecture shift, for manual inspection.
###     architecture_noise_check() - Tests whether a promoter feature associates with a species' own noise (DISP), controlling for mean expression.
###   11. SEURAT CLUSTERING DIAGNOSTICS (data-driven feature count, resolution, and metric choices)
###     hvg_elbow() - Elbow point on the ranked standardized-variance curve from FindVariableFeatures(), for a data-driven nfeatures.
###     plot_hvg_elbow() - Diagnostic plot for hvg_elbow()'s output: full curve, chosen cutoff, and floor marked.
###     sweep_cluster_resolution() - Resolution sweep with a minimum-cluster-size guard; picks the coarsest resolution near-maximal silhouette width.
###     plot_resolution_sweep() - Diagnostic plot for sweep_cluster_resolution()'s output: silhouette vs. resolution, guard-excluded points and the chosen resolution marked.
###     bootstrap_cluster_stability() - Bootstrap resampling of cells at a chosen resolution; mean adjusted Rand index against the original clustering.
###     bootstrap_compare_resolutions() - Bootstrap-compares the chosen resolution against the grid's highest-silhouette plateau; returns the bootstrap-validated final resolution and reclustered object.
###     report_bootstrap_compare() - Prints bootstrap_compare_resolutions()'s table and logs a resolution override when the bootstrap-validated final resolution differs from the one originally chosen.
###   12. CLUSTER-BASED NOISE PARTITIONING (within/between-cluster variance vs. the intrinsic/extrinsic decomposition)
###     within_between_decomp() - Within- and between-cluster variance per gene, in shot-noise-corrected, mean-normalized rate space, from a Seurat cluster partition.
###     plot_within_between_hist() - Genome-wide distribution of within_between_decomp()'s ratio, one dataset, with the within = between line marked.
###     plot_within_between_vs_quantity() - Scatter of within_between_decomp()'s ratio against one burst kinetics quantity (mean, burst frequency, or burst size), with Spearman rho reported.
###     plot_within_between_cross_species() - Cross-species scatter of within_between_decomp()'s ratio, one point per ortholog gene pair, with Spearman rho reported.
###     plot_within_between_by_class() - Boxplot of within_between_decomp()'s ratio split by regulatory or dominance class, with class_anova()'s omnibus test and Tukey pairwise comparisons.
###     report_within_between_by_class() - Runs plot_within_between_by_class() for two datasets side by side against one classification axis; writes the combined figure and prints both datasets' test statistics.
###     intrinsic_extrinsic_components() - Raw intr/extr components underlying intrinsic_fraction(), exposed separately for use as a ratio.
###     plateau_coarsest() - Coarsest resolution in the contiguous same-cluster-count run surrounding the grid's argmax.
###     compare_distance_metrics() - Clusters the same PCs under two annoy.metric choices at the same resolution; adjusted Rand index between them.
###   13. POWER ANALYSIS (simulation and fit, shared with the SLURM job power_grid.R)
###     fit_offset_nb() - Offset NB fit for one simulated gene: exposure-weighted mean plus MLE dispersion via direct log-likelihood optimization.
###     fit_offset_nb_mm() - Closed-form method-of-moments dispersion estimate, used for permutation-null replicates.
###     fit_split_nb() - Splits a pooled simulated sample by a pre-drawn permutation and fits each half with fit_offset_nb().
###     fit_split_nb_mm() - Same split, fit with the method-of-moments estimator.
###     perm_pval() - Two-sided permutation p-value with add-one continuity correction.
###     size_log2_ratio() - Log2 ratio of SIZE (burst frequency) between two fits, NA if either side is non-positive or non-finite.
###     power_grid_row() - Power for one (MEAN.READS, N.CELLS, SIZE) row, across every SIZE.RATIO value at once.
###     open_grid_pdf() - Opens a PDF sized to the panel grid it is about to hold, with tight margins throughout.
###     line_colors() - Color ramp sized to the number of lines drawn in one panel.
###     plot_lines() - Plots one line per column of a matrix against a shared x vector, with reference lines at power 0.05 and 0.9.
###     legend_page() - One legend page mapping each line color to the value it represents.
###     min_detectable_ratio() - Minimum SIZE.RATIO reaching a target power, for the minimum-detectable-ratio summary heatmap.
###   14. EXTERNAL NOISE VALIDATION
###     cor_row() - Spearman rank correlation between two vectors, with pairwise-complete filtering.
###     add_burst_terms() - Implied Fano factor, burst size, and burst frequency algebraically recovered from a reported mean and CV^2.
###############################################################

library(MASS)   # glm.nb, ships with base R

## ckpt_path(n, dir): full path of the section n checkpoint. The default dir
## is read when the function is called, so every section names its checkpoint
## with the same section{N}_checkpoint.rda convention.
ckpt_path <- function(n, dir = CHECKPOINT.DIR) {
  file.path(dir, sprintf("section%d_checkpoint.rda", n))
}

## console_start(n, dir, append): copies everything printed to the console from
## this point on into section{n}_console.txt. split = TRUE keeps the text
## visible in the R console as well. Calling it for the next section closes the
## previous file, so each section gets its own record. append = TRUE adds to an
## existing file, which suits re-running a single subsection.
console_start <- function(n, dir = CONSOLE.DIR, append = FALSE) {
  console_stop()
  sink(file.path(dir, sprintf("section%d_console.txt", n)), append = append, split = TRUE)
}

## console_stop(): closes any open console file.
console_stop <- function() {
  while (sink.number() > 0) sink()
  invisible(NULL)
}

## ============================================================
## 1. OFFSET NEGATIVE-BINOMIAL FIT
## ============================================================
# Offset NB fit for one gene: mu from the exposure-weighted rate,
# disp (dispersion) from glm.nb, with a Poisson-vs-NB pre-check so a
# gene with no detectable overdispersion returns disp = Inf instead of
# an unstable glm.nb fit. init.theta seeds the optimizer when supplied.

neg_binom_fit_offset <- function(y, exposure, init.theta = NULL) {
  if (!is.null(dim(y)))
    stop("neg_binom_fit_offset() takes one gene's counts; use fit_counts_offset() for a matrix")
  na_out <- c(disp = NA_real_, mu = NA_real_, disp_logse = NA_real_, mu_logse = NA_real_)
  if (length(y) < 2 || length(exposure) != length(y)) return(na_out)
  keep <- is.finite(y) & is.finite(exposure) & exposure > 0
  y <- y[keep]; exposure <- exposure[keep]
  if (length(y) < 2) return(na_out)

  rate0 <- sum(y) / sum(exposure)
  if (sum(y) == 0) return(c(disp = NA_real_, mu = 0, disp_logse = NA_real_, mu_logse = NA_real_))

  mu_pois <- exposure * rate0
  pearson <- sum((y - mu_pois)^2 / mu_pois) / (length(y) - 1)
  if (pearson <= 1)
    return(c(disp = Inf, mu = rate0, disp_logse = NA_real_, mu_logse = NA_real_))

  logexp <- log(exposure)
  use_init <- !is.null(init.theta) && is.finite(init.theta) && 
              init.theta > 0 && init.theta < 1e6
  fit <- suppressWarnings(tryCatch(
    if (use_init) glm.nb(y ~ 1 + offset(logexp), init.theta = init.theta)
    else          glm.nb(y ~ 1 + offset(logexp)),
    error = function(e) NULL))
  if (is.null(fit) || isFALSE(fit$converged))
    return(c(disp = NA_real_, mu = rate0, disp_logse = NA_real_, mu_logse = NA_real_))

  theta  <- fit$theta
  lambda <- unname(exp(coef(fit)[1]))
  if (!is.finite(theta) || theta > 1e6)
    return(c(disp = Inf, mu = lambda, disp_logse = NA_real_, mu_logse = NA_real_))

  mu_logse   <- tryCatch(sqrt(vcov(fit)[1, 1]), error = function(e) NA_real_)
  disp_logse <- if (is.finite(fit$SE.theta)) fit$SE.theta / theta else NA_real_
  c(disp = theta, mu = lambda, disp_logse = disp_logse, mu_logse = mu_logse)
}

# Offset NB fit for one gene, used throughout the bootstrap and
# permutation machinery below (boot_contrasts_one, permute_contrasts_one):
# mu from the exposure-weighted rate, disp from direct optimization of
# the NB log-likelihood, with a Poisson-vs-NB pre-check
.fit_one <- function(counts, expo) {
  keep <- is.finite(counts) & is.finite(expo) & expo > 0
  y <- counts[keep]; e <- expo[keep]
  if (length(y) < 2) return(c(mu = NA_real_, disp = NA_real_))

  mu_hat <- sum(y) / sum(e)
  if (sum(y) == 0) return(c(mu = 0, disp = NA_real_))

  mu_i <- mu_hat * e
  pearson <- sum((y - mu_i)^2 / mu_i) / (length(y) - 1)
  if (pearson <= 1) return(c(mu = mu_hat, disp = Inf))

  ll <- function(ltheta) {
    th <- exp(ltheta)
    sum(lgamma(y + th) - lgamma(th) + th*log(th) - (th + y)*log(th + mu_i) + y*log(mu_i))
  }

  opt <- tryCatch(optimize(ll, c(-4, 15), maximum = TRUE), error = function(e) NULL)
  if (is.null(opt)) return(c(mu = mu_hat, disp = NA_real_))

  theta <- exp(opt$maximum)
  c(mu = mu_hat, disp = if (theta > 1e6) Inf else theta)
}

# Splits a pooled sample at position n1 using a pre-drawn permutation
# of indices and fits each half with .fit_one()
.fit_split <- function(counts, expo, perm, n1) {
  g1 <- perm[seq_len(n1)]; g0 <- perm[(n1 + 1):length(perm)]
  list(a = .fit_one(counts[g1], expo[g1]), b = .fit_one(counts[g0], expo[g0]))
}

fit_counts_offset <- function(mat, exposure) {
  stopifnot(ncol(mat) == length(exposure))
  n <- nrow(mat)
  size <- mu <- disp_logse <- mu_logse <- numeric(n)
  for (i in seq_len(n)) {
    f <- neg_binom_fit_offset(mat[i, ], exposure)
    size[i] <- f["disp"]; mu[i] <- f["mu"]
    disp_logse[i] <- f["disp_logse"]; mu_logse[i] <- f["mu_logse"]
  }
  data.frame(
    DISP = size, MU = mu, DISP_LOGSE = disp_logse, MU_LOGSE = mu_logse,
    MEAN_CT = rowSums(mat) / ncol(mat), N_EXPR = rowSums(mat > 0),
    VAR = mu + mu^2 / size, FANO = 1 + mu / size, CV = sqrt(1 / size + 1 / mu),
    BFREQ = size, BSIZE = mu / size,
    row.names = rownames(mat)
  )
}

## Bootstrap dispersion uncertainty over cells (one gene)
boot_disp_logse <- function(y, exposure, B = 200, seed = 1) {
  if (!is.null(dim(y))) stop("boot_disp_logse() takes one gene's counts")
  set.seed(seed); n <- length(y)
  lt <- replicate(B, log(neg_binom_fit_offset(y[sample.int(n, n, TRUE)], exposure[sample.int(n, n, TRUE)])["disp"]))
  finite <- is.finite(lt)
  c(disp_logse_boot = if (sum(finite) > 1) sd(lt[finite]) else NA_real_, boundary_frac   = mean(!finite))
}

## ============================================================
## 1b. INTERNAL-PILOT SPLIT FRACTION (f*)
## ============================================================
## Before the hybrid cells are partitioned into Nc/Nt, bootstrap the
## UNDIVIDED hybrid (HYB.SC, HYB.SE) and the parents (MIX.SC, MIX.SE)
## to get the per-gene variance components that set the correct
## partition. cis draws Nc = f*Nh hybrid cells; trans draws the
## remaining Nt = (1-f)*Nh hybrid cells plus the fixed Np parent
## cells, which is why trans always carries an extra, non-tunable
## variance term (B) that cis does not.
##
## This only uses variance information (A, B), never a cis or trans
## point estimate, so choosing f from it and then re-fitting the same
## data under that split is not the circular reuse that would matter,
## it is the standard internal-pilot / sample-size-re-estimation
## design, where nuisance-parameter information (not the effect being
## tested) sets the design before the real analysis runs

## One gene, bootstrap SE of log(mu) and log(size) for the four
## datasets needed to compute A and B. HYB.SC and HYB.SE are fit from
## the SAME resampled cell indices each replicate, since they are the
## two alleles measured in the same cells.
##
## Seeded per gene from a hash of its own name, not from position in
## a gene list. That matters here specifically because, unlike
## boot_contrasts_one and permute_contrasts_one, this function draws
## its own resamples internally rather than consuming a pre-built
## DRAWS list shared across genes, and it runs via parLapply on a
## FORK cluster, where genes get split across workers in chunks that
## can vary with core count. Seeding by name keeps the result
## reproducible regardless of chunk size or core count; seeding by
## position would not.
pilot_split_se_one <- function(g, mats, expos, B = 200, seed = 1) {
  set.seed(seed + sum(utf8ToInt(g)))
  n.p.sc <- length(expos$MIX.SC); n.p.se <- length(expos$MIX.SE); n.h <- length(expos$HYB)
  logmu <- matrix(NA_real_, B, 4, dimnames = list(NULL, c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")))
  logsz <- logmu
  for (b in seq_len(B)) {
    i.p.sc <- sample.int(n.p.sc, n.p.sc, replace = TRUE)
    i.p.se <- sample.int(n.p.se, n.p.se, replace = TRUE)
    i.h    <- sample.int(n.h,    n.h,    replace = TRUE)
    f.mc <- .fit_one(mats$MIX.SC[g, i.p.sc], expos$MIX.SC[i.p.sc])
    f.me <- .fit_one(mats$MIX.SE[g, i.p.se], expos$MIX.SE[i.p.se])
    f.hc <- .fit_one(mats$HYB.SC[g, i.h],    expos$HYB[i.h])
    f.he <- .fit_one(mats$HYB.SE[g, i.h],    expos$HYB[i.h])
    if (is.finite(f.mc["mu"])   && f.mc["mu"]   > 0) logmu[b, "MIX.SC"] <- log(f.mc["mu"])
    if (is.finite(f.me["mu"])   && f.me["mu"]   > 0) logmu[b, "MIX.SE"] <- log(f.me["mu"])
    if (is.finite(f.hc["mu"])   && f.hc["mu"]   > 0) logmu[b, "HYB.SC"] <- log(f.hc["mu"])
    if (is.finite(f.he["mu"])   && f.he["mu"]   > 0) logmu[b, "HYB.SE"] <- log(f.he["mu"])
    if (is.finite(f.mc["disp"]) && f.mc["disp"] > 0) logsz[b, "MIX.SC"] <- log(f.mc["disp"])
    if (is.finite(f.me["disp"]) && f.me["disp"] > 0) logsz[b, "MIX.SE"] <- log(f.me["disp"])
    if (is.finite(f.hc["disp"]) && f.hc["disp"] > 0) logsz[b, "HYB.SC"] <- log(f.hc["disp"])
    if (is.finite(f.he["disp"]) && f.he["disp"] > 0) logsz[b, "HYB.SE"] <- log(f.he["disp"])
  }
  data.frame(
    gene = g,
    MIX.SC_logmu_se   = sd(logmu[, "MIX.SC"], na.rm = TRUE),
    MIX.SE_logmu_se   = sd(logmu[, "MIX.SE"], na.rm = TRUE),
    HYB.SC_logmu_se   = sd(logmu[, "HYB.SC"], na.rm = TRUE),
    HYB.SE_logmu_se   = sd(logmu[, "HYB.SE"], na.rm = TRUE),
    MIX.SC_logdisp_se = sd(logsz[, "MIX.SC"], na.rm = TRUE),
    MIX.SE_logdisp_se = sd(logsz[, "MIX.SE"], na.rm = TRUE),
    HYB.SC_logdisp_se = sd(logsz[, "HYB.SC"], na.rm = TRUE),
    HYB.SE_logdisp_se = sd(logsz[, "HYB.SE"], na.rm = TRUE),
    row.names = NULL, check.names = FALSE)
}

## Serial convenience wrapper, useful for a quick local check on a
## handful of genes. The main pipeline no longer calls this directly;
## it calls pilot_split_se_one via parLapply in Run_Pilot_Cluster.R
## instead, the same way boot_contrasts_one and permute_contrasts_one
## are called directly rather than through their serial wrappers.
pilot_split_se <- function(genes, mats, expos, B = 200, seed = 1)
  do.call(rbind, lapply(genes, pilot_split_se_one, mats = mats, expos = expos, B = B, seed = seed))

## Closed-form balance point. With A the hybrid-driven variance
## coefficient and B the fixed, non-tunable parent-driven term,
## setting SE_cis(f) = SE_trans(f) gives a quadratic in f whose root
## in (0, 0.5] is f* = [(2+r) - sqrt(r^2+4)] / (2r), r = B*Nh/A.
## r -> 0 (parent term negligible) recovers f* -> 0.5, the even
## split, correctly, since there is then no asymmetry to correct for.
.fstar_from_r <- function(r) {
  out <- rep(0.5, length(r))
  ok <- is.finite(r) & abs(r) > 1e-8
  out[ok] <- ((2 + r[ok]) - sqrt(r[ok]^2 + 4)) / (2 * r[ok])
  out
}

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
estimate_f_star <- function(pilot, n_h) {
  A_mean <- n_h * (pilot$HYB.SC_logmu_se^2   + pilot$HYB.SE_logmu_se^2)
  B_mean <-        pilot$MIX.SC_logmu_se^2   + pilot$MIX.SE_logmu_se^2
  A_disp <- n_h * (pilot$HYB.SC_logdisp_se^2 + pilot$HYB.SE_logdisp_se^2)
  B_disp <-        pilot$MIX.SC_logdisp_se^2 + pilot$MIX.SE_logdisp_se^2

  ok_mean <- is.finite(A_mean) & is.finite(B_mean) & A_mean > 0
  ok_disp <- is.finite(A_disp) & is.finite(B_disp) & A_disp > 0

  r_mean <- median(B_mean[ok_mean]) * n_h / median(A_mean[ok_mean])
  r_disp <- median(B_disp[ok_disp]) * n_h / median(A_disp[ok_disp])

  list(f_mean = .fstar_from_r(r_mean), f_disp = .fstar_from_r(r_disp),
       r_mean = r_mean, r_disp = r_disp,
       n_genes_mean = sum(ok_mean), n_genes_disp = sum(ok_disp))
}

## Depth-matched split of depth-ordered hybrid cells at an arbitrary
## fraction f. Generalizes the old strict every-other-cell
## alternation (the f = 0.5 special case) so both groups stay
## depth-matched at any f, rather than assigning contiguous blocks
## from a depth-sorted order, which would confound the split with
## depth.
split_indices_by_depth <- function(n, f) {
  n_c <- max(1, round(f * n))
  c_pos <- unique(round(seq(1, n, length.out = n_c)))
  list(c = c_pos, t = setdiff(seq_len(n), c_pos))
}

## Single-gene version of fit_counts_offset, for parLapply across the
## cluster. Reassembling per-gene results in row order reproduces
## exactly what fit_counts_offset returns in one serial pass; this
## exists only to spread the same NB fits across cores instead of
## running them one core at a time locally.
## DISP here is the NB dispersion parameter (theta), the same quantity as
## burst frequency (k), NOT burst size. This is the naming trap in this
## pipeline: DISP and BSIZE sound like a pair but are not one, DISP
## equals BFREQ (defined below, a bare copy). BSIZE is the separate,
## derived quantity, MU / DISP (mean per burst), and is not equal to
## DISP under any name.
fit_counts_offset_row <- function(i, mat, exposure) {
  f <- neg_binom_fit_offset(mat[i, ], exposure)
  y <- mat[i, ]
  data.frame(
    DISP = unname(f["disp"]), MU = unname(f["mu"]),
    DISP_LOGSE = unname(f["disp_logse"]), MU_LOGSE = unname(f["mu_logse"]),
    MEAN_CT = sum(y) / length(y), N_EXPR = sum(y > 0),
    VAR = unname(f["mu"]) + unname(f["mu"])^2 / unname(f["disp"]),
    FANO = 1 + unname(f["mu"]) / unname(f["disp"]),
    CV = sqrt(1 / unname(f["disp"]) + 1 / unname(f["mu"])),
    BFREQ = unname(f["disp"]), BSIZE = unname(f["mu"]) / unname(f["disp"]),
    row.names = rownames(mat)[i])
}

## Runs fit_counts_offset_row across all genes on an existing cluster
## `cl`, in place of the serial for-loop inside fit_counts_offset
fit_counts_offset_parallel <- function(mat, exposure, cl) {
  rows <- parLapply(cl, seq_len(nrow(mat)), fit_counts_offset_row, mat = mat, exposure = exposure)
  do.call(rbind, rows)
}

## Fit-free preliminary gene filter, using only raw count summaries
## (mean count per cell, fraction of cells with nonzero counts) never
## the NB dispersion itself. This exists so the pilot step needs zero
## NB fits before it goes to the cluster; requiring a fit here would
## reintroduce the exact cost the pilot is meant to avoid paying
## locally. Used only to keep the pilot's gene pool reasonable, not
## as the final gene filter (build_gene_sets, using real fits, still
## does that after the split-dependent fits come back).
raw_gene_prefilter <- function(mats, min_mean = 0.001, min_expr_frac = 0.10) {
  ok <- Reduce(`&`, lapply(mats, function(m) {
    (rowSums(m) / ncol(m)) >= min_mean & (rowSums(m > 0) / ncol(m)) >= min_expr_frac
  }))
  rownames(mats[[1]])[ok]
}

## ============================================================
## 2. GENE FILTER (marginal information only; never on a contrast)
## ============================================================

gene_pass_group <- function(fit, ncells, min_mean = 1, min_expr_frac = 0.10) {
  is.finite(fit$DISP) & fit$DISP < 1e6 &
    fit$MEAN_CT >= min_mean & (fit$N_EXPR / ncells) >= min_expr_frac
}

build_gene_sets <- function(fits, ncells, min_mean = 1, min_expr_frac = 0.10) {
  groups <- names(fits)
  stopifnot(all(groups %in% names(ncells)))
  genes <- rownames(fits[[1]])
  for (g in groups)
    if (!identical(rownames(fits[[g]]), genes))
      stop("fit frames are not gene-aligned; reorder to a common gene set first")
  pass <- vapply(groups, function(g) gene_pass_group(fits[[g]], ncells[[g]], min_mean, min_expr_frac), logical(length(genes)))
  rownames(pass) <- genes
  contrasts <- list(total = c("MIX.SC","MIX.SE"), cis = c("HYC.SC","HYC.SE"),
                    trans = c("MIX.SC","MIX.SE","HYT.SC","HYT.SE"),
                    dom = c("MIX.SC","MIX.SE","HYB.COMB"),
                    dpar_sc = c("HYB.COMB","MIX.SC"), dpar_se = c("HYB.COMB","MIX.SE"),
                    inh_sc = c("HYB.SC","MIX.SC"), inh_se = c("HYB.SE","MIX.SE"),
                    full = groups)
  sets <- lapply(contrasts, function(gr) genes[rowSums(pass[, gr, drop = FALSE]) == length(gr)])
  list(pass = pass, sets = sets, n = vapply(sets, length, integer(1)), contrasts = contrasts)
}

refine_by_boundary <- function(mats, expos, candidates, max_boundary = 0.10, B = 200) {
  groups <- names(mats)
  keep <- vapply(candidates, function(gn) all(vapply(groups, function(g) {
    bf <- boot_disp_logse(mats[[g]][gn, ], expos[[g]], B = B)["boundary_frac"]
    is.finite(bf) && bf <= max_boundary
  }, logical(1))), logical(1))
  candidates[keep]
}

## Calibration helpers, run on one group's fit frame. chk() bins by
## mean count and reports the fraction with a finite size, the
## identifiability gradient. chk_prec() reports the median size
## standard error per bin; note this is the ASYMPTOTIC SE and is the
## unreliable ruler, the definitive precision gradient comes from the
## bootstrap size_*_se in BURST.CONTRASTS after step 3.
chk <- function(fit) {
  b <- cut(fit$MEAN_CT, c(0, 0.5, 1, 2, 4, 8, Inf))
  round(tapply(is.finite(fit$DISP) & fit$DISP < 1e6, b, mean), 2)
}
chk_prec <- function(fit) {
  b <- cut(fit$MEAN_CT, c(0, 0.5, 1, 2, 4, 8, Inf))
  round(tapply(fit$DISP_LOGSE, b, median, na.rm = TRUE), 2)
}

## ============================================================
## 3. PAIRED PER-GENE BOOTSTRAP OF CONTRASTS
## ============================================================

make_draws <- function(ncells, B, seed = 1) {
  set.seed(seed)
  lapply(seq_len(B), function(b) list(
    MIX.SC = sample.int(ncells[["MIX.SC"]], ncells[["MIX.SC"]], replace = TRUE),
    MIX.SE = sample.int(ncells[["MIX.SE"]], ncells[["MIX.SE"]], replace = TRUE),
    HYC    = sample.int(ncells[["HYC.SC"]], ncells[["HYC.SC"]], replace = TRUE),
    HYT    = sample.int(ncells[["HYT.SC"]], ncells[["HYT.SC"]], replace = TRUE),
    HYC.N  = sample.int(ncells[["HYC.SC.N"]], ncells[["HYC.SC.N"]], replace = TRUE),
    HYT.N  = sample.int(ncells[["HYT.SC.N"]], ncells[["HYT.SC.N"]], replace = TRUE),
    HYB    = sample.int(ncells[["HYB.COMB"]], ncells[["HYB.COMB"]], replace = TRUE)))
}

## Contrast definitions. Each mode lists the groups it needs (for the
## finite-input gate) and a formula on a named list of fitted values.
## Divergence: total, cis, trans. Dominance, midparent form: dom
## (arithmetic midparent for mean, geometric for burst frequency).
## Classic hybrid-vs-each-parent: dpar_sc, dpar_se. Per-allele
## inheritance (Option B): inh_sc, inh_se. Add an axis by adding one
## row here.
##
## cis_n and trans_n are internal only, not reported directly. They
## are cis and trans computed from the HYC.N/HYT.N split (chosen to
## balance noise-axis power) rather than the HYC/HYT split (chosen to
## balance mean-axis power). A single output column pair, bfreq_cis
## and bfreq_trans, draws its value from these instead of from cis and
## trans, since the noise axis needs its own split, not the mean
## axis's split reused. See .BFREQ_SOURCE and boot_contrasts_one.
.MODES <- list(
  total   = c("MIX.SC","MIX.SE"),
  cis     = c("HYC.SC","HYC.SE"),
  trans   = c("MIX.SC","MIX.SE","HYT.SC","HYT.SE"),
  cis_n   = c("HYC.SC.N","HYC.SE.N"),
  trans_n = c("MIX.SC","MIX.SE","HYT.SC.N","HYT.SE.N"),
  dom     = c("MIX.SC","MIX.SE","HYB.COMB"),
  dpar_sc = c("HYB.COMB","MIX.SC"),
  dpar_se = c("HYB.COMB","MIX.SE"),
  inh_sc  = c("HYB.SC","MIX.SC"),
  inh_se  = c("HYB.SE","MIX.SE"))

## Only these eight produce columns in BURST.CONTRASTS; cis_n/trans_n are sourced
## from, never reported under their own name
.OUT_MODES <- c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se")

## Which mode supplies the burst-frequency (noise axis) value for each
## reported mode. Identity for everything split-independent; cis and
## trans redirect to the noise-split version. mean_* always uses the
## mode's own name directly, never remapped, since MU (mean axis) is
## what the HYC/HYT split (f_mean) was chosen to balance.
.BFREQ_SOURCE <- c(total = "total", cis = "cis_n", trans = "trans_n", dom = "dom", dpar_sc = "dpar_sc", dpar_se = "dpar_se", inh_sc = "inh_sc", inh_se = "inh_se")

.contrast_value <- function(mode, gv, q) {
  l2 <- log2
  switch(mode,
    total   = l2(gv[["MIX.SC"]]) - l2(gv[["MIX.SE"]]),
    cis     = l2(gv[["HYC.SC"]]) - l2(gv[["HYC.SE"]]),
    cis_n   = l2(gv[["HYC.SC.N"]]) - l2(gv[["HYC.SE.N"]]),
    trans   = (l2(gv[["MIX.SC"]]) - l2(gv[["MIX.SE"]])) -
              (l2(gv[["HYT.SC"]]) - l2(gv[["HYT.SE"]])),
    trans_n = (l2(gv[["MIX.SC"]]) - l2(gv[["MIX.SE"]])) -
              (l2(gv[["HYT.SC.N"]]) - l2(gv[["HYT.SE.N"]])),
    dom     = l2(gv[["HYB.COMB"]]) -
              (if (q == "MU") l2((gv[["MIX.SC"]] + gv[["MIX.SE"]]) / 2)
               else           0.5 * (l2(gv[["MIX.SC"]]) + l2(gv[["MIX.SE"]]))),
    dpar_sc = l2(gv[["HYB.COMB"]]) - l2(gv[["MIX.SC"]]),
    dpar_se = l2(gv[["HYB.COMB"]]) - l2(gv[["MIX.SE"]]),
    inh_sc  = l2(gv[["HYB.SC"]]) - l2(gv[["MIX.SC"]]),
    inh_se  = l2(gv[["HYB.SE"]]) - l2(gv[["MIX.SE"]]))
}

## CV2 (squared coefficient of variation) for one group's fitted mu/disp:
## 1/mu + 1/disp. NA whenever either input is non-finite or non-positive,
## rather than propagating an Inf or NaN into the bootstrap SD.
.cv2_of <- function(z) {
  m <- unname(z[["mu"]]); k <- unname(z[["disp"]])
  if (is.finite(m) && m > 0 && is.finite(k) && k > 0) 1 / m + 1 / k else NA_real_
}

boot_contrasts_one <- function(g, mats, expos, fits, draws) {
  B <- length(draws); modes <- names(.MODES)
  mcol <- paste0("m.", modes); scol <- paste0("s.", modes); ccol <- paste0("c.", modes)
  M <- matrix(NA_real_, B, 3 * length(modes), dimnames = list(NULL, c(mcol, scol, ccol)))
  for (b in seq_len(B)) {
    d <- draws[[b]]
    fg <- function(mat, ex, idx) .fit_one(mat[g, idx], ex[idx])
    f <- list(
      MIX.SC   = fg(mats$MIX.SC,   expos$MIX.SC, d$MIX.SC),
      MIX.SE   = fg(mats$MIX.SE,   expos$MIX.SE, d$MIX.SE),
      HYC.SC   = fg(mats$HYC.SC,   expos$HYC,    d$HYC),
      HYC.SE   = fg(mats$HYC.SE,   expos$HYC,    d$HYC),
      HYT.SC   = fg(mats$HYT.SC,   expos$HYT,    d$HYT),
      HYT.SE   = fg(mats$HYT.SE,   expos$HYT,    d$HYT),
      HYC.SC.N = fg(mats$HYC.SC.N, expos$HYC.N,  d$HYC.N),
      HYC.SE.N = fg(mats$HYC.SE.N, expos$HYC.N,  d$HYC.N),
      HYT.SC.N = fg(mats$HYT.SC.N, expos$HYT.N,  d$HYT.N),
      HYT.SE.N = fg(mats$HYT.SE.N, expos$HYT.N,  d$HYT.N),
      HYB.COMB = fg(mats$HYB.COMB, expos$HYB,    d$HYB),
      HYB.SC   = fg(mats$HYB.SC,   expos$HYB,    d$HYB),
      HYB.SE   = fg(mats$HYB.SE,   expos$HYB,    d$HYB))
    gv.mu <- lapply(f, function(z) unname(z[["mu"]]))
    gv.bf <- lapply(f, function(z) unname(z[["disp"]]))
    gv.cv <- lapply(f, .cv2_of)
    for (k in seq_along(modes)) {
      grp <- .MODES[[modes[k]]]
      mu <- unlist(gv.mu[grp]); bf <- unlist(gv.bf[grp]); cv <- unlist(gv.cv[grp])
      if (all(is.finite(mu) & mu > 0)) M[b, mcol[k]] <- .contrast_value(modes[k], gv.mu, "MU")
      if (all(is.finite(bf) & bf > 0)) M[b, scol[k]] <- .contrast_value(modes[k], gv.bf, "BFREQ")
      if (all(is.finite(cv) & cv > 0)) M[b, ccol[k]] <- .contrast_value(modes[k], gv.cv, "CV2")
    }
  }
  se   <- apply(M, 2, sd, na.rm = TRUE)
  bdry <- setNames(vapply(modes, function(md) mean(is.na(M[, paste0("s.", md)])), numeric(1)), modes)
  point <- function(q, gv_fun) {
    gv <- setNames(lapply(names(fits), function(grp) gv_fun(fits[[grp]][g, ])), names(fits))
    vapply(modes, function(md) .contrast_value(md, gv, q), numeric(1))
  }
  pm <- point("MU",    function(row) row[["MU"]])
  ps <- point("BFREQ", function(row) row[["DISP"]])
  pc <- point("CV2",   function(row) { m <- row[["MU"]]; k <- row[["DISP"]]; if (is.finite(m) && m > 0 && is.finite(k) && k > 0) 1/m + 1/k else NA_real_ })
  out <- list(gene = g)
  for (md in .OUT_MODES) {
    smd <- .BFREQ_SOURCE[[md]]
    out[[paste0("mean_", md, "_est")]]   <- unname(pm[md])
    out[[paste0("mean_", md, "_se")]]    <- unname(se[paste0("m.", md)])
    out[[paste0("bfreq_", md, "_est")]]  <- unname(ps[smd])
    out[[paste0("bfreq_", md, "_se")]]   <- unname(se[paste0("s.", smd)])
    out[[paste0("bfreq_", md, "_bdry")]] <- unname(bdry[smd])
    out[[paste0("cv2_", md, "_est")]]    <- unname(pc[smd])
    out[[paste0("cv2_", md, "_se")]]     <- unname(se[paste0("c.", smd)])
    ## cor_md, the correlation between the mean and burst-frequency
    ## bootstrap draws, is only meaningful when both come from the same
    ## cell resampling stream. That holds for every mode except cis and
    ## trans, where mean now comes from the f_mean split and burst
    ## frequency from the independent f_disp split; the two streams
    ## share no resampled cells to be correlated through, so this is
    ## reported as NA rather than as a number that looks like the old
    ## quantity but no longer means the same thing.
    out[[paste0("cor_", md)]] <- if (identical(md, smd)) {
      x <- M[, paste0("m.", md)]; y <- M[, paste0("s.", md)]; ok <- is.finite(x) & is.finite(y)
      if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
    } else NA_real_
  }
  .r2 <- function(c1, c2) {
	x <- M[, c1]; y <- M[, c2]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
  }
  out$cor_dpar_mean  <- .r2("m.dpar_sc", "m.dpar_se")
  out$cor_dpar_bfreq <- .r2("s.dpar_sc", "s.dpar_se")
  data.frame(out, row.names = NULL, check.names = FALSE)
}

boot_contrasts <- function(mats, expos, fits, genes, draws)
  do.call(rbind, lapply(genes, boot_contrasts_one, mats = mats, expos = expos, fits = fits, draws = draws))

## Derives bsize (burst size) and kbal (kinetic balance) from the mean
## and burst-frequency contrasts boot_contrasts_one() already wrote into
## df. bfreq IS the raw NB dispersion parameter (theta, called "disp"
## internally in the fitting code, e.g. .fit_one()'s return field)
## reported under its biological name; bsize (mean transcripts per
## firing) is computed here, via log2(bsize) = log2(mean) - log2(bfreq).
## Internal fitting code below this point still uses "disp" as the
## per-group field name (matching R's own dnbinom()/glm.nb() "size"
## convention); "bfreq" is the name used everywhere the divergence
## contrasts themselves are read.
##
## kbal (kinetic balance, bfreq - bsize) is derived alongside bsize
## because bfreq and bsize are not two independent readings on a gene:
## bsize is defined as mean minus bfreq, so bfreq + bsize collapses
## back to the mean contrast already stored in df, and bfreq's own
## variance across genes typically dwarfs the mean's. When that holds,
## bsize is close to a rescaled, sign-flipped copy of bfreq (see the
## worked derivation kept with Section 4 below), so comparing bfreq to
## bsize directly, or classifying and enriching on each separately,
## mostly re-tests that same variance imbalance rather than biology.
## kbal is the other half of the same rotation used in
## plot_burst_kinetics()/kbal_sig(): bfreq + bsize is exactly the mean
## contrast (already a column here), and bfreq - bsize is the part of
## that rotation not already captured by the mean, so kbal and mean
## together carry the same information as bfreq and bsize without the
## structural anti-correlation. Reading kbal against mean (does a gene
## with a bigger overall change lean toward more frequent bursts or
## bigger ones) is the well-posed version of the question "does
## divergence in frequency track divergence in size" that reading bfreq
## against bsize is not.
add_burst_contrasts <- function(df) {
  ## .OUT_MODES, not names(.MODES): this runs on the finished
  ## BOOT.CONTRASTS/BURST.CONTRASTS, which only ever has columns for the
  ## eight reportable modes. cis_n and trans_n are intermediate, internal
  ## working modes inside boot_contrasts_one, collapsed into cis and
  ## trans (via .BFREQ_SOURCE) before this function ever sees the
  ## data; iterating over the full .MODES list here would look for
  ## mean_cis_n_est and find nothing, since that column was never
  ## written out.
  for (ct in .OUT_MODES) {
    m  <- df[[paste0("mean_",ct,"_est")]];  sm <- df[[paste0("mean_",ct,"_se")]]
    s  <- df[[paste0("bfreq_",ct,"_est")]]; ss <- df[[paste0("bfreq_",ct,"_se")]]
    r  <- df[[paste0("cor_",ct)]]
    ## cor_<ct> is NA for cis and trans (see boot_contrasts_one): mean
    ## there is drawn from the f_mean HYC/HYT split and burst frequency
    ## from the independent f_disp split, two disjoint resampled cell
    ## sets that share no randomness to be correlated through. That is
    ## not an unknown correlation to propagate as NA, it is a covariance
    ## of exactly zero by construction, so it is substituted here before
    ## the SE propagation rather than left to poison bsize_se with NA.
    r[!is.finite(r)] <- 0
    bs    <- m - s
    bs_se <- sqrt(pmax(0, sm^2 + ss^2 - 2 * r * sm * ss))
    df[[paste0("bsize_",ct,"_est")]] <- bs
    df[[paste0("bsize_",ct,"_se")]]  <- bs_se
    ## covFB is the error covariance between the bfreq and bsize
    ## estimates (not their biological covariance across genes), needed
    ## to propagate the SE of their difference correctly; the same
    ## quantity plot_burst_kinetics() and kbal_sig() compute inline, kept
    ## here once so every mode gets it, not only "total".
    covFB <- r * sm * ss - ss^2
    df[[paste0("kbal_",ct,"_est")]] <- s - bs
    df[[paste0("kbal_",ct,"_se")]]  <- sqrt(pmax(0, ss^2 + bs_se^2 - 2 * covFB))
  }
  df
}

## ============================================================
## 4. DISATTENUATED MEAN-DISP COUPLING (population summary)
## Only the mean-size coupling is reported here. Burst kinetics
## are read per gene from the plots, not from this summary, since
## any second-moment burst statistic is fixed by the variance
## ratio and carries no independent information.
## ============================================================

eiv_components <- function(BURST.CONTRASTS, mode = c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se")) {
  mode <- match.arg(mode)
  X  <- BURST.CONTRASTS[[paste0("mean_",mode,"_est")]]; sx <- BURST.CONTRASTS[[paste0("mean_",mode,"_se")]]
  Y  <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_est")]]; sy <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_se")]]
  r  <- BURST.CONTRASTS[[paste0("cor_",mode)]]
  ## As in add_burst_contrasts: cor_<mode> is NA for cis/trans because
  ## mean and size are drawn from disjoint resampled cell partitions, a
  ## structural zero covariance rather than a missing value. Requiring
  ## is.finite(r) here would silently drop every cis/trans gene instead
  ## of just using r = 0 for them.
  r[!is.finite(r)] <- 0
  ok <- is.finite(X) & is.finite(Y) & is.finite(sx) & is.finite(sy)
  X <- X[ok]; Y <- Y[ok]; sx <- sx[ok]; sy <- sy[ok]; r <- r[ok]
  Vm  <- var(X)    - mean(sx^2)
  Vs  <- var(Y)    - mean(sy^2)
  Cms <- cov(X, Y) - mean(r * sx * sy)
  ## Vm or Vs can be negative when measurement noise dominates true
  ## variance for this mode, giving sqrt(Vm * Vs) = NaN by design; that
  ## NaN is the informative result (attenuation correction failed, not
  ## an error), so the routine "NaNs produced" warning is suppressed here
  ## rather than left to accumulate over thousands of bootstrap draws
  rho_ms <- suppressWarnings(Cms / sqrt(Vm * Vs))
  c(n = length(X), Vm = Vm, Vs = Vs, Cms = Cms, rho_mean_disp = rho_ms, rho_raw_mean_disp = cor(X, Y))
}

eiv_boot_ci <- function(BURST.CONTRASTS, mode, B = 2000, seed = 1, probs = c(0.025, 0.975), stats = c("rho_mean_disp")) {
  base <- eiv_components(BURST.CONTRASTS, mode); set.seed(seed); n <- nrow(BURST.CONTRASTS)
  draws <- replicate(B, eiv_components(BURST.CONTRASTS[sample.int(n, n, TRUE), , drop = FALSE], mode)[stats])
  if (length(stats) == 1) draws <- matrix(draws, nrow = 1, dimnames = list(stats, NULL))
  na_frac <- vapply(stats, function(s) mean(!is.finite(draws[s, ])), numeric(1))
  ## Vm and Vs are attenuation-corrected variances (raw variance minus
  ## mean squared SE), which go negative whenever measurement noise
  ## dominates the true signal for that mode; Cms / sqrt(Vm * Vs) is then
  ## NaN, both for the point estimate and for individual bootstrap
  ## draws. When most draws fail this way, the handful that survive by
  ## chance are not a representative resample and can land well outside
  ## the [-1, 1] range a correlation should have, so the CI is reported
  ## as NA rather than as an interval built from a small, biased subset.
  ci <- apply(draws, 1, quantile, probs = probs, na.rm = TRUE)
  ci[, na_frac > 0.5] <- NA_real_
  list(estimate = base[stats], ci = ci, na_frac = na_frac)
}

eiv_table <- function(BURST.CONTRASTS, modes = c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se"), B = 2000) {
  do.call(rbind, lapply(modes, function(m) {
    e <- eiv_components(BURST.CONTRASTS, m); b <- eiv_boot_ci(BURST.CONTRASTS, m, B = B)
    data.frame(mode = m, n = e["n"],
      rho_raw = round(e["rho_raw_mean_disp"], 3),
      rho_mean_disp = round(e["rho_mean_disp"], 3),
      ms_lo = round(b$ci[1], 3), ms_hi = round(b$ci[2], 3),
      ms_na = round(b$na_frac, 3), row.names = NULL)
  }))
}

## bfreq_bsize_structural() answers, for one mode, how much of the raw
## correlation between bfreq_est and bsize_est is forced by bsize's own
## definition (bsize = mean - bfreq) rather than telling us anything new
## about the genes. Since bsize is that subtraction, algebra alone gives
## Cov(bfreq, bsize) = Cov(bfreq, mean) - Var(bfreq), exactly, before any
## measurement-noise correction. Substituting Cov(bfreq, mean) = 0, i.e.
## the null case where mean divergence and burst-frequency divergence
## are genuinely uncorrelated across genes, gives the predicted
## correlation the identity alone would produce with no real biology
## contributing at all: rho_null = -Var(bfreq) / sqrt(Var(bfreq)*Var(bsize)).
## Comparing that to the actually observed (attenuation-corrected)
## correlation shows whether the raw bfreq-vs-bsize relationship
## reported anywhere in the pipeline is telling us more than the
## variance imbalance between mean and bfreq already guarantees; Vm,
## Vs, and Cms come straight from eiv_components(), so this adds no new
## estimation, only the algebra that turns those three numbers into a
## null for a bfreq/bsize comparison instead of a mean/bfreq one.
bfreq_bsize_structural <- function(BURST.CONTRASTS, mode) {
  e <- eiv_components(BURST.CONTRASTS, mode)
  Vm <- e["Vm"]; Vf <- e["Vs"]; Cmf <- e["Cms"]           # Vs/Cms here are bfreq's, not bsize's
  Vs_bsize   <- Vm + Vf - 2 * Cmf
  Cov_obs    <- Cmf - Vf
  Cov_null   <- -Vf
  rho_obs    <- suppressWarnings(Cov_obs  / sqrt(Vf * Vs_bsize))
  rho_null   <- suppressWarnings(Cov_null / sqrt(Vf * Vs_bsize))
  c(n = e["n"], Vm = Vm, Vf = Vf, Vs_bsize = Vs_bsize,
    rho_bfreq_bsize_observed = rho_obs, rho_bfreq_bsize_null = rho_null,
    excess_over_null = rho_obs - rho_null)
}

## ============================================================
## 5. PLOTS (gene-level, square symmetric panels, SE bars)
## ============================================================

.sym <- function(v, e) { m <- max(abs(c(v + e, v - e)), na.rm = TRUE); c(-m, m) }

## cis vs trans for one quantity: "mean", "bfreq", "bsize", "kbal", "cv2".
## Both axes share one symmetric range, so the dotted +/-45 lines
## (same-direction and compensatory) are meaningful.
plot_cis_trans <- function(BURST.CONTRASTS, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), main = NULL, lim = NULL, bar_col = "#88888855", pt_col = "#1F4E79") {
  quantity <- match.arg(quantity)
  lab <- c(mean = "mean", bfreq = "burst frequency", bsize = "burst size", kbal = "frequency-size balance", cv2 = "CV2")[quantity]
  cx <- BURST.CONTRASTS[[paste0(quantity, "_cis_est")]];  sx <- BURST.CONTRASTS[[paste0(quantity, "_cis_se")]]
  cy <- BURST.CONTRASTS[[paste0(quantity, "_trans_est")]]; sy <- BURST.CONTRASTS[[paste0(quantity, "_trans_se")]]
  ok <- is.finite(cx) & is.finite(cy) & is.finite(sx) & is.finite(sy)
  cx <- cx[ok]; cy <- cy[ok]; sx <- sx[ok]; sy <- sy[ok]
  if (is.null(lim)) { m <- max(abs(c(cx + sx, cx - sx, cy + sy, cy - sy)), na.rm = TRUE); lim <- c(-m, m) }
  if (is.null(main)) main <- lab
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = lim, ylim = lim, xlab = paste(lab, "cis (log2)"), ylab = paste(lab, "trans (log2)"), main = main)
  abline(0, 1, lty = 3, col = "grey65"); abline(0, -1, lty = 3, col = "grey65")
  abline(h = 0, v = 0, col = "grey45")
  segments(cx - sx, cy, cx + sx, cy, col = bar_col)
  segments(cx, cy - sy, cx, cy + sy, col = bar_col)
  points(cx, cy, pch = 16, cex = 0.5, col = pt_col)
}

## mean vs noise (size) for one mode: "cis", "trans", "total".
## Axes scaled independently, each symmetric about zero, since mean
## and size differ in spread. The visual tilt is not the true slope;
## the disattenuated coupling is in eiv_table. This shows the cloud.
plot_mean_bfreq <- function(BURST.CONTRASTS, mode = c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se"), main = NULL, bar_col = "#88888855", pt_col = "#2E6E4E") {
  mode <- match.arg(mode)
  x <- BURST.CONTRASTS[[paste0("mean_", mode, "_est")]]; sx <- BURST.CONTRASTS[[paste0("mean_", mode, "_se")]]
  y <- BURST.CONTRASTS[[paste0("bfreq_", mode, "_est")]]; sy <- BURST.CONTRASTS[[paste0("bfreq_", mode, "_se")]]
  ok <- is.finite(x) & is.finite(y) & is.finite(sx) & is.finite(sy)
  x <- x[ok]; y <- y[ok]; sx <- sx[ok]; sy <- sy[ok]
  if (is.null(main)) main <- paste0("mean vs noise: ", mode)
  m <- max(abs(c(x + sx, x - sx, y + sy, y - sy)), na.rm = TRUE); lim <- c(-m, m)
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = lim, ylim = lim, xlab = "mean (log2)", ylab = "dispersion (log2)", main = main)
  abline(0, 1, lty = 3, col = "grey65"); abline(0, -1, lty = 3, col = "grey65")
  abline(h = 0, v = 0, col = "grey45")
  segments(x - sx, y, x + sx, y, col = bar_col)
  segments(x, y - sy, x, y + sy, col = bar_col)
  points(x, y, pch = 16, cex = 0.5, col = pt_col)
}

## burst kinetics for one mode, ROTATED so the structural
## frequency/size anti-correlation no longer reads as a trend.
##   x = net mean change (mean contrast; bfreq + bsize collapses back
##       to this exactly, so it is read straight from the mean columns
##       rather than recomputed from bfreq and bsize)
##   y = kinetic balance (bfreq - bsize; right is frequency-led, left
##       is amplitude-led), read from the kbal columns add_burst_contrasts()
##       already derived
## Both columns already carry correctly propagated SEs (see
## add_burst_contrasts()), so this function only has to read and plot
## them; it no longer redoes the covFB algebra locally.
plot_burst_kinetics <- function(BURST.CONTRASTS, mode = c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se"), main = NULL, bar_col = "#88888855", pt_col = "#9C3848") {
  mode <- match.arg(mode)
  x <- BURST.CONTRASTS[[paste0("mean_", mode, "_est")]]; sx <- BURST.CONTRASTS[[paste0("mean_", mode, "_se")]]
  y <- BURST.CONTRASTS[[paste0("kbal_", mode, "_est")]]; sy <- BURST.CONTRASTS[[paste0("kbal_", mode, "_se")]]
  ok <- is.finite(x) & is.finite(y) & is.finite(sx) & is.finite(sy)
  x <- x[ok]; y <- y[ok]; sx <- sx[ok]; sy <- sy[ok]
  if (is.null(main)) main <- paste0("burst kinetics: ", mode)
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = .sym(x, sx), ylim = .sym(y, sy), xlab = "net mean change (log2)", ylab = "frequency - amplitude (kinetic balance, log2)", main = main)
  abline(h = 0, v = 0, col = "grey45")
  segments(x - sx, y, x + sx, y, col = bar_col)
  segments(x, y - sy, x, y + sy, col = bar_col)
  points(x, y, pch = 16, cex = 0.5, col = pt_col)
}

## significance-shaded histogram of a contrast, coloured by direction
sig_hist <- function(x, p, sig = 0.05, brk = 0.1, xlim, ylim, xlab, up = "darkred", dn = "#EE7600FF") {
  ## Breaks must span the full range of x or hist() errors ("some 'x' not
  ## counted; maybe 'breaks' do not span range of 'x'"), even though xlim
  ## only controls which part of that range is drawn -- breaks and the
  ## visible axis are two different things, and a fixed [-8, 8] silently
  ## assumed no contrast would ever exceed that. lo/hi are rounded outward
  ## to the nearest brk (with a small epsilon against floating-point
  ## rounding at the exact boundary) so bin edges stay on a consistent
  ## grid regardless of how far the tails extend; outliers are still
  ## counted correctly, they just fall outside the visible xlim window.
  rng <- range(x, na.rm = TRUE)
  lo  <- floor((rng[1] - 1e-9) / brk) * brk
  hi  <- ceiling((rng[2] + 1e-9) / brk) * brk
  b   <- seq(lo, hi, by = brk)
  hist(x, breaks = b, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "# of Genes", main = "")
  hist(x[p < sig & x > 0], breaks = b, xlim = xlim, ylim = ylim, col = up, add = TRUE)
  hist(x[p < sig & x < 0], breaks = b, xlim = xlim, ylim = ylim, col = dn, add = TRUE)
}

## ============================================================
## 6. PERMUTATION NULL  (per-gene significance for the contrasts)
## ============================================================
## Mirrors the bootstrap but shuffles labels to null each contrast.
## Per mode:
##   total       pool the two parents, reshuffle species labels
##   cis         swap the two alleles within each hybrid cell
##   trans       reshuffle parent vs hybrid context on each allele side
##               (nulls the parental-ratio vs hybrid-ratio difference)
##   dom         additive surrogate cells from random parent pairs,
##               compared to the observed midparent (additivity null)
##   dpar_sc/se  pool combined hybrid with one parent, relabel
##   inh_sc/se   pool one hybrid allele with its parent, relabel
## The first four are centred at zero under the null they test. dpar
## and inh test exchangeability; in V2 the offset normalises depth but
## not gene dose, so read those against the bootstrap interval and
## treat their permutation p as descriptive.

make_perms <- function(ncells, NPERM, seed = 1) {
  set.seed(seed)
  nSC <- ncells[["MIX.SC"]]; nSE <- ncells[["MIX.SE"]]
  nHYC <- ncells[["HYC.SC"]]; nHYT <- ncells[["HYT.SC"]]
  nHYC.N <- ncells[["HYC.SC.N"]]; nHYT.N <- ncells[["HYT.SC.N"]]
  nHYB <- ncells[["HYB.COMB"]]
  lapply(seq_len(NPERM), function(b) list(
    total     = sample.int(nSC + nSE),
    cis       = runif(nHYC) < 0.5,
    cis_n     = runif(nHYC.N) < 0.5,
    transSC   = sample.int(nSC + nHYT),
    transSE   = sample.int(nSE + nHYT),
    transSC_n = sample.int(nSC + nHYT.N),
    transSE_n = sample.int(nSE + nHYT.N),
    dom_i   = sample.int(nSC, nHYB, replace = TRUE),
    dom_j   = sample.int(nSE, nHYB, replace = TRUE),
    dparSC  = sample.int(nHYB + nSC),
    dparSE  = sample.int(nHYB + nSE),
    inhSC   = sample.int(nHYB + nSC),
    inhSE   = sample.int(nHYB + nSE)))
}

permute_contrasts_one <- function(g, mats, expos, fits, perms) {
  l2 <- function(x) log2(x); B <- length(perms); modes <- names(.MODES)
  Nm <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Ns <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Nc <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  gv   <- setNames(lapply(names(fits), function(grp) fits[[grp]][g, ]), names(fits))
  gvmu <- lapply(gv, function(z) z[["MU"]]); gvbf <- lapply(gv, function(z) z[["DISP"]])
  gvcv <- lapply(gv, function(z) { m <- z[["MU"]]; k <- z[["DISP"]]; if (is.finite(m) && m > 0 && is.finite(k) && k > 0) 1/m + 1/k else NA_real_ })
  obs.m <- vapply(modes, function(m) .contrast_value(m, gvmu, "MU"),    numeric(1))
  obs.s <- vapply(modes, function(m) .contrast_value(m, gvbf, "BFREQ"), numeric(1))
  obs.c <- vapply(modes, function(m) .contrast_value(m, gvcv, "CV2"),   numeric(1))
  cMIXsc <- mats$MIX.SC[g, ]; cMIXse <- mats$MIX.SE[g, ]
  cHYCsc <- mats$HYC.SC[g, ]; cHYCse <- mats$HYC.SE[g, ]
  cHYTsc <- mats$HYT.SC[g, ]; cHYTse <- mats$HYT.SE[g, ]
  cHYCsc.N <- mats$HYC.SC.N[g, ]; cHYCse.N <- mats$HYC.SE.N[g, ]
  cHYTsc.N <- mats$HYT.SC.N[g, ]; cHYTse.N <- mats$HYT.SE.N[g, ]
  cHYBsc <- mats$HYB.SC[g, ]; cHYBse <- mats$HYB.SE[g, ]; cHYBc <- mats$HYB.COMB[g, ]
  mp.m <- l2((gvmu[["MIX.SC"]] + gvmu[["MIX.SE"]]) / 2)
  mp.s <- 0.5 * (l2(gvbf[["MIX.SC"]]) + l2(gvbf[["MIX.SE"]]))
  mp.c <- 0.5 * (l2(gvcv[["MIX.SC"]]) + l2(gvcv[["MIX.SE"]]))
  nSC <- length(cMIXsc); nSE <- length(cMIXse); nHYB <- length(cHYBc)
  rat <- function(a, b, q) { x <- a[q]; y <- b[q]
    if (is.finite(x) && x > 0 && is.finite(y) && y > 0) l2(x) - l2(y) else NA_real_ }
  ## Same as rat(), but for CV2 = 1/mu + 1/disp, reconstructed from the
  ## two raw fields a .fit_one()/.fit_split() result carries (mu, disp)
  ## rather than read directly off one field.
  cv2_of <- function(f) { m <- f["mu"]; k <- f["disp"]; if (is.finite(m) && m > 0 && is.finite(k) && k > 0) 1/m + 1/k else NA_real_ }
  ratcv2 <- function(a, b) { x <- cv2_of(a); y <- cv2_of(b)
    if (is.finite(x) && x > 0 && is.finite(y) && y > 0) l2(x) - l2(y) else NA_real_ }
  for (b in seq_len(B)) {
    p <- perms[[b]]
    sp <- .fit_split(c(cMIXsc, cMIXse), c(expos$MIX.SC, expos$MIX.SE), p$total, nSC)
    Nm[b,"total"] <- rat(sp$a, sp$b, "mu"); Ns[b,"total"] <- rat(sp$a, sp$b, "disp"); Nc[b,"total"] <- ratcv2(sp$a, sp$b)

    sc <- ifelse(p$cis, cHYCse, cHYCsc); se <- ifelse(p$cis, cHYCsc, cHYCse)
    fa <- .fit_one(sc, expos$HYC); fb <- .fit_one(se, expos$HYC)
    Nm[b,"cis"] <- rat(fa, fb, "mu"); Ns[b,"cis"] <- rat(fa, fb, "disp"); Nc[b,"cis"] <- ratcv2(fa, fb)

    ## Noise-split (f_disp) version of the cis null, used only for
    ## the bfreq_cis/cv2_cis output columns, mirroring cis_n in .MODES
    sc.n <- ifelse(p$cis_n, cHYCse.N, cHYCsc.N); se.n <- ifelse(p$cis_n, cHYCsc.N, cHYCse.N)
    fa.n <- .fit_one(sc.n, expos$HYC.N); fb.n <- .fit_one(se.n, expos$HYC.N)
    Nm[b,"cis_n"] <- rat(fa.n, fb.n, "mu"); Ns[b,"cis_n"] <- rat(fa.n, fb.n, "disp"); Nc[b,"cis_n"] <- ratcv2(fa.n, fb.n)

    sSC <- .fit_split(c(cMIXsc, cHYTsc), c(expos$MIX.SC, expos$HYT), p$transSC, nSC)
    sSE <- .fit_split(c(cMIXse, cHYTse), c(expos$MIX.SE, expos$HYT), p$transSE, nSE)
    Nm[b,"trans"] <- rat(sSC$a, sSE$a, "mu")   - rat(sSC$b, sSE$b, "mu")
    Ns[b,"trans"] <- rat(sSC$a, sSE$a, "disp") - rat(sSC$b, sSE$b, "disp")
    Nc[b,"trans"] <- ratcv2(sSC$a, sSE$a)      - ratcv2(sSC$b, sSE$b)

    ## Noise-split version of the trans null, used only for the
    ## bfreq_trans/cv2_trans output columns, mirroring trans_n in .MODES
    sSC.n <- .fit_split(c(cMIXsc, cHYTsc.N), c(expos$MIX.SC, expos$HYT.N), p$transSC_n, nSC)
    sSE.n <- .fit_split(c(cMIXse, cHYTse.N), c(expos$MIX.SE, expos$HYT.N), p$transSE_n, nSE)
    Nm[b,"trans_n"] <- rat(sSC.n$a, sSE.n$a, "mu")   - rat(sSC.n$b, sSE.n$b, "mu")
    Ns[b,"trans_n"] <- rat(sSC.n$a, sSE.n$a, "disp") - rat(sSC.n$b, sSE.n$b, "disp")
    Nc[b,"trans_n"] <- ratcv2(sSC.n$a, sSE.n$a)      - ratcv2(sSC.n$b, sSE.n$b)

    synth <- cMIXsc[p$dom_i] + cMIXse[p$dom_j]
    es    <- expos$MIX.SC[p$dom_i] + expos$MIX.SE[p$dom_j]
    fs <- .fit_one(synth, es)
    if (is.finite(fs["mu"])   && fs["mu"]   > 0) Nm[b,"dom"] <- l2(fs["mu"])   - mp.m
    if (is.finite(fs["disp"]) && fs["disp"] > 0) Ns[b,"dom"] <- l2(fs["disp"]) - mp.s
    fs.cv <- cv2_of(fs)
    if (is.finite(fs.cv) && fs.cv > 0) Nc[b,"dom"] <- l2(fs.cv) - mp.c
    d1 <- .fit_split(c(cHYBc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$dparSC, nHYB)
    Nm[b,"dpar_sc"] <- rat(d1$a, d1$b, "mu"); Ns[b,"dpar_sc"] <- rat(d1$a, d1$b, "disp"); Nc[b,"dpar_sc"] <- ratcv2(d1$a, d1$b)
    d2 <- .fit_split(c(cHYBc, cMIXse), c(expos$HYB, expos$MIX.SE), p$dparSE, nHYB)
    Nm[b,"dpar_se"] <- rat(d2$a, d2$b, "mu"); Ns[b,"dpar_se"] <- rat(d2$a, d2$b, "disp"); Nc[b,"dpar_se"] <- ratcv2(d2$a, d2$b)
    i1 <- .fit_split(c(cHYBsc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$inhSC, nHYB)
    Nm[b,"inh_sc"] <- rat(i1$a, i1$b, "mu"); Ns[b,"inh_sc"] <- rat(i1$a, i1$b, "disp"); Nc[b,"inh_sc"] <- ratcv2(i1$a, i1$b)
    i2 <- .fit_split(c(cHYBse, cMIXse), c(expos$HYB, expos$MIX.SE), p$inhSE, nHYB)
    Nm[b,"inh_se"] <- rat(i2$a, i2$b, "mu"); Ns[b,"inh_se"] <- rat(i2$a, i2$b, "disp"); Nc[b,"inh_se"] <- ratcv2(i2$a, i2$b)
  }
  pval <- function(obs, null) { ok <- is.finite(null)
    if (!is.finite(obs) || sum(ok) < 1) return(NA_real_)
    (1 + sum(abs(null[ok]) >= abs(obs))) / (1 + sum(ok)) }
  out <- list(gene = g)
  for (md in .OUT_MODES) {
    smd <- .BFREQ_SOURCE[[md]]
    out[[paste0("mean_", md, "_obs")]]  <- unname(obs.m[md])
    out[[paste0("mean_", md, "_p")]]    <- pval(obs.m[md], Nm[, md])
    out[[paste0("bfreq_", md, "_obs")]] <- unname(obs.s[smd])
    out[[paste0("bfreq_", md, "_p")]]   <- pval(obs.s[smd], Ns[, smd])
    out[[paste0("cv2_", md, "_obs")]]   <- unname(obs.c[smd])
    out[[paste0("cv2_", md, "_p")]]     <- pval(obs.c[smd], Nc[, smd])

    ## Burst size (bsize = mean - bfreq in log2 space) needs its own null,
    ## built by differencing the mean and burst-frequency permutation draws
    ## index for index. For total/dom/dpar/inh, md == smd, so the mean
    ## draw and burst-frequency draw at a given b share the same
    ## relabeling and are correlated; differencing within b carries that
    ## correlation through correctly. For cis/trans, md != smd: the mean
    ## draw (b, md) comes from the f_mean-split permutation vector and the
    ## burst-frequency draw (b, smd) comes from the independent f_disp-
    ## split permutation vector in the same perms[[b]] entry. These two
    ## vectors are independently drawn (see make_perms), so they share no
    ## randomness to correlate through, and differencing them index by
    ## index still gives a valid null for the difference of two
    ## independent quantities -- any pairing across b would do, since the
    ## two sides are independent regardless of b.
    bs_null <- Nm[, md] - Ns[, smd]
    obs.bs  <- unname(obs.m[md]) - unname(obs.s[smd])
    out[[paste0("bsize_", md, "_obs")]] <- obs.bs
    out[[paste0("bsize_", md, "_p")]]   <- pval(obs.bs, bs_null)

    ## Kinetic balance (kbal = bfreq - bsize = 2*bfreq - mean) gets its
    ## own permutation null the same way bsize just did, recombining the
    ## identical Nm/Ns draws rather than resampling again: kbal_null is
    ## just bs_null's null reflected back through the same bfreq null
    ## draw it was built from. This gives kbal a real permutation p-value
    ## at every mode, so classify_reg()/classify_dom() can be run on it
    ## exactly as they already are on mean, bfreq, bsize, and cv2, and
    ## the mean-vs-kbal comparison (does a gene's overall divergence lean
    ## toward frequency-led or size-led) is testable without borrowing
    ## bfreq's own null a second time.
    kbal_null <- Ns[, smd] - bs_null
    obs.kbal  <- unname(obs.s[smd]) - obs.bs
    out[[paste0("kbal_", md, "_obs")]] <- obs.kbal
    out[[paste0("kbal_", md, "_p")]]   <- pval(obs.kbal, kbal_null)
  }
  data.frame(out, row.names = NULL, check.names = FALSE)
}

permute_contrasts <- function(mats, expos, fits, genes, perms)
  do.call(rbind, lapply(genes, permute_contrasts_one, mats = mats, expos = expos, fits = fits, perms = perms))

## ============================================================
## 7. CO-EXPRESSION  (residual co-fluctuation, cis/trans decomposed)
## ============================================================
## Network-level extension of the noise analysis. Correlate each
## gene's offset-NB Pearson residuals across cells, so the signal is
## correlated noise rather than shared mean or depth. Decompose the
## divergence in that correlation structure the same way as single
## genes: parents give total, the two hybrid alleles give cis, and
## trans is the remainder. Correlations are shrunk (Schaefer-Strimmer)
## and uncertainty comes from the same cell resampling as elsewhere.

## Pearson residuals from the offset NB fit, genes x cells
nb_residuals <- function(mat, exposure, fit) {
  fitted <- outer(fit$MU, exposure)
  v <- fitted + ifelse(is.finite(fit$DISP), fitted^2 / fit$DISP, 0)
  r <- (mat - fitted) / sqrt(v)
  r[!is.finite(r)] <- 0
  r
}

## Analytic shrinkage of a correlation matrix toward the identity.
## Z is cells x genes. Off-diagonals are pulled toward zero by the
## Schaefer-Strimmer optimal intensity.
shrink_cor <- function(Z) {
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  n <- nrow(Z); p <- ncol(Z)
  Zs <- scale(Z); Zs[!is.finite(Zs)] <- 0
  R  <- crossprod(Zs) / (n - 1)
  num <- 0; den <- 0
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    w <- Zs[, i] * Zs[, j]
    num <- num + (n / (n - 1)^3) * sum((w - mean(w))^2)
    den <- den + R[i, j]^2
  }
  lam <- if (den > 0) max(0, min(1, num / den)) else 1
  Rs <- R * (1 - lam); diag(Rs) <- 1
  attr(Rs, "lambda") <- lam
  Rs
}

## total / cis / trans divergence of the correlation structure, plus the
## dominance contrast when HYB.COMB is supplied. resid: list of residual
## matrices (genes x cells), same gene order, for MIX.SC, MIX.SE, HYB.SC,
## HYB.SE, and optionally HYB.COMB (both hybrid alleles summed, allele
## identity ignored, the co-expression analog of the single-gene "hybrid
## total" used in dpar_sc/dpar_se). dpar_sc/dpar_se compare the hybrid's
## own pairwise correlation (Rhyb) to each parent's separately, mirroring
## the single-gene dominance contrast at the pair level rather than the
## gene level. HYB.COMB is optional and computed conditionally, since
## coexpr_perm_one()'s permutation null never needs the dominance
## contrast and never supplies it.
coexpr_decompose <- function(resid) {
  Rsc  <- shrink_cor(t(resid$MIX.SC)); Rse  <- shrink_cor(t(resid$MIX.SE))
  Rhsc <- shrink_cor(t(resid$HYB.SC)); Rhse <- shrink_cor(t(resid$HYB.SE))
  total <- Rsc - Rse; cis <- Rhsc - Rhse
  out <- list(total = total, cis = cis, trans = total - cis,
             lambda = c(Sc = attr(Rsc, "lambda"), Se = attr(Rse, "lambda"),
                        HYB.SC = attr(Rhsc, "lambda"), HYB.SE = attr(Rhse, "lambda")))
  if (!is.null(resid$HYB.COMB)) {
    Rhyb <- shrink_cor(t(resid$HYB.COMB))
    out$dpar_sc <- Rhyb - Rsc
    out$dpar_se <- Rhyb - Rse
    out$lambda  <- c(out$lambda, HYB.COMB = attr(Rhyb, "lambda"))
  }
  out
}

## Raw per-dataset shrinkage correlation matrices, kept separately from
## coexpr_decompose() (which only keeps differences). This is for testing
## whether observed correlation MAGNITUDE is attenuated at low
## reliability, a bias question, distinct from whether bootstrap SE
## tracks reliability (a variance question, already tested and found
## flat within CO.GENES). Same shrink_cor() calls, nothing new computed.
coexpr_raw_cor <- function(resid) {
  list(Rsc  = shrink_cor(t(resid$MIX.SC)), Rse  = shrink_cor(t(resid$MIX.SE)), Rhsc = shrink_cor(t(resid$HYB.SC)), Rhse = shrink_cor(t(resid$HYB.SE)))
}

## Exact cis/trans decomposition of a given eigenvector's eigenvalue.
## total_mat = cis_mat + trans_mat elementwise, so for any vector v, the
## quadratic form v'(total)v splits exactly as v'(cis)v + v'(trans)v.
## More precise than correlating separately-computed dominant
## eigenvectors of cis and trans against a total-matrix eigenvector,
## since it's an identity on the same vector rather than a comparison
## between two different decompositions.
## v: a unit eigenvector (e.g. RANK.CHECK$vectors[, k]). PT: the
## coexpr_decompose() point estimate (COEXPR.POINT), with $cis/$trans/$total.
coexpr_axis_cis_trans <- function(v, PT) {
  cis_part   <- as.numeric(t(v) %*% PT$cis   %*% v)
  trans_part <- as.numeric(t(v) %*% PT$trans %*% v)
  total_part <- as.numeric(t(v) %*% PT$total %*% v)
  c(cis = cis_part, trans = trans_part, total = total_part, check_sum = cis_part + trans_part, frac_trans = trans_part / total_part)
}

## Cell-resampling draws for the co-expression bootstrap, generated once so
## a cluster job can parallelize over draws instead of looping serially.
## Mirrors make_draws() for the per-gene bootstrap. Both hybrid alleles
## share one resample per draw (paired), matching the pairing there too.
make_coexpr_draws <- function(nSC, nSE, nH, B, seed = 1) {
  set.seed(seed)
  lapply(seq_len(B), function(b) list(
    SC = sample.int(nSC, nSC, replace = TRUE),
    SE = sample.int(nSE, nSE, replace = TRUE),
    H  = sample.int(nH,  nH,  replace = TRUE)))
}

## One bootstrap draw of the co-expression decomposition. draw: one element
## of make_coexpr_draws(). resid: the RESID list (genes x cells per
## dataset), now including HYB.COMB, so this draw also returns the
## dominance contrasts (dpar_sc, dpar_se) alongside total/cis/trans.
## HYB.COMB is resampled with the same H draw as HYB.SC/HYB.SE, since
## it's the same cells, just the two alleles summed. Returns total/cis/
## trans/dpar_sc/dpar_se at the upper-triangle pair positions for this one
## draw; this is the unit of work a cluster worker does.
coexpr_bootstrap_one <- function(draw, resid) {
  d <- coexpr_decompose(list(
    MIX.SC = resid$MIX.SC[, draw$SC], MIX.SE = resid$MIX.SE[, draw$SE],
    HYB.SC = resid$HYB.SC[, draw$H],  HYB.SE = resid$HYB.SE[, draw$H],
    HYB.COMB = resid$HYB.COMB[, draw$H]))
  p  <- nrow(resid$MIX.SC)
  up <- which(upper.tri(matrix(0, p, p)))
  list(total = d$total[up], cis = d$cis[up], trans = d$trans[up], dpar_sc = d$dpar_sc[up], dpar_se = d$dpar_se[up])
}

## Draws for the joint total/cis/trans/dpar_sc/dpar_se permutation null
## (rank-matched axis testing), one cluster job producing all five null
## spectra at once from the same draw. idx: a random re-split of the
## pooled parent cells into pseudo-Sc/pseudo-Se groups of the original
## sizes (total's null). swap: an independent per-hybrid-cell coin flip
## of which allele is labeled Sc vs Se (cis's null). idx_dpar_sc: a
## random re-split of the pooled Sc-parent and allele-summed hybrid cells
## into pseudo-Sc/pseudo-hybrid groups of the original sizes, the same
## logic as total's null applied to "does condition (Sc vs hybrid)
## matter" instead of "does species matter". idx_dpar_se: the same,
## pooling Se-parent and hybrid cells instead. All four are drawn
## independently per draw, matching how total.null and cis.null are
## drawn independently so that trans.null = total.null - cis.null
## inherits the real total/cis independence correctly; dpar_sc and
## dpar_se have no such identity linking them to total/cis/trans or to
## each other, so there is no dependency to preserve between any of the
## five null spectra.
make_coexpr_perm_draws <- function(nSC, nSE, nH, B, seed = 1) {
  set.seed(seed)
  n_tot <- nSC + nSE
  lapply(seq_len(B), function(b) list(
    idx  = sample(n_tot),
    swap = sample(c(TRUE, FALSE), nH, replace = TRUE),
    idx_dpar_sc = sample(nSC + nH),
    idx_dpar_se = sample(nSE + nH)))
}

## One permutation draw of the total/cis/trans/dpar_sc/dpar_se null
## spectra, for rank-matched testing of every candidate axis (real rank k
## tested against the null's own rank k, not one shared threshold). draw:
## one element of make_coexpr_perm_draws(). resid: the RESID list, now
## including HYB.COMB. nSC, nSE: parent cell counts, used to split each
## pooled, reshuffled pool back into groups of the original sizes. n_keep:
## how many ranks to retain per decomposition (15, matching the top-15
## candidate window used throughout). Returns the top n_keep squared
## eigenvalue magnitudes, descending, for all five decompositions from
## this one draw; this is the unit of work a cluster worker does.
coexpr_perm_one <- function(draw, resid, nSC, nSE, n_keep = 15) {
  pooled   <- cbind(resid$MIX.SC, resid$MIX.SE)
  perm.sc  <- pooled[, draw$idx[seq_len(nSC)]]
  perm.se  <- pooled[, draw$idx[-seq_len(nSC)]]
  perm.hsc <- resid$HYB.SC; perm.hse <- resid$HYB.SE
  perm.hsc[, draw$swap] <- resid$HYB.SE[, draw$swap]
  perm.hse[, draw$swap] <- resid$HYB.SC[, draw$swap]

  d <- coexpr_decompose(list(MIX.SC = perm.sc, MIX.SE = perm.se, HYB.SC = perm.hsc, HYB.SE = perm.hse))
  topk <- function(m) {
    ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    (ev[order(abs(ev), decreasing = TRUE)][seq_len(n_keep)])^2
  }

  ## dpar_sc's null: pool Sc-parent and allele-summed hybrid cells,
  ## reshuffle into pseudo-Sc/pseudo-hybrid groups of the original sizes,
  ## so the null spectrum reflects what the largest eigenvalues of
  ## "condition A minus condition B" look like when condition doesn't
  ## actually distinguish which cells are which. dpar_se's null does the
  ## same with Se-parent and hybrid cells instead.
  pooled.sc <- cbind(resid$MIX.SC, resid$HYB.COMB)
  perm.a.sc <- pooled.sc[, draw$idx_dpar_sc[seq_len(nSC)]]
  perm.b.sc <- pooled.sc[, draw$idx_dpar_sc[-seq_len(nSC)]]
  dpar_sc_null <- topk(shrink_cor(t(perm.b.sc)) - shrink_cor(t(perm.a.sc)))

  pooled.se <- cbind(resid$MIX.SE, resid$HYB.COMB)
  perm.a.se <- pooled.se[, draw$idx_dpar_se[seq_len(nSE)]]
  perm.b.se <- pooled.se[, draw$idx_dpar_se[-seq_len(nSE)]]
  dpar_se_null <- topk(shrink_cor(t(perm.b.se)) - shrink_cor(t(perm.a.se)))

  list(total = topk(d$total), cis = topk(d$cis), trans = topk(d$trans), dpar_sc = dpar_sc_null, dpar_se = dpar_se_null)
}

## Assembles CB (per-pair est/se/z/p) from a point estimate
## (coexpr_decompose(resid)) and a list of per-draw results, one element
## per B from coexpr_bootstrap_one(). Same est/se/z/p structure the old
## inline loop produced, just built from pre-computed draw results so this
## step is cheap regardless of where the draws themselves were run.
assemble_coexpr_bootstrap <- function(resid, pt, draw_results) {
  p  <- nrow(resid$MIX.SC); gn <- rownames(resid$MIX.SC)
  up <- which(upper.tri(matrix(0, p, p)))
  ij <- arrayInd(up, c(p, p))
  base <- data.frame(gene_i = gn[ij[, 1]], gene_j = gn[ij[, 2]])

  acc <- list(total   = do.call(rbind, lapply(draw_results, `[[`, "total")),
              cis     = do.call(rbind, lapply(draw_results, `[[`, "cis")),
              trans   = do.call(rbind, lapply(draw_results, `[[`, "trans")),
              dpar_sc = do.call(rbind, lapply(draw_results, `[[`, "dpar_sc")),
              dpar_se = do.call(rbind, lapply(draw_results, `[[`, "dpar_se")))

  mk <- function(est, boot) { se <- apply(boot, 2, sd, na.rm = TRUE)
    data.frame(base, est = est, se = se, z = est / se, p = 2 * pnorm(-abs(est / se))) }

  list(total   = mk(pt$total[up],   acc$total),
       cis     = mk(pt$cis[up],     acc$cis),
       trans   = mk(pt$trans[up],   acc$trans),
       dpar_sc = mk(pt$dpar_sc[up], acc$dpar_sc),
       dpar_se = mk(pt$dpar_se[up], acc$dpar_se), lambda = pt$lambda)
}

## Streaming accumulator for the co-expression bootstrap SE. Holding all B
## draws at once, a B x pairs matrix for each of total/cis/trans, scales
## memory as O(B x pairs), and pairs scales as p^2, so a larger gene set
## can blow past available memory well before the CPU work finishes.
## Standard deviation only needs sum(x) and sum(x^2), not every individual
## draw, so accumulating those per chunk keeps peak memory at O(chunk_size
## x pairs) instead, and each chunk's raw draws can be discarded once
## folded in. p: number of genes (nrow of a RESID matrix), used to size
## the upper-tri index once.
coexpr_acc_init <- function(p) {
  up <- which(upper.tri(matrix(0, p, p)))
  z  <- numeric(length(up))
  list(up = up, n = 0L,
       sum_total = z, sumsq_total = z,
       sum_cis   = z, sumsq_cis   = z,
       sum_trans = z, sumsq_trans = z,
       sum_dpar_sc = z, sumsq_dpar_sc = z,
       sum_dpar_se = z, sumsq_dpar_se = z)
}

## Folds one chunk's worth of coexpr_bootstrap_one() results into the
## accumulator; the caller can discard chunk_results immediately after
coexpr_acc_update <- function(acc, chunk_results) {
  m_total   <- do.call(rbind, lapply(chunk_results, `[[`, "total"))
  m_cis     <- do.call(rbind, lapply(chunk_results, `[[`, "cis"))
  m_trans   <- do.call(rbind, lapply(chunk_results, `[[`, "trans"))
  m_dpar_sc <- do.call(rbind, lapply(chunk_results, `[[`, "dpar_sc"))
  m_dpar_se <- do.call(rbind, lapply(chunk_results, `[[`, "dpar_se"))
  acc$n <- acc$n + nrow(m_total)
  acc$sum_total   <- acc$sum_total   + colSums(m_total)
  acc$sumsq_total <- acc$sumsq_total + colSums(m_total^2)
  acc$sum_cis     <- acc$sum_cis     + colSums(m_cis)
  acc$sumsq_cis   <- acc$sumsq_cis   + colSums(m_cis^2)
  acc$sum_trans   <- acc$sum_trans   + colSums(m_trans)
  acc$sumsq_trans <- acc$sumsq_trans + colSums(m_trans^2)
  acc$sum_dpar_sc   <- acc$sum_dpar_sc   + colSums(m_dpar_sc)
  acc$sumsq_dpar_sc <- acc$sumsq_dpar_sc + colSums(m_dpar_sc^2)
  acc$sum_dpar_se   <- acc$sum_dpar_se   + colSums(m_dpar_se)
  acc$sumsq_dpar_se <- acc$sumsq_dpar_se + colSums(m_dpar_se^2)
  acc
}

## Finishes the accumulator into the same CB structure
## assemble_coexpr_bootstrap() produces (total/cis/trans/dpar_sc/dpar_se
## data.frames with gene_i/gene_j/est/se/z/p, plus lambda), computed from
## running sums rather than a full B x pairs matrix
coexpr_acc_finalize <- function(resid, pt, acc) {
  p  <- nrow(resid$MIX.SC); gn <- rownames(resid$MIX.SC)
  ij <- arrayInd(acc$up, c(p, p))
  base <- data.frame(gene_i = gn[ij[, 1]], gene_j = gn[ij[, 2]])

  se_of <- function(s, ss) sqrt(pmax(0, (ss - s^2 / acc$n) / (acc$n - 1)))
  mk <- function(est, s, ss) { se <- se_of(s, ss)
    data.frame(base, est = est, se = se, z = est / se, p = 2 * pnorm(-abs(est / se))) }

  list(total   = mk(pt$total[acc$up],   acc$sum_total,   acc$sumsq_total),
       cis     = mk(pt$cis[acc$up],     acc$sum_cis,     acc$sumsq_cis),
       trans   = mk(pt$trans[acc$up],   acc$sum_trans,   acc$sumsq_trans),
       dpar_sc = mk(pt$dpar_sc[acc$up], acc$sum_dpar_sc, acc$sumsq_dpar_sc),
       dpar_se = mk(pt$dpar_se[acc$up], acc$sum_dpar_se, acc$sumsq_dpar_se),
       lambda = pt$lambda)
}

## Serial, local convenience wrapper built from the same primitives the
## cluster job uses. Fine for a quick small-B check; real runs at
## N.COEXPR go through coexpr_boot.R instead, since each draw redoes a
## full shrink_cor per dataset and that cost adds up fast at B in the
## thousands.
coexpr_bootstrap <- function(resid, B = 200, seed = 1) {
  pt    <- coexpr_decompose(resid)
  draws <- make_coexpr_draws(ncol(resid$MIX.SC), ncol(resid$MIX.SE), ncol(resid$HYB.SC), B, seed)
  draw_results <- lapply(draws, coexpr_bootstrap_one, resid = resid)
  assemble_coexpr_bootstrap(resid, pt, draw_results)
}

## Per-gene reliability, the fraction of a gene's total variance that is
## biological rather than Poisson sampling noise. rho = mu / (mu + k), from
## the MU/DISP already fit for every gene, no new model needed. A pairwise
## correlation attenuates relative to its true value by approximately
## sqrt(rho_i * rho_j), so this is the quantity that should predict how
## noisy a pair's residual correlation is, and the natural basis for a
## principled inclusion floor on the co-expression gene set.
## fits: named list of per-dataset FITS.OBS data.frames (MU, DISP columns,
## rownames = gene). Reliability is the minimum across datasets, since a
## gene used in every one of them is only as reliable as its worst context.
gene_reliability <- function(fits, genes) {
  rho_mat <- sapply(fits, function(fr) {
    mu <- fr[genes, "MU"]; k <- fr[genes, "DISP"]
    mu / (mu + k)     # DISP = Inf (Poisson limit) correctly gives rho = 0
  })
  setNames(apply(rho_mat, 1, min, na.rm = TRUE), genes)
}

## Tests whether observed bootstrap SE for co-expression pairs actually
## tracks the NB-predicted attenuation factor sqrt(rho_i * rho_j), before
## committing to a rho floor for gene set inclusion.
## pairs: one of CB$total/cis/trans (data.frame with gene_i, gene_j, se).
## rho: named vector from gene_reliability(), covering all genes in pairs.
## floors: candidate rho cutoffs to report retained pair count and median
## SE for, so the floor can be read off where SE actually stabilizes
## rather than assumed from the attenuation formula alone.
check_coexpr_reliability <- function(pairs, rho, n_bins = 10, floors = seq(0.1, 0.8, by = 0.1)) {
  attn <- sqrt(rho[pairs$gene_i] * rho[pairs$gene_j])
  ok   <- is.finite(attn) & is.finite(pairs$se)
  attn <- attn[ok]; se <- pairs$se[ok]

  bins  <- cut(attn, breaks = quantile(attn, seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  bin_x <- tapply(attn, bins, median)
  bin_y <- tapply(se,   bins, median)
  bin_n <- tapply(se,   bins, length)

  floor_summary <- do.call(rbind, lapply(floors, function(f) {
    keep <- attn >= f
    data.frame(floor = f, n_pairs = sum(keep),
               median_se = if (any(keep)) median(se[keep]) else NA_real_)
  }))

 plot(attn, se, pch = 19, cex = 0.4, col = "grey70", xlab = expression(sqrt(rho[i] * rho[j])), ylab = "bootstrap SE", main = "Co-expression SE vs predicted reliability")
  lines(bin_x, bin_y, type = "b", pch = 19, lwd = 2)

  list(spearman = cor(attn, se, method = "spearman"),
       bins = data.frame(x = bin_x, y = bin_y, n = bin_n),
       floor_summary = floor_summary)
}

## ============================================================
## Intrinsic / extrinsic noise: reliability calibration
## ============================================================
## Before picking a reliability floor for the allele-pair correlation
## analysis (a gene's own within-hybrid-cell correlation between its Sc-
## allele and Se-allele NB Pearson residuals), this checks how the
## correlation's bootstrap SE actually behaves as a function of predicted
## attenuation, the same diagnostic idea as check_coexpr_reliability()
## above. Two differences from that check matter here. First, this is
## one correlation per gene, not a p x p matrix, so there is no
## contamination-of-neighbors argument for pre-filtering the gene set;
## every gene can be checked on its own. Second, the intended use here is
## to disattenuate the correlation (divide by sqrt(rho_Sc*rho_Se)), which
## co-expression never does, so a second, gene-specific failure mode
## exists that the co-expression check never had to look for: at low
## reliability, that division inflates estimation noise on top of the
## attenuation itself. Both plots are produced so a floor can be read off
## each independently rather than assumed to be the same place.

## Row-wise Pearson correlation between two same-shape matrices (here,
## genes x cells), vectorized across rows rather than looped per gene
row_cor <- function(A, B) {
  am <- rowMeans(A); bm <- rowMeans(B)
  Ac <- A - am; Bc <- B - bm
  den <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  ifelse(den > 0, rowSums(Ac * Bc) / den, NA_real_)
}

## A reliability-stratified sample of genes for the calibration check, so
## the low end of the attenuation range, rare among genes overall since
## most fits cluster at moderate-to-high reliability, is still
## represented well enough to see whether the SE actually flattens out
## down there or the check just runs out of genes to show it
stratified_rho_sample <- function(attn, n_per_bin = 40, n_bins = 10, seed = 1) {
  set.seed(seed)
  ok    <- is.finite(attn)
  bins  <- cut(attn[ok], breaks = quantile(attn[ok], seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  genes <- names(attn)[ok]
  unname(unlist(tapply(genes, bins, function(g) sample(g, min(n_per_bin, length(g))))))
}

## Bootstrap SE of each sampled gene's own allele-residual correlation,
## from B hybrid-cell resamples per gene. Looped per gene rather than
## vectorized across genes, since this only runs on the few hundred genes
## in the calibration sample, not genome-wide; a genome-wide production
## version would move this to the cluster the way every other per-gene
## bootstrap in this pipeline does.
allele_cor_boot_se <- function(sc, se, B = 300, seed = 1) {
  set.seed(seed)
  n <- ncol(sc)
  vapply(seq_len(nrow(sc)), function(i) {
    a <- sc[i, ]; b <- se[i, ]
    r <- vapply(seq_len(B), function(k) {
      idx <- sample.int(n, n, replace = TRUE)
      x <- a[idx] - mean(a[idx]); y <- b[idx] - mean(b[idx])
      den <- sqrt(sum(x^2) * sum(y^2))
      if (den > 0) sum(x * y) / den else NA_real_
    }, numeric(1))
    sd(r, na.rm = TRUE)
  }, numeric(1))
}

## sc, se: HYB.SC / HYB.SE residual matrices, genes x cells, same cell
## columns, built genome-wide rather than restricted to CO.GENES, since
## nothing here needs the restriction co-expression's shared correlation
## matrix required.
## rho_sc, rho_se: gene_reliability() applied to the HYB.SC and HYB.SE
## fits separately (not RHO.ALL's cross-dataset minimum, since the two
## alleles' reliabilities enter the attenuation factor sqrt(rho_Sc*rho_Se)
## separately, not as a single combined number).
check_intrinsic_reliability <- function(sc, se, rho_sc, rho_se, n_per_bin = 40, n_bins = 10, B = 300, seed = 1, floors = seq(0.05, 0.6, by = 0.05)) {
  genes <- intersect(rownames(sc), rownames(se))
  attn  <- sqrt(rho_sc[genes] * rho_se[genes])
  samp  <- stratified_rho_sample(attn, n_per_bin, n_bins, seed)

  rho_obs <- row_cor(sc[samp, , drop = FALSE], se[samp, , drop = FALSE])
  se_obs  <- allele_cor_boot_se(sc[samp, , drop = FALSE], se[samp, , drop = FALSE], B, seed)
  a       <- attn[samp]
  rho_true <- rho_obs / a
  ## Exact error propagation for a fixed-scale division: SE(x/a) = SE(x)/a
  se_true  <- se_obs / a

  bins  <- cut(a, breaks = quantile(a, seq(0, 1, length.out = n_bins + 1)), include.lowest = TRUE)
  bin_x <- tapply(a, bins, median)

  floor_summary <- do.call(rbind, lapply(floors, function(f) {
    keep <- a >= f
    data.frame(floor = f, n_genes = sum(keep),
               median_se_obs  = if (any(keep)) median(se_obs[keep])  else NA_real_,
               median_se_true = if (any(keep)) median(se_true[keep]) else NA_real_)
  }))

  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  plot(a, se_obs, pch = 19, cex = 0.4, col = "grey70", xlab = expression(sqrt(rho[Sc] * rho[Se])), ylab = "bootstrap SE, raw correlation", main = "Raw allele correlation")
  lines(bin_x, tapply(se_obs, bins, median), type = "b", pch = 19, lwd = 2)

 plot(a, se_true, pch = 19, cex = 0.4, col = "grey70", xlab = expression(sqrt(rho[Sc] * rho[Se])), ylab = "bootstrap SE, disattenuated", main = "Disattenuated allele correlation")
  lines(bin_x, tapply(se_true, bins, median), type = "b", pch = 19, lwd = 2)

  list(genes = samp, rho_obs = rho_obs, se_obs = se_obs,
       rho_true = rho_true, se_true = se_true, attn = a,
       floor_summary = floor_summary)
}

## Correlation between two residual vectors with a third variable (per-cell
## depth) partialled out. Used to test whether an allele-pair correlation
## reflects shared biological noise or leftover depth structure the NB
## offset failed to remove; if depth is the driver, this collapses toward
## zero relative to the raw correlation, and if not, it tracks the raw
## correlation closely.
partial_cor_depth <- function(x, y, z) {
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

## Mean and variance of a continuous score within each level of a
## classification vector, both indexed identically. Called on more than
## one gene subset (all genes, clipped genes excluded) and more than one
## classification (regulatory, dominance, burst kinetics), so kept as a
## function rather than repeated inline for each combination.
class_mean_var <- function(x, class) {
  ok  <- !is.na(class) & !is.na(x)
  agg <- aggregate(x[ok], by = list(class = class[ok]), FUN = function(v) c(n = length(v), mean = mean(v), var = var(v)))
  data.frame(class = agg$class, n = agg$x[, "n"], mean = agg$x[, "mean"],
             var = agg$x[, "var"], row.names = NULL)
}

## Splits genes into low, average, and high groups by a continuous score,
## using the outer quantiles (probs) as cutoffs. Returns three gene
## vectors, one per group, for enrichment against the shared background.
frac_group_sets <- function(genes, score, probs = c(0.25, 0.75)) {
  cuts <- quantile(score, probs, na.rm = TRUE)
  list(Low     = genes[score <= cuts[1]],
       Average = genes[score >  cuts[1] & score < cuts[2]],
       High    = genes[score >= cuts[2]])
}

## One-way ANOVA testing whether a continuous score differs across levels
## of a classification vector, with Tukey's HSD pairwise comparisons run
## only when the omnibus F-test clears sig, so a non-significant overall
## difference isn't followed by pairwise comparisons that have nothing to
## explain. eta_sq is the between-class sum of squares as a fraction of
## total sum of squares, the size of the class effect independent of
## sample size, meant to sit alongside a significant F/p rather than
## stand in for the Tukey gaps themselves. tukey is NULL when the
## omnibus test doesn't clear sig.
class_anova <- function(x, class, sig = 0.05) {
  ok  <- !is.na(class) & !is.na(x)
  df  <- data.frame(x = x[ok], class = factor(class[ok]))
  fit <- aov(x ~ class, data = df)
  ov  <- summary(fit)[[1]]
  f_stat  <- ov["class", "F value"]; p_val <- ov["class", "Pr(>F)"]
  ss_between <- ov["class", "Sum Sq"]; ss_total <- ss_between + ov["Residuals", "Sum Sq"]
  eta_sq  <- ss_between / ss_total
  tukey   <- if (is.finite(p_val) && p_val < sig) TukeyHSD(fit)$class else NULL
  list(f = f_stat, df1 = ov["class", "Df"], df2 = ov["Residuals", "Df"], p = p_val, eta_sq = eta_sq, tukey = tukey)
}

## total/cis/trans is BH-adjusted across all tested pairs before
## classifying, since a quarter of a million simultaneous single-pair
## tests is exactly the situation BH correction exists for. This is the
## direct way to account for the divergence rate expected from noise
## alone, in place of the rho-based correction that did not hold up
## empirically. CB: coexpr_bootstrap() output (or the cluster-assembled
## equivalent). One row per pair, with the five-class call from cis/trans
## carried alongside so no downstream step reclassifies from raw p.
coexpr_class_table <- function(CB, sig = 0.05) {
  padj_total <- p.adjust(CB$total$p, "BH")
  padj_cis   <- p.adjust(CB$cis$p,   "BH")
  padj_trans <- p.adjust(CB$trans$p, "BH")
  cls <- classify_reg(padj_cis, padj_trans, CB$cis$est, CB$trans$est, sig = sig)$class
  data.frame(gene_i = CB$total$gene_i, gene_j = CB$total$gene_j,
             total_est = CB$total$est, total_padj = padj_total,
             cis_est   = CB$cis$est,   cis_padj   = padj_cis,
             trans_est = CB$trans$est, trans_padj = padj_trans,
             class = cls, stringsAsFactors = FALSE)
}

## Pair-level dominance classification, the co-expression analog of
## classify_dom() at the single-gene level. Compares the hybrid's own
## pairwise correlation (both alleles summed, allele identity ignored,
## dpar_sc/dpar_se from coexpr_decompose()) against each parent's
## pairwise correlation separately: does the hybrid pair's co-expression
## resemble Sc, resemble Se, both (Additive), neither in a consistent
## direction (Over/Underdominant), or match neither parent while the two
## parents also differ from each other. Same two-test, BH-adjusted-first
## pattern as coexpr_class_table(), just fed classify_dom() instead of
## classify_reg().
coexpr_dom_class_table <- function(CB, sig = 0.05) {
  padj_dpar_sc <- p.adjust(CB$dpar_sc$p, "BH")
  padj_dpar_se <- p.adjust(CB$dpar_se$p, "BH")
  cls <- classify_dom(padj_dpar_sc, padj_dpar_se, CB$dpar_sc$est, CB$dpar_se$est, sig = sig)$class
  data.frame(gene_i = CB$dpar_sc$gene_i, gene_j = CB$dpar_sc$gene_j,
             dpar_sc_est = CB$dpar_sc$est, dpar_sc_padj = padj_dpar_sc,
             dpar_se_est = CB$dpar_se$est, dpar_se_padj = padj_dpar_se,
             class = cls, stringsAsFactors = FALSE)
}

## Per-gene degree in the co-expression divergence network, how many
## partners each gene has a BH-significant divergent pair with, split by
## which contrast drives it. degree_any uses the five-class call so it
## matches Figure 11 exactly rather than re-deriving significance. Genes
## with high degree_any are candidates for being actual rewiring hubs
## rather than one member of a single divergent pair.
## class_table: coexpr_class_table() output. genes: CO.GENES, so genes
## with zero hits still appear with degree 0 rather than being dropped.
coexpr_gene_degree <- function(class_table, genes, sig = 0.05) {
  count_for <- function(hit) {
    tab <- table(c(class_table$gene_i[hit], class_table$gene_j[hit]))
    out <- setNames(integer(length(genes)), genes)
    out[names(tab)] <- tab[names(tab)]
    out
  }
  data.frame(gene         = genes,
             degree_total = count_for(class_table$total_padj < sig),
             degree_cis   = count_for(class_table$cis_padj   < sig),
             degree_trans = count_for(class_table$trans_padj < sig),
             degree_any   = count_for(class_table$class != "Conserved"),
             row.names = NULL)
}

## ============================================================
## Fig 11 : co-expression cis vs trans (per gene pair)
## ============================================================
## CT is coexpr_class_table(CB), not CB itself, so the class shown here is
## the same BH-adjusted call used for the pair lists and gene degree, not
## a separate raw-p reclassification
plot_coexpr_cis_trans <- function(CT, main = NULL, lim = NULL) {
  x <- CT$cis_est; y <- CT$trans_est; cls <- CT$class
  ok <- is.finite(x) & is.finite(y) & !is.na(cls)
  x<-x[ok]; y<-y[ok]; cls<-cls[ok]
  col <- COLOR.LIST.1[match(cls,REG.CLASS)]
  if (is.null(lim)) { m<-max(abs(c(x,y)),na.rm=TRUE); lim<-c(-m,m) }
  if (is.null(main)) main<-"co-expression: cis vs trans"
  op<-par(pty="s"); on.exit(par(op))
  plot(NA, xlim=lim, ylim=lim, xlab="cis: hybrid Sc-Se change in correlation", ylab="trans: parents - hybrid (change in correlation)", main=main)
  abline(0,1,lty=3,col="grey65"); abline(0,-1,lty=3,col="grey65")
  abline(h=0,v=0,col="grey45")
  points(x,y,pch=16,cex=0.4,col=col)
  legend("topleft",legend=REG.CLASS,col=COLOR.LIST.1,pch=16,bty="n",cex=0.8)
}

## ============================================================
## Coexpr dominance : dpar_sc vs dpar_se scatter for co-expression pairs
## ============================================================
## DT is coexpr_dom_class_table(CB) output. Same axis pair as the
## single-gene dominance scatter (plot_dom_class, frame = "parent"), just
## for the pairwise correlation change rather than a single burst
## quantity, so the two figures read the same way side by side.
plot_coexpr_dom_class <- function(DT, main = NULL, lim = NULL) {
  x <- DT$dpar_sc_est; y <- DT$dpar_se_est; cls <- DT$class
  ok <- is.finite(x) & is.finite(y) & !is.na(cls)
  x <- x[ok]; y <- y[ok]; cls <- cls[ok]
  col <- COLOR.LIST.2[match(cls, DOM.CLASS)]
  if (is.null(lim)) { m <- max(abs(c(x, y)), na.rm = TRUE); lim <- c(-m, m) }
  if (is.null(main)) main <- "co-expression: dominance"
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = lim, ylim = lim, xlab = "hybrid - Sc parent (correlation change)", ylab = "hybrid - Se parent (correlation change)", main = main)
  abline(0, 1, lty = 3, col = "grey65"); abline(0, -1, lty = 3, col = "grey65")
  abline(h = 0, v = 0, col = "grey45")
  points(x, y, pch = 16, cex = 0.4, col = col)
  legend("topleft", legend = DOM.CLASS, col = COLOR.LIST.2, pch = 16, bty = "n", cex = 0.8)
}

## Residual scatter for one gene pair in each parent dataset, Sc in red
## and Se in blue, with the per-dataset correlation in the legend.
## resid: RESID (genes x cells per dataset). gi, gj: gene names.
plot_coexpr_pair <- function(resid, gi, gj, main = NULL) {
  xsc <- resid$MIX.SC[gi, ]; ysc <- resid$MIX.SC[gj, ]
  xse <- resid$MIX.SE[gi, ]; yse <- resid$MIX.SE[gj, ]
  lim_x <- range(c(xsc, xse), finite = TRUE)
  lim_y <- range(c(ysc, yse), finite = TRUE)
  if (is.null(main)) main <- sprintf("%s vs %s", gi, gj)
  plot(xsc, ysc, pch = 16, cex = 0.5, col = adjustcolor("firebrick", 0.5), xlim = lim_x, ylim = lim_y, xlab = gi, ylab = gj, main = main)
  points(xse, yse, pch = 16, cex = 0.5, col = adjustcolor("steelblue", 0.5))
  legend("topleft", bty = "n", pch = 16, col = c("firebrick", "steelblue"), legend = c(sprintf("Sc  r = %.2f", cor(xsc, ysc, use = "complete.obs")), sprintf("Se  r = %.2f", cor(xse, yse, use = "complete.obs"))))
}

## Two-seed adequacy check for the co-expression bootstrap, same idea as
## the per-gene two-seed comparison. If SE hasn't converged at this B, a
## large fraction of pairs can look spuriously significant simply because
## their SE was underestimated by chance in one particular run, which
## would show up as a widespread, near-uniform high degree rather than a
## sparse set of real hub genes, exactly the pattern under question here.
## cb1, cb2: two coexpr_bootstrap() / cluster outputs at the same B, gene
## set, and everything else, differing only in seed. Rows are assumed to
## be in the same gene_i/gene_j order, true whenever both came from the
## same CO.GENES and upper.tri() call.
## Verifies a loaded coexpr_boot1/2_output.rda-style CB object has the pair
## count CO.GENES currently implies, catching a stale or mismatched file
## before it propagates into CB.CLASS, the pair lists, or a seed
## comparison. One check, applied identically to every CB object loaded,
## rather than one written separately per file.
check_coexpr_pairs <- function(cb, expected_pairs, label) {
  if (nrow(cb$total) != expected_pairs)
  stop(sprintf("%s has %d pairs but CO.GENES expects %d; rerun the matching cluster job", label, nrow(cb$total), expected_pairs))
}

## Shared core for the two-seed bootstrap SE adequacy checks below.
## Given two aligned SE vectors from independent bootstrap seeds,
## reports whether the SE has converged: a scatter against the y = x
## line, the seed-to-seed Spearman correlation, and the median and IQR
## of the SE ratio. Used by both gene_seed_compare() (per-gene) and
## coexpr_seed_compare() (per-pair), which differ only in how they
## extract and align se1/se2 before calling this.
seed_compare_core <- function(se1, se2, main, count_label) {
  ok  <- is.finite(se1) & is.finite(se2)
  se1 <- se1[ok]; se2 <- se2[ok]
  ratio <- se1 / se2

  plot(se1, se2, pch = 16, cex = 0.4, col = adjustcolor("black", 0.3), xlab = "SE, seed 1", ylab = "SE, seed 2", main = main)
  abline(0, 1, col = "firebrick", lty = 2)

  out <- list(cor = cor(se1, se2, method = "spearman"),
              ratio_median = median(ratio), ratio_iqr = IQR(ratio))
  out[[count_label]] <- length(se1)
  out
}

coexpr_seed_compare <- function(cb1, cb2, mode = c("total", "cis", "trans", "dpar_sc", "dpar_se")) {
  mode <- match.arg(mode)
  seed_compare_core(cb1[[mode]]$se, cb2[[mode]]$se, main = sprintf("%s: bootstrap SE, two seeds", mode), count_label = "n_pairs")
}

## Two-seed adequacy check for the per-gene bootstrap (BURST.CONTRASTS), the check
## coexpr_seed_compare() above was itself modeled on. If N.BOOT hasn't
## converged, a gene's SE estimate is unstable from one bootstrap run to
## the next, which shows up here as a low seed-to-seed correlation or a
## wide SE-ratio spread, a failure mode a single bootstrap run alone
## can't reveal since it only ever sees one draw of that instability.
## bc1, bc2: two add_burst_contrasts()-style BURST.CONTRASTS data frames at the same
## N.BOOT, gene set, and everything else, differing only in SEED.BOOT.
## Matched by gene rather than assumed to share row order.
gene_seed_compare <- function(bc1, bc2, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), mode = c("total", "cis", "trans", "dom", "dpar_sc", "dpar_se", "inh_sc", "inh_se")) {
  quantity <- match.arg(quantity); mode <- match.arg(mode)
  col <- paste0(quantity, "_", mode, "_se")
  m   <- match(bc1$gene, bc2$gene)
  seed_compare_core(bc1[[col]], bc2[[col]][m], main = sprintf("%s %s: bootstrap SE, two seeds", quantity, mode), count_label = "n_genes")
}

## Fits a rank-k eigendecomposition to a symmetric gene-by-gene
## divergence matrix (M ~ sum_k lambda_k v_k v_k^T) and reports how well
## that low-rank reconstruction predicts the observed off-diagonal
## entries. total_mat: full p x p point-estimate matrix, e.g.
## COEXPR.POINT$total. k: number of leading components to use for the
## reconstruction. Returns the reconstruction R^2, the full eigenvalues
## and eigenvectors, and loading1 (the first eigenvector, named by
## gene).
coexpr_rank_check <- function(total_mat, k = 1) {
  diag(total_mat) <- 0
  eig  <- eigen(total_mat, symmetric = TRUE)

  V     <- eig$vectors[, seq_len(k), drop = FALSE]
  recon <- V %*% diag(eig$values[seq_len(k)], k, k) %*% t(V)

  ij   <- which(upper.tri(total_mat), arr.ind = TRUE)
  obs  <- total_mat[ij]
  pred <- recon[ij]
  r2   <- cor(obs, pred)^2

 plot(pred, obs, pch = 16, cex = 0.3, col = adjustcolor("black", 0.3), xlab = sprintf("rank-%d reconstruction", k), ylab = "observed total_est", main = sprintf("rank-%d model R^2 = %.3f", k, r2))
  abline(0, 1, col = "firebrick", lty = 2)

  list(r2 = r2, values = eig$values, vectors = eig$vectors, loading1 = setNames(eig$vectors[, 1], rownames(total_mat)))
}

## Takes the eigendecomposition from coexpr_rank_check() and identifies
## candidate axes: the top n_candidate axes by |eigenvalue|, each with
## its share of total variance (axis_var) and its effective number of
## driving genes (axis_pr, the participation ratio 1/sum(loading^4) for
## a unit-length loading vector). sig_axes keeps axes at or above
## var_floor and at or above the eff_genes_min participation-ratio
## floor. label prefixes the console message when an axis is dropped.
coexpr_candidate_axes <- function(rank_check, n_candidate = 40, var_floor = 0.01, eff_genes_min = 10, label = "") {
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
       table = data.frame(axis = axis_order, var = axis_var, eff_genes = axis_pr)[1:15, ])
}

## Fits a two-component Gaussian mixture to each candidate axis's gene
## loadings, splits genes into two poles by posterior probability
## (> 0.5), and runs BP GO enrichment on each pole against co_genes.
## Writes one page per axis to pdf_path and returns a named list, one
## entry per axis, with variance/participation-ratio stats, mixture
## parameters, the two gene sets, and their enrichment results.
coexpr_axis_mixtures <- function(rank_check, mat, sig_axes, axis_var, axis_pr, co_genes, pdf_path, axis_label = "axis") {
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

    enrich_lo <- if (length(genes_lo) >= 5) simplify(enrichGO(gene = genes_lo, universe = co_genes, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP")) else NULL
    enrich_hi <- if (length(genes_hi) >= 5) simplify(enrichGO(gene = genes_hi, universe = co_genes, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP")) else NULL

    out[[paste0("axis", k)]] <- list(
      var_explained = axis_var[as.character(k)], eff_genes = axis_pr[as.character(k)],
      mu = mix_k$mu[c(lo, hi)], sigma = mix_k$sigma[c(lo, hi)], lambda = mix_k$lambda[c(lo, hi)],
      genes_lo = genes_lo, genes_hi = genes_hi, enrich_lo = enrich_lo, enrich_hi = enrich_hi)
  }
  dev.off()
  out
}

## Tests each of the top n_top candidate axes against its own
## permutation null: the k-th largest candidate axis (by variance
## share) against the null's own k-th largest squared eigenvalue,
## rank-matched rather than compared to a single shared threshold.
## Restricts extra_axes (the output of coexpr_axis_mixtures()) to
## whichever axes clear alpha, and returns that restricted list
## alongside the p-value table and the vector of validated axis numbers.
coexpr_axis_validate <- function(rank_check, axis_var, null_ranks, extra_axes, n_top = 15, alpha = 0.05) {
  total_ss       <- sum(rank_check$values^2)
  candidate_axes <- as.integer(names(axis_var)[1:n_top])
  candidate_raw  <- axis_var[as.character(candidate_axes)] * total_ss
  ## Add-one correction (as in class_identity_overlap()'s kappa_p),
  ## since a leading axis routinely beats every one of the permutation
  ## draws: without it, mean(null >= obs) reports an exact 0 that
  ## overstates precision no finite permutation count can support,
  ## rather than the true floor of 1 / (n_perm + 1)
  candidate_p    <- sapply(seq_along(candidate_raw), function(k) (1 + sum(null_ranks[, k] >= candidate_raw[k])) / (1 + nrow(null_ranks)))

  validated_axes <- candidate_axes[candidate_p < alpha]

  list(table = data.frame(axis = candidate_axes, raw = candidate_raw, p_value = candidate_p),
       validated_axes = validated_axes,
       extra_axes = extra_axes[paste0("axis", validated_axes)])
}

## ============================================================
## 8. REGULATORY AND DOMINANCE CLASSIFICATION  (offset level)
## ============================================================
## Per-gene calls read significance from PERM.RESULTS and effect
## direction from BURST.CONTRASTS. One offset level, no NORM index. Classifiers
## return class as a plain character vector so downstream relabeling
## (clean_reg) can reassign freely; ordering is applied at the call
## site with factor(x, levels = REG.CLASS / DOM.CLASS).

## --- regulatory cis/trans, two-test scheme ----------------------------
## Significance from the cis and trans permutation p-values only. The
## parental (total) test is not used, so genes significant in both with
## opposing cis and trans directions cannot be split into true compensatory
## versus cis-by-trans and are reported together as Compensatory. Same-sign
## genes are Cis + Trans. REG.CLASS is defined in the main script.

classify_reg <- function(p_cis, p_trans, est_cis, est_trans, sig = 0.05, colors = COLOR.LIST.1) {
  n   <- length(p_cis)
  cls <- character(n)
  sgn <- sign(est_cis * est_trans)
  cs  <- p_cis   < sig
  ts  <- p_trans < sig
  cls[!cs & !ts]            <- "Conserved"
  cls[ cs & !ts]            <- "Cis"
  cls[!cs &  ts]            <- "Trans"
  cls[ cs &  ts & sgn >= 0] <- "Cis + Trans"
  cls[ cs &  ts & sgn <  0] <- "Compensatory"
  ## Genes where the permutation p-value is significant but the bootstrap
  ## point estimate is NA (sgn undefined) get no label from the five
  ## rules above and are flagged Ambiguous rather than left blank
  cls[cls == ""] <- "Ambiguous"
  data.frame(class = cls,
             color = colors[match(cls, REG.CLASS)],
             stringsAsFactors = FALSE)
}

###Dominance classification (hybrid vs each parent)
## Parents are haploid, the hybrid diploid. The per-sample library offset
## puts every MU on a relative-abundance scale, so the genome-wide ploidy
## factor cancels in the dpar and dom contrasts and the calls below are on
## a common scale. The per-allele inh contrasts measure one allele against
## the full diploid library and so carry a ~1 log2 baseline, which is why
## inh is descriptive and is not used here. DOM.CLASS is defined in the
## main script.

classify_dom <- function(p_dpar_sc, p_dpar_se, est_dpar_sc, est_dpar_se, sig = 0.05, colors = COLOR.LIST.2) {
  n   <- length(p_dpar_sc)
  cls <- character(n)
  sc  <- p_dpar_sc < sig                              # differs from Sc parent
  se  <- p_dpar_se < sig                              # differs from Se parent
  tsg <- sign(sign(est_dpar_sc) + sign(est_dpar_se))
  cls[!sc & !se]            <- "Conserved"
  cls[ sc & !se]            <- "Se.Dominant"
  cls[!sc &  se]            <- "Sc.Dominant"
  cls[ sc &  se & tsg >  0] <- "Overdominant"
  cls[ sc &  se & tsg <  0] <- "Underdominant"
  cls[ sc &  se & tsg == 0] <- "Additive"
  ## Genes where tsg is undefined (bootstrap estimate NA despite a
  ## significant permutation p-value) get no label above and are
  ## flagged Ambiguous rather than left blank
  cls[cls == ""] <- "Ambiguous"
  data.frame(class = cls,
             color = colors[match(cls, DOM.CLASS)],
             stringsAsFactors = FALSE)
}

## --- label cleaner for the relationship panels ------------------------
#Collapse "Cis x Trans" into "Compensatory" and drop "Ambiguous" in a class vector
clean_reg <- function(v) {
  v[v == "Cis x Trans"] <- "Compensatory"
  v[v == "Ambiguous"]   <- NA
  v
}

## ---- per-gene class vectors aligned to BURST.CONTRASTS rows -------------------------
reg_class_vec <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05) {
  quantity <- match.arg(quantity)
  i <- match(BURST.CONTRASTS$gene, PR$gene)
  classify_reg(PR[[paste0(quantity, "_cis_p")]][i], PR[[paste0(quantity, "_trans_p")]][i], BURST.CONTRASTS[[paste0(quantity, "_cis_est")]], BURST.CONTRASTS[[paste0(quantity, "_trans_est")]], sig = sig)$class
}

dom_class_vec <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05) {
  quantity <- match.arg(quantity)
  i <- match(BURST.CONTRASTS$gene, PR$gene)
  classify_dom(PR[[paste0(quantity, "_dpar_sc_p")]][i], PR[[paste0(quantity, "_dpar_se_p")]][i], BURST.CONTRASTS[[paste0(quantity, "_dpar_sc_est")]], BURST.CONTRASTS[[paste0(quantity, "_dpar_se_est")]], sig = sig)$class
}

## --- class association heatmap ----------------------------------------
## Shows the association between two class vectors as a grid of log2
## observed-over-expected counts, pale where a per-cell two-sided
## hypergeometric test is not significant after BH. brk sets colour breaks.
## Vertical color-scale strip just outside the right edge of the current
## plot region, mapping the heatmap's colors to fold-enrichment values.
## Called right after image(), while par("usr") still reflects the
## heatmap's own 0-1 coordinate system.
.heatmap_legend <- function(pal, rng, brk) {
  usr <- par("usr")
  lx <- usr[2] + diff(usr[1:2]) * 0.12
  rx <- usr[2] + diff(usr[1:2]) * 0.24
  ys <- seq(usr[3], usr[4], length.out = brk + 1)
  rect(lx, ys[-length(ys)], rx, ys[-1], col = pal, border = NA, xpd = NA)
  rect(lx, usr[3], rx, usr[4], border = "grey30", xpd = NA)
  text(rx, usr[4],         sprintf("%.2gx", 2^rng),  pos = 4, cex = 0.55, xpd = NA)
  text(rx, mean(usr[3:4]), "1x",                      pos = 4, cex = 0.55, xpd = NA)
  text(rx, usr[3],         sprintf("%.2gx", 2^-rng),  pos = 4, cex = 0.55, xpd = NA)
  text(rx + diff(usr[1:2])*0.16, mean(usr[3:4]), "fold enrichment (obs/exp)", srt = 270, cex = 0.55, xpd = NA)
}

## Shared core for the two class-overlap heatmaps below: a categorical
## contingency table between two classifications, tested cell by cell
## against the hypergeometric null and colored by log2(observed/expected).
## The star in each cell marks BH-corrected significance at fdr; the
## color always encodes the fold enrichment itself, significant or not.
## levels_a/levels_b: if NULL, the table's levels come from whichever
## categories are actually observed (class_heatmap's original
## behavior). If supplied, every level appears in the table even when
## empty, so a fixed axis ordering is preserved across genes
## (plot_class_overlap's original behavior, used for the publication
## Figures 2 and 6). cex_cell, fmt, and star_frac each pick their
## original per-caller default automatically when left NULL, so
## neither existing call site's appearance changes.
## rng, when supplied, fixes the color scale's half-range in log2(obs/exp)
## units instead of letting each call size its own scale to its own
## table. Pass the same rng to a set of related panels (e.g. the three
## regulatory overlap panels of Figure 2) so a given color always means
## the same fold enrichment across all of them.
class_overlap_heatmap <- function(class_a, class_b, levels_a = NULL, levels_b = NULL, brk = length(cols), fdr = 0.05, cols = COLOR.LIST.3, cex_axis = 0.75, cex_cell = NULL, fmt = NULL, star_frac = NULL, xlab = NULL, ylab = NULL, rng = NULL) {
  explicit_levels <- !is.null(levels_a) || !is.null(levels_b)
  if (is.null(cex_cell))  cex_cell  <- if (explicit_levels) 0.7   else 0.65
  if (is.null(fmt))       fmt       <- if (explicit_levels) "%.1f" else "%.2g"
  if (is.null(star_frac)) star_frac <- if (explicit_levels) 0.6   else 0.55

  if (brk %% 2 == 0) brk <- brk + 1   # odd brk puts a true white sample at the center
  keep <- !is.na(class_a) & !is.na(class_b)
  a <- if (is.null(levels_a)) factor(class_a[keep]) else factor(class_a[keep], levels = levels_a)
  b <- if (is.null(levels_b)) factor(class_b[keep]) else factor(class_b[keep], levels = levels_b)
  tab <- table(a, b); n <- sum(tab)
  exp <- outer(rowSums(tab), colSums(tab)) / n
  lor <- log2((tab + 0.5) / (exp + 0.5))
  pv  <- matrix(NA_real_, nrow(tab), ncol(tab))
  for (i in seq_len(nrow(tab))) for (j in seq_len(ncol(tab))) {
    q <- tab[i, j]; m <- rowSums(tab)[i]; k <- colSums(tab)[j]
    po <- phyper(q - 1, m, n - m, k, lower.tail = FALSE)
    pu <- phyper(q,     m, n - m, k, lower.tail = TRUE)
    pv[i, j] <- 2 * min(po, pu, 0.5)
  }
  padj <- matrix(p.adjust(pv, "BH"), nrow(tab))
  ## Color encodes log2(obs/exp) for every cell, not just significant
  ## ones; the * in the cell text is the only significance indicator.
  ## rng falls back to this panel's own max if none was supplied, the
  ## original single-panel behavior.
  rng  <- if (is.null(rng)) max(abs(lor), 1e-6) else rng
  breaks <- seq(-rng, rng, length.out = brk + 1)
  pal  <- colorRampPalette(cols)(brk)

  lab_x <- colnames(tab); lab_y <- rev(rownames(tab))
  bottom_in <- max(strwidth(lab_x, units = "inches", cex = cex_axis)) + 0.5
  left_in   <- max(strwidth(lab_y, units = "inches", cex = cex_axis)) + 0.5
  op <- par(mai = c(bottom_in, left_in, 0.3, 1.0)); on.exit(par(op))

  image(t(lor[nrow(lor):1, , drop = FALSE]), col = pal, breaks = breaks, axes = FALSE)
  xs <- seq(0, 1, length.out = ncol(tab)); ys <- seq(0, 1, length.out = nrow(tab))
  axis(1, at = xs, labels = lab_x, las = 2, cex.axis = cex_axis)
  axis(2, at = ys, labels = lab_y, las = 2, cex.axis = cex_axis)
  box()

  for (i in seq_len(nrow(tab))) for (j in seq_len(ncol(tab))) {
    star <- if (padj[i, j] < fdr) "*" else ""
    lbl  <- paste0(sprintf(fmt, 2^lor[i, j]), star)
    txt_col <- if (abs(lor[i, j]) > rng * star_frac) "white" else "black"
    text(xs[j], ys[nrow(tab) - i + 1], lbl, col = txt_col, cex = cex_cell)
  }

  if (!is.null(xlab)) mtext(xlab, side = 1, line = bottom_in/par("csi") - 1)
  if (!is.null(ylab)) mtext(ylab, side = 2, line = left_in/par("csi") - 1)
  .heatmap_legend(pal, rng, brk)
  invisible(list(table = tab, log2_obs_exp = lor, padj = padj))
}

## Wraps class_overlap_heatmap() with no explicit levels: table levels
## come from whichever categories are observed. Used for the Section
## 3.5 / 4.12 diagnostic heatmaps, where the two classifications being
## compared don't always share the same level set.
class_heatmap <- function(class_a, class_b, brk = length(COLOR.LIST.3), fdr = 0.05, cols = COLOR.LIST.3, cex_axis = 0.75, cex_cell = 0.65, xlab = NULL, ylab = NULL) {
  class_overlap_heatmap(class_a, class_b, levels_a = NULL, levels_b = NULL, brk = brk, fdr = fdr, cols = cols, cex_axis = cex_axis, cex_cell = cex_cell, xlab = xlab, ylab = ylab)
}

## ============================================================
## Gene-identity overlap between two categorical classifications
## ============================================================
## class_heatmap() answers "is any one cell of the mean-class x
## noise-class table over- or under-represented" -- eleven separate
## per-cell hypergeometric tests. The claim that heatmap is usually used
## to support in text ("mean-divergence and noise-divergence
## classification aren't just two labels attached to the same underlying
## gene identity") needs one overall number instead: Cohen's kappa is the
## standard chance-corrected agreement statistic for two categorical
## classifications of the same items, with kappa = 0 meaning the two
## classifications agree no more than their own class-size distributions
## would produce by chance, and kappa = 1 meaning perfect agreement.
##
## po is the observed fraction of genes assigned the same class label by
## both classifications (the table's diagonal). pe is the fraction that
## would agree by chance alone, given each classification's own marginal
## class-size distribution (rowSums(tab)/n and colSums(tab)/n), which is
## exactly the expectation used by class_heatmap's per-cell test, summed
## over the diagonal instead of tested cell by cell.
.cohen_kappa <- function(tab) {
  n  <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  (po - pe) / (1 - pe)
}

## class_identity_overlap: kappa plus two complements. Per-class Jaccard
## is computed directly on gene membership (|A_k intersect B_k| /
## |A_k union B_k|), not from the contingency table margins, so a high
## overall kappa driven mostly by one large class (typically Conserved)
## doesn't get misread as evidence that the classes of biological
## interest (Cis, Trans, ...) overlap too -- each level gets its own
## number. The permutation null reshuffles which gene gets which class_b
## label; that is a label permutation, not a resample, so it holds both
## classifications' own class-size distributions fixed, meaning pe is
## identical for every replicate (rowSums come from the unpermuted a,
## colSums are a permutation of the same multiset of b labels) and is
## computed once rather than recomputed nperm times. Integer coding plus
## tabulate() replaces repeated table()/factor() calls inside the loop,
## since this runs nperm times rather than once.
class_identity_overlap <- function(class_a, class_b, levels = REG.CLASS, nperm = 2000, seed = 1) {
  keep <- !is.na(class_a) & !is.na(class_b)
  a <- factor(class_a[keep], levels = levels)
  b <- factor(class_b[keep], levels = levels)
  ## factor() turns any value not among `levels` (e.g. "Ambiguous", a
  ## data-quality flag rather than a real class) into a fresh NA that the
  ## is.na() check above, run before this conversion, can't have caught.
  ## Re-filter now so nothing downstream has to assume upstream
  ## classification only ever hands back the five real levels.
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]; b <- b[ok]
  n <- length(a); K <- length(levels)
  tab <- table(a, b)

  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  kappa_obs <- (po - pe) / (1 - pe)

  jacc <- setNames(numeric(K), levels)
  for (k in levels) {
    inA <- a == k; inB <- b == k
    uni <- sum(inA | inB)
    jacc[k] <- if (uni > 0) sum(inA & inB) / uni else NA_real_
  }

  set.seed(seed)
  ai <- as.integer(a); bi <- as.integer(b)
  kappa_null <- numeric(nperm)
  for (i in seq_len(nperm)) {
    bp   <- bi[sample.int(n, n)]
    tabp <- matrix(tabulate((bp - 1L) * K + ai, K * K), K, K)
    po_p <- sum(diag(tabp)) / n
    kappa_null[i] <- (po_p - pe) / (1 - pe)
  }
  ## One-sided: does observed agreement exceed what reshuffled labels
  ## produce? That is the direction relevant to "not simply reclassified"
  ## -- a small p here would mean the two classifications share gene
  ## identity more than their class sizes alone predict, which is the
  ## outcome the text claim says is NOT the case
  p_kappa <- (1 + sum(kappa_null >= kappa_obs)) / (1 + nperm)

  list(table = tab, n_genes = n,
       concordance = po, expected_concordance = pe,
       kappa = kappa_obs, kappa_null = kappa_null, kappa_p = p_kappa,
       jaccard = jacc)
}

## One printable row per classification-pair comparison: n, concordance,
## chance concordance, kappa, permutation p, and one Jaccard column per
## class level. Bind rows from multiple class_identity_overlap() calls
## (e.g. mean vs noise, mean vs bfreq, mean vs bsize) into one table.
summarize_class_overlap <- function(ov, label = NULL) {
  data.frame(comparison = label, n_genes = ov$n_genes,
            concordance = ov$concordance, expected = ov$expected_concordance,
            kappa = ov$kappa, kappa_p = ov$kappa_p,
            as.list(ov$jaccard), check.names = FALSE, stringsAsFactors = FALSE)
}

## ============================================================
## Figs 1 & 2 : cis vs trans, coloured by regulatory class
## ============================================================
## Maps a vector of SE values to per-point bar colors so that precise
## (small-SE) estimates get a slightly more visible bar and noisy
## (large-SE) estimates fade toward the background, the reverse of what
## a single fixed bar color does when every bar is drawn at the same
## opacity regardless of how long it is. lo/hi are the SE quantiles
## anchoring the low- and high-opacity ends of the ramp.
se_alpha_col <- function(se, base_col = "grey40", lo = 0.05, hi = 0.9, alpha_range = c(0.08, 0.4)) {
  rng <- quantile(se, c(lo, hi), na.rm = TRUE)
  w   <- 1 - pmin(pmax((se - rng[1]) / (rng[2] - rng[1]), 0), 1)   # w = 1 for precise, 0 for noisy
  a   <- alpha_range[1] + diff(alpha_range) * w
  # adjustcolor(alpha.f = <vector>) doesn't do what this needs: current R
  # builds alpha.f into a single 4x4 transform matrix applied to every
  # color at once, so it only accepts one alpha value, not one per point.
  # col2rgb()/rgb() are plain vectorized functions, so building the RGBA
  # color directly here gives every point its own opacity, one base color
  # recycled against a whole vector of alpha values.
  rgb_base <- grDevices::col2rgb(base_col) / 255
  grDevices::rgb(rgb_base[1], rgb_base[2], rgb_base[3], alpha = a)
}

## All SE bars drawn before any points so bars are never in front of
## coloured points. Bar opacity is scaled by se_alpha_col() so the
## longest (least reliable) bars recede rather than dominating the
## panel. Points are solid (filled) rather than hollow, matching the
## Figure 5 dominance panels.
plot_cis_trans_class <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05, main = NULL, lim = NULL, bar_col = "grey40", colors = setNames(COLOR.LIST.1[seq_along(REG.CLASS)], REG.CLASS)) {
  quantity <- match.arg(quantity)
  lab <- c(mean = "mean", bfreq = "burst frequency", bsize = "burst size", kbal = "frequency-size balance", cv2 = "CV2")[quantity]
  cx <- BURST.CONTRASTS[[paste0(quantity,"_cis_est")]];  sx <- BURST.CONTRASTS[[paste0(quantity,"_cis_se")]]
  cy <- BURST.CONTRASTS[[paste0(quantity,"_trans_est")]]; sy <- BURST.CONTRASTS[[paste0(quantity,"_trans_se")]]
  cls <- reg_class_vec(BURST.CONTRASTS, PR, quantity, sig)
  ok  <- is.finite(cx)&is.finite(cy)&is.finite(sx)&is.finite(sy)&cls %in% REG.CLASS
  cx<-cx[ok]; cy<-cy[ok]; sx<-sx[ok]; sy<-sy[ok]; cls<-cls[ok]
  col <- colors[match(cls, REG.CLASS)]
  if (is.null(lim)) { m <- max(abs(c(cx+sx,cx-sx,cy+sy,cy-sy)),na.rm=TRUE); lim<-c(-m,m) }
  if (is.null(main)) main <- paste0(lab,": cis vs trans")
  op <- par(pty="s"); on.exit(par(op))
  plot(NA, xlim=lim, ylim=lim, xlab="log2(Sc/Se) in hybrid", ylab="parents - hybrid (log2)", main=main)
  abline(0,1,lty=3,col="grey65"); abline(0,-1,lty=3,col="grey65")
  abline(h=0,v=0,col="grey45")
  bar_cols_x <- se_alpha_col(sx, bar_col); bar_cols_y <- se_alpha_col(sy, bar_col)
  segments(cx-sx, cy, cx+sx, cy, col=bar_cols_x)   # all bars first, precision-scaled opacity
  segments(cx, cy-sy, cx, cy+sy, col=bar_cols_y)
  points(cx, cy, pch=16, cex=0.5, col=col)   # solid points, matches Figure 5's style
  legend("topleft", legend=REG.CLASS, col=colors, pch=16, bty="n", cex=0.8)
  invisible(cls)
}

## ============================================================
## Fig 4 : mean vs noise per gene, per-class slopes
## ============================================================
## Points drawn in neutral grey. Regulatory slopes are solid lines;
## dominance slopes are dashed. All lines span the full panel (abline).
plot_mean_bfreq_class <- function(BURST.CONTRASTS, PR, mode = "total", reg_class = NULL, dom_class = NULL, y_quantity = c("bfreq", "bsize", "kbal", "cv2"), class_quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05, main = NULL, reg_colors = setNames(COLOR.LIST.1[seq_along(REG.CLASS)], REG.CLASS), dom_colors = setNames(COLOR.LIST.2[seq_along(DOM.CLASS)], DOM.CLASS)) {
  y_quantity     <- match.arg(y_quantity)
  class_quantity <- match.arg(class_quantity)
  y_lab <- c(bfreq = "burst frequency", bsize = "burst size", kbal = "frequency-size balance", cv2 = "CV2")[y_quantity]
  x <- BURST.CONTRASTS[[paste0("mean_",mode,"_est")]]
  y <- BURST.CONTRASTS[[paste0(y_quantity,"_",mode,"_est")]]
  if (is.null(reg_class)) reg_class <- reg_class_vec(BURST.CONTRASTS, PR, class_quantity, sig)
  if (is.null(dom_class)) dom_class <- dom_class_vec(BURST.CONTRASTS, PR, class_quantity, sig)
  ok <- is.finite(x) & is.finite(y)
  x<-x[ok]; y<-y[ok]; rc<-reg_class[ok]; dc<-dom_class[ok]
  if (is.null(main)) main <- paste0("mean vs ", y_lab, ": ", mode)
  px <- diff(range(x))*0.04; py <- diff(range(y))*0.04
  plot(x, y, pch=16, cex=0.4, col="grey85", xlim=range(x)+c(-px, px), ylim=range(y)+c(-py, py), xlab="mean divergence (log2)", ylab=paste0(y_lab, " divergence (log2)"), main=main)
  abline(h=0,v=0,col="grey75")
  for (k in REG.CLASS) {
    sel <- !is.na(rc) & rc==k
    if (sum(sel)>2) abline(lm(y[sel]~x[sel]), col=reg_colors[match(k,REG.CLASS)], lwd=2, lty=1)
  }
  for (k in DOM.CLASS) {
    sel <- !is.na(dc) & dc==k
    if (sum(sel)>2) abline(lm(y[sel]~x[sel]), col=dom_colors[match(k,DOM.CLASS)], lwd=2, lty=2)
  }
  legend("topleft", legend = c(paste0("Reg — ", REG.CLASS), paste0("Dom – ", DOM.CLASS)), col    = c(reg_colors, dom_colors), lty    = c(rep(1, length(REG.CLASS)), rep(2, length(DOM.CLASS))), lwd=2, bty="n", cex=0.65)
}

## ============================================================
## 8b. GO / KEGG enrichment for the regulatory classification
## ============================================================
## build_reg_go_sets: one gene set per (class, direction) pair for a given
## classification quantity ("mean" or "bfreq"). Direction comes from the
## sign of the parental total contrast (mean_total_est or bfreq_total_est)
## at the same sig threshold used for the class call, not from which
## component (cis or trans) drives the class, so "Sc" and "Se" always mean
## the same thing they mean everywhere else in this pipeline: which parent
## is higher in the direct parent-vs-parent comparison. Genes whose total
## contrast is not itself significant fall into neither direction and are
## dropped from both sets (they still count in class totals elsewhere).
## Set names strip spaces/punctuation from the class label for safe list
## and file names, e.g. "Cis + Trans" -> "CisTrans_Sc".
build_reg_go_sets <- function(BURST.CONTRASTS, PR, universe, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), levels = REG.CLASS, sig = 0.05) {
  quantity <- match.arg(quantity)
  cls  <- reg_class_vec(BURST.CONTRASTS, PR, quantity, sig)
  est  <- BURST.CONTRASTS[[paste0(quantity, "_total_est")]]
  i    <- match(BURST.CONTRASTS$gene, PR$gene)
  p0   <- PR[[paste0(quantity, "_total_p")]][i]
  dirn <- ifelse(!is.finite(p0) | p0 >= sig, NA, ifelse(est > 0, "Sc", "Se"))
  genes <- BURST.CONTRASTS$gene

  sets <- list()
  for (k in levels) {
    tag <- gsub("[^A-Za-z]", "", k)
    for (d in c("Sc", "Se")) {
      sel <- !is.na(cls) & cls == k & !is.na(dirn) & dirn == d
      sets[[paste0(tag, "_", d)]] <- intersect(genes[sel], universe)
    }
  }
  sets
}

## build_component_go_sets: "any-cis" or "any-trans" gene sets, pooled
## across regulatory class and direction-matched on the sign of that
## component's OWN estimate, not the total parental contrast
##
## This is deliberately different from build_reg_go_sets. There, direction
## comes from the total contrast regardless of which component drives the
## class, so e.g. "Cis_Sc" means a gene classified Cis whose overall
## parent-vs-parent difference happens to be Sc-higher, which need not
## track the cis estimate's own sign when trans partially opposes it in a
## borderline call. Here, "any-cis" pools every gene with a significant cis
## effect (p_cis < sig) regardless of trans significance -- i.e. the three
## classes where the two-test scheme's cs flag is TRUE: Cis, Cis + Trans,
## and Compensatory -- and splits that pool by whether the cis component
## itself, cis_est, is Sc-higher or Se-higher. "any-trans" mirrors this
## using p_trans and trans_est, pooling Trans, Cis + Trans, and
## Compensatory. Genes are computed straight from p_cis/p_trans and
## cis_est/trans_est rather than by filtering REG.CLASS labels, since that
## is the exact definition of cs/ts in classify_reg and does not depend on
## REG.CLASS staying in sync with classify_reg's internals.
##
## Returns a two-element named list, "AnyCis_Sc"/"AnyCis_Se" or
## "AnyTrans_Sc"/"AnyTrans_Se" depending on `component`
build_component_go_sets <- function(BURST.CONTRASTS, PR, universe, component = c("cis", "trans"), quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05) {
  component <- match.arg(component)
  quantity  <- match.arg(quantity)
  i     <- match(BURST.CONTRASTS$gene, PR$gene)
  p_c   <- PR[[paste0(quantity, "_", component, "_p")]][i]
  est   <- BURST.CONTRASTS[[paste0(quantity, "_", component, "_est")]]
  sigc  <- !is.na(p_c) & p_c < sig
  genes <- BURST.CONTRASTS$gene
  tag   <- if (component == "cis") "AnyCis" else "AnyTrans"

  by_dir <- list(Sc = intersect(genes[sigc & est > 0], universe), Se = intersect(genes[sigc & est < 0], universe))
  setNames(by_dir, paste0(tag, "_", names(by_dir)))
}

## run_enrichment: GO (BP/CC/MF, simplified) and KEGG enrichment for one
## gene set against a fixed universe. Genes not in the universe are dropped
## silently. Sets with fewer than 2 genes after intersection return all-NULL
## rather than being passed to enrichGO/enrichKEGG. kegg_data (from
## kegg_local()) runs KEGG offline through enricher() with the same
## defaults enrichKEGG() uses (BH, p < 0.05, gene-set size 10 to 500). Sets above that floor
## but still small (roughly under 10-15 genes) will run, but any hit there
## should be treated as exploratory, not confirmatory, since GO enrichment
## on sets that small rarely rises above noise; the class-specific splits
## by direction (e.g. Reinforcing_Sc, Compensatory_Se) are the ones most
## likely to fall in this range and n_genes in the summary table flags it.
run_enrichment <- function(genes, universe, orgdb = org.Sc.sgd.db, keytype = "ORF", kegg_org = "sce", qval = 0.2, kegg_data = NULL) {
  genes <- intersect(genes, universe)
  if (length(genes) < 2) return(list(BP = NULL, CC = NULL, MF = NULL, KEGG = NULL))
  go_one <- function(ont) {
    tryCatch(simplify(enrichGO(gene = genes, universe = universe, OrgDb = orgdb, keyType = keytype, ont = ont, qvalueCutoff = qval)), error = function(e) NULL)
  }
  kegg_one <- function() {
    tryCatch(if (is.null(kegg_data)) enrichKEGG(gene = genes, universe = universe, organism = kegg_org, keyType = "kegg", qvalueCutoff = qval)
             else enricher(gene = genes, universe = universe, TERM2GENE = kegg_data$KEGGPATHID2EXTID,
                           TERM2NAME = kegg_data$KEGGPATHID2NAME, qvalueCutoff = qval),
             error = function(e) NULL)
  }
  list(BP = go_one("BP"), CC = go_one("CC"), MF = go_one("MF"), KEGG = kegg_one())
}

## n_sig_terms: number of terms at qvalue < q in one enrichResult, 0 for a
## NULL or empty result (a set too small to test, or no hits)
n_sig_terms <- function(e, q = 0.2) {
  if (is.null(e) || is.null(e@result) || nrow(e@result) == 0) return(0L)
  sum(e@result$qvalue < q, na.rm = TRUE)
}

## print_enrich_brief: console-friendly view of one enrichResult, since
## the default print() dumps every column (including the full comma-
## separated gene list per term) and R's console truncates that to
## fit width, hiding exactly the columns worth reading. Shows only
## Description, p.adjust, and Count, for the n_top most significant
## terms at qvalue < q. Prints one line and returns invisibly for a
## NULL result or one with nothing significant, rather than an empty
## table with no explanation.
print_enrich_brief <- function(e, q = 0.2, n_top = 10) {
  if (is.null(e) || is.null(e@result) || nrow(e@result) == 0) { cat("  (no terms tested)\n"); return(invisible(NULL)) }
  tab <- e@result[e@result$qvalue < q, c("Description", "p.adjust", "Count"), drop = FALSE]
  if (nrow(tab) == 0) { cat(sprintf("  (0 terms at qvalue < %.2f)\n", q)); return(invisible(NULL)) }
  tab <- tab[order(tab$p.adjust), ][seq_len(min(n_top, nrow(tab))), ]
  tab$p.adjust <- signif(tab$p.adjust, 3)
  print(tab, row.names = FALSE)
}

## summarize_go_sets: one row per set, gene set size plus significant-term
## counts for BP/CC/MF/KEGG. `sets` and `enrich` must be the matched
## gene-set list and run_enrichment output (same names, same order).
summarize_go_sets <- function(sets, enrich, q = 0.2) {
  data.frame(
    set     = names(sets),
    n_genes = vapply(sets, length, integer(1)),
    GO_BP   = vapply(enrich, function(e) n_sig_terms(e$BP,   q), integer(1)),
    GO_CC   = vapply(enrich, function(e) n_sig_terms(e$CC,   q), integer(1)),
    GO_MF   = vapply(enrich, function(e) n_sig_terms(e$MF,   q), integer(1)),
    KEGG    = vapply(enrich, function(e) n_sig_terms(e$KEGG, q), integer(1)),
    row.names = NULL
  )
}

## Barplots each ontology (BP, MF, CC, KEGG) from a pair of
## run_enrichment()-style lists (one for an UP gene set, one for DOWN),
## skipping any ontology that is NULL or has no rows (too small to
## test, or nothing significant at the qvalueCutoff run_enrichment()
## used). Each barplot is wrapped in tryCatch, and critically that
## tryCatch wraps both building the ggplot object AND print()ing it, not
## just the build step: ggplot2 defers most computation (scale training,
## stat transforms, the internal ifelse() that enrichplot's barplot()
## method can throw on for a very-few-row enrichResult) until the plot
## is actually drawn, which happens inside print(), not at `+` layer
## construction. A tryCatch around construction alone lets that error
## through unguarded. Without this, the error propagates out of
## whatever pdf() block called this, skipping dev.off() and leaving a
## broken, incomplete file, plus, if the caller loops over several
## clusters, it silently drops every remaining cluster in that same
## pdf(). A failed plot is now noted on the console and skipped instead.
## label should name the comparison (e.g. "Sc major vs minor cluster")
## since "up"/"down" always means ident.1 vs ident.2 from the FindMarkers()
## call that produced enrich_up/enrich_down, and that pairing differs
## across call sites (species direction at some, cluster-vs-rest at
## others), so it cannot be inferred generically inside this function.
barplot_enrich_pair <- function(enrich_up, enrich_down, show = 10, label = NULL) {
  prefix <- if (is.null(label)) "" else paste0(label, " — ")
  safe_plot <- function(e, tag) {
    if (is.null(e) || is.null(e@result) || nrow(e@result) == 0) return(invisible(NULL))
    ok <- tryCatch({
      print(barplot(e, showCategory = show) + ggtitle(paste0(prefix, tag)))
      TRUE
    }, error = function(err) {
      cat(sprintf("  (%s%s: enrichment barplot failed to render — %s; skipped)\n", prefix, tag, conditionMessage(err)))
      FALSE
    })
    invisible(ok)
  }
  for (ont in c("BP", "MF", "CC", "KEGG")) {
    safe_plot(enrich_up[[ont]],   paste0(ont, ", up"))
    safe_plot(enrich_down[[ont]], paste0(ont, ", down"))
  }
}

## Marker/enrichment comparison of each cluster against the rest of the
## SAME dataset's own cells, generalized to however many clusters that
## dataset actually has, not hardcoded to two ("major vs minor"), the
## same generalization already applied to the hybrid cluster-pair
## comparison. obj must carry its own final, validated clustering
## (Idents already set via sweep_cluster_resolution() /
## bootstrap_compare_resolutions(), e.g. YSC.MIX.SC), not a
## merged/combined object's clusters restricted by cell-ID range: the
## question here is whether the clustering actually validated for THIS
## dataset corresponds to distinguishable biology, not whether some
## other object's clusters happen to look meaningful when relabeled.
## ident.2 is left at FindMarkers()'s default (NULL), which compares
## ident.1 against every other cell already inside obj, i.e. "the rest
## of this dataset".
##
## GSEA (gseGO/gseKEGG) is deliberately not computed here: the prior
## version of this code computed it for every comparison but never
## plotted, printed, or otherwise used the result, so it was pure
## unused compute. The over-representation analysis already plotted via
## barplot_enrich_pair() answers the same question (which processes are
## enriched among cluster-specific markers) and is the one actually
## surfaced.
##
## apply_fun sets how clusters are distributed (lapply locally, a forked
## mclapply in cluster_stability.R), and kegg_data passes through to
## run_enrichment().
##
## Returns NULL (with a message) if obj has fewer than two clusters, or
## a list with markers/up/down/up_enrich/down_enrich/background/
## cluster_ids, one entry per cluster, keyed by cluster ID.
cluster_marker_enrichment <- function(obj, label, kegg_data = NULL, apply_fun = lapply) {
  cluster_ids <- sort(unique(as.character(Idents(obj))))
  if (length(cluster_ids) < 2) {
    cat(sprintf("%s: only one cluster found; skipping the per-cluster marker/enrichment comparison.\n", label))
    return(NULL)
  }

  markers_list <- apply_fun(cluster_ids, function(cc) {
    m <- suppressWarnings(FindMarkers(obj, ident.1 = cc))
    m[order(m$avg_log2FC, decreasing = TRUE), ]
  })
  names(markers_list) <- cluster_ids

  up_list   <- lapply(markers_list, function(m) m[abs(m$avg_log2FC) > log2(1.25) & -log10(m$p_val_adj) > 20 & m$avg_log2FC > 0, ])
  down_list <- lapply(markers_list, function(m) m[abs(m$avg_log2FC) > log2(1.25) & -log10(m$p_val_adj) > 20 & m$avg_log2FC < 0, ])

  gene_vec <- function(m) { v <- m[, 2]; names(v) <- row.names(m); v }
  background_list <- lapply(markers_list, gene_vec)
  up_genes_list    <- lapply(up_list,   gene_vec)
  down_genes_list  <- lapply(down_list, gene_vec)

  enrich <- setNames(apply_fun(cluster_ids, function(cc) list(
    up   = run_enrichment(names(up_genes_list[[cc]]),   names(background_list[[cc]]), kegg_data = kegg_data),
    down = run_enrichment(names(down_genes_list[[cc]]), names(background_list[[cc]]), kegg_data = kegg_data))), cluster_ids)
  up_enrich   <- lapply(enrich, `[[`, "up")
  down_enrich <- lapply(enrich, `[[`, "down")

  list(markers = markers_list, up = up_list, down = down_list,
       up_enrich = up_enrich, down_enrich = down_enrich,
       background = background_list, cluster_ids = cluster_ids)
}

## Writes cluster_marker_enrichment()'s up/down enrichment for every
## cluster to one pdf via barplot_enrich_pair(), one page set per
## cluster, each labeled with which cluster is "up" so pages are
## identifiable when the pdf is opened out of context.
plot_cluster_marker_enrichment <- function(res, label, pdf_path, width = 7, height = 5) {
  if (is.null(res)) return(invisible(NULL))
  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  for (cc in res$cluster_ids) {
    barplot_enrich_pair(res$up_enrich[[cc]], res$down_enrich[[cc]],
                         label = sprintf("%s cluster %s vs rest (up = higher in %s)", label, cc, cc))
  }
  dev.off()
}

## Console-text mirror of plot_cluster_marker_enrichment(): for every
## cluster and direction (up/down), prints each ontology's significant
## terms via print_enrich_brief() (Description, p.adjust, Count), so the
## results are readable without opening the pdf and remain available
## even for a term set whose barplot failed to render (see
## barplot_enrich_pair()'s per-plot tryCatch).
report_cluster_marker_enrichment <- function(res, label, q = 0.2) {
  if (is.null(res)) return(invisible(NULL))
  for (cc in res$cluster_ids) {
    cat(sprintf("\n%s cluster %s vs rest\n", label, cc))
    for (ont in c("BP", "MF", "CC", "KEGG")) {
      cat(sprintf(" %s, up (higher in cluster %s):\n", ont, cc));   print_enrich_brief(res$up_enrich[[cc]][[ont]],   q = q)
      cat(sprintf(" %s, down (lower in cluster %s):\n", ont, cc));  print_enrich_brief(res$down_enrich[[cc]][[ont]], q = q)
    }
  }
  invisible(NULL)
}

## Cell-cycle phase scoring extended to a dataset's own validated
## clustering (Idents already set, e.g. YSC.MIX.SE, YSC.HYB.SC,
## YSC.HYB.SE), rather than only YSC.MERGE. Reuses the same three
## curated yeast regulons (S.GENES/G2M.GENES/MG1.GENES) via Seurat's
## CellCycleScoring()/AddModuleScore(), draws the same violin plot plus
## a phase-composition-by-cluster barplot, and prints the phase table to
## console (the direct quantitative check of a GO-based "this cluster
## looks like G1/S, S, M, ..." reading, rather than resting on term
## names alone). Returns obj with S.Score/G2M.Score/MG1.Score1/Phase
## added, for the caller to reassign (e.g. YSC.MIX.SE <- score_cell_cycle_by_cluster(YSC.MIX.SE, ...)).
score_cell_cycle_by_cluster <- function(obj, label, pdf_path, s_genes = S.GENES, g2m_genes = G2M.GENES, mg1_genes = MG1.GENES, width = 9, height = 8) {
  obj <- suppressWarnings(suppressMessages(CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes)))
  obj <- suppressWarnings(suppressMessages(AddModuleScore(obj, features = list(mg1_genes), name = "MG1.Score")))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  par(mfrow = c(1, 1))
  print(suppressWarnings(VlnPlot(obj, features = c("S.Score", "G2M.Score", "MG1.Score1"), ncol = 1, pt.size = 0)))

  phase_by_cluster <- prop.table(table(Idents(obj), obj$Phase), margin = 1)
  barplot(t(phase_by_cluster), col = c(G1 = "grey70", S = "#1f78b4", G2M = "#e31a1c")[colnames(phase_by_cluster)],
          legend.text = TRUE, args.legend = list(x = "topright", bty = "n"),
          las = 2, ylab = "fraction of cells", main = sprintf("%s: cell-cycle phase composition by cluster", label))
  dev.off()

  cat(sprintf("%s: cell-cycle phase composition by cluster\n", label)); print(round(phase_by_cluster, 3))
  obj
}

## Pulls all genes annotated (including descendant/child terms, via
## GOALL) to a GO biological-process term from org.Sc.sgd.db, the same
## annotation source already used for every enrichment test in this
## pipeline, rather than a hand-typed gene list. This avoids systematic-
## name transcription errors and keeps the same independent-of-these-
## clusters sourcing principle as the curated cell-cycle regulons
## (S.GENES/G2M.GENES/MG1.GENES came from established literature
## regulons, not from any cluster's own markers; the same applies here).
go_gene_set <- function(go_id, orgdb = org.Sc.sgd.db) {
  unique(AnnotationDbi::select(orgdb, keys = go_id, keytype = "GOALL", columns = "ORF")$ORF)
}

## Continuous module-score validation, analogous to
## score_cell_cycle_by_cluster() but for gene sets without a discrete
## phase classification (e.g. metabolic-state regulons like glycolysis,
## oxidative phosphorylation, ribosome biogenesis): scores obj by
## AddModuleScore() for each named gene set in gene_sets, draws a violin
## plot per module, and a per-cluster mean +/- SE barplot for each
## module, so a qualitative GO-term reading of a cluster's markers can
## be checked against continuous, independently-sourced module scores
## rather than resting on term names alone. gene_sets is a named list
## (e.g. list(Glycolysis = go_gene_set("GO:0006096"), ...)); names
## become the AddModuleScore() name prefix and the plot labels.
score_modules_by_cluster <- function(obj, gene_sets, label, pdf_path, width = 9, height = 8) {
  score_names <- paste0(names(gene_sets), "1")
  for (i in seq_along(gene_sets)) obj <- suppressWarnings(suppressMessages(AddModuleScore(obj, features = list(gene_sets[[i]]), name = names(gene_sets)[i])))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  par(mfrow = c(1, 1))
  print(suppressWarnings(VlnPlot(obj, features = score_names, ncol = 1, pt.size = 0)))

  cl    <- Idents(obj)
  means <- sapply(score_names, function(sn) tapply(obj[[sn]][, 1], cl, mean))
  ses   <- sapply(score_names, function(sn) tapply(obj[[sn]][, 1], cl, function(x) sd(x) / sqrt(length(x))))

  par(mfrow = c(1, length(score_names)), mar = c(5, 4.5, 3, 1))
  for (i in seq_along(score_names)) {
    b <- barplot(means[, i], ylim = range(c(means[, i] - ses[, i], means[, i] + ses[, i])),
                 main = sprintf("%s: %s", label, names(gene_sets)[i]), ylab = "mean module score", las = 2)
    arrows(b, means[, i] - ses[, i], b, means[, i] + ses[, i], angle = 90, code = 3, length = 0.05)
  }
  dev.off()

  cat(sprintf("%s: mean module score by cluster\n", label)); print(round(means, 3))
  obj
}

## ============================================================
## 8g. Covariate-noise diagnostic: does per-cell NB noise track
## cell-cycle position or metabolic state?
## ============================================================
## Continuous cell-cycle position from CellCycleScoring()'s two scores
## (S.Score, G2M.Score). These trace a roughly circular path through
## the cycle rather than sitting on a single line, so this takes their
## first principal component as the closest one-dimensional summary,
## rather than using one score alone or an arbitrary difference between
## them. Returns a named numeric vector, one value per cell
cell_cycle_continuum <- function(obj) {
  missing <- setdiff(c("S.Score", "G2M.Score"), colnames(obj[[]]))
  if (length(missing) > 0)
    stop(sprintf("cell_cycle_continuum(): %s not found; run score_cell_cycle_by_cluster() (or CellCycleScoring()) on this object first and reassign its result", paste(missing, collapse = ", ")))
  sc  <- obj[[c("S.Score", "G2M.Score")]]
  pc1 <- prcomp(sc, scale. = TRUE)$x[, 1]
  setNames(pc1, rownames(sc))
}

## Two-object version of cell_cycle_continuum(), for comparing cell-
## cycle position BETWEEN two objects (e.g. species_composition_bound()
## in 8h). prcomp() always centers its own output to exactly mean zero,
## so calling cell_cycle_continuum() separately on two objects gives
## two axes that both have mean zero by construction, regardless of any
## real difference between them, making a between-object mean
## comparison meaningless. This instead fits one PCA on the pooled
## S.Score/G2M.Score from both objects, then projects each object's
## cells onto that shared PC1, so a real difference in cell-cycle
## composition can actually appear in the resulting means. Only use
## this for between-object comparisons; cell_cycle_continuum() remains
## correct for the within-object correlation test in 8g, where the
## per-object fit doesn't matter since correlation is invariant to
## per-dataset linear rescaling or recentering
cell_cycle_continuum_shared <- function(obj1, obj2) {
  sc1    <- obj1[[c("S.Score", "G2M.Score")]]
  sc2    <- obj2[[c("S.Score", "G2M.Score")]]
  scores <- prcomp(rbind(sc1, sc2), scale. = TRUE)$x[, 1]
  list(x1 = setNames(scores[seq_len(nrow(sc1))], rownames(sc1)),
       x2 = setNames(scores[(nrow(sc1) + 1):length(scores)], rownames(sc2)))
}

## Discrete metabolic state from the continuous module scores computed
## by score_modules_by_cluster() (e.g. Glycolysis1/OXPHOS1/RiBi1). A
## cell's metabolic state can move along more than one of these axes at
## once (e.g. simultaneous glycolytic and respiratory use), so this
## clusters cells jointly on all of them via k-means rather than
## thresholding each pathway on its own, which would miss states
## defined by a combination. The number of states k is chosen by mean
## silhouette width over k_range, the same coarsest-preferred logic
## already used for the unsupervised transcriptome clustering earlier
## in this section, so a state count isn't fixed by hand. Returns a
## factor of state labels, one per cell
metabolic_state_cluster <- function(obj, score_names = c("Glycolysis1", "OXPHOS1", "RiBi1"), k_range = 2:4) {
  missing <- setdiff(score_names, colnames(obj[[]]))
  if (length(missing) > 0)
    stop(sprintf("metabolic_state_cluster(): %s not found; run score_modules_by_cluster() on this object first and reassign its result", paste(missing, collapse = ", ")))
  scores <- scale(as.matrix(obj[[score_names]]))
  d      <- dist(scores)
  sil    <- sapply(k_range, function(k) mean(silhouette(kmeans(scores, centers = k, nstart = 10)$cluster, d)[, 3]))
  best_k <- k_range[which.max(sil)]
  km     <- kmeans(scores, centers = best_k, nstart = 10)
  setNames(factor(km$cluster), rownames(scores))
}

## Tests whether per-cell NB Pearson residual noise (nb_residuals(),
## already fit against the exposure offset only) depends on a cell's
## continuous cell-cycle position or discrete metabolic state. This is
## a diagnostic against the existing fit, not a refit with these as
## covariates: for the cell-cycle axis, a per-gene Spearman correlation
## between residual and cell-cycle position; for metabolic state, a
## per-gene Kruskal-Wallis test across states plus eta-squared as the
## fraction of residual variance the state explains. Both are
## summarized by the fraction of genes reaching BH significance and the
## median effect size, the direct check of whether noise divergence is
## being driven by unmodeled cell state rather than genuine regulatory
## divergence
covariate_noise_diagnostic <- function(mat, expo, fit, cc_axis, metab_state, label) {
  cells <- Reduce(intersect, list(colnames(mat), names(cc_axis), names(metab_state)))
  R   <- nb_residuals(mat[, cells], expo[cells], fit)
  cc  <- cc_axis[cells]
  met <- droplevels(metab_state[cells])

  cc_rho <- apply(R, 1, function(r) suppressWarnings(cor(r, cc, method = "spearman", use = "complete.obs")))
  cc_p   <- apply(R, 1, function(r) suppressWarnings(cor.test(r, cc, method = "spearman")$p.value))
  cc_q   <- p.adjust(cc_p, method = "BH")

  eta_sq  <- function(r, g) { ss <- summary(aov(r ~ g))[[1]][["Sum Sq"]]; ss[1] / sum(ss) }
  met_eta <- apply(R, 1, function(r) eta_sq(r, met))
  met_p   <- apply(R, 1, function(r) suppressWarnings(kruskal.test(r, met)$p.value))
  met_q   <- p.adjust(met_p, method = "BH")

  cat(sprintf("%s: cell-cycle axis, %d / %d genes (%.1f%%) at q < 0.05, median |rho| = %.3f\n",
              label, sum(cc_q < 0.05, na.rm = TRUE), length(cc_q), 100 * mean(cc_q < 0.05, na.rm = TRUE), median(abs(cc_rho), na.rm = TRUE)))
  cat(sprintf("%s: metabolic state, %d / %d genes (%.1f%%) at q < 0.05, median eta^2 = %.3f\n",
              label, sum(met_q < 0.05, na.rm = TRUE), length(met_q), 100 * mean(met_q < 0.05, na.rm = TRUE), median(met_eta, na.rm = TRUE)))

  data.frame(gene = rownames(R), cc_rho = cc_rho, cc_q = cc_q, met_eta = met_eta, met_q = met_q, row.names = NULL)
}

## Histogram pair for covariate_noise_diagnostic()'s output: per-gene
## Spearman rho against the cell-cycle axis, and per-gene eta-squared
## against metabolic state, so the effect-size distribution behind the
## significance counts is visible rather than only the summary numbers
plot_covariate_noise_diagnostic <- function(diag_df, label) {
  par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))
  hist(diag_df$cc_rho, breaks = 40, col = "grey70", border = NA,
       xlab = "Spearman rho, residual vs cell-cycle axis", main = sprintf("%s: cell cycle", label))
  hist(diag_df$met_eta, breaks = 40, col = "grey70", border = NA,
       xlab = "eta-squared, residual vs metabolic state", main = sprintf("%s: metabolic state", label))
}

## ============================================================
## 8h. Species composition comparison and confound bound
## ============================================================
## A small covariate-noise relationship (8g) only rules out a confound
## on the cross-species noise comparison if it's paired with a check of
## whether Sc and Se actually differ in composition along that axis in
## the first place; 8g alone shows the relationship is weak within each
## species but never looks at whether the axis differs between them.
## This compares one continuous axis (the cell-cycle axis, or one
## metabolic module score) between two species/allele views via a
## standardized mean difference (Cohen's d) and a rank-sum test, then
## multiplies that composition difference by the matching covariate-
## noise effect size from 8g to bound the expected shift in per-gene
## noise a composition difference of that size could produce. This
## bound is approximate and linear (composition difference in SD units
## x noise-association effect size, the expected standardized-outcome
## shift under a bivariate-normal approximation), not a refit or a
## direct test of any specific contrast; it complements the 8g
## effect-size argument rather than replacing it
species_composition_bound <- function(x1, x2, effect_size, label, axis_label) {
  d     <- (mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) / sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
  wt    <- wilcox.test(x1, x2)
  bound <- abs(d) * effect_size
  cat(sprintf("%s, %s: Cohen's d = %.3f (Wilcoxon p = %.2g), noise-association effect size = %.3f, bound on expected noise shift = %.3f SD\n",
              label, axis_label, d, wt$p.value, effect_size, bound))
  data.frame(label = label, axis = axis_label, cohens_d = d, wilcox_p = wt$p.value, effect_size = effect_size, bound = bound)
}

## Runs species_composition_bound() across the cell-cycle axis and the
## three metabolic module scores for one species/allele pair. The
## noise-association effect size for each axis is the more conservative
## (larger) of the two species' own 8g estimates, so the resulting bound
## is an upper bound rather than an average. For the metabolic module
## scores, sqrt(median eta-squared) is used as the standardized effect
## size: eta-squared is the categorical-predictor analog of a squared
## correlation, so its square root is on the same standardized scale as
## Cohen's d and the cell-cycle rho
species_composition_report <- function(cc1, cc2, met1, met2, diag1, diag2, label) {
  cc_effect  <- max(median(abs(diag1$cc_rho), na.rm = TRUE), median(abs(diag2$cc_rho), na.rm = TRUE))
  met_effect <- sqrt(max(median(diag1$met_eta, na.rm = TRUE), median(diag2$met_eta, na.rm = TRUE)))
  rows <- list(species_composition_bound(cc1, cc2, cc_effect, label, "cell-cycle axis"))
  for (sn in colnames(met1)) rows[[length(rows) + 1]] <- species_composition_bound(met1[, sn], met2[, sn], met_effect, label, sn)
  do.call(rbind, rows)
}

## ============================================================
## plot_geneset_direction_stack: stacked, class-coloured, cross-hatched
## enrichment of one specific gene set (e.g. a GO/KEGG hit's core genes)
## across the regulatory classification
## ============================================================
## For each level of the regulatory classification (Conserved, Cis, Trans,
## Cis + Trans, Compensatory), shows how much of that class's membership in
## a pre-defined gene set comes from genes with significantly higher
## parental noise or mean in Sc (bottom segment) versus Se (top segment).
## Direction is defined exactly as in build_reg_go_sets: the sign of the
## total parental contrast, gated by its own permutation p at `sig`,
## independent of which component (cis or trans) drives the class call.
##
## Each segment's height and its significance test share the same
## denominator, the direction-restricted subgroup size within that class
## (n), not the class total (n_class), so a class with very few
## direction-significant genes gets a short, honestly uncertain segment
## rather than an inflated one. n_class and n are both returned so this is
## checkable. A per-segment two-sided hypergeometric test (gene_set vs
## universe, n draws) is BH-corrected across every level x direction cell;
## cells with padj < fdr get a star. bg_rate (gene_set size / universe
## size) is drawn as a dashed reference line, since a class x direction bar
## is only interpretable relative to how often the gene set turns up by
## chance genome-wide.
##
## Both segments use the class's own colour from cols (default
## COLOR.LIST.1, positionally matched to levels exactly as classify_reg
## does); the Se segment is cross-hatched (two overlaid polygon() calls at
## angle and angle + 90) so it stays visually distinct from Sc without a
## second colour scheme. Sample size n is printed inside each segment,
## since some splits (Reinforcing and Compensatory especially) will be
## thin, and a tall bar built on a handful of genes should not visually
## outweigh one built on many.
##
## Args: BURST.CONTRASTS, PR as elsewhere; gene_set (character vector of gene IDs to
##   test, e.g. a significant GO term's core genes); universe (character
##   vector, the tested background, e.g. BURST.CONTRASTS$gene); quantity ("mean" or
##   "bfreq"); levels (REG.CLASS by default); sig (threshold for the
##   direction call); fdr (BH threshold for stars); cols; hatch_density,
##   hatch_angle (polygon() hatching parameters for the Se segment).
## Plots to the current device. Returns invisibly a data.frame with one row
##   per level x direction: level, direction, n_class, n, n_in_set,
##   frac_within_dir, height, p, padj.
plot_geneset_direction_stack <- function(BURST.CONTRASTS, PR, gene_set, universe, quantity = c("mean", "bfreq", "bsize", "kbal"), levels = REG.CLASS, sig = 0.05, fdr = 0.05, cols = COLOR.LIST.1, hatch_density = 18, hatch_angle = 45, main = NULL, ylab = NULL) {
  quantity <- match.arg(quantity)
  gene_set <- intersect(gene_set, universe)
  bg_rate  <- length(gene_set) / length(universe)

  class_vec <- reg_class_vec(BURST.CONTRASTS, PR, quantity, sig)
  genes <- BURST.CONTRASTS$gene
  i    <- match(genes, PR$gene)
  est  <- BURST.CONTRASTS[[paste0(quantity, "_total_est")]]
  p0   <- PR[[paste0(quantity, "_total_p")]][i]
  dirn <- ifelse(!is.finite(p0) | p0 >= sig, "ns", ifelse(est > 0, "Sc", "Se"))

  calc <- function(k, d) {
    n_class <- length(intersect(genes[!is.na(class_vec) & class_vec == k], universe))
    sel <- !is.na(class_vec) & class_vec == k & dirn == d
    g   <- intersect(genes[sel], universe)
    n   <- length(g); hit <- length(intersect(g, gene_set))
    pp  <- if (n > 0) {
      po <- phyper(hit - 1, length(gene_set), length(universe) - length(gene_set), n, lower.tail = FALSE)
      pu <- phyper(hit,     length(gene_set), length(universe) - length(gene_set), n, lower.tail = TRUE)
      2 * min(po, pu, 0.5)
    } else NA_real_
    data.frame(level = k, direction = d, n_class = n_class, n = n, n_in_set = hit,
               frac_within_dir = if (n > 0) hit / n else NA_real_, p = pp,
               stringsAsFactors = FALSE)
  }

  RES <- do.call(rbind, lapply(levels, function(k) rbind(calc(k, "Sc"), calc(k, "Se"))))
  RES$padj   <- p.adjust(RES$p, "BH")
  RES$height <- ifelse(is.na(RES$frac_within_dir), 0, RES$frac_within_dir)

  pick <- function(v) { x <- RES[[v]][RES$direction == "Sc"]; names(x) <- RES$level[RES$direction == "Sc"]; x[levels] }
  pick_se <- function(v) { x <- RES[[v]][RES$direction == "Se"]; names(x) <- RES$level[RES$direction == "Se"]; x[levels] }
  h_sc <- pick("height");    h_se <- pick_se("height")
  n_sc <- pick("n");         n_se <- pick_se("n")
  padj_sc <- pick("padj");   padj_se <- pick_se("padj")

  if (is.null(ylab)) ylab <- "fraction of direction-specific subgroup in gene set"
  if (is.null(main))  main <- paste0(quantity, ": gene set representation by class")

  n_lev <- length(levels)
  bar_w <- 0.7
  xpos  <- seq_len(n_lev)
  ymax  <- suppressWarnings(max(c(h_sc + h_se, bg_rate), na.rm = TRUE)) * 1.25
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  op <- par(mar = c(8, 5, 3, 1)); on.exit(par(op))
  plot(NA, xlim = c(0.3, n_lev + 0.7), ylim = c(0, ymax), xaxt = "n", xlab = "", ylab = ylab, main = main)
  axis(1, at = xpos, labels = levels, las = 2)
  abline(h = bg_rate, lty = 2, col = "grey40")

  for (idx in seq_len(n_lev)) {
    col_k <- cols[idx]
    x0 <- xpos[idx] - bar_w / 2; x1 <- xpos[idx] + bar_w / 2
    hs <- if (is.na(h_sc[idx])) 0 else h_sc[idx]
    he <- if (is.na(h_se[idx])) 0 else h_se[idx]

    rect(x0, 0, x1, hs, col = col_k, border = "black")
    if (he > 0) {
      rect(x0, hs, x1, hs + he, col = col_k, border = "black")
      polygon(c(x0, x1, x1, x0), c(hs, hs, hs + he, hs + he), density = hatch_density, angle = hatch_angle, col = "black", border = NA)
      polygon(c(x0, x1, x1, x0), c(hs, hs, hs + he, hs + he), density = hatch_density, angle = hatch_angle + 90, col = "black", border = NA)
    }

    if (!is.na(n_sc[idx]) && n_sc[idx] > 0) text(xpos[idx], hs / 2,      paste0("n=", n_sc[idx]), cex = 0.65)
    if (!is.na(n_se[idx]) && n_se[idx] > 0) text(xpos[idx], hs + he / 2, paste0("n=", n_se[idx]), cex = 0.65)

    if (!is.na(padj_sc[idx]) && padj_sc[idx] < fdr) text(xpos[idx], hs + ymax * 0.02,      "*", cex = 1.3, font = 2)
    if (!is.na(padj_se[idx]) && padj_se[idx] < fdr) text(xpos[idx], hs + he + ymax * 0.02, "*", cex = 1.3, font = 2)
  }

 legend("topright", legend = c("Sc-higher", "Se-higher"), fill = "grey70", density = c(NA, hatch_density), angle = hatch_angle, border = "black", bty = "n", cex = 0.8)

  invisible(RES[, c("level", "direction", "n_class", "n", "n_in_set", "frac_within_dir", "height", "p", "padj")])
}

## ============================================================
## 9. PUBLICATION FIGURES
## ============================================================

## ---- figure colour palettes ---------------------------------------------
## COLOR.LIST.1 (regulatory classes) and COLOR.LIST.2 (dominance classes)
## are defined and named in the main script, matching REG.CLASS and
## DOM.CLASS

## ============================================================
## Figs 3 & 8 : mean-class x size-class enrichment heatmap
## ============================================================
## Cell colour encodes log2(obs/exp) for all cells (no masking).
## Cell text shows obs/exp as fold enrichment (2^lor) with a trailing
## "*" where the BH-adjusted two-sided hypergeometric p < fdr.
## A value of 2.0* means the co-occurrence is twice as frequent as
## expected by chance; 0.4* means 0.4x (2.5x depleted).
## Wraps class_overlap_heatmap() with explicit levels, so every class
## appears in the table even when empty. Used for the main-text
## Figures 2 and 6, where the two heatmaps in a row need matching axes.
plot_class_overlap <- function(class_a, class_b, levels_a, levels_b, fdr = 0.05, cols = COLOR.LIST.3, brk = length(cols), xlab = "burst frequency class", ylab = "mean class", cex_axis = 0.75, cex_cell = 0.7, rng = NULL) {
  class_overlap_heatmap(class_a, class_b, levels_a = levels_a, levels_b = levels_b, brk = brk, fdr = fdr, cols = cols, cex_axis = cex_axis, cex_cell = cex_cell, xlab = xlab, ylab = ylab, rng = rng)
}

## Shared log2(obs/exp) half-range across a list of (class_a, class_b)
## pairs, all built against the same levels, so a set of related overlap
## panels (e.g. Figure 2's three regulatory pairings) can be handed one
## common rng and therefore one common color scale.
shared_overlap_rng <- function(pairs, levels_common) {
  lors <- lapply(pairs, function(p) {
    tab <- table(factor(p[[1]], levels = levels_common), factor(p[[2]], levels = levels_common))
    exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    log2((tab + 0.5) / (exp + 0.5))
  })
  max(sapply(lors, function(l) max(abs(l))), 1e-6)
}

## ============================================================
## Fig 5 : rotated burst kinetics, z-test shaded
## ============================================================
## kbal_sig() computes the bootstrap z-test on the kinetic balance y.
## plot_burst_kinetics_sig() returns this data frame invisibly so the
## driver can build the complementary barplot without re-computing.
kbal_sig <- function(BURST.CONTRASTS, mode = "total", sig = 0.05) {
  bf <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_est")]]; bs <- BURST.CONTRASTS[[paste0("bsize_",mode,"_est")]]
  sm <- BURST.CONTRASTS[[paste0("mean_",mode,"_se")]];   ss <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_se")]]
  r  <- BURST.CONTRASTS[[paste0("cor_",mode)]]
  ## As in add_burst_contrasts: cor_<mode> is NA for cis/trans, a
  ## structural zero covariance rather than a missing value
  r[!is.finite(r)] <- 0
  ok <- is.finite(bf)&is.finite(bs)&is.finite(sm)&is.finite(ss)
  bf<-bf[ok]; bs<-bs[ok]; sm<-sm[ok]; ss<-ss[ok]; r<-r[ok]
  covFB <- r*sm*ss - ss^2
  seS   <- sqrt(pmax(0, sm^2+ss^2-2*r*sm*ss))
  y  <- bf - bs
  sy <- sqrt(pmax(0, ss^2+seS^2-2*covFB))
  p  <- 2*pnorm(-abs(y)/sy)
  dir <- ifelse(!is.finite(p)|p>=sig, "ns", ifelse(y>0, "sig_pos", "sig_neg"))
  data.frame(gene=BURST.CONTRASTS$gene[ok], y=y, sy=sy, p=p, direction=dir,
             stringsAsFactors=FALSE)
}

plot_burst_kinetics_sig <- function(BURST.CONTRASTS, mode = "total", sig = 0.05, main = NULL, bar_col = "#88888840", sig_col = "black", ns_col = "grey70") {
  ks <- kbal_sig(BURST.CONTRASTS, mode, sig)
  bf <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_est")]]; bs <- BURST.CONTRASTS[[paste0("bsize_",mode,"_est")]]
  sm <- BURST.CONTRASTS[[paste0("mean_",mode,"_se")]];   ss <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_se")]]
  r  <- BURST.CONTRASTS[[paste0("cor_",mode)]]
  ## As in add_burst_contrasts: cor_<mode> is NA for cis/trans, a
  ## structural zero covariance rather than a missing value
  r[!is.finite(r)] <- 0
  ok <- is.finite(bf)&is.finite(bs)&is.finite(sm)&is.finite(ss)
  bf<-bf[ok]; bs<-bs[ok]; sm<-sm[ok]; ss<-ss[ok]; r<-r[ok]
  covFB <- r*sm*ss - ss^2
  seS   <- sqrt(pmax(0, sm^2+ss^2-2*r*sm*ss))
  x  <- bf + bs; sx <- sqrt(pmax(0, ss^2+seS^2+2*covFB))
  y  <- ks$y;    sy <- ks$sy
  sig_idx <- ks$direction != "ns"
  if (is.null(main)) main <- paste0("burst kinetics: ",mode)
  op <- par(pty="s"); on.exit(par(op))
  plot(NA, xlim=.sym(x, sx), ylim=.sym(y, sy), xlab="net mean change (log2)", ylab="frequency - amplitude (kinetic balance)", main=main)
  abline(h=0,v=0,col="grey45")
  segments(x-sx, y, x+sx, y, col=bar_col)        # bars first
  segments(x, y-sy, x, y+sy, col=bar_col)
  points(x[!sig_idx], y[!sig_idx], pch=16, cex=0.5, col=ns_col)  # ns below
  points(x[ sig_idx], y[ sig_idx], pch=16, cex=0.5, col=sig_col)  # sig on top
  legend("topleft", legend=c(paste0("sig (p<", sig, ")"), "n.s."), col=c(sig_col, ns_col), pch=16, bty="n", cex=0.8)
  invisible(ks)
}

## ============================================================
## Figs 6 & 7 : dominance, parent frame (primary) or A/D rotation (S)
## ============================================================
## frame="parent" (x = hybrid - Sc, y = hybrid - Se) is the primary view
## for both mean and noise: classify_dom() tests exactly these two
## contrasts, with no midparent construction involved (Additive only
## needs the two contrasts to have opposite sign, i.e. hybrid lies
## between the two parents on whichever scale is being tested, a claim
## that needs no separate mean/noise definition since sign is preserved
## under any monotonic rescaling). A rotation preserves straight lines
## through the origin, so the six class regions are exactly as clean in
## this frame as in the A/D rotation, just oriented along the axes and
## the y=-x diagonal instead of along a and d.
##
## frame="AD" is kept as a supplementary view. For noise it is an exact
## rotation of the tested contrasts (the "dom" permutation null's
## midparent is the geometric/log average, identical in log2 space to
## (dpar_sc+dpar_se)/2). For mean it is not: "dom"'s midparent there is
## the arithmetic average of the two parents' point estimates (the
## dosage-additivity model, hybrid = one Sc-like allele plus one Se-like
## allele), a nonlinear function of log2(Sc) and log2(Se) that no
## rotation of dpar_sc/dpar_se can reproduce exactly. When frame="AD" and
## cor_dpar_mean / cor_dpar_bfreq exist in BURST.CONTRASTS (produced by the updated
## boot_contrasts_one), SE bars are exact for the rotation being plotted;
## otherwise they fall back to the independence approximation.
plot_dom_class <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), frame = c("parent","AD"), sig = 0.05, main = NULL, lim = NULL, bar_col = "#88888840") {
  quantity <- match.arg(quantity); frame <- match.arg(frame)
  dsc <- BURST.CONTRASTS[[paste0(quantity,"_dpar_sc_est")]]; ssc <- BURST.CONTRASTS[[paste0(quantity,"_dpar_sc_se")]]
  dse <- BURST.CONTRASTS[[paste0(quantity,"_dpar_se_est")]]; sse <- BURST.CONTRASTS[[paste0(quantity,"_dpar_se_se")]]
  rdp_col <- paste0("cor_dpar_",quantity)
  rdp_all <- if (rdp_col %in% names(BURST.CONTRASTS)) BURST.CONTRASTS[[rdp_col]] else NULL
  cls <- dom_class_vec(BURST.CONTRASTS, PR, quantity, sig)
  ok  <- is.finite(dsc)&is.finite(dse)&is.finite(ssc)&is.finite(sse)&!is.na(cls)
  dsc<-dsc[ok]; dse<-dse[ok]; ssc<-ssc[ok]; sse<-sse[ok]; cls<-cls[ok]
  rdp <- if (!is.null(rdp_all)) { r<-rdp_all[ok]; ifelse(is.finite(r),r,0) } else rep(0,sum(ok))
  if (frame=="parent") {
    x<-dsc; sx<-ssc; y<-dse; sy<-sse
    xl<-"hybrid - Sc parent (log2)"; yl<-"hybrid - Se parent (log2)"
  } else {
    x  <- (dse-dsc)/2; y  <- (dsc+dse)/2
    sx <- sqrt(pmax(0, ssc^2+sse^2-2*rdp*ssc*sse))/2
    sy <- sqrt(pmax(0, ssc^2+sse^2+2*rdp*ssc*sse))/2
    xl<-"additive  (Sc - Se)/2  (log2)"; yl<-"dominance  hybrid - midparent  (log2)"
  }
  col <- COLOR.LIST.2[match(cls,DOM.CLASS)]
  if (is.null(lim)) { m<-max(abs(c(x+sx,x-sx,y+sy,y-sy)),na.rm=TRUE); lim<-c(-m,m) }
  if (is.null(main)) main<-paste0("dominance (",quantity,", ",frame,")")
  op<-par(pty="s"); on.exit(par(op))
  plot(NA,xlim=lim,ylim=lim,xlab=xl,ylab=yl,main=main)
  abline(0,1,lty=3,col="grey65"); abline(0,-1,lty=3,col="grey65")
  abline(h=0,v=0,col="grey45")
  segments(x-sx,y,x+sx,y,col=bar_col); segments(x,y-sy,x,y+sy,col=bar_col)
  points(x,y,pch=16,cex=0.5,col=col)
  legend("topleft",legend=DOM.CLASS,col=COLOR.LIST.2,pch=16,bty="n",cex=0.8)
  invisible(cls)
}

## ============================================================
## Fig 9 : ggplot2 violins of a burst quantity by class
## ============================================================
plot_violins <- function(value, reg_class, dom_class, ylab = "kinetic balance (burst frequency - amplitude)") {
  stopifnot(requireNamespace("ggplot2", quietly=TRUE))
  df <- rbind(
    data.frame(value=value, class=factor(reg_class,levels=REG.CLASS), panel="Regulatory"),
    data.frame(value=value, class=factor(dom_class, levels=DOM.CLASS), panel="Dominance"))
  df <- df[is.finite(df$value) & !is.na(df$class), ]
  df$panel <- factor(df$panel, levels=c("Regulatory","Dominance"))
  cols <- c(COLOR.LIST.1, COLOR.LIST.2[setdiff(names(COLOR.LIST.2),names(COLOR.LIST.1))])
  ggplot2::ggplot(df, ggplot2::aes(class, value, fill=class)) +
    ggplot2::geom_hline(yintercept=0, colour="grey75", linewidth=0.3) +
    ggplot2::geom_violin(scale="width", trim=TRUE, colour="grey30", linewidth=0.3) +
    ggplot2::geom_boxplot(width=0.12, outlier.size=0.3, fill="white", colour="grey30", linewidth=0.3) +
    ggplot2::facet_wrap(~panel, scales="free_x") +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::labs(x=NULL, y=ylab) +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35, hjust=1), legend.position="none", strip.background=ggplot2::element_blank())
}

## ============================================================
## Fig 10 : intrinsic / extrinsic noise from hybrid alleles
## ============================================================
## Two-allele decomposition pooling HYC and HYT cells. Each allele
## is rescaled to its own mean rate to remove the cis mean difference.
## Poisson shot noise is subtracted: Var_Poisson(a'_i) = a_i/(e_i^2 * ma^2).
## Returns the raw intr/extr components (not just their fraction), so
## they can also be used directly as a ratio elsewhere (Section 12)
## rather than only through intrinsic_fraction()'s bounded [0,1] fraction.
intrinsic_extrinsic_components <- function(mats, expos) {
  e_vec <- c(expos$HYC, expos$HYT)
  a_raw <- cbind(mats$HYC.SC, mats$HYT.SC)
  b_raw <- cbind(mats$HYC.SE, mats$HYT.SE)
  a  <- sweep(a_raw, 2, e_vec, "/")
  b  <- sweep(b_raw, 2, e_vec, "/")
  ma <- rowMeans(a); mb <- rowMeans(b)
  good <- is.finite(ma) & is.finite(mb) & ma > 0 & mb > 0
  ap <- a / pmax(ma, 1e-10); bp <- b / pmax(mb, 1e-10)
  intr_raw <- 0.5 * rowMeans((ap-bp)^2)
  extr     <- rowMeans(ap*bp) - 1
  shot_a <- rowMeans(sweep(a_raw, 2, e_vec^2, "/")) / pmax(ma^2, 1e-10)
  shot_b <- rowMeans(sweep(b_raw, 2, e_vec^2, "/")) / pmax(mb^2, 1e-10)
  intr <- intr_raw - 0.5*(shot_a+shot_b)
  data.frame(gene = rownames(mats$HYC.SC), intr = intr, extr = extr, good = good, row.names = NULL)
}

intrinsic_fraction <- function(mats, expos) {
  ie   <- intrinsic_extrinsic_components(mats, expos)
  frac <- ie$intr / (ie$intr + pmax(ie$extr, 0))
  frac[!ie$good | !is.finite(frac)] <- NA_real_
  setNames(frac, ie$gene)
}

## Histogram with regulatory medians (triangles above) and dominance
## medians (inverted triangles below). Legend at topleft avoids the
## right-side bars where the distribution is tallest.
plot_intrinsic_hist <- function(frac, reg_class, dom_class, brk = 30, main = NULL) {
  ok <- is.finite(frac) & frac>=0 & frac<=1
  f  <- frac[ok]; rc <- reg_class[ok]; dc <- dom_class[ok]
  h  <- hist(f, breaks=seq(0,1,length.out=brk+1), plot=FALSE)
  top <- max(h$counts)
  plot(h, col="grey85", border="white", xlim=c(0, 1), ylim=c(-top*0.06, top*1.20), xlab="intrinsic / (intrinsic + extrinsic)", ylab="# of genes", main=if(is.null(main)) "" else main)
  med <- function(cls, lv) vapply(lv, function(k) median(f[cls==k],na.rm=TRUE), numeric(1))
  rm <- med(rc, REG.CLASS); dm <- med(dc, DOM.CLASS)
  yR <- top*1.10; yD <- -top*0.04
  for (i in seq_along(rm))
    if (is.finite(rm[i])) segments(rm[i],0,rm[i],yR, col=COLOR.LIST.1[i], lty=3, lwd=0.7)
  points(rm, rep(yR,length(rm)), pch=17, col=COLOR.LIST.1, cex=1.1, xpd=NA)
  points(dm, rep(yD,length(dm)), pch=25, col=COLOR.LIST.2, bg=COLOR.LIST.2, cex=1.1, xpd=NA)
  legend("topleft", legend=c(REG.CLASS, NA, DOM.CLASS), col=c(COLOR.LIST.1, NA, COLOR.LIST.2), pch=c(rep(17, length(REG.CLASS)), NA, rep(25, length(DOM.CLASS))), pt.bg=c(rep(NA, length(REG.CLASS)), NA, COLOR.LIST.2), bty="n", cex=0.72)
}

##############################################################################
## 10. PROMOTER ARCHITECTURE (TATA box, poly(dA:dT), and predicted
##     nucleosome occupancy)
##############################################################################
## Three sequence features describe promoter architecture here. A TATA
## box consensus match and the longest poly(dA:dT) tract are direct
## sequence heuristics, scored with functions written for this project.
## Predicted nucleosome occupancy comes from NuPoP (Xi et al. 2010,
## Bioinformatics; extending the duration hidden Markov approach of
## Wang et al. 2008), which sits in the same lineage as the Kaplan and
## Segal probabilistic models and carries a real profile trained on
## yeast nucleosome data, rather than a model fit here from scratch.
## All three features depend on sequence alone, so the same functions
## score Sc and Se promoters without needing an experimental occupancy
## map for Se, which does not exist. Neither species has an annotated
## promoter, so promoters are defined directly from a genome FASTA and
## its GFF gene annotation below, rather than read from a pre-extracted
## file.
##
## NuPoP is a Bioconductor package, not part of base R or CRAN. Its
## trained profile (species = 7) is fit on S. cerevisiae; applying the
## same profile to Se scores both species against one shared reference
## rather than two different ones, which keeps the comparison
## symmetric, but it does mean the Se scores carry a real cross-species
## extrapolation, worth noting in the
## methods alongside the Bucher matrix's cross-eukaryote origin.

## Reads a genome FASTA into a named DNAStringSet, one sequence per
## chromosome or scaffold, keyed by the name before any whitespace in
## the FASTA header, matching how chromosomes are referenced in a GFF's
## seqid column
read_genome_fasta <- function(path) {
  genome <- Biostrings::readDNAStringSet(path)
  names(genome) <- sub("\\s.*$", "", names(genome))
  genome
}

## Reads a GFF3 annotation (Gene/exon/intron feature types, no separate
## mRNA level) and returns one row per gene: chromosome, start, end
## (1-based, start < end regardless of strand), and strand, everything
## extract_promoters() needs to locate each gene's promoter
##
## Anchored to exon rather than gene. This annotation's gene feature is
## the transcript span, which can include annotated UTR, and how
## completely UTRs are annotated is exactly the kind of thing that
## differs between two independent genome annotations, not something
## that reflects real promoter divergence; anchoring there mixes
## annotation-extent differences into what should be a promoter
## comparison. exon is the coding region and starts at the ATG, which is
## the reference point Basehoar et al.'s location window also uses.
##
## This Geneious-exported GFF3 has no Parent attribute linking an exon
## to its gene; instead an exon's Name is the gene's own name with a
## literal "_CDS" suffix (confirmed against the actual file: true for
## 6302 of 6306 exon rows), so id_field defaults to "Name" and id_suffix
## strips that suffix to recover the bare gene identifier. A multi-exon
## gene's fragments share the same Name (and so the same stripped
## identifier), which is what lets them collapse to one min(start)/
## max(end) span below, correctly spanning any introns. The few rows
## where Name has no "_CDS" to begin with pass through unchanged, since
## sub() only acts where the suffix is actually present. Set id_suffix
## to "" for an annotation that doesn't use this convention. Errors
## immediately, listing the attribute columns actually present, if
## id_field isn't one of them, rather than letting a wrong field name
## silently turn into a column of NULL and cascade into NA scores with
## no indication of where it went wrong.
read_gff_genes <- function(path, feature_type = "exon", id_field = "Name", id_suffix = "_CDS") {
  gr  <- rtracklayer::import(path)
  gr  <- gr[gr$type == feature_type]
  available <- names(S4Vectors::mcols(gr))
  if (!id_field %in% available)
  stop(sprintf("id_field '%s' not found on '%s' features; available attribute columns are: %s", id_field, feature_type, paste(available, collapse = ", ")))
  ids <- S4Vectors::mcols(gr)[[id_field]]
  if (methods::is(ids, "List")) ids <- vapply(ids, `[`, character(1), 1)
  ids <- as.character(ids)
  if (nzchar(id_suffix)) ids <- sub(paste0(id_suffix, "$"), "", ids)
  df <- data.frame(gene   = ids,
                    seqid  = as.character(GenomeInfoDb::seqnames(gr)),
                    start  = BiocGenerics::start(gr),
                    end    = BiocGenerics::end(gr),
                    strand = as.character(BiocGenerics::strand(gr)),
                    stringsAsFactors = FALSE)
  starts <- tapply(df$start, df$gene, min)
  ends   <- tapply(df$end,   df$gene, max)
  first  <- df[!duplicated(df$gene), c("gene", "seqid", "strand")]
  data.frame(gene      = first$gene,
             seqid     = first$seqid,
             start     = starts[first$gene],
             end       = ends[first$gene],
             strand    = first$strand,
             row.names = NULL,
             stringsAsFactors = FALSE)
}

## Extracts the promoter for every gene in genes: the sequence
## immediately upstream of its start codon, up to max_bp long or the
## distance to the nearest neighboring gene on the same chromosome,
## whichever is shorter, so a promoter never crosses into another
## gene's coding sequence. This also gives a divergently transcribed
## gene pair the full width of its shared intergenic region rather than
## a fixed window, and a naturally short promoter when genes sit close
## together, both of which are real features of the genome rather than
## something to special-case.
##
## min_bp is the one exception to that neighbor boundary. Below min_bp,
## the extracted region is too short to score a TATA box, poly(dA:dT)
## tract, or NuPoP occupancy window meaningfully at all (a 4 or 5 bp
## promoter has essentially nothing for any of the three to find), so
## when the neighbor-capped region falls short of min_bp, extraction
## continues past the neighboring gene's boundary and into its coding
## sequence, up to min_bp total. Real regulatory elements do sit inside
## upstream open reading frames often enough that this isn't just
## padding with noise, though a hit found there could equally be
## transcriptional readthrough rather than a promoter element proper,
## which is not something sequence alone can distinguish. Every gene
## that needed this extension is flagged in the "extended" attribute
## attached to the returned vector, so anything scored from an extended
## region stays identifiable downstream rather than blending in with
## promoters that respected the neighbor boundary throughout.
##
## Every returned sequence is oriented 5' to 3' relative to its own
## gene, with the end of the string nearest the start codon regardless
## of strand; minus-strand promoters are reverse-complemented to make
## this true. That orientation matters for tata_box_match() below,
## since the TATA consensus is not a palindrome, but not for
## poly_at_tract(), since an A/T run has the same length read from
## either strand.
extract_promoters <- function(genome, genes, max_bp = 300, min_bp = 50) {
  out      <- setNames(vector("character", nrow(genes)), genes$gene)
  extended <- setNames(vector("logical",   nrow(genes)), genes$gene)
  coords   <- data.frame(gene = genes$gene, seqid = NA_character_,
                          start = NA_integer_, end = NA_integer_,
                          strand = NA_character_, stringsAsFactors = FALSE)
  rownames(coords) <- genes$gene
  for (sq in unique(genes$seqid)) {
    g <- genes[genes$seqid == sq, ]
    g <- g[order(g$start), ]
    chrom <- genome[[sq]]
    for (i in seq_len(nrow(g))) {
      ext <- FALSE
      if (g$strand[i] == "+") {
        limit   <- if (i > 1) g$end[i - 1] + 1 else 1
        p_end   <- g$start[i] - 1
        p_start <- max(limit, p_end - max_bp + 1, 1)
        if (p_end - p_start + 1 < min_bp) {
          p_start_ext <- max(p_end - min_bp + 1, 1)
          if (p_start_ext < p_start) { p_start <- p_start_ext; ext <- TRUE }
        }
        seq <- if (p_start > p_end) "" else
          as.character(Biostrings::subseq(chrom, p_start, p_end))
      } else {
        limit   <- if (i < nrow(g)) g$start[i + 1] - 1 else length(chrom)
        p_start <- g$end[i] + 1
        p_end   <- min(limit, p_start + max_bp - 1, length(chrom))
        if (p_end - p_start + 1 < min_bp) {
          p_end_ext <- min(p_start + min_bp - 1, length(chrom))
          if (p_end_ext > p_end) { p_end <- p_end_ext; ext <- TRUE }
        }
        seq <- if (p_start > p_end) "" else
          as.character(Biostrings::reverseComplement(Biostrings::subseq(chrom, p_start, p_end)))
      }
      out[g$gene[i]]      <- toupper(seq)
      extended[g$gene[i]] <- ext && nchar(seq) > 0
      ## Recorded only when a real (nonempty) region was extracted;
      ## genes that stayed "" (boxed in on both sides) keep NA
      ## coordinates, the same signal has_seq uses elsewhere to skip
      ## them. start/end here are always the lower/higher genomic
      ## position regardless of strand, matching how p_start/p_end are
      ## already computed above; strand is what tells a downstream
      ## reader whether to read that span forward or to reverse it.
      if (nchar(seq) > 0)
        coords[g$gene[i], c("seqid", "start", "end", "strand")] <-
          list(sq, p_start, p_end, g$strand[i])
    }
  }
  attr(out, "extended") <- extended
  attr(out, "coords")   <- coords
  out
}

## Position frequency matrix for the eukaryotic TATA box (Bucher 1990,
## J Mol Biol 212:563-578), the eight TATAWAWR positions Basehoar's
## consensus also uses, but with real base percentages rather than a
## categorical call. Bucher derived this from 502 unrelated eukaryotic
## RNA polymerase II promoter sequences, so it is a real, quantitative
## matrix, not one measured in Saccharomyces specifically; worth
## flagging as a limitation if this analysis reaches the paper. Values
## are the published percentages divided by 100.
TATA.PWM <- rbind(
  A = c(4.1, 90.5,  0.8, 91.0, 68.9, 92.5, 57.1, 39.8),
  T = c(79.5, 9.0, 96.1,  7.7, 31.1,  1.6, 31.1,  8.5),
  G = c(4.6,  0.5,  0.5,  1.3,  0.0,  5.1, 11.3, 40.4),
  C = c(11.8, 0.0,  2.6,  0.0,  0.0,  0.8,  0.5, 11.3)) / 100

## Basehoar et al. located functional TATA boxes 50 to 200 bp upstream
## of the ATG start codon (their location criterion). TATA.PWM says
## what a TATA box looks like; this window says where to look for one,
## and the two come from separate papers for exactly that reason.
TATA.WINDOW <- c(50, 200)

## Scans the Basehoar location window of a promoter sequence and scores
## every 8-mer against TATA.PWM as a log2 odds ratio relative to a
## uniform 25% base composition, returning the best-scoring window as a
## continuous score rather than a binary hit. pseudocount replaces the
## matrix's exact-zero cells so a single mismatch at an otherwise
## strongly conserved position doesn't force the whole window's score
## to -Inf; 0.001 is a standard small-fraction choice, not a fitted
## value. Directional, since the matrix is not a palindrome, unlike the
## poly(dA:dT) search below; seq must already be oriented 5' to 3'
## relative to its own gene (see extract_promoters), with the
## promoter-proximal end at the end of the string. motif is the actual
## best-scoring 8-mer, kept so a candidate gene's real sequence is
## visible without separately pulling position and prom_len apart to
## find it by hand.
tata_box_score <- function(seq, window = TATA.WINDOW, pwm = TATA.PWM, pseudocount = 0.001) {
  motif_len <- ncol(pwm)
  n <- nchar(seq)
  idx_lo <- max(1, n - window[2] + 1)
  idx_hi <- n - window[1] - motif_len + 2
  if (idx_hi < idx_lo) return(list(score = NA, position = NA, motif = NA))
  starts <- idx_lo:idx_hi
  pwm_adj <- pmax(pwm, pseudocount)
  scores <- sapply(starts, function(i) {
    bases <- strsplit(substring(seq, i, i + motif_len - 1), "")[[1]]
    row_idx <- match(bases, rownames(pwm_adj))
    if (length(bases) < motif_len || any(is.na(row_idx))) return(NA)
    sum(log2(pwm_adj[cbind(row_idx, seq_len(motif_len))] / 0.25))
  })
  if (all(is.na(scores))) return(list(score = NA, position = NA, motif = NA))
  best <- which.max(scores)
  list(score = scores[best], position = starts[best], motif = substring(seq, starts[best], starts[best] + motif_len - 1))
}

## Finds the longest run of A or T bases anywhere in the sequence, the
## feature behind poly(dA:dT)-mediated nucleosome exclusion. Not
## direction-specific, so it is run on the full promoter as extracted,
## unlike tata_box_score(). An empty sequence (extract_promoters()
## returned "" because a neighboring gene left no room for a promoter)
## gives length NA, not 0: a real promoter of any usable length almost
## always contains at least one A or T base, so an actual measured
## length of 0 essentially never happens, and coding "no promoter to
## measure" as 0 would silently read as "no divergence" downstream,
## diluting any real signal in polyat_delta. tract is the actual
## longest A/T run, kept for the same reason motif is above.
poly_at_tract <- function(seq) {
  if (nchar(seq) == 0) return(list(length = NA, position = NA, tract = NA))
  runs <- gregexpr("[AT]+", seq)[[1]]
  lens <- attr(runs, "match.length")
  if (runs[1] == -1) return(list(length = 0, position = NA, tract = NA))
  best <- which.max(lens)
  list(length = lens[best], position = runs[best], tract = substring(seq, runs[best], runs[best] + lens[best] - 1))
}

## Applies both the TATA PWM score and the poly(dA:dT) tract to every
## sequence in a named vector and returns one row per gene, ready to
## merge across species. prom_len is the actual extracted promoter
## length, kept alongside the two scores so a short or empty promoter
## (capped by a nearby neighboring gene, see extract_promoters()) can be
## told apart from a full-length one that simply scored low. prom_
## extended carries through extract_promoters()'s "extended" attribute
## when present (FALSE for every gene otherwise, e.g. if seqs came from
## somewhere other than extract_promoters()), flagging genes whose
## promoter crossed into a neighboring gene's coding sequence to reach
## the min_bp floor; any score built from one of those regions should
## be read with that in mind; see extract_promoters() for why. tata_pos
## and polyat_pos report where each feature's best match starts, in bp
## upstream of the ATG (prom_len - raw index + 1) rather than as a raw
## string index, since raw indices aren't comparable between two
## promoters of different lengths the way a distance from the ATG is.
## Both functions already return that raw index as their "position"
## element; this just converts it. NA propagates automatically wherever
## tata_score or polyat_len is NA, since position is NA in the same
## cases. Useful for telling an indel (the same real element present in
## both species but shifted to a different distance from the ATG) apart
## from a mutation that changed which region scores best, once tata_pos
## or polyat_pos is compared alongside tata_delta or polyat_delta.
score_promoters <- function(seqs, tata_window = TATA.WINDOW) {
  tata     <- lapply(seqs, tata_box_score, window = tata_window)
  polyat   <- lapply(seqs, poly_at_tract)
  len      <- nchar(seqs)
  tata_idx   <- sapply(tata,   `[[`, "position")
  polyat_idx <- sapply(polyat, `[[`, "position")
  extended   <- attr(seqs, "extended")
  if (is.null(extended)) extended <- setNames(rep(FALSE, length(seqs)), names(seqs))
  data.frame(gene          = names(seqs),
             prom_len      = len,
             prom_extended = unname(extended[names(seqs)]),
             tata_score    = sapply(tata, `[[`, "score"),
             tata_pos      = len - tata_idx + 1,
             tata_motif    = sapply(tata, `[[`, "motif"),
             polyat_len    = sapply(polyat, `[[`, "length"),
             polyat_pos    = len - polyat_idx + 1,
             polyat_tract  = sapply(polyat, `[[`, "tract"),
             row.names  = NULL,
             stringsAsFactors = FALSE)
}

## Predicted nucleosome occupancy runs on the Linux cluster, the same way
## the bootstrap and permutation jobs do. Section 6.1 packages each
## species' chromosomes and promoter coordinates with
## nupop_cluster_inputs(), nupop_occupancy.R runs
## nupop_occupancy_cluster() on the cluster, and Section 6.1 scores the
## returned tracks with score_promoters_nupop().
##
## Each chromosome is scored as a series of core windows (window_bp,
## 250 kb by default), and each core is extended by `flank` bp of real
## genomic sequence on both sides before NuPoP sees it. Only the core is
## kept. NuPoP's own Fortran partitions long sequences the same way, with
## 7000 bp of overlap on each side of every partition, and a direct
## comparison on a 400 kb test sequence gave occupancy identical to the
## whole-sequence prediction at every position, to the three decimals
## NuPoP writes.
##
## Every window runs in its own forked process (mclapply with
## mc.preschedule = FALSE), so the job carries on when NuPoP cannot
## process a window. Such a window is split into two halves, each with
## its own flanks, and the halves run in the next round. The split
## repeats down to min_core_bp. A smallest window that still ends early
## runs once more with fallback_flank bp of context (2000 bp by default).
## On the same test sequence, 2000 bp flanks matched the whole-sequence
## prediction to within 0.001 at every position, while 1000 bp flanks
## differed by up to 0.034. The returned list carries two tables as
## attributes, "reduced_flank_regions" for windows scored with fallback
## flanks and "failed_regions" for windows left as NA.
##
## NuPoP keeps occupancy between 0 and 1 and writes -0.05 at every N
## position and at every non-N stretch shorter than one nucleosome
## (148 bp). Those positions are recorded as NA, so a promoter beside an
## assembly gap averages only real predictions.
##
## species = 7 selects NuPoP's S. cerevisiae profile, applied to both
## species so they share one reference. model = 4 selects the
## fourth-order Markov model.

## Packages one species for the cluster job. genome is a DNAStringSet of
## full chromosomes (e.g. GENOME.SC) and seqs is extract_promoters()'s
## return value, which carries the "coords" attribute. Chromosomes become
## plain character strings, so the cluster job needs NuPoP but not
## Biostrings. NuPoP's parser reads A, C, G, T and N, so any other IUPAC
## code in a draft assembly is written as N, the symbol NuPoP uses for an
## unknown base.
nupop_cluster_inputs <- function(genome, seqs) {
  coords <- attr(seqs, "coords")
  if (is.null(coords))
    stop("nupop_cluster_inputs() reads coords from the 'coords' attribute of extract_promoters()'s return value. Pass PROM.SC or PROM.SE directly as returned in Section 6.1.")
  coords <- coords[!is.na(coords$start), c("gene", "seqid", "start", "end", "strand")]
  chroms <- vapply(unique(coords$seqid), function(sq) {
    s <- toupper(as.character(genome[[sq]]))
    n_ambig <- nchar(gsub("[ACGTN]", "", s))
    if (n_ambig > 0) message(sprintf("  %s: %d IUPAC code(s) written as N", sq, n_ambig))
    gsub("[^ACGTN]", "N", s)
  }, character(1))
  list(chroms = chroms, coords = coords)
}

## Predicts occupancy for one sequence in the current process and returns
## positions core_from to core_to. Each call works in its own temporary
## folder, because NuPoP writes its prediction file into the working
## directory. The folder name carries the process ID, which is unique
## among running processes, so forked workers that start from the same
## tempfile() state still receive separate folders. The folders sit under
## work_root, which nupop_occupancy_cluster() places beside R's session
## temp directory rather than inside it, so each worker's files stay
## independent of every other worker.
nupop_predict_window <- function(seg_seq, core_from, core_to, species = 7, model = 4,
                                 work_root = tempdir()) {
  n <- nchar(seg_seq)
  ## NuPoP scores sequences of at least 148 bp, one nucleosome plus one base
  if (n < 148) return(rep(NA_real_, core_to - core_from + 1))
  wd <- tempfile(pattern = sprintf("nupop_pid%d_", Sys.getpid()), tmpdir = work_root)
  dir.create(wd, recursive = TRUE)
  old_wd <- setwd(wd)
  on.exit({ setwd(old_wd); unlink(wd, recursive = TRUE) }, add = TRUE)

  st <- seq(1, n, by = 80)
  writeLines(c(">window", substring(seg_seq, st, pmin(st + 79, n))), "window.fa")
  invisible(utils::capture.output(NuPoP::predNuPoP("window.fa", species = species, model = model)))
  pred <- paste0("window.fa_Prediction", model, ".txt")
  if (!file.exists(pred)) return(NULL)

  ## Columns are Position, P.start, Occup, N/L, Affinity. NULL skips a column.
  tab  <- scan(pred, what = list(0L, NULL, 0, NULL, NULL), skip = 1, quiet = TRUE)
  occ  <- rep(NA_real_, n)
  keep <- tab[[1]] >= 1 & tab[[1]] <= n & tab[[3]] >= 0
  occ[tab[[1]][keep]] <- tab[[3]][keep]
  occ[core_from:core_to]
}

## Runs a table of core windows (seqid, start, end) with `flank` bp of
## context, one forked process per window, and returns one result per
## row. A window whose process ends early returns NULL or a try-error,
## which the caller reads as "split and rerun".
nupop_run_windows <- function(chroms, tasks, flank, species, model, cores, work_root) {
  parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    s  <- chroms[[tasks$seqid[i]]]
    ws <- max(1, tasks$start[i] - flank)
    we <- min(nchar(s), tasks$end[i] + flank)
    nupop_predict_window(substr(s, ws, we), tasks$start[i] - ws + 1, tasks$end[i] - ws + 1,
                         species = species, model = model, work_root = work_root)
  }, mc.cores = cores, mc.preschedule = FALSE)
}

## Scores every chromosome in one species' inputs and returns each gene's
## promoter occupancy track. Rounds of windows run until every window is
## scored or reaches min_core_bp. Minus-strand slices are reversed so the
## promoter-proximal end sits at the end of every vector, matching
## tata_box_score() and poly_at_tract() on the same gene's sequence. The
## coords the tracks came from travel with the result, so
## score_promoters_nupop() can confirm they match the current promoters.
nupop_occupancy_cluster <- function(inputs, species = 7, model = 4,
                                    window_bp = 250000, flank = 7000, fallback_flank = 2000,
                                    min_core_bp = 5000, cores = 1) {
  chroms <- inputs$chroms
  coords <- inputs$coords
  occ    <- lapply(chroms, function(s) rep(NA_real_, nchar(s)))

  ## Node-local scratch space for every window's files, removed at the end
  work_root <- file.path(dirname(tempdir()), sprintf("nupop_work_%d", Sys.getpid()))
  dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_root, recursive = TRUE), add = TRUE)

  tasks <- do.call(rbind, lapply(names(chroms), function(sq) {
    n <- nchar(chroms[[sq]])
    s <- seq(1, n, by = window_bp)
    data.frame(seqid = sq, start = s, end = pmin(s + window_bp - 1, n), stringsAsFactors = FALSE)
  }))
  window_ok <- function(r, t) is.numeric(r) && length(r) == t$end - t$start + 1

  reduced <- NULL
  failed  <- NULL
  round   <- 0
  while (nrow(tasks) > 0) {
    round <- round + 1
    t0    <- Sys.time()
    res   <- nupop_run_windows(chroms, tasks, flank, species, model, cores, work_root)
    ok    <- vapply(seq_len(nrow(tasks)), function(i) window_ok(res[[i]], tasks[i, ]), logical(1))
    for (i in which(ok)) occ[[tasks$seqid[i]]][tasks$start[i]:tasks$end[i]] <- res[[i]]

    bad   <- tasks[!ok, , drop = FALSE]
    small <- bad[bad$end - bad$start + 1 <= min_core_bp, , drop = FALSE]
    big   <- bad[bad$end - bad$start + 1 >  min_core_bp, , drop = FALSE]

    if (nrow(small) > 0) {
      res2 <- nupop_run_windows(chroms, small, fallback_flank, species, model, cores, work_root)
      ok2  <- vapply(seq_len(nrow(small)), function(i) window_ok(res2[[i]], small[i, ]), logical(1))
      for (i in which(ok2)) occ[[small$seqid[i]]][small$start[i]:small$end[i]] <- res2[[i]]
      if (any(ok2))  reduced <- rbind(reduced, cbind(small[ok2, , drop = FALSE], flank = fallback_flank))
      if (any(!ok2)) failed  <- rbind(failed, small[!ok2, , drop = FALSE])
    }

    mid   <- (big$start + big$end) %/% 2
    tasks <- rbind(data.frame(seqid = big$seqid, start = big$start, end = mid,     stringsAsFactors = FALSE),
                   data.frame(seqid = big$seqid, start = mid + 1,   end = big$end, stringsAsFactors = FALSE))
    cat(sprintf("[%s] round %d: %d windows scored, %d split for the next round, %.1f min\n",
                format(Sys.time(), "%H:%M:%S"), round, sum(ok), nrow(big),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }

  out <- setNames(vector("list", nrow(coords)), coords$gene)
  for (i in seq_len(nrow(coords))) {
    s <- occ[[coords$seqid[i]]][coords$start[i]:coords$end[i]]
    out[[coords$gene[i]]] <- if (coords$strand[i] == "-") rev(s) else s
  }
  attr(out, "coords")                <- coords
  attr(out, "reduced_flank_regions") <- reduced
  attr(out, "failed_regions")        <- failed
  out
}
## Reduces each gene's occupancy track to the mean predicted occupancy
## over the Basehoar window used for the TATA score (50 to 200 bp
## upstream of the ATG), so all three architecture features read the
## same stretch of promoter. A window made entirely of NA positions
## returns NA.
nupop_window_score <- function(seqs, occ_list, window = TATA.WINDOW) {
  vapply(names(seqs), function(g) {
    n      <- nchar(seqs[[g]])
    occ    <- occ_list[[g]]
    idx_lo <- max(1, n - window[2] + 1)
    idx_hi <- min(n, n - window[1] + 1)
    if (n == 0 || idx_hi < idx_lo || length(occ) < n) return(NA_real_)
    vals <- occ[idx_lo:idx_hi]
    if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
  }, numeric(1))
}

## Returns one row per gene, the same shape score_promoters() returns,
## from the occupancy tracks the cluster job produced (NUPOP.OCC.SC or
## NUPOP.OCC.SE, loaded from nupop_output.rda). seqs must be
## extract_promoters()'s return value. The coordinate check confirms the
## tracks were computed from the same promoters as seqs, so a change to
## PROM.MAX.BP, PROM.MIN.BP or an annotation is caught here. Genes in an
## unscored region carry NA, and the region tables travel with the
## scores as attributes.
score_promoters_nupop <- function(seqs, occ_list, window = TATA.WINDOW) {
  coords <- attr(seqs, "coords")
  if (is.null(coords))
    stop("score_promoters_nupop() reads coords from the 'coords' attribute of extract_promoters()'s return value. Pass PROM.SC or PROM.SE directly as returned in Section 6.1.")
  coords <- coords[!is.na(coords$start), c("gene", "seqid", "start", "end", "strand")]
  if (!isTRUE(all.equal(coords, attr(occ_list, "coords"), check.attributes = FALSE)))
    stop("score_promoters_nupop(): these promoters differ from the ones the cluster job scored. Save nupop_inputs.rda again in Section 6.1 and rerun nupop_occupancy.R.")
  occ_score <- nupop_window_score(seqs, occ_list, window = window)
  scores    <- data.frame(gene = names(seqs), occ_score = unname(occ_score),
                          row.names = NULL, stringsAsFactors = FALSE)
  attr(scores, "reduced_flank_regions") <- attr(occ_list, "reduced_flank_regions")
  attr(scores, "failed_regions")        <- attr(occ_list, "failed_regions")
  scores
}

## Merges Sc and Se promoter scores on gene identity, keeping only genes
## scored in both species, and adds the between-species shift in each
## feature, the quantity a genome-wide comparison of promoter
## architecture actually needs rather than the per-species scores alone.
## tata_motif and polyat_tract (the actual best-matching sequences,
## carried through unchanged from score_promoters()) merge along with
## everything else, giving tata_motif_sc/_se and polyat_tract_sc/_se
## without any extra handling here.
## tata_pos_delta and polyat_pos_delta give the same shift for where
## each feature's best match sits (bp upstream of the ATG); a large
## score or length shift paired with a small position shift points to a
## mutation at roughly the same site, while a large position shift
## alongside a small score or length shift is more consistent with an
## indel moving an otherwise similar element, or with the match jumping
## to a different site entirely.
##
## occ_delta is the raw shift in predicted occupancy (Se minus Sc), kept
## for anyone who wants the literal NuPoP quantity. occ_access_delta
## flips its sign, since occupancy runs the opposite direction from the
## other two features: a higher TATA score or a longer poly(dA:dT)
## tract both mean a more open, more accessible promoter, while a
## higher occupancy score means a more nucleosome-covered, less
## accessible one. Every downstream test in this section that compares
## a feature's direction against noise divergence (6.2's validation,
## 6.6's concordance, 6.7's candidates) uses occ_access_delta so that
## "a positive value means more accessible" holds for all three
## features at once, without needing a separate sign convention for
## occupancy everywhere it appears.
promoter_divergence <- function(sc_scores, se_scores) {
  m <- merge(sc_scores, se_scores, by = "gene", suffixes = c("_sc", "_se"))
  m$tata_delta        <- m$tata_score_se - m$tata_score_sc
  m$tata_pos_delta    <- m$tata_pos_se   - m$tata_pos_sc
  m$polyat_delta      <- m$polyat_len_se - m$polyat_len_sc
  m$polyat_pos_delta  <- m$polyat_pos_se - m$polyat_pos_sc
  m$occ_delta         <- m$occ_score_se  - m$occ_score_sc
  m$occ_access_delta  <- -m$occ_delta
  m
}

## Tests whether the direction of a promoter feature's between-species
## shift matches the direction NOISE.VALIDATE (Section 6.2) predicts:
## higher TATA score or longer poly(dA:dT) associates with lower DISP
## (burstier, noisier expression) in both species on their own, so the
## species that gained more of a feature should be the noisier allele
## in the hybrid. bfreq_cis_est = log2(DISP_Sc) - log2(DISP_Se), so a
## positive value means Sc has the higher DISP and Se is the noisier
## allele; that predicts Se should also carry the larger feature value,
## i.e. delta (Se - Sc) should be positive too. Concordant is
## sign(delta) == sign(bfreq_cis_est); the two should move together if
## the feature is genuinely acting on noise through the mechanism
## NOISE.VALIDATE found.
##
## Run on the full any-cis gene set (cis_classes, not narrowed by an
## arch_frac magnitude cutoff the way promoter_noise_candidates() is),
## since this is asking whether there is a real direction-of-effect
## signal to begin with, not looking for individual genes worth
## inspecting; restricting to the most extreme deltas first would
## condition on magnitude before testing direction and only cost power.
## A binomial sign test against 50/50, one-sided for an excess of
## concordant genes specifically (alternative = "greater"), is the
## right test: under the null that promoter architecture has nothing to
## do with which allele is noisier, concordant and discordant genes
## should be equally common regardless of how many genes are tested. A
## non-significant result here means the candidate list that follows in
## Section 6.7 should be read with real skepticism, since it would say
## the promoter/noise relationship established per-species in 6.2 isn't
## showing up in the direction it should for the cis divergence itself.
promoter_direction_test <- function(BURST.CONTRASTS, PR, ARCH, reg_class, quantity = c("bfreq", "bsize", "kbal"), cis_classes = c("Cis", "Cis + Trans", "Compensatory")) {
  quantity <- match.arg(quantity)
  est <- BURST.CONTRASTS[[paste0(quantity, "_cis_est")]]
  cis_flag <- reg_class %in% cis_classes
  one_feature <- function(delta) {
    ok <- cis_flag & !is.na(delta) & !is.na(est) &
          sign(delta) != 0 & sign(est) != 0
    n          <- sum(ok)
    concordant <- sign(delta[ok]) == sign(est[ok])
    n_conc     <- sum(concordant)
    bt <- if (n > 0) binom.test(n_conc, n, p = 0.5, alternative = "greater")
          else list(p.value = NA_real_)
    data.frame(n = n, n_concordant = n_conc,
               frac_concordant = if (n > 0) n_conc / n else NA_real_,
               p = bt$p.value)
  }
  rbind(cbind(feature = "TATA", one_feature(ARCH$tata_delta)), cbind(feature = "poly(dA:dT)", one_feature(ARCH$polyat_delta)), cbind(feature = "nucleosome occupancy (NuPoP)", one_feature(ARCH$occ_access_delta)))
}

## Splits the same any-cis gene set promoter_direction_test() tests into
## n_bins equal-count bins by |delta| magnitude and reports the fraction
## concordant with bfreq_cis_est's sign in each bin. A real effect that
## is only detectable in the most divergent promoters would be diluted
## to nothing by the many weakly divergent genes in a whole-set test
## (promoter_direction_test()'s result), but should still show up here
## as a rising concordance rate from the smallest-|delta| bin to the
## largest; if it doesn't, that is further evidence there is no signal
## to find at any magnitude, not just that the whole-set test lacked
## power. Eligibility (cis_classes membership, both values defined,
## neither sign exactly 0) is rebuilt here rather than shared with
## promoter_direction_test(), a small duplication that keeps this
## function usable on its own without depending on internal state from
## another call.
##
## ci_lo/ci_hi are a normal-approximation 95% interval on each bin's
## proportion, adequate for the bin sizes this produces and meant only
## to show roughly how much a bin's estimate should be trusted, not as
## a formal per-bin test. trend_coef and trend_p (attached as
## attributes, since they describe the whole feature rather than any
## one bin) come from a logistic regression of concordance on
## log(|delta|), the formal version of "is the rate rising toward the
## extremes" that the binned table and its plot show informally.
concordance_by_magnitude <- function(delta, bfreq_cis_est, reg_class, cis_classes = c("Cis", "Cis + Trans", "Compensatory"), n_bins = 10) {
  cis_flag <- reg_class %in% cis_classes
  ok <- cis_flag & !is.na(delta) & !is.na(bfreq_cis_est) &
        sign(delta) != 0 & sign(bfreq_cis_est) != 0

  d    <- abs(delta[ok])
  conc <- sign(delta[ok]) == sign(bfreq_cis_est[ok])
  ## Rank-based equal-count binning rather than cut(quantile(...)):
  ## poly(dA:dT) length differences take only a few small integer
  ## values, so quantile() breakpoints collide and cut() fails on
  ## non-unique breaks once ties outnumber the requested bins. Ranking
  ## first and dividing by n avoids that regardless of how many ties
  ## are in the data; ties.method = "first" only decides which side of
  ## a bin boundary a tied value falls on, it doesn't bias which bin
  ## counts as "extreme".
  bin <- ceiling(rank(d, ties.method = "first") / length(d) * n_bins)
  bin <- factor(pmin(bin, n_bins), levels = seq_len(n_bins))

  agg <- data.frame(
    bin          = seq_len(n_bins),
    n            = as.integer(tapply(conc, bin, length)),
    n_concordant = as.integer(tapply(conc, bin, sum)),
    mag_mean     = as.numeric(tapply(d, bin, mean)))
  agg$frac_concordant <- agg$n_concordant / agg$n
  se <- sqrt(agg$frac_concordant * (1 - agg$frac_concordant) / agg$n)
  agg$ci_lo <- pmax(0, agg$frac_concordant - 1.96 * se)
  agg$ci_hi <- pmin(1, agg$frac_concordant + 1.96 * se)

  fit <- glm(conc ~ log(d), family = binomial)
  attr(agg, "trend_coef") <- unname(coef(fit)[2])
  attr(agg, "trend_p")    <- unname(summary(fit)$coefficients[2, 4])
  agg
}

## Bar plot of concordance_by_magnitude()'s output: one bar per bin,
## ordered smallest to largest |delta|, with the 95% interval as an
## error bar and a dashed reference line at the 50/50 null. Bars
## labeled with each bin's mean |delta| rather than a bin index, so the
## x-axis reads directly as a magnitude scale.
plot_concordance_by_magnitude <- function(cb, main = NULL) {
  bp <- barplot(cb$frac_concordant, ylim = c(0, 1), col = "grey75", names.arg = round(cb$mag_mean, 1), las = 2, ylab = "fraction concordant", xlab = "mean |delta| in bin", main = main)
  segments(bp, cb$ci_lo, bp, cb$ci_hi)
  abline(h = 0.5, lty = 2, col = "#D62728")
}

## Flags genes worth pulling up by eye: a significant cis component of
## noise divergence, paired with a large shift in a promoter feature.
## Neither condition alone is enough, a promoter difference with no cis
## noise signature isn't evidence of anything, and cis noise divergence
## with no promoter difference isn't explained by anything measured in
## this section.
##
## Returns three independent tables rather than one combined one, since
## a gene only needs a large shift in one of the three features, not
## all of them, to be a candidate, and keeping every feature's columns
## on every row regardless of which one actually flagged the gene made
## the table wide for no reason. A gene can appear in more than one
## table if it clears more than one cutoff; that is expected, not a
## duplicate to remove.
##
## Gates on regulatory class rather than the total noise contrast.
## Cis, Cis + Trans, and Compensatory all require a significant cis
## permutation p-value by construction (classify_reg(), Section 3.1),
## so requiring reg_class to be one of those three is a real, tested
## cis effect, not just a large total point estimate; Trans-only and
## Conserved genes are excluded because a local promoter difference on
## this gene isn't the expected explanation for either, and a
## Conserved call paired with a large total effect size in particular
## usually means the total estimate isn't trustworthy (large SE, not a
## real large effect that failed to show up in either component test).
## Ambiguous genes fall out of the same filter automatically, since
## "Ambiguous" is not one of the three classes being kept.
##
## There is no null distribution for a single-genome PWM score or
## tract-length difference the way there is for the bootstrap and
## permutation-based noise contrasts, so "large" here is defined by
## rank, the top arch_frac of |delta| for that feature, rather than a
## p-value; each feature's cutoff is computed on its own values, so the
## two tables are not required to end up the same size
##
## bfreq_cis_p ties at a shared minimum value for many genes when this
## runs on the provisional N.PERM = 500 pass (Section 2), since a
## permutation p-value can't resolve past 1/(N.PERM+1); that is a
## resolution limit of the permutation count, not a problem with these
## particular genes, and will spread out once the final N.PERM = 10000
## pass replaces it. It does not affect the ranking below, which sorts
## on effect size rather than p-value.
##
## BURST.CONTRASTS, PR, ARCH, and reg_class must already share the same gene order
## (true for BURST.CONTRASTS, PR, and REG.BFREQ.CLASS by construction, and true for
## ARCH once it has been realigned to BURST.CONTRASTS$gene as Section 6.3 does).
## Each table is sorted by |bfreq_cis_est| rather than the total effect:
## cis compares the two alleles within the hybrid, sharing one trans
## environment, so it is the allele-specific quantity a promoter
## difference should actually track, while total compares separately
## grown parental strains and can carry confounds a cis contrast
## doesn't. cis, trans, and total effect/p-value are all kept in both
## tables, since reconciling why a gene landed in a given reg_class
## needs cis and trans side by side, and total is kept alongside cis
## for a quick magnitude comparison even though it isn't what the
## ranking uses. concordant flags whether this gene's own delta sign
## matches the direction promoter_direction_test() checks for
## enrichment across the whole any-cis gene set; a discordant candidate
## isn't necessarily wrong (a single gene can go against an overall
## trend), but it is a reason for extra scrutiny before treating it as
## a clean example. Within each table, both species' raw values (score or
## length, position, and the actual matched sequence) come before that
## feature's own deltas, so an indel can be told apart from a
## substitution by pulling the actual sequences directly out of the
## table, the way the YML001W and YOR162C-style checks earlier in this
## section were done by hand.
promoter_noise_candidates <- function(BURST.CONTRASTS, PR, ARCH, reg_class, quantity = c("bfreq", "bsize", "kbal"), cis_classes = c("Cis", "Cis + Trans", "Compensatory"), arch_frac = 0.90) {
  quantity <- match.arg(quantity)
  stopifnot(nrow(BURST.CONTRASTS) == nrow(ARCH), nrow(BURST.CONTRASTS) == length(reg_class))

  required_cols <- c("tata_score_sc", "tata_score_se", "tata_pos_sc", "tata_pos_se",
                      "tata_motif_sc", "tata_motif_se", "tata_delta", "tata_pos_delta",
                      "polyat_len_sc", "polyat_len_se", "polyat_pos_sc", "polyat_pos_se",
                      "polyat_tract_sc", "polyat_tract_se", "polyat_delta", "polyat_pos_delta",
                      "occ_score_sc", "occ_score_se", "occ_delta", "occ_access_delta",
                      "prom_len_sc", "prom_len_se", "prom_extended_sc", "prom_extended_se")
  missing_cols <- setdiff(required_cols, names(ARCH))
  if (length(missing_cols) > 0)
  stop(sprintf("ARCH is missing columns: %s. This ARCH was built with an older version of score_promoters()/promoter_divergence(); rerun Section 6.1 through 6.3 to regenerate SCORE.SC, SCORE.SE, and ARCH with the current column set before calling promoter_noise_candidates().", paste(missing_cols, collapse = ", ")))

  est_cis   <- BURST.CONTRASTS[[paste0(quantity, "_cis_est")]]
  est_trans <- BURST.CONTRASTS[[paste0(quantity, "_trans_est")]]
  est_total <- BURST.CONTRASTS[[paste0(quantity, "_total_est")]]
  p_cis     <- PR[[paste0(quantity, "_cis_p")]]
  p_trans   <- PR[[paste0(quantity, "_trans_p")]]
  p_total   <- PR[[paste0(quantity, "_total_p")]]

  cis_flag <- reg_class %in% cis_classes

  tata_cut    <- quantile(abs(ARCH$tata_delta),        arch_frac, na.rm = TRUE)
  polyat_cut  <- quantile(abs(ARCH$polyat_delta),      arch_frac, na.rm = TRUE)
  occ_cut     <- quantile(abs(ARCH$occ_access_delta),  arch_frac, na.rm = TRUE)
  tata_flag   <- !is.na(ARCH$tata_delta)         & abs(ARCH$tata_delta)        >= tata_cut
  polyat_flag <- !is.na(ARCH$polyat_delta)       & abs(ARCH$polyat_delta)      >= polyat_cut
  occ_flag    <- !is.na(ARCH$occ_access_delta)   & abs(ARCH$occ_access_delta)  >= occ_cut

  keep_tata   <- cis_flag & tata_flag;   keep_tata[is.na(keep_tata)]     <- FALSE
  keep_polyat <- cis_flag & polyat_flag; keep_polyat[is.na(keep_polyat)] <- FALSE
  keep_occ    <- cis_flag & occ_flag;    keep_occ[is.na(keep_occ)]       <- FALSE

  ## Shared columns every candidate needs regardless of which feature
  ## flagged it, kept as one helper so the two tables build identically.
  ## Column names carry the quantity prefix (bfreq_ or bsize_) so a
  ## reader can tell which axis a candidate table was built from
  ## without checking the call site.
  context_cols <- function(keep) {
    out <- data.frame(
      gene      = BURST.CONTRASTS$gene[keep],
      reg_class = reg_class[keep],
      stringsAsFactors = FALSE)
    out[[paste0(quantity, "_cis_est")]]   <- est_cis[keep];   out[[paste0(quantity, "_cis_p")]]   <- p_cis[keep]
    out[[paste0(quantity, "_trans_est")]] <- est_trans[keep]; out[[paste0(quantity, "_trans_p")]] <- p_trans[keep]
    out[[paste0(quantity, "_total_est")]] <- est_total[keep]; out[[paste0(quantity, "_total_p")]] <- p_total[keep]
    out
  }
  est_col <- paste0(quantity, "_cis_est")

  tata_out <- cbind(context_cols(keep_tata), data.frame(
    tata_score_sc  = ARCH$tata_score_sc[keep_tata], tata_score_se  = ARCH$tata_score_se[keep_tata],
    tata_pos_sc    = ARCH$tata_pos_sc[keep_tata],   tata_pos_se    = ARCH$tata_pos_se[keep_tata],
    tata_motif_sc  = ARCH$tata_motif_sc[keep_tata], tata_motif_se  = ARCH$tata_motif_se[keep_tata],
    tata_delta     = ARCH$tata_delta[keep_tata],    tata_pos_delta = ARCH$tata_pos_delta[keep_tata],
    concordant     = sign(ARCH$tata_delta[keep_tata]) == sign(est_cis[keep_tata]),
    prom_len_sc    = ARCH$prom_len_sc[keep_tata],   prom_len_se    = ARCH$prom_len_se[keep_tata],
    prom_extended_sc = ARCH$prom_extended_sc[keep_tata], prom_extended_se = ARCH$prom_extended_se[keep_tata],
    row.names = NULL, stringsAsFactors = FALSE))
  tata_out <- tata_out[order(-abs(tata_out[[est_col]])), ]

  polyat_out <- cbind(context_cols(keep_polyat), data.frame(
    polyat_len_sc    = ARCH$polyat_len_sc[keep_polyat],     polyat_len_se    = ARCH$polyat_len_se[keep_polyat],
    polyat_pos_sc    = ARCH$polyat_pos_sc[keep_polyat],     polyat_pos_se    = ARCH$polyat_pos_se[keep_polyat],
    polyat_tract_sc  = ARCH$polyat_tract_sc[keep_polyat],   polyat_tract_se  = ARCH$polyat_tract_se[keep_polyat],
    polyat_delta     = ARCH$polyat_delta[keep_polyat],      polyat_pos_delta = ARCH$polyat_pos_delta[keep_polyat],
    concordant       = sign(ARCH$polyat_delta[keep_polyat]) == sign(est_cis[keep_polyat]),
    prom_len_sc      = ARCH$prom_len_sc[keep_polyat],       prom_len_se      = ARCH$prom_len_se[keep_polyat],
    prom_extended_sc = ARCH$prom_extended_sc[keep_polyat],  prom_extended_se = ARCH$prom_extended_se[keep_polyat],
    row.names = NULL, stringsAsFactors = FALSE))
  polyat_out <- polyat_out[order(-abs(polyat_out[[est_col]])), ]

  ## occ_score_sc/_se are raw predicted occupancy (NuPoP's own scale, a
  ## probability between 0 and 1 that a base is nucleosome-covered),
  ## kept alongside occ_access_delta so a candidate can be checked
  ## against the actual occupancy track directly. concordant still uses
  ## occ_access_delta rather than occ_delta, matching the sign
  ## convention used everywhere else in this section.
  occ_out <- cbind(context_cols(keep_occ), data.frame(
    occ_score_sc      = ARCH$occ_score_sc[keep_occ],      occ_score_se     = ARCH$occ_score_se[keep_occ],
    occ_delta          = ARCH$occ_delta[keep_occ],         occ_access_delta = ARCH$occ_access_delta[keep_occ],
    concordant         = sign(ARCH$occ_access_delta[keep_occ]) == sign(est_cis[keep_occ]),
    prom_len_sc        = ARCH$prom_len_sc[keep_occ],       prom_len_se      = ARCH$prom_len_se[keep_occ],
    prom_extended_sc   = ARCH$prom_extended_sc[keep_occ],  prom_extended_se = ARCH$prom_extended_se[keep_occ],
    row.names = NULL, stringsAsFactors = FALSE))
  occ_out <- occ_out[order(-abs(occ_out[[est_col]])), ]

  list(tata = tata_out, polyat = polyat_out, occ = occ_out)
}

## Tests whether a promoter feature associates with a species' own
## absolute noise level, independent of any cross-species comparison.
## fit is one of SPLIT.FITS' per-species data frames (Section 2);
## response names the column to test, DISP (burst frequency, the
## default and the axis used throughout the rest of this section) or
## BSIZE (MU/DISP, already computed on every fit object by
## fit_counts_offset_row()); score is a named vector of a single
## feature (tata_score or polyat_len) for the same species, keyed the
## same way as fit's row names.
##
## covariate names the column held fixed while testing the feature,
## MU by default. This matters specifically for response = "BSIZE":
## since BSIZE = MU/DISP, testing BSIZE against covariate MU
## is algebraically forced to reproduce the DISP-against-MU result with
## the coefficient sign flipped, log(BSIZE) = log(MU) - log(DISP),
## so conditioning on log(MU) leaves only -log(DISP) for the feature to
## explain, the same one degree of freedom as before under a different
## name, not an independent test. Conditioning on DISP instead asks the
## question that is actually different from the frequency test: among
## genes that fire at the same rate, does the feature predict how much
## comes out per burst. That is not degenerate with the DISP-on-MU test,
## since DISP and MU are two distinct, only partially correlated
## quantities across genes.
##
## For response = "DISP" (covariate = "MU"): lower DISP means burstier,
## noisier expression, so the published relationship (TATA-containing
## and long poly(dA:dT) promoters are noisier: Newman 2006, Blake 2006,
## Hornung 2012, Carey 2013) predicts a negative coefficient on the
## feature once MU is held fixed. For response = "BSIZE"
## (covariate = "DISP"), no direction is asserted here; the molecular
## literature motivating a promoter-architecture effect on burst
## frequency doesn't make an equally clear prediction for burst
## amplitude, so this is run as an open check rather than one with an
## expected sign to confirm.
##
## r2_partial is the R-squared gained by adding the feature to a model
## that already has log(covariate) (full model R-squared minus the
## covariate-only model's), the variance in log(response) the feature
## explains beyond what the covariate alone accounts for. p alone can't
## distinguish a real, large effect from a real, tiny one once n is in
## the thousands, since even a tiny slope clears significance easily at
## that sample size; r2_partial is the number that actually answers
## "how much of the variance is this."
architecture_noise_check <- function(fit, score, feature, response = "DISP", covariate = "MU") {
  ## Confirms the two fit columns and the feature score are numeric before
  ## modelling, and names the column and its class when one is not
  for (col in c(response, covariate)) {
    v <- if (col %in% colnames(fit)) fit[, col] else NULL
    if (!is.numeric(v))
      stop(sprintf("architecture_noise_check(%s): fit column '%s' must be numeric, and it is %s. Fit columns: %s",
                   feature, col, if (is.null(v)) "absent" else paste(class(v), collapse = "/"),
                   paste(colnames(fit), collapse = ", ")), call. = FALSE)
  }
  if (!is.numeric(score) || is.null(names(score)))
    stop(sprintf("architecture_noise_check(%s): score must be a numeric vector named by gene, and it is %s",
                 feature, paste(class(score), collapse = "/")), call. = FALSE)
  g <- intersect(rownames(fit), names(score))
  y <- fit[g, response]; z <- fit[g, covariate]; x <- score[g]
  ok  <- is.finite(y) & y > 0 & is.finite(z) & z > 0 & is.finite(x)
  m_full    <- lm(log(y[ok]) ~ x[ok] + log(z[ok]))
  m_reduced <- lm(log(y[ok]) ~ log(z[ok]))
  s <- summary(m_full)$coefficients
  data.frame(feature    = feature, response = response, covariate = covariate, n = sum(ok),
             estimate   = unname(s["x[ok]", "Estimate"]),
             p          = unname(s["x[ok]", "Pr(>|t|)"]),
             r2_partial = summary(m_full)$r.squared - summary(m_reduced)$r.squared)
}

##############################################################################
## 11. SEURAT CLUSTERING DIAGNOSTICS                                       ##
##############################################################################
## Elbow point on the ranked standardized-variance curve from
## FindVariableFeatures(method = "vst"), the same style as the PC elbow
## used for dimensionality selection: the latest point where the
## rank-to-rank drop exceeds drop_frac of the curve's own range, plus
## one. obj is re-run with nfeatures set to the full gene count first, so
## the standardized-variance ranking the elbow search sees is never
## truncated by a prior, arbitrary nfeatures choice. Returns both the
## chosen count and the full sorted curve, so the curve can be plotted
## with the chosen cutoff marked, the same check the PC elbow plot
## provides.
##
## The search is restricted to genes with standardized variance >=
## floor (default 1, the fitted trend line itself) before the elbow
## logic runs. Below that a gene is UNDER-dispersed relative to the
## mean-variance trend, not a candidate to be a top variable gene at
## all, and a handful of very low-count genes with a poorly fit trend
## can produce a second, sharp drop at the extreme low-variance tail
## that has nothing to do with the real elbow: a rank-to-rank threshold
## has no way to tell that cliff apart from the genuine one, and being
## further out in rank and often steeper, it wins. Filtering to floor
## first removes that artifact at its source. curve returned is the
## full, unfiltered ranking, so the tail is still visible on the
## diagnostic plot even though the search itself ignores it.
hvg_elbow <- function(obj, drop_frac = 0.001, floor = 1) {
  obj   <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj))
  v_all <- sort(HVFInfo(obj, method = "vst")$variance.standardized, decreasing = TRUE)
  v     <- v_all[v_all >= floor]
  if (length(v) < 2) v <- v_all   # fallback if nothing clears the floor
  rng <- diff(range(v, na.rm = TRUE))
  d   <- v[-length(v)] - v[-1]
  hit <- which(d > drop_frac * rng)
  n   <- if (length(hit) == 0) length(v) else max(hit) + 1
  list(n_features = n, curve = v_all, floor = floor)
}

## Diagnostic plot for hvg_elbow()'s output: full standardized-variance
## curve (grey past the cutoff), chosen cutoff, and floor marked. label
## names the output file (extra/S_hvg_elbow_<label>.pdf) and plot title.
plot_hvg_elbow <- function(hvg, label, fig_dir) {
  pdf(file.path(fig_dir, sprintf("extra/S_hvg_elbow_%s.pdf", label)), width = 6, height = 5, useDingbats = FALSE)
  plot(hvg$curve, pch = 16, cex = 0.4, col = ifelse(seq_along(hvg$curve) <= hvg$n_features, "black", "grey75"),
       xlab = "gene rank (by standardized variance)", ylab = "standardized variance", main = sprintf("%s: HVG elbow at %d genes", label, hvg$n_features))
  abline(v = hvg$n_features, lty = 2, col = "red")
  abline(h = hvg$floor, lty = 3, col = "blue")
  legend("topright", legend = c("chosen cutoff", sprintf("floor (%.1f)", hvg$floor)), lty = c(2, 3), col = c("red", "blue"), bty = "n")
  dev.off()
}

## Resolution sweep with an oversplit guard. At each resolution in
## res_grid, clusters obj on the given PCs, discards the resolution
## outright if any resulting cluster falls below min_cells (the power
## floor is enforced before silhouette is ever consulted, not after),
## and computes mean silhouette width (cluster::silhouette) on the same
## PCA-space distance used for clustering for every resolution that
## survives. Rather than the resolution that maximizes silhouette, this
## returns the COARSEST (smallest) resolution whose silhouette is within
## tol of the observed maximum: Louvain silhouette curves are often flat
## near their top with several local maxima, since splitting a real
## cluster into two arbitrary halves has no gap to find and so tends to
## lower silhouette rather than raise it, so taking the coarsest
## near-maximal point only lets resolution increase when doing so buys a
## real, non-marginal improvement in separation.
## Returns the chosen resolution, its silhouette and cluster count, the
## fitted Seurat object at that resolution (for downstream use without
## re-clustering), and the full per-resolution grid for plotting.
sweep_cluster_resolution <- function(obj, dims, res_grid = seq(0.05, 1, by = 0.05), min_cells = 50, tol = 0.02, metric = "manhattan") {
  emb <- Embeddings(obj, "pca")[, dims, drop = FALSE]
  d   <- dist(emb, method = metric)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims, annoy.metric = metric, verbose = FALSE)

  grid <- lapply(res_grid, function(r) {
    o2    <- FindClusters(obj, resolution = r, verbose = FALSE)
    cl    <- Idents(o2)
    sizes <- table(cl)
    # A resolution coarse enough to merge everything into one cluster
    # passes the min_cells check trivially (one big cluster is well
    # above the floor) but leaves silhouette width undefined, since
    # there is no other cluster for any point to be compared against.
    # That case needs its own guard, not just the size floor.
    if (length(sizes) < 2 || min(sizes) < min_cells) return(data.frame(res = r, ok = FALSE, sil = NA_real_, n_clusters = length(sizes), min_size = min(sizes)))
    sil <- mean(silhouette(as.integer(cl), d)[, 3])
    data.frame(res = r, ok = TRUE, sil = sil, n_clusters = length(sizes), min_size = min(sizes))
  })
  grid <- do.call(rbind, grid)

  ok_grid <- grid[grid$ok, ]
  if (nrow(ok_grid) == 0) stop("No resolution in res_grid keeps every cluster >= min_cells; widen res_grid or lower min_cells")
  best_sil  <- max(ok_grid$sil)
  near_max  <- ok_grid[ok_grid$sil >= best_sil - tol, ]
  chosen    <- near_max[which.min(near_max$res), ]

  obj <- FindClusters(obj, resolution = chosen$res, verbose = FALSE)
  list(chosen_res = chosen$res, chosen_sil = chosen$sil, chosen_n_clusters = chosen$n_clusters, obj = obj, grid = grid)
}

## Diagnostic plot for sweep_cluster_resolution()'s output: silhouette
## vs. resolution for every resolution that cleared the min-cluster-size
## guard, resolutions excluded by that guard marked separately along the
## bottom, and the chosen resolution marked. label names the output file
## (extra/S_resolution_sweep_<label>.pdf) and plot title.
plot_resolution_sweep <- function(sweep, label, fig_dir, min_cells) {
  g  <- sweep$grid
  ok <- g$ok
  pdf(file.path(fig_dir, sprintf("extra/S_resolution_sweep_%s.pdf", label)), width = 6, height = 5, useDingbats = FALSE)
  plot(g$res[ok], g$sil[ok], type = "b", pch = 16, xlab = "resolution", ylab = "mean silhouette width",
       main = sprintf("%s: resolution sweep", label), ylim = range(g$sil[ok], na.rm = TRUE))
  if (any(!ok)) points(g$res[!ok], rep(min(g$sil[ok], na.rm = TRUE), sum(!ok)), pch = 4, col = "grey60")
  abline(v = sweep$chosen_res, lty = 2, col = "red")
  legend("bottomright", legend = c("silhouette", sprintf("below %d cells/cluster", min_cells), "chosen"), pch = c(16, 4, NA), lty = c(NA, NA, 2), col = c("black", "grey60", "red"), bty = "n")
  dev.off()
}

## Bootstrap stability of a chosen clustering: resamples cells with
## replacement from the raw count matrix (not by subsetting the fitted
## Seurat object, which requires unique cell names and cannot represent
## a cell drawn twice), rebuilds a fresh Seurat object with uniquified
## barcodes for the resampled columns, and reruns the full
## Normalize/FindVariableFeatures/Scale/PCA/Neighbors/Clusters pipeline
## at the SAME nfeatures/dims/resolution/metric as the original fit, so
## the comparison isolates sampling noise rather than re-deciding any of
## those choices per replicate. Compares the resampled clustering to the
## original clustering restricted to the same resampled cells (in the
## same draw order) via adjusted Rand index. High mean ARI means the
## partition is not an artifact of Louvain being forced to draw a
## boundary through what is really continuous variation; a low or
## widely-spread ARI is itself informative about where the practical
## resolution limit sits for this dataset's cell count.
##
## Gene names are converted from underscore to dash and the matrix
## coerced to CsparseMatrix once, up front, the same fix used everywhere
## else in the pipeline before CreateSeuratObject() (see to_seurat_counts()
## in the driver): without it, CreateSeuratObject() silently does both
## conversions itself but with a warning every time, so B calls means B
## repeats of the same two warnings.
## make_boot_idx: draws every bootstrap resample up front, one column
## per replicate, from a single seeded stream. Fixing the draws before
## any Seurat call runs keeps each replicate an independent resample,
## because RunPCA() reseeds the global RNG (seed.use = 42) inside every
## replicate. The same matrix serves every candidate resolution of a
## dataset, so the resolution comparison is paired on identical draws.
make_boot_idx <- function(n, B, seed = 1) {
  set.seed(seed)
  matrix(replicate(B, sample.int(n, n, replace = TRUE)), nrow = n)
}

## boot_ari_one: one bootstrap replicate. Rebuilds the full
## Normalize/HVG/Scale/PCA/Neighbors/Clusters pipeline on the resampled
## cells at the original settings and scores agreement with the
## reference partition on those same cells by adjusted Rand index.
## Self-contained, so it runs equally well serially or on a forked
## worker (cluster_stability.R).
boot_ari_one <- function(counts, idx, ref_clusters, nfeatures, dims_n, resolution, metric = "manhattan") {
  boot_counts <- counts[, idx]
  rownames(boot_counts) <- gsub("_", "-", rownames(counts), fixed = TRUE)
  colnames(boot_counts) <- make.unique(colnames(counts)[idx])
  boot_counts <- as(boot_counts, "CsparseMatrix")
  boot_obj <- CreateSeuratObject(counts = boot_counts)
  boot_obj <- NormalizeData(boot_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  boot_obj <- FindVariableFeatures(boot_obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  boot_obj <- ScaleData(boot_obj, features = rownames(boot_obj), verbose = FALSE)
  boot_obj <- RunPCA(boot_obj, features = VariableFeatures(boot_obj), verbose = FALSE)
  boot_obj <- FindNeighbors(boot_obj, reduction = "pca", dims = 1:dims_n, annoy.metric = metric, verbose = FALSE)
  boot_obj <- FindClusters(boot_obj, resolution = resolution, verbose = FALSE)
  adjustedRandIndex(as.integer(ref_clusters[idx]), as.integer(Idents(boot_obj)))
}

## Serial wrapper around boot_ari_one(), kept as the local path. Draws
## come from make_boot_idx() unless a matrix is supplied.
bootstrap_cluster_stability <- function(counts, clusters, nfeatures, dims_n, resolution, metric = "manhattan", B = 20, seed = 1,
                                        idx = make_boot_idx(ncol(counts), B, seed)) {
  ari <- vapply(seq_len(ncol(idx)), function(b)
    boot_ari_one(counts, idx[, b], clusters, nfeatures, dims_n, resolution, metric), numeric(1))
  list(mean_ari = mean(ari), ari = ari)
}

## Coarsest resolution within the contiguous plateau of ok that shares
## the same cluster count as the resolution with the single highest
## silhouette. Grouping by cluster count rather than by a silhouette
## value threshold is deliberate: two genuinely distinct partitions (say
## 3 clusters vs. 5) can differ in silhouette by less than a reasonable
## noise tolerance, so any fixed silhouette-based cutoff risks bridging
## a real regime change if the gap happens to be smaller than that
## tolerance, exactly what happened when a 0.02 tolerance failed to
## separate a 3-cluster partition at 0.05 from a 5-cluster partition
## only ~0.01 higher in silhouette. Cluster count is discrete and
## unambiguous where it changes, so it does not have this failure mode.
## Contiguous runs are found by cutting wherever n_clusters changes
## between adjacent (in resolution order) grid points; the run
## containing the grid's argmax is the target plateau, and its coarsest
## resolution is returned.
plateau_coarsest <- function(ok) {
  ok  <- ok[order(ok$res), ]
  grp <- cumsum(c(1, diff(ok$n_clusters) != 0))
  target_grp <- grp[which.max(ok$sil)]
  min(ok$res[grp == target_grp])
}

## Louvain resolution sweeps are not smooth: two nearby resolutions can
## land on distinct community structures, so a lone high-silhouette
## point separated from its neighbors by a crash or by small
## fluctuations could be either a genuinely better partition that the
## coarsest-near-max rule in sweep_cluster_resolution() missed (its tol
## window is centered on the chosen point, not built to also catch a
## detached higher peak elsewhere in the grid), or a small-N fluke that
## bootstrap resampling would not reproduce. This runs
## bootstrap_cluster_stability() at the chosen resolution and, if it
## differs, at the coarsest resolution of the contiguous plateau
## surrounding the grid's single highest silhouette (plateau_coarsest()),
## not the literal argmax itself, so the comparison is a fair contest
## between the coarsest representative of each candidate region rather
## than chosen vs. an arbitrary interior point of an otherwise-flat
## plateau.
##
## The FINAL resolution is then whichever candidate has the higher
## bootstrap mean ARI, the direct evidence of which partition actually
## reproduces under resampling, rather than always keeping the
## originally chosen resolution regardless of how that comparison comes
## out. Re-clusters sweep$obj (whose neighbor graph is already fitted)
## at each resolution being checked, rather than reusing sweep$obj's
## existing Idents, since those reflect only the chosen resolution.
## Returns a list: table (one row per resolution checked, a single row
## if the chosen resolution's own plateau already contains the highest
## silhouette), final_res (the resolution to actually use downstream),
## and final_obj (sweep$obj reclustered at final_res, ready to replace
## the caller's Seurat object for UMAP/DimPlot).
bootstrap_compare_resolutions <- function(sweep, counts, nfeatures, dims_n, metric = "manhattan", B = 20, seed = 1) {
  ok       <- sweep$grid[sweep$grid$ok, ]
  best_res <- plateau_coarsest(ok)
  res_list <- unique(c(sweep$chosen_res, best_res))

  rows <- lapply(res_list, function(r) {
    cl   <- Idents(FindClusters(sweep$obj, resolution = r, verbose = FALSE))
    boot <- bootstrap_cluster_stability(counts, cl, nfeatures = nfeatures, dims_n = dims_n, resolution = r, metric = metric, B = B, seed = seed)
    data.frame(res = r, role = if (r == sweep$chosen_res) "chosen" else "highest silhouette",
               sil = ok$sil[ok$res == r], n_clusters = ok$n_clusters[ok$res == r],
               boot_mean_ari = boot$mean_ari, boot_min_ari = min(boot$ari), boot_max_ari = max(boot$ari))
  })
  table <- do.call(rbind, rows)

  final_res <- table$res[which.max(table$boot_mean_ari)]
  final_obj <- FindClusters(sweep$obj, resolution = final_res, verbose = FALSE)
  list(table = table, final_res = final_res, final_obj = final_obj)
}

## Prints bootstrap_compare_resolutions()'s comparison table and, if the
## bootstrap-validated final resolution differs from the resolution
## sweep_cluster_resolution() originally chose, a note naming the
## override and the bootstrap mean ARI at each so the switch is
## traceable in the log rather than silent.
report_bootstrap_compare <- function(boot, label) {
  cat(sprintf("%s: bootstrap comparison across candidate resolutions\n", label)); print(boot$table)
  chosen_row <- boot$table[boot$table$role == "chosen", ]
  if (boot$final_res != chosen_row$res) {
    final_ari <- boot$table$boot_mean_ari[boot$table$res == boot$final_res]
    cat(sprintf("%s: switching to resolution %.2f (bootstrap mean ARI %.3f vs %.3f at the originally chosen %.2f)\n",
                label, boot$final_res, final_ari, chosen_row$boot_mean_ari, chosen_row$res))
  }
}

## ---- Cluster execution of the stability bootstrap (Section 7.3) ----
## cluster_stability_inputs: packages everything cluster_stability.R
## needs. For each dataset it records the candidate resolutions
## (chosen, plus the coarsest point of the highest-silhouette plateau),
## the reference partition at each candidate (computed here on the
## fitted neighbor graph, so the local and cluster sides share one
## labelling), sparse counts, and the pre-drawn resample matrix. key
## fingerprints the inputs so the returning output can be matched to
## them.
cluster_stability_inputs <- function(sweeps, counts, nfeatures, dims_n, B = 100, seed = 1, metric = "manhattan") {
  ds <- names(sweeps)
  tasks <- do.call(rbind, lapply(ds, function(d) {
    ok <- sweeps[[d]]$grid[sweeps[[d]]$grid$ok, ]
    rl <- unique(c(sweeps[[d]]$chosen_res, plateau_coarsest(ok)))
    data.frame(dataset = d, res = rl,
               role = ifelse(rl == sweeps[[d]]$chosen_res, "chosen", "highest silhouette"),
               sil = ok$sil[match(rl, ok$res)], n_clusters = ok$n_clusters[match(rl, ok$res)],
               stringsAsFactors = FALSE)
  }))
  tasks$task <- sprintf("%s@%.2f", tasks$dataset, tasks$res)
  ref <- setNames(lapply(seq_len(nrow(tasks)), function(k)
    Idents(FindClusters(sweeps[[tasks$dataset[k]]]$obj, resolution = tasks$res[k], verbose = FALSE))), tasks$task)
  data <- setNames(lapply(ds, function(d) {
    stopifnot(ncol(counts[[d]]) == ncol(sweeps[[d]]$obj))
    list(counts = as(counts[[d]], "CsparseMatrix"), nfeatures = nfeatures[[d]], dims_n = dims_n[[d]], metric = metric)
  }), ds)
  idx <- setNames(lapply(ds, function(d) make_boot_idx(ncol(counts[[d]]), B, seed)), ds)
  key <- list(tasks = tasks, cells = lapply(counts, colnames), B = B, seed = seed)
  list(tasks = tasks, ref = ref, data = data, idx = idx, key = key)
}

## diet_for_markers: keeps the counts and log-normalized data layers
## that FindMarkers() reads, which keeps the cluster input file small.
diet_for_markers <- function(obj) {
  if (packageVersion("Seurat") >= "5.0.0") DietSeurat(obj, layers = c("counts", "data")) else DietSeurat(obj)
}

## apply_cluster_labels: installs a partition on a Seurat object as its
## identities, the same state FindClusters() leaves behind.
apply_cluster_labels <- function(obj, labels) {
  stopifnot(identical(names(labels), colnames(obj)))
  Idents(obj) <- labels
  obj$seurat_clusters <- labels
  obj
}

## assemble_cluster_stability: turns the returned ARI vectors for one
## dataset into the same list bootstrap_compare_resolutions() returns
## (table, final_res, final_obj), so report_bootstrap_compare() and all
## downstream code read it unchanged. The final partition is the
## reference labelling at the winning resolution. n_ok counts
## replicates that completed.
assemble_cluster_stability <- function(inputs, ari, dataset, obj) {
  rows <- inputs$tasks[inputs$tasks$dataset == dataset, ]
  a    <- ari[rows$task]
  table <- data.frame(res = rows$res, role = rows$role, sil = rows$sil, n_clusters = rows$n_clusters,
                      boot_mean_ari = vapply(a, mean, numeric(1), na.rm = TRUE),
                      boot_min_ari  = vapply(a, min,  numeric(1), na.rm = TRUE),
                      boot_max_ari  = vapply(a, max,  numeric(1), na.rm = TRUE),
                      n_ok = vapply(a, function(x) sum(is.finite(x)), integer(1)), row.names = NULL)
  win <- which.max(table$boot_mean_ari)
  list(table = table, final_res = rows$res[win], final_obj = apply_cluster_labels(obj, inputs$ref[[rows$task[win]]]), ari = a)
}

## ---- Shared helpers for cluster round trips ----
## kegg_local: downloads the organism's KEGG pathway map once, locally,
## so cluster nodes run KEGG enrichment offline through enricher() with
## the same pathway release used everywhere in the run.
kegg_local <- function(org = "sce") {
  kg <- download_KEGG(org)
  kg$KEGGPATHID2NAME[[2]] <- sub(" - Saccharomyces cerevisiae \\(budding yeast\\)$", "", kg$KEGGPATHID2NAME[[2]])
  kg
}

## load_cluster_output: loads a cluster result into the caller's
## environment, with a message naming the script to run when the file
## is not there yet.
load_cluster_output <- function(path, script, envir = parent.frame()) {
  if (!file.exists(path))
    stop(sprintf("%s not found. Copy the new inputs to the cluster, run %s, copy the output back, then rerun this subsection.",
                 basename(path), script), call. = FALSE)
  invisible(load(path, envir = envir))
}

## check_cluster_key: confirms a cluster output was built from the inputs
## packaged in this session, the same safeguard score_promoters_nupop()
## applies to NuPoP coordinates.
check_cluster_key <- function(out_key, in_key, what, script) {
  if (!identical(out_key, in_key))
    stop(sprintf("%s was built from different inputs than those just packaged. Rerun %s with the new inputs.", what, script), call. = FALSE)
  invisible(TRUE)
}

## go_enrich_one: one unit of work for go_enrich.R. Over-representation
## jobs call run_enrichment(); rank-based jobs call gseGO() with a
## job-specific seed so permutation p-values reproduce.
go_enrich_one <- function(job, kegg_data, qval, seed) {
  if (job$kind == "ora") return(run_enrichment(job$genes, job$universe, qval = qval, kegg_data = kegg_data))
  set.seed(seed)
  suppressWarnings(gseGO(geneList = job$ranks, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = job$ont, nPermSimple = 100000))
}

## Clusters the same retained PCs under two annoy.metric choices at the
## same resolution and returns the adjusted Rand index between the two
## partitions: a one-time robustness check for whether the Manhattan
## vs. Euclidean distance choice materially changes which cells group
## together, rather than a step re-run as part of the normal pipeline.
## High ARI means the metric choice doesn't change the biology and
## either is defensible; low ARI means the choice is load-bearing and
## belongs in Methods as a stated, justified decision.
compare_distance_metrics <- function(obj, dims, resolution, metrics = c("manhattan", "euclidean")) {
  cl <- lapply(metrics, function(m) {
    o2 <- FindNeighbors(obj, reduction = "pca", dims = dims, annoy.metric = m, verbose = FALSE)
    o2 <- FindClusters(o2, resolution = resolution, verbose = FALSE)
    Idents(o2)
  })
  names(cl) <- metrics
  list(clusters = cl, ari = adjustedRandIndex(as.integer(cl[[1]]), as.integer(cl[[2]])))
}

##############################################################################
## 12. CLUSTER-BASED NOISE PARTITIONING                                    ##
##############################################################################
## Within- and between-cluster variance per gene, from a Seurat cluster
## partition, in the same mean-normalized rate space (count / exposure,
## rescaled to mean 1) and with the same Poisson shot-noise correction
## that intrinsic_extrinsic_components() uses for the hybrid allele-pair
## decomposition, so the two are on comparable footing: BETWEEN-cluster
## variance is the cell-state-driven component (variance in per-cluster
## mean expression), the "extrinsic-like" analog; WITHIN-cluster
## variance is the residual among cells sharing a cluster after
## subtracting shot noise, the "intrinsic-like" analog.
##
## mat: genes x cells count matrix. expo: per-cell exposure (same column
## order as mat). clusters: per-cell cluster label (same order).
##
## Shot noise is subtracted from the within-cluster component only:
## between-cluster differences are differences in cluster MEANS, which
## are not directly inflated by per-cell Poisson sampling noise the way
## within-cluster variance is. This ignores the smaller, n_c-dependent
## sampling noise of the group means themselves; worth revisiting if
## small clusters turn out to matter for this specific comparison.
## var_within can come out slightly negative for a gene whose true
## biological within-cluster variance is near zero (shot noise
## subtraction is only unbiased on average, not gene by gene); it is
## returned as-is for transparency and only floored at zero where used
## as a ratio, the same convention intrinsic_extrinsic_components() uses
## for extr.
##
## Returns a list: table (one row per gene: var_within, var_between,
## ratio_within_between) and cluster_n (one row per cluster: cell
## count), so cluster sizes travel with the result rather than being
## dropped after this step.
within_between_decomp <- function(mat, expo, clusters) {
  stopifnot(ncol(mat) == length(expo), ncol(mat) == length(clusters))
  cl        <- as.character(clusters)
  cl_levels <- sort(unique(cl))
  n_c       <- table(cl)[cl_levels]
  N         <- ncol(mat)

  a  <- sweep(mat, 2, expo, "/")
  ma <- rowMeans(a)
  ok_gene <- is.finite(ma) & ma > 0
  ap <- a / pmax(ma, 1e-10)
  shot_cell <- sweep(mat, 2, expo^2, "/") / pmax(ma^2, 1e-10)

  raw_within  <- numeric(nrow(mat))
  var_between <- numeric(nrow(mat))
  for (cc in cl_levels) {
    idx   <- which(cl == cc)
    w     <- length(idx) / N
    cmean <- rowMeans(ap[, idx, drop = FALSE])
    cvar  <- rowMeans((ap[, idx, drop = FALSE] - cmean)^2)
    cshot <- rowMeans(shot_cell[, idx, drop = FALSE])
    raw_within  <- raw_within + w * (cvar - cshot)
    var_between <- var_between + w * (cmean - 1)^2
  }
  var_within <- raw_within
  ratio <- pmax(var_within, 0) / var_between
  ratio[!ok_gene | !is.finite(ratio)] <- NA_real_

  list(
    table = data.frame(gene = rownames(mat), var_within = var_within, var_between = var_between,
                        ratio_within_between = ratio, row.names = NULL),
    cluster_n = data.frame(cluster = cl_levels, n_cells = as.integer(n_c), row.names = NULL)
  )
}

## Genome-wide distribution of within_between_decomp()'s ratio for one
## dataset. Values above 1 (log10 > 0) mean within-cluster variance
## exceeds between-cluster variance for that gene, i.e. residual noise
## among cells sharing a cluster outweighs the variance attributable to
## cluster identity itself. A large fraction of genes above 1 says the
## clustering is only capturing a modest share of total variance, mild
## cell-state separation rather than sharp, well-separated states, which
## fits the weak silhouette values already seen for these datasets in
## Section 8. Returns n and the fraction above 1 invisibly for logging.
plot_within_between_hist <- function(wb, main = NULL, brk = 40) {
  x <- wb$table$ratio_within_between
  x <- x[is.finite(x) & x > 0]
  lx <- log10(x)
  frac_above_1 <- mean(x > 1)
  h <- hist(lx, breaks = brk, plot = FALSE)
  if (is.null(main)) main <- sprintf("n = %d genes, %.0f%% with within > between", length(x), 100 * frac_above_1)
  plot(h, col = "grey85", border = "white", xlab = "log10(within / between)", ylab = "# of genes", main = main)
  abline(v = 0, lty = 2, col = "red")
  invisible(list(n = length(x), frac_above_1 = frac_above_1))
}

## Scatter of within_between_decomp()'s ratio (log10) against one burst
## kinetics quantity, for one dataset. quantity is a named vector (names
## = gene, already on the desired plotting scale, e.g. log2(mean)) so the
## same function covers mean, burst frequency, and burst size without
## three near-duplicate copies. Comparison is by rank (Spearman), the
## same convention used for the now-dropped intrinsic/extrinsic check,
## since the ratio can span orders of magnitude without needing a linear
## relationship to the burst quantity. Returns n and rho invisibly.
plot_within_between_vs_quantity <- function(wb, quantity, xlab, main = NULL) {
  idx <- match(wb$table$gene, names(quantity))
  ok  <- is.finite(wb$table$ratio_within_between) & wb$table$ratio_within_between > 0 & !is.na(idx) & is.finite(quantity[idx])
  x   <- quantity[idx[ok]]; y <- log10(wb$table$ratio_within_between[ok])
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.null(main)) main <- sprintf("n = %d, Spearman rho = %.3f", sum(ok), rho)
  plot(x, y, pch = 16, cex = 0.4, col = "#00000060", xlab = xlab, ylab = "log10(within / between)", main = main)
  abline(lm(y ~ x), col = "red", lty = 2)
  invisible(list(n = sum(ok), rho = rho))
}

## Cross-species scatter of within_between_decomp()'s ratio (log10), one
## point per ortholog gene pair, for two datasets whose gene sets are
## already ortholog-matched by name (e.g. MIX.SC vs MIX.SE, both indexed
## by the shared GENES ortholog-pair set built earlier in the pipeline,
## so no additional ortholog mapping is needed here). Spearman rho
## reported, same rank-based convention used throughout this section.
plot_within_between_cross_species <- function(wb_a, wb_b, lab_a, lab_b, main = NULL) {
  m  <- merge(wb_a$table[, c("gene", "ratio_within_between")], wb_b$table[, c("gene", "ratio_within_between")], by = "gene", suffixes = c("_a", "_b"))
  ok <- is.finite(m$ratio_within_between_a) & m$ratio_within_between_a > 0 & is.finite(m$ratio_within_between_b) & m$ratio_within_between_b > 0
  x   <- log10(m$ratio_within_between_a[ok]); y <- log10(m$ratio_within_between_b[ok])
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  if (is.null(main)) main <- sprintf("n = %d genes, Spearman rho = %.3f", sum(ok), rho)
  plot(x, y, pch = 16, cex = 0.4, col = "#00000060",
       xlab = sprintf("log10(within / between), %s", lab_a), ylab = sprintf("log10(within / between), %s", lab_b), main = main)
  abline(0, 1, lty = 3, col = "grey60")
  abline(lm(y ~ x), col = "red", lty = 2)
  invisible(list(n = sum(ok), rho = rho))
}

## Boxplot of within_between_decomp()'s ratio (log10), split by a
## regulatory or dominance class, for one dataset, with class_anova()'s
## omnibus test run internally and returned invisibly (F, df, p, eta^2,
## and Tukey pairwise comparisons if the omnibus test clears sig). class
## should already be aligned to wb$table$gene (match() against
## BURST.CONTRASTS$gene) with Ambiguous genes set to NA; class_levels
## fixes the plotting/level order (e.g. REG.CLASS or DOM.CLASS) and
## colors should be the matching per-level color vector (e.g.
## COLOR.LIST.1 or COLOR.LIST.2). main is left to the caller (e.g. a
## species label); the test statistics are reported separately rather
## than crowding the plot title, the same convention used for the
## bootstrap comparison reporting in Section 8.
plot_within_between_by_class <- function(wb, class, class_levels, colors, main = NULL) {
  x  <- log10(wb$table$ratio_within_between)
  ok <- is.finite(x) & !is.na(class)
  x  <- x[ok]; cl <- factor(class[ok], levels = class_levels)
  res <- class_anova(x, cl)
  boxplot(x ~ cl, col = colors, las = 2, ylab = "log10(within / between)", main = main, border = "grey30")
  invisible(res)
}

## Runs plot_within_between_by_class() for two datasets side by side
## (e.g. Sc and Se parent) against one classification axis, in one call:
## aligns the raw class vector (indexed like BURST.CONTRASTS$gene, e.g.
## REG.BFREQ.CLASS or DOM.BSIZE.CLASS, not pre-cleaned since Ambiguous
## genes are already excluded automatically by class_levels not
## containing "Ambiguous") to each dataset's own gene table via match(),
## writes both panels to one pdf, and prints the omnibus test and Tukey
## pairwise comparisons (when significant) for both datasets. Returns
## both class_anova() results invisibly, named by lab_a/lab_b.
report_within_between_by_class <- function(wb_a, wb_b, class, gene_ref, class_levels, colors,
                                            lab_a, lab_b, axis_label, pdf_path, width = 11, height = 5.5) {
  cls_a <- class[match(wb_a$table$gene, gene_ref)]
  cls_b <- class[match(wb_b$table$gene, gene_ref)]

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
  res_a <- plot_within_between_by_class(wb_a, cls_a, class_levels, colors, main = lab_a)
  res_b <- plot_within_between_by_class(wb_b, cls_b, class_levels, colors, main = lab_b)
  dev.off()

  cat(sprintf("%s: within/between by %s, F(%d,%d) = %.2f, p = %.3g, eta^2 = %.3f\n",
              lab_a, axis_label, res_a$df1, res_a$df2, res_a$f, res_a$p, res_a$eta_sq))
  cat(sprintf("%s: within/between by %s, F(%d,%d) = %.2f, p = %.3g, eta^2 = %.3f\n",
              lab_b, axis_label, res_b$df1, res_b$df2, res_b$f, res_b$p, res_b$eta_sq))
  if (!is.null(res_a$tukey)) { cat(sprintf("%s: Tukey pairwise (%s)\n", lab_a, axis_label)); print(res_a$tukey) }
  if (!is.null(res_b$tukey)) { cat(sprintf("%s: Tukey pairwise (%s)\n", lab_b, axis_label)); print(res_b$tukey) }

  invisible(setNames(list(res_a, res_b), c(lab_a, lab_b)))
}

## ============================================================
## 13. POWER ANALYSIS (simulation and fit, shared with the SLURM job power_grid.R)
## ============================================================
## Offset NB fit for one simulated gene, matching .fit_one() above: mu is
## the exposure-weighted mean, size is the MLE of the NB dispersion given
## that mu, found by direct optimization of the log-likelihood over
## log(theta). Kept as its own small fit rather than reusing .fit_one()
## so the power simulation has no dependency on glm.nb, which is far
## slower over the many thousands of replicate fits the grid requires.
fit_offset_nb <- function(y, expo) {
  keep <- is.finite(y) & is.finite(expo) & expo > 0
  y <- y[keep]; expo <- expo[keep]
  if (length(y) < 2) return(c(mu = NA_real_, size = NA_real_))

  mu_hat <- sum(y) / sum(expo)
  if (sum(y) == 0) return(c(mu = 0, size = NA_real_))

  mu_i <- mu_hat * expo
  pearson <- sum((y - mu_i)^2 / mu_i) / (length(y) - 1)
  if (pearson <= 1) return(c(mu = mu_hat, size = Inf))

  ll <- function(ltheta) {
    th <- exp(ltheta)
    sum(lgamma(y + th) - lgamma(th) + th*log(th) - (th + y)*log(th + mu_i) + y*log(mu_i))
  }

  opt <- tryCatch(optimize(ll, c(-4, 15), maximum = TRUE), error = function(e) NULL)
  if (is.null(opt)) return(c(mu = mu_hat, size = NA_real_))

  theta <- exp(opt$maximum)
  c(mu = mu_hat, size = if (theta > 1e6) Inf else theta)
}

## Method-of-moments dispersion estimate, closed-form, used for
## permutation-null replicates. Derivation: the Pearson statistic
## sum((y-mu)^2/mu) has expectation (n-1) + sum(mu)/theta under the NB
## model, so solving for theta gives a direct estimate with no
## iteration. This is the fast estimator that makes the NI permutation
## draws per simulated dataset affordable at the scale of the power grid.
fit_offset_nb_mm <- function(y, expo) {
  keep <- is.finite(y) & is.finite(expo) & expo > 0
  y <- y[keep]; expo <- expo[keep]
  n <- length(y)
  if (n < 2) return(c(mu = NA_real_, size = NA_real_))

  mu_hat <- sum(y) / sum(expo)
  if (sum(y) == 0) return(c(mu = 0, size = NA_real_))

  mu_i <- mu_hat * expo
  pearson <- sum((y - mu_i)^2 / mu_i) / (n - 1)
  if (pearson <= 1) return(c(mu = mu_hat, size = Inf))

  theta_mm <- sum(mu_i) / ((n - 1) * (pearson - 1))
  c(mu = mu_hat, size = theta_mm)
}

## Splits a pooled sample at position n1 using a pre-drawn permutation
## of indices and fits each half with the MLE, mirroring .fit_split()
## above for the power simulation's own fit function.
fit_split_nb <- function(y, expo, perm, n1) {
  g1 <- perm[seq_len(n1)]; g0 <- perm[(n1 + 1):length(perm)]
  list(a = fit_offset_nb(y[g1], expo[g1]), b = fit_offset_nb(y[g0], expo[g0]))
}

## Same split, fit with the method-of-moments estimator. Used inside the
## permutation loop, which runs NI times per simulated dataset.
fit_split_nb_mm <- function(y, expo, perm, n1) {
  g1 <- perm[seq_len(n1)]; g0 <- perm[(n1 + 1):length(perm)]
  list(a = fit_offset_nb_mm(y[g1], expo[g1]), b = fit_offset_nb_mm(y[g0], expo[g0]))
}

## Two-sided permutation p-value with add-one continuity correction,
## matching pval() used by the real-data permutation null.
perm_pval <- function(obs, null) {
  ok <- is.finite(null)
  if (!is.finite(obs) || sum(ok) < 1) return(NA_real_)
  (1 + sum(abs(null[ok]) >= abs(obs))) / (1 + sum(ok))
}

## Log2 ratio of SIZE (burst frequency) between two fits, NA if either
## side is non-positive or non-finite. Mirrors rat() used elsewhere for
## real-data contrasts.
size_log2_ratio <- function(fit_a, fit_b) {
  a <- fit_a["size"]; b <- fit_b["size"]
  if (is.finite(a) && a > 0 && is.finite(b) && b > 0) log2(a) - log2(b) else NA_real_
}

## ============================================================
## Power for one (MEAN.READS, N.CELLS, SIZE) row, across every
## SIZE.RATIO value at once. This is the unit of work parLapply
## distributes in power_grid.R: one call per row of GRID, returning a
## vector of power values, one per SIZE.RATIO.
##
## Exposures and permutation index sets are seeded from row$n (the cell
## count) alone, so any two rows sharing a cell-count level draw the
## identical sequence. The reference group's counts and its full-MLE
## fit are drawn once per replicate and reused across the SIZE.RATIO
## sweep below, so the sweep isolates the effect of SIZE.RATIO alone.
##
## N.CELLS.X gives the reference (Sc) group's cell count; the second
## group is drawn at round(N.CELLS.X * CELL.RATIO) cells. CELL.RATIO
## defaults to 1 (equal group sizes) in the driver.
## ============================================================
power_grid_row <- function(i, GRID, MEAN.READS, N.CELLS, SIZE, SIZE.RATIO, EXPOSURE.CV, CELL.RATIO, ALPHA, NJ, NI, SEED.BASE) {
  row <- GRID[i, ]
  MEAN.READS.X <- MEAN.READS[row$m]
  N.SC.X       <- N.CELLS[row$n]
  N.SE.X       <- round(N.SC.X * CELL.RATIO)
  SIZE.1.X     <- SIZE[row$p]

  set.seed(SEED.BASE + 1e6 + row$n)  #shared draws, keyed only by cell count
  sdlog <- sqrt(log(1 + EXPOSURE.CV^2))
  EXPO.X.LIST <- lapply(seq_len(NJ), function(j) rlnorm(N.SC.X, meanlog = -0.5*sdlog^2, sdlog = sdlog))
  EXPO.Y.LIST <- lapply(seq_len(NJ), function(j) rlnorm(N.SE.X, meanlog = -0.5*sdlog^2, sdlog = sdlog))
  PERM.LIST   <- lapply(seq_len(NI), function(k) sample.int(N.SC.X + N.SE.X))

  set.seed(SEED.BASE + i)  #row-specific draws (X depends on m, p, n)
  X.LIST  <- vector("list", NJ)
  FA.LIST <- vector("list", NJ)
  for (j in seq_len(NJ)) {
    X.LIST[[j]]  <- rnbinom(n = N.SC.X, size = SIZE.1.X, mu = MEAN.READS.X*EXPO.X.LIST[[j]])
    FA.LIST[[j]] <- fit_offset_nb(X.LIST[[j]], EXPO.X.LIST[[j]])  #full MLE, shared across SIZE.RATIO below
  }

  POWER.ROW <- numeric(length(SIZE.RATIO))
  for (qi in seq_along(SIZE.RATIO)) {
    SIZE.2.X <- SIZE.1.X / SIZE.RATIO[qi]
    P.VALUE <- numeric(NJ)
    for (j in seq_len(NJ)) {
      EXPO.Y <- EXPO.Y.LIST[[j]]
      Y  <- rnbinom(n = N.SE.X, size = SIZE.2.X, mu = MEAN.READS.X*EXPO.Y)
      fb <- fit_offset_nb(Y, EXPO.Y)   #full MLE for the observed contrast
      OBS <- size_log2_ratio(FA.LIST[[j]], fb)

      XY <- c(X.LIST[[j]], Y); EXPO.XY <- c(EXPO.X.LIST[[j]], EXPO.Y)
      NULL.DIST <- numeric(NI)
      for (k in seq_len(NI)) {
        sp <- fit_split_nb_mm(XY, EXPO.XY, PERM.LIST[[k]], N.SC.X)  #MM estimate for the null
        NULL.DIST[k] <- size_log2_ratio(sp$a, sp$b)
      }
      P.VALUE[j] <- perm_pval(OBS, NULL.DIST)
    }
    POWER.ROW[qi] <- sum(P.VALUE < ALPHA, na.rm = TRUE) / NJ
  }
  POWER.ROW
}

## Opens a PDF sized to the nr x nc panel grid it is about to hold.
## Margins and label spacing (mgp) are set tight throughout, so every
## power-analysis figure shares the same compact panel layout.
open_grid_pdf <- function(file, nr, nc, panel_w = 2.3, panel_h = 2.3, mar = c(3, 3, 2, 1), mgp = c(1.6, 0.5, 0)) {
  pdf(file, width = nc * panel_w, height = nr * panel_h, useDingbats = FALSE)
  par(mfrow = c(nr, nc), mar = mar, mgp = mgp)
}

## Color ramp with n colors, one per line in a panel.
line_colors <- function(n) colorRampPalette(POWER.COLOR)(n)

## Plots one line per column of `mat` against `xv`, using a color ramp
## sized to ncol(mat). Reference lines at power 0.05 (nominal false
## positive rate) and 0.9 (a common target) give every panel the same
## visual anchor regardless of which quantity is on the x-axis.
plot_lines <- function(xv, mat, xlab = "", ylab = "Power", main = "", show_axes = TRUE, lwd = 1.4, x_at = NULL, x_labels = NULL) {
  cols <- line_colors(ncol(mat))
  plot(xv, mat[, 1], type = "l", lwd = lwd, col = cols[1], ylim = c(0, 1), xlab = if (show_axes && is.null(x_at)) xlab else "", ylab = if (show_axes) ylab else "", xaxt = if (!is.null(x_at) || !show_axes) "n" else "s", yaxt = if (show_axes) "s" else "n", main = main, cex.main = 0.8)
  abline(h = 0.05, col = "red"); abline(h = 0.9, col = "blue")
  if (ncol(mat) > 1) for (k in 2:ncol(mat)) lines(xv, mat[,k], lwd = lwd, col = cols[k])
  #x_at/x_labels give exact tick positions and real-valued labels for a
  #log-scaled axis (e.g. log2(SIZE.RATIO) with SIZE.RATIO as labels)
  if (!is.null(x_at)) {
    axis(1, at = x_at, labels = x_labels, cex.axis = 0.55, las = 2)
    title(xlab = xlab, mgp = c(2.6,1,0), cex.lab = 0.75)
  }
  invisible(cols)
}

## One legend page appended after a grid of panels, mapping each line
## color to the value it represents.
legend_page <- function(labels, title) {
  cols <- line_colors(length(labels))
  plot.new()
  legend("center", legend = labels, col = cols, lwd = 2, title = title, ncol = ceiling(length(labels)/15), bty = "n", cex = 0.9)
}

## Minimum SIZE.RATIO detectable at a target power, one number per
## (MEAN.READS, N.CELLS, SIZE) instead of the full power-vs-ratio curve.
## Lower means more sensitive. Returns NA when even the largest tested
## ratio fails to reach the target, which image() then renders as blank.
min_detectable_ratio <- function(power_vec, ratios, target = 0.8) {
  hit <- which(power_vec >= target)
  if (length(hit) == 0) return(NA_real_)
  ratios[min(hit)]
}

## ============================================================
## 14. EXTERNAL NOISE VALIDATION
## ============================================================
## Spearman rank correlation between two vectors, keeping only pairwise-
## complete (finite, finite) observations. Rank correlation is used
## throughout the external validation since both the NB-fitted
## quantities and the external abundance and noise statistics are
## heavy-tailed, so rank correlation avoids any assumption about the
## shape of that tail.
cor_row <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  # exact = FALSE belongs to cor.test() itself: it tells Spearman's test
  # to use the asymptotic t approximation for the p-value rather than
  # the exact permutation null, which is undefined once ranks contain
  # ties and otherwise triggers a warning on every tied comparison here.
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(n = sum(keep), rho = unname(ct$estimate), p = ct$p.value)
}

## Implied Fano factor, burst size, and burst frequency algebraically
## recovered from a reported mean and CV^2, for an external data frame
## with columns Mean and CV2. The same NB identity used for MIX.SC runs
## in reverse: Fano = mean * CV^2, implied burst size = Fano - 1, and
## implied burst frequency = mean / (Fano - 1). This is an algebraic
## transform of the reported mean and CV^2, not a re-fit, and it carries
## an interpretive caveat worth keeping in mind downstream: burst
## frequency and burst size are mRNA-level transcriptional-kinetics
## constructs, while the protein-level external sources fold
## translational and, where relevant, degradation noise on top of
## transcriptional bursting. Fano <= 1 (sub-Poissonian) is not expected
## biologically and is treated as measurement noise, so it is set to NA
## rather than producing a negative burst size.
add_burst_terms <- function(d) {
  fano <- d$Mean * d$CV2
  ok <- fano > 1
  d$Fano  <- fano
  d$BSIZE <- ifelse(ok, fano - 1, NA_real_)
  d$BFREQ <- ifelse(ok, d$Mean / (fano - 1), NA_real_)
  d
}

## read_header_line() parses just the header line of a delimited file
## directly, independent of fread()'s own header detection. Several of
## the scRNA-seq source files write their header with one fewer field
## than the data rows carry, since the row-identifying column has no
## name of its own there, the same implicit-row-names convention
## read.delim(row.names = 1) resolves automatically but fread() does
## not, so reading the header line by itself, and counting on it being
## one field short, avoids the field-count mismatch that otherwise
## confuses fread()'s own header heuristics on these files.
read_header_line <- function(path, sep = "\t", skip = 0) {
  lines  <- readLines(path, n = skip + 1)
  fields <- strsplit(lines[skip + 1], sep, fixed = TRUE)[[1]]
  gsub('^"|"$', "", fields)
}

## get_data_header() returns just the sample/gene-name portion of a
## file's header line, with the row-identifying field removed no matter
## how that file represents it. Two conventions show up across the four
## scRNA-seq sources: Jackson's header omits the row-id field entirely,
## so every data row carries one more field than the header does;
## Gasch's header gives that field an explicit label ("gene") and
## Jariani's header gives it a blank placeholder, so in both of those
## the header and every data row have the same field count. Comparing
## the header's own field count to one data row's field count tells the
## two conventions apart without hard-coding which one a given file
## uses, so fread_matrix() below always ends up with a header that
## describes only the true data columns and is always exactly one field
## shorter than a data row, regardless of which convention the source
## file used.
get_data_header <- function(path, sep = "\t", skip = 0) {
  header.raw <- read_header_line(path, sep = sep, skip = skip)
  data.line  <- readLines(path, n = skip + 2)[skip + 2]
  n.data     <- length(strsplit(data.line, sep, fixed = TRUE)[[1]])
  if (length(header.raw) == n.data - 1) {
    header.raw            # header already omits the row-id field
  } else if (length(header.raw) == n.data) {
    header.raw[-1]         # header carries a label (or blank) for the row-id field; drop it
  } else {
    stop("get_data_header(): header has ", length(header.raw),
         " field(s) but the first data row has ", n.data,
         " field(s) in ", path, "; this file matches neither the ",
         "'header omits row-id' nor the 'header labels row-id' convention")
  }
}

## fread_matrix() is a fast stand-in for read.delim(path, row.names = 1)
## followed by as.matrix(), built on data.table's fread(). keep is a
## logical vector the same length as get_data_header()'s output, marking
## which of those columns to parse; leaving it NULL keeps them all. The
## data rows are read with header = FALSE, so fread() never has to
## reconcile the header's own offset itself; each row's own first field
## becomes the matrix's row names, and keep narrows the parse to only
## the columns a given source actually needs rather than reading every
## column and discarding most of them afterward. Because get_data_header()
## already normalizes away whichever row-id convention the file uses,
## the +1 shift below is always correct, for every source.
fread_matrix <- function(path, keep = NULL, sep = "\t", skip = 0, ...) {
  header <- get_data_header(path, sep = sep, skip = skip)
  if (is.null(keep)) keep <- rep(TRUE, length(header))
  col.idx <- c(1L, which(keep) + 1L)   # +1 shifts past the row-id field, now always absent from header
  dt  <- fread(file = path, sep = sep, skip = skip + 1, header = FALSE, select = col.idx, ...)
  mat <- as.matrix(dt[, -1, with = FALSE])
  colnames(mat) <- header[keep]
  rownames(mat) <- dt[[1]]
  mat
}

## qc_filter_counts() drops cells below a minimum library size and genes
## detected in too few cells, the same two checks applied before any NB
## fit elsewhere in the pipeline. Thresholds are left as arguments rather
## than hard-coded, since the four sources differ enormously in depth and
## cell number and each threshold should be sanity-checked against that
## source's own count distribution before the fit is trusted.
qc_filter_counts <- function(mat, min.cell.count = 1, min.cells.expr = 1) {
  keep.cells <- colSums(mat) >= min.cell.count
  mat <- mat[, keep.cells, drop = FALSE]
  keep.genes <- rowSums(mat > 0) >= min.cells.expr
  mat[keep.genes, , drop = FALSE]
}

## to_numeric_matrix() guards against the common gotcha where a single
## non-numeric column (a leftover metadata field from a Seurat export, a
## stray total-count column, a blank trailing column) makes as.matrix()
## fall back to character for the whole table, which then fails silently
## until something downstream demands numeric input. Any column that does
## not convert cleanly is reported by name and dropped rather than left
## to break rowsum() or the NB fit further down.
to_numeric_matrix <- function(df, label) {
  mat <- as.matrix(df)
  if (is.character(mat)) {
    num <- suppressWarnings(apply(mat, 2, as.numeric))
    # a column is genuinely non-numeric if it produced an NA somewhere
    # that was not already a blank entry in the original character data
    bad <- colSums(is.na(num) & mat != "", na.rm = TRUE) > 0
    if (any(bad))
      message(label, ": dropping ", sum(bad), " non-numeric column(s): ",
              paste(head(colnames(mat)[bad], 10), collapse = ", "))
    mat <- num[, !bad, drop = FALSE]
    rownames(mat) <- rownames(df)
  }
  storage.mode(mat) <- "numeric"
  mat
}

## map_to_orf() takes a vector of gene identifiers that may already be
## systematic ORF names or may be common names, and returns the ORF name
## for each. Anything that already matches the standard systematic-name
## pattern is passed through unchanged; everything else is looked up
## against org.Sc.sgd.db by common name. Identifiers that cannot be
## resolved (SUTs, CUTs, snoRNAs, ERCC spike-ins, and any common name not
## in SGD) come back as NA and are dropped by the caller.
## The dash that marks a dubious-ORF suffix (e.g. YAL047C-A) comes through
## as a period rather than a hyphen in at least the Jackson file, likely
## from an upstream make.names() call before the file was deposited, so
## the pattern accepts either separator.
ORF.PATTERN <- "^Y[A-P][LR][0-9]{3}[CW]([.-][A-Z])?$"
map_to_orf <- function(ids, label = "") {
  is.orf <- grepl(ORF.PATTERN, ids)
  out <- ids
  need <- unique(ids[!is.orf])
  if (length(need) > 0) {
    # mapIds() throws when NONE of the supplied keys exist under that
    # keytype, rather than just returning all-NA the way it does for a
    # partial mismatch. That total-mismatch case almost always means the
    # identifiers in "need" are not actually common names at all (most
    # often because the systematic-name regex above failed to match the
    # real format this source uses), so it is reported rather than left
    # to crash the whole script, with a sample of the unmatched strings
    # printed so the real format can be diagnosed directly.
    looked.up <- tryCatch(
      suppressMessages(mapIds(org.Sc.sgd.db, keys = need, keytype = "GENENAME",
                               column = "ORF", multiVals = "first")),
      error = function(e) {
        message(label, ": ", length(need), " identifier(s) did not match the systematic-ORF ",
                "pattern and also failed as GENENAME lookups. Sample: ",
                paste(head(need, 10), collapse = ", "))
        setNames(rep(NA_character_, length(need)), need)
      })
    out[!is.orf] <- unname(looked.up[match(ids[!is.orf], need)])
  }
  n.unresolved <- sum(is.na(out))
  if (n.unresolved > 0)
    message(label, ": ", n.unresolved, " of ", length(out), " row(s) did not resolve to an ORF and will be dropped")
  out
}

## collapse_to_orf() sums counts across any duplicate rows sharing the
## same resolved ORF (isoform- or TSS-level rows in Nadal-Ribelles collapse
## to one gene-level row this way) and drops rows that never mapped to
## an ORF at all.
collapse_to_orf <- function(mat, raw.ids, label = "") {
  orf <- map_to_orf(raw.ids, label)
  mat <- mat[!is.na(orf), , drop = FALSE]
  orf <- orf[!is.na(orf)]
  rowsum(mat, group = orf)
}

## fit_source() is the common path from a QC'd genes-by-cells count
## matrix to a data.frame keyed by ORF with the same MU/CV2/BFREQ/BSIZE
## columns NB.SC already carries, so every new source merges against
## NB.SC the same way.
##
## Exposure is scaled to the source's own average library size rather
## than left as raw colSums(). Raw colSums() puts MU on a per-read
## scale, where 1/MU is so large that it swamps the 1/DISP term in
## CV2 = 1/MU + 1/DISP, making CV2 track the mean almost exactly for
## every gene. Scaling exposure to the mean library size puts MU in
## counts per typical cell instead, matching the kind of unit MIX.SC
## uses, so 1/MU reflects this dataset's own sequencing depth and the
## 1/DISP term contributes its intended share of CV2.
##
## cl is optional. Leave it NULL to fit on one core. Pass a cluster
## object (from parallel::makeCluster()) to split genes across workers;
## the same cluster is reused across all four sources in Section 9 of
## the driver script, created once and stopped once, which avoids the
## leaked connections that come from making and (sometimes) forgetting
## to close a fresh cluster for each source individually.
fit_source <- function(mat, cl = NULL) {

  exposure <- colSums(mat) / mean(colSums(mat))

  if (is.null(cl)) {

    fit <- fit_counts_offset(mat, exposure)

  } else {

    # Genes are split into one chunk per worker. Splitting by gene
    # rather than by cell is what makes this safe: each gene is fit
    # independently of every other gene, so no information needs to
    # cross chunk boundaries, and the worker's chunk of the matrix
    # carries every cell for its genes.
    n.genes   <- nrow(mat)
    n.workers <- length(cl)
    chunk.id  <- cut(seq_len(n.genes), n.workers, labels = FALSE)
    chunks    <- split(seq_len(n.genes), chunk.id)
    mat.chunks <- lapply(chunks, function(i) mat[i, , drop = FALSE])

    # PSOCK workers start with an empty workspace on each call, so
    # anything the workers need is sent explicitly rather than assumed
    # to already be there. Doing this here, right before use, means a
    # cluster that gets reused across sources does not depend on state
    # left behind by an earlier call.
    parallel::clusterEvalQ(cl, library(MASS))
    parallel::clusterExport(cl, c("fit_counts_offset", "neg_binom_fit_offset"),
                             envir = globalenv())

    fit.chunks <- parallel::parLapply(cl, mat.chunks, function(m, e) {
      fit_counts_offset(m, e)
    }, e = exposure)

    # unname() before rbind() keeps ORF row names clean. rbind() on a
    # named list would otherwise prefix every row name with its chunk
    # number, which breaks the later merge(..., by = "ORF").
    fit <- do.call(rbind, unname(fit.chunks))
  }

  # BFREQ and BSIZE are read straight from fit_counts_offset()'s own
  # BFREQ and BSIZE columns, the same DISP-based definition documented
  # for NB.SC in Section 9.2 of the driver script (burst frequency =
  # DISP, burst size = MU / DISP).
  data.frame(
    ORF   = rownames(fit),
    Mean  = fit$MU,
    CV2   = fit$CV^2,
    Fano  = fit$FANO,
    BFREQ = fit$BFREQ,
    BSIZE = fit$BSIZE,
    stringsAsFactors = FALSE
  )
}

## mean_adjusted_noise() expresses each gene's noise relative to other
## genes of similar abundance in the same dataset, the same idea as
## Newman et al.'s DM statistic. Fitting the trend within each source
## lets every dataset keep its own capture efficiency and units, so the
## residual isolates the gene-specific part of noise that should travel
## across datasets, independent of each dataset's own abundance-noise
## relationship.
mean_adjusted_noise <- function(mean, cv2, span = 0.3) {
  lm_ <- log(mean)
  lc  <- log(cv2)
  ok  <- is.finite(lm_) & is.finite(lc)
  res <- rep(NA_real_, length(mean))
  fit <- loess(lc[ok] ~ lm_[ok], span = span)   # smooth abundance trend, as in Newman's DM
  res[ok] <- lc[ok] - predict(fit)              # residual = noise beyond the expected level at that abundance
  res
}
