###############################################################
### functions.R
### Function library for the cis/trans mean and noise analysis,
### sourced by analysis.R and by the SLURM cluster
### scripts (gene_pilot.R, gene_boot.R, gene_perm.R,
### coexpr_boot.R, coexpr_perm.R, nupop_occupancy.R,
### cluster_stability.R, go_enrich.R, power_grid.R)
###
### Two analyses keep their functions in their own files: the power analysis (functions_power.R, used by
### power_analysis.R and power_grid.R) and the comparison with published data (functions_external.R, used by
### external_validation.R). Both rely on the helpers here (.fit_one, .fit_split, perm_pval, fit_counts_offset,
### fig_pdf, ckpt_path, console_start).
###
### A function lives here when its code runs more than once: it is called from two or more places, or it is
### applied across many items (Section 13b for analysis.R, Section 13c for the cluster scripts). Code that
### runs once is written inline where it runs. analysis.R and the cluster scripts define no functions.
###
### Outline (function name - purpose), grouped by section:
###
###   0. SESSION HELPERS
###     ckpt_path() - Path of a section's checkpoint file (section{N}_checkpoint.rda).
###     console_start() / console_stop() - Per-section console transcript (section{N}_console.txt).
###     with_local_seed() - Evaluates an expression on its own seeded random stream and restores the caller's random number state.
###   1. OFFSET NEGATIVE-BINOMIAL FIT
###     neg_binom_fit_offset() - Offset NB fit for one gene: rate (mu) and NB size (disp) from .fit_one(), with glm.nb asymptotic log-scale SEs; disp = Inf when a Poisson-vs-NB pre-check finds no overdispersion.
###     fit_counts_offset() - Matrix version: one neg_binom_fit_offset() fit per gene (row) against its exposure vector; optionally splits genes across a PSOCK cluster, shipping the functions the fit needs.
###     .fit_one() - The pipeline's single NB estimator (mu = sum(y)/sum(exposure), disp by 1-D likelihood given mu): observed fits, every bootstrap and permutation replicate, and the power grid.
###     .fit_split() - Splits a pooled sample by a pre-drawn permutation and fits each half with .fit_one(); the gene-level null and the power grid share it.
###   1b. INTERNAL-PILOT SPLIT FRACTION (f*)
###     split_indices_by_depth() - Depth-matched split of depth-ordered hybrid cells at an arbitrary fraction.
###   2. GENE FILTER (marginal information only; never on a contrast)
###     qc_gene_keep() - Cross-dataset gene QC with one relative abundance floor anchored to the shallowest dataset and an absolute detection floor, and sets the pilot gene pool.
###   3. PAIRED PER-GENE BOOTSTRAP OF CONTRASTS
###     make_draws() - Builds the resampling index draws used by the paired per-gene bootstrap: parents resampled whole, hybrid cells resampled within the four strata of the mean-split by noise-split overlap and every hybrid dataset assembled from the same resampled cell IDs.
###     add_burst_contrasts() - Derives burst-frequency, burst-size, and kinetic-balance contrasts from the mean and size contrasts already in a data frame.
###     .contrast_value() - Log2 contrast for one mode (total, cis, trans, dom, dpar, inh) from a named list of per-group fitted values.
###     .cv2_of() - CV2 (1/mu + 1/disp) of one fitted group; NA unless mu and disp are finite and positive.
###   4. DISATTENUATED MEAN-DISP COUPLING (population summary)
###     eiv_components() - Errors-in-variables mean-bfreq coupling for one mode: attenuation-corrected variances, covariance and correlation, plus the raw correlation.
###   5. PLOTTING HELPERS AND GENE-LEVEL PLOTS
###     umap_plot() - UMAP with the shared cluster palette and on-plot labels.
###     fig_pdf() - Opens a PDF under FIGURE.DIR with the pipeline's figure settings.
###     class_count_barplot() - Gene (or pair) counts per class, bars colored by class.
###     class_overlap_triptych() - Figure of three class-overlap heatmaps: mean class against burst frequency, burst size and frequency-size balance class, on one shared color range.
###     class_heatmap_grid() - Figure of class-overlap heatmaps for a table of class-vector pairs (regulatory vectors cleaned with clean_reg()).
###     .se_scatter() - Shared SE-bar scatter body of the three gene-level scatters below.
###     plot_contrast_scatter() - Gene-level contrast scatter with SE bars: cis vs trans for one quantity, mean vs burst frequency or rotated burst kinetics for one mode.
###     sig_hist_panel() - Histogram of one contrast (quantity, mode) with the genes significant at q < sig shaded by direction.
###   6. PERMUTATION NULL  (per-gene significance for the contrasts)
###     perm_pval() - Two-sided permutation p-value with add-one continuity correction, shared with the power grid.
###   7. CO-EXPRESSION  (residual co-fluctuation, cis/trans decomposed)
###     nb_residuals() - Pearson residuals from the offset NB fit, genes x cells.
###     shrink_cor() - Analytic shrinkage of a correlation matrix toward the identity matrix.
###     coexpr_decompose() - Decomposes the point-estimate correlation structure into total/cis/trans/dpar_sc/dpar_se matrices.
###     coexpr_axis_cis_trans() - Exact cis/trans decomposition of one eigenvector's eigenvalue, via total = cis + trans.
###     make_coexpr_draws() - Builds the cell-resampling draws used by the co-expression bootstrap.
###     gene_reliability() - Per-gene reliability score: the fraction of a gene's variance that is real signal versus sampling noise.
###     row_cor() - Row-wise Pearson correlation between two same-shape matrices.
###     class_mean_var() - Mean and variance of a continuous score within each level of a class vector.
###     frac_group_sets() - Splits scored genes into low/average/high groups by a continuous score, for downstream GO enrichment by group.
###     class_anova() - One-way ANOVA testing whether a continuous score differs across class levels.
###     plot_coexpr_scatter() - Co-expression pair scatter coloured by pair-level class: cis vs trans by regulatory class (Figure 11) or hybrid-vs-parent dominance.
###     seed_compare_core() - Shared scatter/correlation/ratio core for the two two-seed adequacy checks below.
###     gene_seed_compare() - Two-seed adequacy check for the per-gene bootstrap: compares SE between two independent seeds.
###     coexpr_rank_check() - Rank-k eigendecomposition of a divergence matrix, eigenpairs sorted once by |eigenvalue|, with reconstruction R^2 against the observed values.
###   8. REGULATORY AND DOMINANCE CLASSIFICATION  (offset level)
###     classify_reg() - Five-way regulatory classification (Conserved/Cis/Trans/Cis+Trans/Compensatory) from cis and trans p-values.
###     classify_dom() - Six-way dominance classification (hybrid vs each parent) from dpar_sc and dpar_se p-values.
###     clean_reg() - Maps the "Cis x Trans" label to "Compensatory" and sets "Ambiguous" to NA in a class vector.
###     reg_class_vec() - Vector form of classify_reg(), aligned to BURST.CONTRASTS's rows, for a given quantity.
###     dom_class_vec() - Vector form of classify_dom(), aligned to BURST.CONTRASTS's rows, for a given quantity.
###     class_stats() - Applies a per-class summary (class_mean_var, class_anova) to every regulatory and dominance class vector, named REG.MEAN ... DOM.KBAL.
###     ploidy_shift() - Per-gene expected log2 rise in bfreq when HYB.COMB sums two alleles (from allele weights and the intrinsic fraction).
###     .overlap_lor() - Log2((observed + 0.5) / (expected + 0.5)) of a class contingency table; one definition for the heatmap colours and shared_overlap_rng().
###     class_overlap_heatmap() - Log2(observed/expected) association heatmap between two class vectors with BH-adjusted significance; levels auto-derived from the data or fixed so panels share axes (Figures 2 and 6, diagnostics).
###     class_identity_overlap() - Cohen's kappa plus per-class Jaccard overlap between two class vectors on the same genes.
###     summarize_class_overlap() - One printable summary row (n, concordance, kappa, permutation p) per class-overlap comparison.
###     plot_cis_trans_class() - Cis vs trans scatter for one quantity, coloured by class, with SE bars drawn behind points.
###     plot_mean_bfreq_class() - Mean divergence vs a burst-parameter divergence (bfreq/bsize/kbal/cv2) scatter with per-class regulatory and dominance slopes.
###   8b. GO / KEGG enrichment for the regulatory classification
###     build_reg_go_sets() - Builds one GO gene set per (class, direction) pair for a given classification quantity.
###     build_component_go_sets() - Builds pooled any-cis or any-trans GO gene sets, direction-matched on that component's own sign.
###     run_enrichment() - Runs GO (BP/MF/CC, simplified) and KEGG enrichment for one gene set against a fixed universe.
###     axis_pole_enrichment() - GO (BP/MF/CC) or KEGG enrichment of one pole of a co-expression axis against the co-expressed universe.
###     n_sig_terms() - Counts significant terms (q < threshold) in one enrichResult, 0 for a NULL or empty result.
###     print_enrich_brief() - Console view of one enrichResult: Description, p.adjust and Count for the n_top most significant terms.
###     plot_cluster_marker_enrichment() - Writes cluster_marker_enrichment()'s per-cluster up/down enrichment to one pdf (its bar-pair helper is nested inside).
###     score_cell_cycle_by_cluster() - Cell-cycle phase scoring (CellCycleScoring()/AddModuleScore()) extended to a dataset's own validated clustering; violin plot, phase-composition-by-cluster barplot, and console table.
###     go_gene_set() - Pulls all genes (including descendant terms) annotated to a GO term from org.Sc.sgd.db, for building an independently-sourced module-score gene set.
###     score_modules_by_cluster() - Continuous module-score validation (AddModuleScore()) of a cluster's marker-based identity against named, independently-sourced gene sets; violin plot, per-cluster mean +/- SE barplot, and console table.
###     cell_cycle_continuum() - Single continuous cell-cycle axis, first PC of S.Score/G2M.Score.
###     cell_cycle_continuum_shared() - Two-object version of cell_cycle_continuum(), pooled PCA so between-object means are comparable.
###     metabolic_state_cluster() - Discrete metabolic state from joint k-means clustering on the Glycolysis/OXPHOS/RiBi module scores, k chosen by silhouette.
###     covariate_noise_diagnostic() - Per-gene test of whether NB Pearson residual noise depends on cell-cycle position (Spearman) or metabolic state (Kruskal-Wallis, eta-squared).
###     plot_covariate_noise_diagnostic() - Histogram pair of covariate_noise_diagnostic()'s per-gene effect sizes.
###     species_composition_report() - Runs the Cohen's d / rank-sum bound across the cell-cycle axis and the three metabolic module scores for one species/allele pair.
###   9. PUBLICATION FIGURES
###     shared_overlap_rng() - One common log2(obs/exp) colour half-range for a set of overlap heatmaps.
###     plot_dom_class() - Dominance scatter in the parent frame (or A/D rotation), coloured by dominance class.
###     intrinsic_extrinsic_components() - Raw intrinsic and extrinsic noise components (Poisson-corrected) from the two hybrid alleles, per gene.
###   10. PROMOTER ARCHITECTURE (TATA box, poly(dA:dT), nucleosome occupancy)
###     read_genome_fasta() - Reads a genome FASTA into a named DNAStringSet, one sequence per chromosome.
###     read_gff_genes() - Reads a GFF3 annotation into a per-gene coordinate table.
###     extract_promoters() - Extracts each gene's promoter sequence, bounded by the nearest upstream gene end (running maximum over all genes that start before it).
###     score_promoters() - Applies both the TATA PWM score and poly(dA:dT) tract length to every promoter in a set.
###     nupop_cluster_inputs() - Packages one species' chromosomes and promoter coordinates for the NuPoP cluster job.
###     score_promoters_nupop() - Scores each promoter from the cluster occupancy tracks after confirming they match the current promoters.
###     .concordance_eligible() - Genes eligible for the concordance tests: cis-class, both values defined and nonzero (with the .CIS_CLASSES and .PROMOTER_PREDICTED_SIGN constants).
###     promoter_direction_test() - Tests whether a promoter feature's between-species direction matches the predicted direction of cis divergence (burst frequency primary, burst size with the sign flipped).
###     concordance_by_magnitude() - Splits an any-cis gene set into magnitude bins and computes concordance within each bin.
###     plot_concordance_by_magnitude() - Bar plot of concordance_by_magnitude()'s output, one bar per magnitude bin.
###     promoter_noise_candidates() - Lists genes with a cis noise component and a top-decile promoter-feature shift, with each gene's concordance flag, for manual inspection.
###     architecture_noise_check() - Tests whether a promoter feature associates with a species' own DISP or BSIZE, controlling for a second fit column (MU by default).
###   11. SEURAT CLUSTERING DIAGNOSTICS (data-driven feature count, resolution, and metric choices)
###     hvg_elbow() - Elbow point on the ranked standardized-variance curve from FindVariableFeatures(), for a data-driven nfeatures.
###     plot_hvg_elbow() - Diagnostic plot for hvg_elbow()'s output: full curve, chosen cutoff, and floor marked.
###     sweep_cluster_resolution() - Resolution sweep with a minimum-cluster-size guard; picks the coarsest resolution with near-maximal silhouette width.
###     plot_resolution_sweep() - Diagnostic plot for sweep_cluster_resolution()'s output: silhouette vs. resolution, guard-excluded points and the chosen resolution marked.
###     prepare_dataset() - Log-normalization, elbow-selected variable features, scaling and PCA for one Seurat dataset.
###     elbow_pcs() - Percent variance per PC and the elbow PC count to retain.
###     cluster_dataset() - Resolution sweep, plot and log line for one dataset on its retained PCs.
###     umap_dataset() - UMAP on the retained PCs of one clustered dataset, drawn on the open device.
###     assemble_cluster_stability() - Summarizes the returned ARI vectors per candidate resolution, picks the final resolution and relabels the object.
###     kegg_local() - Downloads the KEGG pathway map once, locally, for offline enrichment on cluster nodes.
###     load_cluster_output() - Loads a cluster result into the caller's environment, naming the script to run when the file is missing.
###     check_cluster_key() - Confirms a cluster output was built from the inputs packaged in this session.
###     load_replicates() - Loads the per-replicate outputs of one cluster job array and checks each has the expected size.
###   12. CLUSTER-BASED NOISE PARTITIONING (within/between-cluster variance vs. the intrinsic/extrinsic decomposition)
###     within_between_decomp() - Within- and between-cluster variance per gene, in shot-noise-corrected, mean-normalized rate space, from a Seurat cluster partition.
###     plot_within_between_hist() - Genome-wide distribution of within_between_decomp()'s ratio, one dataset, with the within = between line marked.
###     plot_within_between_vs_quantity() - Scatter of within_between_decomp()'s ratio against one burst kinetics quantity (mean, burst frequency, or burst size), with Spearman rho reported.
###     report_within_between_by_class() - Plots the within/between ratio by class (with class_anova()'s omnibus test and Tukey comparisons) for two datasets side by side against one classification axis; writes the combined figure and prints both datasets' test statistics.
###   13b. FUNCTIONS APPLIED ACROSS ITEMS BY analysis.R
###     qc_cell_cutoff() - library-size cell filter for one dataset s (a list with lib, the per-cell library.
###     fstar_from_r() - closed-form split fraction f* = [(2+r) - sqrt(r^2+4)] / (2r) for r = B*Nh/A.
###     gene_pass_group() - genes passing the filters in one group g: a finite NB dispersion, a mean count at or.
###     perm_label_draw() - one permutation's relabelings for every mode.
###     eiv_mode_ci_row() - errors-in-variables correlation for one mode m with a gene-resampling bootstrap CI.
###     seed_check_rows() - Summary rows (correlation, SE ratio median and IQR, count) of a two-seed bootstrap SE comparison, for the gene-level and co-expression checks.
###     se_floor_row() - count and median bootstrap SE (one column per SE vector) among items whose attenuation is at least floor f.
###     coexpr_seed_check() - two-seed adequacy check for the co-expression bootstrap.
###     coexpr_perm_draw() - one permutation draw for the co-expression null: the pooled-cell order (n_tot cells),.
###     axis_mixture_summary() - variance explained, effective genes, mixture means and sigmas, and pole sizes of.
###     pole_trans_fractions() - for one gene group g, its size and the fraction classified any-trans (Trans or.
###     allele_cor_boot_se_row() - bootstrap SE of row i's allele-residual correlation between sc and se: B.
###     partial_cor_depth() - correlation of gene g's Sc- and Se-allele residuals with cell depth partialled out.
###     to_seurat_counts() - converts a count matrix to a Seurat-ready sparse matrix with dash-delimited feature names.
###     stability_task_rows() - the candidate resolutions to bootstrap for dataset d: the chosen resolution plus the.
###     stability_dataset_inputs() - sparse counts and the fit settings (nfeatures, dims_n, metric) the cluster job.
###     boot_resample_matrix() - draws every bootstrap resample of dataset d up front, one column per replicate,.
###     metric_clusters() - clusters of a Seurat object at one resolution using the annoy distance metric m on the.
###     consistent_pair_row() - keeps Sc-hybrid cluster sc as a pair only when its best Se-hybrid partner picks it.
###     diet_for_markers() - strips a Seurat object to the counts and data layers FindMarkers needs, keeping the cluster job's input small.
###   13c. FUNCTIONS THE CLUSTER SCRIPTS RUN
###     pilot_split_se_one() - one gene's bootstrap SE of log(mu) and log(size) for the four datasets needed to compute A and B, plus the bootstrap.
###     boot_contrasts_one() - one gene's bootstrap contrasts across all modes (the unit of work a cluster worker does).
###     permute_contrasts_one() - one gene's permutation null: refits the relabeled groups for every mode and permutation, compares the observed.
###     coexpr_bootstrap_one() - one bootstrap draw of the co-expression decomposition.
###     coexpr_part_table() - estimate, bootstrap SE (from the running sums), z and two-sided normal p for decomposition nm, one row per gene pair.
###     coexpr_perm_one() - one permutation draw of the five null spectra, for rank-matched testing of every candidate axis.
###     nupop_occupancy_cluster() - Scores every chromosome in one species' inputs and returns each gene's.
###     report_regions() - Summarizes the regions that needed fallback flanks or stayed unscored,.
###     bootstrap_ari_job() - one bootstrap replicate (job j of the tasks x replicates table).
###     dataset_marker_enrichment() - for each cluster of dataset d at its final resolution, marker genes against the rest of the SAME dataset's.
###     go_enrich_job() - one enrichment job k.
###     ora_jobs() - one over-representation job per gene set in sets (named group::set), all against the same universe.
###     gse_job_set() - The rank-based enrichment jobs (BP, MF, CC) for one rank list.
###############################################################

library(MASS)   # glm.nb, ships with base R

## ckpt_path(n, dir): full path of the section n checkpoint. The default dir
## is read when the function is called, so every section names its checkpoint
## with the same section{N}_checkpoint.rda convention.
ckpt_path <- function(n, dir = CHECKPOINT.DIR) {
  file.path(dir, sprintf("section%d_checkpoint.rda", n))
}

## console_start(n, dir, append): copies printed (stdout) output from this point on
## into section{n}_console.txt; split = TRUE keeps the text visible in the R console
## as well. Starting the next section closes the previous file, so each section gets
## its own record. append = TRUE adds to an existing file, which suits re-running a
## single subsection. Any graphics device a failed block left open is closed here, so a
## section starts with none.
console_start <- function(n, dir = CONSOLE.DIR, append = FALSE) {
  console_stop()
  graphics.off()
  sink(file.path(dir, sprintf("section%d_console.txt", n)), append = append, split = TRUE)
}

## with_local_seed(seed, expr): evaluates expr after set.seed(seed) and gives the caller's random number state
## back afterwards (or removes the state when none existed). A function that needs a reproducible stream
## draws it through here, so it neither depends on nor moves any other random draw. The values drawn
## are those of set.seed(seed) followed by expr.
with_local_seed <- function(seed, expr) {
  stopifnot(is.numeric(seed), length(seed) == 1, is.finite(seed))
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()),
          add = TRUE)
  set.seed(seed)
  expr
}

## console_stop(): closes any open console file.
console_stop <- function() {
  while (sink.number() > 0) sink()
  invisible(NULL)
}

## ============================================================
## 1. OFFSET NEGATIVE-BINOMIAL FIT
## ============================================================
# Offset NB fit for one gene. Counts follow y ~ NB(mean = exposure * mu, size = disp),
# so mu is a per-unit-exposure rate and disp is the NB size theta (larger = less noise).
# The point estimates come from .fit_one(): mu = sum(y)/sum(exposure), and disp is the
# NB likelihood maximizer given that mu (Inf when Pearson chi-square / df <= 1, i.e. no
# detectable overdispersion). Every bootstrap, permutation and power-grid replicate fits
# with this same estimator, so an observed value and the distribution it is compared
# against are measured identically. glm.nb supplies only the asymptotic log-scale SEs.
# Returns c(disp, mu, disp_logse, mu_logse).
neg_binom_fit_offset <- function(y, exposure) {
  if (!is.null(dim(y)))
    stop("neg_binom_fit_offset() takes one gene's counts; use fit_counts_offset() for a matrix")
  na_out <- c(disp = NA_real_, mu = NA_real_, disp_logse = NA_real_, mu_logse = NA_real_)
  if (length(y) < 2 || length(exposure) != length(y)) return(na_out)
  keep <- is.finite(y) & is.finite(exposure) & exposure > 0
  y <- y[keep]; exposure <- exposure[keep]
  if (length(y) < 2) return(na_out)

  core <- .fit_one(y, exposure)
  out  <- c(disp = unname(core["disp"]), mu = unname(core["mu"]), disp_logse = NA_real_, mu_logse = NA_real_)
  if (!is.finite(out[["disp"]])) return(out)

  logexp <- log(exposure)
  fit <- suppressWarnings(tryCatch(glm.nb(y ~ 1 + offset(logexp), init.theta = out[["disp"]]),
                                   error = function(e) NULL))
  if (!is.null(fit) && !isFALSE(fit$converged)) {
    out["mu_logse"] <- tryCatch(sqrt(vcov(fit)[1, 1]), error = function(e) NA_real_)
    if (is.finite(fit$SE.theta) && is.finite(fit$theta) && fit$theta > 0)
      out["disp_logse"] <- fit$SE.theta / fit$theta
  }
  out
}

# Bounds of the NB size (theta) estimate. The 1-D likelihood search runs over log(theta) in
# LOG.THETA.RANGE. A fitted theta above THETA.CAP is indistinguishable from Poisson at any count the data
# can resolve, so it is reported as Inf (the Poisson limit), and gene filters treat DISP >= THETA.CAP as
# not identifiable.
THETA.CAP       <- 1e6
LOG.THETA.RANGE <- c(-4, 15)

# The pipeline's single NB estimator: mu = sum(y)/sum(exposure), and disp = the NB size
# from a 1-D likelihood maximization over log(theta) given that mu, after a Poisson-vs-NB
# pre-check (Pearson chi-square / df <= 1 gives disp = Inf). neg_binom_fit_offset() uses it
# for the observed fits, gene_boot.R and gene_perm.R for every resample,
# and power_grid.R for the simulated observed contrast and its null. Returns c(mu, disp);
# base R only and RNG-free, so fork-safe.
.fit_one <- function(counts, expo) {
  keep <- is.finite(counts) & is.finite(expo) & expo > 0
  y <- counts[keep]; e <- expo[keep]
  if (length(y) < 2) return(c(mu = NA_real_, disp = NA_real_))

  mu_hat <- sum(y) / sum(e)
  if (sum(y) == 0) return(c(mu = 0, disp = NA_real_))

  mu_i <- mu_hat * e
  pearson <- sum((y - mu_i)^2 / mu_i) / (length(y) - 1)
  if (pearson <= 1) return(c(mu = mu_hat, disp = Inf))

  # NB log-likelihood in log(theta), with the means fixed at mu_i
  opt <- tryCatch(optimize(function(ltheta) {
    th <- exp(ltheta)
    sum(lgamma(y + th) - lgamma(th) + th*log(th) - (th + y)*log(th + mu_i) + y*log(mu_i))
  }, LOG.THETA.RANGE, maximum = TRUE), error = function(e) NULL)
  if (is.null(opt)) return(c(mu = mu_hat, disp = NA_real_))

  theta <- exp(opt$maximum)
  c(mu = mu_hat, disp = if (theta > THETA.CAP) Inf else theta)
}

# Splits a pooled sample at position n1 of a pre-drawn index permutation and fits each
# half with .fit_one(): the permutation null for a two-group contrast (labels
# are exchangeable under H0). Used by the gene-level null and the power grid. Returns
# list(a = first n1 cells, b = remaining cells).
.fit_split <- function(counts, expo, perm, n1) {
  stopifnot(n1 >= 1, n1 <= length(perm))
  g1 <- perm[seq_len(n1)]; g0 <- perm[seq_along(perm) > n1]   # empty second group when n1 == length(perm)
  list(a = .fit_one(counts[g1], expo[g1]), b = .fit_one(counts[g0], expo[g0]))
}

## Fits every gene (row) of a genes x cells matrix against one shared exposure vector
## with neg_binom_fit_offset() and returns a per-gene table: DISP (NB size), MU, their log
## SEs, MEAN_CT, N_EXPR, and derived VAR, FANO, CV, BFREQ (= DISP) and BSIZE (= MU / DISP).
## cl is optional: NULL fits on one core; a PSOCK cluster from parallel::makeCluster() splits
## the genes into one chunk per worker, each chunk carrying all cells for its genes. A PSOCK
## worker starts with an empty workspace, so the functions the fit needs are sent to it here:
## fit_counts_offset() and everything it reaches are found from the code itself
## (codetools::findGlobals), so a new dependency of the fit reaches the workers without
## editing a list. Genes are fit independently, so chunks need no communication.
fit_counts_offset <- function(mat, exposure, cl = NULL) {
  stopifnot(ncol(mat) == length(exposure))
  if (is.null(cl))
    return(do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
      ## One NB fit per gene. DISP is the NB size parameter (theta), which this pipeline also reports
      ## as burst frequency (BFREQ = DISP). BSIZE is the separate derived quantity MU / DISP (mean
      ## count per burst); DISP and BSIZE are distinct quantities.
      y <- mat[i, ]
      f <- neg_binom_fit_offset(y, exposure)
      disp <- f[["disp"]]; mu <- f[["mu"]]
      data.frame(
        DISP = disp, MU = mu, DISP_LOGSE = f[["disp_logse"]], MU_LOGSE = f[["mu_logse"]],
        MEAN_CT = sum(y) / length(y), N_EXPR = sum(y > 0),
        VAR = mu + mu^2 / disp, FANO = 1 + mu / disp, CV = sqrt(1 / disp + 1 / mu),
        BFREQ = disp, BSIZE = mu / disp,
        row.names = rownames(mat)[i])
    })))

  needed <- character(); todo <- "fit_counts_offset"
  while (length(todo)) {
    f <- todo[1]; todo <- todo[-1]
    if (f %in% needed) next
    needed <- c(needed, f)
    todo   <- c(todo, intersect(codetools::findGlobals(get(f, envir = globalenv()), merge = TRUE), ls(globalenv(), all.names = TRUE)))
  }
  chunk_id   <- cut(seq_len(nrow(mat)), length(cl), labels = FALSE)
  mat_chunks <- split.data.frame(mat, chunk_id)
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS)))
  parallel::clusterExport(cl, needed, envir = globalenv())
  fit_chunks <- parallel::parLapply(cl, mat_chunks, fit_counts_offset, exposure = exposure)
  ## unname() before rbind() keeps the gene row names unprefixed by the chunk labels.
  do.call(rbind, unname(fit_chunks))
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

## Depth-matched split of depth-ordered hybrid cells at fraction f (0 < f <= 0.5).
## The c group takes round(f * n) evenly spaced positions along the depth order and
## the t group takes the rest, so both groups span the full depth range (f = 0.5 is
## strict every-other-cell alternation). Contiguous blocks of a depth-sorted order
## would confound the split with depth. Returns list(c, t) of cell positions.
split_indices_by_depth <- function(n, f) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 2, is.numeric(f), length(f) == 1, f > 0, f < 1)
  n_c <- max(1, round(f * n))
  c_pos <- unique(round(seq(1, n, length.out = n_c)))
  list(c = c_pos, t = setdiff(seq_len(n), c_pos))
}

## ============================================================
## 2. GENE FILTER (marginal information only; never on a contrast)
## ============================================================

## ---- Adaptive QC thresholds ----
## Cutoffs are defined in relative units so that datasets with different
## depth and cell number lose the same kind of cell or gene, which keeps
## the retained data comparable.

## Gene QC across all datasets at once. Abundance is a fraction of each
## dataset's total reads, so every dataset is judged on the same relative
## scale. lambda0 is the mean count per cell required at the shallowest
## dataset; deeper datasets face a proportionally higher mean-count floor.
## The detection floor is an absolute number of cells because the precision
## of a dispersion fit depends on informative cells, not on their fraction.
qc_gene_keep <- function(mats, lambda0 = 0.20, cell_frac = 0.10) {
  stopifnot(is.list(mats), length(mats) >= 1, is.numeric(lambda0), length(lambda0) == 1, lambda0 > 0,
            is.numeric(cell_frac), length(cell_frac) == 1, cell_frac > 0, cell_frac <= 1)
  genes <- Reduce(intersect, lapply(mats, rownames))
  mats  <- lapply(mats, `[`, genes, , drop = FALSE)
  depth <- vapply(mats, sum, numeric(1)) / vapply(mats, ncol, integer(1))
  p_min <- lambda0 / min(depth)
  n_min <- ceiling(cell_frac * min(vapply(mats, ncol, integer(1))))
  pass  <- vector("list", length(mats)); names(pass) <- names(mats)
  for (i in seq_along(mats)) pass[[i]] <- rowSums(mats[[i]]) / sum(mats[[i]]) >= p_min & rowSums(mats[[i]] > 0) >= n_min
  list(genes = sort(genes[Reduce(`&`, pass)]), p_min = p_min, n_min = n_min,
       depth = depth, exp_mean = p_min * depth)
}

## ============================================================
## 3. PAIRED PER-GENE BOOTSTRAP OF CONTRASTS
## ============================================================

## Draws all B sets of bootstrap resampling indices up front, from one seeded stream, so the
## replicates are independent of any RNG use inside the fitting code. The two parents are
## resampled as whole datasets (MIX.SC, MIX.SE).
##
## The hybrid cells are cut into four strata, the 2 x 2 overlap of the mean split (HYC / HYT, cells
## hyc and the rest) with the noise split (HYC.N / HYT.N, cells hyc_n and the rest). Each replicate
## resamples the cells of every stratum with replacement and keeps the stratum size, then builds every
## hybrid dataset from those same resampled cell IDs: HYC = strata (C, C.N) + (C, T.N), HYT =
## (T, C.N) + (T, T.N), HYC.N = (C, C.N) + (T, C.N), HYT.N = (C, T.N) + (T, T.N), and HYB = all four.
## Group sizes stay fixed, the SC and SE alleles of a cell stay paired (one index vector serves both),
## and every overlap between the splits is preserved, so the mean-split and noise-split contrasts of a
## mode are resampled from the same cells and their bootstrap correlation is measured. The strata are
## depth matched (split_indices_by_depth), so the marginal SEs differ little from resampling each
## group on its own.
##
## hyc and hyc_n are the hybrid-cell positions (columns of HYB.SC) of the HYC and HYC.N groups.
## Each returned hybrid index vector gives positions within its own dataset's columns, the layout
## boot_contrasts_one() expects; HYB indexes the full hybrid column order.
make_draws <- function(ncells, B, seed = 1, hyc, hyc_n) {
  n_h <- ncells[["HYB.COMB"]]
  stopifnot(is.numeric(B), length(B) == 1, B >= 1,
            is.numeric(hyc), is.numeric(hyc_n),
            !anyDuplicated(hyc), !anyDuplicated(hyc_n),
            all(hyc %in% seq_len(n_h)), all(hyc_n %in% seq_len(n_h)),
            length(hyc) < n_h, length(hyc_n) < n_h)
  ## Group membership in the dataset's own column order (the order the splits were cut in)
  members <- list(HYC = sort(hyc), HYT = setdiff(seq_len(n_h), hyc),
                  HYC.N = sort(hyc_n), HYT.N = setdiff(seq_len(n_h), hyc_n))
  ## Stratum code: 1 = (C, C.N), 2 = (T, C.N), 3 = (C, T.N), 4 = (T, T.N)
  code   <- 1L + as.integer(!(seq_len(n_h) %in% hyc)) + 2L * as.integer(!(seq_len(n_h) %in% hyc_n))
  strata <- split(seq_len(n_h), factor(code, levels = 1:4))
  of     <- list(HYC = c(1, 3), HYT = c(2, 4), HYC.N = c(1, 2), HYT.N = c(3, 4))
  with_local_seed(seed, lapply(seq_len(B), function(b) {
    ## Resampled cell IDs per stratum; an empty stratum stays empty
    res <- lapply(strata, function(s) s[sample.int(length(s), length(s), replace = TRUE)])
    d <- list(
      MIX.SC = sample.int(ncells[["MIX.SC"]], ncells[["MIX.SC"]], replace = TRUE),
      MIX.SE = sample.int(ncells[["MIX.SE"]], ncells[["MIX.SE"]], replace = TRUE))
    for (k in names(of)) d[[k]] <- match(unlist(res[of[[k]]], use.names = FALSE), members[[k]])
    d$HYB <- unlist(res, use.names = FALSE)
    d
  }))
}

## Contrast definitions. Each mode lists the groups it needs (for the
## finite-input gate) and a formula on a named list of fitted values.
## Divergence: total, cis, trans. Dominance, midparent form: dom
## (arithmetic midparent for mean, geometric for burst frequency).
## Classic hybrid-vs-each-parent: dpar_sc, dpar_se. Per-allele
## inheritance: inh_sc, inh_se. Add an axis by adding one row here.
##
## HAPLOID PARENTS AND THE DIPLOID HYBRID
## Every fit offsets by the cell's total library (EXPO), so MU reports a
## share of that library. A hybrid cell carries both genomes, so its library
## and its per-gene counts scale together and the share stays comparable to
## a haploid parent's share. HYB.COMB sums the two alleles. HYB.SC and
## HYB.SE keep them apart.
## Mean axis
##   total, cis  Ratios inside one ploidy state, so a shared factor cancels.
##   trans       (parent ratio) - (hybrid allele ratio). The factor cancels
##               inside each ratio, so only a species-specific ploidy effect
##               survives the subtraction.
##   dpar, dom   An additive gene holds the midparent share of the hybrid
##               library. dom therefore uses the arithmetic midparent, and
##               its null sums parent cells and their libraries to mirror
##               the hybrid.
##   inh         One allele against the full diploid library, which places
##               pure inheritance near -1 log2. inh serves as a descriptive
##               contrast and no classifier reads it.
## Noise axis (DISP is the NB size and has no summing rule)
##   total, cis, trans  Ratios again, so a ploidy effect shared by both
##               species (cell size, cell-cycle timing) cancels.
##   dom         The midparent is geometric, the mean of the two log2 values.
##               Its observed noise contrast and null keep the raw HYB.COMB
##               scale, so the noise axis of dom stays descriptive.
##   dpar        HYB.COMB sums two alleles per cell, which averages their
##               private (intrinsic) noise and raises DISP in the hybrid.
##               For alleles with mean weights w1 and w2 and intrinsic
##               fraction phi, the summed pair keeps the shared noise and
##               scales the private noise by w1^2 + w2^2. The expected rise
##               in log2 DISP is s = -log2(1 - (1 - w1^2 - w2^2) * phi),
##               between 0 and 1 (ploidy_shift). ploidy_adjust_dpar removes
##               s from the dpar estimates as follows. bfreq falls by s,
##               bsize (mean - bfreq) rises by s, kbal (bfreq - bsize)
##               falls by 2s, and cv2 recomputes with the latent term scaled
##               by 2^s. The mean axis needs no adjustment. The bootstrap SEs
##               stay as fitted, so the shift enters as a fixed per-gene
##               constant. The permutation job reads the adjusted contrast
##               against the zero-centred exchangeability null and reports
##               _p_ploidy (own shift) and _p_ind (shift of 1, fully
##               independent equal alleles) beside the raw _p.
##   inh         One allele against one haploid genome, so no averaging
##               enters and the noise axis of inh compares like with like.
##
## cis_n and trans_n are internal only, not reported directly. They
## are cis and trans computed from the HYC.N/HYT.N split (chosen to
## balance noise-axis power) rather than the HYC/HYT split (chosen to
## balance mean-axis power). A single output column pair, bfreq_cis
## and bfreq_trans, draws its value from these instead of from cis and
## trans, since the noise axis needs its own split, not the mean
## axis's split reused. See .BFREQ_SOURCE and boot_contrasts_one() in gene_boot.R.
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

## Maps each dataset to the draws-list element (see make_draws()) that resamples it. The four hybrid
## keys (HYC, HYT, HYC.N, HYT.N) and HYB are built from the same resampled hybrid cell IDs, so every
## mode shares its resampled cells with every other mode that uses the same cells.
.DRAW.KEY <- c(MIX.SC = "MIX.SC", MIX.SE = "MIX.SE",
               HYC.SC = "HYC", HYC.SE = "HYC", HYT.SC = "HYT", HYT.SE = "HYT",
               HYC.SC.N = "HYC.N", HYC.SE.N = "HYC.N", HYT.SC.N = "HYT.N", HYT.SE.N = "HYT.N",
               HYB.COMB = "HYB", HYB.SC = "HYB", HYB.SE = "HYB")

## Log2 contrast for one mode from a named list gv of per-group fitted values. q selects the
## quantity: "MU" (library share; arithmetic midparent for dom), "BFREQ" (NB size; geometric
## midparent) or "CV2". total and cis are ratios within one ploidy state; trans is the parental
## ratio minus the hybrid allele ratio.
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

## CV2 (squared coefficient of variation) of one fitted group: 1/mu + 1/disp, the Poisson plus
## overdispersion terms. Returns NA unless mu and disp are finite and positive, so the bootstrap
## SD and the permutation null see only valid values. The single CV2 formula for the bootstrap,
## the permutation null and the observed fits.
.cv2_of <- function(mu, disp) {
  if (is.finite(mu) && mu > 0 && is.finite(disp) && disp > 0) 1 / mu + 1 / disp else NA_real_
}

## Adds burst size (bsize = mean - bfreq) and kinetic balance (kbal = bfreq - bsize = 2*bfreq - mean),
## both on the log2 scale, with propagated SEs, to the finished bootstrap table. bfreq is the NB
## dispersion contrast under its biological name. Because bsize is defined by that subtraction,
## bfreq and bsize are structurally anti-correlated whenever bfreq varies more than the mean
## (see the structural bfreq-bsize check in analysis.R Section 3.6). kbal together with the mean carries the same information as
## bfreq and bsize without that dependence, so kbal against the mean asks whether a gene's overall
## change leans toward frequency or size. Fitting code keeps the field name "disp".
add_burst_contrasts <- function(df) {
  ## Loops over the reportable modes; BOOT.CONTRASTS has columns only for .OUT_MODES (cis_n and
  ## trans_n are folded into cis and trans through .BFREQ_SOURCE).
  for (ct in .OUT_MODES) {
    m  <- df[[paste0("mean_",ct,"_est")]];  sm <- df[[paste0("mean_",ct,"_se")]]
    s  <- df[[paste0("bfreq_",ct,"_est")]]; ss <- df[[paste0("bfreq_",ct,"_se")]]
    r  <- df[[paste0("cor_",ct)]]
    ## cor_<ct> is measured from the bootstrap draws for every mode; a gene with too few finite draws
    ## (NA) enters the SE propagation as a zero covariance.
    r[!is.finite(r)] <- 0
    bs    <- m - s
    bs_se <- sqrt(pmax(0, sm^2 + ss^2 - 2 * r * sm * ss))
    df[[paste0("bsize_",ct,"_est")]] <- bs
    df[[paste0("bsize_",ct,"_se")]]  <- bs_se
    ## covFB = Cov(bfreq_est, bsize_est) = r*sm*ss - ss^2, the covariance of the two estimation
    ## errors. It gives kbal = bfreq - bsize its SE for every mode.
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

## Errors-in-variables summary of mean vs bfreq divergence across genes for one mode. Subtracting the
## mean squared bootstrap SE from the raw variances and covariance removes the measurement-noise
## attenuation, so rho_mean_disp estimates the correlation of the underlying contrasts.
eiv_components <- function(BURST.CONTRASTS, mode = .OUT_MODES) {
  mode <- match.arg(mode)
  X  <- BURST.CONTRASTS[[paste0("mean_",mode,"_est")]]; sx <- BURST.CONTRASTS[[paste0("mean_",mode,"_se")]]
  Y  <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_est")]]; sy <- BURST.CONTRASTS[[paste0("bfreq_",mode,"_se")]]
  r  <- BURST.CONTRASTS[[paste0("cor_",mode)]]
  ## cor_<mode> is measured from the bootstrap draws; a gene with too few finite draws (NA) enters as a
  ## zero covariance so it stays in the summary.
  r[!is.finite(r)] <- 0
  ok <- is.finite(X) & is.finite(Y) & is.finite(sx) & is.finite(sy)
  X <- X[ok]; Y <- Y[ok]; sx <- sx[ok]; sy <- sy[ok]; r <- r[ok]
  Vm  <- var(X)    - mean(sx^2)
  Vs  <- var(Y)    - mean(sy^2)
  Cms <- cov(X, Y) - mean(r * sx * sy)
  ## Vm or Vs below zero (measurement noise exceeds signal) makes sqrt(Vm * Vs) NaN. The NaN is the
  ## result for that mode (attenuation correction not identifiable), so the warning is suppressed.
  rho_ms <- suppressWarnings(Cms / sqrt(Vm * Vs))
  c(n = length(X), Vm = Vm, Vs = Vs, Cms = Cms, rho_mean_disp = rho_ms, rho_raw_mean_disp = cor(X, Y))
}

## ============================================================
## 5. PLOTTING HELPERS AND GENE-LEVEL PLOTS
## ============================================================

## ---- Shared plotting helpers: color ramps and UMAP ----
## The multi-panel layout and power-curve helpers (open_grid_pdf, plot_lines, legend_page) live in
## functions_power.R with the rest of the power analysis.

## UMAP with the shared cluster palette and on-plot labels as a second cue,
## so identity never depends on color alone. group.by = NULL uses Idents(obj).
umap_plot <- function(obj, title = NULL, group.by = NULL) {
  ids  <- if (is.null(group.by)) Idents(obj) else obj[[group.by, drop = TRUE]]
  lv   <- if (is.factor(ids)) levels(droplevels(ids)) else sort(unique(as.character(ids)))
  cols <- setNames(colorRampPalette(COLOR.CLUSTER)(length(lv)), lv)
  DimPlot(obj, reduction = "umap", group.by = group.by, cols = cols, label = TRUE, repel = TRUE) + ggtitle(title)
}

## Labels of the burst quantities, shared by the figure functions below.
QUANTITY.LABEL <- c(mean = "mean", bfreq = "burst frequency", bsize = "burst size", kbal = "frequency-size balance", cv2 = "CV2")

## fig_pdf: opens a PDF at path (relative to FIGURE.DIR) with the pipeline's figure settings.
fig_pdf <- function(path, width, height, useDingbats = FALSE) {
  pdf(file.path(FIGURE.DIR, path), width = width, height = height, useDingbats = useDingbats)
}

## class_count_barplot: number of genes (or pairs, via ylab) in each class of cls, bars colored by class.
class_count_barplot <- function(cls, levels, cols, main, ylab = "# of genes") {
  barplot(table(factor(cls, levels = levels)), col = cols, las = 2, ylab = ylab, border = NA, main = main)
}

## class_overlap_triptych: figure of three class-overlap heatmaps, mean class (rows) against burst frequency,
## burst size and frequency-size balance class (columns). classes is a list of class vectors named by quantity
## (REG.VEC or DOM.VEC), levels the class levels. The three heatmaps share one color range
## (shared_overlap_rng), so enrichment strength reads the same across panels. Cell text is the fold
## enrichment (obs/exp), * marks BH q < OVERLAP.FDR, and color is log2(obs/exp), diverging through white.
class_overlap_triptych <- function(path, classes, levels, width, height, mar) {
  others <- c("bfreq", "bsize", "kbal")
  pairs <- vector("list", length(others))
  for (i in seq_along(others)) pairs[[i]] <- list(classes$mean, classes[[others[i]]])
  rng <- shared_overlap_rng(pairs, levels)
  fig_pdf(path, width, height)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 3), mar = mar)
  for (o in others)
    class_overlap_heatmap(classes$mean, classes[[o]], levels_a = levels, levels_b = levels,
                          xlab = paste(QUANTITY.LABEL[[o]], "class"), ylab = "mean class", rng = rng)
}

## class_heatmap_grid: figure of class-overlap heatmaps for the pairs of class vectors listed in panels
## (columns y_kind, y_q, x_kind, x_q: kind is the name in classes, e.g. REG or DOM, and q the quantity). Regulatory
## vectors pass through clean_reg() so "Cis x Trans" counts as Compensatory and Ambiguous drops. Axis labels read
## "<quantity> <kind word> class<suffix>"; a kind word of "" leaves it out.
class_heatmap_grid <- function(path, panels, classes, fdr, width, height, mfrow,
                               kind_words = c(REG = "regulatory", DOM = "dominance"), suffix = "") {
  fig_pdf(path, width, height)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = mfrow)
  for (i in seq_len(nrow(panels))) {
    pn <- panels[i, ]
    ## Class vector and axis label of the y and x side of this panel.
    side <- list(y = list(), x = list())
    for (a in names(side)) {
      kind <- pn[[paste0(a, "_kind")]]; q <- pn[[paste0(a, "_q")]]
      v <- classes[[kind]][[q]]
      side[[a]]$vec <- if (kind == "REG") clean_reg(v) else v
      side[[a]]$lab <- paste0(QUANTITY.LABEL[[q]], if (nzchar(kind_words[[kind]])) paste0(" ", kind_words[[kind]]), " class", suffix)
    }
    class_overlap_heatmap(side$y$vec, side$x$vec, fdr = fdr, ylab = side$y$lab, xlab = side$x$lab)
  }
}

## ---- Gene-level scatters and histograms (square symmetric panels, SE bars) ----

## .se_scatter(): square scatter of y against x, the shared body of plot_contrast_scatter(), plot_cis_trans_class(),
## plot_dom_class() and plot_coexpr_scatter(). sx, sy are SEs drawn as bars behind the points (NULL draws none).
## bar_col is one colour or list(x, y) with one colour per bar on each axis (already restricted to the plotted points); pt_col is one colour or one per point. legend_labels with legend_cols adds a top-left legend. With diagonals = TRUE both axes
## share one symmetric range (unless lim is given) and the dotted +/-45 degree lines are drawn, so
## same-direction and opposite-direction changes are readable. With diagonals = FALSE each axis gets its
## own symmetric range. Only points with a finite x, y and (when given) SEs are drawn.
.se_scatter <- function(x, sx = NULL, y, sy = NULL, xlab, ylab, main, lim = NULL, diagonals = TRUE, bar_col = NULL, pt_col,
                        pt_cex = 0.5, legend_labels = NULL, legend_cols = NULL) {
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(sx)) ok <- ok & is.finite(sx) & is.finite(sy)
  x <- x[ok]; y <- y[ok]
  if (length(pt_col) == length(ok)) pt_col <- pt_col[ok]
  bars <- !is.null(sx)
  if (bars) { sx <- sx[ok]; sy <- sy[ok] } else sx <- sy <- numeric(length(x))
  if (diagonals) {
    if (is.null(lim)) { m <- max(abs(c(x + sx, x - sx, y + sy, y - sy)), na.rm = TRUE); lim <- c(-m, m) }
    xlim <- ylim <- lim
  } else {
    xlim <- c(-1, 1) * max(abs(c(x + sx, x - sx)), na.rm = TRUE)   # symmetric about zero, covering value +/- error
    ylim <- c(-1, 1) * max(abs(c(y + sy, y - sy)), na.rm = TRUE)
  }
  op <- par(pty = "s"); on.exit(par(op))
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  if (diagonals) { abline(0, 1, lty = 3, col = COLOR.GREY[["mid"]]); abline(0, -1, lty = 3, col = COLOR.GREY[["mid"]]) }
  abline(h = 0, v = 0, col = COLOR.GREY[["dark"]])
  if (bars) {
    segments(x - sx, y, x + sx, y, col = if (is.list(bar_col)) bar_col$x else bar_col)
    segments(x, y - sy, x, y + sy, col = if (is.list(bar_col)) bar_col$y else bar_col)
  }
  points(x, y, pch = 16, cex = pt_cex, col = pt_col)
  if (!is.null(legend_labels)) legend("topleft", legend = legend_labels, col = legend_cols, pch = 16, bty = "n", cex = 0.8)
}

## plot_contrast_scatter: one gene-level contrast scatter with SE bars, drawn by .se_scatter(). type picks
## the pairing and what its argument:
##   "cis_trans"      what is a quantity ("mean", "bfreq", "bsize", "kbal", "cv2"): cis (x) vs trans (y) on one
##                    symmetric range, so the dotted +/-45 degree lines (same direction, compensatory) are
##                    meaningful. lim overrides that range.
##   "mean_bfreq"     what is a mode: mean (x) vs burst-frequency (NB dispersion, y) contrast on one symmetric
##                    range. The tilt of the cloud is not the true slope; the attenuation-corrected coupling is
##                    in the errors-in-variables table of analysis.R Section 3.
##   "burst_kinetics" what is a mode: rotated coordinates that separate net mean change from the
##                    frequency/size balance: x = net mean change (mean contrast; bfreq + bsize equals it
##                    exactly), y = kinetic balance (bfreq - bsize; right of zero is frequency-led, below is
##                    amplitude-led). Each axis gets its own symmetric range. Both columns carry SEs
##                    propagated in add_burst_contrasts().
plot_contrast_scatter <- function(BURST.CONTRASTS, type = c("cis_trans", "mean_bfreq", "burst_kinetics"), what, main = NULL, lim = NULL,
                                  bar_col = adjustcolor(COLOR.GREY[["dark"]], 0.33), pt_col = "black") {
  type <- match.arg(type)
  diagonals <- type != "burst_kinetics"
  if (type == "cis_trans") {
    quantity <- match.arg(what, c("mean", "bfreq", "bsize", "kbal", "cv2"))
    lab <- QUANTITY.LABEL[[quantity]]
    if (is.null(main)) main <- lab
    cols <- paste0(quantity, c("_cis_est", "_cis_se", "_trans_est", "_trans_se"))
    xlab <- paste(lab, "cis (log2)"); ylab <- paste(lab, "trans (log2)")
  } else {
    mode <- match.arg(what, .OUT_MODES)
    lim <- NULL
    if (type == "mean_bfreq") {
      if (is.null(main)) main <- paste0("mean vs noise: ", mode)
      cols <- paste0(c("mean_", "mean_", "bfreq_", "bfreq_"), mode, c("_est", "_se", "_est", "_se"))
      xlab <- "mean (log2)"; ylab <- "dispersion (log2)"
    } else {
      if (is.null(main)) main <- paste0("burst kinetics: ", mode)
      cols <- paste0(c("mean_", "mean_", "kbal_", "kbal_"), mode, c("_est", "_se", "_est", "_se"))
      xlab <- "net mean change (log2)"; ylab <- "frequency - amplitude (kinetic balance, log2)"
    }
  }
  v <- BURST.CONTRASTS[cols]
  .se_scatter(v[[1]], v[[2]], v[[3]], v[[4]], xlab = xlab, ylab = ylab, main = main, lim = lim, diagonals = diagonals,
              bar_col = bar_col, pt_col = pt_col)
}

## sig_hist_panel: histogram of one contrast (quantity, mode) with the genes significant at q < sig shaded by
## direction (Sc-higher up, Se-higher dn). Breaks span the full range of the estimate, rounded outward to a
## multiple of brk, so every gene is counted and bin edges stay on one grid; the x window is only the visible
## range. Burst frequency spans a narrower range and uses finer bins. ymax is the top of the y axis.
sig_hist_panel <- function(contrasts, pr, mode, quantity, ymax, sig = 0.05, up = SPECIES.COLOR[["Sc"]], dn = SPECIES.COLOR[["Se"]]) {
  x <- contrasts[[paste0(quantity, "_", mode, "_est")]]
  p <- pr[[paste0(quantity, "_", mode, "_q")]]
  fine <- quantity == "bfreq"
  brk  <- if (fine) 0.05 else 0.1
  xlim <- if (fine) c(-2.5, 2.5) else c(-5, 5)
  ylim <- c(0, ymax)
  xlab <- paste0(c(total = "parents", cis = "cis", trans = "trans")[[mode]], " log2(Sc/Se) ", QUANTITY.LABEL[[quantity]])
  stopifnot(any(is.finite(x)))
  rng <- range(x[is.finite(x)])    # an Inf contrast must not set the bin range
  lo  <- floor((rng[1] - 1e-9) / brk) * brk
  hi  <- ceiling((rng[2] + 1e-9) / brk) * brk
  b   <- seq(lo, hi, by = brk)
  hist(x, breaks = b, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "# of Genes", main = "", col = COLOR.GREY[["light"]], border = "white")
  hist(x[p < sig & x > 0], breaks = b, xlim = xlim, ylim = ylim, col = up, border = "white", add = TRUE)
  hist(x[p < sig & x < 0], breaks = b, xlim = xlim, ylim = ylim, col = dn, border = "white", add = TRUE)
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
## Centering of each null
##   total, cis, trans  Center at zero on both axes.
##   dpar        Pools hybrid and parent cells, so it tests exchangeability
##               and centers at zero. On the mean axis the offset keeps an
##               additive gene at zero. On the noise axis the observed
##               contrast carries the allele averaging shift, so the
##               permutation reads the ploidy-adjusted contrast against the
##               same null and reports _p_ploidy and _p_ind (see .MODES).
##   dom         The surrogate sums independent parent cells. The mean axis
##               centers near zero. The noise axis centers near the
##               averaging shift (up to +1 log2 for equal parents), so the
##               dom p-value is descriptive and the classifiers read dpar.
##   inh         Pools a half-share hybrid allele with a full-share parent,
##               so the mean axis carries the -1 log2 baseline and the
##               p-value is descriptive.
## Downstream classifiers read the BH-adjusted _q columns (fdr_columns).

## perm_pval: two-sided permutation p-value on |statistic| with an add-one correction, so a
## finite null never yields p = 0; non-finite null draws are dropped. Shared by the gene-level
## permutation null (permute_contrasts_one() in gene_perm.R) and the power grid (power_grid_row() in power_grid.R).
perm_pval <- function(obs, null) {
  ok <- is.finite(null)
  if (!is.finite(obs) || sum(ok) < 1) return(NA_real_)
  (1 + sum(abs(null[ok]) >= abs(obs))) / (1 + sum(ok))
}

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

## Pearson residuals (y - mu*e) / sqrt(mu*e + (mu*e)^2 / k) from the offset NB fit, genes x cells.
## The exposure offset removes depth and the NB SD puts genes on a common variance scale, so
## correlations between residuals reflect co-fluctuation. Non-finite entries (zero fitted mean) are set to 0.
nb_residuals <- function(mat, exposure, fit) {
  fitted <- outer(fit$MU, exposure)
  v <- fitted + ifelse(is.finite(fit$DISP), fitted^2 / fit$DISP, 0)
  r <- (mat - fitted) / sqrt(v)
  r[!is.finite(r)] <- 0
  r
}

## Analytic shrinkage of a correlation matrix toward the identity. Z is cells x genes.
## Off-diagonals are scaled by (1 - lambda), where lambda = sum Var(r_ij) / sum r_ij^2 is the
## Schaefer-Strimmer optimal intensity clipped to [0, 1]: noisier correlations contract more.
## The diagonal stays 1 and lambda is returned as attr(, "lambda").
shrink_cor <- function(Z) {
  stopifnot(is.matrix(Z), is.numeric(Z))
  Z <- Z[stats::complete.cases(Z), , drop = FALSE]
  n <- nrow(Z); p <- ncol(Z)
  if (p < 2 || n < 3) stop(sprintf("shrink_cor(): needs at least 2 genes and 3 complete cells, got %d genes and %d cells", p, n))
  Zs <- scale(Z); Zs[!is.finite(Zs)] <- 0
  R  <- crossprod(Zs) / (n - 1)
  ## Var(r_ij) summed over the pairs i < j, in closed form. For w = z_i * z_j over the n cells,
  ## sum(w^2) is crossprod(Zs^2)[i, j] and mean(w) = R_ij * (n - 1) / n, so
  ## sum((w - mean(w))^2) = crossprod(Zs^2)[i, j] - (n - 1)^2 / n * R_ij^2, with no loop over pairs.
  up  <- upper.tri(R)
  num <- sum(n / (n - 1)^3 * (crossprod(Zs^2)[up] - (n - 1)^2 / n * R[up]^2))
  den <- sum(R[up]^2)
  lam <- if (den > 0) max(0, min(1, num / den)) else 1
  Rs <- R * (1 - lam); diag(Rs) <- 1
  attr(Rs, "lambda") <- lam
  Rs
}

## Total, cis and trans divergence of the pairwise residual-correlation structure, plus the
## dominance contrasts when HYB.COMB is supplied. resid: list of residual matrices (genes x cells,
## same gene order) for MIX.SC, MIX.SE, HYB.SC, HYB.SE and optionally HYB.COMB (both hybrid alleles
## summed, allele identity ignored: the pair-level analog of the single-gene hybrid total).
## total = Rsc - Rse (parents), cis = Rhsc - Rhse (hybrid alleles), trans = total - cis.
## dpar_sc / dpar_se contrast the hybrid's own correlation (Rhyb) with each parent separately,
## mirroring the single-gene dominance contrasts at the pair level. The permutation null supplies
## no HYB.COMB, so the dpar terms are optional.
## Ploidy. Summing two alleles lowers each gene's latent variance and leaves the between-gene
## covariance unchanged, which raises the hybrid's residual correlations above a haploid genome's.
## ploidy_f (one named factor per gene, PLOIDY.F in analysis.R) rescales the shrunken Rhyb to
## R_ij * f_i * f_j, the per-genome scale of the haploid parents, before the dpar contrasts. The
## default reads attr(resid, "ploidy_f"), which RESID carries into the cluster jobs. Total, cis
## and trans never use HYB.COMB and are unaffected.
coexpr_decompose <- function(resid, ploidy_f = attr(resid, "ploidy_f")) {
  Rsc  <- shrink_cor(t(resid$MIX.SC)); Rse  <- shrink_cor(t(resid$MIX.SE))
  Rhsc <- shrink_cor(t(resid$HYB.SC)); Rhse <- shrink_cor(t(resid$HYB.SE))
  total <- Rsc - Rse; cis <- Rhsc - Rhse
  out <- list(total = total, cis = cis, trans = total - cis,
             lambda = c(Sc = attr(Rsc, "lambda"), Se = attr(Rse, "lambda"),
                        HYB.SC = attr(Rhsc, "lambda"), HYB.SE = attr(Rhse, "lambda")))
  if (!is.null(resid$HYB.COMB)) {
    Rhyb <- shrink_cor(t(resid$HYB.COMB))
    lam.h <- attr(Rhyb, "lambda")
    if (!is.null(ploidy_f)) {
      fh <- unname(ploidy_f[rownames(resid$HYB.COMB)])
      stopifnot(length(fh) == nrow(Rhyb), all(is.finite(fh)))
      Rhyb <- Rhyb * outer(fh, fh); diag(Rhyb) <- 1   # per-genome scale of a haploid parent
    }
    out$dpar_sc <- Rhyb - Rsc
    out$dpar_se <- Rhyb - Rse
    out$lambda  <- c(out$lambda, HYB.COMB = lam.h)
  }
  out
}

## Cis/trans split of one axis's eigenvalue. Because total = cis + trans elementwise,
## v'(total)v = v'(cis)v + v'(trans)v holds exactly for any vector v, so the axis's eigenvalue
## partitions into cis and trans parts on the same vector.
## v: a unit eigenvector (e.g. RANK.CHECK$vectors[, k], axis k by |eigenvalue|). PT: the coexpr_decompose() point estimate
## (COEXPR.POINT) with $cis/$trans/$total. Returns the three quadratic forms, their cis + trans sum
## (check_sum, equal to total) and frac_trans = trans / total.
coexpr_axis_cis_trans <- function(v, PT) {
  cis_part   <- as.numeric(t(v) %*% PT$cis   %*% v)
  trans_part <- as.numeric(t(v) %*% PT$trans %*% v)
  total_part <- as.numeric(t(v) %*% PT$total %*% v)
  c(cis = cis_part, trans = trans_part, total = total_part, check_sum = cis_part + trans_part, frac_trans = trans_part / total_part)
}

## Cell-resampling draws for the co-expression bootstrap. All B index sets are drawn up front from
## one seed, so draws are independent of any later RNG use and a cluster job can parallelize over them.
## SC, SE and H resample the Sc-parent, Se-parent and hybrid cells; HYB.SC, HYB.SE and HYB.COMB share
## the one H resample per draw (paired alleles of the same cells), as in make_draws().
make_coexpr_draws <- function(nSC, nSE, nH, B, seed = 1) {
  stopifnot(all(c(nSC, nSE, nH) >= 1), is.numeric(B), length(B) == 1, B >= 1)
  with_local_seed(seed, lapply(seq_len(B), function(b) list(
    SC = sample.int(nSC, nSC, replace = TRUE),
    SE = sample.int(nSE, nSE, replace = TRUE),
    H  = sample.int(nH,  nH,  replace = TRUE))))
}

## Per-gene reliability: the fraction of a gene's total NB variance (mu + mu^2/k) that is biological
## rather than Poisson sampling noise, rho = mu / (mu + k), from the MU/DISP already fitted for
## every gene. A pairwise correlation is attenuated by about sqrt(rho_i * rho_j), so rho predicts how
## noisy a pair's residual correlation is and is a natural basis for an inclusion floor.
## fits: named list of per-dataset fit data.frames (MU, DISP columns, rownames = gene). Reliability
## is the minimum across datasets, since a gene used in every dataset is only as reliable as its
## worst context. Returns a named vector over genes.
gene_reliability <- function(fits, genes) {
  stopifnot(is.list(fits), length(fits) >= 1, is.character(genes), length(genes) >= 1,
            all(vapply(fits, function(fr) all(c("MU", "DISP") %in% names(fr)), logical(1))))
  ## One column per dataset (a matrix also for a single gene). A gene absent from a fit gets NA there.
  rho_mat <- do.call(cbind, lapply(fits, function(fr) {
    mu <- fr[genes, "MU"]; k <- fr[genes, "DISP"]
    mu / (mu + k)     # DISP = Inf (Poisson limit) correctly gives rho = 0
  }))
  ## A gene with no value in any dataset is NA rather than Inf from min() of nothing
  setNames(apply(rho_mat, 1, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE)), genes)
}

## ============================================================
## Intrinsic / extrinsic noise: reliability calibration
## ============================================================
## Before choosing a reliability floor for the allele-pair correlation (a gene's within-hybrid-cell
## correlation between its Sc- and Se-allele NB Pearson residuals), this checks how the correlation's
## bootstrap SE varies with the predicted attenuation sqrt(rho_Sc * rho_Se), the same diagnostic as
## the co-expression reliability check in analysis.R Section 4. There is one correlation per gene rather than a p x p matrix, so every
## gene can be checked on its own. Because the intended use divides the correlation by
## sqrt(rho_Sc * rho_Se) (disattenuation), a second effect is examined: at low reliability the
## division inflates estimation noise on top of the attenuation. Raw and disattenuated SE are plotted
## side by side so a floor can be read off each.

## Row-wise Pearson correlation between two same-shape matrices (here genes x cells), vectorized
## across rows; rows with zero variance give NA.
row_cor <- function(A, B) {
  am <- rowMeans(A); bm <- rowMeans(B)
  Ac <- A - am; Bc <- B - bm
  den <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  ifelse(den > 0, rowSums(Ac * Bc) / den, NA_real_)
}

## Count, mean and variance of a continuous score within each level of a classification vector
## (rows with NA class or score are dropped). Used for several gene subsets and classifications
## (regulatory, dominance, burst kinetics).
class_mean_var <- function(x, class) {
  ok  <- !is.na(class) & !is.na(x)
  by  <- list(class = class[ok])
  agg <- aggregate(x[ok], by = by, FUN = length)
  data.frame(class = agg$class, n = as.numeric(agg$x), mean = aggregate(x[ok], by = by, FUN = mean)$x,
             var = aggregate(x[ok], by = by, FUN = var)$x, row.names = NULL)
}

## Splits genes into low, average and high groups by a continuous score, using the outer quantiles
## (probs) as cutoffs. Genes with a missing score or a missing ID belong to no group, so the three
## gene vectors hold only valid IDs for enrichment against a shared background.
frac_group_sets <- function(genes, score, probs = c(0.25, 0.75)) {
  scored <- !is.na(genes) & !is.na(score)
  genes <- genes[scored]; score <- score[scored]
  cuts <- quantile(score, probs)
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
  ## One class, or no residual degrees of freedom (every class holds a single value): no test to run
  if (nlevels(df$class) < 2 || nrow(df) <= nlevels(df$class))
    return(list(f = NA_real_, df1 = max(nlevels(df$class) - 1L, 0L), df2 = max(nrow(df) - nlevels(df$class), 0L),
                p = NA_real_, eta_sq = NA_real_, tukey = NULL))
  fit <- aov(x ~ class, data = df)
  ov  <- summary(fit)[[1]]
  f_stat  <- ov["class", "F value"]; p_val <- ov["class", "Pr(>F)"]
  ss_between <- ov["class", "Sum Sq"]; ss_total <- ss_between + ov["Residuals", "Sum Sq"]
  eta_sq  <- ss_between / ss_total
  tukey   <- if (is.finite(p_val) && p_val < sig) TukeyHSD(fit)$class else NULL
  list(f = f_stat, df1 = ov["class", "Df"], df2 = ov["Residuals", "Df"], p = p_val, eta_sq = eta_sq, tukey = tukey)
}

## ============================================================
## Figure 11 : co-expression cis vs trans (per gene pair); supplement: co-expression dominance
## ============================================================
## plot_coexpr_scatter(): scatter of co-expression pairs coloured by pair-level class. For type
## "cis_trans" CT is coexpr_class_table(CB): cis against trans change in correlation, coloured by the
## BH-adjusted regulatory call used for the pair lists and gene degree (Figure 11). For type
## "dominance" CT is coexpr_dom_class_table(CB): hybrid minus Sc parent against hybrid minus Se
## parent, the pairwise analogue of plot_dom_class(frame = "parent"), so the two read the same way
## side by side. Both axes share one symmetric range (unless lim is given).
plot_coexpr_scatter <- function(CT, type = c("cis_trans", "dominance"), main = NULL, lim = NULL) {
  type <- match.arg(type)
  spec <- switch(type,
    cis_trans = list(x = "cis_est", y = "trans_est", levels = REG.CLASS, cols = COLOR.LIST.1,
                     xlab = "cis: hybrid Sc-Se change in correlation", ylab = "trans: parents - hybrid (change in correlation)",
                     main = "co-expression: cis vs trans"),
    dominance = list(x = "dpar_sc_est", y = "dpar_se_est", levels = DOM.CLASS, cols = COLOR.LIST.2,
                     xlab = "hybrid - Sc parent (correlation change)", ylab = "hybrid - Se parent (correlation change)",
                     main = "co-expression: dominance"))
  x <- CT[[spec$x]]; y <- CT[[spec$y]]; cls <- CT$class
  ok <- is.finite(x) & is.finite(y) & !is.na(cls)
  if (is.null(main)) main <- spec$main
  .se_scatter(x[ok], y = y[ok], xlab = spec$xlab, ylab = spec$ylab, main = main, lim = lim, pt_col = spec$cols[match(cls[ok], spec$levels)],
              pt_cex = 0.4, legend_labels = spec$levels, legend_cols = spec$cols)
}

## Shared core for the two-seed bootstrap SE adequacy checks below.
## Given two aligned SE vectors from independent bootstrap seeds,
## reports whether the SE has converged: a scatter against the y = x
## line, the seed-to-seed Spearman correlation, and the median and IQR
## of the SE ratio. Used by both gene_seed_compare() (per-gene) and
## coexpr_seed_check() (per-pair), which differ only in how they
## extract and align se1/se2 before calling this.
seed_compare_core <- function(se1, se2, main, count_label) {
  ok  <- is.finite(se1) & is.finite(se2)
  se1 <- se1[ok]; se2 <- se2[ok]
  ratio <- se1 / se2

  plot(se1, se2, pch = 16, cex = 0.4, col = adjustcolor(COLOR.GREY[["dark"]], 0.3), xlab = "SE, seed 1", ylab = "SE, seed 2", main = main)
  abline(0, 1, col = COLOR.ACCENT, lty = 2)

  out <- list(cor = cor(se1, se2, method = "spearman"),
              ratio_median = median(ratio), ratio_iqr = IQR(ratio))
  out[[count_label]] <- length(se1)
  out
}

## Two-seed adequacy check for the per-gene bootstrap (BURST.CONTRASTS), the per-gene counterpart of
## coexpr_seed_check(). If N.BOOT has not converged, a gene's SE changes between bootstrap runs,
## seen as a low seed-to-seed correlation or a wide SE-ratio spread. bc1, bc2: two
## add_burst_contrasts()-style data frames at the same N.BOOT and gene set that differ only in
## SEED.BOOT. Genes are matched by name, so row order need not agree.
gene_seed_compare <- function(bc1, bc2, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), mode = c("total", "cis", "trans", "dom", "dpar_sc", "dpar_se", "inh_sc", "inh_se")) {
  quantity <- match.arg(quantity); mode <- match.arg(mode)
  col <- paste0(quantity, "_", mode, "_se")
  m   <- match(bc1$gene, bc2$gene)
  seed_compare_core(bc1[[col]], bc2[[col]][m], main = sprintf("%s %s: bootstrap SE, two seeds", quantity, mode), count_label = "n_genes")
}

## Rank-k eigendecomposition of a symmetric gene-by-gene divergence matrix (M ~ sum_k lambda_k v_k
## v_k^T) and how well that low-rank reconstruction predicts the observed off-diagonal entries.
## total_mat: full p x p point-estimate matrix, e.g. COEXPR.POINT$total. k: number of leading
## components. A divergence matrix has eigenvalues of both signs, so the eigenpairs are sorted once here by
## |eigenvalue|, largest first, and every consumer (the rank-k reconstruction, the candidate axes, the
## axis tests, loading1 and vectors[, j]) reads the same ordering: axis j is the j-th largest in
## magnitude. Plots reconstruction vs observed on the open device. Returns the reconstruction R^2 (squared
## correlation), all eigenvalues and eigenvectors in that order, and loading1 (the first eigenvector, named by gene).
coexpr_rank_check <- function(total_mat, k = 1) {
  stopifnot(is.numeric(k), length(k) == 1, k >= 1, k <= nrow(total_mat))
  diag(total_mat) <- 0
  eig  <- eigen(total_mat, symmetric = TRUE)
  by_magnitude <- order(abs(eig$values), decreasing = TRUE)
  eig$values   <- eig$values[by_magnitude]
  eig$vectors  <- eig$vectors[, by_magnitude, drop = FALSE]

  V     <- eig$vectors[, seq_len(k), drop = FALSE]
  recon <- V %*% diag(eig$values[seq_len(k)], k, k) %*% t(V)

  ij   <- which(upper.tri(total_mat), arr.ind = TRUE)
  obs  <- total_mat[ij]
  pred <- recon[ij]
  r2   <- cor(obs, pred)^2

 plot(pred, obs, pch = 16, cex = 0.3, col = adjustcolor(COLOR.GREY[["dark"]], 0.3), xlab = sprintf("rank-%d reconstruction", k), ylab = "observed total_est", main = sprintf("rank-%d model R^2 = %.3f", k, r2))
  abline(0, 1, col = COLOR.ACCENT, lty = 2)

  list(r2 = r2, values = eig$values, vectors = eig$vectors, loading1 = setNames(eig$vectors[, 1], rownames(total_mat)))
}

## ============================================================
## 8. REGULATORY AND DOMINANCE CLASSIFICATION  (offset level)
## ============================================================
## Per-gene calls take significance from the permutation q-values in PR and effect
## direction from the bootstrap estimates in BURST.CONTRASTS (offset level). Classifiers
## return class as a plain character vector, so clean_reg() can relabel freely; the
## call site applies ordering with factor(x, levels = REG.CLASS / DOM.CLASS).

## --- regulatory cis/trans, two-test scheme ----------------------------
## Significance comes from the cis and trans permutation p-values only (the parental
## total test is not used). Cis-only and trans-only genes get Cis and Trans; genes
## significant in both get Cis + Trans when the effects share a sign and Compensatory
## when they oppose, because the two-test scheme cannot separate compensatory from
## cis-by-trans. REG.CLASS is defined in the main script.

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
  ## A significant p-value with an NA point estimate leaves the sign undefined, so
  ## the gene matches none of the five rules and is labelled Ambiguous.
  cls[cls == ""] <- "Ambiguous"
  data.frame(class = cls,
             color = colors[match(cls, REG.CLASS)],
             stringsAsFactors = FALSE)
}

### Dominance classification (hybrid vs each parent)
## Inputs are the dpar q-values and estimates (dpar = log2 hybrid minus log2 parent). The
## parents are haploid and the hybrid is diploid. On the mean axis the library offset
## places an additive gene at the midparent, so the calls need no ploidy correction. On
## the noise axis HYB.COMB averages the intrinsic noise of two alleles, which raises DISP
## in the hybrid relative to both parents; noise-axis calls therefore use the
## ploidy-adjusted estimates (the BURST.CONTRASTS block in analysis.R) and the _q_ploidy q-values
## (bfreq, bsize, kbal, cv2; dom_class_vec selects them), while the mean axis uses _q.
## The inh contrasts carry a -1 log2 mean-axis baseline, so they stay descriptive and
## outside this classifier. DOM.CLASS is defined in the main script.

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
  ## A significant p-value with an NA bootstrap estimate leaves tsg undefined, so the
  ## gene matches none of the rules above and is labelled Ambiguous.
  cls[cls == ""] <- "Ambiguous"
  data.frame(class = cls,
             color = colors[match(cls, DOM.CLASS)],
             stringsAsFactors = FALSE)
}

## --- label cleaner for the relationship panels ------------------------
## clean_reg: maps the label "Cis x Trans" to "Compensatory" and sets "Ambiguous" to NA,
## so overlap statistics and tables see only the five REG.CLASS levels.
clean_reg <- function(v) {
  v[v == "Cis x Trans"] <- "Compensatory"
  v[v == "Ambiguous"]   <- NA
  v
}

## ---- per-gene class vectors aligned to BURST.CONTRASTS rows -------------------------
## reg_class_vec: classify_reg() for one quantity; q-values come from PR (matched by gene),
## estimates from BURST.CONTRASTS, and the returned class vector follows BURST.CONTRASTS row order.
reg_class_vec <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05) {
  quantity <- match.arg(quantity)
  i <- match(BURST.CONTRASTS$gene, PR$gene)
  classify_reg(PR[[paste0(quantity, "_cis_q")]][i], PR[[paste0(quantity, "_trans_q")]][i], BURST.CONTRASTS[[paste0(quantity, "_cis_est")]], BURST.CONTRASTS[[paste0(quantity, "_trans_est")]], sig = sig)$class
}

## dom_class_vec: classify_dom() for one quantity, aligned to BURST.CONTRASTS rows. basis picks
## the noise-axis q-values: "own" reads _q_ploidy (this gene's own ploidy shift) and "ind"
## reads _q_ind (a shift of 1). The mean axis always reads _q. A missing column stops with a
## message, so a permutation run without PLOIDY.SHIFT cannot pass through unnoticed.
dom_class_vec <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05, basis = c("own", "ind")) {
  quantity <- match.arg(quantity); basis <- match.arg(basis)
  i   <- match(BURST.CONTRASTS$gene, PR$gene)
  sfx <- if (quantity == "mean") "_q" else if (basis == "own") "_q_ploidy" else "_q_ind"
  cs <- paste0(quantity, "_dpar_sc", sfx); ce <- paste0(quantity, "_dpar_se", sfx)
  if (!all(c(cs, ce) %in% names(PR))) stop(sprintf("PR lacks %s. Rerun gene_perm.R with PLOIDY.SHIFT in the inputs.", cs))
  classify_dom(PR[[cs]][i], PR[[ce]][i], BURST.CONTRASTS[[paste0(quantity, "_dpar_sc_est")]], BURST.CONTRASTS[[paste0(quantity, "_dpar_se_est")]], sig = sig)$class
}

## class_stats: applies fun(values, class vector) to every class vector in classes (a list by kind, e.g.
## list(REG = REG.VEC, DOM = DOM.VEC), each a list by quantity), restricted to the genes in idx (positions in the
## class vectors). The result is named "<KIND>.<QUANTITY>", e.g. REG.MEAN, DOM.KBAL.
class_stats <- function(fun, values, classes, idx) {
  out <- list()
  for (kind in names(classes))
    for (q in names(classes[[kind]]))
      out[[paste0(kind, ".", toupper(q))]] <- fun(values, classes[[kind]][[q]][idx])
  out
}

## .overlap_lor(tab): log2((observed + 0.5) / (expected + 0.5)) for a class-by-class contingency
## table, expected from the marginals under independence. The 0.5 pseudocount keeps empty cells finite.
## One definition for the heatmap colours and for shared_overlap_rng().
.overlap_lor <- function(tab) {
  exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  log2((tab + 0.5) / (exp + 0.5))
}

## Class-overlap heatmaps (Figures 2 and 6 and the Section 3.5 / 4.12 diagnostics): a categorical contingency table between two
## classifications, tested cell by cell against the hypergeometric null and colored by
## log2((observed + 0.5) / (expected + 0.5)). The star in each cell marks BH-corrected
## significance at fdr; the color always encodes the fold enrichment itself.
## levels_a/levels_b: if NULL, table levels are the categories actually observed (the diagnostic
## heatmaps, where the two classifications do not always share a level set). If supplied, every level
## appears even when empty, which fixes the axis order across panels (Figures 2 and 6, where each
## panel in a row shares axes). cex_cell, fmt and star_frac default to values suited to each case
## when left NULL. Cell text shows obs/exp as fold enrichment (2^lor) with a trailing "*" where the
## BH-adjusted two-sided hypergeometric p < fdr.
## rng, when supplied, fixes the color scale's half-range in log2(obs/exp) units. Passing the
## same rng to related panels (e.g. the three regulatory overlap panels of Figure 2) makes a
## given color mean the same fold enrichment in all of them.
class_overlap_heatmap <- function(class_a, class_b, levels_a = NULL, levels_b = NULL, brk = length(cols), fdr = 0.01, cols = COLOR.LIST.3, cex_axis = 0.75, cex_cell = NULL, fmt = NULL, star_frac = NULL, xlab = NULL, ylab = NULL, rng = NULL) {
  stopifnot(is.numeric(fdr), length(fdr) == 1, fdr > 0, fdr < 1, is.numeric(brk), length(brk) == 1, brk >= 1,
            is.null(rng) || (is.numeric(rng) && length(rng) == 1 && rng > 0))
  explicit_levels <- !is.null(levels_a) || !is.null(levels_b)
  if (is.null(cex_cell))  cex_cell  <- if (explicit_levels) 0.7   else 0.65
  if (is.null(fmt))       fmt       <- if (explicit_levels) "%.1f" else "%.2g"
  if (is.null(star_frac)) star_frac <- if (explicit_levels) 0.6   else 0.55

  if (brk %% 2 == 0) brk <- brk + 1   # odd brk puts a true white sample at the center
  keep <- !is.na(class_a) & !is.na(class_b)
  a <- if (is.null(levels_a)) factor(class_a[keep]) else factor(class_a[keep], levels = levels_a)
  b <- if (is.null(levels_b)) factor(class_b[keep]) else factor(class_b[keep], levels = levels_b)
  tab <- table(a, b); n <- sum(tab)
  if (n == 0 || nrow(tab) < 2 || ncol(tab) < 2)
    stop(sprintf("class_overlap_heatmap(): needs genes classified in both vectors and at least 2 levels on each axis (got %d genes, %d x %d levels)", n, nrow(tab), ncol(tab)))
  lor <- .overlap_lor(tab)
  pv  <- matrix(NA_real_, nrow(tab), ncol(tab))
  for (i in seq_len(nrow(tab))) for (j in seq_len(ncol(tab))) {
    q <- tab[i, j]; m <- rowSums(tab)[i]; k <- colSums(tab)[j]
    po <- phyper(q - 1, m, n - m, k, lower.tail = FALSE)
    pu <- phyper(q,     m, n - m, k, lower.tail = TRUE)
    pv[i, j] <- 2 * min(po, pu, 0.5)
  }
  padj <- matrix(p.adjust(pv, "BH"), nrow(tab))
  ## Color encodes log2(obs/exp) for every cell, not just significant ones; the * in the
  ## cell text is the only significance indicator. Without a supplied rng the scale spans
  ## this panel's own largest |log2(obs/exp)|.
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
  ## --- class association heatmap ----------------------------------------
  ## The heatmap shows the association between two class vectors as a grid of log2
  ## observed-over-expected counts. Every cell is colored by its fold enrichment; a star marks
  ## cells whose two-sided hypergeometric test is significant after BH. brk sets the number of
  ## color bins.
  ## .heatmap_legend draws a vertical color strip just outside the right edge of the plot
  ## region, mapping the heatmap colors to fold enrichment. It runs right after image(), while
  ## par("usr") still describes the heatmap's own 0-1 coordinate system.
  local({
    usr <- par("usr")
    lx <- usr[2] + diff(usr[1:2]) * 0.12
    rx <- usr[2] + diff(usr[1:2]) * 0.24
    ys <- seq(usr[3], usr[4], length.out = brk + 1)
    rect(lx, ys[-length(ys)], rx, ys[-1], col = pal, border = NA, xpd = NA)
    rect(lx, usr[3], rx, usr[4], border = COLOR.GREY[["dark"]], xpd = NA)
    text(rx, usr[4],         sprintf("%.2gx", 2^rng),  pos = 4, cex = 0.55, xpd = NA)
    text(rx, mean(usr[3:4]), "1x",                      pos = 4, cex = 0.55, xpd = NA)
    text(rx, usr[3],         sprintf("%.2gx", 2^-rng),  pos = 4, cex = 0.55, xpd = NA)
    text(rx + diff(usr[1:2])*0.16, mean(usr[3:4]), "fold enrichment (obs/exp)", srt = 270, cex = 0.55, xpd = NA)
  })
  invisible(list(table = tab, log2_obs_exp = lor, padj = padj))
}

## ============================================================
## Gene-identity overlap between two categorical classifications
## ============================================================
## class_identity_overlap: Cohen's kappa plus per-class Jaccard overlap between two
## classifications of the same genes. Per-class Jaccard is computed on gene membership
## (|A_k intersect B_k| / |A_k union B_k|), so a high overall kappa driven by one large class
## (typically Conserved) is not read as overlap of the classes of biological interest
## (Cis, Trans, ...); each level gets its own number. The null reshuffles which gene carries
## which class_b label. A label permutation holds both class-size distributions fixed, so the
## chance concordance pe is identical in every replicate and is computed once. Integer
## coding plus tabulate() keeps the nperm-fold loop cheap. Concordance po is the diagonal fraction of the
## class table, chance concordance pe = sum(rowSums(tab) * colSums(tab)) / n^2 comes from the two
## marginal class-size distributions, and kappa = (po - pe) / (1 - pe) is 0 at chance agreement and 1
## for identical labelling.
class_identity_overlap <- function(class_a, class_b, levels = REG.CLASS, nperm = 2000, seed = 1) {
  stopifnot(length(class_a) == length(class_b), is.numeric(nperm), length(nperm) == 1, nperm >= 1)
  keep <- !is.na(class_a) & !is.na(class_b)
  a <- factor(class_a[keep], levels = levels)
  b <- factor(class_b[keep], levels = levels)
  ## factor() maps any value outside `levels` (e.g. "Ambiguous", a data-quality flag rather
  ## than a class) to NA after the is.na() filter above; filtering again here leaves only
  ## genes with a real level in both classifications.
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

  ai <- as.integer(a); bi <- as.integer(b)
  kappa_null <- numeric(nperm)
  with_local_seed(seed, for (i in seq_len(nperm)) {
    bp   <- bi[sample.int(n, n)]
    tabp <- matrix(tabulate((bp - 1L) * K + ai, K * K), K, K)
    po_p <- sum(diag(tabp)) / n
    kappa_null[i] <- (po_p - pe) / (1 - pe)
  })
  ## One-sided test: does observed agreement exceed what reshuffled labels produce? A small
  ## p means the two classifications share gene identity beyond what their class sizes
  ## alone predict. The add-one correction keeps p > 0.
  p_kappa <- (1 + sum(kappa_null >= kappa_obs)) / (1 + nperm)

  list(table = tab, n_genes = n,
       concordance = po, expected_concordance = pe,
       kappa = kappa_obs, kappa_null = kappa_null, kappa_p = p_kappa,
       jaccard = jacc)
}

## summarize_class_overlap: one printable row per classification-pair comparison (n,
## concordance, chance concordance, kappa, permutation p, and one Jaccard column per class
## level). Bind rows from several class_identity_overlap() results into one table.
summarize_class_overlap <- function(ov, label = NULL) {
  data.frame(comparison = label, n_genes = ov$n_genes,
            concordance = ov$concordance, expected = ov$expected_concordance,
            kappa = ov$kappa, kappa_p = ov$kappa_p,
            as.list(ov$jaccard), check.names = FALSE, stringsAsFactors = FALSE)
}

## se_alpha_col: maps a vector of SEs to per-point bar colors. Precise (small-SE) estimates get
## a more opaque bar and imprecise (large-SE) ones fade toward the background. lo/hi are the
## SE quantiles anchoring the high- and low-opacity ends of the ramp.
se_alpha_col <- function(se, base_col = COLOR.GREY[["dark"]], lo = 0.05, hi = 0.9, alpha_range = c(0.08, 0.4)) {
  stopifnot(is.numeric(lo), is.numeric(hi), lo >= 0, hi <= 1, lo < hi,
            is.numeric(alpha_range), length(alpha_range) == 2, all(alpha_range >= 0 & alpha_range <= 1))
  rng <- if (any(is.finite(se))) quantile(se[is.finite(se)], c(lo, hi)) else c(0, 0)
  width <- rng[2] - rng[1]
  ## SEs with no spread (all equal) are all equally precise; a missing SE gets the faintest bar
  w   <- if (width > 0) 1 - pmin(pmax((se - rng[1]) / width, 0), 1) else rep(1, length(se))   # w = 1 for precise, 0 for noisy
  w[!is.finite(w)] <- 0
  a   <- alpha_range[1] + diff(alpha_range) * w
  # col2rgb()/rgb() are vectorized, so building the RGBA color directly gives every point
  # its own opacity (adjustcolor() accepts a single alpha value).
  rgb_base <- grDevices::col2rgb(base_col) / 255
  grDevices::rgb(rgb_base[1], rgb_base[2], rgb_base[3], alpha = a)
}

## ============================================================
## Figure 1 : cis vs trans, coloured by regulatory class
## ============================================================
## plot_cis_trans_class: cis (x) vs trans (y) divergence for one quantity, points colored by
## regulatory class. All SE bars are drawn before the points so bars sit behind them; bar
## opacity follows se_alpha_col() so imprecise bars recede. Points are solid, matching the
## Figure 5 dominance panels.
plot_cis_trans_class <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), sig = 0.05, main = NULL, lim = NULL, bar_col = COLOR.GREY[["dark"]], colors = setNames(COLOR.LIST.1[seq_along(REG.CLASS)], REG.CLASS)) {
  quantity <- match.arg(quantity)
  lab <- QUANTITY.LABEL[[quantity]]
  cx <- BURST.CONTRASTS[[paste0(quantity,"_cis_est")]];  sx <- BURST.CONTRASTS[[paste0(quantity,"_cis_se")]]
  cy <- BURST.CONTRASTS[[paste0(quantity,"_trans_est")]]; sy <- BURST.CONTRASTS[[paste0(quantity,"_trans_se")]]
  cls <- reg_class_vec(BURST.CONTRASTS, PR, quantity, sig)
  ok  <- is.finite(cx)&is.finite(cy)&is.finite(sx)&is.finite(sy)&cls %in% REG.CLASS
  if (is.null(main)) main <- paste0(lab,": cis vs trans")
  .se_scatter(cx[ok], sx[ok], cy[ok], sy[ok], xlab = "log2(Sc/Se) in hybrid", ylab = "parents - hybrid (log2)", main = main, lim = lim,
              bar_col = list(x = se_alpha_col(sx[ok], bar_col), y = se_alpha_col(sy[ok], bar_col)), pt_col = colors[match(cls[ok], REG.CLASS)],
              legend_labels = REG.CLASS, legend_cols = colors)
  invisible(cls[ok])
}

## ============================================================
## Figure 3 (supplement) : mean vs noise per gene, per-class slopes
## class_slope_lines: one least-squares line of y on x per class level (levels with more than two genes), drawn
## across the whole panel in that level's colour.
class_slope_lines <- function(x, y, cls, levels, cols, lty) {
  for (k in levels) {
    sel <- !is.na(cls) & cls == k
    if (sum(sel) > 2) abline(lm(y[sel] ~ x[sel]), col = cols[match(k, levels)], lwd = 2, lty = lty)
  }
}

## ============================================================
## plot_mean_bfreq_class: mean divergence (x) vs a burst-parameter divergence (y) at one mode,
## with per-class regression lines. Points are neutral grey; regulatory slopes are solid and
## dominance slopes dashed. Lines span the full panel (abline).
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
  plot(x, y, pch=16, cex=0.4, col=COLOR.GREY[["light"]], xlim=range(x)+c(-px, px), ylim=range(y)+c(-py, py), xlab="mean divergence (log2)", ylab=paste0(y_lab, " divergence (log2)"), main=main)
  abline(h=0,v=0,col=COLOR.GREY[["mid"]])
  class_slope_lines(x, y, rc, REG.CLASS, reg_colors, lty = 1)
  class_slope_lines(x, y, dc, DOM.CLASS, dom_colors, lty = 2)
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
  p0   <- PR[[paste0(quantity, "_total_q")]][i]
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
  p_c   <- PR[[paste0(quantity, "_", component, "_q")]][i]
  est   <- BURST.CONTRASTS[[paste0(quantity, "_", component, "_est")]]
  sigc  <- !is.na(p_c) & p_c < sig
  genes <- BURST.CONTRASTS$gene
  tag   <- if (component == "cis") "AnyCis" else "AnyTrans"

  by_dir <- list(Sc = intersect(genes[sigc & est > 0], universe), Se = intersect(genes[sigc & est < 0], universe))
  setNames(by_dir, paste0(tag, "_", names(by_dir)))
}

## run_enrichment: GO (BP/CC/MF, simplified) and KEGG enrichment for one gene set against a
## fixed universe. Genes outside the universe are dropped; sets with fewer than 2 genes after
## intersection return all-NULL without calling enrichGO/enrichKEGG. kegg_data (from
## kegg_local()) runs KEGG offline through enricher() with enrichKEGG()'s defaults (BH,
## gene-set size 10 to 500). Every term is judged on its Benjamini-Hochberg adjusted p-value
## (p.adjust) alone: fdr is the BH cutoff (passed as pvalueCutoff, which clusterProfiler applies to
## p.adjust) and the Storey q-value cutoff is left open (qvalueCutoff = 1), so the qvalue column
## never filters terms and is not used. Results for small sets (roughly under 10-15 genes) are
## exploratory; the direction-split class sets (e.g. Compensatory_Se) are the most likely to
## be that small, and n_genes in the summary table flags it.
run_enrichment <- function(genes, universe, orgdb = org.Sc.sgd.db, keytype = "ORF", kegg_org = "sce", fdr = 0.2, kegg_data = NULL) {
  stopifnot(is.numeric(fdr), length(fdr) == 1, fdr > 0, fdr <= 1)
  genes    <- unique(genes[!is.na(genes)])
  universe <- unique(universe[!is.na(universe)])
  genes    <- intersect(genes, universe)
  if (length(genes) < 2) return(list(BP = NULL, CC = NULL, MF = NULL, KEGG = NULL))
  go_one <- function(ont) {
    tryCatch(simplify(enrichGO(gene = genes, universe = universe, OrgDb = orgdb, keyType = keytype, ont = ont, pvalueCutoff = fdr, qvalueCutoff = 1)), error = function(e) NULL)
  }
  kegg <- tryCatch(if (is.null(kegg_data)) enrichKEGG(gene = genes, universe = universe, organism = kegg_org, keyType = "kegg", pvalueCutoff = fdr, qvalueCutoff = 1)
                   else enricher(gene = genes, universe = universe, TERM2GENE = kegg_data$KEGGPATHID2EXTID,
                                 TERM2NAME = kegg_data$KEGGPATHID2NAME, pvalueCutoff = fdr, qvalueCutoff = 1),
                   error = function(e) NULL)
  list(BP = go_one("BP"), CC = go_one("CC"), MF = go_one("MF"), KEGG = kegg)
}

## axis_pole_enrichment: GO (BP, MF or CC, simplified) or KEGG enrichment of one pole of a co-expression axis
## (genes) against the co-expressed universe, at a Benjamini-Hochberg cutoff fdr on p.adjust (the Storey q-value
## cutoff is left open). A pole with fewer than min_genes genes returns NULL without testing.
axis_pole_enrichment <- function(genes, universe, ont = c("BP", "MF", "CC", "KEGG"), min_genes = 0, fdr = 0.05) {
  ont <- match.arg(ont)
  stopifnot(is.numeric(fdr), length(fdr) == 1, fdr > 0, fdr <= 1)
  if (length(genes) < min_genes) return(NULL)
  if (ont == "KEGG") enrichKEGG(gene = genes, universe = universe, organism = "sce", pvalueCutoff = fdr, qvalueCutoff = 1)
  else simplify(enrichGO(gene = genes, universe = universe, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = ont, pvalueCutoff = fdr, qvalueCutoff = 1))
}

## n_sig_terms: number of terms at BH q (p.adjust) < q in one enrichResult, 0 for a
## NULL or empty result (a set too small to test, or no hits)
n_sig_terms <- function(e, q = 0.2) {
  if (is.null(e) || is.null(e@result) || nrow(e@result) == 0) return(0L)
  sum(e@result$p.adjust < q, na.rm = TRUE)
}

## print_enrich_brief: console-friendly view of one enrichResult, since
## the default print() dumps every column (including the full comma-
## separated gene list per term) and R's console truncates that to
## fit width, hiding exactly the columns worth reading. Shows only
## Description, p.adjust, and Count, for the n_top most significant
## terms at BH q (p.adjust) < q. Prints one line and returns invisibly for a
## NULL result or one with nothing significant, rather than an empty
## table with no explanation.
print_enrich_brief <- function(e, q = 0.2, n_top = 10) {
  if (is.null(e) || is.null(e@result) || nrow(e@result) == 0) { cat("  (no terms tested)\n"); return(invisible(NULL)) }
  tab <- e@result[e@result$p.adjust < q, c("Description", "p.adjust", "Count"), drop = FALSE]
  if (nrow(tab) == 0) { cat(sprintf("  (0 terms at BH q < %g)\n", q)); return(invisible(NULL)) }
  tab <- tab[order(tab$p.adjust), ][seq_len(min(n_top, nrow(tab))), ]
  tab$p.adjust <- signif(tab$p.adjust, 3)
  print(tab, row.names = FALSE)
}

## plot_cluster_marker_enrichment: writes the up/down marker enrichment of every cluster (from
## cluster_stability.R's dataset_marker_enrichment()) to one pdf; each page set is labeled with its
## "up" cluster so pages are identifiable out of context.
plot_cluster_marker_enrichment <- function(res, label, pdf_path, width = 7, height = 5) {
  if (is.null(res)) return(invisible(NULL))
  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  for (cc in res$cluster_ids) {
    ## barplot_enrich_pair: barplots each ontology (BP, MF, CC, KEGG) from a pair of
    ## run_enrichment()-style lists (UP and DOWN gene sets), skipping any ontology that is NULL
    ## or empty. Each plot is built and print()ed inside one tryCatch, because ggplot2 defers
    ## scale training and stat transforms (where enrichplot's barplot() can fail on a
    ## very-few-row enrichResult) until print(). A failed plot is reported to the console and
    ## skipped, so the enclosing pdf() block stays open until the function closes it on exit and later clusters are drawn.
    ## label names the comparison (e.g. "Sc major vs minor cluster"), since "up"/"down" means
    ## ident.1 vs ident.2 of the FindMarkers() call and that pairing differs across call sites.
    local({
      enrich_up <- res$up_enrich[[cc]]
      enrich_down <- res$down_enrich[[cc]]
      show <- 10
      label <- sprintf("%s cluster %s vs rest (up = higher in %s)", label, cc, cc)
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
    })
  }
}

## score_cell_cycle_by_cluster: cell-cycle phase scoring on a dataset's own validated
## clustering (Idents already set, e.g. YSC$MIX.SE, YSC$HYB.SC, YSC$HYB.SE). Scores the
## curated yeast regulons (S.GENES/G2M.GENES/MG1.GENES) with CellCycleScoring()/
## AddModuleScore(), draws a violin plot and a phase-composition-by-cluster barplot, and
## prints the phase table as a quantitative check of a GO-based reading of each cluster.
## Returns obj with S.Score/G2M.Score/MG1.Score1/Phase added, for the caller to reassign
## (e.g. YSC$MIX.SE <- score_cell_cycle_by_cluster(YSC$MIX.SE, ...)).
score_cell_cycle_by_cluster <- function(obj, label, pdf_path, s_genes = S.GENES, g2m_genes = G2M.GENES, mg1_genes = MG1.GENES, width = 9, height = 8) {
  obj <- suppressWarnings(suppressMessages(CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes)))
  obj <- suppressWarnings(suppressMessages(AddModuleScore(obj, features = list(mg1_genes), name = "MG1.Score")))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 1))
  print(suppressWarnings(VlnPlot(obj, features = c("S.Score", "G2M.Score", "MG1.Score1"), ncol = 1, pt.size = 0, cols = colorRampPalette(COLOR.CLUSTER)(nlevels(Idents(obj))))))

  phase_by_cluster <- prop.table(table(Idents(obj), obj$Phase), margin = 1)
  barplot(t(phase_by_cluster), col = COLOR.PHASE[colnames(phase_by_cluster)],
          legend.text = TRUE, args.legend = list(x = "topright", bty = "n"),
          las = 2, ylab = "fraction of cells", main = sprintf("%s: cell-cycle phase composition by cluster", label))

  cat(sprintf("%s: cell-cycle phase composition by cluster\n", label)); print(round(phase_by_cluster, 3))
  obj
}

## go_gene_set: all genes annotated to a GO biological-process term, including descendant
## terms (GOALL), from org.Sc.sgd.db, the annotation source used for every enrichment test in
## this pipeline. Gene sets drawn from the annotation are independent of any cluster's own
## markers, like the curated cell-cycle regulons (S.GENES/G2M.GENES/MG1.GENES).
go_gene_set <- function(go_id, orgdb = org.Sc.sgd.db) {
  unique(AnnotationDbi::select(orgdb, keys = go_id, keytype = "GOALL", columns = "ORF")$ORF)
}

## score_modules_by_cluster: continuous module-score validation, the counterpart of
## score_cell_cycle_by_cluster() for gene sets without discrete phases (e.g. glycolysis,
## oxidative phosphorylation, ribosome biogenesis). Scores obj with AddModuleScore() for each
## named set in gene_sets, draws a violin plot per module and a per-cluster mean +/- SE
## barplot, so a qualitative GO-term reading of a cluster's markers can be checked against
## independently sourced scores. gene_sets is a named list (e.g. list(Glycolysis =
## go_gene_set("GO:0006096"), ...)); names become the AddModuleScore() name prefix and the
## plot labels.
score_modules_by_cluster <- function(obj, gene_sets, label, pdf_path, width = 9, height = 8) {
  score_names <- paste0(names(gene_sets), "1")
  for (i in seq_along(gene_sets)) obj <- suppressWarnings(suppressMessages(AddModuleScore(obj, features = list(gene_sets[[i]]), name = names(gene_sets)[i])))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 1))
  print(suppressWarnings(VlnPlot(obj, features = score_names, ncol = 1, pt.size = 0, cols = colorRampPalette(COLOR.CLUSTER)(nlevels(Idents(obj))))))

  cl    <- Idents(obj)
  sc    <- as.data.frame(obj[[score_names]])
  means <- sapply(sc, tapply, cl, mean)
  ses   <- sapply(sc, tapply, cl, sd) / sqrt(sapply(sc, tapply, cl, length))

  par(mfrow = c(1, length(score_names)), mar = c(5, 4.5, 3, 1))
  for (i in seq_along(score_names)) {
    b <- barplot(means[, i], ylim = range(c(means[, i] - ses[, i], means[, i] + ses[, i])),
                 main = sprintf("%s: %s", label, names(gene_sets)[i]), ylab = "mean module score", las = 2)
    arrows(b, means[, i] - ses[, i], b, means[, i] + ses[, i], angle = 90, code = 3, length = 0.05)
  }

  cat(sprintf("%s: mean module score by cluster\n", label)); print(round(means, 3))
  obj
}

## ============================================================
## 8g. Covariate-noise diagnostic: does per-cell NB noise track
## cell-cycle position or metabolic state?
## ============================================================
## cell_cycle_continuum(obj): one continuous cell-cycle position per cell.
## S.Score and G2M.Score from CellCycleScoring() trace a roughly circular
## path through the cycle, so their first principal component is the closest
## one-dimensional summary. A rank correlation against this axis is
## unchanged by per-object recentering or rescaling. Returns a named numeric
## vector, one value per cell
cell_cycle_continuum <- function(obj) {
  missing <- setdiff(c("S.Score", "G2M.Score"), colnames(obj[[]]))
  if (length(missing) > 0)
    stop(sprintf("cell_cycle_continuum(): %s not found; run score_cell_cycle_by_cluster() (or CellCycleScoring()) on this object first and reassign its result", paste(missing, collapse = ", ")))
  sc  <- obj[[c("S.Score", "G2M.Score")]]
  pc1 <- prcomp(sc, scale. = TRUE)$x[, 1]
  setNames(pc1, rownames(sc))
}

## cell_cycle_continuum_shared(obj1, obj2): cell-cycle axis for comparing two
## objects (used by species_composition_bound() in 8h). prcomp() centers its
## scores to mean zero, so each object's separately fitted axis has mean zero
## by construction. One PCA on the pooled S.Score/G2M.Score of both objects,
## with each object's cells projected onto the shared PC1, puts both on one
## scale so a difference in cell-cycle composition shows up as a difference
## in means. Use cell_cycle_continuum() for the within-object correlation
## test in 8g. Returns list(x1, x2), named by cell
cell_cycle_continuum_shared <- function(obj1, obj2) {
  sc1    <- obj1[[c("S.Score", "G2M.Score")]]
  sc2    <- obj2[[c("S.Score", "G2M.Score")]]
  scores <- prcomp(rbind(sc1, sc2), scale. = TRUE)$x[, 1]
  list(x1 = setNames(scores[seq_len(nrow(sc1))], rownames(sc1)),
       x2 = setNames(scores[(nrow(sc1) + 1):length(scores)], rownames(sc2)))
}

## metabolic_state_cluster(obj): discrete metabolic state per cell from the
## continuous module scores of score_modules_by_cluster() (e.g.
## Glycolysis1/OXPHOS1/RiBi1). A cell can sit high on several pathways at
## once, so k-means on all scores jointly captures states defined by a
## combination that per-pathway thresholds would miss. k-means is fit once for every k in
## k_range (nstart = 10 random starts each) and the fit with the highest mean silhouette
## width is returned as it was scored (ties go to the smaller k), so the labels are the
## clustering the silhouette judged. The starts come from their own seeded stream (seed),
## and the caller's random number state is restored on exit, so the states are the same on
## every run and the function does not move any other draw.
## Returns a factor of state labels, one per cell
metabolic_state_cluster <- function(obj, score_names = c("Glycolysis1", "OXPHOS1", "RiBi1"), k_range = 2:4, seed = 1) {
  missing <- setdiff(score_names, colnames(obj[[]]))
  if (length(missing) > 0)
    stop(sprintf("metabolic_state_cluster(): %s not found; run score_modules_by_cluster() on this object first and reassign its result", paste(missing, collapse = ", ")))
  scores <- scale(as.matrix(obj[[score_names]]))
  stopifnot(is.numeric(seed), length(seed) == 1, is.finite(seed),
            is.numeric(k_range), length(k_range) >= 1, all(k_range >= 2), all(k_range < nrow(scores)))
  d    <- dist(scores)
  fits <- with_local_seed(seed, lapply(k_range, function(k) kmeans(scores, centers = k, nstart = 10)))
  sil  <- vapply(fits, function(fit) mean(silhouette(fit$cluster, d)[, 3]), numeric(1))
  km   <- fits[[which.max(sil)]]
  setNames(factor(km$cluster), rownames(scores))
}

## covariate_noise_diagnostic(): tests whether per-cell NB Pearson residuals
## (nb_residuals(), from the exposure-offset fit) depend on cell-cycle
## position or metabolic state. The squared correlation of a gene's residual
## with a covariate is the share of its residual variance that covariate
## explains, so a small effect size means unmodelled cell state contributes
## little to the measured noise. Cell-cycle axis: per-gene Spearman rho and
## BH q. Metabolic state: per-gene Kruskal-Wallis BH q and eta-squared
## (fraction of residual variance explained by state). The fit is not
## refitted with these covariates. Prints the fraction of genes at q < 0.05
## and the median effect size, and returns one row per gene
covariate_noise_diagnostic <- function(mat, expo, fit, cc_axis, metab_state, label) {
  cells <- Reduce(intersect, list(colnames(mat), names(cc_axis), names(metab_state)))
  R   <- nb_residuals(mat[, cells], expo[cells], fit)
  cc  <- cc_axis[cells]
  met <- droplevels(metab_state[cells])

  ## One pass over genes: Spearman correlation with the cell-cycle axis, and eta^2 (share of residual
  ## variance explained by metabolic state) with its Kruskal-Wallis p-value.
  cc_rho <- cc_p <- met_eta <- met_p <- setNames(rep(NA_real_, nrow(R)), rownames(R))
  for (i in seq_len(nrow(R))) {
    r <- R[i, ]
    cc_rho[i] <- suppressWarnings(cor(r, cc, method = "spearman", use = "complete.obs"))
    cc_p[i]   <- suppressWarnings(cor.test(r, cc, method = "spearman")$p.value)
    ss <- summary(aov(r ~ met))[[1]][["Sum Sq"]]
    met_eta[i] <- ss[1] / sum(ss)
    met_p[i]   <- suppressWarnings(kruskal.test(r, met)$p.value)
  }
  cc_q  <- p.adjust(cc_p, method = "BH")
  met_q <- p.adjust(met_p, method = "BH")

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
  hist(diag_df$cc_rho, breaks = 40, col = COLOR.GREY[["mid"]], border = NA,
       xlab = "Spearman rho, residual vs cell-cycle axis", main = sprintf("%s: cell cycle", label))
  hist(diag_df$met_eta, breaks = 40, col = COLOR.GREY[["mid"]], border = NA,
       xlab = "eta-squared, residual vs metabolic state", main = sprintf("%s: metabolic state", label))
}

## ============================================================
## 8h. Species composition comparison and confound bound
## ============================================================
## species_composition_report(): runs species_composition_bound() on the
## cell-cycle axis and on each metabolic module score for one species/allele
## pair. The effect size for each axis is the larger of the two views' median
## per-gene 8g effect sizes, so the bound reflects the stronger covariate
## dependence of the two. The result's bound_sd is an expected shift in residual-SD units (|Cohen's d| x
## effect size). For the metabolic scores sqrt(median eta-squared)
## is the effect size: eta-squared is the categorical analog of a squared
## correlation, so its square root is on the same standardized scale as
## Cohen's d and the cell-cycle rho
species_composition_report <- function(cc1, cc2, met1, met2, diag1, diag2, label, paired = FALSE) {
  ## A small covariate-noise effect (8g) rules out a cell-state confound on the
  ## cross-species comparison only if the two species also differ little in
  ## composition along that axis. Each axis (the cell-cycle axis, or one
  ## metabolic module score) is compared between two species/allele views by
  ## Cohen's d (pooled SD, so d is the composition shift in SD units of the axis) and a rank test.
  ## Multiplying |d| by the matching 8g effect size (a correlation) gives the expected shift, in
  ## residual-SD units, that a composition difference of that size can produce (a linear, bivariate-
  ## normal approximation: E[Y | X shifted by d SD] = rho * d SD). That shift (bound_sd) is not a share
  ## of residual variance; the variance share would be roughly its square. The bound sizes the
  ## confound; it is not a refit or a test of any specific contrast. With thousands of cells any
  ## test rejects, so the p-value is descriptive and the bound is the informative output.
  ## paired = FALSE compares different cells (the parents) with the rank-sum test. paired = TRUE
  ## compares two views of the same cells (the Sc and Se alleles of the hybrid cells): the values
  ## are matched by cell name and the signed-rank test is used, since the rank-sum test would
  ## treat the 2n paired values as independent samples. Cohen's d keeps the pooled SD in both cases
  ## (a paired d_z would inflate the bound with the allele correlation).
  species_composition_bound <- function(x1, x2, effect_size, label, axis_label, paired) {
    if (paired) {
      stopifnot(!is.null(names(x1)), !is.null(names(x2)), setequal(names(x1), names(x2)))
      x2 <- x2[names(x1)]
    }
    d     <- (mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) / sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
    ok    <- is.finite(x1) & is.finite(x2)
    test  <- if (paired) "signed-rank" else "rank-sum"
    wt    <- if (paired) wilcox.test(x1[ok], x2[ok], paired = TRUE, exact = FALSE) else wilcox.test(x1, x2, exact = FALSE)
    bound_sd <- abs(d) * effect_size
    cat(sprintf("%s, %s: Cohen's d = %.3f (%s p = %.2g, descriptive), noise-association effect size = %.3f, bound on expected noise shift = %.3f residual SD\n",
                label, axis_label, d, test, wt$p.value, effect_size, bound_sd))
    data.frame(label = label, axis = axis_label, cohens_d = d, test = test, test_p = wt$p.value, effect_size = effect_size, bound_sd = bound_sd)
  }

  cc_effect  <- max(median(abs(diag1$cc_rho), na.rm = TRUE), median(abs(diag2$cc_rho), na.rm = TRUE))
  met_effect <- sqrt(max(median(diag1$met_eta, na.rm = TRUE), median(diag2$met_eta, na.rm = TRUE)))
  rows <- list(species_composition_bound(cc1, cc2, cc_effect, label, "cell-cycle axis", paired))
  for (sn in colnames(met1)) rows[[length(rows) + 1]] <- species_composition_bound(met1[, sn], met2[, sn], met_effect, label, sn, paired)
  do.call(rbind, rows)
}

## ============================================================
## 9. PUBLICATION FIGURES
## ============================================================

## ---- figure colour palettes ---------------------------------------------
## COLOR.LIST.1 (regulatory classes) and COLOR.LIST.2 (dominance classes)
## are defined and named in the main script, matching REG.CLASS and
## DOM.CLASS

## ============================================================
## Figures 2 and 6 : mean-class x size-class enrichment heatmaps (regulatory classes; dominance classes)
## ============================================================
## shared_overlap_rng(pairs, levels_common): common log2(obs/exp) half-range
## across a list of (class_a, class_b) pairs built on the same levels, using
## the same 0.5 pseudocount as class_overlap_heatmap(). Passing it as rng to
## a set of related panels gives them one colour scale.
shared_overlap_rng <- function(pairs, levels_common) {
  lors <- lapply(pairs, function(p) {
    .overlap_lor(table(factor(p[[1]], levels = levels_common), factor(p[[2]], levels = levels_common)))
  })
  max(abs(unlist(lors)), 1e-6)
}

## ============================================================
## Dominance scatter: parent frame (primary) or A/D rotation (supplement)
## ============================================================
## plot_dom_class(): dominance scatter of the hybrid against each parent (Figure 5, frame = "parent").
## frame="parent" (x = hybrid - Sc, y = hybrid - Se) plots exactly the two
## contrasts classify_dom() tests, with no midparent construction. Additive
## needs only that the two contrasts have opposite sign (the hybrid lies
## between the parents), which holds on any monotonic rescaling. A rotation
## keeps straight lines through the origin, so the six class regions are
## equally clean in either frame.
##
## frame="AD" plots additive (Sc - Se)/2 against dominance (hybrid -
## midparent). For noise this is an exact rotation of the tested contrasts,
## because the noise midparent is the geometric (log2) average, (dpar_sc +
## dpar_se)/2. For mean it is a different construction: the mean midparent
## is the arithmetic average of the parents' estimates (counts add across
## alleles, hybrid = one Sc-like plus one Se-like allele), a nonlinear
## function of log2(Sc) and log2(Se). In the AD frame SE bars use
## cor_dpar_<quantity> when present (exact for the rotation) and fall back
## to independence otherwise. The noise quantities plot the ploidy-adjusted
## estimates, so a cloud on the zero lines means hybrid noise matches the
## parents at the per-genome scale. Returns the class vector invisibly
## (finite genes only).
plot_dom_class <- function(BURST.CONTRASTS, PR, quantity = c("mean", "bfreq", "bsize", "kbal", "cv2"), frame = c("parent","AD"), sig = 0.05, main = NULL, lim = NULL, bar_col = adjustcolor(COLOR.GREY[["dark"]], 0.25)) {
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
  if (is.null(main)) main<-paste0("dominance (",quantity,", ",frame,")")
  .se_scatter(x, sx, y, sy, xlab = xl, ylab = yl, main = main, lim = lim, bar_col = bar_col, pt_col = COLOR.LIST.2[match(cls,DOM.CLASS)],
              legend_labels = DOM.CLASS, legend_cols = COLOR.LIST.2)
  invisible(cls)
}

## ============================================================
## Intrinsic / extrinsic noise from the two hybrid alleles
## ============================================================
## intrinsic_extrinsic_components(mats, expos): two-allele decomposition,
## pooling HYC and HYT cells. Each allele's depth-normalised rate is divided
## by its own mean, which removes the cis mean difference. Within a cell the
## two alleles share extrinsic fluctuations and differ by intrinsic ones, so
## extr = mean(a'b') - 1 (allele covariance) and intr = 0.5 mean((a'-b')^2).
## Poisson shot noise is subtracted from intr: Var_Poisson(a'_i) =
## a_i/(e_i^2 ma^2), and intr is clipped at zero because a variance estimate that falls
## below zero is sampling noise. Returns gene, intr, extr and a `good` flag (both allele
## means positive), so the components can also be used directly as a ratio
## (Section 12) as well as through the intrinsic fraction in analysis.R.
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
  ## A variance cannot be negative. Where the Poisson term exceeds the measured allele difference the
  ## estimate is sampling noise, so intr is set to zero.
  intr <- pmax(intr_raw - 0.5*(shot_a+shot_b), 0)
  data.frame(gene = rownames(mats$HYC.SC), intr = intr, extr = extr, good = good, row.names = NULL)
}

## Expected log2 rise in DISP (bfreq) when a diploid hybrid's two alleles are
## summed inside one cell and compared with a haploid parent.
## Derivation. Write allele k of a gene in one cell as its mean times
## (1 + e + i_k). The term e is the cell's extrinsic deviation and both alleles
## carry it, because both sit in the same cell and share its transcription
## factors, cell size and cell-cycle state. The term i_k is the allele's own
## private deviation from stochastic bursting. Private terms of the two
## alleles are independent, and the terms have mean zero. Let E = var(e) and
## I = mean of var(i_1), var(i_2). One allele then has latent CV2 = E + I,
## which is the noise a haploid genome with the same regulation shows.
## The sum of the alleles weights each by its share of the summed mean,
## w1 = m1 / (m1 + m2) and w2 = 1 - w1, so the summed deviation is
## e + w1 * i_1 + w2 * i_2. The shared term e passes through at full weight
## because the weights sum to one. The private terms are independent, so their
## variances add after squaring the weights, which gives
## var = E + (w1^2 + w2^2) * I. Because w1^2 + w2^2 is at most 1, the sum is
## quieter than either allele. It equals 1/2 for equal alleles and approaches 1
## when one allele dominates, since one allele leaves nothing to average.
## DISP is the inverse of latent CV2, so the rise in log2 DISP is
## s = -log2((E + (w1^2 + w2^2) * I) / (E + I))
##   = -log2(1 - (1 - w1^2 - w2^2) * phi), with phi = I / (I + E).
## Equal alleles with fully private noise give s = 1. Fully shared noise gives 0.
## E and I come from intrinsic_extrinsic_components(), which removes Poisson
## shot noise. The average I treats the two alleles' private variances as equal.
## Returns one value per gene, named by gene, in log2 units.
ploidy_shift <- function(mats, expos) {
  ie    <- intrinsic_extrinsic_components(mats, expos)
  e_vec <- c(expos$HYC, expos$HYT)
  ma <- rowMeans(sweep(cbind(mats$HYC.SC, mats$HYT.SC), 2, e_vec, "/"))   # Sc allele mean rate
  mb <- rowMeans(sweep(cbind(mats$HYC.SE, mats$HYT.SE), 2, e_vec, "/"))   # Se allele mean rate
  w1   <- ma / (ma + mb)
  wsq  <- w1^2 + (1 - w1)^2                          # weight on the private noise of the sum
  intr <- ie$intr; extr <- pmax(ie$extr, 0)           # intr is already clipped at zero
  phi  <- intr / (intr + extr)
  s    <- -log2(1 - (1 - wsq) * phi)
  s[!ie$good | !is.finite(s)] <- NA_real_
  setNames(s, ie$gene)
}

##############################################################################
## 10. PROMOTER ARCHITECTURE (TATA box, poly(dA:dT), and predicted
##     nucleosome occupancy)
##############################################################################
## Three sequence features describe promoter architecture. A TATA box
## consensus match and the longest poly(dA:dT) tract are direct sequence
## heuristics scored by functions in this file. Predicted nucleosome
## occupancy comes from NuPoP (Xi et al. 2010, Bioinformatics; extending
## the duration hidden Markov approach of Wang et al. 2008), whose profile
## is trained on yeast nucleosome data.
## All three features depend on sequence alone, so the same functions score
## Sc and Se promoters and no experimental occupancy map for Se is needed.
## Promoters are defined from each species' genome FASTA and GFF gene
## annotation, because neither species has an annotated promoter set.
##
## NuPoP is a Bioconductor package. Its trained profile (species = 7) is
## fit on S. cerevisiae and is applied to both species, so Sc and Se are
## scored against one shared reference. The Se scores therefore carry a
## cross-species extrapolation, to be stated in the methods alongside the
## Bucher matrix's cross-eukaryote origin.

## Reads a genome FASTA into a named DNAStringSet, one sequence per
## chromosome or scaffold. Names are the header text before any whitespace,
## which is how a GFF's seqid column references chromosomes.
read_genome_fasta <- function(path) {
  genome <- Biostrings::readDNAStringSet(path)
  names(genome) <- sub("\\s.*$", "", names(genome))
  genome
}

## Reads a GFF3 annotation and returns one row per gene: gene, seqid,
## start, end (1-based, start < end regardless of strand) and strand, which
## is what extract_promoters() needs to locate each promoter.
##
## Genes are anchored on their exon (coding) features. The exon starts at
## the ATG, the reference point of Basehoar et al.'s location window, and
## excludes annotated UTRs, whose completeness varies between independent
## annotations and is unrelated to promoter divergence.
##
## This Geneious-exported GFF3 has no Parent attribute linking an exon to
## its gene. An exon's Name is the gene name plus a literal "_CDS" suffix
## (true for 6302 of 6306 exon rows), so id_field = "Name" and
## id_suffix = "_CDS" recover the bare gene identifier; sub() leaves names
## without the suffix unchanged. Fragments of a multi-exon gene share one
## identifier and collapse to a single min(start)/max(end) span that
## includes any introns. Use id_suffix = "" for annotations without this
## convention. If id_field is not an attribute column, the function stops
## and lists the available columns.
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

## Extracts the promoter of every gene in genes: the sequence immediately
## upstream of its start codon, up to max_bp long or up to the nearest
## neighboring gene on the same chromosome, whichever is shorter, so a
## promoter stops at the neighbor's boundary. The boundary on the + strand is
## the largest end among all genes that start before the gene (a running
## maximum over the genes sorted by start, so a long gene upstream still
## bounds the promoter when a shorter gene lies between them); on the - strand, the mirror
## image, it is the smallest start among all genes that end after the gene (a gene nested
## inside it does not bound it, while a host gene that contains it does). A gene that lies
## inside another gene or overlaps its neighbor has no region left, which the min_bp rule below
## handles. A divergent gene pair gets
## its full shared intergenic region and tightly spaced genes get short
## promoters, both genuine features of the genome.
##
## min_bp is the one exception to the neighbor boundary. A region shorter
## than min_bp is too short to score a TATA box, poly(dA:dT) tract or
## NuPoP window, so extraction continues into the neighboring gene's
## coding sequence until min_bp is reached (or the chromosome ends). Such
## genes are flagged in the "extended" attribute, so scores from extended
## regions stay identifiable downstream. Genes with no usable region keep
## "" and NA coordinates.
##
## Every returned sequence reads 5' to 3' on its own gene's strand, with
## the end of the string nearest the start codon; minus-strand promoters
## are reverse-complemented. Orientation matters for score_promoters() (the
## consensus is not a palindrome) but not for the poly(dA:dT) run length (an A/T run
## has the same length on either strand). The "coords" attribute holds
## each promoter's genomic seqid, start, end (start < end) and strand.
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
    ## Neighbor boundaries over all genes, not only the adjacent one. prev_end[i] is the largest end among
    ## genes whose start is smaller than gene i's (0 when none; findInterval counts the genes that start
    ## strictly before), and next_start[i] the smallest start among genes whose end is larger than gene i's
    ## (Inf when none), a suffix minimum over the genes sorted by end. A value at or inside the gene
    ## leaves the promoter region empty.
    n_before   <- findInterval(g$start - 1, g$start)
    prev_end   <- c(0, cummax(g$end))[n_before + 1]
    by_end     <- order(g$end)
    n_le_end   <- findInterval(g$end, g$end[by_end])
    next_start <- c(rev(cummin(rev(g$start[by_end]))), Inf)[n_le_end + 1]
    for (i in seq_len(nrow(g))) {
      ext <- FALSE
      if (g$strand[i] == "+") {
        limit   <- prev_end[i] + 1
        p_end   <- g$start[i] - 1
        p_start <- max(limit, p_end - max_bp + 1, 1)
        if (p_end - p_start + 1 < min_bp) {
          p_start_ext <- max(p_end - min_bp + 1, 1)
          if (p_start_ext < p_start) { p_start <- p_start_ext; ext <- TRUE }
        }
        seq <- if (p_start > p_end) "" else
          as.character(Biostrings::subseq(chrom, p_start, p_end))
      } else {
        limit   <- if (is.finite(next_start[i])) next_start[i] - 1 else length(chrom)
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
## Coordinates are recorded only for nonempty promoters; genes with "" keep
## NA coordinates, which is how downstream code skips them. start/end are
## the lower/higher genomic position regardless of strand; strand tells a
## reader whether to read the span forward or reverse-complement it.
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
## consensus also uses, giving base percentages per position (divided by
## 100 here) for a continuous score. Bucher derived it from 502 eukaryotic
## RNA polymerase II promoters, so it is quantitative but not
## Saccharomyces-specific, a limitation to state in the methods.
TATA.PWM <- rbind(
  A = c(4.1, 90.5,  0.8, 91.0, 68.9, 92.5, 57.1, 39.8),
  T = c(79.5, 9.0, 96.1,  7.7, 31.1,  1.6, 31.1,  8.5),
  G = c(4.6,  0.5,  0.5,  1.3,  0.0,  5.1, 11.3, 40.4),
  C = c(11.8, 0.0,  2.6,  0.0,  0.0,  0.8,  0.5, 11.3)) / 100

## Basehoar et al. located functional TATA boxes 50 to 200 bp upstream of
## the ATG (their location criterion). TATA.PWM says what a TATA box looks
## like; this window says where to look for one.
TATA.WINDOW <- c(50, 200)

## Scores the TATA box and the longest poly(dA)/poly(dT) run of every sequence in a
## named vector and returns one row per gene, ready to merge across
## species. prom_len is the extracted promoter length, so short or empty
## promoters can be told apart from full-length ones that scored low.
## prom_extended carries extract_promoters()'s "extended" attribute
## (FALSE for all genes when absent), flagging promoters that reached the
## min_bp floor by crossing into a neighboring gene's coding sequence.
## tata_pos and polyat_pos give each best match's first base in bp
## upstream of the ATG (prom_len - string index + 1), a distance that is
## comparable between promoters of different lengths; NA propagates from
## the score. Comparing tata_pos or polyat_pos alongside tata_delta or
## polyat_delta separates an indel (same element at a shifted distance)
## from a substitution that changed which site scores best.
score_promoters <- function(seqs, tata_window = TATA.WINDOW) {
  ## Scans the Basehoar location window of a promoter sequence and scores
  ## every 8-mer against TATA.PWM as a log2 odds ratio relative to a uniform
  ## 25% base composition, returning the best-scoring 8-mer as a continuous
  ## score. A pseudocount of 0.001 floors the matrix's zero cells, so one mismatch at an
  ## otherwise conserved position lowers the score without forcing -Inf.
  ## Each 8-mer lies entirely within the window. Directional: seq must be
  ## oriented 5' to 3' relative to its own gene, with the promoter-proximal
  ## end at the end of the string (see extract_promoters()). position is the
  ## string index of the best 8-mer's first base and motif is its sequence.
  ## Windows containing N or other non-ACGT bases score NA; ties go to the
  ## most upstream 8-mer.
  tata     <- lapply(seqs, function(seq) {
    motif_len <- ncol(TATA.PWM)
    n <- nchar(seq)
    idx_lo <- max(1, n - tata_window[2] + 1)
    idx_hi <- n - tata_window[1] - motif_len + 2
    if (idx_hi < idx_lo) return(list(score = NA, position = NA, motif = NA))
    starts <- idx_lo:idx_hi
    pwm_adj <- pmax(TATA.PWM, 0.001)
    scores <- sapply(starts, function(i) {
      bases <- strsplit(substring(seq, i, i + motif_len - 1), "")[[1]]
      row_idx <- match(bases, rownames(pwm_adj))
      if (length(bases) < motif_len || any(is.na(row_idx))) return(NA)
      sum(log2(pwm_adj[cbind(row_idx, seq_len(motif_len))] / 0.25))
      })
      if (all(is.na(scores))) return(list(score = NA, position = NA, motif = NA))
      best <- which.max(scores)
      list(score = scores[best], position = starts[best], motif = substring(seq, starts[best], starts[best] + motif_len - 1))
  })
  ## Finds the longest poly(dA) or poly(dT) run in the sequence: a run is one repeated base,
  ## so A and T runs are measured separately and the longer wins, and a mixed stretch such as
  ## ATATAT does not count as a tract. These homopolymer runs, the two strands of a
  ## poly(dA:dT) tract, are the feature behind poly(dA:dT)-mediated nucleosome exclusion. The
  ## search is strand-symmetric, so it runs on the promoter as extracted. An empty sequence
  ## gives NA (no promoter to measure), so a missing promoter is not read as a measured length
  ## of 0 in polyat_delta. position is the string index of the tract's first base and tract is
  ## its sequence (its first base says whether it is a dA or a dT run).
  polyat   <- lapply(seqs, function(seq) {
    if (nchar(seq) == 0) return(list(length = NA, position = NA, tract = NA))
    runs <- gregexpr("A+|T+", seq)[[1]]
    lens <- attr(runs, "match.length")
    if (runs[1] == -1) return(list(length = 0, position = NA, tract = NA))
    best <- which.max(lens)
    list(length = lens[best], position = runs[best], tract = substring(seq, runs[best], runs[best] + lens[best] - 1))
  })
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

## Predicted nucleosome occupancy runs on the Linux cluster, in the same
## package-compute-load pattern as the bootstrap and permutation jobs.
## analysis.R Section 6.1 packages each species' chromosomes and promoter
## coordinates with nupop_cluster_inputs(), cluster/scripts/nupop_occupancy.R
## runs nupop_occupancy_cluster(), and Section 6.1 scores the returned
## tracks with score_promoters_nupop().
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

## Returns one row per gene (gene, occ_score) from the occupancy tracks
## the cluster job produced (NUPOP.OCC.SC or NUPOP.OCC.SE, loaded from
## nupop_output.rda); merge it onto score_promoters()'s table by gene.
## seqs must be extract_promoters()'s return value. The coordinate check
## confirms the tracks were computed from the same promoters as seqs, so a
## change to PROM.MAX.BP, PROM.MIN.BP or the annotation is caught here.
## Genes in an unscored region carry NA, and the region tables travel with
## the scores as attributes.
score_promoters_nupop <- function(seqs, occ_list, window = TATA.WINDOW) {
  coords <- attr(seqs, "coords")
  if (is.null(coords))
    stop("score_promoters_nupop() reads coords from the 'coords' attribute of extract_promoters()'s return value. Pass PROM.SC or PROM.SE directly as returned in Section 6.1.")
  coords <- coords[!is.na(coords$start), c("gene", "seqid", "start", "end", "strand")]
  if (!isTRUE(all.equal(coords, attr(occ_list, "coords"), check.attributes = FALSE)))
    stop("score_promoters_nupop(): these promoters differ from the ones the cluster job scored. Save nupop_inputs.rda again in Section 6.1 and rerun nupop_occupancy.R.")
  ## Reduces each gene's occupancy track to the mean predicted occupancy
  ## over the Basehoar window used for the TATA score (50 to 200 bp
  ## upstream of the ATG), so all three architecture features read the
  ## same stretch of promoter. A window made entirely of NA positions
  ## returns NA.
  occ_score <- local({
    vapply(names(seqs), function(g) {
      n      <- nchar(seqs[[g]])
      occ    <- occ_list[[g]]
      idx_lo <- max(1, n - window[2] + 1)
      idx_hi <- min(n, n - window[1] + 1)
      if (n == 0 || idx_hi < idx_lo || length(occ) < n) return(NA_real_)
      vals <- occ[idx_lo:idx_hi]
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }, numeric(1))
  })
  scores    <- data.frame(gene = names(seqs), occ_score = unname(occ_score),
                          row.names = NULL, stringsAsFactors = FALSE)
  attr(scores, "reduced_flank_regions") <- attr(occ_list, "reduced_flank_regions")
  attr(scores, "failed_regions")        <- attr(occ_list, "failed_regions")
  scores
}

## Regulatory classes with a cis component. Cis, Cis + Trans and Compensatory all require a
## significant cis permutation p-value by construction (classify_reg()).
.CIS_CLASSES <- c("Cis", "Cis + Trans", "Compensatory")

## Predicted sign of each tested quantity's cis contrast relative to the promoter shift delta (Se - Sc), the one
## prespecified direction the concordance tests use. bfreq is the primary test: a higher TATA score or a longer
## poly(dA:dT) tract goes with lower DISP (Section 6.2), and bfreq_cis_est = log2(DISP_Sc) - log2(DISP_Se), so
## sign(delta) should equal sign(bfreq_cis_est) (+1). bsize is read from the identity log BSIZE = log MU - log DISP:
## a feature that acts through DISP at fixed mean moves log BSIZE the opposite way, so sign(delta) should equal
## -sign(bsize_cis_est) (-1). The bsize test follows from the bfreq prediction and the mean contrast enters it, so
## it is a secondary check that carries no independent mechanism. kbal has no independent prediction and is not tested.
.PROMOTER_PREDICTED_SIGN <- c(bfreq = 1, bsize = -1)

## Genes eligible for the promoter concordance tests: in a class with a cis component, with the
## promoter shift (delta) and the cis estimate both defined, and neither exactly 0 (a zero has no
## sign to match). Shared by promoter_direction_test() and concordance_by_magnitude().
.concordance_eligible <- function(delta, est, reg_class, cis_classes = .CIS_CLASSES) {
  reg_class %in% cis_classes & !is.na(delta) & !is.na(est) & sign(delta) != 0 & sign(est) != 0
}

## Tests whether the direction of each promoter feature's between-species
## shift matches the direction predicted from the per-species relationship
## in NOISE.VALIDATE (Section 6.2): a higher TATA score or longer
## poly(dA:dT) tract goes with lower DISP (burstier, noisier expression).
## bfreq_cis_est = log2(DISP_Sc) - log2(DISP_Se), so a positive value means
## Se has the lower DISP and is the noisier allele, and the prediction is
## that Se carries the larger feature value, delta (Se - Sc) > 0. A gene is
## concordant when sign(delta) == .PROMOTER_PREDICTED_SIGN[[quantity]] *
## sign(<quantity>_cis_est). bfreq is the primary, prespecified test. For
## quantity "bsize" the predicted sign is flipped (log BSIZE = log MU - log DISP),
## so concordance there asks whether delta has the sign opposite to bsize_cis_est;
## it is a secondary test that follows from the bfreq prediction (see
## .PROMOTER_PREDICTED_SIGN). kbal is not tested. The result carries the
## quantity and the predicted sign beside each feature.
##
## The test runs on the full any-cis gene set (cis_classes), without the
## arch_frac magnitude cutoff used by promoter_noise_candidates(), because
## selecting on magnitude first would cost power. It is a one-sided
## binomial sign test against 50/50 (alternative = "greater"): if promoter
## architecture is unrelated to which allele is noisier, concordant and
## discordant genes are equally likely. Genes with delta or estimate equal
## to 0 are dropped. A non-significant result means the candidate tables
## in Section 6.5 should be read with caution.
promoter_direction_test <- function(BURST.CONTRASTS, PR, ARCH, reg_class, quantity = c("bfreq", "bsize"), cis_classes = .CIS_CLASSES) {
  quantity <- match.arg(quantity)
  est <- BURST.CONTRASTS[[paste0(quantity, "_cis_est")]]
  predicted_sign <- .PROMOTER_PREDICTED_SIGN[[quantity]]
  one_feature <- function(delta) {
    ok <- .concordance_eligible(delta, est, reg_class, cis_classes)
    n          <- sum(ok)
    concordant <- sign(delta[ok]) == predicted_sign * sign(est[ok])
    n_conc     <- sum(concordant)
    bt <- if (n > 0) binom.test(n_conc, n, p = 0.5, alternative = "greater")
          else list(p.value = NA_real_)
    data.frame(quantity = quantity, predicted_sign = predicted_sign, n = n, n_concordant = n_conc,
               frac_concordant = if (n > 0) n_conc / n else NA_real_,
               p = bt$p.value)
  }
  rbind(cbind(feature = "TATA", one_feature(ARCH$tata_delta)), cbind(feature = "poly(dA:dT)", one_feature(ARCH$polyat_delta)), cbind(feature = "nucleosome occupancy (NuPoP)", one_feature(ARCH$occ_access_delta)))
}

## Splits the any-cis gene set used by promoter_direction_test() into
## n_bins equal-count bins by |delta| and reports the fraction concordant
## with the predicted sign of the cis estimate in each bin (.PROMOTER_PREDICTED_SIGN[[quantity]]). An effect confined to the most
## divergent promoters is diluted in the whole-set test but appears here as
## concordance rising from the smallest to the largest |delta| bin; no
## rise suggests there is no signal at any magnitude. Eligibility
## (cis_classes membership, both values defined, neither exactly 0) is
## rebuilt here so the function works on its own. est is the cis estimate of
## the tested quantity (bfreq_cis_est or bsize_cis_est).
##
## ci_lo/ci_hi are a normal-approximation 95% interval on each bin's
## proportion, adequate for the bin sizes this produces and meant only
## to show roughly how much a bin's estimate should be trusted, not as
## a formal per-bin test. trend_coef and trend_p (attached as
## attributes, since they describe the whole feature rather than any
## one bin) come from a logistic regression of concordance on
## log(|delta|), the formal version of "is the rate rising toward the
## extremes" that the binned table and its plot show informally.
concordance_by_magnitude <- function(delta, est, reg_class, quantity = c("bfreq", "bsize"), cis_classes = .CIS_CLASSES, n_bins = 10) {
  quantity <- match.arg(quantity)
  ok <- .concordance_eligible(delta, est, reg_class, cis_classes)

  d    <- abs(delta[ok])
  conc <- sign(delta[ok]) == .PROMOTER_PREDICTED_SIGN[[quantity]] * sign(est[ok])
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
  attr(agg, "quantity")   <- quantity
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
  bp <- barplot(cb$frac_concordant, ylim = c(0, 1), col = COLOR.GREY[["mid"]], names.arg = round(cb$mag_mean, 1), las = 2, ylab = "fraction concordant", xlab = "mean |delta| in bin", main = main)
  segments(bp, cb$ci_lo, bp, cb$ci_hi)
  abline(h = 0.5, lty = 2, col = COLOR.ACCENT)
}

## Lists genes worth inspecting by eye: a cis component of noise
## divergence together with a large shift in a promoter feature. Neither
## alone is informative: a promoter difference without a cis noise
## signature is unexplained by this section, and cis divergence without a
## promoter difference is not explained by anything measured here.
##
## Returns three tables (tata, polyat, occ). A gene needs a large shift in
## only one feature to qualify and may appear in several tables.
##
## Gating is on regulatory class. Cis, Cis + Trans and Compensatory all
## require a significant cis permutation p-value by construction
## (classify_reg(), Section 3.1), so each candidate has a tested cis
## effect rather than only a large total point estimate. Trans-only,
## Conserved and Ambiguous genes are excluded.
##
## "Large" is defined by rank, because a single-genome PWM score or
## tract-length difference has no null distribution: a gene qualifies when
## its |delta| is at or above the arch_frac quantile of |delta| over all
## genes (arch_frac = 0.90 keeps roughly the top 10%). Each feature's
## cutoff is computed on its own values, so the tables differ in size, and
## integer-valued polyat_delta can keep more than 10% through ties.
##
## Permutation p-values cannot resolve below 1/(N.PERM + 1), so many genes
## share the minimum value; the ranking below sorts on effect size and is
## unaffected.
##
## BURST.CONTRASTS, PR, ARCH and reg_class must share one gene order
## (true for BURST.CONTRASTS, PR and the REG.*.CLASS vectors by
## construction, and for ARCH after realignment to BURST.CONTRASTS$gene in
## Section 6.3). Each table is sorted by |<quantity>_cis_est|: cis compares
## the two alleles within the hybrid in one trans environment, the
## allele-specific quantity a promoter difference should track. Cis, trans
## and total estimates and p/q-values are all kept, so the reg_class call
## can be reconciled. concordant flags whether the gene's delta sign
## matches the predicted direction tested in promoter_direction_test() (the predicted sign of
## the quantity, .PROMOTER_PREDICTED_SIGN); a discordant
## gene can be a valid example but merits extra scrutiny. Raw values for
## both species (score or length, position, matched sequence) precede each
## feature's deltas, so an indel can be told apart from a substitution by
## reading the sequences in the table.
promoter_noise_candidates <- function(BURST.CONTRASTS, PR, ARCH, reg_class, quantity = c("bfreq", "bsize"), cis_classes = .CIS_CLASSES, arch_frac = 0.90) {
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
  q_cis     <- PR[[paste0(quantity, "_cis_q")]]
  q_trans   <- PR[[paste0(quantity, "_trans_q")]]
  q_total   <- PR[[paste0(quantity, "_total_q")]]

  cis_flag <- reg_class %in% cis_classes

  est_col <- paste0(quantity, "_cis_est")

  ## One table per promoter feature. delta names the shift that gates membership (at or above
  ## the arch_frac quantile of |delta| over all genes, computed on that feature's own values) and
  ## sets `concordant`; cols are the raw species values and deltas that precede it. occ_score_sc/_se
  ## are NuPoP's own occupancy scale (a 0-1 probability that a base is nucleosome-covered), kept
  ## beside occ_access_delta so a candidate can be checked against the occupancy track; concordant
  ## uses occ_access_delta, matching the sign convention used everywhere else in this section.
  feature_spec <- list(
    tata   = list(delta = "tata_delta",
                  cols = c("tata_score_sc", "tata_score_se", "tata_pos_sc", "tata_pos_se",
                           "tata_motif_sc", "tata_motif_se", "tata_delta", "tata_pos_delta")),
    polyat = list(delta = "polyat_delta",
                  cols = c("polyat_len_sc", "polyat_len_se", "polyat_pos_sc", "polyat_pos_se",
                           "polyat_tract_sc", "polyat_tract_se", "polyat_delta", "polyat_pos_delta")),
    occ    = list(delta = "occ_access_delta",
                  cols = c("occ_score_sc", "occ_score_se", "occ_delta", "occ_access_delta")))
  prom_cols <- c("prom_len_sc", "prom_len_se", "prom_extended_sc", "prom_extended_se")

  lapply(feature_spec, function(spec) {
    delta <- ARCH[[spec$delta]]
    cut   <- quantile(abs(delta), arch_frac, na.rm = TRUE)
    keep  <- cis_flag & !is.na(delta) & abs(delta) >= cut
    keep[is.na(keep)] <- FALSE
    ## Columns shared by all three tables: gene, regulatory class and the cis, trans and total
    ## estimates, p-values and q-values. Names carry the quantity prefix (bfreq_ or bsize_) so a
    ## table states which axis it was built from.
    ctx <- data.frame(
      gene      = BURST.CONTRASTS$gene[keep],
      reg_class = reg_class[keep],
      stringsAsFactors = FALSE)
    ctx[[paste0(quantity, "_cis_est")]]   <- est_cis[keep];   ctx[[paste0(quantity, "_cis_p")]]   <- p_cis[keep]
    ctx[[paste0(quantity, "_trans_est")]] <- est_trans[keep]; ctx[[paste0(quantity, "_trans_p")]] <- p_trans[keep]
    ctx[[paste0(quantity, "_total_est")]] <- est_total[keep]; ctx[[paste0(quantity, "_total_p")]] <- p_total[keep]
    ## FDR-adjusted values sit beside the raw p-values so a reader can see both
    ctx[[paste0(quantity, "_cis_q")]] <- q_cis[keep]; ctx[[paste0(quantity, "_trans_q")]] <- q_trans[keep]; ctx[[paste0(quantity, "_total_q")]] <- q_total[keep]
    arch_feature <- as.data.frame(ARCH[spec$cols], stringsAsFactors = FALSE)[keep, , drop = FALSE]
    arch_prom    <- as.data.frame(ARCH[prom_cols], stringsAsFactors = FALSE)[keep, , drop = FALSE]
    row.names(arch_feature) <- row.names(arch_prom) <- NULL
    out <- cbind(ctx, arch_feature,
                 data.frame(concordant = sign(delta[keep]) == .PROMOTER_PREDICTED_SIGN[[quantity]] * sign(est_cis[keep])),
                 arch_prom)
    out[order(-abs(out[[est_col]])), ]
  })
}

## Tests whether a promoter feature associates with a species' own
## absolute noise level, independent of any cross-species comparison.
## fit is one of SPLIT.FITS' per-species data frames (Section 2);
## response names the column to test, DISP (burst frequency, the
## default and the axis used throughout the rest of this section) or
## BSIZE (MU/DISP, already computed on every fit object by
## fit_counts_offset()); score is a named vector of a single
## feature (tata_score or polyat_len) for the same species, keyed the
## same way as fit's row names.
##
## covariate names the column held fixed while testing the feature, MU by
## default. With response = "BSIZE" pass covariate = "DISP": since
## BSIZE = MU/DISP, testing BSIZE against covariate MU is algebraically
## forced to reproduce the DISP-against-MU result with the coefficient sign
## flipped, log(BSIZE) = log(MU) - log(DISP),
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
## hvg_elbow: data-driven nfeatures for FindVariableFeatures(). The ranked
## standardized-variance (vst) curve falls steeply over the few strongly
## variable genes and flattens into the bulk; the elbow is the last
## rank-to-rank drop larger than drop_frac of the curve's range, plus one,
## the same style as the PC elbow used for dimensionality selection. obj is
## refit with nfeatures = all genes so the search sees the complete ranking.
##
## The search covers genes with standardized variance >= floor (default 1,
## the fitted mean-variance trend itself); genes below the trend are
## under-dispersed and are not candidates for top variable genes, so the
## elbow is located within the above-trend part of the curve.
## Returns n_features (the chosen count), the full sorted curve (so
## plot_hvg_elbow() can show the cutoff on the whole ranking) and floor.
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
  on.exit(dev.off(), add = TRUE)
  plot(hvg$curve, pch = 16, cex = 0.4, col = ifelse(seq_along(hvg$curve) <= hvg$n_features, "black", COLOR.GREY[["mid"]]),
       xlab = "gene rank (by standardized variance)", ylab = "standardized variance", main = sprintf("%s: HVG elbow at %d genes", label, hvg$n_features))
  abline(v = hvg$n_features, lty = 2, col = COLOR.ACCENT)
  abline(h = hvg$floor, lty = 3, col = COLOR.ACCENT)
  legend("topright", legend = c("chosen cutoff", sprintf("floor (%.1f)", hvg$floor)), lty = c(2, 3), col = c(COLOR.ACCENT, COLOR.ACCENT), bty = "n")
}

## sweep_cluster_resolution: chooses the Louvain resolution for one dataset.
## At each resolution in res_grid it clusters obj on the given PCs, skips
## any resolution that gives a single cluster or a cluster below min_cells
## (the size floor is applied before silhouette is consulted), and scores
## the rest by mean silhouette width (cluster::silhouette) on the same
## PCA-space distance used for clustering. It returns the COARSEST
## resolution whose silhouette is within tol of the maximum: silhouette
## curves are often flat near their top, and splitting a real cluster into
## arbitrary halves does not raise silhouette, so resolution increases only
## when it buys a real gain in separation.
## Returns chosen_res, chosen_sil, chosen_n_clusters, obj (the clustered
## Seurat object at that resolution, reusable without reclustering) and
## grid (one row per resolution, for plotting).
sweep_cluster_resolution <- function(obj, dims, res_grid = seq(0.05, 1, by = 0.05), min_cells = 50, tol = 0.02, metric = "manhattan") {
  emb <- Embeddings(obj, "pca")[, dims, drop = FALSE]
  d   <- dist(emb, method = metric)
  obj <- FindNeighbors(obj, reduction = "pca", dims = dims, annoy.metric = metric, verbose = FALSE)

  grid <- lapply(res_grid, function(r) {
    o2    <- FindClusters(obj, resolution = r, verbose = FALSE)
    cl    <- Idents(o2)
    sizes <- table(cl)
    # A single-cluster result passes the size floor but leaves silhouette
    # undefined (no second cluster to compare against), so it is excluded
    # together with results that break min_cells.
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
  on.exit(dev.off(), add = TRUE)
  plot(g$res[ok], g$sil[ok], type = "b", pch = 16, xlab = "resolution", ylab = "mean silhouette width",
       main = sprintf("%s: resolution sweep", label), ylim = range(g$sil[ok], na.rm = TRUE))
  if (any(!ok)) points(g$res[!ok], rep(min(g$sil[ok], na.rm = TRUE), sum(!ok)), pch = 4, col = COLOR.GREY[["mid"]])
  abline(v = sweep$chosen_res, lty = 2, col = COLOR.ACCENT)
  legend("bottomright", legend = c("silhouette", sprintf("below %d cells/cluster", min_cells), "chosen"), pch = c(16, 4, NA), lty = c(NA, NA, 2), col = c("black", COLOR.GREY[["mid"]], COLOR.ACCENT), bty = "n")
}

## prepare_dataset: normalization, variable features and PCA for one Seurat dataset. Counts are
## log-normalized; the number of variable features comes from the elbow of the ranked standardized-variance
## curve (hvg_elbow, plotted by plot_hvg_elbow); every gene in `genes` is scaled; PCA runs on the variable
## features. Returns the prepared object and the elbow result.
prepare_dataset <- function(obj, name, genes, figure_dir) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  hvg <- hvg_elbow(obj)
  plot_hvg_elbow(hvg, name, figure_dir)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvg$n_features, verbose = FALSE)
  obj <- ScaleData(obj, features = genes, verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = FALSE)
  list(obj = obj, hvg = hvg)
}

## elbow_pcs: number of PCs to retain, the elbow of the ranked percent-variance curve: the latest PC among
## consecutive PCs whose percent-variance drop exceeds 0.05 percentage points, plus one. Returns the percent
## variance (pct), its cumulative sum (cumu) and the PC count (pcs).
elbow_pcs <- function(obj) {
  pct <- 100 * obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev)
  list(pct = pct, cumu = cumsum(pct), pcs = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = TRUE)[1] + 1)
}

## cluster_dataset: resolution sweep for dataset d on its retained PCs (pcs[[d]]), the sweep plot, and a log
## line naming the chosen resolution (labels[[d]] is the dataset's readable name). Returns the
## sweep_cluster_resolution() result, whose obj carries the chosen clustering.
cluster_dataset <- function(d, ysc, pcs, labels, min_cells, figure_dir) {
  sweep <- sweep_cluster_resolution(ysc[[d]], dims = 1:pcs[[d]], min_cells = min_cells)
  plot_resolution_sweep(sweep, d, figure_dir, min_cells)
  cat(sprintf("%s: chosen resolution = %.2f, %d clusters, mean silhouette = %.3f\n",
              labels[[d]], sweep$chosen_res, sweep$chosen_n_clusters, sweep$chosen_sil))
  sweep
}

## umap_dataset: UMAP on the retained PCs of one clustered dataset, drawn with the shared cluster palette on
## the open device. Returns the object with its UMAP reduction.
umap_dataset <- function(obj, pcs, title) {
  obj <- suppressWarnings(RunUMAP(obj, dims = 1:pcs, verbose = FALSE))
  print(umap_plot(obj, title, group.by = "seurat_clusters"))
  obj
}

## assemble_cluster_stability: summarizes one dataset's returned ARI vectors.
## table has one row per candidate resolution (sil, n_clusters, mean/min/max
## ARI, n_ok = completed replicates); final_res is the candidate with the
## highest mean ARI (ties go to the first row, the sweep's chosen
## resolution); final_obj is obj relabelled with the reference partition at
## final_res; ari keeps the raw vectors. The comparison report in analysis.R
## Section 7.3 and the downstream code read this list.
assemble_cluster_stability <- function(inputs, ari, dataset, obj) {
  rows <- inputs$tasks[inputs$tasks$dataset == dataset, ]
  a    <- ari[rows$task]
  table <- data.frame(res = rows$res, role = rows$role, sil = rows$sil, n_clusters = rows$n_clusters,
                      boot_mean_ari = vapply(a, mean, numeric(1), na.rm = TRUE),
                      boot_min_ari  = vapply(a, min,  numeric(1), na.rm = TRUE),
                      boot_max_ari  = vapply(a, max,  numeric(1), na.rm = TRUE),
                      n_ok = vapply(lapply(a, is.finite), sum, integer(1)), row.names = NULL)
  win <- which.max(table$boot_mean_ari)
  ## apply_cluster_labels: installs a partition on a Seurat object as its
  ## identities, the same state FindClusters() leaves behind.
  list(table = table, final_res = rows$res[win], final_obj = local({
    labels <- inputs$ref[[rows$task[win]]]
    stopifnot(identical(names(labels), colnames(obj)))
    Idents(obj) <- labels
    obj$seurat_clusters <- labels
    obj
  }), ari = a)
}

## ---- Shared helpers for cluster round trips ----
## kegg_local: downloads the organism's KEGG pathway map once, locally,
## so cluster nodes run KEGG enrichment offline through enricher() with
## the same pathway release used everywhere in the run.
kegg_local <- function(org = "sce") {
  kg <- download_KEGG(org)
  if (!"to" %in% names(kg$KEGGPATHID2NAME))
    stop("kegg_local(): KEGGPATHID2NAME has no 'to' (pathway name) column; found: ", paste(names(kg$KEGGPATHID2NAME), collapse = ", "))
  kg$KEGGPATHID2NAME[["to"]] <- sub(" - Saccharomyces cerevisiae \\(budding yeast\\)$", "", kg$KEGGPATHID2NAME[["to"]])
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

## load_replicates: loads the outputs of a replicate bootstrap, one file <stem>_output_<k>of<K>.rda per array
## task k (replicate k, seeded SEED + k - 1), and returns the object `object` from each as a list in replicate
## order. Every replicate must have size(x) == expected, so a stale or mismatched output stops here, before it
## reaches downstream sections. Replicate 1 is the reported result; the others check seed adequacy.
load_replicates <- function(stem, object, n_rep, expected, size = nrow, script = paste0(stem, ".R"), dir = OUTPUT.DIR) {
  if (!is.numeric(n_rep) || length(n_rep) != 1 || n_rep < 1) stop("n_rep must be a positive integer")
  reps <- vector("list", n_rep)
  for (k in seq_len(n_rep)) {
    file <- sprintf("%s_output_%dof%d.rda", stem, k, n_rep)
    env <- new.env()
    load_cluster_output(file.path(dir, file), script, envir = env)
    reps[[k]] <- env[[object]]
    if (size(reps[[k]]) != expected)
      stop(sprintf("%s has %d rows but %d are expected; rerun %s with the current inputs", file, size(reps[[k]]), expected, script), call. = FALSE)
  }
  reps
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
## order as mat). clusters: per-cell cluster label (same order, no NA).
## Clusters with fewer than min_cells cells (at least 2, the least a variance
## needs) are left out with a message, and the weights and the gene means
## are computed over the remaining cells.
##
## Within-cluster variance uses the sample divisor n_c - 1, so small clusters do not
## understate it. Shot noise is subtracted from the within-cluster component:
## per-cell Poisson sampling noise inflates within-cluster variance. The
## between-cluster component is the variance of the cluster means, and each
## mean carries its own sampling noise (the cluster's total within variance
## divided by n_c), which is subtracted so small clusters do not inflate the
## between-cluster variance. Both corrections make the pair an unbiased
## estimate rather than an exact partition of the total variance.
## var_within can come out slightly negative for a gene whose true
## biological within-cluster variance is near zero (shot noise
## subtraction is only unbiased on average, not gene by gene); it is
## returned as-is for transparency and only floored at zero where used
## as a ratio, the same convention intrinsic_extrinsic_components() uses
## for extr. var_between can likewise come out at or below zero after its
## correction; the ratio is then NA.
##
## Returns a list: table (one row per gene: var_within, var_between,
## ratio_within_between) and cluster_n (one row per cluster: cell
## count and whether the cluster entered the decomposition), so cluster
## sizes travel with the result rather than being dropped after this step.
within_between_decomp <- function(mat, expo, clusters, min_cells = 2) {
  stopifnot(ncol(mat) == length(expo), ncol(mat) == length(clusters), !anyNA(clusters),
            is.numeric(min_cells), length(min_cells) == 1, min_cells >= 2)
  cl        <- as.character(clusters)
  n_all     <- table(cl)
  small     <- names(n_all)[n_all < min_cells]
  if (length(small) > 0) {
    message(sprintf("within_between_decomp(): leaving out %d cluster(s) with fewer than %d cells (%s)", length(small), min_cells, paste(small, collapse = ", ")))
    use  <- !(cl %in% small)
    mat  <- mat[, use, drop = FALSE]; expo <- expo[use]; cl <- cl[use]
  }
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
    n     <- length(idx)
    w     <- n / N
    cmean <- rowMeans(ap[, idx, drop = FALSE])
    cvar  <- rowSums((ap[, idx, drop = FALSE] - cmean)^2) / (n - 1)
    cshot <- rowMeans(shot_cell[, idx, drop = FALSE])
    raw_within  <- raw_within + w * (cvar - cshot)
    var_between <- var_between + w * ((cmean - 1)^2 - cvar / n)
  }
  var_within <- raw_within
  ratio <- pmax(var_within, 0) / var_between
  ratio[!ok_gene | !is.finite(ratio) | var_between <= 0] <- NA_real_

  list(
    table = data.frame(gene = rownames(mat), var_within = var_within, var_between = var_between,
                        ratio_within_between = ratio, row.names = NULL),
    cluster_n = data.frame(cluster = names(n_all), n_cells = as.integer(n_all), used = !(names(n_all) %in% small), row.names = NULL)
  )
}

## Genome-wide distribution of within_between_decomp()'s ratio for one
## dataset. Values above 1 (log10 > 0) mean within-cluster variance
## exceeds between-cluster variance for that gene, i.e. residual noise
## among cells sharing a cluster outweighs the variance attributable to
## cluster identity. A large fraction of genes above 1 means the clustering
## captures a modest share of total variance (graded rather than sharply
## separated cell states). The fraction above 1 is taken over every gene with
## a finite ratio. Genes whose within-cluster variance is zero after the shot-noise
## subtraction (ratio 0, the most between-dominated end) cannot be placed on a log
## axis, so they are counted in the title and left out of the histogram.
## Returns n (genes with a finite ratio), n_zero and the fraction above 1 invisibly for logging.
plot_within_between_hist <- function(wb, main = NULL, brk = 40) {
  r <- wb$table$ratio_within_between
  r <- r[is.finite(r)]
  n_zero <- sum(r == 0)
  frac_above_1 <- mean(r > 1)
  lx <- log10(r[r > 0])
  h <- hist(lx, breaks = brk, plot = FALSE)
  if (is.null(main)) main <- sprintf("n = %d genes (%d with within = 0 not shown), %.0f%% with within > between", length(r), n_zero, 100 * frac_above_1)
  plot(h, col = COLOR.GREY[["light"]], border = "white", xlab = "log10(within / between)", ylab = "# of genes", main = main)
  abline(v = 0, lty = 2, col = COLOR.ACCENT)
  invisible(list(n = length(r), n_zero = n_zero, frac_above_1 = frac_above_1))
}

## Scatter of within_between_decomp()'s ratio (log10) against one burst
## kinetics quantity, for one dataset. quantity is a named vector (names
## = gene, already on the desired plotting scale, e.g. log2(mean)) so the
## same function covers mean, burst frequency, and burst size without
## three near-duplicate copies. The comparison is by rank (Spearman),
## since the ratio can span orders of magnitude without a linear
## relationship to the burst quantity; the dashed line is a least-squares
## guide on the plotted scale. rho uses every gene with a finite ratio, zeros
## included (ranks need no log); genes with ratio 0 cannot be placed on the log
## axis and are left out of the points and the guide line, and counted in the title.
## Returns n (genes in rho), n_zero and rho invisibly.
plot_within_between_vs_quantity <- function(wb, quantity, xlab, main = NULL) {
  idx <- match(wb$table$gene, names(quantity))
  ok  <- is.finite(wb$table$ratio_within_between) & !is.na(idx) & is.finite(quantity[idx])
  r   <- wb$table$ratio_within_between[ok]; x <- quantity[idx[ok]]
  rho <- suppressWarnings(cor(x, r, method = "spearman"))
  pos <- r > 0
  n_zero <- sum(!pos)
  if (is.null(main)) main <- sprintf("n = %d (%d with within = 0 not shown), Spearman rho = %.3f", length(r), n_zero, rho)
  y <- log10(r[pos])
  plot(x[pos], y, pch = 16, cex = 0.4, col = adjustcolor(COLOR.GREY[["dark"]], 0.38), xlab = xlab, ylab = "log10(within / between)", main = main)
  abline(lm(y ~ x[pos]), col = COLOR.ACCENT, lty = 2)
  invisible(list(n = length(r), n_zero = n_zero, rho = rho))
}

## Runs plot_within_between_by_class() for two datasets side by side
## (e.g. Sc and Se parent) against one classification axis, in one call:
## aligns the raw class vector (indexed like BURST.CONTRASTS$gene, e.g.
## REG.VEC$bfreq or DOM.VEC$bsize, not pre-cleaned since Ambiguous
## genes are already excluded automatically by class_levels not
## containing "Ambiguous") to each dataset's own gene table via match(),
## writes both panels to one pdf, and prints the omnibus test and Tukey
## pairwise comparisons (when significant) for both datasets. Returns
## both class_anova() results invisibly, named by lab_a/lab_b.
report_within_between_by_class <- function(wb_a, wb_b, class, gene_ref, class_levels, colors,
                                            lab_a, lab_b, axis_label, pdf_path, width = 11, height = 5.5) {
  ## Boxplot of within_between_decomp()'s ratio (log10), split by a
  ## regulatory or dominance class, for one dataset, with class_anova()'s
  ## omnibus test run internally and returned invisibly (F, df, p, eta^2,
  ## and Tukey pairwise comparisons if the omnibus test clears sig). class
  ## should already be aligned to wb$table$gene (match() against
  ## BURST.CONTRASTS$gene) with Ambiguous genes set to NA; class_levels
  ## fixes the plotting/level order (e.g. REG.CLASS or DOM.CLASS) and
  ## colors should be the matching per-level color vector (e.g.
  ## COLOR.LIST.1 or COLOR.LIST.2). main is left to the caller (e.g. a
  ## species label); the test statistics are returned for separate
  ## reporting rather than crowding the plot title.
  plot_within_between_by_class <- function(wb, class, class_levels, colors, main = NULL) {
    x  <- log10(wb$table$ratio_within_between)
    ok <- is.finite(x) & !is.na(class)
    ## genes with ratio 0 have no log value; they are counted on the console, not dropped silently
    n_zero <- sum(wb$table$ratio_within_between == 0 & !is.na(class), na.rm = TRUE)
    if (n_zero > 0) message(sprintf("within/between by class: %d gene(s) with within = 0 left out of the log10 comparison", n_zero))
    x  <- x[ok]; cl <- factor(class[ok], levels = class_levels)
    res <- class_anova(x, cl)
    boxplot(x ~ cl, col = colors, las = 2, ylab = "log10(within / between)", main = main, border = COLOR.GREY[["dark"]])
    invisible(res)
  }

  cls_a <- class[match(wb_a$table$gene, gene_ref)]
  cls_b <- class[match(wb_b$table$gene, gene_ref)]

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(7, 4.5, 3, 1))
  res_a <- plot_within_between_by_class(wb_a, cls_a, class_levels, colors, main = lab_a)
  res_b <- plot_within_between_by_class(wb_b, cls_b, class_levels, colors, main = lab_b)

  cat(sprintf("%s: within/between by %s, F(%d,%d) = %.2f, p = %.3g, eta^2 = %.3f\n",
              lab_a, axis_label, res_a$df1, res_a$df2, res_a$f, res_a$p, res_a$eta_sq))
  cat(sprintf("%s: within/between by %s, F(%d,%d) = %.2f, p = %.3g, eta^2 = %.3f\n",
              lab_b, axis_label, res_b$df1, res_b$df2, res_b$f, res_b$p, res_b$eta_sq))
  if (!is.null(res_a$tukey)) { cat(sprintf("%s: Tukey pairwise (%s)\n", lab_a, axis_label)); print(res_a$tukey) }
  if (!is.null(res_b$tukey)) { cat(sprintf("%s: Tukey pairwise (%s)\n", lab_b, axis_label)); print(res_b$tukey) }

  invisible(setNames(list(res_a, res_b), c(lab_a, lab_b)))
}

## map_to_orf() returns the systematic ORF name for each identifier in ids.
## Identifiers that match ORF.PATTERN are returned as given except for the
## suffix separator: the pattern accepts '-' or '.' before the dubious-ORF
## suffix (YAL047C-A / YAL047C.A), since the Jackson file writes it with a
## period, and the result always uses the hyphen that NB.SC and org.Sc.sgd.db
## use, so these genes merge. All other identifiers are
## looked up as common names in org.Sc.sgd.db. Identifiers that do not resolve
## (SUTs, CUTs, snoRNAs, ERCC spike-ins, names absent from SGD) return NA and
## are dropped by the caller. When no identifier at all matches as a common
## name, a message with a sample of them is printed and all are returned as NA.
ORF.PATTERN <- "^Y[A-P][LR][0-9]{3}[CW]([.-][A-Z])?$"
## ============================================================
## 13b. FUNCTIONS APPLIED ACROSS ITEMS BY analysis.R
## ============================================================
## Each function below runs the same code once per element of a vector or list in analysis.R (through
## lapply, sapply, vapply, mapply or apply). It takes the element first and the objects it reads as named
## arguments. Listed roughly in the order analysis.R uses them.

## qc_cell_cutoff: library-size cell filter for one dataset s (a list with lib, the per-cell library
## sizes). Cells below the larger of min_reads and the log10-scale median - k*MAD cutoff are dropped.
## Working in logs makes the rule scale-free, so one setting adapts to any sequencing depth.
qc_cell_cutoff <- function(s, k, min_reads = 500) {
  stopifnot(is.list(s), is.numeric(s$lib), any(s$lib > 0), is.numeric(k), length(k) == 1, k >= 0,
            is.numeric(min_reads), length(min_reads) == 1, min_reads >= 0)
  lib <- s$lib
  lx <- log10(lib[lib > 0])
  lib_cut <- max(min_reads, 10^(median(lx) - k * mad(lx)))
  list(keep = lib >= lib_cut, lib_cut = lib_cut)
}

## fstar_from_r: closed-form split fraction f* = [(2+r) - sqrt(r^2+4)] / (2r) for r = B*Nh/A. r -> 0
## (parent term negligible) recovers the even split f* = 0.5, since there is then no asymmetry to correct for.
fstar_from_r <- function(r) {
  out <- rep(0.5, length(r))
  ok <- is.finite(r) & abs(r) > 1e-8
  out[ok] <- ((2 + r[ok]) - sqrt(r[ok]^2 + 4)) / (2 * r[ok])
  out
}

## gene_pass_group: genes passing the filters in one group g: a finite NB dispersion, a mean count at or
## above the group's depth-scaled floor, and expression in at least n_min cells. Returns a logical vector
## over the rows of fits[[g]].
gene_pass_group <- function(g, fits, floor_mean, n_min) {
  fit <- fits[[g]]
  is.finite(fit$DISP) & fit$DISP < THETA.CAP & fit$MEAN_CT >= floor_mean[[g]] & fit$N_EXPR >= n_min
}

## perm_label_draw: one permutation's relabelings for every mode. total, dpar and inh shuffle pooled
## cell labels, dom pairs random parent cells, and cis and trans are drawn so that the mean-split and
## noise-split versions of a contrast (cis and cis_n, trans and trans_n) are relabeled from the same
## cells. Each version keeps its exact marginal null, so only the coupling between them is new.
## - cis swaps the alleles of a hybrid cell with probability 0.5. One swap flag is drawn per hybrid
##   cell, and the cis and cis_n groups read the flags of their own cells (hyc, hyc_n: positions of the
##   HYC and HYC.N cells among the nHYB hybrid cells), so a cell in both groups swaps in both.
## - trans pools the nP parent cells of one allele with the hybrid cells of the trans group. One random
##   ordering of the parent cells plus every hybrid cell of either trans group (hyt, hyt_n: positions
##   of the HYT and HYT.N cells) serves both versions; each version keeps the members of its own pool
##   in that order. Restricting a uniform ordering to a subset leaves a uniform ordering, so each
##   pooled shuffle is exact, and the parent cells sent to the hybrid-sized group follow the same
##   ordering in both versions. The hybrid positions are in column order of their datasets (ascending),
##   as split_indices_by_depth() returns them. The SC and SE allele pools are ordered separately.
## nSC, nSE are the parent cell counts and b the permutation index (not used).
perm_label_draw <- function(b, hyc, hyc_n, hyt, hyt_n, nHYB, nSC, nSE) {
  stopifnot(is.numeric(nHYB), length(nHYB) == 1, nHYB >= 1,
            all(vapply(list(hyc, hyc_n, hyt, hyt_n), function(x) !is.unsorted(x, strictly = TRUE) && all(x >= 1 & x <= nHYB), logical(1))))
  sc <- perm_trans_pool(nSC, hyt, hyt_n); se <- perm_trans_pool(nSE, hyt, hyt_n)
  flip <- runif(nHYB) < 0.5
  list(
    total     = sample.int(nSC + nSE),
    cis       = flip[hyc],
    cis_n     = flip[hyc_n],
    transSC   = sc$m,
    transSE   = se$m,
    transSC_n = sc$n,
    transSE_n = se$n,
    dom_i   = sample.int(nSC, nHYB, replace = TRUE),
    dom_j   = sample.int(nSE, nHYB, replace = TRUE),
    dparSC  = sample.int(nHYB + nSC),
    dparSE  = sample.int(nHYB + nSE),
    inhSC   = sample.int(nHYB + nSC),
    inhSE   = sample.int(nHYB + nSE))
}

## perm_trans_pool: the pooled relabeling of the trans null for one allele's parent dataset of nP cells.
## Pooled index 1..nP are the parent cells and nP + j the j-th cell of a trans group; the group's first nP
## entries of the returned permutation form the parent-sized group and the rest the hybrid-sized group
## (.fit_split()). One random ordering of the parent cells plus every hybrid cell of either trans group
## (hyt, hyt_n: ascending positions of the HYT and HYT.N cells) serves both versions, each keeping the
## members of its own pool in that order, so each version is an exact uniform pooled shuffle and the
## two share their parental assignment. Returns list(m = HYT pool, n = HYT.N pool). The grid passes
## hyt_n = hyt and reads m. Used by perm_label_draw() and the power grid.
perm_trans_pool <- function(nP, hyt, hyt_n) {
  u    <- sort(union(hyt, hyt_n))
  ord  <- sample.int(nP + length(u))
  id_m <- c(seq_len(nP), nP + match(u, hyt))[ord]
  id_n <- c(seq_len(nP), nP + match(u, hyt_n))[ord]
  list(m = id_m[!is.na(id_m)], n = id_n[!is.na(id_n)])
}

## perm_trans_pool_paired: pairing-preserving version of the trans relabeling for one hybrid group of nHYT
## cells. The Sc side is shuffled as in perm_trans_pool(). The Se side sends the same hybrid cells to its
## parent-sized group, with both alleles of a cell moving together as in the cis swap, and fills the rest of that
## group with random Se parent cells. Returns list(sc, se) permutations for .fit_split(). Each side's
## marginal null is the pooled shuffle exactly when nSC == nSE; for unequal parent counts the number of
## hybrid cells in the parent-sized group follows the Sc side's hypergeometric law.
perm_trans_pool_paired <- function(nSC, nSE, nHYT) {
  stopifnot(nSC >= 1, nSE >= 1, nHYT >= 1)
  ord  <- sample.int(nSC + nHYT)
  a_h  <- ord[seq_len(nSC)]; a_h <- a_h[a_h > nSC] - nSC   # hybrid cells in the parent-sized group
  stopifnot(length(a_h) <= nSE)
  par_a <- sample.int(nSE, nSE - length(a_h))
  list(sc = ord,
       se = c(par_a, nSE + a_h, setdiff(seq_len(nSE), par_a), nSE + setdiff(seq_len(nHYT), a_h)))
}

## fit_ratio_axes: the log2 ratio of two .fit_one() results on each axis, as c(mu, disp, cv2); an axis is NA
## unless both sides are finite and positive. The one definition of a contrast between two groups for the
## gene-level permutation null and the power grid.
fit_ratio_axes <- function(fa, fb) {
  ratio <- function(x, y) if (is.finite(x) && x > 0 && is.finite(y) && y > 0) log2(x) - log2(y) else NA_real_
  c(mu   = ratio(fa[["mu"]], fb[["mu"]]),
    disp = ratio(fa[["disp"]], fb[["disp"]]),
    cv2  = ratio(.cv2_of(fa[["mu"]], fa[["disp"]]), .cv2_of(fb[["mu"]], fb[["disp"]])))
}

## null_cis_axes: one draw of the cis null for hybrid cells with allele counts sc, se and exposure expo. swap
## (logical per cell) exchanges the two alleles within a cell, which keeps the allele pairing; the contrast is the
## swapped Sc allele against the swapped Se allele (fit_ratio_axes()). Shared by permute_contrasts_one()
## and the power grid.
null_cis_axes <- function(sc, se, expo, swap) {
  fit_ratio_axes(.fit_one(ifelse(swap, se, sc), expo), .fit_one(ifelse(swap, sc, se), expo))
}

## null_trans_axes: one draw of the trans null. For each allele the parent cells and the hybrid cells of the
## trans group are pooled and relabeled by perm_sc / perm_se (pooled indices: parents first, then hybrid
## cells; the first length(par) entries form the parent-sized group, .fit_split()). The contrast is the
## parental ratio minus the hybrid allele ratio of the relabeled groups. Shared by permute_contrasts_one()
## and the power grid.
null_trans_axes <- function(par_sc, par_se, hyb_sc, hyb_se, expo_par_sc, expo_par_se, expo_hyb, perm_sc, perm_se) {
  sSC <- .fit_split(c(par_sc, hyb_sc), c(expo_par_sc, expo_hyb), perm_sc, length(par_sc))
  sSE <- .fit_split(c(par_se, hyb_se), c(expo_par_se, expo_hyb), perm_se, length(par_se))
  fit_ratio_axes(sSC$a, sSE$a) - fit_ratio_axes(sSC$b, sSE$b)
}

## eiv_mode_ci_row: errors-in-variables correlation for one mode m with a gene-resampling bootstrap CI
## (B resamples of the rows of contrasts). Draws with negative Vm or Vs give NaN; when more than half of
## the draws are non-finite the survivors are a biased subset that can fall outside [-1, 1], so the CI is
## reported as NA.
eiv_mode_ci_row <- function(m, contrasts, B) {
  e <- eiv_components(contrasts, m)
  n <- nrow(contrasts)
  draws <- with_local_seed(1, replicate(B, eiv_components(contrasts[sample.int(n, n, TRUE), , drop = FALSE], m)["rho_mean_disp"]))
  na_frac <- mean(!is.finite(draws))
  ci <- quantile(draws, probs = c(0.025, 0.975), na.rm = TRUE)
  if (na_frac > 0.5) ci[] <- NA_real_
  b <- list(ci = ci, na_frac = na_frac)
  data.frame(mode = m, n = e["n"],
    rho_raw = round(e["rho_raw_mean_disp"], 3),
    rho_mean_disp = round(e["rho_mean_disp"], 3),
    ms_lo = round(b$ci[1], 3), ms_hi = round(b$ci[2], 3),
    ms_na = round(b$na_frac, 3), row.names = NULL)
}

## se_floor_row: for attenuation floor f, the number of items (pairs or genes, named count_label) whose attenuation attn is at
## least f and the median bootstrap SE among them. se is a named list of SE vectors aligned with attn; each gives
## a median_<name> column.
se_floor_row <- function(f, attn, se, count_label = "n_pairs") {
  keep <- attn >= f
  out <- data.frame(floor = f, n = sum(keep))
  names(out)[2] <- count_label
  for (nm in names(se)) out[[paste0("median_", nm)]] <- if (any(keep)) median(se[[nm]][keep]) else NA_real_
  out
}

## coexpr_seed_check: two-seed adequacy check for the co-expression bootstrap. If SE has not converged at
## this B, SEs underestimated by chance in one run can make many pairs look spuriously significant. cb1,
## cb2: two CB objects at the same B and gene set that differ only in seed, with rows in the same
## gene_i/gene_j order. mode: the contrast whose SE is compared.
coexpr_seed_check <- function(mode, cb1, cb2) {
  mode <- match.arg(mode, c("total", "cis", "trans", "dpar_sc", "dpar_se"))
  seed_compare_core(cb1[[mode]]$se, cb2[[mode]]$se, main = sprintf("%s: bootstrap SE, two seeds", mode), count_label = "n_pairs")
}

## coexpr_perm_draw: one permutation draw for the co-expression null of total and cis: the pooled-cell order
## (n_tot cells) and the per-cell allele swaps for the nH hybrid cells. b is the draw index and is not used.
coexpr_perm_draw <- function(b, nH, nSC, nSE, n_tot) list(
  idx  = sample(n_tot),
  swap = sample(c(TRUE, FALSE), nH, replace = TRUE))

## coexpr_null_draw: one draw of the within-group resamples for the dpar nulls (coexpr_null_dpar_one()):
## cells of the Sc parent, the Se parent and the hybrid, each resampled with replacement within its own
## group at the original size. The hybrid resample is shared by the dpar_sc and dpar_se nulls of the draw,
## as the observed dpar_sc and dpar_se share one hybrid correlation matrix. b is the draw index and is not
## used.
coexpr_null_draw <- function(b, nSC, nSE, nH) list(
  sc = sample.int(nSC, nSC, replace = TRUE),
  se = sample.int(nSE, nSE, replace = TRUE),
  h  = sample.int(nH,  nH,  replace = TRUE))

## axis_mixture_summary: variance explained, effective genes, mixture means and sigmas, and pole sizes of
## one candidate axis a.
axis_mixture_summary <- function(a) c(
  var = a$var_explained, eff_genes = a$eff_genes,
  mu_lo = a$mu[1], mu_hi = a$mu[2],
  sigma_lo = a$sigma[1], sigma_hi = a$sigma[2],
  n_lo = length(a$genes_lo), n_hi = length(a$genes_hi))

## pole_trans_fractions: for one gene group g, its size and the fraction classified any-trans (Trans or
## Cis + Trans) on each burst-kinetics axis.
pole_trans_fractions <- function(g, contrasts, bfreq_class, bsize_class, kbal_class) {
  i <- match(g, contrasts$gene)
  c(n              = length(g),
    pct_bfreq_trans = mean(bfreq_class[i] %in% c("Trans", "Cis + Trans"), na.rm = TRUE),
    pct_bsize_trans = mean(bsize_class[i] %in% c("Trans", "Cis + Trans"), na.rm = TRUE),
    pct_kbal_trans  = mean(kbal_class[i]  %in% c("Trans", "Cis + Trans"), na.rm = TRUE))
}

## allele_cor_boot_se_row: bootstrap SE of row i's allele-residual correlation between sc and se: B
## resamples of the n hybrid cells, with the SD of the resampled correlations returned.
allele_cor_boot_se_row <- function(i, sc, se, n, B) {
  a <- sc[i, ]; b <- se[i, ]
  r <- vapply(seq_len(B), function(k) {
    idx <- sample.int(n, n, replace = TRUE)
    x <- a[idx] - mean(a[idx]); y <- b[idx] - mean(b[idx])
    den <- sqrt(sum(x^2) * sum(y^2))
    if (den > 0) sum(x * y) / den else NA_real_
  }, numeric(1))
  sd(r, na.rm = TRUE)
}


## partial_cor_depth: correlation of gene g's Sc- and Se-allele residuals with cell depth partialled out.
## If depth drove the allele-pair correlation this collapses toward zero relative to the raw correlation;
## otherwise it tracks the raw correlation closely.
partial_cor_depth <- function(g, resid, depth) {
  x <- resid$HYB.SC[g, ]
  y <- resid$HYB.SE[g, ]
  z <- depth
  rxy <- cor(x, y); rxz <- cor(x, z); ryz <- cor(y, z)
  (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
}

## to_seurat_counts: converts a count matrix to a Seurat-ready sparse matrix with dash-delimited feature names.
to_seurat_counts <- function(mat) {
  rownames(mat) <- gsub("_", "-", rownames(mat), fixed = TRUE)
  as(mat, "CsparseMatrix")
}

## stability_task_rows: the candidate resolutions to bootstrap for dataset d: the chosen resolution plus the
## coarsest resolution on the plateau around the highest silhouette. Grid points are grouped into contiguous
## runs (in resolution order) of equal n_clusters, and the run holding the silhouette argmax is the plateau,
## so genuinely different partitions (e.g. 3 vs 5 clusters at nearly equal silhouette) stay on separate
## plateaus, which a silhouette tolerance cannot guarantee. ok is the grid restricted to resolutions that
## passed the size guard.
stability_task_rows <- function(d, sweeps) {
  ok <- sweeps[[d]]$grid[sweeps[[d]]$grid$ok, ]
  rl <- unique(c(sweeps[[d]]$chosen_res, local({
    ok  <- ok[order(ok$res), ]
    grp <- cumsum(c(1, diff(ok$n_clusters) != 0))
    target_grp <- grp[which.max(ok$sil)]
    min(ok$res[grp == target_grp])
  })))
  data.frame(dataset = d, res = rl,
             role = ifelse(rl == sweeps[[d]]$chosen_res, "chosen", "highest silhouette"),
             sil = ok$sil[match(rl, ok$res)], n_clusters = ok$n_clusters[match(rl, ok$res)],
             stringsAsFactors = FALSE)
}

## stability_dataset_inputs: sparse counts and the fit settings (nfeatures, dims_n, metric) the cluster job
## reuses for dataset d.
stability_dataset_inputs <- function(d, counts, dims_n, metric, nfeatures, sweeps) {
  stopifnot(ncol(counts[[d]]) == ncol(sweeps[[d]]$obj))
  list(counts = as(counts[[d]], "CsparseMatrix"), nfeatures = nfeatures[[d]], dims_n = dims_n[[d]], metric = metric)
}

## boot_resample_matrix: draws every bootstrap resample of dataset d up front, one column per replicate,
## from a single seeded stream. Fixing the draws before any Seurat call runs keeps each replicate an
## independent resample, because RunPCA() reseeds the global RNG (seed.use = 42) inside every replicate.
## The same matrix serves every candidate resolution of a dataset, so the resolution comparison is paired
## on identical draws.
boot_resample_matrix <- function(d, counts, B, seed) {
  n <- ncol(counts[[d]])
  with_local_seed(seed, matrix(replicate(B, sample.int(n, n, replace = TRUE)), nrow = n))
}

## metric_clusters: clusters of a Seurat object at one resolution using the annoy distance metric m on the
## first PCs (dims).
metric_clusters <- function(m, obj, dims, resolution) {
  o2 <- FindNeighbors(obj, reduction = "pca", dims = dims, annoy.metric = m, verbose = FALSE)
  o2 <- FindClusters(o2, resolution = resolution, verbose = FALSE)
  Idents(o2)
}

## consistent_pair_row: keeps Sc-hybrid cluster sc as a pair only when its best Se-hybrid partner picks it
## back. sc: cluster id (character); sc_best_se and se_best_sc: best-partner lookups from the hybrid crosstab.
consistent_pair_row <- function(sc, sc_best_se, se_best_sc) {
  sc_id <- as.integer(sc); se_id <- sc_best_se[[sc]]
  if (!is.na(se_id) && identical(se_best_sc[[as.character(se_id)]], sc_id))
    data.frame(sc = sc_id, se = se_id) else NULL
}

## diet_for_markers: strips a Seurat object to the counts and data layers FindMarkers needs, keeping the cluster job's input small.
diet_for_markers <- function(obj) {
  if (packageVersion("Seurat") >= "5.0.0") DietSeurat(obj, layers = c("counts", "data")) else DietSeurat(obj)
}

## seed_check_rows: summary rows of the two-seed bootstrap SE comparison for mode m, read from checks[[m]]:
## correlation, SE ratio median and IQR, and the count (count_label: n_pairs for the co-expression check,
## n_genes for the gene-level check). The co-expression check stores one comparison per mode; the gene-level
## check stores one per quantity within a mode, giving one row per quantity.
seed_check_rows <- function(m, checks, count_label) {
  x <- checks[[m]]
  one <- function(r, ...) {
    out <- data.frame(mode = m, ..., cor = round(r$cor, 3), ratio_median = round(r$ratio_median, 3), ratio_iqr = round(r$ratio_iqr, 3))
    out[[count_label]] <- r[[count_label]]
    out
  }
  if ("cor" %in% names(x)) one(x) else do.call(rbind, unname(Map(one, x, quantity = names(x))))
}

## ============================================================
## 13c. FUNCTIONS THE CLUSTER SCRIPTS RUN
## ============================================================
## The per-item work that the cluster scripts distribute with parLapply or mclapply, and the helpers they share. Each function takes the
## item first and the loaded inputs as named arguments, so the scripts hold only loading, chunking, progress and saving.

## pilot_split_se_one: one gene's bootstrap SE of log(mu) and log(size) for the four datasets needed to compute A and B, plus the bootstrap
## covariance between the two hybrid alleles. HYB.SC and HYB.SE are refit on the SAME resampled cell indices in each replicate because the
## two alleles are measured in the same cells. The RNG is seeded from a hash of the gene's own name, so each gene's draws depend only on its
## name and results are identical for any gene chunking or core count. The hash is a rolling polynomial over the character codes
## (h <- (131 * h + code) mod 2^31 - 1, a prime), so the order of the characters matters and the values spread over the 31-bit seed
## space; it runs in double precision, where every intermediate stays below 2^53 and is exact.
pilot_split_se_one <- function(g, B, expos, mats, seed) {
  stopifnot(is.character(g), length(g) == 1, nzchar(g), is.numeric(seed), length(seed) == 1, is.finite(seed))
  h <- 0
  for (code in utf8ToInt(g)) h <- (131 * h + code) %% 2147483647
  seed_g <- (seed + h) %% 2147483647
  n.p.sc <- length(expos$MIX.SC); n.p.se <- length(expos$MIX.SE); n.h <- length(expos$HYB)
  logmu <- matrix(NA_real_, B, 4, dimnames = list(NULL, c("MIX.SC", "MIX.SE", "HYB.SC", "HYB.SE")))
  logsz <- logmu
  with_local_seed(seed_g, for (b in seq_len(B)) {
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
  })
  allele_cov <- function(m) {
    ok <- is.finite(m[, "HYB.SC"]) & is.finite(m[, "HYB.SE"])
    if (sum(ok) > 2) cov(m[ok, "HYB.SC"], m[ok, "HYB.SE"]) else NA_real_
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
    HYB_logmu_cov     = allele_cov(logmu),
    HYB_logdisp_cov   = allele_cov(logsz),
    row.names = NULL, check.names = FALSE)
}

## fit_group_values: gene g's fitted values in every group of fits (a list of per-group data frames with MU and
## DISP columns): lists named by group holding mu, disp and the squared coefficient of variation cv2.
fit_group_values <- function(fits, g) {
  rows <- vector("list", length(fits)); names(rows) <- names(fits)
  for (grp in names(fits)) rows[[grp]] <- fits[[grp]][g, ]
  mu   <- lapply(rows, `[[`, "MU")
  disp <- lapply(rows, `[[`, "DISP")
  list(mu = mu, disp = disp, cv2 = Map(.cv2_of, mu, disp))
}

## boot_contrasts_one: one gene's bootstrap contrasts across all modes (the unit of work a cluster worker does). draws holds every resample
## index set, drawn before any fit so replicates are independent.
boot_contrasts_one <- function(g, expos, fits, mats, draws) {
  B <- length(draws); modes <- names(.MODES)
  mcol <- paste0("m.", modes); scol <- paste0("s.", modes); ccol <- paste0("c.", modes)
  M <- matrix(NA_real_, B, 3 * length(modes), dimnames = list(NULL, c(mcol, scol, ccol)))
  for (b in seq_len(B)) {
    d <- draws[[b]]
    ## Each dataset is refit on the resample of its draw key (.DRAW.KEY); exposures follow the same key.
    f <- vector("list", length(.DRAW.KEY)); names(f) <- names(.DRAW.KEY)
    for (grp in names(f)) {
      key <- .DRAW.KEY[[grp]]
      f[[grp]] <- .fit_one(mats[[grp]][g, d[[key]]], expos[[key]][d[[key]]])
    }
    mu_all <- lapply(f, `[[`, "mu"); disp_all <- lapply(f, `[[`, "disp")
    gv.mu <- lapply(mu_all, unname)
    gv.bf <- lapply(disp_all, unname)
    gv.cv <- Map(.cv2_of, mu_all, disp_all)
    for (k in seq_along(modes)) {
      grp <- .MODES[[modes[k]]]
      mu <- unlist(gv.mu[grp]); bf <- unlist(gv.bf[grp]); cv <- unlist(gv.cv[grp])
      if (all(is.finite(mu) & mu > 0)) M[b, mcol[k]] <- .contrast_value(modes[k], gv.mu, "MU")
      if (all(is.finite(bf) & bf > 0)) M[b, scol[k]] <- .contrast_value(modes[k], gv.bf, "BFREQ")
      if (all(is.finite(cv) & cv > 0)) M[b, ccol[k]] <- .contrast_value(modes[k], gv.cv, "CV2")
    }
  }
  se   <- apply(M, 2, sd, na.rm = TRUE)
  bdry <- setNames(colMeans(is.na(M[, scol, drop = FALSE])), modes)
  fv <- fit_group_values(fits, g)
  pm <- vapply(modes, .contrast_value, numeric(1), fv$mu,   "MU")
  ps <- vapply(modes, .contrast_value, numeric(1), fv$disp, "BFREQ")
  pc <- vapply(modes, .contrast_value, numeric(1), fv$cv2,  "CV2")
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
    ## cor_<md>: correlation of the mean and bfreq bootstrap draws. Every hybrid dataset of a replicate
    ## is built from the same resampled cell IDs (make_draws()), and the parents are shared, so the
    ## mean contrast (m.<md>) and the bfreq contrast of its source mode (s.<smd>; the noise split for cis
    ## and trans) are measured on the same cells for every mode. NA only when fewer than three draws
    ## give finite values for both.
    x <- M[, paste0("m.", md)]; y <- M[, paste0("s.", smd)]; ok <- is.finite(x) & is.finite(y)
    out[[paste0("cor_", md)]] <- if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
  }
  .r2 <- function(c1, c2) {
	x <- M[, c1]; y <- M[, c2]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) > 2) cor(x[ok], y[ok]) else NA_real_
  }
  out$cor_dpar_mean  <- .r2("m.dpar_sc", "m.dpar_se")
  out$cor_dpar_bfreq <- .r2("s.dpar_sc", "s.dpar_se")
  data.frame(out, row.names = NULL, check.names = FALSE)
}

## permute_contrasts_one: one gene's permutation null: refits the relabeled groups for every mode and permutation, compares the observed
## contrasts (mean, bfreq, CV2, bsize, kbal) with the null by a two-sided permutation p, and adds ploidy-adjusted p-values for the dpar
## contrasts. bsize and kbal nulls recombine the mean and bfreq null draws of the same permutation.
permute_contrasts_one <- function(g, expos, fits, mats, perms, ploidy_shift) {
  l2 <- log2; B <- length(perms); modes <- names(.MODES)
  Nm <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Ns <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  Nc <- matrix(NA_real_, B, length(modes), dimnames = list(NULL, modes))
  fv   <- fit_group_values(fits, g)
  gvmu <- fv$mu; gvbf <- fv$disp; gvcv <- fv$cv2
  obs.m <- vapply(modes, .contrast_value, numeric(1), gvmu, "MU")
  obs.s <- vapply(modes, .contrast_value, numeric(1), gvbf, "BFREQ")
  obs.c <- vapply(modes, .contrast_value, numeric(1), gvcv, "CV2")
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
  ## store: one row of the Nm / Ns / Nc null matrices (mean, bfreq and CV2 axes) for mode md, from the
  ## c(mu, disp, cv2) contrast that fit_ratio_axes() returns
  store <- function(b, md, ax) { Nm[b, md] <<- ax[["mu"]]; Ns[b, md] <<- ax[["disp"]]; Nc[b, md] <<- ax[["cv2"]] }
  for (b in seq_len(B)) {
    p <- perms[[b]]
    sp <- .fit_split(c(cMIXsc, cMIXse), c(expos$MIX.SC, expos$MIX.SE), p$total, nSC)
    store(b, "total", fit_ratio_axes(sp$a, sp$b))

    store(b, "cis", null_cis_axes(cHYCsc, cHYCse, expos$HYC, p$cis))
    ## Noise-split (f_disp) version of the cis null, used only for
    ## the bfreq_cis/cv2_cis output columns, mirroring cis_n in .MODES
    store(b, "cis_n", null_cis_axes(cHYCsc.N, cHYCse.N, expos$HYC.N, p$cis_n))

    store(b, "trans", null_trans_axes(cMIXsc, cMIXse, cHYTsc, cHYTse, expos$MIX.SC, expos$MIX.SE, expos$HYT, p$transSC, p$transSE))
    ## Noise-split version of the trans null, used only for the
    ## bfreq_trans/cv2_trans output columns, mirroring trans_n in .MODES
    store(b, "trans_n", null_trans_axes(cMIXsc, cMIXse, cHYTsc.N, cHYTse.N, expos$MIX.SC, expos$MIX.SE, expos$HYT.N, p$transSC_n, p$transSE_n))

    synth <- cMIXsc[p$dom_i] + cMIXse[p$dom_j]
    es    <- expos$MIX.SC[p$dom_i] + expos$MIX.SE[p$dom_j]
    fs <- .fit_one(synth, es)
    if (is.finite(fs["mu"])   && fs["mu"]   > 0) Nm[b,"dom"] <- l2(fs["mu"])   - mp.m
    if (is.finite(fs["disp"]) && fs["disp"] > 0) Ns[b,"dom"] <- l2(fs["disp"]) - mp.s
    fs.cv <- .cv2_of(fs[["mu"]], fs[["disp"]])
    if (is.finite(fs.cv) && fs.cv > 0) Nc[b,"dom"] <- l2(fs.cv) - mp.c
    d1 <- .fit_split(c(cHYBc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$dparSC, nHYB)
    store(b, "dpar_sc", fit_ratio_axes(d1$a, d1$b))
    d2 <- .fit_split(c(cHYBc, cMIXse), c(expos$HYB, expos$MIX.SE), p$dparSE, nHYB)
    store(b, "dpar_se", fit_ratio_axes(d2$a, d2$b))
    i1 <- .fit_split(c(cHYBsc, cMIXsc), c(expos$HYB, expos$MIX.SC), p$inhSC, nHYB)
    store(b, "inh_sc", fit_ratio_axes(i1$a, i1$b))
    i2 <- .fit_split(c(cHYBse, cMIXse), c(expos$HYB, expos$MIX.SE), p$inhSE, nHYB)
    store(b, "inh_se", fit_ratio_axes(i2$a, i2$b))
  }
  out <- list(gene = g)
  for (md in .OUT_MODES) {
    smd <- .BFREQ_SOURCE[[md]]
    out[[paste0("mean_", md, "_obs")]]  <- unname(obs.m[md])
    out[[paste0("mean_", md, "_p")]]    <- perm_pval(obs.m[md], Nm[, md])
    out[[paste0("bfreq_", md, "_obs")]] <- unname(obs.s[smd])
    out[[paste0("bfreq_", md, "_p")]]   <- perm_pval(obs.s[smd], Ns[, smd])
    out[[paste0("cv2_", md, "_obs")]]   <- unname(obs.c[smd])
    out[[paste0("cv2_", md, "_p")]]     <- perm_pval(obs.c[smd], Nc[, smd])

    ## bsize null: mean null minus bfreq null, permutation by permutation. Where md == smd the two
    ## draws share one relabeling, so their correlation carries through. For cis and trans they come
    ## from the mean-split and noise-split relabelings of the same permutation, which perm_label_draw()
    ## builds from the same cells (shared swap flags for cis, one shared pool ordering for trans).
    bs_null <- Nm[, md] - Ns[, smd]
    obs.bs  <- unname(obs.m[md]) - unname(obs.s[smd])
    out[[paste0("bsize_", md, "_obs")]] <- obs.bs
    out[[paste0("bsize_", md, "_p")]]   <- perm_pval(obs.bs, bs_null)

    ## kbal = bfreq - bsize; its null is the bfreq null minus the bsize null from the same draws, so
    ## kbal gets a permutation p-value for every mode, usable by classify_reg()/classify_dom().
    kbal_null <- Ns[, smd] - bs_null
    obs.kbal  <- unname(obs.s[smd]) - obs.bs
    out[[paste0("kbal_", md, "_obs")]] <- obs.kbal
    out[[paste0("kbal_", md, "_p")]]   <- perm_pval(obs.kbal, kbal_null)

    ## Ploidy-adjusted p-values for the hybrid-versus-parent noise contrasts.
    ## HYB.COMB sums two alleles, which raises bfreq by the shift s (see
    ## ploidy_shift). Removing s from the observed contrast and reading it
    ## against the same zero-centred exchangeability null tests hybrid noise
    ## at the per-genome scale of the parents. _p_ploidy uses this gene's own
    ## s and _p_ind uses s = 1, the fully independent equal-allele bound.
    ## bfreq falls by s, bsize rises by s, kbal falls by 2s, and cv2 rebuilds
    ## with its latent term scaled by 2^s. The mean axis needs no adjustment.
    if (md %in% c("dpar_sc", "dpar_se")) {
      par.nm <- if (md == "dpar_sc") "MIX.SC" else "MIX.SE"
      mh <- gvmu[["HYB.COMB"]]; kh <- gvbf[["HYB.COMB"]]; cvp <- gvcv[[par.nm]]
      nulls <- list(bfreq = Ns[, smd], bsize = bs_null, kbal = kbal_null, cv2 = Nc[, smd])
      s.own <- unname(ploidy_shift[g])
      for (tag in c("ploidy", "ind")) {
        s <- if (tag == "ploidy") s.own else 1
        cvh <- if (is.finite(s) && is.finite(mh) && mh > 0 && is.finite(kh) && kh > 0) 1 / mh + 2^s / kh else NA_real_
        cv2 <- if (is.finite(cvh) && cvh > 0 && is.finite(cvp) && cvp > 0) l2(cvh) - l2(cvp) else NA_real_
        ob <- list(bfreq = unname(obs.s[smd]) - s, bsize = obs.bs + s, kbal = obs.kbal - 2 * s, cv2 = cv2)
        for (q in names(nulls)) out[[paste0(q, "_", md, "_p_", tag)]] <- perm_pval(ob[[q]], nulls[[q]])
      }
    }
  }
  out$ploidy_shift <- unname(ploidy_shift[g])
  data.frame(out, row.names = NULL, check.names = FALSE)
}

## coexpr_bootstrap_one: one bootstrap draw of the co-expression decomposition. draw (from make_coexpr_draws()) resamples the columns of each
## resid dataset. HYB.SC, HYB.SE and HYB.COMB take the same H draw because they are the same cells, and the per-gene ploidy factors are
## fixed, so every draw rescales Rhyb identically. Returns total/cis/trans/dpar_sc/dpar_se at the upper-triangle pair positions.
coexpr_bootstrap_one <- function(draw, resid) {
  d <- coexpr_decompose(list(
    MIX.SC = resid$MIX.SC[, draw$SC], MIX.SE = resid$MIX.SE[, draw$SE],
    HYB.SC = resid$HYB.SC[, draw$H],  HYB.SE = resid$HYB.SE[, draw$H],
    HYB.COMB = resid$HYB.COMB[, draw$H]),
    ploidy_f = attr(resid, "ploidy_f"))
  up <- which(upper.tri(matrix(0, nrow(resid$MIX.SC), nrow(resid$MIX.SC))))
  list(total = d$total[up], cis = d$cis[up], trans = d$trans[up], dpar_sc = d$dpar_sc[up], dpar_se = d$dpar_se[up])
}

## coexpr_part_table: estimate, bootstrap SE (from the running sums), z and two-sided normal p for decomposition nm, one row per gene pair.
coexpr_part_table <- function(nm, point, acc, base) {
  est <- point[[nm]][acc$up]
  s   <- acc[[paste0("sum_", nm)]]; ss <- acc[[paste0("sumsq_", nm)]]
  se  <- sqrt(pmax(0, (ss - s^2 / acc$n) / (acc$n - 1)))      # SD of the draws from the running sums
  data.frame(base, est = est, se = se, z = est / se, p = 2 * pnorm(-abs(est / se)))
}

## top_sq_eigen: the n_keep largest squared eigenvalues of symmetric matrix m by magnitude, descending
## (the rank-matched statistic of the axis tests in Section 4.9).
top_sq_eigen <- function(m, n_keep) {
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  (ev[order(abs(ev), decreasing = TRUE)][seq_len(n_keep)])^2
}

## sym_mat_power: M^power for a symmetric positive definite matrix, from its eigendecomposition (power -0.5
## whitens, 0.5 colors). Whitening needs a full-rank M, so the smallest eigenvalue must exceed tol; a smaller
## one gets a ridge of size `ridge` on the diagonal first, reported on the console, and the check then
## has to pass.
sym_mat_power <- function(M, power, tol = 1e-8, ridge = 1e-6) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M), is.numeric(power), length(power) == 1, tol > 0, ridge > 0)
  e <- eigen(M, symmetric = TRUE)
  if (min(e$values) <= tol) {
    cat(sprintf("sym_mat_power: smallest eigenvalue %.3g, adding a ridge of %.3g\n", min(e$values), ridge))
    e <- eigen(M + ridge * diag(nrow(M)), symmetric = TRUE)
  }
  stopifnot(min(e$values) > 0)
  e$vectors %*% (e$values^power * t(e$vectors))
}

## coexpr_null_dpar_setup: the null of one dpar contrast, built once. The contrast is the leading spectrum of
## Rhyb_f - Rpar, where Rhyb_f is the hybrid's shrunken correlation on the per-genome scale of a haploid
## parent (R_ij f_i f_j, unit diagonal; coexpr_decompose()) and Rpar the parent's shrunken correlation. A
## pooled-label permutation forces both pseudo groups to share one noise level, which neither the
## per-genome rescale nor the groups' different cell counts and correlation structure reproduce. This null
## instead follows Beran and Srivastava (1985): each group keeps its own cells, and the cells are
## recolored so that H0 (equal correlation structure on the per-genome scale) holds exactly in the data
## the recipe sees. The recipe then runs unchanged on within-group resamples.
##   Sigma0   pooled correlation of the two groups on the per-genome scale, weighted by cell count.
##   parent   cells whitened by their own sample correlation and recolored by Sigma0^(1/2), so the recolored
##            cells have sample correlation Sigma0 exactly.
##   hybrid   cells whitened by their own sample correlation and recolored by T_H = C_H^(1/2), where
##            C_H = (Sigma0 - diag(1 - f^2)) / (f f') is the hybrid correlation whose rescaled image is
##            Sigma0 (the rescale is R_ij f_i f_j off the diagonal and 1 on it). When C_H is not
##            positive semidefinite its negative eigenvalues are clipped and the diagonal restored.
## zp, zh: parent and hybrid (HYB.COMB) Pearson residuals, genes x cells in one gene order; ploidy_f: named
## per-gene factors (PLOIDY.F). Returns the recolored cells x genes matrices, f, and an identity check:
## the recipe on the full recolored matrices should give a contrast near zero, and its remainder is the
## difference in shrinkage intensity between the two groups (lambda_p - lambda_h), which scales Sigma0's
## off-diagonal (predicted_rms_remainder). Whitening needs a full-rank sample correlation, so the groups
## need more cells than genes; sym_mat_power() stops on a rank-deficient matrix after a small ridge.
coexpr_null_dpar_setup <- function(zp, zh, ploidy_f) {
  stopifnot(is.matrix(zp), is.matrix(zh), identical(rownames(zp), rownames(zh)), !is.null(rownames(zp)))
  f <- unname(ploidy_f[rownames(zh)])
  stopifnot(length(f) == nrow(zh), all(is.finite(f)), all(f > 0))
  p <- nrow(zh); nP <- ncol(zp); nH <- ncol(zh)
  stdz <- function(z) { s <- scale(t(z)); s[!is.finite(s)] <- 0; s }
  Zp <- stdz(zp); Zh <- stdz(zh)
  ff <- outer(f, f)
  rescale <- function(R) { R <- R * ff; diag(R) <- 1; R }
  Rp <- shrink_cor(t(zp)); Rh <- shrink_cor(t(zh))
  Sigma0 <- (nP * Rp + nH * rescale(Rh)) / (nP + nH)
  ## sample correlations of the standardized cells (unit diagonal also for a constant gene)
  samp_cor <- function(Z) { R <- crossprod(Z) / (nrow(Z) - 1); diag(R) <- 1; R }
  ## hybrid correlation whose rescaled image is Sigma0, projected onto the positive semidefinite cone
  CH <- (Sigma0 - diag(1 - f^2, p)) / ff
  e  <- eigen((CH + t(CH)) / 2, symmetric = TRUE)
  n_clipped <- sum(e$values < 1e-6)
  CH <- e$vectors %*% (pmax(e$values, 1e-6) * t(e$vectors))
  d  <- sqrt(diag(CH)); CH <- CH / outer(d, d)
  Zp_new <- Zp %*% sym_mat_power(samp_cor(Zp), -0.5) %*% sym_mat_power(Sigma0, 0.5)
  Zh_new <- Zh %*% sym_mat_power(samp_cor(Zh), -0.5) %*% sym_mat_power(CH, 0.5)
  colnames(Zp_new) <- colnames(Zh_new) <- rownames(zh)
  ## identity check: recipe on the full recolored matrices
  Rp_img <- shrink_cor(Zp_new); Rh_img <- shrink_cor(Zh_new)
  contrast <- rescale(Rh_img) - Rp_img
  up <- upper.tri(contrast)
  lam_p <- attr(Rp_img, "lambda"); lam_h <- attr(Rh_img, "lambda")
  list(Zp = Zp_new, Zh = Zh_new, f = f, Sigma0 = Sigma0,
       check = c(max_abs_contrast = max(abs(contrast[up])), rms_contrast = sqrt(mean(contrast[up]^2)),
                 rms_sigma0_offdiag = sqrt(mean(Sigma0[up]^2)), lambda_parent = lam_p, lambda_hybrid = lam_h,
                 predicted_rms_remainder = abs(lam_h - lam_p) * sqrt(mean(Sigma0[up]^2)),
                 n_clipped = n_clipped))
}

## coexpr_null_dpar_one: one draw of the dpar null spectrum. setup: coexpr_null_dpar_setup() for the
## contrast. idx_p, idx_h: within-group cell resamples (coexpr_null_draw()). Applies the production recipe to
## the resampled recolored cells (shrunken correlations, hybrid rescaled to the per-genome scale, hybrid minus
## parent) and returns the top n_keep squared eigenvalues of that contrast.
coexpr_null_dpar_one <- function(setup, idx_p, idx_h, n_keep) {
  Rh <- shrink_cor(setup$Zh[idx_h, , drop = FALSE]) * outer(setup$f, setup$f); diag(Rh) <- 1
  top_sq_eigen(Rh - shrink_cor(setup$Zp[idx_p, , drop = FALSE]), n_keep)
}

## coexpr_perm_one: one draw of the five null spectra, for rank-matched testing of every candidate axis. draw:
## list(perm = element of DRAWS.PERM.COEXPR, null = element of DRAWS.NULL.COEXPR). resid: the RESID list.
## nSC, nSE: parent cell counts, used to split each pooled, reshuffled pool back into groups of the original
## sizes. n_keep: ranks retained per decomposition. dpar_setup: list(sc, se) of coexpr_null_dpar_setup()
## results, built once by the cluster script. total and cis are permutation nulls (pooled parents, per-cell
## allele swaps) and trans is their difference. dpar_sc and dpar_se are bootstrap nulls with H0 imposed
## by recoloring (coexpr_null_dpar_setup()), on their own within-group resamples. Returns the top n_keep
## squared eigenvalues by magnitude, descending, for all five decompositions.
coexpr_perm_one <- function(draw, n_keep, nSC, nSE, resid, dpar_setup) {
  pooled   <- cbind(resid$MIX.SC, resid$MIX.SE)
  perm.sc  <- pooled[, draw$perm$idx[seq_len(nSC)]]
  perm.se  <- pooled[, draw$perm$idx[-seq_len(nSC)]]
  perm.hsc <- resid$HYB.SC; perm.hse <- resid$HYB.SE
  perm.hsc[, draw$perm$swap] <- resid$HYB.SE[, draw$perm$swap]
  perm.hse[, draw$perm$swap] <- resid$HYB.SC[, draw$perm$swap]

  d <- coexpr_decompose(list(MIX.SC = perm.sc, MIX.SE = perm.se, HYB.SC = perm.hsc, HYB.SE = perm.hse))
  list(total = top_sq_eigen(d$total, n_keep), cis = top_sq_eigen(d$cis, n_keep), trans = top_sq_eigen(d$trans, n_keep),
       dpar_sc = coexpr_null_dpar_one(dpar_setup$sc, draw$null$sc, draw$null$h, n_keep),
       dpar_se = coexpr_null_dpar_one(dpar_setup$se, draw$null$se, draw$null$h, n_keep))
}

## Scores every chromosome in one species' inputs and returns each gene's
## promoter occupancy track. Rounds of windows run until every window is
## scored or reaches min_core_bp. Minus-strand slices are reversed so the
## promoter-proximal end sits at the end of every vector, matching
## the TATA and poly(dA:dT) scores in score_promoters() on the same gene's sequence. The
## coords the tracks came from travel with the result, so
## score_promoters_nupop() can confirm they match the current promoters.
nupop_occupancy_cluster <- function(inputs, species = 7, model = 4,
                                    window_bp = 250000, flank = 7000, fallback_flank = 2000,
                                    min_core_bp = 5000, cores = 1) {
  stopifnot(is.list(inputs), is.numeric(cores), length(cores) == 1, cores >= 1,
            is.numeric(window_bp), window_bp >= 1, is.numeric(flank), flank >= 0,
            is.numeric(fallback_flank), fallback_flank >= 0, is.numeric(min_core_bp), min_core_bp >= 1)
  ## Runs a table of core windows (seqid, start, end) with `flank` bp of
  ## context, one forked process per window, and returns one result per
  ## row. A window whose process ends early returns NULL or a try-error,
  ## which the caller reads as "split and rerun". NuPoP's compiled code can abort the
  ## process that hosts it, so every window runs in a forked child even with cores = 1
  ## (mclapply() with one core would run in this process): a crash then costs one window,
  ## not the whole job.
  nupop_run_windows <- function(chroms, tasks, flank, species, model, cores, work_root) {
    ## Predicts occupancy for one sequence in the current process and returns
    ## positions core_from to core_to. Each call works in its own temporary
    ## folder, because NuPoP writes its prediction file into the working
    ## directory. The folder name carries the process ID, which is unique
    ## among running processes, so forked workers that start from the same
    ## tempfile() state still receive separate folders. The folders sit under
    ## work_root, which nupop_occupancy_cluster() places beside R's session
    ## temp directory rather than inside it, so each worker's files stay
    ## independent of every other worker.
    one_window <- function(i) {
      s  <- chroms[[tasks$seqid[i]]]
      ws <- max(1, tasks$start[i] - flank)
      we <- min(nchar(s), tasks$end[i] + flank)
      local({
        seg_seq <- substr(s, ws, we)
        core_from <- tasks$start[i] - ws + 1
        core_to <- tasks$end[i] - ws + 1
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
      })
    }
    if (cores >= 2) {
      parallel::mclapply(seq_len(nrow(tasks)), one_window, mc.cores = cores, mc.preschedule = FALSE)
    } else {
      lapply(seq_len(nrow(tasks)), function(i) {
        res <- parallel::mccollect(parallel::mcparallel(one_window(i)))
        if (length(res) == 0) NULL else res[[1]]   # a child that died returns nothing
      })
    }
  }

  chroms <- inputs$chroms
  coords <- inputs$coords
  occ    <- vector("list", length(chroms)); names(occ) <- names(chroms)
  for (sq in names(chroms)) occ[[sq]] <- rep(NA_real_, nchar(chroms[[sq]]))

  ## Node-local scratch space for every window's files, removed at the end
  work_root <- file.path(dirname(tempdir()), sprintf("nupop_work_%d", Sys.getpid()))
  dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_root, recursive = TRUE), add = TRUE)

  tasks <- do.call(rbind, lapply(names(chroms), function(sq) {
    n <- nchar(chroms[[sq]])
    s <- seq(1, n, by = window_bp)
    data.frame(seqid = sq, start = s, end = pmin(s + window_bp - 1, n), stringsAsFactors = FALSE)
  }))

  reduced <- NULL
  failed  <- NULL
  round   <- 0
  while (nrow(tasks) > 0) {
    round <- round + 1
    t0    <- Sys.time()
    res   <- nupop_run_windows(chroms, tasks, flank, species, model, cores, work_root)
    ## A window is scored when its result is numeric and covers every base of the window.
    ok    <- vapply(res, is.numeric, logical(1)) & lengths(res) == tasks$end - tasks$start + 1
    for (i in which(ok)) occ[[tasks$seqid[i]]][tasks$start[i]:tasks$end[i]] <- res[[i]]

    bad   <- tasks[!ok, , drop = FALSE]
    small <- bad[bad$end - bad$start + 1 <= min_core_bp, , drop = FALSE]
    big   <- bad[bad$end - bad$start + 1 >  min_core_bp, , drop = FALSE]

    if (nrow(small) > 0) {
      res2 <- nupop_run_windows(chroms, small, fallback_flank, species, model, cores, work_root)
      ok2  <- vapply(res2, is.numeric, logical(1)) & lengths(res2) == small$end - small$start + 1
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

## Summarizes the regions that needed fallback flanks or stayed unscored,
## so the job log records them next to the progress lines
report_regions <- function(occ, label) {
  red <- attr(occ, "reduced_flank_regions")
  bad <- attr(occ, "failed_regions")
  cat(sprintf("%s: %d promoters, %d region(s) scored with fallback flanks, %d region(s) unscored (%d bp)\n",
              label, length(occ), if (is.null(red)) 0L else nrow(red),
              if (is.null(bad)) 0L else nrow(bad),
              if (is.null(bad)) 0L else sum(bad$end - bad$start + 1)))
  if (!is.null(bad)) print(bad, row.names = FALSE)
}

## bootstrap_ari_job: one bootstrap replicate (job j of the tasks x replicates table). It resamples the cells of the dataset by the job's
## index vector, applies the same preparation as to_seurat_counts() in analysis.R (underscore to dash in gene names, CsparseMatrix),
## reruns Normalize/HVG/Scale/PCA/Neighbors/Clusters at the original settings, and scores the adjusted Rand index against the reference
## clusters on the same resampled cells. A failed replicate returns NA and the job log names it.
bootstrap_ari_job <- function(j, inputs, jobs, tasks) {
  tk  <- tasks[jobs$k[j], ]
  d   <- inputs$data[[tk$dataset]]
  tryCatch({
    idx <- inputs$idx[[tk$dataset]][, jobs$b[j]]
    boot_counts <- d$counts[, idx]
    rownames(boot_counts) <- gsub("_", "-", rownames(d$counts), fixed = TRUE)
    colnames(boot_counts) <- make.unique(colnames(d$counts)[idx])
    boot_counts <- as(boot_counts, "CsparseMatrix")
    boot_obj <- CreateSeuratObject(counts = boot_counts)
    boot_obj <- NormalizeData(boot_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    boot_obj <- FindVariableFeatures(boot_obj, selection.method = "vst", nfeatures = d$nfeatures, verbose = FALSE)
    boot_obj <- ScaleData(boot_obj, features = rownames(boot_obj), verbose = FALSE)
    boot_obj <- RunPCA(boot_obj, features = VariableFeatures(boot_obj), verbose = FALSE)
    boot_obj <- FindNeighbors(boot_obj, reduction = "pca", dims = 1:d$dims_n, annoy.metric = d$metric, verbose = FALSE)
    boot_obj <- FindClusters(boot_obj, resolution = tk$res, verbose = FALSE)
    adjustedRandIndex(as.integer(inputs$ref[[tk$task]][idx]), as.integer(Idents(boot_obj)))
  }, error = function(e) { message(sprintf("%s replicate %d: %s", tk$task, jobs$b[j], conditionMessage(e))); NA_real_ })
}

## dataset_marker_enrichment: for each cluster of dataset d at its final resolution, marker genes against the rest of the SAME dataset's
## cells (FindMarkers with ident.2 left at its default) and GO/KEGG over-representation of the up and down markers (run_enrichment).
## The object carries its own final, validated clustering in Idents, so the question answered is whether the clustering validated for
## THIS dataset corresponds to distinguishable biology; no GSEA is computed. A dataset with fewer than two clusters gives NULL (with a
## message); otherwise the entry is a list of markers/up/down/up_enrich/down_enrich/background/cluster_ids, one element per cluster.
dataset_marker_enrichment <- function(d, ari, inputs, marker_objs, labels, kegg_data, cores) {
  obj <- assemble_cluster_stability(inputs, ari, d, marker_objs[[d]])$final_obj
  cat(sprintf("marker enrichment: %s, %d clusters\n", labels[[d]], length(levels(Idents(obj))))); flush.console()
  cluster_ids <- sort(unique(as.character(Idents(obj))))
  if (length(cluster_ids) < 2) {
    cat(sprintf("%s: only one cluster found; skipping the per-cluster marker/enrichment comparison.\n", labels[[d]]))
    return(NULL)
  }

  markers_list <- mclapply(cluster_ids, function(cc) {
    m <- suppressWarnings(FindMarkers(obj, ident.1 = cc))
    m[order(m$avg_log2FC, decreasing = TRUE), ]
  }, mc.cores = cores, mc.preschedule = FALSE)
  names(markers_list) <- cluster_ids

  ## Strong markers (fold change above 1.25x and adjusted p below 1e-20) split by direction.
  up_list <- down_list <- markers_list
  for (cc in cluster_ids) {
    m <- markers_list[[cc]]
    strong <- abs(m$avg_log2FC) > log2(1.25) & -log10(m$p_val_adj) > 20
    up_list[[cc]]   <- m[strong & m$avg_log2FC > 0, ]
    down_list[[cc]] <- m[strong & m$avg_log2FC < 0, ]
  }

  ## Named fold-change vectors for the background and the up and down sets. avg_log2FC is selected by name:
  ## FindMarkers column order differs across Seurat versions.
  fc <- list(background = markers_list, up = up_list, down = down_list)
  for (kind in names(fc))
    for (cc in cluster_ids) fc[[kind]][[cc]] <- setNames(fc[[kind]][[cc]][["avg_log2FC"]], row.names(fc[[kind]][[cc]]))
  background_list <- fc$background; up_genes_list <- fc$up; down_genes_list <- fc$down

  enrich <- setNames(mclapply(cluster_ids, function(cc) list(
    up   = run_enrichment(names(up_genes_list[[cc]]),   names(background_list[[cc]]), kegg_data = kegg_data),
    down = run_enrichment(names(down_genes_list[[cc]]), names(background_list[[cc]]), kegg_data = kegg_data)), mc.cores = cores, mc.preschedule = FALSE), cluster_ids)
  up_enrich   <- lapply(enrich, `[[`, "up")
  down_enrich <- lapply(enrich, `[[`, "down")

  list(markers = markers_list, up = up_list, down = down_list,
       up_enrich = up_enrich, down_enrich = down_enrich,
       background = background_list, cluster_ids = cluster_ids)
}

## go_enrich_job: one enrichment job k. Over-representation jobs call run_enrichment(); rank-based jobs call gseGO() with a job-specific
## seed so permutation p-values reproduce. A failed job returns NULL and the job log names it.
go_enrich_job <- function(k, go_inputs, jobs, kegg_data) {
  tryCatch({
    job <- jobs[[k]]
    if (job$kind == "ora") {
      run_enrichment(job$genes, job$universe, fdr = go_inputs$GO.FDR, kegg_data = kegg_data)
    } else {
      set.seed(go_inputs$SEED.GO + k)
      suppressWarnings(gseGO(geneList = job$ranks, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = job$ont, nPermSimple = 100000))
    }
  }, error = function(e) { message(sprintf("%s: %s", names(jobs)[k], conditionMessage(e))); NULL })
}

## ora_jobs: one over-representation job per gene set in sets (named group::set), all against the same universe.
ora_jobs <- function(sets, universe, group) {
  jobs <- vector("list", length(sets)); names(jobs) <- paste0(group, "::", names(sets))
  for (i in seq_along(sets)) jobs[[i]] <- list(kind = "ora", genes = sets[[i]], universe = universe)
  jobs
}

## gse_job_set: the rank-based enrichment jobs (BP, MF, CC) for rank list q of ranks_lists, named GSE::GO.GSE.<q>.<ontology>.
gse_job_set <- function(q, ranks_lists) {
  onts <- c("BP", "MF", "CC")
  jobs <- vector("list", length(onts)); names(jobs) <- sprintf("GSE::GO.GSE.%s.%s", q, onts)
  for (i in seq_along(onts)) jobs[[i]] <- list(kind = "gse", ranks = ranks_lists[[q]], ont = onts[i])
  jobs
}
