# Open issues

Source: `docs/function_review.md` (item numbers refer to its section 3), re-checked against the code on `claude/refactor`. Items fixed since the review are not repeated. "Rerun" uses the CLAUDE.md protocol.

## A. Decisions needed (change main results, force cluster reruns)

| Review # | Where | Issue | Rerun if fixed |
|---|---|---|---|
| 4 | `permute_contrasts_one` (13c), `perm_label_draw`, the `PERMS` block in analysis.R | The trans null shuffles the SC and SE allele pools independently, which breaks the within-cell pairing of HYT.SC and HYT.SE that the cis null keeps. This bears on the trans versus cis power asymmetry (project_context.md). The mean-split and noise-split draws of cis and trans are now coupled (shared swap flags, shared pool ordering); only the SC/SE pairing remains. The modes grid runs the independent shuffles and a pairing-preserving swap (`perm_trans_pool_paired`) on the same simulated data; its type I error at ratio 1 is the readout for adopting the swap in `perm_label_draw`. | 2.3 permutation (`gene_perm`) |
| 7 | `power_grid_row` (functions_power.R) | `power_grid_row` simulates one two-group contrast and stays as the total result. The cis (paired alleles) and trans (difference of ratios) contrasts are now simulated by `power_modes_row` with the pipeline's own nulls (Sections 10.9 to 10.14, `power_modes.sub`). Open until that job has run and the MDE ratio, type I and null-SD tables (10.12 to 10.14) are reviewed. | Section 10 modes grid (`power_modes`), additive |
| 8 | `pilot_split_se_one` (13c) | Per-gene seed `seed + sum(utf8ToInt(gene))` yields about 37 distinct seeds for 8000 genes. Marginal SEs unaffected. | Same cascade as the pilot; fix with any pilot rerun |
| 9 | `extract_promoters` | Neighbour bound uses the adjacent gene, not the running maximum end. Rare in yeast. | 6.1 NuPoP |

## B. Local-only changes that alter reported numbers

| Review # | Where | Issue |
|---|---|---|
| 12 | `promoter_direction_test`, `concordance_by_magnitude`, `promoter_noise_candidates` | "Concordant" means the same sign as cis for every quantity, but that prediction was derived for burst frequency. bsize and kbal tails need a decision. |
| 16 | `coexpr_rank_check` | `eigen()` orders by signed value, `coexpr_candidate_axes` and `coexpr_axis_validate` by absolute value. |
| 17 | `kbal_sig`, `plot_burst_kinetics_sig` | Gene-level test uses nominal `p < sig`; the rest of the pipeline uses BH q. |
| 18 | `species_composition_bound` | Label says residual SD, comment says variance fraction; the Wilcoxon p-value is invalid for paired allele cells. |
| 19 | `within_between_decomp` and plots | No NA check on `clusters`; divisor n instead of n - 1; plots drop genes with ratio <= 0. |
| 20 | `map_to_orf` | `YAL047C.A` is passed through unchanged while `NB.SC` uses the hyphen form, so those genes drop out of the Jackson merge. |
| 21 | `add_burst_terms` | Mean x CV2 is a Fano factor only on a count scale; the fluorescence sources depend on instrument units. Also `d$Fano` overwrites the source's own column. |
| 22 | `metabolic_state_cluster` | k-means unseeded; the silhouette-scored fit is discarded and refit. |
| 23 | `run_enrichment` | `pvalueCutoff` stays 0.05 while `qval` feeds only `qvalueCutoff`. |
| 24 | `coexpr_axis_mixtures` | `normalmixEM` unseeded. |

## C. Robustness, no numeric change

- `sig_hist` errors if a contrast is Inf (restrict `range()` to finite values).
- `kegg_local` selects the pathway-name column by position (`KEGGPATHID2NAME[[2]]`), against the select-by-name rule.
- `refine_by_boundary` (unused) draws two independent index vectors for counts and exposures in its bootstrap, pairing one cell's count with another cell's library size.
- `.fit_split`, `fit_split_nb_mm`, `fit_split_nb`: `perm[(n1 + 1):length(perm)]` misbehaves when `n1 == length(perm)`.
- Edge cases in `class_overlap_heatmap`, `se_alpha_col`, `shrink_cor`, `class_anova`, `gene_reliability`, `cor_row`, `to_numeric_matrix`, `read_header_line`.
- `set.seed()` on the global RNG in `class_identity_overlap`, `check_intrinsic_reliability`, `boot_resample_matrix`, `eiv_mode_ci_row` and `pilot_split_se_one`.
- Missing entry validation for most numeric arguments (`fit_source`, `qc_filter_counts`, `mean_adjusted_noise`, `nupop_occupancy_cluster` (13c) with `cores = 1` has no crash isolation).
- `NOISE.CL` leaks if an error occurs between creation and `stopCluster`; several `pdf()` blocks lack `on.exit(dev.off())`.

## D. Project-level open work (CLAUDE.md)

1. `&&` / `||` scalar scan: partial evidence only (reviewers saw no vector operands); no full scan.
2. `pkg_versions()` / `check_pkg_versions()` not finalized.
3. Full pipeline rerun after the R upgrade. Outstanding cluster resubmissions from the earlier estimator, paired-covariance and power-grid changes: `gene_pilot`, `gene_boot` (array of 2), `gene_perm`, `coexpr_boot` (array of 2), `coexpr_perm`, `go_enrich`, `power`, `power_modes` (new, reduced grid first) (and `cluster_stability` only if its inputs change).
4. Trans versus cis power asymmetry (items 4 and 7 above).

## E. Structure and housekeeping

- Section 7.1 to 7.4 now run as list-valued pipelines over `DS.NAMES` (`YSC`, `HVG`, `PCS`, `RES.SWEEP`, `BOOT`) through `prepare_dataset`, `elbow_pcs`, `cluster_dataset` and `umap_dataset`; they were checked with Seurat stubs only, not a real Seurat run. Sections 7.5 to 7.8 (within/between, cell-cycle, metabolic steps) still repeat per dataset, as do the paired Sc/Se blocks, the ploidy block and the power heatmap blocks 10.7 and 10.8; these were left on request.
- `.se_scatter` now serves `plot_contrast_scatter`, `plot_cis_trans_class`, `plot_dom_class` and `plot_coexpr_scatter`. `plot_mean_bfreq_class` keeps its own axes (grey points, per-class regression lines, rectangular panel) and shares only `class_slope_lines`.
- Some `analysis.R` comments quote the algebra of a step in prose that now sits beside inline code; keep them in step with the code when it changes.
- `R/functions.R` holds no one-line function literals except `tryCatch(error = function(e) ...)` handlers; longer anonymous functions inside apply calls (about 34) remain there.
