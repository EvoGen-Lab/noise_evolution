# Open issues after the single-caller and inlining pass

Source: `docs/function_review.md` (item numbers refer to its section 3), re-checked against the code on `claude/refactor`. Items fixed since the review are not repeated. "Rerun" uses the CLAUDE.md protocol.

## A. Decisions needed (change main results, force cluster reruns)

| Review # | Where | Issue | Rerun if fixed |
|---|---|---|---|
| 3 | `boot_contrasts_one`, `add_burst_contrasts`, `eiv_components` | `cor_<mode>` is NA for cis and trans on the premise that the mean and bfreq contrasts share no resampled cells. The trans contrasts share parental draws (`.MODES`). NA enters as r = 0, biasing the trans burst-size and kinetic-balance SEs. | 2.3 bootstrap (`gene_boot`) |
| 4 | `permute_contrasts_one`, `make_perms` | The trans null shuffles the SC and SE sides independently and breaks the within-cell pairing that the cis null keeps. This bears on the trans versus cis power asymmetry (project_context.md). | 2.3 permutation (`gene_perm`) |
| 5 | `coexpr_perm_one` (now in `coexpr_perm.R`) | dpar nulls pool raw hybrid cells while the observed dpar matrices carry the ploidy rescale. Bias direction not quantified. | Section 4 (`coexpr_perm`) |
| 7 | `power_grid_row` (now in `power_grid.R`) | The grid simulates one two-group contrast. It has no paired-allele (cis) or difference-of-ratios (trans) contrast, so it cannot show the asymmetry (open work item 4). | Section 10, additive |
| 8 | `pilot_split_se_one` | Per-gene seed `seed + sum(utf8ToInt(gene))` yields about 37 distinct seeds for 8000 genes. Marginal SEs unaffected. | Same cascade as the pilot; fix with any pilot rerun |
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
- `boot_disp_logse` (unused) draws two independent index vectors for counts and exposures; `refine_by_boundary` inherits it.
- `.fit_split`, `fit_split_nb_mm`, `fit_split_nb`: `perm[(n1 + 1):length(perm)]` misbehaves when `n1 == length(perm)`.
- Edge cases in `class_overlap_heatmap`, `se_alpha_col`, `shrink_cor`, `class_anova`, `gene_reliability`, `cor_row`, `to_numeric_matrix`, `read_header_line`.
- `set.seed()` on the global RNG in `class_identity_overlap`, `make_boot_idx`, `check_intrinsic_reliability`, the `eiv_boot_ci` block and `pilot_split_se_one` (serial use). `make_coexpr_draws` and `make_coexpr_perm_draws` share `SEED.COEXPR`.
- Missing entry validation for most numeric arguments (`fit_source`, `qc_filter_counts`, `mean_adjusted_noise`, `nupop_occupancy_cluster` with `cores = 1` has no crash isolation).
- `NOISE.CL` leaks if an error occurs between creation and `stopCluster`; several `pdf()` blocks lack `on.exit(dev.off())`.

## D. Project-level open work (CLAUDE.md)

1. `&&` / `||` scalar scan: partial evidence only (reviewers saw no vector operands); no full scan.
2. `pkg_versions()` / `check_pkg_versions()` not finalized.
3. Full pipeline rerun after the R upgrade. Outstanding cluster resubmissions from the earlier estimator, paired-covariance and power-grid changes: `gene_pilot`, `gene_boot1/2`, `gene_perm`, `coexpr_boot1/2`, `coexpr_perm`, `go_enrich`, `power` (and `cluster_stability` only if its inputs change).
4. Trans versus cis power asymmetry (items 4 and 7 above).

## E. Structure and housekeeping

- The Section 15 serial wrappers were deleted and the eight cluster-only functions now live in their cluster scripts (`gene_pilot.R`, `gene_boot.R`, `gene_perm.R`, `coexpr_boot.R`, `cluster_stability.R`). `functions.R` Section 15 now holds only the unused calibration helpers, diagnostics and plots, and `fit_offset_nb_mm`.
- Remaining merge candidates: plot bodies shared by `plot_cis_trans`, `plot_mean_bfreq`, `plot_burst_kinetics`; `seed_compare_core` wrappers; `.cohen_kappa` versus the inline kappa; the per-dataset blocks in `analysis.R` around the external-source correlations that could be a loop.
- Single-use helpers nested in other functions (for example `mad_lower`, `.fstar_from_r`) could be inlined into their one caller.
- Hoist repeated constants (theta cap `1e6`, `optimize` interval `c(-4, 15)`, axis count 15).
- Group the shared plotting helpers (`line_colors`, `cluster_cols`, `umap_plot`, `plot_lines`, `legend_page`, `open_grid_pdf`) into their own section.
- Doc comments on some inlined blocks in `analysis.R` still read as the former function signature; reword them to describe the step.
