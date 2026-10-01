# Function review of `R/functions.R`

Reviewed 2026-10-01 on branch `claude/refactor`. No function body was changed.
Eight reviewers each read one slice of `functions.R` against its call sites in
`analysis/analysis.R` and `cluster/scripts/`. The full reports are appendices A1 to A8.

**Line numbers.** The appendices cite line numbers in `functions.R` at commit `8c4f6d8`,
before the comment pass shortened the file. Use the function names to locate code.

**Status after the follow-up pass.** Resolved in the working tree (item numbers are those in section 3):
- Item 1 (MU estimator): the observed fit now takes mu and disp from `.fit_one()`, the same estimator every bootstrap and permutation replicate uses; `glm.nb` supplies only the asymptotic SEs.
- Item 2 (f\*): `pilot_split_se_one()` records the hybrid-allele covariance from a paired cell bootstrap, and `estimate_f_star()` uses the full contrast variance.
- Item 6 (power grid estimators): the observed contrast and the permutation null both use `fit_offset_nb()`. The method-of-moments pair moved to section 15d.
- Item 10 (redundant bsize rows), item 11 (`poly_at_tract` now measures poly(dA) and poly(dT) runs separately), item 14 (`frac_group_sets` and `run_enrichment` drop NA IDs), item 15 (`coexpr_axis_validate`), and the by-name column selection in item 27.
- Section 3.4: all `.sub` files load `r/4.5.2`, and `CLAUDE.md` lists `cluster_stability.R` and `go_enrich.R`.

Still open: items 3 to 5, 7 to 9, 12, 13, 16 to 26, 28 to 31 and the remaining section 3.4 items, including the `gene_split.R` wiring question.

**What I verified myself.**
- The call graph (every function checked for a caller in `analysis.R`, a cluster script, or a reachable function).
- Every call to a library function in `analysis.R`, `functions.R` and the cluster scripts against the function's signature: no mismatches. Calls through `FUN=` or `...` are not covered.
- The headline findings marked **[confirmed]** below, by reading the code and, for R-semantics claims, running a toy example.

Findings marked **[reviewer sim]** rest on a reviewer's simulation and I have not reproduced the numbers.
Everything else is the reviewer's reading of the code, spot-checked.

## 1. What changed in this pass

| Commit | Change | Code changed? |
|---|---|---|
| `6e56959` | 11 more uncalled functions moved into Section 15 (UNUSED), regrouped as 15a/15b/15c | No: all 205 top-level definitions deparse identically |
| `8c4f6d8` | Orphaned banner of a moved function folded into its comment | No |
| `cb28928` | 216 comment edits: purpose-first, positive wording, history narration removed, missing comments added, header outline corrected | No: definitions identical, file parses and sources |

Rerun protocol for all three commits: no sections or subsections need a rerun, and no cluster jobs need resubmission.

Checks run: `parse()` on `functions.R`, `analysis.R` and all 10 cluster scripts; `source('R/functions.R')`; definition-by-definition comparison against the previous commit; call-signature check. R here is 4.3.3 without Seurat, so no analysis code was executed.

## 2. Unused functions (Section 15)

Nineteen functions have no caller in `analysis.R`, in `cluster/scripts/`, or in any function those reach.
- 15a, serial local counterparts of cluster-run functions: `pilot_split_se`, `boot_contrasts`, `permute_contrasts`, `coexpr_bootstrap`, `assemble_coexpr_bootstrap`, `bootstrap_cluster_stability`, `bootstrap_compare_resolutions`.
- 15b, gene-fit calibration: `refine_by_boundary`, `boot_disp_logse`, `chk`, `chk_prec`.
- 15c, interactive diagnostics and plots: `coexpr_raw_cor`, `coexpr_gene_degree`, `plot_coexpr_pair`, `.cohen_kappa`, `report_cluster_marker_enrichment`, `plot_geneset_direction_stack`, `fit_split_nb`, `plot_palette_swatches`.

Functions that look unused from `analysis.R` but are needed by a cluster script, so they stay in the working sections: `.fit_one`, `.fit_split`, `pilot_split_se_one`, `.mode_draw_keys`, `.contrast_value`, `.cv2_of`, `boot_contrasts_one`, `permute_contrasts_one`, `coexpr_bootstrap_one`, `coexpr_perm_one`, `coexpr_acc_init/update/finalize`, `run_enrichment`, `cluster_marker_enrichment`, `nupop_predict_window`, `nupop_run_windows`, `nupop_occupancy_cluster`, `boot_ari_one`, `go_enrich_one`, `fit_offset_nb`, `fit_offset_nb_mm`, `fit_split_nb_mm`, `perm_pval`, `size_log2_ratio`, `power_grid_row`.

Only `boot_disp_logse` is broken if called today (see 3.4); every other Section 15 wrapper matches the current signature of what it wraps. `plot_palette_swatches` needs the colour constants that `analysis.R` defines, so it works only after `analysis.R` has been sourced.

## 3. Discrepancies for your review

Nothing here has been changed. Each item says whether fixing it would change numbers and what it would force you to rerun. Items are grouped by consequence, then roughly by importance.

### 3.1 Would change the main results and force cluster reruns (decisions needed)

| # | Function | Issue | Rerun if fixed |
|---|---|---|---|
| 1 | `neg_binom_fit_offset` vs `.fit_one` | **[confirmed in code; magnitude reviewer sim]** The observed fit takes MU from `glm.nb` (`exp(intercept)`); every bootstrap and permutation replicate takes `sum(y)/sum(exposure)`. With varying exposure the two differ (simulated ratio 0.94 to 1.07, about one SE at n = 300). Observed contrasts are therefore compared with SEs and nulls built from a different estimator. DISP agrees. | 2.3 permutation (`gene_perm`): observed values enter the p-values. `gene_boot` only if `.fit_one` is changed. Local path: observed `_est` columns could be patched in the checkpoint. |
| 2 | `estimate_f_star` | **[reviewer]** Adds the HYB.SC and HYB.SE variances as if independent, but cis resamples the same cells for both alleles. The variance of the difference is overstated, so f* is biased toward 0.5. The pilot stores only SEs, so the covariance cannot be recovered. | 2.1 pilot, then everything downstream of f*: 2.2, 2.3, 4, 10. No local path. |
| 3 | `boot_contrasts_one` | **[reviewer]** `cor_trans` is set to NA on the premise that the mean and burst-frequency trans contrasts share no resampled cells. They share the parental draws (see `.MODES`). The NA becomes r = 0 in `add_burst_contrasts` and `eiv_components`, biasing the trans burst-size and kinetic-balance SEs and the trans EIV row. | 2.3 bootstrap (`gene_boot`). The draw matrix is not saved, so no local path. |
| 4 | `permute_contrasts_one` / `make_perms` | **[reviewer]** The trans null shuffles the SC and SE sides independently, which breaks the within-cell pairing that the cis null keeps. Relevant to the trans-versus-cis power asymmetry in `project_context.md`. | 2.3 permutation (`gene_perm`). |
| 5 | `coexpr_perm_one` | **[reviewer]** The dpar nulls pool raw hybrid cells, but the observed dpar matrices carry the ploidy rescale. The direction of the resulting bias is not quantified. | `coexpr_perm.R` (Section 4). |
| 6 | `power_grid_row` | **[confirmed in code; magnitude reviewer sim]** The observed statistic uses the full MLE, the permutation null uses the method-of-moments estimate. The null SD was 18 to 31% wider in a toy run, so power is underestimated. | `power_grid.R` (Section 10). |
| 7 | `power_grid_row` | **[reviewer]** The grid simulates one contrast between two independent groups. It has no paired-allele (cis) or difference-of-ratios (trans) contrast and no `var_trans = var_cis + 4/N_p` term, so it cannot show the asymmetry that `project_context.md` says any trans-dominant noise claim must address. | `power_grid.R` (Section 10), additive. |
| 8 | `pilot_split_se_one` | **[reviewer sim]** The per-gene seed is `seed + sum(utf8ToInt(gene))`. With 8000 systematic names only 37 distinct seeds came out, so genes sharing a seed share bootstrap streams. Marginal SEs are unaffected. | Same cascade as item 2, so fix together or leave. |
| 9 | `extract_promoters` | **[reviewer]** The neighbour bound uses only the adjacent gene, not the running maximum end, so a gene nested in a longer upstream gene gets a promoter that starts inside it. Rare in yeast. | 6.1 NuPoP: coordinates change, which forces resubmission. |

### 3.2 Local-only fixes that change reported results

| # | Function | Issue |
|---|---|---|
| 10 | `architecture_noise_check` (call sites `analysis.R:1438-1443`) | **[confirmed]** The six `response = "BSIZE"` calls leave `covariate` at its default `"MU"`. The function's own comment says this is the DISP-on-MU test with the sign flipped, so the six bsize rows of `NOISE.VALIDATE` repeat the DISP rows. They need `covariate = "DISP"`. |
| 11 | `poly_at_tract` | **[confirmed]** The pattern `[AT]+` counts mixed runs: `ATATATAT` is a tract of 8, and a toy sequence with a true A-run of 5 returned 8. A poly(dA:dT) tract is a homopolymer run. Changes `polyat_len` and `polyat_delta` and the 6.2 to 6.5 results. If the mixed-run definition is deliberate, the name and comment should say so. |
| 12 | `promoter_direction_test`, `concordance_by_magnitude`, `promoter_noise_candidates` | **[reviewer]** "Concordant" means the same sign as the cis estimate for every quantity. That prediction was derived for burst frequency only. `bsize_cis_est` is mean minus bfreq, so the one-sided test may be on the wrong tail for bsize, and kbal has no derived direction. |
| 13 | `intrinsic_fraction` | **[confirmed]** `intr` is not clipped after the Poisson subtraction. With `intr = -0.1` and `extr = 0` the fraction is exactly 1, and other negative cases fall outside [0, 1] and are silently dropped by `plot_intrinsic_hist`. `ploidy_shift` does clip, so the two treat the same quantity differently. |
| 14 | `frac_group_sets` | **[confirmed]** An NA in `score` puts NA gene IDs into the Low, Average and High sets, which `analysis.R:2144-2145` passes to enrichment. |
| 15 | `coexpr_axis_validate` | **[confirmed]** `paste0("axis", integer(0))` returns `"axis"`, so when no axis validates `extra_axes` becomes `list(<NA> = NULL)`. A validated axis that is not in `extra_axes` gives the same NULL element, and `analysis.R` iterates over those names. Use `intersect()`. |
| 16 | `coexpr_rank_check` | **[reviewer]** `eigen()` orders by signed value, so the rank-1 R², `loading1` and "axis 1" use the largest positive eigenvalue, while `coexpr_candidate_axes` and `coexpr_axis_validate` rank by absolute eigenvalue. |
| 17 | `kbal_sig`, `plot_burst_kinetics_sig` | **[reviewer]** The gene-level test uses nominal `p < sig` while the rest of the pipeline uses BH q-values. Both functions also re-derive columns that `add_burst_contrasts` already stores. |
| 18 | `species_composition_bound` | **[reviewer]** The bound is an SD shift (`|d| x effect`), but `analysis.R` calls it a fraction of residual variance and the plot label says residual SD. For the hybrid Sc-versus-Se allele pair both views are the same cells, so the unpaired Wilcoxon p-value is invalid. |
| 19 | `within_between_decomp` and its plots | **[reviewer]** No NA check on `clusters` or check on `expo`; the within variance uses divisor n rather than n - 1; plots drop genes with ratio <= 0, which conditions the reported rho and ANOVA on positive within-variance. |
| 20 | `map_to_orf` | **[reviewer]** The pattern accepts `YAL047C.A` but passes it through unchanged, while `NB.SC` uses the hyphen form, so those genes silently drop out of the Jackson merge. |
| 21 | `add_burst_terms` | **[reviewer]** Mean x CV² is a Fano factor only on a count scale. For the fluorescence sources (arbitrary units) the `Fano > 1` gate and the burst values depend on instrument units. |
| 22 | `metabolic_state_cluster` | **[reviewer]** k-means has no seed, and the fit that was scored by silhouette is discarded and redone. |
| 23 | `run_enrichment` | **[reviewer]** `pvalueCutoff` stays at 0.05 while `qval` feeds only `qvalueCutoff`, so a q cutoff of 0.2 is first limited by p < 0.05. |
| 24 | `coexpr_axis_mixtures` | **[reviewer]** `normalmixEM` is unseeded, so poles depend on the global RNG state at call time. |

### 3.3 Bugs that change no numbers

| # | Function | Issue |
|---|---|---|
| 25 | `bfreq_bsize_structural` | **[confirmed]** `c(n = e["n"], ...)` keeps the inner names, so the result is named `n.n`, `rho_bfreq_bsize_observed.Cms` and so on. `analysis.R:747` looks up the unsuffixed names and prints NA for all three. Fix with `unname()`. A related reviewer point: `rho_null` divides by the observed Var(bsize), not the null value `Vm + Vf`. |
| 26 | `sig_hist` | **[reviewer]** Errors if a contrast is Inf (`range()` is not restricted to finite values). |
| 27 | `cluster_marker_enrichment` | **[confirmed]** `gene_vec` takes `m[, 2]` by position, against the repo rule on Seurat columns. Also, a failed forked worker returns a try-error that is never checked, so the failure surfaces later as an unrelated error. The function is called by `cluster_stability.R`. |
| 28 | `kegg_local` | **[reviewer]** Selects the pathway-name column by position. |
| 29 | `boot_disp_logse` (unused) | **[reviewer]** Resamples counts and exposures with two independent `sample.int()` calls, pairing one cell's count with another cell's library size. `refine_by_boundary` inherits it. No current result is affected. |
| 30 | `.fit_split`, `fit_split_nb_mm`, `fit_split_nb` | **[reviewer]** `perm[(n1 + 1):length(perm)]` returns `c(n + 1, n)` when `n1 == length(perm)`. No current caller triggers it. |
| 31 | `class_overlap_heatmap`, `se_alpha_col`, `class_identity_overlap`, `shrink_cor`, `class_anova`, `gene_reliability`, `cor_row` and others | **[reviewer]** Edge cases that fail or give NaN on degenerate input (one-level axis, equal SEs, `pe == 1`, p = 1, one-gene input, fewer than 3 pairs). None is reachable with the current data. |

### 3.4 Rule and consistency issues

- **R module version.** All 12 files in `cluster/submit/*.sub` load `r/4.2.2`, but `CLAUDE.md` specifies `r/4.5.2`. **[confirmed]** This matters for open-work item 3 (full rerun after the R upgrade): as written, the cluster jobs run on the old R.
- **`gene_split.R` is orphaned.** **[confirmed]** Nothing in `analysis.R` writes `gene_split_inputs.rda` or loads `gene_split_output.rda`; the split fits run locally. `gene_split.R`, `gene_split.sub` and the cluster-side use of `fit_counts_offset_parallel` are dead unless you want them wired in.
- **Rerun list in `CLAUDE.md`.** `cluster_stability.R` and `go_enrich.R` also source `functions.R` and are not on the list. Changes to `boot_ari_one`, `make_boot_idx`, `assemble_cluster_stability`, `cluster_marker_enrichment` or `go_enrich_one` would need those jobs resubmitted.
- **Section numbering.** `functions.R` section numbers differ from both `analysis.R` and `CLAUDE.md`: its Section 10 is promoter architecture (`analysis.R` Section 6) and its Section 13 is the power analysis (`CLAUDE.md` Section 10).
- **Figure numbers.** The old figure numbers in function comments disagreed with `analysis.R`'s numbering, so I dropped most of them. Eight comment lines still cite figure numbers (for example "Fig 11", "Figs 1 & 2", "Fig 4") and are inconsistent with `analysis.R`; I left them rather than guess which numbering is current.
- **Seeding.** `class_identity_overlap`, `make_boot_idx`, `eiv_boot_ci`, `check_intrinsic_reliability` and `pilot_split_se_one` (serial use) call `set.seed()` on the global RNG. The collapse rule is respected, but later unseeded steps inherit the state. `make_coexpr_draws` and `make_coexpr_perm_draws` are called with the same `SEED.COEXPR`.
- **Validation at entry.** Missing for most numeric arguments. The reviewers list them per function in the appendices. `nupop_occupancy_cluster` with `cores = 1` gives no crash isolation, because `mclapply` falls back to `lapply`.
- **Scalar `&&`/`||` scan (open-work item 1).** Every reviewer checked the conditions in their slice and found no vector-valued operand. That is partial evidence only, not a completed scan.

## 4. Merge, collapse and inline candidates

Functions also called from `cluster/scripts/` cannot move into `analysis.R`, because the cluster jobs source `functions.R` only. They are marked (cluster).

**Duplicate logic that could collapse into one function**
- `class_heatmap` and `plot_class_overlap` are thin wrappers over `class_overlap_heatmap`. `shared_overlap_rng` repeats its log-ratio computation.
- `fit_offset_nb`, `fit_offset_nb_mm` and `.fit_one` are near-duplicates (cluster). `.fit_split`, `fit_split_nb_mm` and the unused `fit_split_nb` are three copies of one body (cluster). `perm_pval` and the inner `pval()` in `permute_contrasts_one` are identical.
- The CV² formula appears three times in Section 3 and 6. `.cv2_of` already exists (cluster).
- `plot_cis_trans`, `plot_mean_bfreq` and `plot_burst_kinetics` share one scatter-with-SE-bars body, as do `plot_coexpr_cis_trans` and `plot_coexpr_dom_class`.
- `coexpr_acc_init`, `coexpr_acc_update` and `coexpr_acc_finalize` repeat five identical blocks (cluster).
- `.cohen_kappa` (unused) duplicates the kappa formula inline in `class_identity_overlap`.
- `fit_counts_offset` and `fit_counts_offset_row` write the same derived columns twice.
- `raw_gene_prefilter` is algebraically `qc_gene_keep` with different arguments. `build_gene_sets` builds nine gene sets but only `$sets$full` is used.
- The three `promoter_noise_candidates` feature blocks and the eligibility test shared by three concordance functions.

**Called once or twice and short enough to inline** (not cluster-called)
- Into their only caller: `mad_lower`, `gene_pass_group`, `eiv_boot_ci`, `species_composition_bound`, `kbal_sig`, `barplot_enrich_pair`, `se_alpha_col`, `.heatmap_legend`, `nupop_window_score`, `apply_cluster_labels`.
- Into `analysis.R`: `check_coexpr_pairs`, `diet_for_markers` (its `else` branch is dead under Seurat v5), `compare_distance_metrics`, `min_detectable_ratio`, `partial_cor_depth`, `bfreq_bsize_structural`, `summarize_go_sets`, `ploidy_coexpr_factor`.
- The `coexpr_seed_compare` and `gene_seed_compare` wrappers over `seed_compare_core`.

Keep as they are: functions with three or more callers and the cluster-called set listed in section 2.

## 5. Other suggestions

1. Fix item 25 and items 10 and 13 to 15 first. They need no cluster rerun; items 10 and 13 to 15 change reported tables or sets, and item 25 changes only a printed line.
2. Decide items 1 to 4 together. They all feed the 2.3 bootstrap and permutation, so one resubmission can carry several of them. Item 2 forces the longest cascade, so decide it before the next full rerun, which open-work item 3 already plans.
3. Update the `.sub` files to `r/4.5.2` when you do that rerun, and add `cluster_stability.R` and `go_enrich.R` to the rerun list.
4. Hoist magic numbers into named constants: the theta cap `1e6` and the `optimize` interval `c(-4, 15)` appear in three places, and the 15 used for axis counts appears in four.
5. Wrap `pdf()` blocks in `on.exit(dev.off())`; several plotting and scoring functions leave a device open if they error partway.
6. Group the shared helpers by purpose: the cluster round-trip helpers (`kegg_local`, `load_cluster_output`, `check_cluster_key`, `go_enrich_one`) and the shared plotting helpers (`line_colors`, `cluster_cols`, `umap_plot`, `plot_lines`, `legend_page`, `open_grid_pdf`) currently sit in Sections 11 and 13.
7. The per-dataset blocks in `analysis.R` around lines 1596 to 1745 repeat the same calls and could be one `for (d in DS.NAMES)` loop.

---



---

# Appendix A1

# Review 1: R/functions.R lines 1-569 (git 8c4f6d8)

Checks run: call-site greps (analysis.R, cluster/scripts), base-R simulations of `neg_binom_fit_offset` vs `.fit_one`, `.fstar_from_r`, `.fit_split`, `split_indices_by_depth`, seed hashing. Scalar `&&`/`||` scan of the range: all operands are length one (L231, L245 assume scalar `init.theta`, L362-369 index single named values). No Seurat calls and no column-by-position selection in this range. Fit functions are RNG-free; `pilot_split_se_one` seeds per gene.

## Discrepancies

| # | Severity | Function, line | Issue and evidence | Direction | Changes numbers? Sections |
|---|---|---|---|---|---|
| 1 | Possible bug | `neg_binom_fit_offset` L254-255 vs `.fit_one` L273 | Observed point estimates use glm.nb `exp(intercept)` for MU (NB-weighted MLE), but every bootstrap/permutation replicate uses `sum(y)/sum(exposure)`. With varying exposure these differ: sim (n=300, lognormal exposure sdlog 0.6, mu=2, theta=1.5, 300 reps) gave mu_glm/mu_fit1 range 0.94-1.07, IQR 0.986-1.020 (about 1 SE at that n). DISP agrees (ratio 0.996-1.002). Observed contrasts (`CONTRAST.FITS`, analysis.R:430; read in `permute_contrasts_one` L1102-1106, `boot_contrasts_one` L762-770) are therefore compared with SEs and nulls from a different estimator. | Use one MU estimator in both (e.g. have the observed fit use the `.fit_one` estimate, or both glm.nb). | Yes. Observed `_est` columns and permutation p-values (2.3 perm job; boot SE unchanged, boot `_est` could be patched locally). Section 4 residuals use MU/DISP from the same fits so they shift too (local). Cluster: gene_perm resubmit (p-values need the observed values inside the job); gene_boot only if `.fit_one` is changed. |
| 2 | Possible bug (statistical) | `estimate_f_star` L410-413 | A = n_h*(se_HYB.SC^2 + se_HYB.SE^2) adds allele variances as if independent. The cis contrast (HYC.SC vs HYC.SE) uses the same resampled cells for both alleles (pilot L355-361; boot `d$HYC` for both), so the variance of the difference is var_sc + var_se - 2cov, and cov > 0 under extrinsic noise. A is overstated, r = B*Nh/A understated, f* biased toward 0.5. PILOT.SE stores only SEs, so the covariance is not recoverable. | Have `pilot_split_se_one` also return the bootstrap SD of log(mu_SC) - log(mu_SE) (and for size) and build A from that. | Yes: f*, then every split-dependent object. Pilot 2.1 rerun, then 2.2, 2.3 (gene_boot, gene_perm), 4 (coexpr), 10. No local path (needs raw bootstrap draws). Decide before the next full rerun. |
| 3 | Possible bug (low impact) | `pilot_split_se_one` L350 | `seed + sum(utf8ToInt(g))` collides for any names with equal character-code sums (anagrams, YAL010C vs YAL001C, etc.). Sim with 8000 Y?L###W names: 37 distinct seeds. Genes sharing a seed share identical index streams, so their bootstrap SEs are correlated; marginal SEs and the pooled median are unaffected. | Positional-weighted hash, e.g. `sum(utf8ToInt(g) * seq_along(utf8ToInt(g)))` or a `digest` hash mod 2^31. | Changes PILOT.SE draws, hence f* slightly and everything downstream (same cascade as #2). Fix together with #2, or leave. |
| 4 | Possible bug (latent) | `.fit_split` L295 | `perm[(n1 + 1):length(perm)]` with n1 == length(perm) evaluates to `c(n+1, n)`; verified returns `b = c(mu=NA, disp=NA)` silently. n1 == 0 also misbehaves. All 9 current callers (L1127-1167) pass n1 < length, so no live effect. | `perm[seq_len(length(perm)) > n1]` or guard `n1`. | No. |
| 5 | Inconsistency | `.fstar_from_r` L392; `estimate_f_star` L418-419 | A non-finite `r` (e.g. NaN/NA from `median` of an empty `ok_*` set, or all-NA pilot) silently yields f* = 0.5 with no warning. Verified `.fstar_from_r(c(NA, Inf)) = 0.5`. | `warning()`/`stop` when `n_genes_* == 0` or r is non-finite. | No. |
| 6 | Inconsistency | `qc_gene_keep` L526 vs analysis.R:256 | Function default `lambda0 = 0.20`; the only caller passes `GENE.LAMBDA0 = 0.10`. `build_gene_sets` and `raw_gene_prefilter` defaults (0.001, 0.10) are repeated as literals at analysis.R:345 and :414. | Make defaults equal the used values or drop defaults; name the constants once in analysis.R. | No. |
| 7 | Note | `build_gene_sets` L559-566 | Only `GS$sets$full` is consumed (analysis.R:415). The other eight sets, `pass`, `n`, `contrasts` are unused. The `cis`/`trans` sets omit the `.N` datasets that the noise-axis cis/trans contrasts also use, so they would be too permissive if someone used them. | Reduce to the full-set filter (see merge section). | No. |
| 8 | Note | `gene_pass_group` L540 | `is.finite(DISP)` drops every gene whose NB fit returned Inf (Poisson-like) or NA in any of the 13 datasets. The filter is marginal (not on a contrast) but conditions on low overdispersion being absent, so the low-noise end is truncated. Worth a sentence in the methods/power discussion. | Document. | No. |
| 9 | Note | `pilot_split_se_one` L362-369 | Inf/NA replicates are dropped before `sd(..., na.rm=TRUE)`, so near-Poisson genes get an SE from surviving draws only (truncation). Pooled medians in `estimate_f_star` blunt this. | Document or report the dropped fraction. | No. |
| 10 | Note | `neg_binom_fit_offset` L245, L247-251 | `init.theta` is never passed by any caller (fit_counts_offset L304, row L450, boot_disp_logse L5188): dead parameter. glm.nb warnings are suppressed and only `converged == FALSE` is trapped, so an alternation-limit fit is accepted silently. | Remove parameter or document; optionally check `fit$th.warn`. | No. |
| 11 | Note | return order | `neg_binom_fit_offset` returns `c(disp, mu, ...)`, `.fit_one` returns `c(mu, disp)`. All callers index by name (verified L362-369, L762, permute `f["mu"]`), so no live bug. | Keep names-only access; consider one order. | No. |
| 12 | Note | `gene_split.R` / `fit_counts_offset_parallel` | analysis.R:386-399 fits the 13 split datasets locally with PSOCK. Nothing saves `gene_split_inputs.rda` or loads `gene_split_output.rda` (grep of analysis/ and cluster/), so cluster script and `.sub` are orphaned and the cluster caller of `fit_counts_offset_parallel` (gene_split.R:41) is dead. CLAUDE.md does not list the split fits as cluster-dependent. | Delete the script or wire it in. | No. |
| 13 | Note | outline / numbering | Outline lists `refine_by_boundary`, `chk`, `chk_prec` under Section 2 and `pilot_split_se` under 1b, but they are defined in Section 15 (L5040-5062 and later). Header says `Functions.R`/`Analysis.R` and omits `power_grid.R`. functions.R "Section 2 GENE FILTER" is used in analysis.R Sections 1 and 2.3, and functions.R section numbers differ from analysis.R numbers. Header and several outline lines fixed in the JSON. | Move those outline lines to 15. | No. |
| 14 | Note (outside chunk) | `boot_disp_logse` L5188 | Calls `sample.int(n, n, TRUE)` twice, so `y` and `exposure` are resampled independently (unpaired). Unused function (15b). | Draw one index vector. | No. |

Validation gaps (repo rule: validate numeric inputs at entry), none are present:
`split_indices_by_depth(n, f)` (f in (0, 0.5], n >= 2), `estimate_f_star(n_h)`, `pilot_split_se_one(B, seed)` (B >= 2), `mad_lower(k)` (empty/all-zero `lib` gives NaN and `qc_cell_keep` then returns NA keep), `qc_cell_keep(k, min_reads)`, `qc_gene_keep(lambda0, cell_frac)`, `build_gene_sets(min_mean, min_expr_frac)` (also: a single-gene input makes `vapply` return a vector so `rownames(pass) <-` fails), `raw_gene_prefilter` (assumes all matrices share row order; combines with positional `&`), `fit_counts_offset_row` / `_parallel` (no `ncol(mat) == length(exposure)` check, which `fit_counts_offset` has), `ckpt_path(n)`.

Magic numbers: theta cap `1e6` (L256, L289, L540), `optimize` interval `c(-4, 15)` (L285), `1e-8` (L392), `min_reads = 500` (L515). Suggest one `THETA.MAX` constant.

## Merge/collapse candidates

| Function | Lines | Call sites | Recommendation |
|---|---|---|---|
| `mad_lower` | 4 | `qc_cell_keep` L516 only | Inline into `qc_cell_keep`. |
| `gene_pass_group` | 4 | `build_gene_sets` L557 only | Inline into the `vapply`. |
| `build_gene_sets` | 20 | analysis.R:414 (uses `$sets$full` only) | Reduce to `full` filter (about 8 lines); drop the eight unused contrast sets. Stays in functions.R (or inline in analysis.R 2.3). |
| `raw_gene_prefilter` | 13 | analysis.R:345 | Near-duplicate of `qc_gene_keep`: `rowSums(m)/sum(m) >= lambda0/min(depth)` is algebraically `mean >= lambda0*depth/min(depth)`, the same rule with the same absolute cell floor. Replace by `qc_gene_keep(PILOT.MATS, lambda0 = 0.001, cell_frac = 0.10)$genes` (or share a `depth_scaled_floor()` helper with `build_gene_sets`). Order differs (sorted); check PILOT.SE row order is not used positionally. Not used by cluster scripts, so inlinable. |
| `fit_counts_offset` vs `fit_counts_offset_row` | 17 / 13 | former: `fit_source` only (L4964, L4989, plus clusterExport string); latter: via `fit_counts_offset_parallel` | Same derived columns (VAR, FANO, CV, BFREQ, BSIZE) written twice. Factor one table-building helper or define `fit_counts_offset` as `do.call(rbind, lapply(seq_len(nrow(mat)), fit_counts_offset_row, ...))`. `fit_counts_offset_row` is repeated `unname(f[...])` seven times. |
| `fit_counts_offset_row` + `fit_counts_offset_parallel` | 13 + 4 | analysis.R:397; gene_split.R:41 (orphaned, see #12) | Keep in functions.R while gene_split.R exists (cluster scripts source functions.R only). If gene_split.R is deleted, the pair is analysis-only and could be one function using `parLapply` with an anonymous closure. |
| `.fstar_from_r` | 6 | `estimate_f_star` L421 (x2) | Keep; closed form is worth a named function. |
| `estimate_f_star` | 16 | analysis.R:354 | Keep (documented pooling logic). |
| `split_indices_by_depth` | 5 | analysis.R:360-361 | Keep. |
| `qc_cell_keep` | 4 | analysis.R:271 | Keep (after absorbing `mad_lower`). |
| `qc_gene_keep` | 10 | analysis.R:284 | Keep; absorbs `raw_gene_prefilter`. |
| `neg_binom_fit_offset` | 36 | `fit_counts_offset` L304, `fit_counts_offset_row` L450 (and unused `boot_disp_logse`) | Keep. |
| `.fit_split` | 4 | `permute_contrasts_one` x9 | Keep. |
| `ckpt_path`, `console_start`, `console_stop` | 3, 4, 4 | analysis.R x10, x10, x1 | Keep. |
| `pilot_split_se_one` | 34 | gene_pilot.R:39 (cluster) and unused `pilot_split_se` | Cannot be inlined into analysis.R. |

Cluster-called functions in this chunk (must stay in functions.R): `pilot_split_se_one`, `fit_counts_offset_parallel`, `fit_counts_offset_row`, `neg_binom_fit_offset`, `.fit_one`, `.fit_split`.

## Other suggestions

- `fit_counts_offset` return: `size[i] <- f["disp"]` carries names; wrap in `unname` for consistency with the row version.
- Section 2 comment block says the filter is "marginal information only"; `gene_pass_group` keys on DISP (an outcome of interest); the comment should say "never on a contrast value" to stay accurate (see #8).
- `pilot_split_se_one` calls `set.seed` on the caller's global RNG when used serially (`pilot_split_se`); acceptable on workers, document or restore `.Random.seed`.
- `console_start` records stdout only (messages and warnings on stderr are not in the transcript); comment now says so.
- `neg_binom_fit_offset` Poisson pre-check uses Pearson/df <= 1, which classifies about half of truly Poisson genes as Inf and some weakly overdispersed genes too; both fit functions share it, so it is consistent, but it makes the DISP = Inf rate depend on N cells.

## Functions reviewed with no issues
`ckpt_path`, `console_start`, `console_stop`, `fit_counts_offset` (formulas VAR = mu + mu^2/size, FANO = 1 + mu/size, CV = sqrt(1/size + 1/mu) verified), `mad_lower`, `qc_cell_keep`, `qc_gene_keep` (apart from #6), `split_indices_by_depth` (f <= 0.5), `fit_counts_offset_row`, `fit_counts_offset_parallel`, `raw_gene_prefilter`. `.fstar_from_r` algebra verified: SE_cis = SE_trans gives r f^2 - (2+r) f + 1 = 0 with root (2+r - sqrt(r^2+4))/(2r); r -> 0 gives 0.5.


---

# Appendix A2

# Review 2: R/functions.R lines 570-1245 (Sections 3-6)

## Discrepancies

| # | Severity | Function, line | Finding |
|---|---|---|---|
| 1 | Likely bug | `bfreq_bsize_structural` 933-935; caller analysis.R:747 | `c(n = e["n"], Vm = Vm, ...)` keeps the inner names, so the result is named `n.n, Vm.Vm, Vf.Vs, Vs_bsize.Vm, rho_bfreq_bsize_observed.Cms, rho_bfreq_bsize_null.Vs, excess_over_null.Cms` (verified by running the function on simulated input). analysis.R:747 reads `STRUCT.BFREQ.BSIZE["rho_bfreq_bsize_observed"]`, `["rho_bfreq_bsize_null"]`, `["excess_over_null"]`; all three lookups return NA, so the printed summary shows NA. Fix direction: `unname()` the elements (or `c(n = unname(e["n"]), ...)`). Numerical results unchanged; only the printed line and the saved named vector (analysis.R:768) change. No cluster section affected. |
| 2 | Possible bug | `boot_contrasts_one` 765-768 and `.MODES` 645/647; consumers `add_burst_contrasts` 829, `eiv_components` 864 | `cor_trans` is set to NA on the premise that the mean (`trans`) and bfreq (`trans_n`) draws "share no resampled cells". `.MODES` shows `trans` uses draws MIX.SC, MIX.SE, HYT and `trans_n` uses MIX.SC, MIX.SE, HYT.N. The two contrasts share the parental MIX.SC/MIX.SE resamples, so the draws are correlated. The NA is converted to r = 0, which biases `bsize_trans_se`, `kbal_trans_se` and the trans row of `eiv_components` (Cms correction). For `cis` the draws are independent (HYC vs HYC.N), but HYC and HYC.N are depth-split subsets of the same hybrid cells (analysis.R:360-361), so the real sampling correlation is not zero either; independent resampling of the two index vectors makes the bootstrap blind to it. Fix direction: compute `cor(M[,m.md], M[,s.smd])` for every mode (it is ~0 where draws are independent) instead of NA. Would change trans (and slightly cis) bsize/kbal SEs and EIV trans estimates. Needs a rerun of the 2.3 bootstrap job (gene_boot.R); M is not saved, so no local path. The permutation side (item 3) is related. |
| 3 | Note | `permute_contrasts_one` 1140-1141, `make_perms` 1084-1087 | Trans null permutes SC and SE sides with independent vectors (`transSC`, `transSE`), so the within-cell pairing of HYT.SC/HYT.SE (same hybrid cells) is not preserved in the null, while cis swaps alleles within cells and keeps it. The null variance of the hybrid allele ratio can therefore differ from the observed one (likely conservative for mean). Relevant to the trans-vs-cis power asymmetry in project_context.md. Changing it alters `_p` for trans/trans_n and every downstream `_q`; needs a 2.3 permutation rerun. |
| 4 | Possible bug | `bfreq_bsize_structural` 928-932 | The null correlation divides by `sqrt(Vf * Vs_bsize)` with `Vs_bsize` the observed `Vm + Vf - 2*Cmf`. Under the stated null (Cov(mean,bfreq) = 0) Var(bsize) is `Vm + Vf`. Observed and null use different denominators, so `excess_over_null` mixes two scales. Using the null variance for `rho_null` makes it internally consistent. Changes the reported rho_null and excess only (analysis.R:747); no cluster rerun. |
| 5 | Possible bug | `sig_hist` 1023-1026 | `range(x, na.rm=TRUE)` with an Inf/-Inf contrast (reachable: `.contrast_value` returns +/-Inf for mu = 0, and `pm` in `boot_contrasts_one` is not gated) makes `seq(lo, hi, by=brk)` error. Callers analysis.R:567-575 pass `*_est` columns. Use `is.finite(x)` in the range. No numeric change otherwise. |
| 6 | Inconsistency | `boot_contrasts_one` 746 vs `.cv2_of` 701 | The CV2 point estimate re-implements `.cv2_of` inline (also inline at 1103 and as `cv2_of` at 1122 with `f["mu"]`). Three copies of one formula; keep one helper. No numeric change. |
| 7 | Inconsistency | `plot_mean_bfreq` 965-968 vs code | Comment says "Axes scaled independently" and "mean vs noise (size)"; code uses one shared symmetric range and the y axis is the bfreq (NB dispersion) contrast, labelled "dispersion (log2)". Also comment lists three modes but eight are accepted. Comment fixed in the JSON; consider the ylab "burst frequency (log2)" for consistency with the other plots. |
| 8 | Note | `eiv_boot_ci` 880 | `set.seed(seed)` is applied to the global RNG on every call (default seed 1), so each mode sees the same gene resample (paired across modes, fine) but the caller's RNG stream is altered. Safe here because no Seurat/RunPCA call follows inside the loop; note only. |
| 9 | Note | `eiv_boot_ci` 883, 892 | `na_frac` counts `!is.finite`, but `quantile(na.rm=TRUE)` drops only NA/NaN; an Inf draw (Vm*Vs = 0) would enter the quantiles. Unlikely; use `draws[is.infinite(draws)] <- NA`. |
| 10 | Note | `eiv_components` 867-876 | `cor(X, Y)`/`var` fail or warn when fewer than 2 genes pass `ok` or when a column is constant (e.g. a bootstrap resample of a tiny set); no minimum-n guard. Mode `dom` etc. uses `match.arg` but `eiv_table` passes `m` unvalidated apart from that. |
| 11 | Note | `permute_contrasts_one` 1234 | `ploidy_shift[g]` indexes by gene name; if `g` is absent or PLOIDY.SHIFT is unnamed it silently yields NA (then `_p_ploidy` NA). analysis.R:431 builds `PLOIDY.SHIFT <- ploidy_shift(...)[GENES]`, so names are present; add `stopifnot(g %in% names(ploidy_shift))` at entry. |
| 12 | Note | `permute_contrasts_one` 1138, 1150 | `Nm[,"cis_n"]`, `Nm[,"trans_n"]` are computed but never read (outputs read `Nm[,md]` for the mean and `Ns/Nc[,smd]` for bfreq/CV2); likewise `Ns/Nc[,"cis"]`, `Ns/Nc[,"trans"]`. Harmless: one fit gives both. |
| 13 | Note | `make_draws` 573-583, `make_perms` 1074-1094 | No entry validation of `B`/`NPERM`, `ncells` names, or positivity (repo rule: validate numeric inputs). Both use `set.seed(seed)`; analysis.R:438-439 uses SEED.BOOT and SEED.PERM; if these are ever equal, boot and perm index streams coincide. All index draws are made before any fit (rule satisfied). |
| 14 | Note | `fdr_columns` 1069 | Pattern `_p(_ploidy|_ind)?$` matches any column ending `_p`; only permutation columns currently do. OK; `p.adjust` leaves NA unadjusted and counts only non-NA tests (comment correct). |
| 15 | Note | `add_burst_contrasts` 829 | `r[!is.finite(r)] <- 0` also converts NA from "fewer than 3 finite draws" (genes with failed fits) to 0, not only the structural cis/trans NA. Minor SE bias for sparse genes. |

Scalar `&&`/`||`: all operands checked (lines 703, 746, 1117-1124, 1157-1158, 1171, 1225, 1229-1230) are length one (named length-1 subsets of fit vectors). No Seurat column selection by position in this chunk. Algebra of `add_burst_contrasts` (bsize SE, kbal SE = sqrt(4ss^2 + sm^2 - 4 r sm ss), covFB) and the ploidy adjustment (bfreq -s, bsize +s, kbal -2s, cv2 latent 2^s) was re-derived and is correct. Bootstrap and permutation functions are pure (no RNG inside), so fork/parLapply/mclapply use is safe.

## Merge/collapse candidates

| Function | Calls | Lines | Recommendation |
|---|---|---|---|
| `.mode_draw_keys` 677 | boot_contrasts_one:765 (twice) | 1 | Keep, or inline as one `identical()` test on `.MODES`/`.DRAW.KEY`; used only in cluster-called function, so it must stay in functions.R. If item 2 is adopted (cor always from M), `.mode_draw_keys` and `.DRAW.KEY` become unnecessary (delete both). |
| `.cv2_of` 701 | boot_contrasts_one:729 only | 4 | Reuse at lines 746, 1103 and replace `cv2_of` (1122); then 3 duplicate definitions collapse to 1. Must stay in functions.R (cluster). |
| `.OUT_MODES`/`.BFREQ_SOURCE` | multiple | - | Keep. |
| `eiv_boot_ci` 879 | eiv_table:899 only | 17 | Fold into `eiv_table` (or keep as internal); not cluster-called. The `stats` argument and the 1-stat matrix special case (882) are unused generality: eiv_table uses the default only. |
| `eiv_table` 897 | analysis.R:471 (once) | 10 | Keep; it is the only analysis entry. Note the return value is auto-printed only if top level (analysis.R:471 is top level). |
| `bfreq_bsize_structural` 925 | analysis.R:738 (once) | 12 | Could be inlined at analysis.R 3.6 after `eiv_components`; keep if reused. Not cluster-called. |
| `plot_cis_trans`, `plot_mean_bfreq`, `plot_burst_kinetics` | analysis.R:479-481 (3 calls each, one panel) | 17/16/14 | The three share the same scatter+SE-bar body (ok filter, abline, segments, points) and also `*_class` and `*_sig` versions elsewhere (plot_cis_trans_class 103, plot_mean_bfreq_class 104, plot_burst_kinetics_sig). Merge into one `plot_se_scatter(x, sx, y, sy, xlab, ylab, main, sym)` helper; each wrapper then becomes 3 lines. Not cluster-called. |
| `.sym` 942 | plot_burst_kinetics:1005, plus Section 15 line 2891 | 1 | `plot_cis_trans` (954) and `plot_mean_bfreq` (976) recompute the same limit inline; use one helper. |
| `sig_hist` 1013 | analysis.R:567-575 (9 calls) | 18 | Keep. |
| `fdr_columns` 1068 | analysis.R:470 (once) | 5 | Short; keep (documents the q-value family rule). |
| `make_draws`/`make_perms` | analysis.R:438,448 / 439 | 11/21 | Keep (outputs saved for cluster). make_perms and make_draws could share a seed/validation helper. |
| `permute_contrasts_one` | gene_perm.R (mclapply), `permute_contrasts` (Section 15) | 148 | Keep (cluster). Long; the nine near-identical `Nm/Ns/Nc` assignment lines (1128-1168) could be one helper `set_row(b, mode, fa, fb)` that writes the three rat/ratcv2 calls. |
| `boot_contrasts_one` | gene_boot.R:50 | 72 | Keep (cluster). |

Functions whose only other callers are Section 15: `permute_contrasts_one`/`boot_contrasts_one` also feed `permute_contrasts`/`boot_contrasts` (5067, 5074); the cluster scripts remain the real callers.

## Other suggestions

- Magic numbers: `B = 2000` in `eiv_boot_ci`/`eiv_table`, `0.5` NA-fraction CI cutoff (893), `0.33` bar alpha repeated in three plot defaults, `by = brk` 0.1; hoist to named constants or arguments.
- `eiv_table` has no `B`/`seed` pass-through for `seed`; add for symmetry.
- Mode vector `c("total","cis","trans","dom","dpar_sc","dpar_se","inh_sc","inh_se")` is typed in four signatures (854, 897, 969, 997); use `.OUT_MODES` as the default, with `match.arg(mode, .OUT_MODES)`.
- Line 771 is indented with a tab inside `.r2`; normalise.
- `plot_*` functions call `max(abs(...), na.rm=TRUE)` on possibly empty vectors (all NA); returns -Inf with a warning and then `plot` fails. Guard with `if (!any(ok)) return(invisible(NULL))`.
- `fdr_columns` could print the number of families adjusted as a sanity check.
- Parameter `ploidy_shift = NULL` at 1096 is optional but gene_perm.R:36-41 branches on `formals()`; since the argument now exists, that compatibility shim (HAS.PLOIDY.ARG) could be removed.

## Functions reviewed with no issues

make_draws (logic), .contrast_value, .cv2_of, add_burst_contrasts (formulas), eiv_components (formulas), eiv_table, plot_cis_trans, plot_burst_kinetics, fdr_columns, make_perms (logic), permute_contrasts_one (p-value, bsize/kbal nulls, ploidy adjustment logic).


---

# Appendix A3

# Review 3: R/functions.R lines 1246-1976 (Section 7 co-expression)

## Discrepancies

| # | Severity | Function (file:line) | Issue / evidence | Direction | Numerical impact / sections |
|---|---|---|---|---|---|
| 1 | Likely bug | `coexpr_axis_validate` (1954-1974) | `extra_axes[paste0("axis", validated_axes)]`: `paste0("axis", integer(0))` returns `"axis"` (verified in R), so when no axis validates the result is `list(<NA> = NULL)` instead of an empty list. Also a validated axis that is not in `extra_axes` (extra_axes only holds `sig_axes`, i.e. axes passing the 1% variance and participation-ratio floors) yields a `NULL` element named NA. Caller analysis.R:977-990 / 1000+ then iterates `names(EXTRA.AXES)` and reads `$enrich_lo` from NULL. | Use `intersect(paste0("axis", validated_axes), names(extra_axes))` (or `if (length(validated_axes))`). | No change to the p/q table. Local only. |
| 2 | Likely bug | `frac_group_sets` (1690-1695) | `genes[score <= cuts[1]]` with NA in `score` returns `NA` elements (logical-NA indexing), so Low/Average/High gain NA gene IDs. Callers analysis.R:2144-2145 pass `intrinsic_frac`, which can be NA. NA IDs then go to enrichGO. | Use `which(...)` or `!is.na(score) &`. | Only if NAs exist; changes INTR.SETS and downstream GO. Not cluster-dependent. |
| 3 | Possible bug | `coexpr_rank_check` (1872-1888) | `eigen(symmetric=TRUE)` orders by signed value (verified: eigenvalues of a zero-diagonal difference matrix come back as `3, -3`). The reconstruction, `loading1` and "axis 1" (analysis.R:1057 `RANK.CHECK$loading1`, 1079 `vectors[,2]`) use the largest *positive* eigenvalues, while `coexpr_candidate_axes`/`coexpr_axis_validate` rank by `abs(eigenvalue)`. If the dominant axis is negative (trace is 0, so negative eigenvalues always exist), the rank-1 R^2 and "axis 1" do not describe the top-|lambda| axis. `R^2 = cor^2` is also scale-free, so a rank-k fit with a wrong-magnitude lambda still reports a high R^2. | Order components by `abs(values)` for the reconstruction and `loading1`; or document. Consider R^2 = 1 - SSE/SST. | Could change RANK.CHECK R^2 / AXIS*.CT (local only). Section 4 downstream of cluster output, no resubmission. |
| 4 | Possible bug | `coexpr_axis_mixtures` (1919-1945) | `normalmixEM` uses random starts and the function never seeds, so poles and `sig_axes` enrichment depend on the global RNG state at call time (analysis.R:948). `pdf()` is opened and `dev.off()` is not guarded with `on.exit`; an error in `enrichGO`/mixture leaves the device open. | `set.seed()` (once, before the loop, no Seurat call inside) and `on.exit(dev.off())`. | Mixture membership could change run to run; local only. |
| 5 | Inconsistency | `coexpr_perm_one` (1410-1445) vs `coexpr_decompose` | The dpar nulls pool raw (un-rescaled) hybrid cells, while observed `dpar_sc/dpar_se` (COEXPR.POINT) carry the `ploidy_f` rescale (R_ij f_i f_j). The rank-matched test therefore compares ploidy-adjusted observed spectra with a null built on unadjusted hybrid structure (the code comment states this). The adjusted Rhyb has smaller magnitudes if f < 1, so the test is conservative/anti-conservative in a direction not quantified. | State direction in analysis.R, or apply `ploidy_f` to the pseudo-hybrid group in the null (needs `resid` attribute and a cluster rerun). | Fix would change dpar null spectra: resubmit `coexpr_perm.R` (Section 4); 4.9 dpar validation changes. |
| 6 | Inconsistency | `coexpr_decompose` comment (old 1300) | Claimed "The scaling commutes with shrinkage". Not exact: lambda is data-driven, so shrinking `D R D` would give a different lambda than shrinking `R`. The code shrinks first then rescales, which is a design choice, not an identity. Comment rewritten in the JSON to avoid the claim. | Documentation only. | None. |
| 7 | Possible bug | `shrink_cor` (1268-1283) | `for (i in 1:(p - 1))` with `p == 1` iterates `1:0` and indexes out of range; `n < 3` gives division by zero. All current callers have p in the hundreds, so latent. The double loop is O(p^2) R-level iterations (~245k at p=700) executed 4-5 times per bootstrap draw and 4 times per permutation draw; this is the stated hot spot (coexpr_boot.R header). | Add `stopifnot(p >= 2, n >= 3)`. A vectorized form gives the identical lambda (checked numerically, 1.00906 vs 1.00906): `S2 <- crossprod(Zs^2); V <- n/(n-1)^3 * (S2 - (n-1)^2 * R^2 / n); lam <- sum(V[ut]) / sum(R[ut]^2)`. | Vectorization is numerically identical up to rounding; it would speed Sections 4.2/4.6 cluster jobs only if resubmitted (not required). |
| 8 | Possible bug | `class_anova` (1706-1716) | `aov` errors (`contrasts can be applied only to factors with 2 or more levels`, verified) when fewer than 2 levels remain after NA removal. Row lookups `ov["class", ...]` rely on partial matching of padded row names ("class      "); works but fragile. 16 call sites in analysis.R:1308-1325. | Guard `nlevels(df$class) < 2` returning NA list. | None unless a class is empty. |
| 9 | Possible bug | `gene_reliability` (1520-1526) | If `fits` has several datasets but `genes` has length 1, `sapply` returns a vector and `apply(rho_mat, 1, ...)` errors. A row with all NA gives `Inf` with a warning. Callers (analysis.R:793, 845, 1186-1187) are safe. | `vapply` + `do.call(cbind, ...)`; return NA when all NA. | None. |
| 10 | Possible bug | `coexpr_candidate_axes` (1897-1911) | `table = ...[1:15, ]` hard-codes 15; if `n_candidate < 15` the table contains NA rows. `n_candidate <- min(n_candidate, length(values) - 1)` has an unexplained `- 1`. `%d` in `sprintf` is fed `eff_genes_min` (double); errors if non-integer. 15 also appears in `n_top` (coexpr_axis_validate), `n_keep` (coexpr_perm_one) and `N.KEEP` (analysis.R:886). `coexpr_axis_validate` indexes `axis_var[1:n_top]` with no length check (NA if fewer than n_top axes). | One shared constant; `stopifnot(n_top <= length(axis_var))`; `%g` in sprintf. | None. |
| 11 | Possible bug | `check_coexpr_reliability`, `stratified_rho_sample`, `check_intrinsic_reliability` | `cut(x, quantile(x, ...))` fails with "breaks are not unique" when quantiles tie (heavy ties in attenuation, e.g. many rho = 0 for DISP = Inf). Also these functions plot as a side effect (needs an open device; they are called inside `pdf()` in analysis.R:845-848) and mix compute with plotting. | `unique(quantile(...))`. | None. |
| 12 | Note | `coexpr_bootstrap_one` (1362-1371) | Recomputes `which(upper.tri(matrix(0, p, p)))` every draw (p^2 allocation) although `coexpr_acc_init` already builds it. `resid$HYB.COMB` NULL gives an obscure error (the function requires HYB.COMB). The estimate in `coexpr_acc_finalize` is `pt`, SE from bootstrap SD about the bootstrap mean. | Take `up` as an argument or check `HYB.COMB`. | None. |
| 13 | Note | `coexpr_acc_finalize` (1493-1508) | `se = 0` gives `z = Inf` / `NaN` (est = 0), `p` NaN; `classify_reg` turns NaN p into Ambiguous only via the `cls == ""` rule. Variance via `ss - s^2/n` (cancellation negligible at n = 1000 but Welford is safer). `resid` only supplies gene names and p. | Treat se = 0 explicitly. | Rare. |
| 14 | Note | `make_coexpr_draws` / `make_coexpr_perm_draws` | Both use `set.seed(seed)`; analysis.R:834 and 898 call them with the same `SEED.COEXPR`, so the bootstrap index stream and permutation stream start from the same RNG state (not harmful, but correlated). No Seurat call occurs between seed and draws, so the collapse rule is respected. | Offset the permutation seed (`SEED.COEXPR + 2`). | Would change permutation draws: resubmit `coexpr_perm.R`. |
| 15 | Note | `coexpr_seed_compare` (1843-1846) | Assumes cb1/cb2 rows align and does not check; a mismatched order would give a meaningless correlation. | `stopifnot(identical(cb1[[mode]]$gene_i, cb2[[mode]]$gene_i), ...)`. | None. |
| 16 | Note | `gene_seed_compare` (1857-1862) | `match.arg` mode list includes modes whose column (e.g. `mean_dom_se`) may be absent; `bc1[[col]]` is then NULL and `cor` errors with an unhelpful message. | `stopifnot(col %in% names(bc1))`. | None. |
| 17 | Note | `partial_cor_depth` (1670-1673) | No NA handling and no guard for |r| = 1 (division by zero). | `use = "complete.obs"`. | None. |
| 18 | Note | `nb_residuals` (1257-1263) | No check that `rownames(mat)` equals `rownames(fit)`; mismatch silently pairs wrong MU/DISP. Callers subset both by the same gene vector. | `stopifnot(identical(rownames(mat), rownames(fit)))`. | None. |
| 19 | Note | `check_coexpr_pairs` / `coexpr_seed_compare` comments | The `coexpr_seed_compare` description was glued onto `check_coexpr_pairs` and `coexpr_seed_compare` had none. Fixed in the JSON. | Applied via comments. | None. |
| 20 | Note | `check_intrinsic_reliability` | `rho_true = rho_obs / a` is not clipped to [-1, 1]; at low reliability values exceed 1. `set.seed(seed)` is called twice (sample, then bootstrap) and resets the caller's global RNG. | Clip or flag; use `withr`-style local seed. | None. |

## Merge/collapse candidates

| Function | Lines | Call sites | Recommendation |
|---|---|---|---|
| `plot_coexpr_cis_trans` + `plot_coexpr_dom_class` | 14 + 14 | analysis.R:911 and 878 (one each) | Near-duplicate bodies (scatter with y = +/-x guides, legend). Merge into one `plot_coexpr_scatter(x, y, cls, levels, cols, xlab, ylab, main, lim)` with the two wrappers deleted, or one function with a `type` argument. Analysis-only, so safe. |
| `coexpr_seed_compare` + `gene_seed_compare` + `seed_compare_core` | 4 + 6 + 13 | analysis.R:861 (coexpr), analysis.R:490 (gene, multiple) | Wrappers are 4-6 lines. Could be one `seed_compare(se1, se2, main, count_label)` called directly from analysis.R with `cb1[[m]]$se` / `bc1[[col]]` / matched `bc2[[col]][m]`; or keep `seed_compare_core` and inline the two wrappers. Analysis-only. |
| `check_coexpr_pairs` | 4 | analysis.R:857 (once; CB loaded without check at 846) | Inline as a `stopifnot(nrow(CB2$total) == EXPECTED.PAIRS)` in analysis.R Section 4.4, or apply to both CB objects in a loop. |
| `partial_cor_depth` | 4 | analysis.R:1211 (once, via sapply) | Inline at Section 4.x intrinsic check (formula is one line) or keep for readability; analysis-only. |
| `frac_group_sets` | 6 | analysis.R:2144-2145 (twice) | Keep (two calls) but fix NA handling. |
| `coexpr_axis_cis_trans` | 6 | analysis.R:1057, 1079 | Keep (two calls, the algebraic identity deserves a name). Could use `crossprod`. |
| `allele_cor_boot_se`, `stratified_rho_sample`, `row_cor` | 14, 7, 6 | only inside `check_intrinsic_reliability` (row_cor also analysis.R:1235) | `stratified_rho_sample` and `allele_cor_boot_se` have a single caller, but `check_intrinsic_reliability` is already 33 lines; keep as separate helpers or fold into it (all analysis-only). `row_cor` keep (two callers). |
| `check_coexpr_reliability` and `check_intrinsic_reliability` | 23 + 33 | analysis.R:847, 1190 (once each) | Share the same bin/floor-summary/plot pattern (cut by quantile, `tapply` median, floor table). Factor a `.attn_floor_summary()` helper if desired. |
| `coexpr_acc_init/update/finalize` | 10 + 19 + 16 | cluster/scripts/coexpr_boot.R:63,69,80 only | Cluster-called, cannot be inlined into analysis.R. Could be collapsed with a generic accumulator over a vector of names (`c("total","cis","trans","dpar_sc","dpar_se")`) to remove the 5x repeated lines. |
| `coexpr_bootstrap_one`, `coexpr_perm_one` | 10, 36 | cluster scripts coexpr_boot.R:67, coexpr_perm.R:58 (also `coexpr_bootstrap_one` in unused Section 15 line 5085) | Cluster-called, keep. |
| `make_coexpr_draws` | 7 | analysis.R:834, 837; Section 15 line 5084 | Keep. |
| `make_coexpr_perm_draws` | 9 | analysis.R:898 (once) | Not cluster-called (draws are saved to the input .rda), so inlinable into analysis.R Section 4.6, but keeping preserves the draws-before-compute pattern. Keep. |
| `coexpr_class_table` / `coexpr_dom_class_table` | 11 / 9 | analysis.R:873, 874 (once each) | Nearly parallel (adjust p, classify, build data.frame). Could be one function with `type = c("reg","dom")`. Low priority. |
| `class_mean_var`, `class_anova` | 6, 11 | 16 each in analysis.R; `class_anova` also Section 12 line 4488 | Keep. The 16-call blocks in analysis.R:1269-1325 are loops candidates (`lapply` over the 8 class vectors) but not a function issue. |
| `coexpr_candidate_axes`, `coexpr_axis_mixtures`, `coexpr_axis_validate`, `coexpr_rank_check` | 15, 27, 21, 17 | analysis.R 4.7-4.9 (loop over `COEXPR.AXES`, once each) | Keep: each is called in a loop across five matrices. |
| Section 15 only callers | | `coexpr_bootstrap_one` also used by `coexpr_bootstrap` (5085); `shrink_cor` (5202); `coexpr_decompose` (5083) | These are used by the main path too, so not dead. |

## Other suggestions

- Magic numbers: 15 (table rows, `n_top`, `n_keep`, `N.KEEP`), 40 (`n_candidate`), 0.01 / 10 floors, 5 (minimum pole size), 0.5 posterior, 300 (B), 40 per bin. Define once in analysis.R and pass.
- Naming: parameters `sc`, `se` in `allele_cor_boot_se` / `check_intrinsic_reliability` read as standard error; `se` also shadows `se_obs`. Rename `res_sc`, `res_se`. Mixed code style in `plot_coexpr_cis_trans` (compact `x<-x[ok]`) vs `plot_coexpr_dom_class`.
- Indentation slips at 1552, 1656 and 1884 (` plot(` one space).
- `coexpr_acc_*` repeat five blocks of identical code; loop over `names`.
- Validation at entry is missing for most numeric arguments: `B`, `n_keep`, `n_per_bin`, `n_bins`, `floors`, `sig`, `probs` (not in (0,1) or unsorted), `alpha`, `k` in `coexpr_rank_check` (k < p), `n_top`.
- `coexpr_class_table` returns `total_padj` but the class ignores total; `coexpr_dom_class_table` ignores pair order alignment between dpar_sc and dpar_se (assumed equal).
- Plotting functions rely on globals (`COLOR.GREY`, `COLOR.LIST.1/2`, `REG.CLASS`, `DOM.CLASS`, `COLOR.ACCENT`, `org.Sc.sgd.db`, `normalmixEM`, `enrichGO`); `coexpr_axis_mixtures` could reuse `run_enrichment()` (2411) but that returns four ontologies, so a BP-only variant would be needed.
- `plot_coexpr_cis_trans` y-axis label "parents - hybrid" matches `trans = total - cis` only if read as (Sc-Se parents) minus (Sc-Se alleles); consider "trans: parental minus allelic change in correlation".

## Functions reviewed with no issues

`coexpr_decompose` (logic; comments updated), `make_coexpr_draws`, `coexpr_acc_init`, `coexpr_acc_update`, `row_cor`, `class_mean_var`, `coexpr_class_table`, `coexpr_dom_class_table`, `plot_coexpr_dom_class`, `seed_compare_core`, `coexpr_axis_cis_trans`, `make_coexpr_perm_draws`, `allele_cor_boot_se`.


---

# Appendix A4

# Review 4: R/functions.R lines 1977-2656 (Sections 8, 8b)

Sign conventions checked against `.OUT_MODES` builders (~l.683-694): dpar = log2(hybrid) - log2(parent); trans = parental ratio - hybrid ratio; cis = HYC.SC - HYC.SE. classify_dom's Over/Underdominant/Additive logic and build_reg_go_sets' "Sc if est > 0" are consistent with these.

## Discrepancies

| # | Severity | Function | Line | Issue / evidence | Direction | Numerical change? / sections |
|---|---|---|---|---|---|---|
| 1 | Possible bug | cluster_marker_enrichment | 2551 | `gene_vec` takes `m[, 2]` (positional) as the logFC vector. Repo rule: Seurat columns by name (FindMarkers column order changed in v5). Works only while avg_log2FC is column 2. Used for the enrichment background, up and down lists (values unused downstream except names, but wrong column silently). | `m[["avg_log2FC"]]` | No change if col 2 is avg_log2FC today. Called by cluster/scripts/cluster_stability.R:65 -> would need job resubmission only if edited and rerun; not required for results. |
| 2 | Possible bug | cluster_marker_enrichment | 2539-2542, 2528 | With `apply_fun = par_apply` (forked `mclapply`, cluster_stability.R:61), a failed worker returns a `try-error` string, not an error. `names(markers_list) <- ...` succeeds and the failure surfaces later as an obscure `$`/`[` error (or, for the enrichment mclapply, as `[[` on a try-error). Not caught per cluster. | Check `inherits(x, "try-error")` after each apply_fun and stop with the cluster ID. | No numerical change. Cluster job cluster_stability. |
| 3 | Possible bug | cluster_marker_enrichment | 2546-2547 | Marker thresholds are hard-coded magic numbers (`log2(1.25)`, `-log10(p_val_adj) > 20`). Also `FindMarkers` is not given `min.pct`/`logfc.threshold`, so the "background" is only genes passing Seurat defaults, not the full expressed set. | Hoist thresholds to arguments (defaults unchanged); document that the universe = tested genes. | No change if defaults kept. |
| 4 | Possible bug | class_overlap_heatmap | 2133 | Explicit `rng` smaller than the cell's |log2(obs/exp)| puts that value outside `breaks`; `image()` leaves such cells blank (unfilled), not clipped. analysis.R:556/660 use `shared_overlap_rng()` (same tables) so it is safe there; any other caller passing a smaller rng loses cells silently. | Clip `lor` to [-rng, rng] for `image()` (text keeps true value). | Plot only. |
| 5 | Possible bug | class_overlap_heatmap | 2150-2151 | `seq(0, 1, length.out = n)` for a 1-level axis gives `0` only; `image()` of a 1-column matrix and axis positions then mismatch (cell drawn across whole panel, label at 0). Also `n == 0` after NA filtering gives NaN `exp`. | Guard `nrow/ncol < 2` and `n == 0` (stop with message). | Plot only. |
| 6 | Possible bug | se_alpha_col | 2259-2261 | If `quantile(se, c(lo,hi))` gives equal values (all SE equal, or fewer than two finite SEs) the ratio is 0/0 = NaN, `pmin/pmax` keep NaN, `alpha` NaN and `rgb()` errors. | `w <- if (diff(rng) > 0) ... else 1`. | Plot only. |
| 7 | Possible bug | class_identity_overlap | 2203-2212 | `pe == 1` (a single level present in both vectors) gives 0/0 kappa = NaN and `kappa_p` = NA; `n == 0` makes `sample.int(0,0)` fine but `po` NaN. No validation of `nperm`, `n`. | Early return NA kappa with a message when `pe >= 1` or `n == 0`; validate `nperm` is a positive integer. | No change for current data. |
| 8 | Inconsistency (RNG) | class_identity_overlap | 2218 | `set.seed(seed)` resets the global RNG on every call (6 calls: analysis.R:734-737, 1109-1110). Any later code that relies on the global stream continues from a fixed state, and the RNG kind is process-global. Reproducible, but a hidden side effect. | Save/restore `.Random.seed` (on.exit) or use a local stream. | Would change downstream draws only if later code uses the global stream without its own seed. Affects no cluster section (all local). |
| 9 | Inconsistency | reg_class_vec | 2055-2059 | Unlike dom_class_vec, no check that `<q>_cis_q`/`_trans_q` columns exist. A missing column gives `NULL`, `classify_reg` sees `length(NULL) == 0` and returns a zero-row frame; callers then fail later or misalign silently. | Mirror dom_class_vec's `stop()` check. | None. |
| 10 | Inconsistency | classify_reg / classify_dom | 1992, 2025 | Arguments are named `p_*` and the comment says "permutation p-values", but every caller passes BH q-values (`PR$*_q`). Also no length check that p/est vectors agree (`length(p_cis)` drives `n`; a shorter `est` recycles silently via `sign(est_cis*est_trans)`). | Rename to `q_*` or comment "adjusted p"; add length check. | None. |
| 11 | Inconsistency | classify_reg | 1998 | `sgn >= 0` sends a gene with an exact-zero cis or trans estimate to "Cis + Trans" (rather than Ambiguous); only NA is Ambiguous. Almost never occurs with bootstrap estimates, but the rule is implicit. | Note in comment or use `> 0` and treat 0 as Ambiguous. | Negligible. |
| 12 | Inconsistency (stale label) | clean_reg | 2048-2052; header l.92 | Maps "Cis x Trans" -> "Compensatory", but classify_reg never emits "Cis x Trans" (only Compensatory). Dead branch kept from an older 6-class scheme. | Drop the first line (keep NA mapping) if no saved objects carry the label. | None. |
| 13 | Inconsistency | build_reg_go_sets | 2352 | `ifelse(est > 0, "Sc", "Se")`: `est == 0` is labeled "Se". Only matters for exact ties; NA est already gives NA. | `est > 0 -> "Sc"`, `est < 0 -> "Se"`, else NA. | Negligible. |
| 14 | Inconsistency | run_enrichment | 2411-2424 | `qval` feeds `qvalueCutoff` but `pvalueCutoff` stays at the clusterProfiler default 0.05 for GO/KEGG, so `qval = 0.2` is effectively limited by p < 0.05 first; n_sig_terms/print_enrich_brief then count `qvalue < 0.2` among terms that already passed p < 0.05. Comment mentions defaults but not this interplay. Also `simplify()` is applied to the GO result without guarding the NULL return (caught only by the tryCatch, hiding real errors such as a missing OrgDb). | Pass `pvalueCutoff = 1` or document; log the tryCatch error message. | Would change which terms appear (more terms). Re-run of 8b local only; ORA jobs at l.4322 share this function. |
| 15 | Note | run_enrichment comment | 2409 (old) | Example set name `Reinforcing_Sc` does not exist (REG.CLASS has no Reinforcing). Fixed in the comments JSON. | | None. |
| 16 | Note | plot_mean_bfreq_class | 2305-2328 | Name says "bfreq" but y can be bsize/kbal/cv2; `mode` is not `match.arg`-checked (a typo yields `NULL` column and a confusing error); `basis` for dom classes is not exposed; legend uses two different dash characters ("Reg — " vs "Dom – ") and non-ASCII text in code (Windows/CRLF encoding risk, also in `barplot_enrich_pair` l.2492,2498). `abline(lm(...))` on constant x gives NA coefficients and a warning. | Validate `mode`; use ASCII in strings. | Plot only. |
| 17 | Note | score_cell_cycle_by_cluster / score_modules_by_cluster | 2591-2607, 2632-2654 | `pdf()` is opened with no `on.exit(dev.off())`; any error in VlnPlot/barplot leaves the device open and a partial file (the very failure mode barplot_enrich_pair guards against). `score_modules_by_cluster`: with exactly one cluster `sapply(...)` returns a vector and `means[, i]` errors; with many modules `par(mfrow = c(1, k))` makes narrow panels; names containing `_` or spaces are altered by AddModuleScore column naming, so `paste0(names, "1")` may not match. | `on.exit`; `drop = FALSE`; validate names. | Plot only. |
| 18 | Note | go_gene_set | 2617-2619 | `select()` can return `NA` ORFs and prints the "1:many mapping" message; `unique(...)` keeps NA. | `na.omit`, `suppressMessages`. | AddModuleScore drops unknown features, so results unchanged. |
| 19 | Note | cluster_marker_enrichment | 2556-2558 | `up`/`down` calls pass `names(...)` as gene lists; ensure FindMarkers rownames are ORF IDs (Seurat v5 may convert `_` to `-` in feature names; yeast ORFs contain `-` rarely but some contain none). OK today; flag only if gene IDs change. | | |

No `&&`/`||` on possibly-vector operands found in this chunk (class_overlap_heatmap l.2115 uses `||` on two `is.null()` scalars; run_enrichment/n_sig_terms `||` on scalars). Dominance/regulatory semantics match project_context.md (hybrid diploid caveat is encoded in the ploidy-adjusted noise-axis calls). Not part of the cluster-dependent section list: classification uses PR from Section 2.3 output; fixing items 1-13 does not require resubmitting any cluster job except the optional cluster_stability rerun if item 1-3 are edited and results regenerated.

## Merge/collapse candidates

| Function | Lines | Callers | Recommendation |
|---|---|---|---|
| class_heatmap | 3 | analysis.R x11 (no cluster) | Thin wrapper over class_overlap_heatmap. Fold in by making class_overlap_heatmap's defaults (`levels_* = NULL`) the single entry point, and replace calls in analysis.R with `class_overlap_heatmap(...)`; or keep as an alias. Do together with plot_class_overlap (l.2829-2833, outside chunk), also a thin wrapper over the same core. |
| .heatmap_legend | 12 | only class_overlap_heatmap l.2162 | Keep as a helper or nest inside class_overlap_heatmap; not shared. |
| se_alpha_col | 13 | only plot_cis_trans_class l.2292 (x2) | Single caller; nest as a local helper inside plot_cis_trans_class. Not used by cluster scripts. |
| barplot_enrich_pair | 18 | only plot_cluster_marker_enrichment l.2575 (plus a comment at l.5271) | Single caller; merge into plot_cluster_marker_enrichment as an inner function. Not used by cluster scripts. |
| plot_cluster_marker_enrichment | 9 | analysis.R x4 (1825-1834) | Keep (analysis-level wrapper over a loop); absorbs barplot_enrich_pair. |
| summarize_go_sets | 11 | analysis.R:2160 only | Could inline at Section 8b of analysis.R as a data.frame call; it depends on n_sig_terms (kept). Low value either way. |
| n_sig_terms | 4 | analysis.R x2, summarize_go_sets x4 | Keep. |
| run_enrichment | 14 | l.2557-2558 (cluster_marker_enrichment), job runner l.4322; not in analysis.R | Keep; reached from cluster_stability via cluster_marker_enrichment, so it must stay in functions.R. |
| cluster_marker_enrichment | 31 | cluster/scripts/cluster_stability.R:65 only | Cannot be inlined into analysis.R (cluster job sources functions.R). Keep. |
| reg_class_vec | 5 | functions only (2282, 2346, 5329, + dom plot) | Keep; used by 3 other functions. |
| classify_reg / classify_dom | 19 / 20 | analysis.R x4 each; also coexpr tables (l.1730, 1752) and reg/dom_class_vec | Keep. classify_reg and classify_dom are near-duplicates (two p, two signs -> labels); could share a table-driven core, but the label rules differ, so a merge adds little. |
| build_reg_go_sets, build_component_go_sets | 19 / 13 | analysis.R x4 / x8 | Keep. analysis.R repeats the cis/trans pair per quantity (l.2125-2132): a `component = c("cis","trans")` loop or a `build_component_go_sets` that returns both components would remove 4 duplicate pairs. |
| go_gene_set | 3 | analysis.R x3 | Keep (shared by 3 gene sets). |
| score_cell_cycle_by_cluster / score_modules_by_cluster | 17 / 23 | analysis.R x4 each | Keep; both share the pdf+VlnPlot+mean-by-cluster pattern; a common `open_pdf_safely()` would fix item 17 once. |
| clean_reg | 5 | analysis.R x25 | Keep; drop the dead "Cis x Trans" line (item 12). |
| Section 15-only callers | | none found in this chunk's functions: all have real callers in analysis.R or functions.R | |

## Other suggestions

- Magic numbers: default `fdr = 0.01` for heatmaps vs `sig = 0.05` for classes; marker thresholds (item 3); `0.5` pseudocount in `lor` (state it in the heatmap comment; the displayed fold `2^lor` is the pseudocount-smoothed ratio, not raw obs/exp, so a zero cell prints about 0.3x-1x rather than 0); `n_top = 10`; `alpha_range`; `sum(sel) > 2` for slopes.
- class_overlap_heatmap's two-sided p is `2*min(P(X>=q), P(X<=q), 0.5)`: a doubled-tail test, slightly conservative vs the exact two-sided hypergeometric; fine, but state it. Multiple-testing BH is applied over all K x K cells including empty-margin cells (p = 1), which inflates the denominator slightly.
- Validate numeric inputs at entry per repo rule: `sig` in (0,1] for classify_*, `nperm`, `fdr`, `brk >= 2`, `q`, `n_top`.
- classify_reg/classify_dom: `colors[match(cls, REG.CLASS)]` returns NA color for "Ambiguous"; callers using `$color` should handle it (none found reading `$color` in analysis.R besides plotting by class names).
- `plot_cis_trans_class` legend lists all of `REG.CLASS` using `colors` in their given order; if a user-supplied `colors` is not in REG.CLASS order the legend is mis-colored. Index by `colors[REG.CLASS]`.
- reg_class_vec/dom_class_vec/build_*_go_sets repeat `i <- match(BURST.CONTRASTS$gene, PR$gene)`; a small `pr_col(PR, BC, name)` helper would deduplicate and could carry the missing-column check (item 9).
- The `levels` argument default `REG.CLASS` (global) in class_identity_overlap is evaluated lazily; all current callers pass it explicitly, so the default can be dropped.

## Functions reviewed with no issues

ckpt-independent helpers: n_sig_terms, print_enrich_brief, summarize_go_sets, summarize_class_overlap, plot_cluster_marker_enrichment, build_component_go_sets, classify_dom (logic; only naming note in item 10), dom_class_vec, class_heatmap, .heatmap_legend, build_reg_go_sets (apart from item 13).


---

# Appendix A5

# Review 5: functions.R lines 2657-3126 (8g, 8h, 9, intrinsic/extrinsic, ploidy)

Call sites verified with Grep. All 19 functions are reachable from analysis.R. None is called from cluster/scripts/*.R, except `ploidy_shift()`, which is only named in a comment in gene_perm.R; the object `PLOIDY.SHIFT` is saved by analysis.R:444. A scan of the `&&`/`||` conditions in this range found no vector-valued operand. The only conditions are `length(missing) > 0` and `is.null(...)`.

## Discrepancies

| # | Severity | Function, line | Issue and evidence | Direction | Numerical change / cluster sections |
|---|---|---|---|---|---|
| 1 | Likely bug | `intrinsic_fraction`, 3012 | `frac <- intr/(intr + pmax(extr,0))` does not clip `intr`. Poisson subtraction (3006) can make `intr` negative. With intr = -0.1 and extr = 0 the result is exactly 1 (checked in R), so a gene with no measurable intrinsic noise is scored fully intrinsic. With intr < 0 and extr > 0 the result is negative, or above 1 when intr + extr < 0. `plot_intrinsic_hist` (3111) then drops everything outside [0,1], so the histogram and the class medians are biased. `ploidy_shift` (3049) clips `intr` and `extr` at 0, so the two functions treat the same quantity differently. | Clip `intr` at 0, as in `ploidy_shift`, or set the fraction to NA when `intr + extr <= 0`. | Changes INTR.FRAC and Fig 8 (analysis.R:1337-1341). The `ploidy_shift` path is unaffected. No cluster sections (2.1/2.3/4/6.1/10). |
| 2 | Possible bug | `species_composition_bound`, 2789; analysis.R:2051 | The bound is `|d| * effect`, which is a shift in mean residual in SD units (E[Y given X shifted by d SD] = rho*d). The 8g effect sizes are medians over genes, so the bound is a typical-gene value, not an upper bound. The old comment (2797) said "upper bound". analysis.R:2034-2036 says "fraction of residual variance", but the code computes an SD shift, and the plot label (2051) says "fraction of residual SD". The text, label and statistic are three different things. A composition difference moves the residual variance by about (d*rho)^2, not d*rho. | Choose one interpretation. For a variance claim, use `d^2 * effect^2`. Align the analysis.R:2034-2036 comment with the plot label. | Affects COMP.BOUND.* and S_species_composition_bound.pdf only. No cluster sections. |
| 3 | Possible bug | `species_composition_bound` / report, analysis.R:2046 | For "Hybrid, Sc vs Se allele", HYB.SC and HYB.SE are the same hybrid cells scored on allele-specific counts. Both use `CONTRAST.EXPOS$HYB` (analysis.R:2021-2022). `wilcox.test(x1, x2)` is unpaired, which treats 2N observations as independent. Pooled Cohen's d is fine, but the p-value is not valid. In `cell_cycle_continuum_shared`, the pooled PCA also counts each cell twice. | Use `paired = TRUE` (by shared barcode) for the hybrid pair, or report d only and drop the p-value. | Presentation only. No cluster sections. |
| 4 | Possible bug | `metabolic_state_cluster`, 2714-2716 | (a) `kmeans(..., nstart=10)` has no `set.seed`. The k selection and the final labels depend on global RNG state, so the labels are not reproducible. `RunPCA` and others reseed the global RNG earlier, so the state is incidental. (b) The silhouette k-means fits are discarded and the final fit is redone, so the scored partition is not the returned one. (c) `dist(scores)` is n x n for every dataset, which is memory-heavy at 10^4+ cells. (d) `silhouette` is called unqualified (works only because analysis.R:134 attaches `cluster`). The old comment ("coarsest-preferred logic already used earlier") misdescribed `which.max`. | Save the k-means result per k and return the best one. Take a local seed argument or `withr`-style save and restore. Use `cluster::silhouette`. Optionally compute the silhouette on a subsample. | Changes MET.STATE.* labels and COV.DIAG met_* columns slightly. No cluster sections. |
| 5 | Possible bug | `covariate_noise_diagnostic`, 2738-2744 | (a) `eta_sq` returns 1 if `droplevels` leaves a single state (aov has then no factor row), and `kruskal.test` errors. The current call sites always use k >= 2, but the function does not check. (b) `cor()` and `cor.test()` each rank every gene, which doubles the cost. `cor.test(...)$estimate` carries rho. (c) `aov` per gene (thousands of genes x 4 datasets) is slow. Eta-squared has a closed form: `sum(tapply(r,g,sum)^2/n_g)-sum(r)^2/n`, divided by the total SS. (d) `mat[, cells]`, `expo[cells]` and `fit` assume `rownames(mat) == rownames(fit)`. (e) No validation that `cells` is non-empty. | Validate at entry: `length(levels(met))>=2`, `length(cells)>0`, and `identical(rownames(mat), rownames(fit))`. Vectorise the statistics. | Results unchanged (c, b). No cluster sections. |
| 6 | Inconsistency | `cell_cycle_continuum_shared`, 2689-2695 | Has no missing-column check (`cell_cycle_continuum` has one at 2667-2669). A missing S.Score gives an obscure Seurat error. | Share one validator. | None. |
| 7 | Inconsistency | `kbal_sig`, 2862-2872; `plot_burst_kinetics_sig`, 2876-2887 | Both re-derive quantities that `add_burst_contrasts` already stores: `y`/`sy` equal `kbal_<mode>_est`/`_se`. In `plot_burst_kinetics_sig`, `x`/`sx` reduce algebraically to `mean_<mode>_est`/`_se` (checked in R: max difference ~1e-16). The covFB/seS block is copy-pasted three times (add_burst_contrasts, kbal_sig, plot), and `plot_burst_kinetics_sig` relies on its own `ok` mask matching `kbal_sig`'s. | Read the stored columns. `kbal_sig` becomes about 5 lines and the plot function drops about 10. | None beyond ~1e-16. No cluster sections. |
| 8 | Inconsistency | `kbal_sig`, 2868-2869 | The gene-level threshold is nominal `p < sig` on thousands of genes. Everywhere else the pipeline classifies on BH `_q` columns (analysis.R:467-469). Also, `sy == 0` with `y != 0` gives p = 0 (significant), and `sy == 0` with `y == 0` gives NaN, which is "ns". | Apply `p.adjust(p, "BH")` (or label the panel "nominal p"). Guard `sy <= 0`. | Changes Fig 4 counts. No cluster sections. |
| 9 | Inconsistency | `ploidy_shift` vs `ploidy_coexpr_factor` / `ploidy_adjust_dpar` | Genes with a non-finite shift (`!good`, 3052) get NA in every adjusted dpar column (bfreq, bsize, kbal, cv2), so they drop out of classification (`is.finite` filters). `ploidy_coexpr_factor` gives the same genes f = 1, "unadjusted" (3073), so they stay in the co-expression analysis. | State the intent. Either exclude these genes from co-expression or accept the asymmetry. | The count is printed at analysis.R:816. |
| 10 | Inconsistency | `ploidy_coexpr_factor` comment, 3064-3067 | The comment said haploid correlation = hybrid x `sd_c/sd_h` and that g "lowers theta". The code (`f^2 = g(mu+th)/(mu+g th)`, f <= 1) is f = sd_h/sd_c with theta raised by 1/g. The derivation in the code is correct; only the comment was reversed. Fixed in the comments JSON. | Comment only. | None. |
| 11 | Note | `ploidy_coexpr_factor`, 3076 | f uses the rate `MU`, not per-cell `mu*exposure`. It is exact only when exposure is about 1. | Document as an approximation. | None. |
| 12 | Note | `ploidy_adjust_dpar`, 3092-3102 | `cv2_*_est` is recomputed from full-data fitted MU/DISP, while the other columns are bootstrap point estimates shifted by s. The adjusted cv2 and its `_se` come from different estimators, and `_raw` and adjusted cv2 differ by more than s. The raw twin of cv2 exists only if cv2 columns were already present. | Note it in the function comment, or shift the bootstrap cv2 by the same factor. | Small. 2.3 outputs are not read. |
| 13 | Note | `plot_dom_class`, 2950, 2959 | `lim` and the return value: if no finite genes remain, `max()` of an empty vector gives -Inf and a warning. The invisible return value is the class vector subset to the finite genes, so it is not aligned with the rows of BURST.CONTRASTS (callers ignore it). `sig` and `quantity` are not validated. | Return a full-length vector (NA outside `ok`) or document it. | None. |
| 14 | Note | `intrinsic_extrinsic_components`, 2995-3000 | Assumes `HYC.SC` and `HYC.SE` (and `HYT.*`) have identical row order and cell order matching `expos`; this is not checked. Magic number 1e-10 is repeated 4 times. | `stopifnot(identical(rownames(...)), ncol(...) == length(expos$HYC))`. | None. |
| 15 | Note | `species_composition_bound`, 2787 | Module scores (AddModuleScore) and S/G2M scores come from separately normalised objects (Sc and Se genomes), so between-object means are only approximately on a common scale. The shared PCA fixes this for the cell-cycle axis only, not for the three metabolic scores. | State the caveat in the comment. | None. |
| 16 | Note | Naming | `missing` (2667, 2709) shadows `base::missing`. functions.R labels these blocks 8g/8h; analysis.R calls them 7.7/7.8 and checkpoints them in section7. | Rename to `absent`. Cross-reference the numbering. | None. |

## Merge/collapse candidates

| Function | Lines | Call sites | Recommendation |
|---|---|---|---|
| `species_composition_bound` | 8 | only inside `species_composition_report` (2807, 2808) | Fold into `species_composition_report` as a local helper. Not used in cluster scripts. |
| `cell_cycle_continuum` + `_shared` | 8 + 7 | analysis.R:2009-2012 (4), analysis.R:2042-2043 (2) | Could be one function taking one or two objects (the single-object case is the pooled PCA of one set). Keep both if you want the "within vs between" contrast explicit. |
| `plot_class_overlap` | 3 | analysis.R:559-561, 663-665 | Thin wrapper of `class_overlap_heatmap` (which already defaults `cex_cell` etc. when levels are supplied). `class_heatmap` (f.txt:2170) is the other wrapper. Call `class_overlap_heatmap` directly and delete both wrappers. |
| `shared_overlap_rng` | 8 | analysis.R:556, 660 | Duplicates the `lor` computation at 2124-2125. Factor out a small `.overlap_lor(a, b)` helper used by both, or let `class_overlap_heatmap` return `lor` invisibly. |
| `kbal_sig` | 18 | only `plot_burst_kinetics_sig` (2875) | Reduce to 5 lines reading the stored columns (see #7), or merge into the plot function. |
| `intrinsic_fraction` | 6 | analysis.R:1337 | Thin wrapper over `intrinsic_extrinsic_components` (called only by it and `ploidy_shift`). Keep, since the components are used elsewhere; or add a `frac` column to the components table. |
| `ploidy_coexpr_factor` | 8 | analysis.R:814 | One caller. Could move inline into analysis.R section 7, but the derivation comment justifies keeping it. Not used in cluster scripts, so inlining is possible. |
| `plot_violins` | 18 | analysis.R:684 | One caller; keep, it is a self-contained ggplot. |
| `plot_burst_kinetics_sig` | 26 | analysis.R:671 | One caller; keep. |
| `ploidy_shift`, `ploidy_adjust_dpar` | 13, 17 | analysis.R:436, 465 | Keep. `ploidy_shift` is named in gene_perm.R, but only as a comment. |

None of these functions is exclusively in the unused Section 15.

## Other suggestions

- `plot_dom_class` and `plot_burst_kinetics_sig` both repeat the pattern "read est/se, set non-finite correlation to 0, mask, subset"; a shared reader returning the masked columns would remove the repeated 6-line blocks (also in `eiv_components` and `add_burst_contrasts`).
- Magic numbers: 1e-10 floors (`intrinsic_extrinsic_components`), `breaks = 40` (diagnostic hist), `nstart = 10`, `k_range = 2:4`, `0.05` thresholds repeated as literals in `covariate_noise_diagnostic` (2747-2750).
- Validation gaps at entry: `mode` in `kbal_sig` and `plot_burst_kinetics_sig` (a bad mode gives NULL columns and an empty plot), `sig` in (0,1), `brk`/`fdr` in `plot_class_overlap`, `k_range >= 2` in `metabolic_state_cluster`.
- `plot_covariate_noise_diagnostic` changes `par(mfrow)` without restoring it.
- `species_composition_bound` should guard against `NA` `effect_size` and an empty `x1`/`x2`.
- Header outline: `shared_overlap_rng` is missing, `plot_burst_kinetics_sig` is described as a barplot (it is a scatter), and `intrinsic_extrinsic_components` is listed under Section 12 although it is defined in Section 9. All three are fixed in the comments JSON.

## Functions reviewed with no issues

`plot_covariate_noise_diagnostic`, `species_composition_report`, `ploidy_shift` (formulas verified: s = 1 for equal alleles with fully private noise, 0 for fully shared; consistent sign with `ploidy_adjust_dpar`: bfreq - s, bsize + s, kbal - 2s), `plot_class_overlap`, `shared_overlap_rng`, `plot_violins`, `plot_intrinsic_hist`.


---

# Appendix A6

# Review 6: functions.R lines 3127-4019 (Section 10, promoter architecture)

Note: this is analysis.R Section 6 (6.1-6.5); comments inside the chunk refer to those numbers. Numbering of "Section 10" in functions.R differs from CLAUDE.md's "10 = power analysis" (functions.R Section 13). Worth a one-line cross-reference.

## Discrepancies

| # | Severity | Function | Line | Issue / evidence | Direction | Result change / sections |
|---|---|---|---|---|---|---|
| 1 | Likely bug | architecture_noise_check (call site) | fn 3994; analysis.R:1438-1443 | The six `response = "BSIZE"` calls omit `covariate`, so they use the default `covariate = "MU"`. The function's own comment (3962-3973) says this is algebraically the DISP-on-MU test with the coefficient sign flipped (log BSIZE = log MU - log DISP). The p-value and slope are identical to the DISP rows (r2_partial differs only through a different total SS). The 6 BSIZE rows of NOISE.VALIDATE are therefore redundant. | Pass `covariate = "DISP"` at the call sites, or make the default `covariate = if (response == "BSIZE") "DISP" else "MU"`. | Changes the 6 bsize rows of NOISE.VALIDATE. Local only; no cluster section affected. |
| 2 | Possible bug | poly_at_tract | 3374 (regex at 3372) | `gregexpr("[AT]+")` finds runs of any mixture of A and T (e.g. ATATATAT counts as 8). Poly(dA:dT) tracts are homopolymeric A (or T) runs. Toy check: `poly_at_tract("GGATATATATGGAAAAAGG")` returns length 8 at position 3 (the true A-run is 5). In AT-rich yeast promoters this inflates length and adds noise. | If homopolymer is intended use `"A+|T+"` (or `(A+|T+)`). If mixed A/T is deliberate, say so in the name/comment (proposed comment states the current behavior). | Changes polyat_len/pos/tract, polyat_delta, 6.2-6.5 results. No cluster section (NuPoP unaffected). |
| 3 | Possible bug | promoter_direction_test, concordance_by_magnitude, promoter_noise_candidates | 3702, 3744, 3855 | Concordant is defined as `sign(delta) == sign(est)` for every `quantity`. The derivation (3674-3701) holds only for bfreq_cis_est = log2 DISP_Sc - log2 DISP_Se. bsize_cis_est = mean - bfreq (line ~829), so under the burst-frequency mechanism the predicted sign for bsize is mostly opposite; the one-sided "greater" test would then test the wrong tail. For kbal (bfreq - bsize) no direction is derived at all. architecture_noise_check says explicitly that no direction is asserted for BSIZE. analysis.R:1451-1463 acknowledges the bsize mirror but not the sign. | Decide the predicted sign per quantity (e.g. a `sign` argument) or use a two-sided test for bsize/kbal. | Changes bsize/kbal p-values, bins, concordant column. Local only. |
| 4 | Possible bug (cosmetic) | plot_concordance_by_magnitude | 3784 | `names.arg = round(cb$mag_mean, 1)`. For occ_access_delta (|delta| roughly 0.01-0.3) labels collapse to 0, 0.1, 0.2 and are uninformative. | `signif(cb$mag_mean, 2)`. | None numerically. |
| 5 | Inconsistency | promoter_noise_candidates | 3820, 3855 | Comment said "top arch_frac of |delta|", but `arch_frac = 0.90` is used as a quantile probability (keeps about the top 10%). Name suggests the opposite. Cutoffs are also computed over all genes, not the cis subset. | Rename to `arch_quantile` or document (comment fix proposed). | None. |
| 6 | Inconsistency | promoter_noise_candidates | 3832-3834, 3857 | The gene-order requirement is documented but only row counts are checked (`nrow`). PR is not checked. Misaligned inputs of equal length would silently mislabel. Likewise promoter_direction_test and concordance_by_magnitude check no lengths, and `&` recycles silently on a mismatch. | Add `stopifnot(identical(BURST.CONTRASTS$gene, ARCH$gene), nrow(PR) == nrow(BURST.CONTRASTS))`; length checks on delta/est/reg_class. | None unless inputs are misaligned. |
| 7 | Inconsistency | promoter_divergence | 3663-3672 | `merge()` is an inner join on gene id: genes with NA scores stay. Comment said "scored in both species". It also needs `occ_score_*` to exist (else `m$occ_delta <- NULL - NULL` errors with a 0-row replacement). The data flow SCORE <- merge(OCC) is only in analysis.R:1416-1417. | Comment fixed; optionally check `c("occ_score_sc","occ_score_se") %in% names(m)`. | None. |
| 8 | Possible bug | extract_promoters | 3266 | `g$strand[i] == "+"` else-branch treats everything else ("-", "*", ".", NA) as minus (NA gives an `if` error). read_gff_genes does not validate strand. | Validate `strand %in% c("+","-")` at entry (or in read_gff_genes). | None for current data if all genes are +/-. |
| 9 | Note | extract_promoters | 3267, 3277 | The neighbor bound uses only the adjacent gene in start order (`g$end[i-1]`), not the running maximum end. A gene nested inside a longer upstream gene sees the nested gene's end, so the promoter can start inside the longer gene. Rare in yeast. | `cummax(g$end)` shifted by one. | Could shift a few promoters; affects 6.x and requires re-packaging NuPoP inputs (6.1 cluster job; coords check forces resubmission). |
| 10 | Note | tata_box_score | 3352 | Log-odds background is uniform 0.25, but yeast promoters are ~62% AT. Scores are inflated for AT-rich sequence, and the same background is used in both species, so the delta is somewhat protected. | State in methods, or use the promoter-set base composition. | Would change tata scores/deltas. |
| 11 | Note | poly_at_tract / score_promoters | 3370 | Tract length is searched over the whole promoter (variable, 50-300 bp), while TATA and NuPoP use the fixed 50-200 window. polyat_delta mixes real tract differences with prom_len differences between species (longer promoters have a longer expected maximum run). | Consider restricting to the same window or adding prom_len covariates in the tests. | Changes polyat deltas. |
| 12 | Note | nupop_occupancy_cluster | 3552 | Crash isolation relies on forking. `parallel::mclapply` with `mc.cores < 2` falls back to plain `lapply` (verified), so `cores = 1` (the function default) gives no protection from a NuPoP abort. Fine under nupop.sub (24 CPUs). | Use `max(2, cores)` or document. | None. |
| 13 | Note | nupop_occupancy_cluster | 3537 | No validation of numerics (`window_bp`, `flank`, `fallback_flank`, `min_core_bp`, `cores`, `species`, `model`) at entry (repo rule). `window_bp <= 0` would make `seq(by=)` error or loop. | Add `stopifnot` checks. Also in tata_box_score (`window` length 2, ordered) and extract_promoters (`max_bp >= min_bp >= 0`). | None. |
| 14 | Note | nupop.sub (outside chunk) | cluster/submit/nupop.sub | `module load r/4.2.2` while CLAUDE.md specifies r/4.5.2. Also cluster/scripts/nupop_occupancy.R header says chromosomes "that carry a promoter in both species"; nupop_cluster_inputs actually keeps chromosomes with a promoter in that species. | Update. | Would change nothing numerically. |
| 15 | Note | promoter_noise_candidates | 3867 | Error message narrates an "older version" of score_promoters; no `&&` hazard found. Message should state the required columns. Also a missing PR column (e.g. `kbal_cis_q`) silently yields a column-less table (`out[[nm]] <- NULL`). | Check `PR` names too. | None. |
| 16 | Note | concordance_by_magnitude | 3744 | Parameter named `bfreq_cis_est` but called with bsize/kbal estimates (analysis.R:1478-1483). Empty `ok` set (length 0) makes `glm` error. Wald CI degenerate at 0 or 1 (se = 0); Wilson would be better. `trend_p` is two-sided while the headline test is one-sided. | Rename to `cis_est`; guard `length(d) < 2`. | None. |

Checked and fine: 0/1-based and strand handling in extract_promoters (verified by hand on both strands, boundary clipping, empty cases); TATA window arithmetic (8-mer fully inside 50-200 bp upstream; tata_pos and nupop_window_score cover the same stretch); PWM columns each sum to 100; occ sign flip; `&&`/`||` in chunk (only 3288 `ext && nchar(seq) > 0`, both scalar); all.equal coordinate check in score_promoters_nupop does detect small shifts (verified with a 40-gene x 10 bp toy, all.equal averages over differing elements only); no Seurat column-by-position use; no RNG use.

## Merge/collapse candidates

| Function | Lines | Callers | Recommendation |
|---|---|---|---|
| tata_box_score | 18 | only score_promoters (3403) | Keep; `pwm`/`pseudocount` args never varied. Could be inlined into score_promoters as a local, but standalone is testable. |
| poly_at_tract | 8 | only score_promoters (3404) | Same: keep or fold in. |
| nupop_run_windows | 9 | nupop_occupancy_cluster 3560, 3569 | Keep (used twice, avoids duplicate mclapply). Cluster-only. |
| nupop_predict_window | 23 | nupop_run_windows only | Keep (isolation of working directory). |
| nupop_window_score | 11 | score_promoters_nupop 3627 only | Fold into score_promoters_nupop (no other user); low value. |
| nupop_cluster_inputs | 13 | analysis.R:1408-1409 (2) | Keep: must run locally (needs Biostrings genome). |
| nupop_occupancy_cluster | 60 | cluster/scripts/nupop_occupancy.R:44,48 | Cannot be inlined (cluster). Keep. |
| score_promoters_nupop | 14 | analysis.R:1414-1415 | Keep. |
| read_genome_fasta | 5 | analysis.R:1387-1388 | Thin wrapper; reasonable to keep (name cleanup). Could be inlined into 6.1. |
| promoter_divergence | 10 | analysis.R:1448 (1 call) | Single call, short; could inline into 6.3, but names the sign convention; keep. |
| plot_concordance_by_magnitude | 5 | analysis.R:1505-1513 (9) | Keep. |
| Eligibility logic (cis_flag & !is.na & sign != 0) | - | promoter_direction_test one_feature, concordance_by_magnitude, promoter_noise_candidates | Factor a helper `concordance_set(delta, est, reg_class, cis_classes)` returning logical `ok`; default `cis_classes` is repeated in 3 signatures (3702, 3744, 3855) - use one constant. |
| tata/polyat/occ blocks in promoter_noise_candidates | 3911-3946 | - | Three near-identical builders; could loop over a spec list (feature columns), cutting about 30 lines. |
| promoter_noise_candidates / promoter_direction_test / concordance_by_magnitude | - | analysis.R only (3, 3, 9 calls) | None callers in Section 15 or cluster. |
| Section 15 | - | none of this chunk's functions are Section-15-only | n/a |

## Other suggestions

- Magic numbers: 148 (NuPoP minimum length; used at 3494 only, name it), 80 (FASTA line width), 0.25 background, 1.96, 10 bins, `arch_frac = 0.90`.
- `score_promoters(tata_window=)` exposes only the window; `pwm`/`pseudocount` are unreachable. Either expose or drop the arguments.
- `ARCH` realignment `ARCH[match(BURST.CONTRASTS$gene, ARCH$gene), ]` (analysis.R:1449) is a precondition of three functions; consider doing it in a helper or asserting inside.
- architecture_noise_check uses coefficient name `"x[ok]"` (3014); renaming via a data.frame (`d <- data.frame(...)`) is more robust.
- nupop_occupancy_cluster scores whole chromosomes though only ~1.8 Mb of promoters are needed; restricting windows to those covering promoters would cut cluster time (result identical; would need resubmission only if changed).
- A blank line is missing between nupop_occupancy_cluster's closing brace (3594) and the nupop_window_score comment (3595); not a comment edit, so not in the JSON.
- Mean occupancy over the window accepts any number of non-NA positions (no minimum coverage), so a window with 1 valid base scores as if full.
- `read_gff_genes` takes the seqid/strand of the first exon only; a gene id appearing on two seqids or strands is merged silently. Add a uniqueness check.

## Functions reviewed with no issues

read_genome_fasta, nupop_cluster_inputs, nupop_predict_window, nupop_run_windows, nupop_window_score, score_promoters, score_promoters_nupop (comment fix only), read_gff_genes (comment fix only)


---

# Appendix A7

# Review 7: R/functions.R lines 4020-4740 (Sections 11, 12, 13)

Line numbers = commit 8c4f6d8 (CR-stripped copy f.txt). Cluster-dependent sections per CLAUDE.md: 2.1/2.3/4/6.1/10. Note that `cluster_stability.R` (Sections 7.3/7.4) and `go_enrich.R` also source functions.R but are NOT in the CLAUDE.md rerun list; changes to boot_ari_one/make_boot_idx/assemble_cluster_stability/go_enrich_one would need those jobs resubmitted too.

## Discrepancies

| # | Severity | Function (file:line) | Issue / evidence | Direction | Numeric change? Affected |
|---|---|---|---|---|---|
| 1 | Possible bug | `power_grid_row` 4646-4664 with `fit_offset_nb` / `fit_offset_nb_mm` | Observed statistic uses the full MLE size (`fit_offset_nb`, 4655-4656) but the permutation null uses the method-of-moments size (`fit_split_nb_mm`, 4661). The real-data test uses the MLE for both (`.fit_split` -> `.fit_one`, 268-296). The MM log2-ratio has a wider sampling distribution than the MLE one. Toy check (400+400 cells, EXPOSURE.CV 0.8, 200 datasets, 200 perms, SIZE.RATIO=1): SD of log2 ratio MLE vs MM = 0.49 vs 0.57 (mean 2, size 4) and 0.64 vs 0.84 (mean 0.5, size 1); FPR at p<0.05 came out 0.025 and 0.045 (SE ~0.015), i.e. conservative in the low-depth/low-size corner. Power is therefore underestimated and not an exact image of the gene-level test. | Use the same estimator for observed and null (MLE both, if cost allows, or MM both, then say the real test differs). The SIZE.RATIO=1 FPR output (already returned) should be reviewed per grid cell. | Yes, power values change. Section 10 + `power_grid.R` job resubmission. |
| 2 | Note (open issue) | `power_grid_row` 4628-4678 | Simulates ONE contrast: two independent groups, equal means, size ratio. It has no hybrid/paired-allele structure, no parental-vs-hybrid ratio contrast, and no `4/N_p` extra variance term, so it cannot show the trans-vs-cis asymmetry (var_trans = var_cis + 4/N_p) documented in project_context.md. `CELL.RATIO` (analysis.R:2511, set to 1) is the only handle on parent/hybrid depth imbalance. Also the second group's exposures are drawn independently of the first, whereas cis alleles share a cell (and thus exposure/extrinsic noise). | Add a cis-like (paired alleles, shared exposure) and a trans-like (difference of ratios) contrast to the simulation, or state in Section 10 that the grid is a parental-contrast power only. | Additive; Section 10 + `power_grid.R`. |
| 3 | Possible bug (robustness) | `power_grid_row` 4668-4672 | `i0 <- which(SIZE.RATIO == 1)`; `if (qi == i0)` errors when `SIZE.RATIO` lacks 1 (length 0) or has 1 twice (length 2), and `P.NULL <- P.ALL[, i0]` would be wrong. analysis.R:2506 includes 1, so currently fine. | `stopifnot(sum(SIZE.RATIO == 1) == 1)` at entry. Also validate `PI1` in (0,1), `NJ`, `NI`, `N.MIX` >= 1, `CELL.RATIO > 0` (N.SE.X >= 1). | No. |
| 4 | Possible bug (edge) | `fit_split_nb_mm` 4581-4584 (and `fit_split_nb` 5405) | `perm[(n1 + 1):length(perm)]` returns `c(n+1, n)` (NA plus one element) when `n1 == length(perm)`, instead of an empty group. Reached only if `round(N.CELLS*CELL.RATIO) == 0`. | `perm[-seq_len(n1)]`. | No. |
| 5 | Possible bug | `within_between_decomp` 4378-4411 | (a) `clusters` is not checked for NA; cells with NA are dropped by `which(cl == cc)` but `N <- ncol(mat)` still includes them, so weights `w` sum to < 1 and both components are biased. analysis.R:1906 builds `Idents(...)[colnames(CONTRAST.MATS$MIX.SC)]`, which gives NA for any contrast cell absent from the Seurat object. (b) `expo` is not validated (`expo <= 0`/non-finite give Inf/NaN in `a` and `shot_cell`). (c) For genes failing `ok_gene` (mean 0) only `ratio` is set NA; `var_within` / `var_between` remain numbers (`var_between` = 1 for an all-zero gene). | `stopifnot(!anyNA(clusters), all(is.finite(expo) & expo > 0))`; set `var_within`/`var_between` to NA when `!ok_gene`; use `N <- length(cl)`. | (a),(b) no change if inputs are clean; (c) changes only table rows for zero-mean genes. No cluster job. |
| 6 | Possible bug (bias) | `within_between_decomp` 4393-4400 | Within variance uses `rowMeans((x-cmean)^2)` (divisor n_c, not n_c-1) while the Poisson shot term is unbiased, so `var_within` is biased low by factor (n_c-1)/n_c (2% at the min_cells=50 floor). Between variance is not corrected for the sampling noise of the cluster means (shot/n_c), which inflates `var_between` for small clusters; comment 4361-4366 acknowledges the second. | Use n_c/(n_c-1) correction, and optionally subtract mean shot noise / n_c from the between term. | Yes, small; Section 7.5 only. |
| 7 | Possible bug (selection bias) | `within_between_decomp` 4403, `plot_within_between_hist` 4424, `plot_within_between_vs_quantity` 4444, `plot_within_between_cross_species` 4461, `plot_within_between_by_class` 4485-4486 | `ratio = pmax(var_within,0)/var_between` becomes 0 for genes with `var_within <= 0` (exactly the genes with the lowest within/between ratio); every plot then drops `ratio <= 0` (`x > 0` or `is.finite(log10(0))`). Reported fraction above 1, Spearman rho, ANOVA and class boxplots are conditioned on positive within-variance. | Report the number dropped (`n` and `n_dropped`), or plot on a pseudo-log scale; at minimum return the dropped count. | Yes (reported n, rho, p). Section 7.5 only. |
| 8 | Inconsistency | `plot_lines` 4711 vs `min_detectable_ratio` 4734 | Reference line drawn at power 0.9 ("common target"), but `min_detectable_ratio(target = 0.8)` default and analysis.R:2750 comment ("80% power") use 0.8; the figure line and the heatmap use different targets. | Pick one target (pass `target` to `plot_lines` or draw 0.8). | No (figures only). |
| 9 | Inconsistency | `cluster_stability_inputs` 4252 | `key` fingerprints tasks (res, sil, n_clusters), cell names, `B`, `seed`, but not `nfeatures`, `dims_n`, `metric`, or the count values. If only those change, `check_cluster_key()` still passes and stale ARI is accepted. Also `stopifnot(ncol(counts)==ncol(obj))` (4248) does not check cell-name order; `boot_ari_one` indexes `ref_clusters[idx]` positionally and `assemble_cluster_stability` relies on `apply_cluster_labels` name check only for the final object. | Add `nfeatures`, `dims_n`, `metric` to `key` and `stopifnot(identical(colnames(counts[[d]]), colnames(sweeps[[d]]$obj)))`. | No numerics. Forces cluster_stability.R resubmit once (key changes). |
| 10 | Inconsistency | `kegg_local` 4295 | Selects the pathway-name column by position (`KEGGPATHID2NAME[[2]]`); repo rule says select by name. | Use the column name (`$to` in clusterProfiler's download_KEGG output; confirm on cluster-side version). | No. |
| 11 | Note | `compare_distance_metrics` 4342 | Uses `cl[[1]]`, `cl[[2]]`, silently ignoring any metric beyond two (or erroring with one). | `stopifnot(length(metrics) == 2)`. | No. |
| 12 | Note | `make_boot_idx` 4163 | `set.seed(seed)` overwrites the caller's global RNG state; analysis.R:1758 calls it (via `cluster_stability_inputs`) mid-script, so later unseeded stochastic steps depend on it. Draws themselves are correct: all indices are drawn before any Seurat call, so the RunPCA reseed (seed.use 42) cannot collapse them; each replicate gets a distinct column and `cluster_stability.R:40` passes `idx[, b]` (verified). Same seed for all datasets is fine (different n). | Optional: save/restore `.Random.seed`. | No (Section 7.3 unaffected unless RNG consumers downstream change). |
| 13 | Note | `assemble_cluster_stability` 4285-4286 | If every replicate of every candidate is NA, `mean` gives NaN, `which.max` returns `integer(0)` and the function errors with an unclear message; `n_ok` is not checked against `B`. | `stopifnot(any(is.finite(table$boot_mean_ari)))` with message; warn when `n_ok < B`. | No. |
| 14 | Note | `cluster_stability_inputs` 4240 | `role = "highest silhouette"` labels the coarsest resolution on the highest-silhouette cluster-count plateau (`plateau_coarsest`), not necessarily the argmax. Name is misleading in reports (`report_bootstrap_compare`). | Rename role to "plateau coarsest" (update `report_bootstrap_compare` and analysis.R consumers if any read role). | No. |
| 15 | Note | `plateau_coarsest` 4203-4208 | Runs are contiguous in the guard-filtered grid `ok`; excluded resolutions between two same-count points are skipped, so a run can bridge a gap. Harmless with the 0.05 grid but undocumented. | Mention or group on the unfiltered grid. | No. |
| 16 | Note | `hvg_elbow` 4045-4055 | `drop_frac = 0.001` and `floor = 1` are magic numbers; the "last" drop above 0.1% of the range is typically deep in the tail, so `n_features` is large and sensitive to `drop_frac`. Behaviour is documented, not wrong. | Document the chosen values in analysis.R. | No. |
| 17 | Note | `plot_within_between_by_class` 4486-4489 | `ok` filters on raw `class` NA, so "Ambiguous" genes pass and only become NA after `factor(levels=class_levels)`; `boxplot(x ~ cl, col = colors)` drops empty levels, so with any empty class `colors` shift off their levels. | `cl <- droplevels(cl)` and `colors[match(levels(cl), class_levels)]`. | Figure only. |
| 18 | Note | `plot_within_between_vs_quantity` 4449 / `_cross_species` 4468 | Reports Spearman rho but draws an OLS line (`lm`); on log10 ratios heavy-tailed points can make the line disagree with rho. | Use a rank-based or loess guide, or label as OLS. | Figure only. |
| 19 | Note | `size_log2_ratio` 4599 | Returns a *named* numeric (name `size`) on the finite branch, unnamed `NA_real_` otherwise. Harmless today. | `unname(...)`. | No. |
| 20 | Note | `fit_offset_nb_mm` 4575 vs `fit_offset_nb` 4553 | MLE caps theta > 1e6 to `Inf` (-> NA ratio, dropped from the null) whereas MM returns any finite theta; the two estimators handle near-Poisson groups differently, which feeds #1. | Cap MM identically or document. | Part of #1. |
| 21 | Note | Stale cross-references | 4152 "to_seurat_counts() in the driver" ok; 4421 and 4483 cite "Section 8" for silhouette/bootstrap reporting, which live in Section 7.3 of analysis.R; 4526 cites ".fit_one() above" (it is ~4,200 lines above, Section 1); outline lines 155, 164-166 list functions in the wrong section (bootstrap_compare_resolutions is Section 15 line 5125; intrinsic_extrinsic_components is at 2993). Fixed in the comments JSON except line 164 (outside this chunk). | - | No. |

Verified OK (no issue): ARI direction (`as.integer(ref_clusters[idx])` vs `Idents(boot_obj)`, same draw order, make.unique barcodes), fork safety of `boot_ari_one` (no global state, errors caught per job in `cluster_stability.R:40`), `power_grid_row` seeding (`SEED.BASE + i` / `+ 1e6 + row$n`: results independent of how rows are split across array tasks and workers), BH mixture indexing (`q[NJ + seq_len(n_alt)]`, `n_alt/(NJ+n_alt) = PI1`), NB MM formula (E[Pearson] = (n-1) + sum(mu)/theta), `perm_pval` add-one two-sided form, `&&`/`||` uses in this chunk (all scalar: `length(sizes) < 2 || min(sizes) < min_cells`, `is.finite(a) && a > 0 && ...`, `!is.finite(obs) || sum(ok) < 1`).

Relation to the trans/cis power asymmetry: only #2 (grid has no trans/cis contrast; `var_trans = var_cis + 4/N_p` not represented). Fano-after-thinning compression is not modeled either (simulation draws NB at the target mean directly).

## Merge/collapse candidates

| Function | Lines | Call sites | Recommendation |
|---|---|---|---|
| `fit_offset_nb` + `fit_offset_nb_mm` | 23 + 16 | `power_grid_row` 4646, 4655; `fit_split_nb_mm` 4583; `fit_split_nb` 5406 (Sec 15) | Near-duplicates (identical validation, mu, Pearson pre-check). Merge into one `fit_offset_nb(y, expo, method = c("mle","mm"))`. Also `fit_offset_nb` duplicates `.fit_one` (268-292) apart from the result name `size` vs `disp`; a single function with a `name` argument (or letting power code read `disp`) removes it. Cluster-sourced (power_grid.R), so must stay in functions.R. |
| `fit_split_nb_mm` | 4 | `power_grid_row` 4661 only; sibling `fit_split_nb` (5405, Section 15 only) and `.fit_split` (294) | Three copies of the same split. Collapse to one `fit_split(y, expo, perm, n1, fit)` or inline the 2-line split into `power_grid_row`. Stays in functions.R (cluster). |
| `perm_pval` | 4 | `power_grid_row` 4664; inner `pval()` at 1170 is identical | Have the inner `pval` at 1170 call `perm_pval` (or delete it). Keep (cluster-sourced). |
| `size_log2_ratio` | 3 | `power_grid_row` 4656, 4662; mirrors inner `rat()` 1117 | Keep (2 calls, cluster path). |
| `cluster_cols` / `line_colors` | 1 each | `cluster_cols`: functions.R 2597, 2638, 4697, 5420; `line_colors`: 4709, 4725, 5419 | Trivial wrappers over `colorRampPalette(<palette>)(n)`; could be one `ramp_cols(palette, n)`. Low value; keep. Used only in local code paths. |
| `diet_for_markers` | 3 | analysis.R:1765 only | Thin wrapper; the `else DietSeurat(obj)` branch is dead under the Seurat v5 requirement. Inline in analysis.R Section 7.3 as `DietSeurat(obj, layers = c("counts","data"))`. Not used by cluster scripts, so inlining is allowed. |
| `apply_cluster_labels` | 5 | `assemble_cluster_stability` 4286 only | Inline there. Note `assemble_cluster_stability` is called from cluster_stability.R:63, but `apply_cluster_labels` is only reached through it, so it can be folded into it. |
| `report_bootstrap_compare` | 9 | analysis.R:1773 only (inside the dataset loop) | Could be inlined into the 7.3 loop; keep if the logging should stay testable. |
| `compare_distance_metrics` | 8 | analysis.R:1778 only | Single use on MIX.SC; inline into 7.3 (cluster scripts do not call it) or keep as a documented one-off. |
| `min_detectable_ratio` | 5 | analysis.R:2753 only (inside a triple loop) | Inline as `ratios[which(power >= 0.8)[1]]` in 10.8, or keep; not called by power_grid.R, so inlining is allowed. |
| `kegg_local`, `load_cluster_output`, `check_cluster_key` | 4 / 6 / 4 | analysis.R:1766, 2151 / 1769, 2154 / 1770, 2155 | Used twice each, generic cluster round-trip helpers. Keep, but move out of Section 11 into a "cluster round-trip helpers" block (they are not clustering diagnostics). |
| `go_enrich_one` | 5 | `cluster/scripts/go_enrich.R:43` only | Cannot be inlined into analysis.R (cluster script sources functions.R). Keep; belongs with the enrichment functions (~2411). |
| `make_boot_idx`, `boot_ari_one` | 4 / 14 | `boot_ari_one`: cluster_stability.R:40, functions.R 5178 (Sec 15); `make_boot_idx`: 4251, 5176 (Sec 15) | Keep; cluster-sourced. Sec 15 `bootstrap_cluster_stability` is the only other caller. |
| `plateau_coarsest` | 6 | `cluster_stability_inputs` 4238; Sec 15 `bootstrap_compare_resolutions` 5127 | Keep (single logic shared with the unused serial path). |
| `legend_page`, `plot_lines`, `open_grid_pdf`, `umap_plot` | - | analysis.R (12, 12, 15, 7 sites) | Keep. |
| `hvg_elbow`, `plot_hvg_elbow`, `sweep_cluster_resolution`, `plot_resolution_sweep` | - | 7 sites each in analysis.R | Keep; the seven near-identical call blocks (analysis.R:1596-1633, 1707-1745) are better as a `for (d in DS.NAMES)` loop, not in functions.R. |
| `plot_resolution_sweep` | 10 | takes `min_cells` separately from `sweep` | Store `min_cells` in the `sweep_cluster_resolution` return list and drop the argument. |

Only-Section-15 callers inside this chunk: none of the chunk's functions has only Section 15 callers; `make_boot_idx`, `boot_ari_one`, `plateau_coarsest`, `fit_offset_nb` also have Section 15 callers in addition to live ones.

## Other suggestions

- Section 11/12/13 header styles differ (`####` banner vs `## ====`); unify.
- Section 13 mixes simulation code and generic plotting helpers (`cluster_cols`, `umap_plot`, `line_colors`); `cluster_cols`/`umap_plot` belong with Section 11 (clustering), and a cluster-round-trip helper block would hold `kegg_local`/`load_cluster_output`/`check_cluster_key`/`go_enrich_one`.
- Outline (lines ~150-180) is stale for Section 11 (missing 10 functions, lists `bootstrap_compare_resolutions` and `intrinsic_extrinsic_components` in the wrong section); replacement lines are in the comments JSON (line 164 is left to the chunk that owns 2993).
- Magic numbers: `drop_frac = 0.001`, `floor = 1`, `min_cells = 50`, `tol = 0.02`, `res_grid`, `nPermSimple = 100000` (`go_enrich_one`), `1e-10` floors in `within_between_decomp`, optimize range `c(-4, 15)` and cap `1e6` (duplicated in `.fit_one`), `ncol = ceiling(n/15)` in `legend_page`. Hoist the shared theta bounds to constants used by `.fit_one`, `fit_offset_nb`, `fit_offset_nb_mm`.
- Input validation at function entry is missing in: `hvg_elbow` (`drop_frac` in (0,1), `floor` numeric), `sweep_cluster_resolution` (`dims` subset of PCs, `res_grid` > 0, `min_cells`, `tol`), `make_boot_idx` (`n`, `B` positive integers), `boot_ari_one` (`idx` within 1:ncol, `length(ref_clusters) == ncol(counts)`), `fit_offset_nb*` (numeric `y`/`expo`, equal length), `power_grid_row` (see #3), `within_between_decomp` (see #5; also `nrow(mat) > 0`).
- `plot_hvg_elbow` / `plot_resolution_sweep` open `pdf()` and rely on `dev.off()` at the end; use `on.exit(dev.off())` so an error does not leak a device. `report_within_between_by_class` has the same pattern.
- `plot_resolution_sweep` and `plot_hvg_elbow` assume `fig_dir/extra/` exists; `dir.create(recursive = TRUE)` or fail with a clear message.
- `power_grid_row` runtime: the inner loop makes `NJ * length(SIZE.RATIO) * NI` calls to `fit_split_nb_mm`, each re-running `is.finite` and `keep` filtering on both halves; a vectorised MM estimator over a permutation matrix would cut the cluster cost, at the price of changing nothing numerically.
- `within_between_decomp` allocates three dense genes x cells copies (`a`, `ap`, `shot_cell`); fine for current sizes, but note if run on sparse input (`sweep` densifies).
- `cluster_stability_inputs` `role` strings and `table$role == "chosen"` in `report_bootstrap_compare` are string-coupled; use constants.
- The CLAUDE.md rerun protocol lists cluster-dependent sections but omits the `cluster_stability.R` and `go_enrich.R` jobs; consider adding them since functions in this chunk feed both.

## Functions reviewed with no issues

plot_hvg_elbow, sweep_cluster_resolution (logic; see #14/#15 for neighbors), plot_resolution_sweep, report_bootstrap_compare, diet_for_markers, apply_cluster_labels, load_cluster_output, check_cluster_key, go_enrich_one, plot_within_between_hist (apart from #7), open_grid_pdf, cluster_cols, umap_plot, line_colors, legend_page, min_detectable_ratio, perm_pval, fit_offset_nb (logic correct), fit_offset_nb_mm (formula correct; see #20).


---

# Appendix A8

# Review 8: R/functions.R lines 4741-5436 (Section 14 external validation, Section 15 unused)

Line numbers refer to commit 8c4f6d8. No code was changed. Comment edits are in `review_8_comments.json` (32 edits).

## Discrepancies

| # | Severity | Function | Line | Finding |
|---|---|---|---|---|
| 1 | Likely bug (Section 15, unused) | `boot_disp_logse` | 5188 | `neg_binom_fit_offset(y[sample.int(n, n, TRUE)], exposure[sample.int(n, n, TRUE)])` draws two independent index vectors. Counts and exposures are resampled separately, so each replicate pairs cell i's count with cell j's library size. This inflates the apparent dispersion and its SE. Fix: draw one `idx <- sample.int(n, n, TRUE)` and use it for both. `refine_by_boundary()` (5093) inherits this through `boundary_frac`. Neither is called, so no current result changes and no cluster section (2.1/2.3/4/6.1/10) is affected. The pipeline's own bootstraps (`make_boot_ky()` / `make_draws()`) resample pairs. |
| 2 | Possible bug | `map_to_orf` / `ORF.PATTERN` | 4896-4925 | The pattern accepts a `.` before the dubious-ORF suffix (`YAL047C.A`), but ids that match are passed through unchanged. The period form therefore stays `YAL047C.A`. `NB.SC$ORF` uses the SGD hyphen form, so `merge(NB.SC, JACKSON, by = "ORF")` (analysis.R:2342 and `NEW.MERGE`) drops those genes. `JACKSON.ORF.KEEP` (analysis.R:2318) counts them as resolved, so they stay in the matrix but never merge. Fix: `out[is.orf] <- sub("\\.", "-", ids[is.orf])`. The comment now states the actual behavior. Effect: a few dozen `-A` genes enter the Jackson rows of `NEW.CORR` and `ALL.CORR.PAIRWISE`. This is a local Section 9 change only, with no cluster sections affected. |
| 3 | Possible bug | `add_burst_terms` | 4772-4779 | Fano = Mean x CV2 is a Fano factor only when Mean is on a count scale. Newman (GFP, AU), Stewart-Ornstein (analysis.R:2238, "Mean Expression (AU)") and the Keren YFP data are fluorescence in arbitrary units. For those sources the `fano > 1` gate, BSIZE = Fano - 1 and BFREQ = Mean/(Fano - 1) depend on the instrument's units, because the "-1" is not scale-free. The NB identity holds for MIX.SC and the count-based scRNA sources only. The BSIZE and BFREQ rows of `EXT.CORR` and `ALL.CORR.PAIRWISE` for these three sources are therefore not on a common scale with MIX.SC. Spearman rho on `Mean` and `CV2` is unaffected. The comment now states the caveat. Direction: treat the protein-source burst rows as descriptive, or drop them. Section 9 only. |
| 4 | Inconsistency | `add_burst_terms` | 4775 | `d$Fano <- fano` overwrites a pre-existing column. `KEREN` carries the source's own `Fano` (analysis.R:2236), which is replaced silently by Mean x CV2. Harmless if the two agree. Check once and say so. Section 9 only. |
| 5 | Possible bug (edge) | `cor_row` | 4749-4756 | With fewer than 2 finite pairs `cor.test` stops ("not enough finite observations"). With exactly 2 pairs it returns p = NaN with a warning that `suppressWarnings` hides. With a constant vector rho = NA. An error inside the `do.call(rbind, lapply(...))` loops (analysis.R:2352, 2375, 2395) aborts Section 9. No input validation (equal lengths, n >= 3). Direction: return `data.frame(n, NA, NA)` when `sum(keep) < 3`. Reachable only for a nearly empty ORF overlap. No numeric change today. |
| 6 | Note | `fit_source` | 4958 | No input validation (CLAUDE.md "validate numeric inputs at function entry"). An empty matrix after QC, or `mean(colSums(mat)) == 0`, gives NaN exposure and every gene NA, silently. `cut(seq_len(n.genes), n.workers)` on a very small matrix gives fewer chunks than workers (harmless). Add `stopifnot(is.matrix(mat), is.numeric(mat), ncol(mat) >= 2, nrow(mat) >= 1)`. |
| 7 | Note | `qc_filter_counts` | 4854 | Thresholds not validated (numeric, length 1, >= 0). Cells are filtered first and genes second, so the gene criterion applies to the surviving cells. The old comment claimed these are "the same two checks applied before any NB fit elsewhere in the pipeline"; this was not verified and is removed from the comment. |
| 8 | Note | `mean_adjusted_noise` | 5020 | The argument `mean` shadows `base::mean` inside the function (harmless here). `span` is not validated. `loess` on fewer than about 10 finite points errors, and there is no guard for `sum(ok) < 10`. Both callers pass thousands of genes. Section 9.5 only. |
| 9 | Note | `to_numeric_matrix` | 4868-4883 | With a one-row character matrix `apply(mat, 2, as.numeric)` returns a vector and `num[, !bad, drop = FALSE]` fails. Not reachable with the current four sources. |
| 10 | Note | `read_header_line` / `get_data_header` | 4790, 4810 | `strsplit` drops a trailing empty field, so a header ending in a tab is one field short and `get_data_header` may stop with its mismatch error. This fails loudly rather than silently. The Jackson header is cleaned with `trimws` in analysis.R:2317 only for the ORF test, while `fread_matrix` names its columns with the raw header. A header field with stray whitespace would pass `map_to_orf` in the first call and fail it in `collapse_to_orf`. |
| 11 | Note | `plot_geneset_direction_stack` comment | 5288-5323 | The comment said "gated by its own permutation p" (the code reads `PR[[quantity_total_q]]`, the BH q), described `quantity` as "mean" or "bfreq" (the code allows mean/bfreq/bsize/kbal), and cited a "Reinforcing" class that is not in `REG.CLASS`. Corrected in the comment edits. The code needs no change. Not broken if called. |
| 12 | Note | `bootstrap_compare_resolutions` vs `assemble_cluster_stability` | 5125, 4277 | The local table has 7 columns and the cluster table 8 (`n_ok`). The local path uses `mean(ari)` / `min(ari)` without `na.rm`, while the cluster path uses `na.rm = TRUE`. Both feed `report_bootstrap_compare()` safely (it reads only `res`, `role`, `boot_mean_ari`). Cosmetic divergence in an unused function. |
| 13 | Note | `.cohen_kappa` | 5259 | Duplicates the kappa formula in `class_identity_overlap()` (line 2207) and is never called. The old comment described a rationale that belongs to `class_identity_overlap` and referred to "class_heatmap's per-cell test". Rewritten. See Merge/collapse. |

### Section 15 interactive-call check (signatures against wrapped functions)

I checked each against the definitions and the cluster scripts.
- `pilot_split_se`, `boot_contrasts`, `permute_contrasts` (with `...` for `ploidy_shift`), `coexpr_bootstrap` (`coexpr_decompose`, `make_coexpr_draws(nSC, nSE, nH, B, seed)`, `coexpr_bootstrap_one(draw, resid)`, `assemble_coexpr_bootstrap(resid, pt, draw_results)`) and `bootstrap_cluster_stability` (`boot_ari_one` positional order) all match their targets' current signatures and return shapes.
- `assemble_coexpr_bootstrap` yields the same CB layout as `coexpr_acc_finalize()`: est/se/z/p, plus `lambda`.
- `bootstrap_compare_resolutions` uses `sweep$grid` (`ok`, `res`, `sil`, `n_clusters`), `sweep$chosen_res` and `sweep$obj`, all of which `sweep_cluster_resolution()` returns.
- `chk` and `chk_prec` use the `MEAN_CT`, `DISP` and `DISP_LOGSE` columns of `fit_counts_offset()`, and they exist.
- `coexpr_raw_cor` uses `shrink_cor(t(...))`, the same orientation as `coexpr_decompose`.
- `coexpr_gene_degree` reads the `coexpr_class_table()` columns `gene_i`, `gene_j`, `total_padj`, `cis_padj`, `trans_padj` and `class`, and they exist.
- `plot_coexpr_pair` uses `SPECIES.COLOR[["Sc"/"Se"]]`, defined in analysis.R:150 before the function is called.
- `report_cluster_marker_enrichment` matches the `cluster_marker_enrichment()` return (`up_enrich`, `down_enrich`, `cluster_ids`) and `print_enrich_brief(e, q)`.
- `plot_geneset_direction_stack` matches `reg_class_vec()` and uses the same `_total_est` / `_total_q` columns as `build_reg_go_sets`. `COLOR.GREY[["dark"]]` and `[["mid"]]` exist.
- `fit_split_nb` calls `fit_offset_nb()`, which exists with the matching `c(mu, size)` return.
- `plot_palette_swatches` uses `COLOR.LIST.1/2/3`, `COLOR.SEQ`, `COLOR.PHASE`, `COLOR.GREY`, `COLOR.ACCENT`, `line_colors`, `cluster_cols`. All are defined, but only in analysis.R (lines 150-168), so it works interactively after sourcing analysis.R and not from `functions.R` alone.
- The only broken function is `boot_disp_logse` (finding 1), and then silently, not by error.
- `permute_contrasts`, `boot_contrasts` and `pilot_split_se` add no error handling. The cluster scripts check for try-errors, the serial wrappers do not.

### Algebra check (Section 14)

NB gives CV^2 = 1/mu + 1/theta, so Fano = mu x CV^2 = 1 + mu/theta, burst size = Fano - 1 = mu/theta, and burst frequency = mu/(Fano - 1) = theta. This matches `NB.SC$BFREQ <- DISP` and `BSIZE <- MU/DISP` in analysis.R:2227-2228 and `fit_counts_offset()` (BFREQ = size, BSIZE = mu/size), so the algebra is correct.
- `fit_source` returns `CV2 = fit$CV^2` with `CV = sqrt(1/size + 1/mu)`, which is right.
- DISP = Inf gives BSIZE 0 and BFREQ Inf. Both drop out of `cor_row` through `is.finite` after `log`.
- Exposure = colSums/mean(colSums) makes MU count per average cell, as documented.
- `mean_adjusted_noise` returns a loess residual on log scale. `predict(fit)` returns one value per training point, aligned with `ok`. This is correct.

## Merge/collapse candidates

| Function | Call sites | Lines | Recommendation |
|---|---|---|---|
| `read_header_line` | analysis.R:2301 (1), `get_data_header` 4811 | 5 | Keep as is. Its only direct use is the Nadal file, whose header is not one field short. Alternative: fold into `get_data_header(raw = TRUE)`. Low value. |
| `get_data_header` | analysis.R:2287, 2316, `fread_matrix` 4838 | 15 | Keep. Three callers. |
| `collapse_to_orf` | analysis.R:2293, 2305, 2330, 2340 | 6 | Keep. Four callers, and it wraps `map_to_orf`. |
| `map_to_orf` | analysis.R:2318 (1), `collapse_to_orf` | 28 | Keep. Used directly for `JACKSON.ORF.KEEP`. |
| `add_burst_terms` | analysis.R:2259-2261 (3) | 8 | Could be inlined as a 3-line loop in analysis.R Section 9.2. Keep: a single named transform is clearer than three copies. Not called from cluster scripts, so inlining is possible. |
| `mean_adjusted_noise` | analysis.R:2418, 2420 (2) | 9 | Keep (short, two uses). |
| `fit_source` | analysis.R:2295, 2307, 2332, 2342 (4) | 54 | Keep. The pipeline's only per-source fit. |
| `cor_row` | 14 call sites in analysis.R | 8 | Keep. Many uses. The repeated `mean/cv2/bfreq/bsize` blocks in analysis.R 9.4 (three near-identical `rbind(cor_row(...))` blocks) could share a helper, but that is in analysis.R. |
| `.cohen_kappa` | none | 6 | Delete, or have `class_identity_overlap()` (line 2207) call it for `kappa_obs`. It is a duplicate of the formula inlined there. |
| `fit_split_nb` | none | 4 | Delete, or merge with `.fit_split` (294) and `fit_split_nb_mm` (4581) into one `fit_split(y, expo, perm, n1, fit_fun)`. Three copies of the same 3-line body. `.fit_split` and `fit_split_nb_mm` are called from cluster scripts, so a merged function must stay in functions.R. |
| `chk` / `chk_prec` | none | 4 each | Merge into one `chk(fit, what = c("finite", "prec"))`, or drop. Same `cut` call in both. |
| `pilot_split_se`, `boot_contrasts`, `permute_contrasts` | none | 2 each | Thin wrappers over `lapply` plus `do.call(rbind)`. Keep only if interactive use is wanted. Otherwise delete. |
| `assemble_coexpr_bootstrap` | `coexpr_bootstrap` only | 21 | Duplicates `coexpr_acc_finalize()` (1493) in structure. A single "finalize" taking either a stack of draws or running sums would remove one of them. `coexpr_acc_finalize` is cluster-used, so keep it. |
| `bootstrap_compare_resolutions` | none | 18 | Overlaps `cluster_stability_inputs` + `assemble_cluster_stability`. Could be deleted together with `bootstrap_cluster_stability` if no local fallback is wanted. |

## Other suggestions

- Header outline: Section 14 lists only 2 of its 11 functions; the edit adds the other nine, and Section 15 now says "kept for interactive use".
- Section 15 layout: the preamble promises 15a/15b/15c. In the file, `refine_by_boundary`, `chk` and `chk_prec` (15b) come right after the 15a functions without a marker, and `bootstrap_compare_resolutions` (15a) follows them. An `## ---- 15b ----` marker is inserted before `refine_by_boundary`. The "(continued)" markers then read correctly. Reordering the code would make the layout cleaner and was not proposed.
- Section numbering: functions.R calls external validation Section 14, but analysis.R calls it Section 9 (9.2-9.5). The comments in Section 14 refer to "Section 9 of the driver script" and are correct. The preamble now says so. A reader going from `functions.R` Section 14 to analysis.R Section 14 would find nothing.
- Misplaced helpers: `line_colors`, `cluster_cols`, `umap_plot`, `plot_lines`, `legend_page`, `open_grid_pdf` sit in Section 13 (power analysis) though they are used throughout (`plot_palette_swatches` and `score_cell_cycle_by_cluster` call them). Consider a "shared plotting helpers" section.
- Out of chunk, noticed on the way: `cluster_marker_enrichment()` line 2548 reads `m[, 2]` by position (`gene_vec`). That is the avg_log2FC column in Seurat v5 today, but it violates the "select Seurat columns by name" rule. Suggest `m$avg_log2FC`. Section numbering also skips 8c-8f (8b is followed by 8g).
- `NOISE.CL` (analysis.R:2269) is created before Section 9.3 and stopped after; an error between them leaks the cluster. Wrap in `tryCatch(..., finally = parallel::stopCluster(NOISE.CL))` or use `on.exit` inside a helper.
- `fit_source` has no seed dependence; the fit is deterministic, so no RNG concerns. `&&`/`||` in the range: `refine_by_boundary` (`is.finite(bf) && bf <= max_boundary`) is scalar-safe; `plot_geneset_direction_stack` (`!is.na(n_sc[idx]) && n_sc[idx] > 0`) is scalar-safe; `fit_counts_offset` use is scalar. No violations found in this chunk.
- No Seurat column selection by position in this chunk.

## Functions reviewed with no issues

`collapse_to_orf`, `fread_matrix`, `get_data_header` (apart from edge notes), `read_header_line` (edge note), `qc_filter_counts` (validation note), `fit_source` (validation note), `mean_adjusted_noise`, `pilot_split_se`, `boot_contrasts`, `permute_contrasts`, `coexpr_bootstrap`, `assemble_coexpr_bootstrap`, `bootstrap_cluster_stability`, `bootstrap_compare_resolutions` (note 12), `chk`, `chk_prec`, `coexpr_raw_cor`, `coexpr_gene_degree`, `plot_coexpr_pair`, `report_cluster_marker_enrichment`, `plot_geneset_direction_stack` (comment only), `fit_split_nb`, `plot_palette_swatches`, `.cohen_kappa` (correct formula).
