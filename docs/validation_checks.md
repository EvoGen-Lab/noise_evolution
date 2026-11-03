# Validation checks of the methods document

Covers `docs/methods.md` (revision 17945) and the validation handoff of Nov 3. The scripts are in `analysis/checks/`, their functions in `R/functions_checks.R`, and the calibration array in `cluster/scripts/null_calibration_grid.R` and `cluster/submit/null_calibration.sub`. Results are written to `results/checks/`, which git ignores. This file keeps the code and document mismatches, the status of the code changes K1 to K8, the rerun note and the latest numbers.

## Status of the code changes K1 to K8 (claude/refactor, commit e45e4e9)

| ID | Status on the branch | Where |
|---|---|---|
| K1 | Done. The Pearson pre-check, `PHI.FLOOR` and `THETA.CAP` are gone. The Cox–Reid adjusted fit uses the slope rule. Truncation at zero comes after the share correction, and `hit_range` records a maximum at either end of `LOG.THETA.RANGE`. | `fit_dispersion()`, `.cr_loglik()` |
| K2 | Done. Each fit has a class from the boundary LRT and the upper bound times M̄, assigned after the correction. There are gates per contrast and scale, gene sets per contrast, and the matched set `in_matched_set`. | `dispersion_class_stats()`, `GENE.TABLE`, `GATES` |
| K3 | Done. A zero φ maps to −∞ in replicates. Intervals and SEs come from percentiles, and genes above 2.5% zero replicates are flagged. | `.boot_summary()`, `BOOT.ZERO.FLAGGED` |
| K4 | Done. There is a third partition at f_κ. γ is pooled over genes that are class 1 in the parents and both alleles. `kbalsup` is the supplementary κ. | `SPLIT.FRAC$f_kbal`, `.KBAL_SOURCE` |
| K5 | Done. Each allele mean is Σy/ΣO. | `fit_dispersion()` |
| K6 | Done. The means are divided by 2c̄ and 2(1 − c̄). Appendix B variances are corrected before truncation. `share_sub` equals the document's summed-count term exactly (1000 random cases, 0 differences above 1e-12). | `.correct_means()`, `ploidy_components()` |
| K7 | Done. Each bootstrap replicate refits the share terms on its resampled hybrid cells. Every null refit uses the observed c̄ and η, and k̂ is held fixed. | `boot_contrasts_one()`, `null_contrasts_one()` |
| K8 | Done. Mean nulls are permutations. The midparent mean null and every dispersion and kinetic null are simulations with projected null values. | `null_contrasts_one()`, `gene_perm.R` |

## Code and document mismatches

1. **Appendix A estimator.** The document's moment estimator has weights 1/(1 + φM)², the corrected denominator B̃ and the closed form φ̂ = (1 + 1/Σy) φ̌ / (1 − h φ̌). The code's `fit_dispersion(method = "weighted")` differs. Its weights are 1/[(1 + φM)(1 + φ + φM(1 + 3φ))]. It solves the weighted moment equation Σ w (r² − A − φB) = 0 with `uniroot()`, and it has no (1 + 1/Σy)/(1 − hφ̌) factor. The checks implement the documented estimator (`appendixA_phi()`) and report both: F1, B4 and D10 (the D10 nulls use the code's estimator).
2. **Upper likelihood bound.** The document says "upper 95% likelihood bound" and does not name the critical value. The code uses 2.71, the 0.90 quantile of χ²₁, which gives a one-sided 95% bound. A two-sided 95% bound would use 3.84. E1(b) reports the coverage of the code's bound.
3. **CV in the variance cap.** N Σ O²/(Σ O)² = 1 + CV²(O) holds exactly only with the population CV (divisor N). The document does not state the divisor. G3 shows the identity at rounding error with the population CV and not with the sample CV. The difference is of order 1/N, so it is negligible at the observed group sizes.
4. **Gene set for the quoted median dispersions.** The text quotes median φ of Hc 0.150 and He 0.152 without naming the genes. G1 reports both the class-1 median and the median of all fits.
5. **Handoff B4, allele weight π.** For a single allele fit, the handoff's "allele weight π" is read as the allele's mean share, so t = c/π with η = σ²_c/π². This is an interpretation, not a statement of the document.
6. **Numerical precision of `.cr_loglik()` (not a document mismatch).** The lgamma form takes the difference of two terms of size log Γ(1/φ). At φ = 1e-8 its finite-difference slope is dominated by rounding (A2). The fit searches φ ≥ 3.06e-7 and A2(b) checks the fitted maxima against a grid. The class LRT uses the same function.
7. **C6 call stability** uses a Wald approximation, |estimate / bootstrap SE| > 1.96, instead of rerunning the permutation tests for each of the 50 random partitions.

## Rerun note

These commits change no code in `R/functions.R`, `analysis/analysis.R` or the pipeline's cluster scripts, so no section or pipeline cluster job needs a rerun. The checks run as follows:

- `checks_sim.R` runs locally from the generator. It uses the observed design when the Section 2 checkpoint exists and stand-in offsets otherwise.
- The calibration array needs one new cluster job, `null_calibration.sub` (16 array tasks). Its inputs come from `checks_sim.R` block `CAL`, and block `CALMERGE` reads its outputs.
- `checks_data.R` needs the Section 1 and 2 checkpoints. Their cluster-dependent parts (C3, D9, E2, E3) read the `gene_boot` and `gene_perm` outputs already merged into the Section 2 checkpoint. A4, A6, C7 and D8 run a local bootstrap on a gene subsample instead of a cluster rerun.
- `checks_text.R` G2 needs network access to api.crossref.org.

## Not testable with these data

These cannot be tested: species-specific ploidy effects (they need diploid parental strains), capture differences between parental and hybrid cells, and the bursty-limit reading of θ and ζ.

## Latest numbers (local run, Nov 3)

These come from the simulation checks at full size (`CHECK.SIZE = 1`). The Section 2 checkpoint is not in the repository, so they ran on stand-in offsets with the reported offset CVs (Pc 0.72, Pe 0.87, hybrid 0.44), 1000 cells per parent, 2000 hybrid cells, c̄ = 0.5024 and σ²_c = 1.9e-5. Rerun `checks_sim.R` once the checkpoint exists to use the observed design. `checks_data.R` and `checks_text.R` were tested on synthetic checkpoints only. They have not run on the data, and G2 could not reach Crossref from this container.

Rows marked **outside** have numbers outside the handoff's stated criterion.

| ID | Headline numbers | Monte Carlo error |
|---|---|---|
| A2 | (a) Closed-form slope at φ = 0 against the finite difference of the stable likelihood: 0 sign disagreements in 9873 datasets (median relative error 4.6e-7). (b) The grid finds no better interior maximum. In 19 datasets it finds a gain of at most 4.7e-5 at the upper end φ = 54.6, where the fit stopped at 54.59. The lgamma form `.cr_loglik()` differenced at φ = 1e-8 disagrees in sign in 4910 datasets (rounding; see mismatch 6). | exact counts |
| A1 | **outside**: the bias is not well below the SD in sparse cells. Largest Cox–Reid \|bias\|/SD is 0.61 (N 800, M̄ 0.2, φ 0.01), and \|bias\| exceeds 2 MC SE in 21 of 60 cells (unadjusted profile ML: 19 of 60). In the sparsest cells the Cox–Reid bias is no smaller than the unadjusted one, for example N 800, M̄ 0.05, φ 0.01: 0.383 against 0.357 (SD 0.64). Relative bias × N M̄ φ (cells with M̄φ ≥ 0.05) has median 3.07 for Cox–Reid and 1.43 for unadjusted, so the 1/(N M̄ φ) pattern is not removed. | per cell; 500 genes |
| A7 sim | Including the gene's own counts in the library biases φ by −0.13% at a library share of 1e-4, −2.1% at 0.01 and −9.5% at 0.05 (φ = 0.1); −0.05% to −8.1% at φ = 0.5. Read against the observed shares from `checks_data.R`. | below 2e-4 |
| E1 | (a) LRT size at φ = 0: median 0.044, range 0.029 to 0.059 over 36 cells. (b) Coverage of the upper bound: median 0.992, lowest 0.934 (N 800, M̄ 0.2, φ 0.5). (c) At φM̄ = 0.1, 5 to 6% of genes are called class 2 and 22 to 31% class 3. At φM̄ = 0.02, 11 to 13% are called class 1. | binomial SE 0.007 |
| B2 sim | Bias of σ̂²_c is within 2 MC SE in all 8 cells. SD at σ²_c = 1.9e-5 is 1.7e-6 (N 800) and 1.2e-6 (N 1600). | per row |
| B3 | Median \|bias\|/SD 0.033, max 0.61 (φ⁰ 0.005, M̄ 0.1). The correction truncates up to 7% of fits. Relative bias of the corrected means is at most 1.2%. | per row |
| B4 | **outside** for Cox–Reid with non-gamma factors. Lognormal: median relative bias −6.3% (range −32% to +22%), median bias/SD −0.60. Bimodal: +4.2% (up to +37%), bias/SD 0.35. Gamma: −0.1%. The moment estimators (Appendix A and the code's) have median relative bias near 0 for all three families, max +29% at M̄ 0.5, φ 0.05. | per row |
| F1 sim | Appendix A as documented: max \|bias\| 0.022 (relative 2.6%); \|bias\| > 2 MC SE in 4 of 48 cells; converged in every fit. | per row |
| F2 | Star limits: error 0.0011 (σ²*_d) and 0.0012 (σ²*_q) against the target 0.0005. Corrected σ²_d and σ²_q: up to 0.0011 and 0.0012. These are within about 2 MC SE, and the MC SE (0.0006 to 0.0013 at 8e6 cells) exceeds the target, so this run cannot resolve 0.0005. Summed-factor variance against the formula: 6e-5 (target 0.0002). | MC SE 0.0006 to 0.0013 |
| F3 | δ̂ bias/SD up to 0.58; \|bias\| > 2 MC SE in 9 of 30 cells. Up to 47% of σ̂²_q are truncated at zero (ω 0.1, φ 0.1, M̄ 0.5). | 300 genes per cell |
| F4 | δ̂ error up to 0.014 (MC SE 0.005) at an intrinsic variance ratio of 4 (ω 0.75). The other ratios are within 2 MC SE. | per row |
| C4 | f estimated from the same cells: cis size 0.0505 (SE 0.0035), trans 0.043 (0.0032). f fixed at 0.4: cis 0.050 (0.0034), trans 0.041 (0.0031, KS p 0.003). | binomial SE |
| A4 sim | **outside**. Within a group on the same cells, cor(μ̂, φ̂) reaches 0.129 (M̄ 5, φ 1, ω 0.75) and \|z\| > 2 in 6 of 18 cells. cor(Δμ, Δφ) of the cis contrast on the same cells reaches 0.233 (M̄ 5, φ 1, ω 0.25), \|z\| > 2 in 7 of 18. The disjoint and half-overlapping cell sets are near zero (\|z\| > 2 in 0 and 1 of 18). The departures grow with φ and M̄, which bears on the derivation in Methods: Separate partitions. | Fisher SE 0.032 |
| C7 sim | cor(μ̂ of one allele, φ̂ of the other): \|cor\| ≤ 0.065, \|z\| > 2 in 1 of 36 cells. | Fisher SE 0.032 |
| E3 sim | Coverage of the 95% percentile intervals: log₂ φ 0.915 to 0.965, burst size 0.905 to 0.965, κ 0.910 to 0.965. | binomial SE 0.015 |
| D5 | (a) 10000 projections: 0 violations. (b) Raw total null, reused log₂ draws: size 0.033 to 0.060. Separate arithmetic centre: 0.000 to 0.067. | binomial SE 0.015 to 0.02 |
| D6 | 0 of 10000 cases violate μ_H = μ_MP or the preservation of the cis contrast. | exact |
| D1 local | Smoke test only, 20 genes per test (N.SIM 200, N.PERM 400): the pipeline's null generator runs end to end on the generator. Mean rejection at 0.05 by scale is 0.017 to 0.092; too few genes to judge size. The array (`null_calibration.sub`) gives the D1 to D4, D7 and D10 numbers. | binomial SE about 0.05 |
| G3 | The ratio lies outside [1, cap] in 0 of 10000 cases. The cap identity holds with the population CV (error 1.4e-14) and not with the sample CV. The quadratic residual is 4.5e-12, and `fstar_from_r()` equals the formula. The γ → 0 remainder is γ³/128, so the document's f ≈ ½ − γ/8 is correct to that order. The large-γ form is γf − 1 ≈ −1/γ. | exact |

Not run here, needing the data or the cluster: A3, A4 data, A5, A6, A7 data, A8, B1, B2 data, B5, B6, C1 to C3, C5 (`power_grid`), C6, C7 data, D1 to D4, D7, D8 to D10, E2, E3 data, E4, F1 (b, c), F5, G1, G2.
