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

1. **Appendix A estimator (resolved in commit 1ffbccb).** The code's weighted estimator used other weights and had no bias-correction factor. `fit_dispersion(method = "weighted")` now follows Appendix A. It uses weights 1/(1 + φM)², φ̂ = (1 + 1/Σy) φ̌ / (1 − h φ̌), and iteration from φ = 0 with the convergence status recorded. It matches the reference implementation `appendixA_phi()` to 1.6e-13 on 2000 random genes, and 2 of 2000 reached the 100-step limit.
2. **Upper likelihood bound (being added to the document).** The document says "upper 95% likelihood bound" and does not name the critical value. The code uses 2.71, the 0.90 quantile of χ²₁, which gives a one-sided 95% bound. A two-sided 95% bound would use 3.84. E1(b) reports the coverage of the code's bound.
3. **CV in the variance cap (being added to the document).** N Σ O²/(Σ O)² = 1 + CV²(O) holds exactly only with the population CV (divisor N). The document does not state the divisor. G3 shows the identity at rounding error with the population CV and not with the sample CV. The difference is of order 1/N, so it is negligible at the observed group sizes.
4. **Gene set for the quoted median dispersions.** The text quotes median φ of Hc 0.150 and He 0.152 without naming the genes. G1 reports both the class-1 median and the median of all fits.
5. **Handoff B4, allele weight π.** The allele weight π is read as the allele's mean share, so t = c/π with η = σ²_c/π². This reading is confirmed.
6. **Numerical precision of `.cr_loglik()` (not a document mismatch).** The lgamma form takes the difference of two terms of size log Γ(1/φ). At φ = 1e-8 its finite-difference slope is dominated by rounding (A2). The fit searches φ ≥ 3.06e-7 and A2(b) checks the fitted maxima against a grid. The class LRT uses the same function.
7. **Burst size standard errors (the document now says the SEs include the covariance; resolved).** Previously the document said the burst size SEs treat Δ(μ) and Δ(φ) as uncorrelated. The code takes them from bootstrap replicates that resample the same cells for both partitions, so they include the covariance. The bsize and kbalsup nulls simulate the mean groups independently of the dispersion groups, so the nulls do treat them as uncorrelated. See the section on overlapping cell sets below.
8. **C6 call stability** (approach confirmed) uses a Wald approximation, |estimate / bootstrap SE| > 1.96, instead of rerunning the permutation tests for each of the 50 random partitions.

## Overlapping cell sets: the derivation in "Separate partitions for each quantity"

**Claim in the document.** Within a group, fit the mean on cells A and the dispersion on cells B. To first order the covariance of the two errors is a sum over the shared cells of Cov[Y_i, s_i]. That covariance is zero because Σ_y y ∂Pr(y)/∂φ = ∂M_i/∂φ = 0. The conclusion is that overlapping or identical cell sets do not correlate the mean and dispersion errors.

**Where it fails.** The identity Cov[Y_i, s_i] = ∂M_i/∂φ holds only when Pr(y) is the true distribution of the count. Here s_i is the slope of the NB log-likelihood, and Y_i follows the NB only when the expression factor is gamma. The document's own model gives a non-gamma factor to every hybrid allele: t·a with the share multiplier (Transcriptome share corrections), and d·q with shared and private parts (Appendix B). The Notes for revision raise the same point for the fits ("ML fit of non-gamma expression factors"). For any other factor the covariance is not zero. For the moment estimator the first-order term is exact and has a closed form. With equal offsets, φ̂ = (S² − Ȳ)/Ȳ² linearizes to δφ ≈ (δS² − δȲ)/μ² − 2(σ² − μ)δȲ/μ³. For a Poisson mixture with factor a (mean 1, variance φ), the third central moment of the count is κ₃ = μ + 3φμ² + μ³E[(a − 1)³]. Then

Cov[μ̂_A, φ̂_B] = |A ∩ B| / (N_A N_B) · [κ₃ − σ² − 2φμσ²] / μ² = |A ∩ B| / (N_A N_B) · μ (E[(a − 1)³] − 2φ²),

because σ² + 2φμσ² = μ + 3φμ² + 2φ²μ³ is the third central moment of the NB. The covariance is zero for disjoint sets, as the document says. For overlapping sets it is zero only when the factor has the gamma's third central moment 2φ². The correlation it implies does not shrink with N, because the covariance and the product of the SDs are both of order |A ∩ B|/(N_A N_B).

**Numerical confirmation** (check A4_theory in `checks_sim.R`; 16 cells of family × μ {1, 5} × φ {0.3, 1}, N = 2000, equal offsets, 4000 replicates per cell). The formula matches the simulated covariance to within about 2 MC SE in every cell. The table shows μ = 5, φ = 1:

| Factor | E[(a − 1)³] (gamma: 2) | Cov(μ̂, φ̂_mom), formula | simulated (MC SE) | cor, moment | cor, Cox–Reid |
|---|---|---|---|---|---|
| gamma | 1.98 | −4.0e-5 | 1.1e-4 (1.0e-4) | 0.017 | 0.015 |
| lognormal | 3.91 | 4.8e-3 | 5.1e-3 (2.8e-4) | 0.327 | 0.377 |
| bimodal | 0.94 | −2.7e-3 | −2.7e-3 (0.9e-4) | −0.541 | −0.362 |
| product of shared and private gamma (Var d = φ/4) | 2.63 | 1.6e-3 | 1.5e-3 (1.4e-4) | 0.177 | 0.172 |

The correlation grows with μ and φ: at μ = 1, φ = 0.3 every family lies within ±0.05. The product factor is the one the generator and `sim_alleles()` use. This explains the A4 simulation numbers (up to 0.13 within a group and 0.23 for the cis contrast on the same cells). The cross-allele term C7 stays near zero: cor(μ̂_Hc, φ̂_He) was 0.008 and 0.007 at shared variances of 0.25 and 0.75. The Cox–Reid estimator follows the same sign and size as the moment formula, but has no closed form here.

**What this affects.**
- **Bootstrap SEs and intervals:** not affected. They resample the observed cells, so they carry the data's covariance, including the shared cells of the mean, noise and kinetic partitions.
- **The document's statements:** not supported. These are the sentence on uncorrelated errors and the statement that burst size SEs treat Δ(μ) and Δ(φ) as uncorrelated (mismatch 7).
- **The burst size and κ nulls:** affected. They hold the fitted means fixed and simulate the dispersion groups, so within a null draw the mean and dispersion errors have the generator's third moment, not the data's. bsize and kbalsup also draw the mean groups independently of the dispersion groups. A positive covariance raises the variance of both burst size (Δμ + Δφ) and κ (−Δμ − 2Δφ), which makes the nulls too narrow, and a negative one makes them too wide. D4a (lognormal and bimodal truths) and D7 in the calibration array measure the resulting size. The data arms A4 and C7 (`checks_data.R`) measure the correlation in the data.

## Handoff section 9 (Nov 3 corrections): status

| ID | Status |
|---|---|
| H1 | Done: `share_terms()` sets c̄̂ = mean(ĉ_i), unweighted. |
| H2 | Done: `.mode_contrasts()` divides every allele mean by 2c̄̂ or 2(1 − c̄̂) (`.correct_means()`) before the mid and par contrasts average them. |
| H3 | Done: route A, geometric mean of the corrected allele dispersions on the log₂ scale and arithmetic mean on the raw scale; δ̂ and δ̂_raw are reported only (`APPENDIX.B`). |
| H4 | Done: contrasts are named mid and par. |
| H5 | Not done. This is a full-pipeline simulation, beyond the synthetic harness used to test the code. |
| H6 | Not done. The Oct 23 scripts are not in the repository or in this session. Checks A4_theory, A10, F2 and F3 cover the same ground with new code. |
| H7 | Done in code (commit 4d53e03): `STRUCT.MEAN.KBAL` in Section 3.6. It needs the data. |
| H8 | Not done. It needs a modeling choice (the dropout model) and the data. The suggesting reference stays uncited until its DOI is supplied. |
| H9 | Not done. It needs a modeling choice (clusters from Section 7.3, a shared φ across clusters) and the data. |
| H10 | Done: `PHI.MU.COR` in Section 2.6. |
| K2 addendum | Done: the class threshold is `qchisq(1 − 2 CLASS.ALPHA, 1)`, which is 2.7055 at α = 0.05. `cox_reid_check.R` asserts it. |
| K9 | Done. The burst size and κ bootstrap draws are sums of the mean and dispersion draws on the same resampled cells, so their SEs and intervals carry the covariance. `eiv_components()` no longer enters a missing correlation as zero (commit 4d53e03). The bsize and kbalsup simulation nulls still draw the mean groups independently (see the overlap section). |
| A4, A10 | Code: A4_theory and A10 (overlaps 0, 0.5, 1; equal and observed offsets; moment, Appendix A, unadjusted ML, Cox–Reid). |
| A5 | Code: A5 in `checks_data.R` gives each gene's realized ratio against its group's cap. Replacing the CV of 0.53 in the text waits for the numbers. |
| A9 | Code: Cox–Reid revisit with the untruncated moment estimator. My factor third-moment data check is renamed A11. |
| B4 addendum | Code: 20,000 cells with equal offsets, and the observed hybrid depths. |
| F2, F3 | Redone with new code (F2 replicated; see above). |
| F5, B5 | Code: F5_routes adds the bootstrap SD ratio of the two routes with Cox–Reid fits, next to ROUTE.CHECK's Spearman correlation. |

## Rerun note

Commit 1ffbccb changes `fit_dispersion(method = "weighted")` (Appendix A). It needs a rerun of Section 2.3 (`SPLIT.FITS.W` and the weighted columns of `GENE.TABLE`) and of the weighted columns of `BOOT.CONTRASTS` and `PR` with every table read from them, and resubmission of `gene_boot` and `gene_perm`. `SPLIT.FITS.W` can be refitted locally; the bootstrap and null columns need the cluster. The NB results are unchanged. The other check commits change no code in `R/functions.R`, `analysis/analysis.R` or the pipeline's cluster scripts, so no section or pipeline cluster job needs a rerun. The checks run as follows:

- `checks_sim.R` runs locally from the generator. It uses the observed design when the Section 2 checkpoint exists and stand-in offsets otherwise.
- The calibration array needs one new cluster job, `null_calibration.sub` (16 array tasks). Its inputs come from `checks_sim.R` block `CAL`, and block `CALMERGE` reads its outputs.
- `checks_data.R` needs the Section 1 and 2 checkpoints. Their cluster-dependent parts (C3, D9, E2, E3) read the `gene_boot` and `gene_perm` outputs already merged into the Section 2 checkpoint. A4, A6, C7 and D8 run a local bootstrap on a gene subsample instead of a cluster rerun.
- `checks_text.R` G2 needs network access to api.crossref.org.

## Open items

- **F2 statement for Appendix B (for the document).** The earlier "within 0.0005" could not be resolved with one run of 8e6 cells, because its MC SE was 0.0006 to 0.0013. The redone check gives a statement that holds. With inflated σ²_c (0.005 and 0.01), 20 replicates of 8e6 cells reproduce both star limits and recover σ²_d and σ²_q with no error beyond 1.1 MC SE (errors ≤ 1.8e-4, MC SE ≤ 2e-4). At the observed σ²_c the share terms (8.7e-5 and 1.8e-4) are below the Monte Carlo resolution and about 400 times smaller than the sampling SD of the estimators at the observed number of hybrid cells. So they cannot change any reported quantity.
- **A9, the third moment of the expression factor (new data check).** Factorial moments give each gene's E[(a − 1)³] − 2φ², which sets both the mean–dispersion covariance on shared cells and the NB bias under non-gamma factors (B4). Each tested gene is compared with NB data simulated at its own fitted mean and variance on the same offsets, so the null carries the estimator's finite-sample bias. The group statistic is the mean standardized deviation z over genes.
  - In simulation with 100 genes and 2000 cells, mean z is −0.07 to −0.02 for gamma factors (SE 0.10). It is 0.18 to 1.10 for lognormal, −0.10 to −0.85 for bimodal, and −0.07 to 0.34 for the shared-times-private gamma product.
  - One gene's test has little power: 13% of lognormal genes at φ 0.5 reach p < 0.05.
  - On the synthetic NB checkpoint, p < 0.05 in 0 to 6.8% of 44 genes per group, and mean z ranges from −0.19 to 0.03 (SE 0.15).
- **Cox–Reid conclusions (A1, B4):** revisit when the Section 2 checkpoint is available, by rerunning `checks_sim.R` on the observed design and `checks_data.R`. The current conclusions rest on stand-in offsets.

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
| B4 | (Rerun after commit 1ffbccb; the code's weighted estimator now equals Appendix A.) **outside** for Cox–Reid with non-gamma factors. Lognormal: median relative bias −6.3% (range −32% to +22%), median bias/SD −0.60. Bimodal: +4.2% (up to +37%), bias/SD 0.35. Gamma: −0.1%. The moment estimators (Appendix A and the code's) have median relative bias near 0 for all three families, max +29% at M̄ 0.5, φ 0.05. | per row |
| F1 sim | Appendix A (documented and code, identical after commit 1ffbccb): max \|bias\| 0.022 (relative 2.6%); \|bias\| > 2 MC SE in 4 of 48 cells; converged in every fit. | per row |
| A4 theory | The first-order covariance μ(E[(a − 1)³] − 2φ²)/N matches simulation within about 2 MC SE in 16 cells. The correlation of μ̂ and Cox–Reid φ̂ on the same cells ranges from −0.36 (bimodal) to 0.38 (lognormal) at μ 5, φ 1, and is −0.02 to 0.03 for gamma. | Fisher SE 0.016 |
| F2 | Redone with 20 replicates of 8e6 cells per setting, reading each error against its Monte Carlo SE. The star limits and the corrected estimators are confirmed: \|z\| < 1.1 in all 15 comparisons, errors at most 1.8e-4 with MC SE 1e-4 to 3e-4. At the inflated σ²_c (0.005, 0.01) the share terms the corrections remove (−0.022 to 0.098) are 100 to 900 MC SE. At the observed σ²_c (1.9e-5) they are 8.7e-5 (σ²*_d) and 1.8e-4 (σ²*_q), below the MC SE of 3e-4 and about 400 times below the sampling SD of the estimators at N_H = 2000 (0.07 to 0.08). The summed-factor variance matches its formula to within 1 MC SE (3e-5 to 5e-5). | 20 × 8e6 cells per setting |
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
