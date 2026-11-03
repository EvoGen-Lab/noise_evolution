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

## Latest numbers

See below; filled from the local run of `checks_sim.R`.
