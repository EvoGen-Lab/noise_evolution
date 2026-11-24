<!-- Exported from the Claude Docs document "Testing contrasts near a parameter boundary"
     (https://claude.ai/artifact/WZTnkkh5yM9XoYak4z19ec), revision 7, retrieved 2026-11-20. Display equations are LaTeX blocks.
     This file is a snapshot for the repository; the document remains the specification and is edited there.
     Implementation: sign_boot_p() and boot_contrasts_one() in R/functions.R (Section 2.3 of analysis/analysis.R); checks P1 (analysis/checks/sign_boot_test.R). -->

# Testing contrasts near a parameter boundary

Nov 19, 2026 · Brian Patrick Ha Metzger

## The problem

Many estimators live on a bounded parameter space, and the most common bound is zero. A variance component, an overdispersion such as the negative binomial dispersion φ, a heritability or a rate of extra variation cannot fall below zero, and when the true value is at or near zero the fitted value lands exactly on zero in a large share of datasets. Near the boundary these estimators behave badly. Their sampling distribution is a mixture of a point mass at zero and a continuous part, the usual normal approximations fail, the likelihood ratio statistic no longer follows χ² (Self and Liang 1987), and the ordinary bootstrap is inconsistent (Andrews 2000).

The trouble compounds when the question is a contrast between groups, such as whether one group is noisier than another. On a log scale a zero estimate has no value, so a log ratio involving it is ±∞. Bootstrap replicates that land on zero then sit at infinity, intervals lose one or both ends, and null distributions built by shifting or rescaling estimates cannot move a value that is already zero. The common responses are to gate out contrasts with too many zero replicates, to move each group onto a null value and rescale its replicates, or to switch to likelihood bounds near the boundary. Each adds machinery and choices that are hard to justify.

## The idea

The boundary is a problem for estimating a parameter, but not necessarily for deciding the direction of a difference between parameters. A bootstrap replicate in which group A refits to zero and group B to a positive value still says, unambiguously, that A is below B. On a log scale that replicate is −∞, which is useless as a magnitude but perfectly informative as a sign.

So test the contrast through its sign rather than its value. Resample the data, refit every group, compute the contrast in each replicate, and ask how often the replicates fall on the null side of zero. Zeros and infinities enter as ordinary values, because only their side of zero is used. This is the percentile bootstrap test applied to the contrast. Nothing has to be gated, projected or shifted, and for a comparison of two groups the test gives the same answer on any monotone scale of the parameter.

## Procedure

**1. Check the observed fits.** If the fit of any group in the contrast failed, report no contrast. If an observed estimate is zero, report no contrast on any scale where zero is infinite or undefined, such as the log scale or a reciprocal like 1/φ, because there is no value to report. Scales where zero is finite, such as the raw parameter, proceed.

**2. Bootstrap.** Resample the sampling units, cells for single-cell data, within each group, keeping paired units together. Refit every group with the same estimator and any corrections, and count a failed refit as zero, since fits usually fail where there is no signal above the null model.

**3. Compute the contrast in every replicate.** Write the contrast as a weighted sum over groups, Δ = Σ λ_r x_r. On a log or reciprocal scale, a zero in group r places the replicate at −∞ or +∞, on the side set by the sign of λ_r and of the transformation.

**4. Count the null side.** For an observed contrast above zero, the null side is at or below zero, and for one below zero it is at or above zero. A replicate exactly at zero, or undefined because it holds both +∞ and −∞, counts on the null side. With n₀ the number of the B replicates on the null side, the two-sided p-value is:

```latex
p = \min\left\{1,\ \frac{2\left(1 + n_0\right)}{1 + B}\right\}
```

**5. Report the interval.** Use the percentile interval of the replicates. On a scale where zeros are infinite it can be bounded on one side only. Report it that way, as a one-sided bound.

The add-one count keeps the p-value valid for a finite number of replicates and never zero (Phipson and Smyth 2010).

## Why it works, and what it gives up

**Direction survives the boundary.** When at least one group is clearly away from zero, the replicates of a group near zero, whether exactly zero or small and positive, almost always fall on the same side of the other group. How many of them are exactly zero, which is what the bootstrap gets wrong at the boundary, then barely affects the count of replicates on the null side. The infinities that break magnitude-based methods are, for this test, just very confident signs.

**Scale invariance for two groups.** For a two-group contrast the sign of the difference is the same on every monotone scale, so log, raw and reciprocal tests give one p-value. For contrasts of three or more groups, such as a difference of differences, the null hypothesis itself depends on the scale, and each scale needs its own test.

**Size is not protected.** The test decides whether and in which direction groups differ. It does not repair the estimate or interval of the size of a contrast near the boundary. On the log scale those can be infinite or one-sided, and on the raw scale the end of the interval nearest the boundary still inherits the bootstrap's trouble there.

**The weak case is every group near zero.** When the true values of all groups in a contrast are at or near zero, the null is often true and replicates are dominated by ties and zeros. The bootstrap's inconsistency at the boundary then matters most. Counting ties and undefined replicates on the null side pushes the test toward conservative there, but its calibration should be checked by simulation. A modest loss of power follows for a group at zero against one only slightly above it, a comparison the data can barely resolve anyway.

**Ordinary percentile-test limits remain.** The contrast is not studentized, so the accuracy is first order, the same as other unstudentized bootstrap tests.

## What it replaces

| Approach | What it does | Why the sign test is preferred |
| --- | --- | --- |
| Boundary gate | Withholds tests when more than a set share of replicates are at zero | Zeros are informative about direction, so the gate discards usable contrasts, and its threshold interacts awkwardly with contrasts that combine several groups |
| Shifted bootstrap null | Moves the bootstrap distribution onto the null, by subtracting the observed contrast or by rescaling each group to a null value | On a log scale it reduces to subtracting the observed contrast. Rescaling needs per-group null values chosen by convention, which change the answer on a raw scale, and replicates at zero cannot be moved. Recentering with absolute values also counts replicates far from the null as extreme, which makes tests needlessly conservative |
| Likelihood bounds | Compares one-sided profile-likelihood bounds of two groups when one is near zero | It is descriptive, lies outside multiple-testing control, needs a joint profile for contrasts of more than two groups, and relies on the likelihood model near the boundary |
| Pooled permutation | Reassigns units among groups | Under a variance null the groups can still differ in their means, so pooled units mix means and inflate the refitted variance |

Where a test with formal boundary theory is preferred, the established route is a likelihood ratio test with a null distribution from a parametric bootstrap under the fitted null model, which simulates a common zero exactly. Its cost is reliance on the model being correct.

## Checks before use

1. **Calibration where every group is near zero.** Simulate data whose groups all sit at the boundary, or share a small equal value, at the sample sizes and depths of the real data. Confirm that the false-positive rate is at or below nominal and that the p-value histogram is flat or piled toward one.
2. **Calibration away from the boundary.** Simulate equal values with differing nuisance parameters, such as different means for a dispersion contrast, and under departures from the fitted model, and confirm uniform p-values.
3. **Power at the edge.** Simulate one group at zero against groups slightly above it to measure what counting ties on the null side costs.
4. **Failed fits.** Check that failures occur where the signal is null-like, so that counting them as zero is reasonable for the data at hand.

## References

Andrews DWK (2000) Inconsistency of the bootstrap when a parameter is on the boundary of the parameter space. Econometrica 68:399–405. https://doi.org/10.1111/1468-0262.00114

Phipson B, Smyth GK (2010) Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology 9:39. https://doi.org/10.2202/1544-6115.1585

Self SG, Liang KY (1987) Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard conditions. Journal of the American Statistical Association 82:605–610. https://doi.org/10.1080/01621459.1987.10478472
