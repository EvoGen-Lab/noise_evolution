# Project context and settled decisions

Synthesized from 29 working sessions (March to August 2026). Code changes must preserve the scientific logic below.

## Central claim
Noise and mean expression diverge through distinct regulatory mechanisms between S. cerevisiae (YPS1000) and S. eubayanus. The hybrid design decomposes divergence into cis and trans components. Target journal is Nature Ecology and Evolution.

## Five questions the pipeline answers
1. How much does noise diverge relative to mean expression?
2. Does burst size or burst frequency drive noise changes?
3. What are the cis and trans contributions to noise divergence?
4. How much noise is intrinsic versus extrinsic?
5. What are the dominance and additivity patterns in the hybrid?

## Biological logic the code relies on
- Burst frequency changes move mean and noise in opposite directions. Burst size changes move them in the same direction. The direction of joint change therefore indicates the kinetic parameter.
- Cis effects come from allele-specific differences within the hybrid. Trans effects come from the difference between the parental ratio and the hybrid allelic ratio.
- Extrinsic noise produces correlated fluctuations between alleles in the same cell. Intrinsic noise does not.
- Hybrids are diploid and parents are haploid. This can inflate parental noise and contribute to the trans signal. Keep it stated as a caveat.

## Open analytical issue: power asymmetry between trans and cis
The trans contrast carries an extra variance term that scales with 1/N_parental (var_trans = var_cis + 4/N_p). Trans is therefore easier to detect than cis at finite depth, and the gap grows with the parent to hybrid depth ratio.
Depth-equalizing downsampling does not remove this for noise. Thinning a negative binomial compresses noise toward the Poisson floor unevenly across genes (Fano after thinning = 1 + p x (Fano before - 1)).
Any code that reports the trans-dominant noise result must include a power analysis, a simulation, or an explicit statement of the asymmetry and its direction.

## Open analytical question: RP and Ribi genes
Ribosomal protein and Ribi genes show higher mean and higher noise in S. eubayanus. This fits burst size driven divergence. A planned check is to run the burst decomposition on this gene set alone. If burst size dominates, one trans regulatory difference in this regulon links the mean, noise, and kinetic results.

## Planned figures
1. Noise versus mean divergence, regulatory categories (cis, trans, compensatory, reinforcing), intrinsic versus extrinsic.
2. Burst size versus burst frequency contributions.
3. Genome-wide cis and trans noise divergence with power analysis.
4. Dominance and additivity of noise differences.
5. GO enrichment by regulatory category, including RP and Ribi genes.

## What the analysis does not claim
- It does not claim selection drove noise divergence. A two species design cannot separate selection from mutational input.
- It does not claim the cis and trans test finds the same variants as QTL studies.
- It does not claim dedicated noise regulating trans factors. The trans signal describes regulatory architecture in hybrids.
- It does not claim noise evolves faster or slower than mean without stating the metric and scale. Use relative values, not absolute values, where comparing the two.

## Methodological notes
- Noise is estimated by fitting a negative binomial across cells. The Fano factor and the coefficient of variation both give consistent gene sets.
- Significance of noise differences comes from a permutation test developed in the power analysis.
- Noise is also computed with and without adjustment for mean expression. Keep both code paths.
