<!-- Draft appendix for "Statistical models for allele-specific expression, noise, and ploidy" (docs/methods.md, revision 18037).
     Written against analysis.R Section 2.5, cluster/scripts/gene_cond.R and R/functions.R (fit_dispersion(z = ), .pois_loglinear(), cond_contrasts_one()).
     Symbols were checked against the notation table and every section of methods.md: Z, β, ψ and the set 𝒫 are new and appear nowhere else in the
     document. W, w, s, t, a and T were avoided (Appendix A weights, the per-cell slope, the share multiplier, the expression factor). -->

## Appendix C. State-conditioned estimates

This appendix gives an estimate of every contrast that is conditioned on one cell state. It flags the significant contrasts that the state, rather than the genotype, carries.

**Why the culture matters.** We grew the two parents together in one culture and the hybrid in a second culture. Every contrast between parental and hybrid cells, which includes the *trans*, midparent and parental contrasts, therefore also carries any difference between the two cultures, and the design cannot separate that difference from the genotype. The total contrast compares the two parents within one culture, and the *cis* contrast compares two alleles within the same cells, so neither carries the culture difference. One culture difference is visible in the cells. A program of glucose derepression and of the induced stress response (Gasch et al. 2000) varies among cells. Its stress-response score, as a depth-aware z-score, exceeds 2.5 in 4.0% of *S. cerevisiae* parental cells and 5.9% of *S. eubayanus* parental cells, but in only 0.25% of hybrid cells, and its derepression score exceeds 4 in 3 of 834 and 6 of 1240 parental cells and in none of 2791 hybrid cells (check ST1, `analysis/checks/state_groups.R`). In simulations of a state present in 0.1% to 0.3% of cells with a 50-fold expression change, the type I error of the dispersion test at the 5% level averaged 0.12 for total and 0.08 to 0.10 for *trans*, and reached 0.155, because a group that samples no state cells cannot show them in its bootstrap. The *cis* test averaged at most 0.06 across all settings and reached 0.09, because both alleles of a cell share its state (check M3, 200 datasets per setting). We therefore kept the unconditioned estimates and tests as the main analysis, refitted every contrast with the program in the mean, and flagged the calls that the program carries.

**Program score.** Twenty marker genes, the set 𝒫, define the program. Ten mark glucose derepression (ICL1, MLS1, CIT2, MDH2, ACS1, FOX2, ADH2, IDP2, JEN1, ALD4) and ten the induced stress response (HSP12, HSP26, HSP42, HSP104, CTT1, TPS1, TPS2, GSY2, GRE2, SSA4). We chose them before any noise estimate, as the first members of ordered candidate lists with a one-to-one ortholog in both species. For gene g in cell i, the score is the log₂ of the summed marker counts per unit offset, with gene g's own count left out so that no gene enters its own covariate. In a hybrid cell the sum runs over both alleles, so both alleles of a cell share one score. We centered the score at its median within each culture, separately for the two parents, because the marker genes differ in expression between the species:

$$
Z_i = \log_2\frac{\tfrac12 + \sum_{g' \in \mathcal{P},\, g' \neq g} y_{g'i}}{O_i} - \mathrm{median}_{i'}\left(\log_2\frac{\tfrac12 + \sum_{g' \in \mathcal{P},\, g' \neq g} y_{g'i'}}{O_{i'}}\right)
$$

The median runs over the cells i′ of the same parent, or over all hybrid cells. The pseudocount ½ keeps the score finite in a cell without marker counts. A hybrid allele's offset is half the cell's, and this constant factor cancels in the centering. Z\_i = 0 marks the typical cell of its culture, and each unit of Z\_i doubles the program's expression. We computed Z\_i once from the observed counts and held it fixed.

**Conditioned mean.** The mean count of a cell now depends on its score through a log-linear term with slope β\_r, the log₂ change in the gene's mean for each log₂ unit of the program:

$$
\mathrm{E}\left[Y_{ri} \mid Z_i\right] = M_i = O_i\,\mu_r\,2^{\beta_r Z_i}, \qquad \mathrm{Var}\left[Y_{ri} \mid Z_i\right] = M_i + \varphi_r\,M_i^2
$$

Here μ\_r is the mean rate of a cell at Z\_i = 0, the typical state of its culture, and φ\_r is the dispersion around the conditioned mean. We estimated μ\_r and β\_r from the two Poisson score equations, which extend the estimating equation of μ̂\_r in *Estimating μ* by one equation for β\_r:

$$
\begin{aligned}
\sum_i \left(y_{ri} - \hat{M}_i\right) &= 0, \qquad \sum_i Z_i\left(y_{ri} - \hat{M}_i\right) = 0, \qquad \hat{M}_i = O_i\,\hat{\mu}_r\,2^{\hat{\beta}_r Z_i} \\
\hat{\beta}_r = 0 \quad &\Rightarrow\quad \sum_i \left(y_{ri} - O_i\,\hat{\mu}_r\right) = 0 \quad\Rightarrow\quad \hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i}
\end{aligned}
$$

These equations assume only the form of the mean, so they are the quasi-Poisson estimating equations (Wedderburn 1974, McCullagh and Nelder 1989), and their solution contains no φ. Without the covariate they give the ratio estimate of the main text. We solved them by Newton steps from the ratio estimate and β\_r = 0, halving each step that lowered the Poisson log-likelihood. The corrected allele means, the per-genome hybrid mean and the transcriptome share corrections then follow from μ̂\_r as in the main text. In 300 simulated groups of 500 cells with β\_r = 0.8 and φ = 0.5, β̂\_r averaged 0.798 (SD 0.080).

**Conditioned dispersion.** The mean now uses two parameters, so the Cox–Reid adjustment uses the information matrix of both (Cox and Reid 1987, McCarthy et al. 2012). We write the three weighted sums of powers of the score as ψ\_0, ψ\_1 and ψ\_2. Up to a constant factor that does not move the maximum, the determinant of the information matrix of log μ\_r and β\_r is ψ\_0ψ\_2 − ψ\_1²:

$$
\begin{aligned}
\ell^{\mathrm{CR}}_r(\varphi) &\equiv \ell_r(\varphi) - \tfrac12\log\left[\psi_0(\varphi)\,\psi_2(\varphi) - \psi_1(\varphi)^2\right], \qquad \psi_n(\varphi) = \sum_i \frac{Z_i^{\,n}\,\hat{M}_i}{1 + \varphi\,\hat{M}_i}, \quad n = 0, 1, 2 \\
\hat{\beta}_r = 0 \ \text{fixed} \quad &\Rightarrow\quad \ell^{\mathrm{CR}}_r(\varphi) = \ell_r(\varphi) - \tfrac12\log\psi_0(\varphi) = \ell_r(\varphi) - \tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}
\end{aligned}
$$

where ℓ\_r(φ) is the NB log-likelihood with every cell's mean fixed at M̂\_i. With the slope fixed at zero, the adjustment returns the one-parameter form of *Estimating φ*. A score with no variance among the cells of a group carries no information, so we fitted that group without it.

Whether the maximum lies at zero again depends on the adjusted slope at zero. The slope of ℓ\_r at zero is the sum of the per-cell slopes of *Estimating φ*, and the first score equation sets Σ\_i y\_ri equal to Σ\_i M̂\_i. Each weighted sum ψ\_n falls with φ at the rate Σ\_i Z\_i^n M̂\_i² at zero, and ψ\_n(0) = Σ\_i Z\_i^n M̂\_i:

$$
\begin{aligned}
\left.\frac{\partial \ell_r}{\partial \varphi}\right|_{\varphi = 0} &= \tfrac12\left[\sum_i \left(y_{ri} - \hat{M}_i\right)^2 - \sum_i y_{ri}\right] = \tfrac12\left[\sum_i \left(y_{ri} - \hat{M}_i\right)^2 - \sum_i \hat{M}_i\right] \\
\left.\frac{\partial \psi_n}{\partial \varphi}\right|_{\varphi = 0} &= -\sum_i Z_i^{\,n}\,\hat{M}_i^2 \\
\left.\frac{\partial}{\partial \varphi}\left(-\tfrac12\log\left[\psi_0\psi_2 - \psi_1^2\right]\right)\right|_{\varphi = 0} &= \frac{\psi_2(0)\sum_i \hat{M}_i^2 + \psi_0(0)\sum_i Z_i^2\hat{M}_i^2 - 2\,\psi_1(0)\sum_i Z_i\hat{M}_i^2}{2\left[\psi_0(0)\,\psi_2(0) - \psi_1(0)^2\right]}
\end{aligned}
$$

The estimator therefore keeps the two cases of the main text, with the adjustment's slope in place of Σ\_i M̂\_i²/Σ\_i M̂\_i:

$$
\hat{\varphi}_r = \begin{cases} 0 & \sum_i \left(y_{ri} - \hat{M}_i\right)^2 \le \sum_i \hat{M}_i - \dfrac{\psi_2(0)\sum_i \hat{M}_i^2 + \psi_0(0)\sum_i Z_i^2\hat{M}_i^2 - 2\,\psi_1(0)\sum_i Z_i\hat{M}_i^2}{\psi_0(0)\,\psi_2(0) - \psi_1(0)^2} \\ \arg\max_{\varphi > 0}\ \ell^{\mathrm{CR}}_r(\varphi) & \text{otherwise} \end{cases}
$$

We applied the transcriptome share corrections of the hybrid alleles to the conditioned dispersions unchanged, which assumes that the transcriptome share of a cell is independent of its program score. In simulations with a gamma expression factor and a log-linear score effect (600 cells, φ = 0.3, 200 replicates), the conditioned estimate averaged 0.300, while the unconditioned estimate averaged 0.82 because it absorbed the program's variance. With Poisson counts and a score effect (300 replicates), 48% of conditioned estimates fell on the boundary, close to the 50% expected for a parameter on the boundary.

**Contrasts, bootstrap and flags.** We built every contrast of every quantity from the conditioned fits with the weights of *Contrasts*, so the log₂ and raw dispersion, burst size and kinetic balance contrasts follow from μ̂\_r and φ̂\_r exactly as before. Each bootstrap replicate used the resamples of the reported bootstrap and refitted every group twice, with and without the score, with the transcriptome share terms refit on the resampled hybrid cells. A resampled cell kept its score. This gave each contrast a conditioned estimate Δ^Z with a 95% percentile interval, and the paired change Δ − Δ^Z between the unconditioned and conditioned contrasts on the same cells, with its own percentile interval. We built no null for the conditioned contrasts, so they carry intervals but no p-values. We flagged a contrast with a significant unconditioned test (BH-adjusted q below 0.05) when its conditioned interval contains zero and the interval of its change excludes zero. For the *trans*, midparent and parental contrasts a flag marks an effect confounded with the culture. For the total and *cis* contrasts, whose groups share a culture, it marks an effect of the program itself. The figures draw flagged genes as open triangles.

**What conditioning removes.** Conditioning removes the program's contribution to each group's dispersion, and it evaluates each group's mean at the typical state of its culture, so the share of cells in the program no longer moves the mean contrasts. It does not remove a difference between the cultures in their typical program level, which the centering assigns to each culture's own median. It also leaves every culture difference that the score does not measure. The score is measured with counting noise, so the conditioning is partial and the conditioned dispersion stays above the dispersion given the true state. The marker genes keep their tests, and their conditioned dispersion measures the noise beyond what the other markers explain.

**Rows for Model assumptions.**

| Assumption | Statement | Used in |
| --- | --- | --- |
| Separate cultures | Parental and hybrid cells come from different cultures, so the *trans*, midparent and parental contrasts carry any culture difference. Conditioning on the program score removes only the measured part | Appendix C |
| Log-linear state effect | A gene's mean count changes with the program score as 2^(β\_r Z\_i), with the dispersion constant around that mean | Appendix C |
| Score independent of share | The transcriptome share of a hybrid cell is independent of its program score, so the share corrections apply to the conditioned dispersions unchanged | Appendix C |

**Rows for Notation.**

| Symbol | Definition |
| --- | --- |
| 𝒫 | Set of the program's marker genes |
| Z | Program score of a cell, the log₂ marker count per unit offset centered at its culture's median |
| β | Log₂ change in a gene's mean per log₂ unit of the program score |
| ψ | Weighted sum of powers of the program score in the Cox–Reid adjustment (Appendix C) |

**Pointer sentence for Contrasts.** *Appendix C* refits every contrast with a cell-state program in the mean and flags the significant contrasts that the program carries.

**Reference to add.** Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown PO (2000) Genomic expression programs in the response of yeast cells to environmental changes. Molecular Biology of the Cell 11:4241–4257. https://doi.org/10.1091/mbc.11.12.4241 (verify with `checks_text.R` G2; Crossref was not reachable when this draft was written).
