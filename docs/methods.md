<!-- Exported from the Claude Docs document "Statistical Methods, allele-specific offsets"
     (https://claude.ai/artifact/76UUgMeP9KaXjjcBGwMQa7), retrieved 2026-11-20: Methods tab revision 5638, Appendices tab revision 118, References tab revision 6.
     The Notes for revision tab (revision 50) stays in the document. Display equations are LaTeX blocks. This file is a snapshot for the repository;
     the document remains the specification and is edited there. It replaces the earlier snapshot of "Statistical Methods" (79JHPV3XCH2K8VS7RJZQHS, revision 18040). -->

# Statistical Methods

Nov 5, 2026 · @Brian Patrick Ha Metzger

## Notation

The experiment contains two samples. The first is a mixture of the two parental species, *S. cerevisiae* (Pc) and *S. eubayanus* (Pe). The second is the hybrid of these species, containing the *S. cerevisiae* allele (Hc), and the *S. eubayanus* allele (He). H is used for the hybrid without distinguishg alleles. Cells are indexed by i and cell groups are indexed by r. Each gene is fit separately and we omit gene indexes (g) unless needed for clarity. Lowercase letters are used for latent variables and realized observations. We use matching uppercase letters for the corresponding random variables, for fixed values calculated from the data, and for expectations built from parameters. Greek letters are used for model parameters and the quantities derived from them. Uppercase expectations are the expected values of observable counts, such as M and F, while other quantities derived from parameters, such as the split bias δ, are Greek. Functions, such as ℓ, Ψ and the per-cell slope ℓ′\_i, and fixed design constants, such as the contrast weights λ, take their conventional symbols, and the degree of dominance k follows the convention of quantitative genetics. A bar marks the mean of a variable over cells. A hat marks an estimate of a latent variable or parameter. A tilde marks an intermediate or comparison estimator and a star marks a quantity that carries error from the library split. Parameters of the summed allele counts of a hybrid cell carry the subscript sum.&#32;

| Symbol | Definition |
| --- | --- |
| y, Y | UMI count of one gene in one cell |
| L | Library size of a cell, and of each allele library (L\_Hc, L\_He) in a hybrid cell |
| D | Reference library depth |
| O | Offset of a cell, L/D for its own library |
| N | Number of cells in a group |
| μ | Mean rate, i.e. the expected count at the reference depth. |
| M | Expected count of a cell at its library depth&#32; |
| φ | Dispersion, the variability in counts beyond Poisson counting |
| F | Fano factor, the variance of a count over its mean |
| ℓ | Log-likelihood of a group |
| Ψ | Gamma term of the log probability of one cell |
| θ | Burst frequency, 1/φ |
| β | Burst size, μφ |
| κ | Kinetic balance, θ/β |
| Δ | Contrast between groups, Δ(log₂ x) on the log₂ scale and Δ(x) on the raw scale |
| ℓ′\_i | Slope in φ of the log probability of one cell |
| γ | Ratio of the parental to the hybrid contrast variance |
| χ | Fraction of hybrid cells assigned to the cis estimate |
| σ² | Variance of the variable in its subscript |
| s | Transcriptome share of the S. cerevisiae genome in a hybrid cell |
| ξ | Split error, the relative error of the observed transcriptome share |
| η | Split variance, the variance of the split error of a cell (η\_ri) |
| δ | Split bias, the rise in the fitted allele dispersion caused by the split error (δ\_r), which the split correction subtracts |
| z | Expression factor |
| α | Second moment of the expression factors |
| v | Extrinsic deviation, shared by both alleles of a cell |
| q | Intrinsic deviation, private to one allele |
| ω | Intrinsic fraction of the dispersion |
| π | Count share of Hc in the summed allele counts |
| Λ | Ploidy shift, the fall in dispersion from summing two alleles |
| λ | Weight of a group in a contrast |
| x | Generic quantity on the scale of a contrast, such as log₂ μ, log₂ φ or φ |
| Δ\_inh, Δ\_H,P | Allelic inheritance contrasts, each hybrid allele against its own parent, and parental comparisons, the per-genome hybrid against one parent |
| Δ\_dom | Dominance contrast, the per-genome hybrid minus the midparent, the average of the two parents on the scale of the contrast |
| k | Degree of dominance, k = 2Δ\_dom/Δ\_total, and of trans dominance, k\_trans = 2Δ\_dom/Δ\_trans |
| Ξ | Shared shift, the loess curve of Δ\_dom(log₂ φ) against log₂ μ\_H fitted on calibrator genes and subtracted from every dominance contrast |
| ρ | Correlation of the Pearson residuals of two genes across the cells of a group |
| B, n₀ | Number of bootstrap replicates, and the number of them on the null side of zero |
| N\_sel, N\_test | Numbers of genes selected by both tests and tested in a family, which set the level of the coverage-adjusted intervals |
| G | Set of the marker genes of the cell-state program&#32; |
| S | Program score of a cell, the log₂ marker count per unit offset centered at its culture's median&#32; |
| τ | Log₂ change in a gene's mean per log₂ unit of the program score |
| ℐ | Weighted sums of the zeroth, first and second powers of the program score, ℐ\_0, ℐ\_1 and ℐ\_2, the entries of the information matrix in the conditioned Cox–Reid adjustment |

## Negative binomial model of single-cell counts

**Sequencing depth offset.** Sequencing depth and capture efficiency vary among cells and groups. This causes two problems. First, differences in sequencing depth between groups shift their contrasts, leading to differences in the estimated means. Second, cell-to-cell spread in sequencing depth inflates estimates of dispersion, so genes appear noisier than they are. We solved both problems with cell-specific offsets. The offset of each cell is its library size, L, divided by a reference depth, D. We set D to the median library size across all cells, calculate each offset once, and then hold it fixed.&#32;

```latex
\begin{gathered}
L_{ri} = \sum_g y_{gri}, \qquad D = \mathrm{median}\left(L_{ri}\right) \quad O_{ri} = \frac{L_{ri}}{D}
\end{gathered}
```

**Negative binomial model.** We modeled the counts of each gene in each group as negative binomial (NB). The NB distribution arises from a model of bursty transcription in which a gene switches between an inactive state and a short-lived active state that releases several transcripts at once.&#32;

```latex
Y_{ri} \sim \mathrm{NB}\left(O_i\,\mu_r,\ \varphi_r\right) \\

\mathrm{E}[Y_{ri}] = O_i\,\mu_r \equiv M_i, \qquad \mathrm{Var}[Y_{ri}] = M_i + \varphi_r\,M_i^2
```

where the mean rate of group r is μ\_r, the expected count of a cell sequenced at the reference depth. The expected count of cell i at its own depth is then M\_i. The dispersion, φ, scales how far the variance exceeds the Poisson variance. Small dispersions indicate little excess noise, with zero dispersion the Poisson distribution. When comparing estimates to prior work, we use the squared coefficient of variation:

```latex
\mathrm{CV}^2_{ri} = \frac{\mathrm{Var}[Y_{ri}]}{\mathrm{E}[Y_{ri}]^2} = \frac{O_i\,\mu_r + \varphi_r\,(O_i\,\mu_r)^2}{(O_i\,\mu_r)^2} = \frac{1}{O_i\,\mu_r} + \varphi_r
```

This separates two sources of variance. The first term is the Poisson contribution. It falls as depth or expression rises because counting noise shrinks as more molecules are sampled. The second term does not depend on depth and is the excess variability that bursting adds.

Evaluated at the reference depth, O\_i = 1, these give the moments of a cell group:

```latex
\mathrm{E}_r = \mu_r, \qquad \mathrm{Var}_r = \mu_r + \varphi_r\,\mu_r^2, \qquad \mathrm{CV}^2_r = \frac{1}{\mu_r} + \varphi_r
```

**Estimating μ.** The NB likelihood of group r is the product over its cells of the NB probability of each count. Setting the derivative of its log with respect to μ to zero gives the maximum likelihood estimate μ̃\_r (Lawless 1987). When the offsets differ among cells, μ̃\_r has no closed form, so it is usually found iteratively and jointly with φ:

```latex
\begin{aligned}
\Pr\left(y_{ri} \mid O_i\mu, \varphi\right) &= \frac{\Gamma\left(y_{ri} + \varphi^{-1}\right)}{\Gamma\left(\varphi^{-1}\right)\,y_{ri}!}\left(\frac{1}{1 + \varphi O_i\mu}\right)^{\varphi^{-1}}\left(\frac{\varphi O_i\mu}{1 + \varphi O_i\mu}\right)^{y_{ri}} \\
\ell_r(\mu, \varphi) &= \sum_i \left[\log\Gamma\left(y_{ri} + \varphi^{-1}\right) - \log\Gamma\left(\varphi^{-1}\right) - \log\left(y_{ri}!\right) - \varphi^{-1}\log\left(1 + \varphi O_i\mu\right) + y_{ri}\log\frac{\varphi O_i\mu}{1 + \varphi O_i\mu}\right] \\
&= \sum_i \left[\Psi_{ri}(\varphi) + y_{ri}\log\left(O_i\mu\right) - \log\left(y_{ri}!\right) - \left(\varphi^{-1} + y_{ri}\right)\log\left(1 + \varphi O_i\mu\right)\right], \qquad \Psi_{ri}(\varphi) \equiv \log\frac{\Gamma\left(y_{ri} + \varphi^{-1}\right)}{\Gamma\left(\varphi^{-1}\right)} + y_{ri}\log\varphi = \sum_{j=0}^{y_{ri}-1}\log\left(1 + j\varphi\right) \\
\frac{\partial \ell_r}{\partial \mu} &= \frac{1}{\mu}\sum_i \frac{y_{ri} - O_i\mu}{1 + \varphi O_i\mu}= 0 \\
\tilde{\mu}_r &= \frac{\sum_i y_{ri}/\left(1 + \varphi O_i\tilde{\mu}_r\right)}{\sum_i O_i/\left(1 + \varphi O_i\tilde{\mu}_r\right)}
\end{aligned}
```

Here Ψ\_ri is the gamma term of cell i. It is the part of the likelihood that describes bursting and does not depend on the mean, so it does not enter the derivative in μ (Appendix A).&#32;

μ̃\_r can be rewritten as the sum of a φ independent and φ dependent terms:

```latex
\tilde{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i} - \frac{\varphi}{\sum_i O_i}\sum_i \frac{O_i\tilde{\mu}_r\left(y_{ri} - O_i\tilde{\mu}_r\right)}{1 + \varphi O_i\tilde{\mu}_r}
```

At the true mean, the φ dependent term has expectation zero:

```latex
\mathrm{E}\left[\sum_i \frac{M_i\left(Y_{ri} - M_i\right)}{1 + \varphi M_i}\right] = \sum_i \frac{M_i\left(\mathrm{E}[Y_{ri}] - M_i\right)}{1 + \varphi M_i} = 0
```

This term therefore adds no systematic shift to μ̃\_r. Instead, it changes the relative contribution of cells when estimating the mean. Dropping this term gives a closed form estimator, μ̂\_r, that is independent of φ.

```latex
\hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i}
```

This is the Poisson score with a log link and an offset. It assumes only that the mean count is O\_i μ, which makes it the quasi-Poisson estimating equation (Wedderburn 1974, McCullagh and Nelder 1989). When φ = 0, every weight equals one and μ̃\_r and μ̂\_r are equal. Similarly, when every cell has the same offset, every weight is the same and cancels from the equation.&#32;

Using μ̂\_r over μ̃\_r has three advantages. First, it keeps the mean separate from the dispersion. In a joint fit, μ̃\_r depends on φ̂\_r, so sampling error in φ̂\_r reaches the mean and every mean contrast built from it. Under the NB, orthogonality of μ and φ makes this dependence second order in large samples, but μ̂\_r removes it at every sample size. Second, a joint fit using μ̃\_r weights each rate by O\_i for genes whose φ̂\_r is zero and by O\_i/(1 + φ̂\_r M\_i) for the rest, so the mean estimator would differ among genes according to whether the dispersion was zero. μ̂\_r treats every gene the same way. Third, μ̂\_r is computed once from the counts and offsets, before any dispersion fit, so the dispersion fit runs with every cell mean fixed rather than alternating between μ and φ. This produces a computational savings which is repeated for every gene in every observed, bootstrap, null, and power-simulation fit.

The cost of using  μ̂\_r over μ̃\_r is precision. For μ̂\_r, the sample variance follows from the NB definition. For μ̃\_r , the large-sample variance is the inverse of the Fisher information. The ratio of these variances is:

```latex
\begin{aligned}
\mathrm{Var}[\hat{\mu}_r] &= \frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i O_i\right)^2} = \mu_r^2\,\frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2} \\
\mathrm{Var}\left[\frac{\partial \ell_r}{\partial \mu}\right] &= \frac{1}{\mu_r^2}\sum_i \frac{\mathrm{Var}[Y_{ri}]}{\left(1 + \varphi\,M_i\right)^2} = \frac{1}{\mu_r^2}\sum_i \frac{M_i}{1 + \varphi\,M_i} \quad\Rightarrow\quad \mathrm{Var}[\tilde{\mu}_r] \approx \frac{\mu_r^2}{\sum_i M_i/\left(1 + \varphi\,M_i\right)} \\
\frac{\mathrm{Var}[\hat{\mu}_r]}{\mathrm{Var}[\tilde{\mu}_r]} &\approx \mu_r^2\,\frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2} \cdot \frac{\sum_i M_i/\left(1 + \varphi\,M_i\right)}{\mu_r^2} = \frac{\sum_i M_i\left(1 + \varphi\,M_i\right)\,\sum_i M_i/\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2}
\end{aligned}
```

By the Cauchy–Schwarz inequality, the numerator is always larger or equal to the denominator, so the ratio of variances is at least one and μ̂\_r is never more precise than μ̃\_r. The weights of μ̃\_r gives cells whose counts vary the most (those with large φM\_i) less weight relative to their offsets and therefore produces a more precise estimate.

This ratio does have an upper bound that depends only on the variance of the offsets within a group. In the numerator, both the quotient M\_i/(1 + φM\_i) and the product M\_i(1 + φM\_i) increase with M\_i, so across cells they are positively correlated. The product of their averages is then at most the average of their product:

```latex
\begin{aligned}
\frac{\mathrm{Var}[\hat{\mu}_r]}{\mathrm{Var}[\tilde{\mu}_r]} &= \frac{\frac{1}{N_r}\sum_i M_i\left(1 + \varphi\,M_i\right) \cdot \frac{1}{N_r}\sum_i \frac{M_i}{1 + \varphi\,M_i}}{\left(\frac{1}{N_r}\sum_i M_i\right)^2} \le \frac{\frac{1}{N_r}\sum_i M_i\left(1 + \varphi\,M_i\right)\frac{M_i}{1 + \varphi\,M_i}}{\left(\frac{1}{N_r}\sum_i M_i\right)^2} \\
&= \frac{N_r\sum_i M_i^2}{\left(\sum_i M_i\right)^2} = \frac{N_r\sum_i O_i^2}{\left(\sum_i O_i\right)^2} = \frac{N_r}{N_{\mathrm{eff}}} , \qquad N_{\mathrm{eff}} \equiv \frac{\left(\sum_i O_i\right)^2}{\sum_i O_i^2}
\end{aligned}
```

The cap, N\_r/N\_eff, depends only on the variability of the offsets, with the estimator μ̂\_r essentially acting like there are fewer cells than the estimator for μ̃\_r. The position of a gene between one and the cap depends on φM\_i. When φM\_i is small, as for lowly expressed genes or genes with φ near zero, the variances of the two estimators nearly coincide.

**Estimating φ.** Because μ̂\_r does not depend on φ, we estimated φ over φ ≥ 0 with each cell's mean fixed at its fitted value. Fitting the mean uses part of the information in the counts, so maximizing this log-likelihood underestimates φ. Since μ and φ are orthogonal under the NB model, we corrected for this by maximizing the Cox–Reid adjusted log-likelihood (Cox and Reid 1987, McCarthy et al. 2012).&#32;

```latex
\begin{aligned}
\ell_r(\varphi) &\equiv \ell_r\left(\hat{\mu}_r, \varphi\right) \\
\ell^{\mathrm{CR}}_r(\varphi) &\equiv \ell_r(\varphi) - \tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}, \qquad \hat{M}_i \equiv O_i\,\hat{\mu}_r
\end{aligned}
```

This adjustment removes the first-order bias that fitting μ adds to φ̂\_r. It was originally derived assuming that μ is refit by maximum likelihood at each φ. This assumption also holds with the μ̂\_r estimator and O\_i offsets (Appendix B).

The maximum of this likelihood has no closed form, so we located it numerically over log φ. Since the dispersion cannot be negative, we first determined if the estimate would be located at the zero boundary by determining the slope of the likelihood at zero. We take φ → 0 in the per-cell first derivative and sum over cells (Appendix A):

```latex
\begin{aligned}
\ell_r'(0) &= \lim_{\varphi \to 0}\sum_i\left[\Psi_{ri}'(\varphi) + \frac{\log\hat{F}_i}{\varphi^2} - \frac{\left(\varphi^{-1} + y_{ri}\right)\hat{M}_i}{\hat{F}_i}\right] = \sum_i\left[\tfrac12\,y_{ri}\left(y_{ri} - 1\right) + \tfrac12\,\hat{M}_i^2 - y_{ri}\,\hat{M}_i\right] \\
&= \tfrac12\sum_i\left[\left(y_{ri} - \hat{M}_i\right)^2 - y_{ri}\right] = \tfrac12\left[\sum_i \left(y_{ri} - \hat{M}_i\right)^2 - \sum_i \hat{M}_i\right]
\end{aligned}
```

The Cox-Reid adjustment adds its own slope at zero:

```latex
\left.\frac{\partial}{\partial \varphi}\left(-\tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}\right)\right|_{\varphi = 0} = \left.\tfrac12\,\frac{\sum_i \hat{M}_i^2/\left(1 + \varphi\,\hat{M}_i\right)^2}{\sum_i \hat{M}_i/\left(1 + \varphi\,\hat{M}_i\right)}\right|_{\varphi = 0} = \tfrac12\,\frac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i}
```

An adjusted slope that is not positive at the boundary makes φ = 0 a maximum. We therefore define the dispersion estimator as:

```latex
\hat{\varphi}_r = \begin{cases} 0 & \sum_i (y_{ri} - \hat{M}_i)^2 \le \sum_i \hat{M}_i - \dfrac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i} \\ \arg\max_{\varphi > 0}\ \ell^{\mathrm{CR}}_r(\varphi) & \text{otherwise} \end{cases}
```

We used this estimator for the observed fits and for every bootstrap, null, and power-simulation replicate. *Appendix C* provides a moment estimator for the dispersion that does not assume an NB distribution and uses only the first two moments of the counts. This estimator is less precise then the Cox-Reid fit, but it is unbiased when the NB assumption is violated (unlike Cox-Reid). We therefore calculate both estimators and flag any gene in which the moment estimator falls outside the 95% percentile interval of its Cox–Reid estimate.&#32;

**Burst kinetics.** In the bursty limit of this two-state model, the burst frequency is the inverse of the dispersion, θ = 1/φ, and burst size is then β = μφ (Raj et al. 2006). Burst frequency and burst size are therefore not independent.&#32;

```latex
\log_2\hat{\theta}_r = -\log_2\hat{\varphi}_r, \qquad \log_2\hat{\beta}_r = \log_2\hat{\mu}_r - \log_2\hat{\theta}_r = \log_2\hat{\mu}_r + \log_2\hat{\varphi}_r
```

Because β̂\_r is built from θ̂\_r, any sampling error in θ̂\_r enters β̂\_r with the opposite sign. Expanding their covariance and using the orthogonality of the NB mean and dispersion gives:

```latex
\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\beta}_r\right] = \mathrm{Cov}\left[\log_2\hat{\theta}_r,\ \log_2\hat{\mu}_r - \log_2\hat{\theta}_r\right] = \mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] - \mathrm{Var}\left[\log_2\hat{\theta}_r\right] \approx -\mathrm{Var}\left[\log_2\hat{\theta}_r\right]
```

This covariance is large when the dispersion is poorly estimated. For an individual gene, an imprecise dispersion estimate thus appears as burst frequency traded against burst size. Across genes, estimated burst frequency and burst size will be negatively correlated even when the underlying kinetics are unrelated.&#32;

To address this issue, we also summarized burst kinetics with the kinetic balance κ = θ/β. The pair μ̂\_r and κ̂\_r carries the same information as θ̂\_r and β̂\_r, but are less coupled.

```latex
\log_2\hat{\kappa}_r = \log_2\hat{\theta}_r - \log_2\hat{\beta}_r = 2\log_2\hat{\theta}_r - \log_2\hat{\mu}_r, \qquad \log_2\hat{\theta}_r = \tfrac12\left(\log_2\hat{\mu}_r + \log_2\hat{\kappa}_r\right), \qquad \log_2\hat{\beta}_r = \tfrac12\left(\log_2\hat{\mu}_r - \log_2\hat{\kappa}_r\right) \\
\mathrm{Cov}\left[\log_2\hat{\kappa}_r, \log_2\hat{\mu}_r\right] = 2\,\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] - \mathrm{Var}\left[\log_2\hat{\mu}_r\right] \approx -\mathrm{Var}\left[\log_2\hat{\mu}_r\right]
```

Because the mean is estimated far more precisely than the dispersion, kinetic balance is only weakly coupled to the mean, and typically much less than the coupling between θ̂\_r and β̂\_r. For example, with the mean and dispersion estimates uncorrelated, log₂ β̂\_r has variance Var\[log₂ μ̂\_r\] + Var\[log₂ θ̂\_r\]. Converting each covariance to a correlation so they are the same scale then gives:&#32;

```latex
\begin{aligned}
\rho\left(\log_2\hat{\theta}_r, \log_2\hat{\beta}_r\right) &\approx \frac{-\mathrm{Var}[\log_2\hat{\theta}_r]}{\sqrt{\mathrm{Var}[\log_2\hat{\theta}_r]\left[\mathrm{Var}[\log_2\hat{\mu}_r] + \mathrm{Var}[\log_2\hat{\theta}_r]\right]}} = -\left(1 + \frac{\mathrm{Var}[\log_2\hat{\mu}_r]}{\mathrm{Var}[\log_2\hat{\theta}_r]}\right)^{-1/2} \\
\rho\left(\log_2\hat{\kappa}_r, \log_2\hat{\mu}_r\right) &\approx \frac{-\mathrm{Var}[\log_2\hat{\mu}_r]}{\sqrt{\mathrm{Var}[\log_2\hat{\kappa}_r]\,\mathrm{Var}[\log_2\hat{\mu}_r]}} = -\frac{\mathrm{SE}(\log_2\hat{\mu}_r)}{\mathrm{SE}(\log_2\hat{\kappa}_r)}
\end{aligned}
```

When the dispersion is the less precise estimate, burst frequency and burst size are almost perfectly coupled, while kinetic balance and the mean are coupled only by the ratio of two standard errors, which is small. The cost is precision in estimating kinetic balance:&#32;

```latex
\mathrm{Var}\left[\log_2\hat{\kappa}_r\right] = 4\,\mathrm{Var}\left[\log_2\hat{\theta}_r\right] + \mathrm{Var}\left[\log_2\hat{\mu}_r\right] - 4\,\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] \approx 4\,\mathrm{Var}\left[\log_2\hat{\theta}_r\right]
```

Kinetic balance therefore carries about twice the standard error of burst frequency. It gives up precision to weaken the coupling among estimators.

**Fano factor.** The Fano factor of a count, its variance divided by its mean, is a common summary of bursting. Under the NB, the Fano factor is one plus the dispersion times the expected count, so it rises with the sequencing depth of each cell. We therefore do not use it as an estimator and instead retain it only for derivations.

```latex
F =\frac{\mathrm{Var}[Y_{ri}]}{\mathrm{E}[Y_{ri}]} = \frac{O_i\,\mu_r + \varphi_r\,(O_i\,\mu_r)^2}{O_i\,\mu_r} = 1 + \varphi_r\,O_i\,\mu_r = 1 + \varphi_r\,M_i
```

## Library split correction

A hybrid cell carries two genomes that share one library.  Each allele uses an allele offset, the size of its allele library divided by D, which correctly handles differences in sequencing depth among cells. However, the library split also carries binomial sampling noise, which the allele offsets pass to the dispersion estimates. We define the transcriptome share, s\_i, as the fraction of UMIs that come from the *S. cerevisiae* genome. It is a latent property of the cell and reflects the fraction of mRNA transcribed from each genome weighted by how efficiently each genome’s molecules are captured and assigned. Here we estimate and correct for this effect.

**Library split.** Given the hybrid library size (L\_Hi = L\_Hci + L\_Hei) and transcriptome share, the number of *S. cerevisiae* UMIs is binomial. We write ŝ\_i for the observed share and ξ\_Hci for its relative error, the fraction by which the observed share misses the true share. ξ\_Hci has a mean of zero with binomial variance η\_Hci (Appendix D):

```latex
\begin{gathered}
L_{\mathrm{Hc}i} \mid L_{\mathrm{H}i}, s_i \sim \mathrm{Binomial}\left(L_{\mathrm{H}i}, s_i\right), \qquad \hat{s}_i = \frac{L_{\mathrm{Hc}i}}{L_{\mathrm{H}i}}, \qquad \xi_{\mathrm{Hc}i} \equiv \frac{\hat{s}_i}{s_i} - 1 \\
\mathrm{E}[\xi_{\mathrm{Hc}i}] = 0, \qquad \mathrm{Var}[\xi_{\mathrm{Hc}i}] = \frac{\mathrm{Var}[\hat{s}_i]}{s_i^2} = \frac{s_i(1-s_i)}{L_{\mathrm{H}i}\,s_i^2} = \frac{1-s_i}{s_i\,L_{\mathrm{H}i}} \equiv \eta_{\mathrm{Hc}i}
\end{gathered}
```

Using this, we can rewrite the allele offset of Hc in cell i:

```latex
O_{\mathrm{Hc}i} = \frac{L_{\mathrm{Hc}i}}{D} = \frac{L_{\mathrm{H}i}\,\hat{s}_i}{D} = \frac{L_{\mathrm{H}i}\,s_i}{D}\left(1 + \xi_{\mathrm{Hc}i}\right)
```

A single gene contributes a negligible part of the library, so we take the split error as independent of the count of the gene being fitted.

**Mean.** To determine the impact of the library split on μ̂\_r, we expand its ratio to second order and simplify:&#32;

```latex
\begin{aligned}
\mathrm{E}\left[\hat{\mu}_{\mathrm{Hc}}\right] = \mathrm{E}\left[\frac{\sum_i Y_{\mathrm{Hc}i}}{\sum_i O_{\mathrm{Hc}i}}\right] &\approx \frac{\mathrm{E}\left[\sum_i Y_{\mathrm{Hc}i}\right]}{\mathrm{E}\left[\sum_i O_{\mathrm{Hc}i}\right]}\left(1 + \frac{\mathrm{Var}\left[\sum_i O_{\mathrm{Hc}i}\right]}{\mathrm{E}\left[\sum_i O_{\mathrm{Hc}i}\right]^2}\right) \\
&= \frac{\mu_{\mathrm{Hc}}\sum_i L_{\mathrm{H}i}\,s_i/D}{\sum_i L_{\mathrm{H}i}\,s_i/D}\left(1 + \frac{\sum_i L_{\mathrm{H}i}^2\,s_i^2\,\eta_{\mathrm{Hc}i}}{\left(\sum_i L_{\mathrm{H}i}\,s_i\right)^2}\right) \\ &= \mu_{\mathrm{Hc}}\left(1 + \frac{\sum_i L_{\mathrm{H}i}^2\,s_i^2\,\eta_{\mathrm{Hc}i}}{\left(\sum_i L_{\mathrm{H}i}\,s_i\right)^2}\right)
\end{aligned}
```

Since the split error is independent of the counts, the only bias comes from the variance of the total allele offset. Because η\_Hci is of order 1/L\_Hi and the sums run over N\_H cells, this bias is of order 1/(N\_H L\_H). Every group has 800 cells or more, so this term is small and we drop it. While this makes the fitted means per-genome rates comparable to the parental rates and removes any difference in the capture efficiency of the genomes, it also removes differences in genome-wide mRNA output.&#32;

**Dispersion.** The effect of the library split on the dispersion estimator depends on what the NB dispersion measures. The NB arises when a Poisson count is scaled by an expression factor that has mean one and a gamma distributed variance φ:

```latex

Y_{ri} \mid z_{ri} \sim \mathrm{Poisson}\left(O_i\,\mu_r\,z_{ri}\right), \qquad z_{ri} \sim \mathrm{Gamma}\left(\text{shape} = 1/\varphi_r,\ \text{rate} = 1/\varphi_r\right) \\
\mathrm{E}[z_{ri}] = 1, \qquad \mathrm{Var}[z_{ri}] = \varphi_r

```

Mixing the Poisson over the gamma factor, with θ = 1/φ, gives the NB with mean M\_i and dispersion φ (Greenwood and Yule 1920):

```latex
\Pr(y) = \frac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)}\left(\frac{\theta}{M_i+\theta}\right)^{\theta}\left(\frac{M_i}{M_i+\theta}\right)^{y}
```

An NB fit of the dispersion therefore estimates the variance of the mean-one factor that scales the Poisson mean. Since the offset of each hybrid allele carries the library split, the fit overestimates the allele dispersions. To determine the size of this effect, we determined how much the split error moves the fitted dispersion away from its true value. The fit dispersion is determined by where the log-likelihood slope is zero. At the true mean and dispersion, this slope averages to zero. Using a star for the fitted dispersion before any correction and expanding the slope to first order gives the error of one fit. The expected shift, δ\_Hc, is then approximately the expectation of the ratios:

```latex
\begin{aligned}
\ell'_{\mathrm{Hc}}\left(\hat{\varphi}^*_{\mathrm{Hc}}\right) \approx \ell'_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}}) + \left(\hat{\varphi}^*_{\mathrm{Hc}} - \varphi_{\mathrm{Hc}}\right)\ell''_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}}) = 0 \quad&\Longrightarrow\quad \hat{\varphi}^*_{\mathrm{Hc}} - \varphi_{\mathrm{Hc}} \approx -\frac{\ell'_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}})}{\ell''_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}})} \\
\delta_{\mathrm{Hc}} \equiv \mathrm{E}\left[\hat{\varphi}^*_{\mathrm{Hc}} - \varphi_{\mathrm{Hc}}\right] &\approx \frac{\mathrm{E}\left[\ell'_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}})\right]}{-\mathrm{E}\left[\ell''_{\mathrm{Hc}}(\varphi_{\mathrm{Hc}})\right]}
\end{aligned}
```

The first and second derivatives of ℓ\_Hc in φ are given in Appendix A. We evaluate these at the fitted means, M̂\_i and F̂\_i, and then expand to first order in η\_Hci. For the first derivative, both the slope and the linear split error term average to zero at the true means. Since ξ\_Hci and Y are independent, only the quadratic term, which is of first order in η, is left:

```latex
\begin{aligned}
\mathrm{E}\left[\ell'_{\mathrm{Hc}}(\varphi)\right] &= \sum_i \mathrm{E}\left[\Psi_{\mathrm{Hc}i}'(\varphi) + \frac{\log\hat{F}_i}{\varphi^2} - \frac{\left(\varphi^{-1}+Y_{\mathrm{Hc}i}\right)\hat{M}_i}{\hat{F}_i}\right], \qquad \hat{M}_i = M_i\left(1 + \xi_{\mathrm{Hc}i}\right), \quad \hat{F}_i = 1 + \varphi\hat{M}_i \\
&\approx \sum_i \mathrm{E}\left[\Psi_{\mathrm{Hc}i}'(\varphi) + \frac{\log F_i}{\varphi^2} - \frac{\left(\varphi^{-1}+Y_{\mathrm{Hc}i}\right)M_i}{F_i} + M_i\,\xi_{\mathrm{Hc}i}\,\frac{M_i - Y_{\mathrm{Hc}i}}{F_i^2} + \tfrac12 M_i^2\,\xi_{\mathrm{Hc}i}^2\left(\frac{1}{F_i^2} - \frac{2\varphi\left(M_i - Y_{\mathrm{Hc}i}\right)}{F_i^3}\right)\right] \\
&= \sum_i\left[\mathrm{E}\left[\Psi_{\mathrm{Hc}i}'(\varphi) + \frac{\log F_i}{\varphi^2} - \frac{\left(\varphi^{-1}+Y_{\mathrm{Hc}i}\right)M_i}{F_i}\right] + \frac{M_i\left(M_i\,\mathrm{E}\left[\xi_{\mathrm{Hc}i}\right] - \mathrm{E}\left[\xi_{\mathrm{Hc}i}Y_{\mathrm{Hc}i}\right]\right)}{F_i^2} + \tfrac12 M_i^2\left(\frac{\mathrm{E}\left[\xi_{\mathrm{Hc}i}^2\right]}{F_i^2} - \frac{2\varphi\left(M_i\,\mathrm{E}\left[\xi_{\mathrm{Hc}i}^2\right] - \mathrm{E}\left[\xi_{\mathrm{Hc}i}^2Y_{\mathrm{Hc}i}\right]\right)}{F_i^3}\right)\right] \\
&= \sum_i \frac{M_i^2\,\mathrm{E}\left[\xi_{\mathrm{Hc}i}^2\right]}{2F_i^2} \\
&= \sum_i \frac{\eta_{\mathrm{Hc}i}\,M_i^2}{2F_i^2}
\end{aligned}
```

We do the same for the second derivative. In this case, the first order term provides each cell's information about φ, and is typically not zero (appendix A). The higher order terms then simply adjust this effect and are of order η\_Hci², so we drop them.&#32;

```latex
\begin{aligned}
-\mathrm{E}\left[\ell''_{\mathrm{Hc}}(\varphi)\right] &= \sum_i \mathrm{E}\left[-\Psi_{\mathrm{Hc}i}''(\varphi) + \frac{2\log\hat{F}_i}{\varphi^3} - \frac{2\hat{M}_i}{\varphi^2\hat{F}_i} - \frac{\left(\varphi^{-1}+Y_{\mathrm{Hc}i}\right)\hat{M}_i^2}{\hat{F}_i^2}\right] \\
&\approx \sum_i \mathrm{E}\left[-\Psi_{\mathrm{Hc}i}''(\varphi) + \frac{2\log F_i}{\varphi^3} - \frac{2M_i}{\varphi^2 F_i} - \frac{\left(\varphi^{-1}+Y_{\mathrm{Hc}i}\right)M_i^2}{F_i^2} + M_i\,\xi_{\mathrm{Hc}i}\,\frac{2M_i\left(M_i - Y_{\mathrm{Hc}i}\right)}{F_i^3}\right] \\
&= \sum_i\left[-\mathrm{E}\left[\Psi_{\mathrm{Hc}i}''(\varphi)\right] + \frac{2\log F_i}{\varphi^3} - \frac{2M_i}{\varphi^2 F_i} - \frac{\left(\varphi^{-1}+M_i\right)M_i^2}{F_i^2} + \frac{2M_i^2\left(M_i\,\mathrm{E}\left[\xi_{\mathrm{Hc}i}\right] - \mathrm{E}\left[\xi_{\mathrm{Hc}i}Y_{\mathrm{Hc}i}\right]\right)}{F_i^3}\right] \\
&= \sum_i\left[-\mathrm{E}\left[\Psi_{\mathrm{Hc}i}''(\varphi)\right] + \frac{2\log F_i - F_i + F_i^{-1}}{\varphi^3}\right]
\end{aligned}
```

Dividing the expected first derivative by the expected second derivative gives δ\_Hc:

```latex
\begin{aligned}
\delta_{\mathrm{Hc}} &\approx \frac{\sum_i \eta_{\mathrm{Hc}i}\,M_i^2/\left(2F_i^2\right)}{\sum_i\left[-\mathrm{E}\left[\Psi_{\mathrm{Hc}i}''(\varphi)\right] + \left(2\log F_i - F_i + F_i^{-1}\right)/\varphi^3\right]}
\end{aligned}
```

The fitted dispersion maximizes the Cox–Reid adjusted likelihood rather than ℓ\_Hc itself. The slope of the adjustment does not grow with the number of cells, while the slope of ℓ\_Hc is a sum over the N\_H hybrid cells. The split error therefore acts on the adjustment about N\_H times less than on ℓ\_Hc, and we ignore its effect there. Replacing each quantity in δ\_Hc by its value from the data produces a corrected allele dispersion, φ̂:

```latex
\hat{\varphi}_{\mathrm{Hc}} = \hat{\varphi}^*_{\mathrm{Hc}} - \hat{\delta}_{\mathrm{Hc}}, \qquad \hat{\delta}_{\mathrm{Hc}} = \frac{\sum_i \hat{\eta}_{\mathrm{Hc}i}\,\hat{M}_i^2/\hat{F}_i^2}{2\sum_i\left[-\hat{\mathrm{E}}\left[\Psi_{\mathrm{Hc}i}''\left(\hat{\varphi}^*_{\mathrm{Hc}}\right)\right] + \left(2\log\hat{F}_i - \hat{F}_i + \hat{F}_i^{-1}\right)/\left(\hat{\varphi}^*_{\mathrm{Hc}}\right)^3\right]} \\ \hat{F}_i = 1 + \hat{\varphi}^*_{\mathrm{Hc}}\hat{M}_i, \qquad \hat{\eta}_{\mathrm{Hc}i} = \frac{1-\hat{s}_i}{\hat{s}_i\,L_{\mathrm{H}i}}
```

Replacing s\_i by ŝ\_i changes η\_Hci by an amount of order 1/L\_Hi, which we ignored. In the denominator, −E\[Ψ\_Hci″(φ)\] depends on a cell only through M̂\_i and φ̂\*. We therefore compute it once on a grid and interpolate it for each cell. We corrected every hybrid allele dispersion in the same way and then truncated negative values at zero.&#32;

Each bootstrap replicate refitted the allele means and dispersions, with their split corrections, on its resampled hybrid cells, so the standard errors of every quantity built from them include the variation of the corrections.

## Contrasts

Significant differences between cell groups were determined using contrasts. For the mean, we used a log₂ scale.&#32;

```latex
\Delta(\log_2\mu) = \log_2\mu_1 - \log_2\mu_2
```

For dispersion we used a log₂ scale for relative differences and a raw scale for absolute differences.&#32;

```latex
\Delta(\log_2\varphi) = \log_2\varphi_1 - \log_2\varphi_2, \qquad \Delta(\varphi) = \varphi_1 - \varphi_2
```

Log₂ dispersion contrasts compare genes on a common footing so a doubling of φ counts the same whether the gene starts quiet or noisy. The log₂ scale therefore suits questions about the mechanism behind a difference between groups. In addition, these contrasts can be compared directly with log₂ mean contrasts.

On the raw scale, a contrast is the difference of two dispersions. A raw difference in φ is the difference in cell-to-cell variability beyond Poisson sampling, so it measures how much noise a difference adds. The raw scale suits questions about the size of noise and its likely consequences for phenotype. In addition, the raw scale can be used for genes in which the estimated dispersion is zero.

We use both the log₂ and raw scales for contrasts of burst size, burst frequency, and kinetic balance:

```latex
\begin{aligned}
\Delta(\log_2\beta) &= \log_2\beta_1 - \log_2\beta_2, &\qquad \Delta(\beta) &= \beta_1 - \beta_2 \\
\Delta(\log_2\theta) &= \log_2\theta_1 - \log_2\theta_2, &\qquad \Delta(\theta) &= \theta_1 - \theta_2 \\
\Delta(\log_2\kappa) &= \log_2\kappa_1 - \log_2\kappa_2, &\qquad \Delta(\kappa) &= \kappa_1 - \kappa_2
\end{aligned}
```

**Informative estimates.** A dispersion fit can fail, particularly for lowly expressed genes near the Poisson limit. In such cases, the contrast receives no estimate, interval, or p-value on either scale. The dispersion fit can also return an estimate of zero (the Poisson limit). In such cases, every log₂ contrast of the dispersion, burst frequency, burst size and kinetic balance are infinite. Similarly, the burst frequency and the kinetic balance are infinite on the raw scale when φ is zero. In such cases, the contrast, interval, and p-value are not calculated.&#32;&#32;

**Contrasts between groups.** For each quantity a gene has four values, the two parents and the two hybrid alleles. Contrasts see only differences among them, so three independent contrasts describe them completely, and we used the three that each correspond to a regulatory process. The *cis* contrast compares the two hybrid alleles, which share one *trans* environment. The *trans* contrast is the part of the parental difference that the *cis* contrast leaves unexplained. The dominance contrast Δ\_dom compares the per-genome hybrid with the midparent. Each allele carries its own *cis* sequence in both its parent and the hybrid, so on the log₂ scale the *cis* effects cancel in Δ\_dom, and Δ\_dom measures how the two parental *trans* environments combine in the hybrid (*Dominance*, below).

The per-genome hybrid and the midparent are averages on the scale of each contrast. Expression is measured relative to library size, and each hybrid allele relative to the library of its own genome, so each rate is known only up to a scale factor of its library. Under such factors the average of the log₂ allele rates shifts by the same constant in every gene, which the normalization absorbs, whereas the log₂ of the arithmetic average shifts by an amount that depends on each gene's allele ratio. We therefore averaged the log₂ values, so that log₂ μ\_H = ½(log₂ μ\_Hc + log₂ μ\_He), and likewise for the midparent. For the mean:

```latex
\begin{aligned}
\Delta_{\mathrm{cis}}(\log_2 \mu) &= \log_2\mu_{\mathrm{Hc}} - \log_2\mu_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\log_2 \mu) &= \left(\log_2\mu_{\mathrm{Pc}} - \log_2\mu_{\mathrm{Pe}}\right) - \left(\log_2\mu_{\mathrm{Hc}} - \log_2\mu_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\log_2 \mu) &= \tfrac12\left(\log_2\mu_{\mathrm{Hc}} + \log_2\mu_{\mathrm{He}}\right) - \tfrac12\left(\log_2\mu_{\mathrm{Pc}} + \log_2\mu_{\mathrm{Pe}}\right)
\end{aligned}
```

The dispersion contrasts take the same form. For the log₂ contrasts:

```latex
\begin{aligned}
\Delta_{\mathrm{cis}}(\log_2 \varphi) &= \log_2\varphi_{\mathrm{Hc}} - \log_2\varphi_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\log_2 \varphi) &= \left(\log_2\varphi_{\mathrm{Pc}} - \log_2\varphi_{\mathrm{Pe}}\right) - \left(\log_2\varphi_{\mathrm{Hc}} - \log_2\varphi_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\log_2 \varphi) &= \tfrac12\left(\log_2\varphi_{\mathrm{Hc}} + \log_2\varphi_{\mathrm{He}}\right) - \tfrac12\left(\log_2\varphi_{\mathrm{Pc}} + \log_2\varphi_{\mathrm{Pe}}\right)
\end{aligned}
```

&#32;And the raw contrasts:

```latex
\begin{aligned}
\Delta_{\mathrm{cis}}(\varphi) &= \varphi_{\mathrm{Hc}} - \varphi_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\varphi) &= \left(\varphi_{\mathrm{Pc}} - \varphi_{\mathrm{Pe}}\right) - \left(\varphi_{\mathrm{Hc}} - \varphi_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\varphi) &= \tfrac12\left(\varphi_{\mathrm{Hc}} + \varphi_{\mathrm{He}}\right) - \tfrac12\left(\varphi_{\mathrm{Pc}} + \varphi_{\mathrm{Pe}}\right)
\end{aligned}
```

We computed θ = 1/φ, β = μφ and κ = θ/β for each group and each hybrid allele. On the log₂ scale each is a sum of log₂ μ and log₂ φ:

```latex
\begin{aligned}
\log_2 \theta &= -\log_2 \varphi \\
\log_2 \beta &= \log_2 \mu + \log_2 \varphi \\
\log_2 \kappa &= \log_2 \theta - \log_2 \beta = -\log_2 \mu - 2\log_2 \varphi
\end{aligned}
```

Every log₂ burst frequency, burst size and kinetic balance contrast therefore simplifies to a combination of the mean and dispersion contrasts:

```latex
\begin{aligned}
\Delta(\log_2 \theta) &= -\Delta(\log_2 \varphi) \\
\Delta(\log_2 \beta) &= \Delta(\log_2 \mu) - \Delta(\log_2 \theta) \\
\Delta(\log_2 \kappa) &= \Delta(\log_2 \theta) - \Delta(\log_2 \beta) \\
&= 2\,\Delta(\log_2 \theta) - \Delta(\log_2 \mu)
\end{aligned}
```

The raw contrasts are:

```latex
\begin{aligned}
\Delta_{\mathrm{cis}}(\theta) &= \theta_{\mathrm{Hc}} - \theta_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\theta) &= \left(\theta_{\mathrm{Pc}} - \theta_{\mathrm{Pe}}\right) - \left(\theta_{\mathrm{Hc}} - \theta_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\theta) &= \tfrac12\left(\theta_{\mathrm{Hc}} + \theta_{\mathrm{He}}\right) - \tfrac12\left(\theta_{\mathrm{Pc}} + \theta_{\mathrm{Pe}}\right) \\
\Delta_{\mathrm{cis}}(\beta) &= \beta_{\mathrm{Hc}} - \beta_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\beta) &= \left(\beta_{\mathrm{Pc}} - \beta_{\mathrm{Pe}}\right) - \left(\beta_{\mathrm{Hc}} - \beta_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\beta) &= \tfrac12\left(\beta_{\mathrm{Hc}} + \beta_{\mathrm{He}}\right) - \tfrac12\left(\beta_{\mathrm{Pc}} + \beta_{\mathrm{Pe}}\right) \\
\Delta_{\mathrm{cis}}(\kappa) &= \kappa_{\mathrm{Hc}} - \kappa_{\mathrm{He}} \\
\Delta_{\mathrm{trans}}(\kappa) &= \left(\kappa_{\mathrm{Pc}} - \kappa_{\mathrm{Pe}}\right) - \left(\kappa_{\mathrm{Hc}} - \kappa_{\mathrm{He}}\right) \\
\Delta_{\mathrm{dom}}(\kappa) &= \tfrac12\left(\kappa_{\mathrm{Hc}} + \kappa_{\mathrm{He}}\right) - \tfrac12\left(\kappa_{\mathrm{Pc}} + \kappa_{\mathrm{Pe}}\right)
\end{aligned}
```

**Derived comparisons.** Every other comparison among the four values is a combination of these three. For a quantity x on the scale of a contrast, with x\_H the average of the two hybrid alleles, the total contrast is the sum of the *cis* and *trans* contrasts. The parental comparisons set the per-genome hybrid against each parent. The allelic inheritance contrasts set each hybrid allele against its own parent, so each compares one *cis* sequence in two *trans* environments, that of its parent and that of the hybrid:

```latex
\begin{aligned}
\Delta_{\mathrm{total}}(x) &= x_{\mathrm{Pc}} - x_{\mathrm{Pe}} = \Delta_{\mathrm{cis}}(x) + \Delta_{\mathrm{trans}}(x) \\
\Delta_{\mathrm{H,Pc}}(x) &= x_{\mathrm{H}} - x_{\mathrm{Pc}} = \Delta_{\mathrm{dom}}(x) - \tfrac12\,\Delta_{\mathrm{total}}(x), \qquad \Delta_{\mathrm{H,Pe}}(x) = x_{\mathrm{H}} - x_{\mathrm{Pe}} = \Delta_{\mathrm{dom}}(x) + \tfrac12\,\Delta_{\mathrm{total}}(x) \\
\Delta_{\mathrm{inh,c}}(x) &= x_{\mathrm{Hc}} - x_{\mathrm{Pc}} = \Delta_{\mathrm{dom}}(x) - \tfrac12\,\Delta_{\mathrm{trans}}(x), \qquad \Delta_{\mathrm{inh,e}}(x) = x_{\mathrm{He}} - x_{\mathrm{Pe}} = \Delta_{\mathrm{dom}}(x) + \tfrac12\,\Delta_{\mathrm{trans}}(x) \\
k &= 2\,\Delta_{\mathrm{dom}}(x)/\Delta_{\mathrm{total}}(x), \qquad k_{\mathrm{trans}} = 2\,\Delta_{\mathrm{dom}}(x)/\Delta_{\mathrm{trans}}(x)
\end{aligned}
```

Two ratios place the hybrid between its parents. The degree of dominance k scales Δ\_dom by half the total contrast and describes the inheritance of expression. The degree of *trans* dominance k\_trans scales Δ\_dom by half the *trans* contrast and describes the inheritance of the *trans* environment. It equals 1 when the hybrid *trans* environment matches that of Pc, so that Δ\_inh,c = 0, −1 when it matches that of Pe, so that Δ\_inh,e = 0, and 0 when it lies at the log₂ midpoint of the two. The same identities hold on the raw scale, but the *cis* effects cancel in Δ\_dom only on the log₂ scale.

Each estimated contrast replaces every quantity by its corrected estimate, so it converges to its true contrast within 𝒪(1/N) because each estimator does.

**Kinetic balance.** For a contrast between independent groups, the variances of the two groups add and the covariances between groups are zero, so the covariances, correlations and precision cost derived in *Burst kinetics* carry over under the NB, with each variance replaced by its sum over the two groups. For the *cis* contrasts, whose two alleles share cells, we assumed that the covariance between the mean estimate of one allele and the dispersion estimate of the other is negligible, so the same results hold. A zero dispersion makes kinetic balance infinite, so its contrasts follow the rules for zero estimates in *Informative estimates*.

Together, the kinetic balance and mean contrasts separate two biological questions. Δ(log₂ μ) asks how much expression differs. Δ(log₂ κ) then asks how that difference divides between differences in burst frequency and burst size. Rewriting the relationship between μ, θ, and β, we have:

```latex
\Delta(\log_2 \mu) = \Delta(\log_2 \theta) + \Delta(\log_2 \beta)
```

Setting either Δ(log₂ θ) or Δ(log₂ β) to zero, or setting the two equal, and substituting produces:

```latex
\begin{aligned}
\text{frequency only, } \Delta(\log_2 \beta) = 0 &: \quad \Delta(\log_2 \theta) = \Delta(\log_2 \mu) \ \Rightarrow\ \Delta(\log_2 \kappa) = 2\,\Delta(\log_2 \mu) - \Delta(\log_2 \mu) = \Delta(\log_2 \mu) \\
\text{size only, } \Delta(\log_2 \theta) = 0 &: \quad \Delta(\log_2 \kappa) = 0 - \Delta(\log_2 \mu) = -\Delta(\log_2 \mu) \\
\text{equal, } \Delta(\log_2 \theta) = \Delta(\log_2 \beta) &: \quad \Delta(\log_2 \theta) = \tfrac{1}{2}\Delta(\log_2 \mu) \ \Rightarrow\ \Delta(\log_2 \kappa) = \Delta(\log_2 \mu) - \Delta(\log_2 \mu) = 0
\end{aligned}
```

A difference carried entirely by burst frequency therefore gives Δ(log₂ κ) = Δ(log₂ μ), one carried entirely by burst size gives Δ(log₂ κ) = −Δ(log₂ μ), and equal log₂ differences in both give Δ(log₂ κ) = 0.

**Tests.** To identify significant contrasts, we tested each contrast against its null hypothesis of no difference with the bootstrap distribution of the contrast itself, for the mean and for the dispersion, burst frequency, burst size and kinetic balance alike. Each bootstrap replicate was refitted with the same estimators and corrections as the observed data, and recomputed the split corrections δ̂\_r on its resampled hybrid cells.

Each replicate resampled the cells of each parent with replacement and resampled the undivided hybrid cells once, and each resampled hybrid cell kept the *cis*, *trans* and partition memberships of the cell it copied. The expression of each allele varies through two parts. An extrinsic part is shared by both alleles of a cell and comes from cell-wide sources such as cell size, cell-cycle state and the levels of *trans*-acting factors. An intrinsic part is private to each allele and comes from events at that copy, such as the bursting of its own promoter. Resampling a hybrid cell whole keeps the extrinsic noise that its two alleles share, which lowers the variance of the *cis* contrast. Each replicate refitted every group with the same estimators and corrections and recomputed every contrast on both scales. We write a contrast as a weighted sum over the groups r it uses, with weights λ\_r of ±1 or ±½, so that a group with a zero mean or dispersion places a log₂ contrast at +∞ or −∞ according to the sign of its weight. A refit that failed counted as zero, because fits fail mainly where counts are low and the noise is near the Poisson limit. These infinite replicates stay in the bootstrap distribution, so a percentile interval can be bounded on one side only.

For an observed contrast above zero, the replicates on the null side are those at or below zero, and for an observed contrast below zero, those at or above zero. A replicate exactly at zero, or undefined because it holds both +∞ and −∞, counts on the null side. With n₀ the number of the B replicates on the null side, the two-sided p-value is:

```latex
p = \min\left\{1,\ \frac{2\left(1 + n_0\right)}{1 + B}\right\}
```

The test uses only the side of zero on which each replicate falls, not its distance from zero. A replicate in which one group is at zero still shows which group is the noisier, so zero and failed fits carry information about the direction of a contrast even though they carry none about its size. For the two-group contrasts, total and *cis*, the log₂ and raw tests give the same p-value, because a ratio of two dispersions exceeds one exactly when their difference exceeds zero. On the log₂ scale burst frequency needs no separate test, because Δ(log₂ θ) = −Δ(log₂ φ) in every replicate. The raw burst frequency, burst size and kinetic balance contrasts were tested from their own replicates.

This is the percentile bootstrap test applied to the contrast. The contrast is not studentized, so the test is accurate to first order. The bootstrap is also not consistent for a parameter on the boundary (Andrews 2000). Because the test uses only the direction of a contrast, this matters little when at least one group is clearly away from zero, and it matters most when every group is near zero, as for a gene that is Poisson in every group. Ties and undefined replicates arise mainly there, and counting them on the null side makes the test conservative in that case. In simulated genes with φ = 0 or a mean count of 0.2, the average false-positive rate at the 0.05 level was 0.008, and no simulated condition exceeded 0.05 by more than two standard errors. Where φ was well estimated, zero and failed replicates were rare and the rules left power unchanged.

**Dominance.** We described the inheritance of each quantity by the dominance contrast and the total contrast. Δ\_dom is the dominance effect d of quantitative genetics, and half the total contrast is the additive effect a, so the degree of dominance is the familiar ratio d/a:

```latex
k = \frac{2\,\Delta_{\mathrm{dom}}(x)}{\Delta_{\mathrm{total}}(x)}
```

The degree of dominance is 0 when the hybrid lies at the midparent, 1 when it equals Pc and −1 when it equals Pe. Because the midparent is the average of the parents on the scale of each contrast, the same definition holds for the mean on the log₂ scale and for every noise quantity on both scales.

**Shared shift.** The hybrid alleles are noisier than the parents by a shift shared by most genes. The design cannot separate this shift into effects of ploidy, culture or a genome-wide property of the hybrid, so we reported the shift on its own and classified genes by their departure from it. The shift changes with expression, so we fitted it as a smooth curve, Ξ(log₂ μ\_H), of Δ\_dom(log₂ φ) against the observed per-genome hybrid mean, log₂ μ\_H, by loess with a span of 0.3. We fitted the curve on calibrator genes, those whose total and *trans* dispersion contrasts were not significant. In these genes the two parental *trans* environments act alike, so a departure of the hybrid from the midparent reflects the shared shift rather than *trans* dominance, and fitting on them keeps any real link between noise dominance and expression in the adjusted contrasts. We subtracted Ξ, evaluated at each gene's log₂ μ\_H, from its Δ\_dom(log₂ φ). The median of Δ\_dom(log₂ μ) is close to zero across the expression range, so the mean carries no shift, and the identities of the burst quantities give the shifts of the others, −Ξ for burst frequency, Ξ for burst size and −2Ξ for kinetic balance. The same Ξ is subtracted from the two parental comparisons, Δ\_H,Pc and Δ\_H,Pe, in the observed data and in every replicate, which keeps them consistent with the adjusted Δ\_dom and with the allelic inheritance contrasts through the identities of *Derived comparisons*. On the raw scale we divided each hybrid allele value by 2 raised to the shift of its quantity at its log₂ μ\_H, 2^Ξ for the dispersion, before forming every raw contrast that uses the hybrid, the *cis* and *trans* contrasts as well as Δ\_dom. The shift multiplies both alleles alike, so it cancels in every log₂ ratio, but it also scales the raw difference between the alleles. Without this step the raw *trans* contrast would compare an unshifted parental difference with an inflated hybrid difference. The calibrator set was fixed from the observed data and Ξ was refitted in every bootstrap replicate, so the tests and intervals of the adjusted Δ\_dom include its uncertainty. We report Ξ by its median within deciles of log₂ μ\_H.

We tested the total contrast and the adjusted Δ\_dom against zero, each with the bootstrap test above, and assigned each gene from the two outcomes:

| Total contrast | Dominance contrast (Δ\_dom) | Class |
| --- | --- | --- |
| Not significant | Not significant | Conserved |
| Significant | Not significant | Additive |
| Not significant | Significant | Overdominant when Δ\_dom > 0, underdominant when Δ\_dom < 0 |
| Significant | Significant | Placed by the interval of k |

When both effects are significant, we placed the gene by the percentile interval of k over the bootstrap replicates, k^(b) = 2Δ\_dom^(b)/Δ\_total^(b), so that each class rests on the range of dominance the data allow rather than on the point estimate alone. A gene is Pc dominant when the interval contains 1 and Pe dominant when it contains −1. It is overdominant or underdominant, by the sign of Δ\_dom, when the interval lies entirely beyond 1 or −1, and partially dominant toward a parent when the interval lies entirely between 0 and that parent's value. An interval that contains both 1 and −1 arises when the parents barely differ, so the hybrid departs from both alike, and we classed the gene as overdominant or underdominant by the sign of Δ\_dom, as when the total contrast is not significant. An interval that contains 0 and the value of one parent places the hybrid between the midparent and that parent, and we classed the gene as partially dominant toward that parent. These intervals are built only for the genes that both tests select, so we set their level by the false coverage-statement rate, using marginal intervals of level 1 − 0.05 N\_sel/N\_test, with N\_sel the number of genes selected and N\_test the number tested in the family (Benjamini and Yekutieli 2005).

**Trans dominance.** The same scheme, with the *trans* contrast in place of the total contrast, classifies how the two parental *trans* environments combine in the hybrid. We tested Δ\_trans = 0 and Δ\_dom = 0 with the bootstrap test. A gene with neither significant has no detected *trans* divergence. A gene with only the *trans* contrast significant has a hybrid *trans* environment consistent with the log₂ midpoint of the parental environments, which is log-additive *trans* action. A gene with only Δ\_dom significant has a hybrid *trans* environment beyond those of both parents, above them when Δ\_dom > 0 and below them when Δ\_dom < 0. When both are significant, we placed the gene by the percentile interval of k\_trans with the rules and the coverage level used for k. Its *trans* environment is Pc dominant when the interval contains 1, Pe dominant when it contains −1, beyond both parents when the interval lies entirely beyond 1 or −1, and partially dominant toward a parent when the interval lies between 0 and that parent's value.

*Appendix E* refits every contrast with a cell-state program in the mean and flags the significant contrasts that the program carries.

## Statistical power of *cis* and *trans* estimates

In a hybrid, two alleles share one nucleus and one set of *trans*-acting factors. Differences in expression between the alleles are therefore a measures *cis*-regulatory divergence alone (Wittkopp et al. 2004). A comparison of the hybrid difference to the difference between the parents then reveals *trans*-regulatory divergence. Since the parental ratio for this *trans *contrast comes from two independent sets of cells, its variance adds to that of the hybrid difference. Thus, if the same hybrid samples were used for both the *cis *and *trans *contrasts, *trans* would always be the less precise estimate, with lower power than *cis* to detect differences. In addition, the estimation error for *trans* would correlate negatively with that of *cis* (Fraser 2019). We resolved both problems by using a portion of the hybrid cells to estimate *cis* and the remaining hybrid cells to estimate *trans*.

```latex
N_{\mathrm{cis}} = \chi\,N_{\mathrm{H}}, \qquad N_{\mathrm{trans}} = (1-\chi)\,N_{\mathrm{H}}
```

where N\_H is the number of hybrid cells and χ the fraction assigned to *cis*.

Since the two subsets share no hybrid cells, the *cis* and *trans* estimates are independent and estimation errors in one do not affect the other. We also chose the partition such that the two contrasts have equal sampling variance, and therefore equal power to detect an effect of a given size. The dominance contrast Δ\_dom, the allelic inheritance contrasts and k\_trans are estimated on the trans cells as well, so every contrast that draws on the hybrid alleles, apart from the cis contrast, is independent of the cis estimate, and the partition adds no correlation between cis and any of them.

**Variance of the *cis* contrast.** Both hybrid alleles are fit on the same *cis* cells, so the variance of the difference of their log₂ estimates includes their covariance. In large samples each variance and the covariance are inversely proportional to the number of cells used. The whole sum is therefore a constant divided by χN\_H. We call this constant σ²\_Δ,H, the per-cell variance of the hybrid allele ratio.

```latex
\begin{aligned}
\mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Hc}} - \log_2 \hat{x}_{\mathrm{He}}\right] &= \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Hc}}\right] + \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{He}}\right] - 2\,\mathrm{Cov}\left[\log_2 \hat{x}_{\mathrm{Hc}},\ \log_2 \hat{x}_{\mathrm{He}}\right] \\
&= \frac{\sigma^2_{\Delta,\mathrm{H}}}{\chi\,N_{\mathrm{H}}}
\end{aligned}
```

**Variance of the *trans* contrast.** The *trans* contrast subtracts the hybrid allele ratio, measured on the remaining hybrid cells, from the parental ratio. The two parents are measured in different cells, so their estimates have no covariance, and the variance of the parental ratio is the sum of two terms.

```latex
\sigma^2_{\Delta,P} = \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Pc}}\right] + \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Pe}}\right]
```

The parental cells and the *trans* hybrid cells are also disjoint, so the two parts of the *trans* contrast are independent and their variances add with no covariance. The fitted difference for the *trans* contrast is the parental log ratio minus the hybrid allele log ratio, with the two hybrid alleles fit only on the cells of the *trans* subset. Substituting these gives:

```latex
\mathrm{Var}\left\{\left[\log_2 \hat{x}_{\mathrm{Pc}} - \log_2 \hat{x}_{\mathrm{Pe}}\right] - \left[\log_2 \hat{x}_{\mathrm{Hc}} - \log_2 \hat{x}_{\mathrm{He}}\right]\right\} = \sigma^2_{\Delta,P} + \frac{\sigma^2_{\Delta,\mathrm{H}}}{(1-\chi)\,N_{\mathrm{H}}}
```

**B****alancing**** statistical power.** To balance statistical power between *cis* and *trans* estimates, we set the *cis* variance equal to the *trans* variance and solved for χ. Multiplying both sides by the number of hybrid cells and dividing by the per-cell variance gives:

```latex
\begin{aligned}
\frac{\sigma^2_{\Delta,\mathrm{H}}}{\chi\,N_{\mathrm{H}}} &= \sigma^2_{\Delta,P} + \frac{\sigma^2_{\Delta,\mathrm{H}}}{(1-\chi)\,N_{\mathrm{H}}} \\
\frac{1}{\chi} &= \frac{N_{\mathrm{H}}\,\sigma^2_{\Delta,P}}{\sigma^2_{\Delta,\mathrm{H}}} + \frac{1}{1-\chi}
\end{aligned}
```

The partition now depends on the data through a single ratio, γ. It divides the parental variance by the variance a hybrid allele ratio would have if all hybrid cells were used, so it measures the fixed parental variance in units of the best precision the hybrid can achieve.

```latex
\gamma = \frac{N_{\mathrm{H}}\,\sigma^2_{\Delta,P}}{\sigma^2_{\Delta,\mathrm{H}}} \quad\Longrightarrow\quad \frac{1}{\chi} = \gamma + \frac{1}{1-\chi}
```

Multiplying through by both denominators and collecting terms gives a quadratic in the *cis* fraction with two roots.

```latex
\begin{aligned}
1 - \chi &= \gamma\,\chi\,(1-\chi) + \chi \quad\Longrightarrow\quad \gamma\,\chi^2 - (2 + \gamma)\,\chi + 1 = 0 \\
\Longrightarrow\quad \chi &= \frac{(2 + \gamma) \pm \sqrt{(2 + \gamma)^2 - 4\gamma}}{2\gamma} = \frac{(2 + \gamma) \pm \sqrt{\gamma^2 + 4}}{2\gamma}
\end{aligned}
```

The root with the minus sign is the balanced partition (the root with the plus sign exceeds one).

When the parental variance is negligible compared to the hybrid variance (γ → 0), the optimum is an even partition. When the parental variance is much larger than the hybrid variance, the partition falls as the inverse of the ratio. Power is therefore equalized by moving hybrid cells from estimating *cis* to estimating *trans* as the parental variance grows.

```latex
\begin{aligned}
\gamma \to 0&: \quad \sqrt{\gamma^2 + 4} \approx 2 + \frac{\gamma^2}{4} \;\Longrightarrow\; \chi \approx \frac12 - \frac{\gamma}{8} \\
\gamma \to \infty&: \quad \chi \approx \frac{1}{\gamma}
\end{aligned}
```

**Separate partitions for each quantity.** The optimal partition depends on the parental variance relative to the hybrid variance, which differs among quantities. Estimates of φ carry larger and more skewed sampling error than estimates of μ, and kinetic balance and burst size carry additional sampling error beyond this. We therefore estimated γ separately for each. We set each partition from the log₂ variances, using the same partitions for the raw contrasts. Since genes already differ in their ratio of parental to hybrid variance, the partition pooled across genes makes the power of the *cis *and *trans *contrasts simply comparable rather than identical, on either scale.

For the kinetic balance and burst size partitions, the fits of the mean and the dispersion need to be on the same cells. When the counts follow the NB, the estimation error of the fitted mean and dispersion are uncorrelated. However, when the counts do not follow the NB, the covariance between these errors depends on the third central moment of the factor. This shrinks as 1/N, so the correlation does not fall as the sample grows (Appendix F). The bootstrap standard errors of the kinetic balance and burst size contrasts include this covariance.

**Estimating partition fractions.** To estimate χ for each quantity, we first estimated σ²\_Δ,H and σ²\_Δ,P. We bootstrapped cells from the undivided hybrid and from both parents, refitting the hybrid alleles on the same resampled cells in each replicate to preserve their covariance. Because per-gene ratios of squared bootstrap standard errors are noisy and right-skewed, we pooled the variance components across genes by their medians. The pooled γ is the ratio of their medians across genes. For the dispersion, burst size and kinetic balance the medians ran over genes with no zero or failed bootstrap replicates in the groups compared.&#32;

```latex
\begin{aligned}
\hat{\sigma}^2_{\Delta,\mathrm{H},g} &= N_{\mathrm{H}}\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Hc},g} - \log_2 \hat{x}_{\mathrm{He},g}\right], \qquad \hat{\sigma}^2_{\Delta,P,g} = \widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Pc},g} - \log_2 \hat{x}_{\mathrm{Pe},g}\right] \\
\hat{\gamma} &= \frac{N_{\mathrm{H}}\,\mathrm{median}_g\,\hat{\sigma}^2_{\Delta,P,g}}{\mathrm{median}_g\,\hat{\sigma}^2_{\Delta,\mathrm{H},g}} = \frac{\mathrm{median}_g\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Pc},g} - \log_2 \hat{x}_{\mathrm{Pe},g}\right]}{\mathrm{median}_g\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Hc},g} - \log_2 \hat{x}_{\mathrm{He},g}\right]}
\end{aligned}
```

Since these estimates do not use *cis* or *trans* point estimates (only variance), they do not bias the *cis* and *trans* tests (Wittes and Brittain 1990). Within each partition, the hybrid cells used for *cis* and *trans* span the full range of sequencing depth.

## Multiple testing

We tested every contrast in every gene, so we controlled the false discovery rate across genes with the Benjamini–Hochberg procedure (Benjamini and Hochberg 1995). Each combination of contrast and quantity formed its own family. The contrasts are the total, *cis*, *trans* and dominance contrasts. The quantities are the log₂ mean, the log₂ and raw dispersion, the raw burst frequency, and the log₂ and raw burst size and kinetic balance. Each family was adjusted over the genes that it tested. The moment estimates of *Appendix C* only flag discrepant contrasts, so they were not tested. Co-expression pairs formed families of their own (Co-expression).

A contrast whose observed fit failed, or a log₂, raw burst frequency or raw kinetic balance contrast whose observed estimate is zero, has no p-value (*Contrasts*), so it lies outside its family. These rules change the number of tests in a family but not the p-values of the genes that remain in it. The p-values themselves come from the bootstrap tests of *Contrasts*, and the add-one count keeps every p-value above zero, so that each remains a valid p-value for a finite number of resamples (Phipson and Smyth 2010).

We called a contrast significant when its adjusted p-value was below 0.05. Every call in the analysis reads these adjusted values. These include the classification of each gene by its *cis* and *trans* contrasts, the classification of inheritance (*Dominance*), and the flags of *Appendix E*.

## Co-expression

**Residual correlations.** To ask whether the connections among genes diverged, we correlated the noise of pairs of genes across cells. For gene g in cell i of group r we took the Pearson residual of the offset NB fit, with the split-corrected dispersion for a hybrid allele. The residual removes the mean and the library depth and puts every gene on a common variance scale, so the correlation of two residuals measures how the genes fluctuate together. A residual with a fitted mean of zero was set to zero:

```latex
\rho_{gg',r} = \mathrm{cor}_i\left(\frac{y_{gri} - \hat{M}_{gri}}{\sqrt{\hat{M}_{gri} + \hat{\varphi}_{gr}\hat{M}_{gri}^2}},\ \frac{y_{g'ri} - \hat{M}_{g'ri}}{\sqrt{\hat{M}_{g'ri} + \hat{\varphi}_{g'r}\hat{M}_{g'ri}^2}}\right), \qquad \hat{M}_{gri} = O_{ri}\,\hat{\mu}_{gr}
```

**Genes.** A correlation between two genes whose counts are dominated by Poisson noise mostly measures that noise. At its mean, the share of a gene's count variance that is biological is φμ/(1 + φμ). We kept the 250 most highly expressed genes of the *S. cerevisiae* parent and every other gene whose biological share exceeded 0.3 in all four groups, Pc, Pe, Hc and He.

**Shrinkage.** With many genes and fewer cells, sample correlations are noisy, so we shrank each correlation matrix toward the identity (Schäfer and Strimmer 2005). Every correlation off the diagonal was multiplied by one minus the analytic shrinkage intensity, the summed estimated variance of the pairwise correlations divided by their summed squares, clipped to lie between zero and one. Noisier matrices therefore contract more.

**Decomposition.** We decomposed the divergence of each pair's correlation into the contrasts used for single genes (*Contrasts*). The *cis* contrast compares the correlations of the pair within each genome in the same hybrid cells. The *trans* contrast is the part of the parental difference that the *cis* contrast leaves unexplained. The dominance contrast compares the average of the two within-genome correlations in the hybrid with the average of the two parental correlations, so the hybrid is compared with its haploid parents genome by genome. Correlations were compared on their own scale:

```latex
\begin{aligned}
\Delta^{\mathrm{cis}}_{gg'} &= \rho_{gg',\mathrm{Hc}} - \rho_{gg',\mathrm{He}} \\
\Delta^{\mathrm{trans}}_{gg'} &= \left(\rho_{gg',\mathrm{Pc}} - \rho_{gg',\mathrm{Pe}}\right) - \left(\rho_{gg',\mathrm{Hc}} - \rho_{gg',\mathrm{He}}\right) \\
\Delta^{\mathrm{dom}}_{gg'} &= \tfrac12\left(\rho_{gg',\mathrm{Hc}} + \rho_{gg',\mathrm{He}}\right) - \tfrac12\left(\rho_{gg',\mathrm{Pc}} + \rho_{gg',\mathrm{Pe}}\right)
\end{aligned}
```

The total contrast, the allelic inheritance contrasts, k and k\_trans follow from these three as for single genes (*Derived comparisons*). The allelic inheritance contrasts of a pair compare its correlation within one genome in the hybrid with its correlation in that genome's parent, so they show whether the connection between two genes moves when its genome enters the hybrid *trans* environment. Because the hybrid is noisier, its within-genome correlations can differ from the parental correlations by a shift shared by most pairs. We grouped pairs into deciles of the average observed log₂ μ\_H of their two genes and subtracted from Δ^dom\_gg′ the median of Δ^dom\_gg′ over the pairs of its decile, recomputing the median in every bootstrap replicate. A shift shared by both genomes of the hybrid cancels in the cis and trans contrasts.

**Inference for pairs.** We drew 5000 bootstrap replicates. Each replicate resampled the cells of each parent, and resampled the hybrid cells once for both alleles, so the two alleles of a cell stayed paired. The residuals were held at the observed fits, and each replicate recomputed the correlations, the shrinkage and the contrasts. The standard error of a contrast is the standard deviation of its replicates, and its p-value is two-sided under a normal approximation to the ratio of the estimate to this standard error. We adjusted these p-values across all pairs within each contrast (Multiple testing), called a contrast significant when its adjusted p-value was below 0.05, and classified each pair by its *cis* and *trans* contrasts, and its inheritance by the total, *trans* and dominance contrasts and the percentile intervals of k and k\_trans, as for single genes (*Dominance*).

**Axes of divergence.** We summarized each contrast matrix by its eigenvectors, ordered by the magnitude of their eigenvalues. Because the total contrast is the sum of the *cis* and *trans* contrasts for every pair, the eigenvalue of any total axis u divides exactly into u′Δ^cis u and u′Δ^trans u. We tested the squared eigenvalue at each of the fifteen leading ranks against the squared eigenvalue at the same rank in null draws. The total null pools the parental cells and reassigns them to two groups of the original sizes, the *cis* null swaps the residuals of the two alleles in a random subset of hybrid cells, and the *trans* null is the difference of the two. The dominance null works within each genome. It pools the cells of the *S. cerevisiae* parent with the *S. cerevisiae* residuals of the hybrid cells, and the cells of the *S. eubayanus* parent with the *S. eubayanus* residuals, reassigns each pool to parental and hybrid groups of the original sizes, and recomputes Δ^dom with the same shared-shift adjustment as the observed matrix.

## Model assumptions

| Assumption | Statement | Used in |
| --- | --- | --- |
| Independent cells | Counts of different cells are independent given their offsets | Negative binomial model, Statistical power, Contrasts |
| Bursty limit | The NB arises from bursty transcription, so θ = 1/φ is the burst frequency and β = μφ the burst size | Negative binomial model, Burst kinetics |
| Fixed offsets | Offsets are library sizes divided by D and are held fixed, although a gene's own counts enter its library size. Each hybrid allele takes the library of its own genome | Sequencing depth offset, Library split correction |
| Fitted dispersion | The NB fit converges to the variance of the expression factor only when that factor is gamma | Negative binomial model, Library split correction |
| Large samples | Every group has about 800 cells or more, so terms of order 1/N are small and are dropped, except the plug-in bias of φ̂, which the Cox–Reid adjustment removes to first order. Partition sizes are reported in Statistical power | Estimating φ, Library split correction, Contrasts |
| Capture efficiency | UMI capture thins counts binomially, which leaves φ and θ unchanged but scales μ and β. A capture difference between parental cells, or between parental and hybrid cells, enters the total, trans and dominance contrasts in μ and β. A genome-wide capture difference between the hybrid alleles is removed by the allele offsets | Burst kinetics, Contrasts |
| Hybrid allele counts | Allele counts follow the parental NB form with the true allele library as offset, and the observed allele library differs from it by the library split | Library split correction, Appendix D |
| Independent split error | The split error of a hybrid library is independent of the expression factors and of the count of any one gene | Library split correction, Appendix D |
| Binomial library split | Genome-assigned UMIs are binomial given the hybrid library size and s\_i | Library split correction, Appendix D |
| Split correction | The split corrections of the hybrid dispersions and variance estimates are computed from the observed libraries and hold to first order in the split variance | Library split correction, Contrasts, Appendix D |
| Zero estimates | A zero dispersion estimate is a valid boundary value and is kept, and a failed fit counts as zero in the bootstrap. A failed observed fit gives no contrast, and a zero observed estimate gives no log₂, raw burst frequency or raw kinetic balance contrast | Contrasts |
| Ploidy effects | Ploidy effects act equally on both alleles, so they cancel in the cis ratio and within each ratio of the trans contrast. A species-specific ploidy effect would remain in the trans divergence, and separating it requires diploid parental strains. A shift shared by both hybrid alleles does not cancel in Δ\_dom and is removed by the shared-shift adjustment | Contrasts |
| Relative measurement | Rates are measured relative to library size, so each group and each hybrid allele is known up to a scale factor of its library. Averages of log₂ values turn these factors into a shift shared by every gene | Contrasts |
| Cis effects independent of background | A cis sequence has the same effect in its parent and in the hybrid. A cis effect that changes with genetic background appears in the trans contrast, and Δ\_dom measures trans dominance only to the extent that this holds | Contrasts, Dominance |
| Cross-allele covariance | In the cis kinetic contrasts, the covariance between the mean estimate of one allele and the dispersion estimate of the other is negligible | Contrasts, Kinetic balance |
| Common partition ratio | One ratio γ of parental to hybrid variance, pooled by medians across genes, sets the partition for every gene of a quantity | Statistical power |
| Bootstrap test | Every contrast is tested by the side of zero on which its bootstrap replicates fall. Resampling hybrid cells whole keeps the extrinsic noise that the two alleles share | Contrasts |
| Extrinsic and intrinsic parts | Allele expression factors divide into a shared extrinsic part and equal private intrinsic parts | Contrasts, Appendix D |
| Working weights | The weighted moment estimator uses second moments alone and needs no gamma assumption. Its estimating equation is unbiased only with true means and fixed weights | Appendix C |
| Expression factors and offset | Expression factors are independent of the offset | Appendix D |
| Separate cultures | Parental and hybrid cells come from different cultures, so the *trans* and dominance contrasts carry any culture difference. Conditioning on the program score removes only the measured part | Appendix E |
| Log-linear state effect | A gene's mean count changes with the program score as 2^(τ\_r S\_i), with the dispersion constant around that mean | Appendix E |
| Score independent of share | The transcriptome share of a hybrid cell is independent of its program score, so the split correction applies to the conditioned dispersions unchanged | Appendix E |

## Summary of estimators

The table gives where each estimator is defined and where it is used, and the display gives each estimator as it is computed. A corrected estimate below zero is set to zero after its correction.

| Quantity | Defined in | Used in |
| --- | --- | --- |
| Depth offset | Sequencing depth offset | Every fit |
| Mean rate | Estimating μ | Every mean contrast, fitted means M̂\_i |
| Dispersion before any correction | Estimating φ | Parental dispersions directly, and φ̂\*\_Hc and φ̂\*\_He for the hybrid alleles |
| Split variance | Library split correction | Corrected allele dispersions, Appendix D |
| Allele mean | Library split correction | Cis, trans and dominance mean contrasts |
| Allele dispersion | Library split correction | Cis and trans dispersion contrasts, per-genome hybrid dispersion, kinetic quantities |
| Per-genome hybrid mean | Contrasts, Contrasts between groups | Dominance mean contrasts and kinetic quantities |
| Per-genome hybrid dispersion | Contrasts, Contrasts between groups | Dominance dispersion contrasts and kinetic quantities |
| Partition fraction | Statistical power | Cis and trans partitions of the hybrid cells |
| Kinetic quantities | NB model, Burst kinetics | Burst frequency, burst size and kinetic balance contrasts |
| Ploidy shift | Appendix D | Reported only |

```latex
\begin{aligned}
&\text{Depth offset} && O_{ri} = L_{ri}/D, \qquad L_{\mathrm{Hc}i} \ \text{and} \ L_{\mathrm{He}i} \ \text{for the hybrid alleles} \\
&\text{Mean rate} && \hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_{ri}}, \qquad \hat{M}_i = O_{ri}\,\hat{\mu}_r \\
&\text{Dispersion} && \hat{\varphi}_r = \arg\max_{\varphi > 0} \ell^{\mathrm{CR}}_r(\varphi), \qquad \hat{\varphi}_r = 0 \ \text{if} \ \sum_i (y_{ri} - \hat{M}_i)^2 \le \sum_i \hat{M}_i - \frac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i} \\
&\text{Contrast p-value} && p = \min\left\{1,\ 2\left(1 + n_0\right)/\left(1 + B\right)\right\}, \qquad n_0 = \text{replicates on the null side of zero, ties and undefined replicates included} \\
&\text{Split variance} && \hat{s}_i = \frac{L_{\mathrm{Hc}i}}{L_{\mathrm{H}i}}, \qquad \hat{\eta}_{\mathrm{Hc}i} = \frac{1-\hat{s}_i}{\hat{s}_i\,L_{\mathrm{H}i}}, \qquad \hat{\eta}_{\mathrm{He}i} = \frac{\hat{s}_i}{(1-\hat{s}_i)\,L_{\mathrm{H}i}}, \qquad \hat{\delta}_r = \frac{\sum_i \hat{\eta}_{ri}\,\hat{M}_i^2/(1+\hat{\varphi}^*_r\hat{M}_i)^2}{-2\,\hat{\mathrm{E}}\left[\ell''_r\left(\hat{\varphi}^*_r\right)\right]} \\
&\text{Allele dispersion} && \hat{\varphi}_r = \hat{\varphi}^*_r - \hat{\delta}_r, \qquad r = \mathrm{Hc}, \mathrm{He} \\
&\text{Allele mean} && \hat{\mu}_{\mathrm{Hc}}, \ \hat{\mu}_{\mathrm{He}} \ \text{from the mean rate, with no rescaling} \\
&\text{Per-genome hybrid mean} && \log_2\hat{\mu}_{\mathrm{H}} = \tfrac12\left(\log_2\hat{\mu}_{\mathrm{Hc}} + \log_2\hat{\mu}_{\mathrm{He}}\right) \\
&\text{Per-genome hybrid dispersion} && \log_2\hat{\varphi}_{\mathrm{H}} = \tfrac12\left(\log_2\hat{\varphi}_{\mathrm{Hc}} + \log_2\hat{\varphi}_{\mathrm{He}}\right) \\
&\text{Partition fraction} && \chi = \frac{(2 + \hat{\gamma}) - \sqrt{\hat{\gamma}^2 + 4}}{2\hat{\gamma}}, \ \ \text{one each for } \mu, \varphi, \kappa, \beta \\
&\text{Kinetic quantities} && \hat{\theta} = 1/\hat{\varphi}, \qquad \hat{\beta} = \hat{\mu}\,\hat{\varphi}, \qquad \hat{\kappa} = \hat{\theta}/\hat{\beta} \\
&\text{Ploidy shift (reported)} && \hat{\Lambda} = -\log_2\left[1 - 2\hat{\pi}(1-\hat{\pi})\,\hat{\omega}\right], \qquad \hat{\Lambda}_{\mathrm{raw}} = 2\hat{\pi}(1-\hat{\pi})\,\hat{\sigma}^2_q \ \ \text{(Appendix D)}
\end{aligned}
```

## Appendix A. Derivatives of the NB log-likelihood

This appendix collects the derivatives and expectations used in the main text and in the other appendices. For one cell with count y, expected count M, dispersion φ and Fano factor F = 1 + φM, the log probability is written with the gamma term Ψ (Estimating μ), and its derivatives in the mean are:

```latex
\begin{aligned}
\log\Pr(y \mid M, \varphi) &= \Psi(\varphi) + y\log M - \log y! - \left(\varphi^{-1} + y\right)\log F, \qquad \Psi(\varphi) = \sum_{j=0}^{y-1}\log\left(1 + j\varphi\right) \\
\frac{\partial}{\partial M}\log\Pr(y \mid M, \varphi) &= \frac{y - M}{M F}, \qquad \frac{\partial}{\partial \mu}\log\Pr(y \mid O\mu, \varphi) = \frac{y - M}{\mu F}
\end{aligned}
```

Ψ does not depend on the mean, so it does not enter these derivatives. Primes mark derivatives in φ at fixed M. Only Ψ carries a sum over the count, and the other terms have closed forms:

```latex
\begin{aligned}
\Psi'(\varphi) &= \sum_{j=0}^{y-1}\frac{j}{1 + j\varphi}, \qquad \Psi''(\varphi) = -\sum_{j=0}^{y-1}\frac{j^2}{(1 + j\varphi)^2} \\
\ell'(\varphi) &= \Psi'(\varphi) + \frac{\log F}{\varphi^2} - \frac{\left(\varphi^{-1} + y\right)M}{F} \\
\ell''(\varphi) &= \Psi''(\varphi) - \frac{2\log F}{\varphi^3} + \frac{2M}{\varphi^2 F} + \frac{\left(\varphi^{-1} + y\right)M^2}{F^2}
\end{aligned}
```

Here ℓ′ and ℓ″ are the terms of one cell, and the derivatives for a group sum them over its cells. The split correction and Appendix B replace M by a fitted value and expand in its error. Only the terms of ℓ′ and ℓ″ that contain M change, and their first and second derivatives in M are:

```latex
\begin{aligned}
\frac{\partial \ell'}{\partial M} &= \frac{M - y}{F^2}, \qquad \frac{\partial^2 \ell'}{\partial M^2} = \frac{1}{F^2} - \frac{2\varphi\left(M - y\right)}{F^3} \\
\frac{\partial\left(-\ell''\right)}{\partial M} &= \frac{2M\left(M - y\right)}{F^3}, \qquad \frac{\partial^2\left(-\ell''\right)}{\partial M^2} = \frac{2\left(2M - y\right)}{F^3} - \frac{6\varphi M\left(M - y\right)}{F^4}
\end{aligned}
```

The count Y is NB with mean M, and the split error ξ has mean zero and variance η and is independent of Y. Their moments, and the expectations of ℓ′ and −ℓ″ at the true mean, are:

```latex
\begin{aligned}
\mathrm{E}[Y] &= M, \qquad \mathrm{Var}[Y] = M F, \qquad \mathrm{E}\left[(Y - M)^3\right] = M F\left(1 + 2\varphi M\right), \qquad \mathrm{E}[\xi] = 0, \qquad \mathrm{E}\left[\xi^2\right] = \eta \\
\mathrm{E}\left[\ell'(\varphi)\right] &= 0, \qquad -\mathrm{E}\left[\ell''(\varphi)\right] = -\mathrm{E}\left[\Psi''(\varphi)\right] + \frac{2\log F - F + F^{-1}}{\varphi^3} \\
-\mathrm{E}\left[\Psi''(\varphi)\right] &= \mathrm{E}\left[\sum_{j=0}^{\infty}\frac{j^2}{(1 + j\varphi)^2}\,\mathbf{1}\left[Y > j\right]\right] = \sum_{j=0}^{\infty}\frac{j^2\Pr\left(Y > j\right)}{(1 + j\varphi)^2}
\end{aligned}
```

The closed-form part of −E\[ℓ″\] uses E\[Y\] = M and φM = F − 1. In −E\[Ψ″\], the sum up to Y − 1 equals a sum over every j of the same term times an indicator that j lies below the count, and the expectation of that indicator is Pr(Y > j). As φ → 0, −E\[Ψ″\] tends to M²/2 + M³/3 and the closed-form part to −M³/3, so −E\[ℓ″\] for one cell tends to M²/2.

## Appendix B. The Cox–Reid adjustment with the ratio mean

The Cox–Reid adjustment removes the first-order bias of φ̂\_r when μ is refit by maximum likelihood at each φ. We derive the bias that remains when μ is instead the ratio estimator μ̂\_r. To first order, the bias of φ̂\_r is the expected slope of the adjusted log-likelihood at the true φ divided by the expected information, where A′\_r is the slope of the Cox–Reid term:

```latex
\mathrm{E}\left[\hat{\varphi}_r\right] - \varphi \approx \frac{\mathrm{E}\left[\ell'_r\left(\varphi; \hat{\mu}_r\right)\right] + A'_r(\varphi)}{-\mathrm{E}\left[\ell''_r(\varphi)\right]}, \qquad A'_r(\varphi) = \frac{\partial}{\partial\varphi}\left(-\tfrac12\log\sum_i \frac{M_i}{F_i}\right) = \frac{\tfrac12\sum_i M_i^2/F_i^2}{\sum_i M_i/F_i}
```

The ratio estimator multiplies every fitted mean by the same relative error, ξ̄. Because ξ̄ is built from the counts of every cell, it has mean zero, its variance follows from the independence of cells, and it covaries with each count through that count's own variance, M\_i F\_i:

```latex
\hat{M}_i = M_i\left(1 + \bar{\xi}\right), \qquad \bar{\xi} = \frac{\sum_j \left(Y_j - M_j\right)}{M_{\mathrm{tot}}}, \qquad \mathrm{E}\left[\bar{\xi}\right] = 0, \qquad \mathrm{E}\left[\bar{\xi}^2\right] = \frac{\sum_j M_j F_j}{M_{\mathrm{tot}}^2}, \qquad \mathrm{Cov}\left[\bar{\xi}, Y_i\right] = \frac{M_i F_i}{M_{\mathrm{tot}}}
```

As in the split correction, only the terms of ℓ′\_r that contain the fitted mean change, and their expectation at the true means is zero. Their first and second derivatives in M̂\_i are given in Appendix A. We expand to second order in ξ̄ and distribute the expectation. Unlike the split error, ξ̄ is not independent of the counts, so the linear term keeps its covariance with Y\_i. The term E\[ξ̄²(Y\_i − M\_i)\] is a third moment of order 1/N² per cell, and we drop it:

```latex
\begin{aligned}
\mathrm{E}\left[\ell'_r\left(\varphi; \hat{\mu}_r\right)\right] &\approx \sum_i \mathrm{E}\left[\frac{M_i\,\bar{\xi}\left(M_i - Y_i\right)}{F_i^2} + \tfrac12 M_i^2\,\bar{\xi}^2\left(\frac{1}{F_i^2} - \frac{2\varphi\left(M_i - Y_i\right)}{F_i^3}\right)\right] \\
&= \sum_i\left[-\frac{M_i\,\mathrm{Cov}\left[\bar{\xi}, Y_i\right]}{F_i^2} + \frac{M_i^2\,\mathrm{E}\left[\bar{\xi}^2\right]}{2F_i^2}\right] \\
&= -\frac{1}{M_{\mathrm{tot}}}\sum_i \frac{M_i^2}{F_i} + \frac{\sum_j M_j F_j}{2M_{\mathrm{tot}}^2}\sum_i \frac{M_i^2}{F_i^2}
\end{aligned}
```

Adding the Cox–Reid slope gives R\_r, the expected slope that the adjustment leaves, and the residual bias of φ̂\_r:

```latex
R_r = -\frac{1}{M_{\mathrm{tot}}}\sum_i \frac{M_i^2}{F_i} + \frac{\sum_j M_j F_j}{2M_{\mathrm{tot}}^2}\sum_i \frac{M_i^2}{F_i^2} + \frac{\sum_i M_i^2/F_i^2}{2\sum_i M_i/F_i}, \qquad \mathrm{E}\left[\hat{\varphi}_r\right] - \varphi \approx \frac{R_r}{-\mathrm{E}\left[\ell''_r(\varphi)\right]}
```

R\_r vanishes in two limits. With equal offsets every M\_i equals M, and the three terms are −M/F, M/(2F) and M/(2F), which sum to zero, as expected because μ̂\_r is then the maximum likelihood mean. As φ → 0, every F\_i → 1 and the three terms are −ΣM\_i²/M\_tot, ΣM\_i²/(2M\_tot) and ΣM\_i²/(2M\_tot), which also sum to zero. When φM\_i is large in every cell, F\_i ≈ φM\_i, and the residual grows with the spread of the offsets:

```latex
R_r \approx \frac{1}{2\varphi}\left(\frac{N_r\sum_i M_i^2}{M_{\mathrm{tot}}^2} - 1\right) = \frac{1}{2\varphi}\left(\frac{N_r}{N_{\mathrm{eff}}} - 1\right)
```

A Monte Carlo check with 2 × 10⁵ replicates of 40 to 80 cells, offset CVs of 0 to 1, mean counts of 0.5 to 5 and φ of 0.3 to 1 matched R\_r within 1.5 standard errors in every setting. At 800 cells with lognormal offsets of CV 0.72 and 0.87, as stand-ins for the parental offsets, mean counts of 0.2 to 5 and φ of 0.1 to 1, the residual bias was below 0.3% of φ and below 0.01 standard deviations of φ̂\_r. The bias that the adjustment removes in the same settings ranged from 0.2% to 7% of φ. The adjustment therefore removes almost all of the first-order bias with our offsets, and what remains is negligible next to the sampling error of φ̂\_r.

## Appendix C. Weighted method of moments for the dispersion

The Cox–Reid adjusted NB likelihood estimate is the primary estimator of φ. To flag contrasts where a departure from the gamma assumption may bias the Cox–Reid fit, we also estimated φ from the first two moments of the counts, without requiring them to follow an NB distribution. For the count Y\_i of cell i, the mean is M\_i = O\_i μ. The variance can exceed the Poisson variance M\_i by a term proportional to M\_i².

```latex
\mathrm{E}[Y_i] = O_i\,\mu = M_i, \qquad \mathrm{Var}[Y_i] = \mathrm{E}\left[(Y_i - M_i)^2\right] = M_i + \varphi\,M_i^2
```

Subtracting the Poisson variance M\_i from the squared deviation leaves a random variable whose expectation is this excess variance:

```latex
\mathrm{E}\left[(Y_i - M_i)^2 - M_i\right] = \mathrm{E}\left[(Y_i - M_i)^2\right] - M_i = \left(M_i + \varphi\,M_i^2\right) - M_i = \varphi\,M_i^2
```

If the means were known, we could produce a method of moments estimate of φ by replacing the expectation with the observed values. Each cell provides one such moment, and we averaged them with non-negative weights w\_i so that more precise cells count more:

```latex
\mathrm{E}\left[\sum_i w_i\left((Y_i - M_i)^2 - M_i\right)\right] = \varphi\sum_i w_i\,M_i^2 \\
\sum_i w_i\left[(y_i - M_i)^2 - M_i\right] = \tilde{\varphi}\sum_i w_i\,M_i^2 \\
\tilde{\varphi} = \frac{\sum_i w_i\left[(y_i - M_i)^2 - M_i\right]}{\sum_i w_i\,M_i^2}

```

The means are not known, so we replaced M\_i by its fitted value M̂\_i in φ̃ to give a naive estimator φ̃(M̂).

```latex
\tilde{\varphi}(\hat{M}) = \frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{M}_i\right]}{\sum_i w_i\,\hat{M}_i^2}
```

The offsets are fixed, so the expectation of the fitted mean is the offset times the sum of the expected counts over the sum of the offsets, which equals M\_i. The residual therefore has expectation zero, and the expectation of its square is its variance.

```latex
\begin{aligned}
\mathrm{E}[\hat{M}_i] &= O_i\,\mathrm{E}[\hat{\mu}] = O_i\,\frac{\sum_j \mathrm{E}[Y_j]}{\sum_j O_j} = O_i\,\frac{\sum_j O_j\,\mu}{\sum_j O_j} = O_i\,\mu = M_i \\
\mathrm{E}[Y_i - \hat{M}_i] &= \mathrm{E}[Y_i] - \mathrm{E}[\hat{M}_i] = M_i - M_i = 0 \\
\mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] &= \left(\mathrm{E}[Y_i - \hat{M}_i]\right)^2 + \mathrm{Var}[Y_i - \hat{M}_i] = \mathrm{Var}[Y_i - \hat{M}_i]
\end{aligned}
```

To find this variance, we write O\_tot for the total offset, the sum of the offsets over cells, so that cell i carries the fraction O\_i/O\_tot of it.

```latex
\hat{M}_i = O_i\,\hat{\mu} = O_i\,\frac{\sum_j y_j}{\sum_j O_j} = \frac{O_i}{O_{\mathrm{tot}}}\sum_j y_j, \qquad O_{\mathrm{tot}} \equiv \sum_j O_j
```

Pulling out the part that depends on Y\_i from the residual produces:

```latex
\begin{aligned}
Y_i - \hat{M}_i &= Y_i - \frac{O_i}{O_{\mathrm{tot}}}\sum_j Y_j = Y_i - \frac{O_i}{O_{\mathrm{tot}}}\,Y_i - \frac{O_i}{O_{\mathrm{tot}}}\sum_{j\neq i} Y_j = \left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right)Y_i - \frac{O_i}{O_{\mathrm{tot}}}\sum_{j\neq i} Y_j \\
\mathrm{Var}[Y_i - \hat{M}_i] &= \left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right)^2\mathrm{Var}[Y_i] + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_{j\neq i}\mathrm{Var}[Y_j] = \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\mathrm{Var}[Y_i] + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \mathrm{Var}[Y_j]
\end{aligned}
```

Substituting in the variance and using Σ\_j M\_j = M\_i O\_tot/O\_i separates a Poisson part from a part proportional to φ.

```latex
\begin{aligned}
\mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] &= \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\left(M_i + \varphi M_i^2\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\left(\sum_j M_j + \varphi\sum_j M_j^2\right) \\
&= \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\left(M_i + \varphi M_i^2\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\left(\frac{M_i\,O_{\mathrm{tot}}}{O_i} + \varphi\sum_j M_j^2\right) \\
&= \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)M_i + \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\varphi M_i^2 + \frac{O_i}{O_{\mathrm{tot}}}\,M_i + \varphi\,\frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j M_j^2 \\
&= M_i\left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right) + \varphi\left[M_i^2\left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j M_j^2\right] = A_i + \varphi\,C_i \\
A_i &\equiv M_i\left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right), \qquad C_i \equiv M_i^2\left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j M_j^2
\end{aligned}
```

For the numerator of the naive estimator we then have:

```latex
\begin{aligned}
\mathrm{E}\left[(Y_i - \hat{M}_i)^2 - \hat{M}_i\right] &= \mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] - \mathrm{E}[\hat{M}_i] = A_i + \varphi\,C_i - M_i \\
&= M_i\left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right) + \varphi\,C_i - M_i = \varphi\,C_i - \frac{O_i}{O_{\mathrm{tot}}}\,M_i
\end{aligned}
```

With known means the numerator would have expectation φM\_i². With fitted means it does not, so this estimator is biased. The bias term (O\_i/O\_tot)M\_i does not scale with φ. We therefore defined a corrected numerator that adds its estimate, (O\_i/O\_tot)M̂\_i, whose expectation is (O\_i/O\_tot)M\_i because the fitted mean has expectation M\_i. The expectation of the corrected numerator is then:

```latex
\mathrm{E}\left[(Y_i - \hat{M}_i)^2 - \hat{M}_i + \frac{O_i}{O_{\mathrm{tot}}}\,\hat{M}_i\right] = \left(\varphi\,C_i - \frac{O_i}{O_{\mathrm{tot}}}\,M_i\right) + \frac{O_i}{O_{\mathrm{tot}}}\,M_i = \varphi\,C_i
```

The weighted sum of these terms minus φ Σ w\_i C\_i therefore has expectation zero at the true φ. Solving this for φ gives an unbiased estimator (φ̌) if C\_i is known.

```latex
\mathrm{E}\left[\sum_i w_i\left((Y_i - \hat{M}_i)^2 - \hat{M}_i + \frac{O_i}{O_{\mathrm{tot}}}\,\hat{M}_i\right) - \varphi\sum_i w_i\,C_i\right] = 0 \\
\check{\varphi} \equiv \frac{\sum_i w_i\left(\left[(y_i - \hat{M}_i)^2 - \hat{M}_i\right] + \frac{O_i}{O_{\mathrm{tot}}}\,\hat{M}_i\right)}{\sum_i w_i\,\hat{C}_i} = \frac{\sum_i w_i\left((y_i - \hat{M}_i)^2 - \hat{A}_i\right)}{\sum_i w_i\,\hat{C}_i}, \\ \hat{A}_i \equiv \hat{M}_i\left(1 - \frac{O_i}{O_{\mathrm{tot}}}\right), \quad \hat{C}_i \equiv \hat{M}_i^2\left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \hat{M}_j^2
```

A check marks this second intermediate estimate. It uses the fitted means and the plug-in Ĉ\_i, which keeps it distinct from φ̃, written with the true means, and from the final estimate φ̂, which also corrects the bias of Ĉ\_i.

In practice we evaluated A\_i, C\_i and the weights at M̂\_i. Â\_i is linear in the means when the offsets are known and is therefore unbiased for A\_i. However, Ĉ\_i is quadratic in the means, so it has to be corrected.

```latex
\begin{aligned}
\mathrm{E}\left[\hat{C}_i\right] &= \mathrm{E}\left[\left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\hat{M}_i^2 + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \hat{M}_j^2\right] = \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\mathrm{E}\left[\hat{M}_i^2\right] + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \mathrm{E}\left[\hat{M}_j^2\right] \\
&= \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\left(M_i^2 + \mathrm{Var}[\hat{M}_i]\right) + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j\left(M_j^2 + \mathrm{Var}[\hat{M}_j]\right) \\
&= \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)M_i^2 + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j M_j^2 + \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\mathrm{Var}[\hat{M}_i] + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \mathrm{Var}[\hat{M}_j] \\
&= C_i + \left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)\mathrm{Var}[\hat{M}_i] + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j \mathrm{Var}[\hat{M}_j]
\end{aligned}
```

The variances of the fitted means remain. Since M̂\_j = O\_j μ̂, its variance is O\_j² Var\[μ̂\]. Substituting it collects the terms into a multiple of C\_i.

```latex
\begin{aligned}
\mathrm{Var}[\hat{M}_j] &= \mathrm{Var}[O_j\,\hat{\mu}] = O_j^2\,\mathrm{Var}[\hat{\mu}] = \frac{M_j^2}{\mu^2}\,\mathrm{Var}[\hat{\mu}] = M_j^2\,\mathrm{CV}^2[\hat{\mu}], \qquad \mathrm{CV}^2[\hat{\mu}] \equiv \frac{\mathrm{Var}[\hat{\mu}]}{\mu^2} \\
\mathrm{E}\left[\hat{C}_i\right] &= C_i + \mathrm{CV}^2[\hat{\mu}]\left[\left(1 - \frac{2O_i}{O_{\mathrm{tot}}}\right)M_i^2 + \frac{O_i^2}{O_{\mathrm{tot}}^2}\sum_j M_j^2\right] = \left(1 + \mathrm{CV}^2[\hat{\mu}]\right)C_i
\end{aligned}
```

Because μ̂ is the total count divided by the total offset, and the counts are independent across cells, its variance follows from the variance relationship for each count.

```latex
\begin{aligned}
\mathrm{CV}^2[\hat{\mu}] &= \frac{\mathrm{Var}[\hat{\mu}]}{\mu^2} = \frac{\mathrm{Var}\left[\sum_j Y_j\right]}{\mu^2\left(\sum_j O_j\right)^2} = \frac{\sum_j\left(M_j + \varphi\,M_j^2\right)}{M_{\mathrm{tot}}^2} = \frac{1}{M_{\mathrm{tot}}} + \varphi\,h \\ 
M_{\mathrm{tot}} &\equiv \sum_j M_j = \mu\sum_j O_j, 
\qquad h \equiv \sum_j \frac{O_j^2}{O_{\mathrm{tot}}^2} = \frac{\sum_j M_j^2}{M_{\mathrm{tot}}^2} \\
\mathrm{E}\left[\hat{C}_i\right] &= \left(1 + \mathrm{CV}^2[\hat{\mu}]\right)C_i = \left(1 + \frac{1}{M_{\mathrm{tot}}} + \varphi\,h\right) C_i
\end{aligned}
```

Dividing Ĉ\_i by the coefficient out front then gives an unbiased estimate of C\_i, if the coefficient is known.

```latex
\tilde{C}_i \equiv \frac{\hat{C}_i}{1 + 1/\sum_j y_j + \hat{\varphi}\,h}
```

Plugging C̃\_i into the denominator of φ̌ gives φ̂.

```latex
\begin{aligned}
\hat{\varphi} &= \frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{A}_i\right]}{\sum_i w_i\,\tilde{C}_i} = \left(1 + \frac{1}{\sum_j y_j} + \hat{\varphi}\,h\right)\frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{A}_i\right]}{\sum_i w_i\,\hat{C}_i} \\
&= \left(1 + \frac{1}{\sum_j y_j} + \hat{\varphi}\,h\right)\check{\varphi} \\
&= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} + h\,\check{\varphi}\,\hat{\varphi}
\end{aligned}
```

Collecting the terms in φ̂ and solving gives the corrected estimator.

```latex
\begin{aligned}
\hat{\varphi} - h\,\check{\varphi}\,\hat{\varphi} &= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} \\
\hat{\varphi}\left(1 - h\,\check{\varphi}\right) &= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} \\
\hat{\varphi} &= \frac{\left(1 + 1/\sum_j y_j\right)\check{\varphi}}{1 - h\,\check{\varphi}}
\end{aligned}
```

Finally, we chose the weights w\_i. φ̂ is approximately unbiased for fixed weights, so the weights mainly affect the precision of φ̂ and only slightly its remaining bias. We minimized the variance of φ̃ instead, because with M\_i and w\_i known it can be written exactly.

```latex
\begin{aligned}
\mathrm{Var}[\tilde{\varphi}] &= \mathrm{Var}\left[\frac{\sum_i w_i\left[(Y_i - M_i)^2 - M_i\right]}{\sum_i w_i\,M_i^2}\right] = \frac{\mathrm{Var}\left[\sum_i w_i\left[(Y_i - M_i)^2 - M_i\right]\right]}{\left(\sum_i w_i\,M_i^2\right)^2} \\
&= \frac{\sum_i \mathrm{Var}\left[w_i\left[(Y_i - M_i)^2 - M_i\right]\right]}{\left(\sum_i w_i\,M_i^2\right)^2} = \frac{\sum_i \mathrm{Var}\left[w_i\,(Y_i - M_i)^2 - w_i\,M_i\right]}{\left(\sum_i w_i\,M_i^2\right)^2} \\
&= \frac{\sum_i \mathrm{Var}\left[w_i\,(Y_i - M_i)^2\right]}{\left(\sum_i w_i\,M_i^2\right)^2} = \frac{\sum_i w_i^2\,\mathrm{Var}\left[(Y_i - M_i)^2\right]}{\left(\sum_i w_i\,M_i^2\right)^2}
\end{aligned}
```

What remains unknown is V\_i, the variance of the squared deviation of cell i. The variance of any random variable is the mean of its square minus the square of its mean. The deviation has expectation zero and we write E\[(Y\_i − M\_i)⁴\] as Kurt\[Y\_i\] Var\[Y\_i\]², where Kurt\[Y\_i\] is the kurtosis of the counts of cell i. We assumed that Kurt\[Y\_i\] is the same for every cell.

```latex
\begin{aligned}
V_i\equiv\mathrm{Var}\left[(Y_i - M_i)^2\right] &= \mathrm{E}\left[\left((Y_i - M_i)^2\right)^2\right] - \left(\mathrm{E}\left[(Y_i - M_i)^2\right]\right)^2 \\
&= \mathrm{E}\left[(Y_i - M_i)^4\right] - \mathrm{Var}[Y_i]^2 = \mathrm{Kurt}[Y_i]\,\mathrm{Var}[Y_i]^2 - \mathrm{Var}[Y_i]^2 \\
&= \left(\mathrm{Kurt}[Y_i] - 1\right)\mathrm{Var}[Y_i]^2 = \left(\mathrm{Kurt}[Y_i] - 1\right)M_i^2(1 + \varphi M_i)^2 \\
\mathrm{Kurt}[Y_i] &\equiv \frac{\mathrm{E}\left[(Y_i - M_i)^4\right]}{\mathrm{Var}[Y_i]^2}
\end{aligned}
```

To minimize the variance, we set the derivative with respect to the weights to zero and solved.

```latex
\begin{aligned}
\frac{\partial}{\partial w_i}\,\frac{\sum_j w_j^2\,V_j}{\left(\sum_j w_j\,M_j^2\right)^2} &= \frac{2\,w_i\,V_i}{\left(\sum_j w_j\,M_j^2\right)^2} - \frac{2\,M_i^2\sum_j w_j^2\,V_j}{\left(\sum_j w_j\,M_j^2\right)^3} = 0 \\
\Longrightarrow\quad w_i &= \frac{\sum_j w_j^2\,V_j}{\sum_j w_j\,M_j^2}\cdot\frac{M_i^2}{V_i}
\end{aligned}
```

The first factor in the solution is a sum over all cells, so it is the same for every cell and only scales the weights. Similarly, Kurt\[Y\_i\] − 1 is the same for every cell, so we dropped both factors.

```latex
\begin{aligned}
w_i &= \frac{\sum_j w_j^2\,V_j}{\sum_j w_j\,M_j^2}\cdot\frac{M_i^2}{V_i} \propto \frac{M_i^2}{V_i} = \frac{M_i^2}{\left(\mathrm{Kurt}[Y_i] - 1\right)M_i^2\,(1 + \varphi M_i)^2} \\
&= \frac{1}{\left(\mathrm{Kurt}[Y_i] - 1\right)(1 + \varphi M_i)^2} \propto \frac{1}{(1 + \varphi M_i)^2}
\end{aligned}
```

These weights are the optimum for φ̃ when the kurtosis is the same for every cell. We used them for φ̂, evaluated at M̂\_i, which departs from the derivation in three ways. First, the fitted means lower the variance, because they remove the part of the numerator that moves with the total count. Second, the kurtosis does not have to be the same for every cell, for example for NB counts, where it falls as the mean rises. Third, evaluating the weights at fitted values makes them random, which adds bias.

These assumptions have little consequence at the cell numbers in this study. We tested this by simulation at N = 800 and N = 1600, with M of 0.5, 2 or 5, φ of 0.3 or 1, offset CV of 0.5 or 1, and either NB or Poisson-lognormal counts. We compared the weights from this derivation with more complete weights that use the true fourth moment of each cell and include the effect of the fitted means. Both used the true M, φ and count distribution. The weights from this derivation were nearly as precise. The median standard deviation of φ̂ was 0.5% higher at N = 800 and 0.3% higher at N = 1600. Evaluating them at the fitted M̂\_i and an iterated φ̂ added at most 0.3% more at the median. The loss was concentrated when the true counts were NB with high mean and dispersion. Almost all of this loss comes from leaving out the effect of the fitted means. Bias differed from that of the more complete weights by at most 0.006, and it was at most 0.06 of the standard deviation for every choice.

The weights depend on φ, so we started from φ = 0, calculated the weights at M̂\_i, updated φ̂ with the estimating equation, restricted it to nonnegative values, and repeated until convergence or a maximum number of iterations. We recorded the convergence status, whether the estimate sat at the boundary, and the untruncated estimate.

## Appendix D. The ploidy shift

Summing the two alleles of a hybrid cell averages their private noise, so the summed counts are less dispersed than a single genome. We derive this ploidy shift, Λ, and its estimator, including the correction of the estimator for the library split. We reported Λ as a description of how much noise summing the two alleles hides.

**Summed counts.** Summing two alleles within a cell averages their private fluctuations, so the summed hybrid appears quieter than either allele alone. We derive the size of this effect from the decomposition of dispersion into intrinsic and extrinsic parts (Elowitz et al. 2002). Setting aside variation in the transcriptome share, the Hc count is Poisson around its library size over D times s̄, the mean transcriptome share weighted by library size, its expression per genome, and its expression factor z. The two allele counts of a cell are independent given these expression factors, and the sum of two independent Poisson counts is Poisson with the summed mean. The summed mean is therefore L\_Hi/D, the sum of the two allele offsets, times the sum of the two allele terms.

```latex
Y_{\mathrm{H}i} \mid z_{\mathrm{Hc}i}, z_{\mathrm{He}i} \sim \mathrm{Poisson}\left(\frac{L_{\mathrm{H}i}}{D}\left[\bar{s}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i} + (1-\bar{s})\,\mu_{\mathrm{He}}\,z_{\mathrm{He}i}\right]\right)
```

We write μ\_sum = s̄μ\_Hc + (1 − s̄)μ\_He for the mean rate of the summed counts and define π as the fraction of it contributed by the Hc allele. Rearranging gives a weighted sum for the expression factor of the summed hybrid, z\_Hi. Each expression factor has mean one and the weights sum to one, so z\_Hi has mean one.

```latex
\begin{aligned}
Y_{\mathrm{H}i} \mid z_{\mathrm{Hc}i}, z_{\mathrm{He}i} &\sim \mathrm{Poisson}\left(\frac{L_{\mathrm{H}i}}{D}\,\mu_{\mathrm{sum}}\left[\pi\,z_{\mathrm{Hc}i} + (1-\pi)\,z_{\mathrm{He}i}\right]\right), \qquad \mu_{\mathrm{sum}} \equiv \bar{s}\,\mu_{\mathrm{Hc}} + (1-\bar{s})\,\mu_{\mathrm{He}}, \qquad \pi = \frac{\bar{s}\,\mu_{\mathrm{Hc}}}{\mu_{\mathrm{sum}}} \\
&= \mathrm{Poisson}\left(\frac{L_{\mathrm{H}i}}{D}\,\mu_{\mathrm{sum}}\,z_{\mathrm{H}i}\right), \qquad z_{\mathrm{H}i} \equiv \pi\,z_{\mathrm{Hc}i} + (1-\pi)\,z_{\mathrm{He}i} \\
\mathrm{E}[z_{\mathrm{H}i}] &= \pi\,\mathrm{E}[z_{\mathrm{Hc}i}] + (1-\pi)\,\mathrm{E}[z_{\mathrm{He}i}] = \pi + (1-\pi) = 1
\end{aligned}
```

An NB fit of the summed counts would treat z\_Hi as its expression factor, and we call the variance of that factor φ\_sum.

**Intrinsic and extrinsic noise.** In the hybrid, each allele has an intrinsic deviation (q) that arises from its own stochastic bursting, and the two alleles share an extrinsic deviation (v) due to shared cell size, cell-cycle state, and pool of transcription factors. The deviations are mutually independent and both deviations have a mean of zero. Since the dispersion is the variance of the expression factor, we have, assuming the two alleles have the same intrinsic variance:

```latex
z_{\mathrm{Hc}i} = 1 + v_i + q_{\mathrm{Hc}i}, \qquad z_{\mathrm{He}i} = 1 + v_i + q_{\mathrm{He}i} \\
\mathrm{Var}[v_i] = \sigma^2_v, \qquad \mathrm{Var}[q_{\mathrm{Hc}i}] = \mathrm{Var}[q_{\mathrm{He}i}] = \sigma^2_q \\
\varphi_{\mathrm{Hc}} = \mathrm{Var}[z_{\mathrm{Hc}i}] = \mathrm{Var}[v_i] + \mathrm{Var}[q_{\mathrm{Hc}i}] = \sigma^2_v + \sigma^2_q = \varphi_{\mathrm{He}}
```

The summed expression factor is a weighted average of the two allele expression factors. The weights sum to one, so the extrinsic deviation keeps a coefficient of one, and the independence of the deviations lets the variances add. This gives the dispersion of the summed hybrid.

```latex
\begin{aligned}
\varphi_{\mathrm{sum}} &= \mathrm{Var}\left[\pi\,z_{\mathrm{Hc}i} + (1-\pi)\,z_{\mathrm{He}i}\right] = \mathrm{Var}\left[1 + v_i + \pi\,q_{\mathrm{Hc}i} + (1-\pi)\,q_{\mathrm{He}i}\right] \\
&= \sigma^2_v + \pi^2\sigma^2_q + (1-\pi)^2\sigma^2_q = \sigma^2_v + \left[\pi^2 + (1-\pi)^2\right]\sigma^2_q
\end{aligned}
```

The extrinsic variance passes through in full, whereas the intrinsic variance enters with weight π² + (1 − π)². Expanding the square and completing it around π = 1/2 gives:

```latex

\pi^2 + (1-\pi)^2 = 1 - 2\pi + 2\pi^2 = \tfrac12 + 2\left(\pi - \tfrac12\right)^2 \\
\tfrac12 \le \pi^2 + (1-\pi)^2 < 1 \qquad \text{for } 0 < \pi < 1

```

Whenever both alleles contribute counts, the weight on the intrinsic variance is between one half and one. It equals one half at equal counts and approaches one as a single allele dominates. The summed hybrid therefore carries less intrinsic variance than a single allele, and this averaging can be estimated for each gene.

**Ploidy shift.** A parent carries one genome, so its expression factor has the same form as that of one allele, with the same extrinsic and intrinsic variances. We write φ\_H for the dispersion of a single genome. The difference between φ\_H and the dispersion of the summed hybrid is the averaging component. Dividing this difference by φ\_H gives the fraction of dispersion lost in the summed hybrid. This fraction depends only on a single factor ω, the intrinsic variance share of φ\_H.

```latex
\begin{aligned}
\varphi_{\mathrm{H}} &= \mathrm{Var}\left[1 + v_i + q_{\mathrm{Hc}i}\right] = \sigma^2_v + \sigma^2_q \\
\varphi_{\mathrm{H}} - \varphi_{\mathrm{sum}} &= \left(\sigma^2_v + \sigma^2_q\right) - \left(\sigma^2_v + \left[\pi^2 + (1-\pi)^2\right]\sigma^2_q\right) = \left[1 - \pi^2 - (1-\pi)^2\right]\sigma^2_q \\
\frac{\varphi_{\mathrm{sum}}}{\varphi_{\mathrm{H}}} &= 1 - \frac{\varphi_{\mathrm{H}} - \varphi_{\mathrm{sum}}}{\varphi_{\mathrm{H}}} = 1 - \left[1 - \pi^2 - (1-\pi)^2\right]\omega, \qquad \omega \equiv \frac{\sigma^2_q}{\sigma^2_v + \sigma^2_q}
\end{aligned}
```

The log₂ of the inverse of this ratio is the expected fall in log₂ dispersion from summing.

```latex
\Lambda \equiv \log_2\frac{\varphi_{\mathrm{H}}}{\varphi_{\mathrm{sum}}} = -\log_2\frac{\varphi_{\mathrm{sum}}}{\varphi_{\mathrm{H}}} = -\log_2\left(1 - \left[1 - \pi^2 - (1-\pi)^2\right]\omega\right)
```

The weight 1 − π² − (1 − π)² equals 2π(1 − π) and lies between zero and one half, and ω lies between zero and one, so the argument of the logarithm lies between one half and one and Λ lies between 0 and 1. The shift equals 1 for two equally expressed alleles with purely intrinsic noise, the case in which averaging halves the noise. It falls to 0 when noise is entirely extrinsic, or when one allele dominates expression and leaves nothing to average.

On the raw scale, the same averaging lowers the dispersion of the summed hybrid below that of a single genome by the difference of the two dispersions. It is additive and proportional to the intrinsic variance, so it is larger in absolute terms for genes with more intrinsic noise.

```latex
\Lambda_{\mathrm{raw}} \equiv \varphi_{\mathrm{H}} - \varphi_{\mathrm{sum}} = \left[1 - \pi^2 - (1-\pi)^2\right]\sigma^2_q
```

**Estimating the ploidy shift.** We estimated the ploidy shift from estimators of π, σ²\_v and σ²\_q. For the weight π, we used the count share of Hc in the summed allele counts, whose expectation is π because the summed counts carry the library-weighted mean share.

```latex
\hat{\pi} = \frac{\sum_i y_{\mathrm{Hc}i}}{\sum_i \left(y_{\mathrm{Hc}i} + y_{\mathrm{He}i}\right)}
```

This estimator is a ratio of sums, so it is consistent but biased. The bias comes from the joint fluctuation of the two sums around their means. Each sum runs over N\_H cells, so the bias is 𝒪(1/N\_H), and we dropped it.

To estimate the extrinsic variance, we note that it is also the covariance of the two alleles' expression factors, because the extrinsic deviation is the only part that they share. Each expression factor has mean one, so expanding the covariance gives:

```latex
\begin{aligned}
\sigma^2_v &= \mathrm{Var}[v_i] = \mathrm{Cov}[z_{\mathrm{Hc}i}, z_{\mathrm{He}i}] = \mathrm{E}[z_{\mathrm{Hc}i}\,z_{\mathrm{He}i}] - \mathrm{E}[z_{\mathrm{Hc}i}]\,\mathrm{E}[z_{\mathrm{He}i}] = \mathrm{E}[z_{\mathrm{Hc}i}\,z_{\mathrm{He}i}] - 1
\end{aligned}
```

The intrinsic variance is the same for both alleles and the expected deviations have mean zero, so expanding the variances gives:

```latex
\begin{aligned}
\sigma^2_q &= \mathrm{Var}[q_{\mathrm{Hc}i}] = \tfrac12\left[\mathrm{Var}[q_{\mathrm{Hc}i}] + \mathrm{Var}[q_{\mathrm{He}i}]\right] \\
&= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] - \mathrm{E}[q_{\mathrm{Hc}i}]^2 + \mathrm{E}[q_{\mathrm{He}i}^2] - \mathrm{E}[q_{\mathrm{He}i}]^2\right) 
= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] + \mathrm{E}[q_{\mathrm{He}i}^2]\right) 
\end{aligned}
```

The two intrinsic deviations are independent, so E\[q\_Hc q\_He\] is zero, and subtracting twice this cross moment completes the square. Within the square, we add and subtract one plus the extrinsic deviation, and then substitute with the allele expression factors. Expanding and simplifying gives:

```latex
\begin{aligned}
\sigma^2_q &= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] + \mathrm{E}[q_{\mathrm{He}i}^2] - 2\,\mathrm{E}[q_{\mathrm{Hc}i}\,q_{\mathrm{He}i}]\right) = \tfrac12\,\mathrm{E}\left[\left(q_{\mathrm{Hc}i} - q_{\mathrm{He}i}\right)^2\right]\\
&= \tfrac12\,\mathrm{E}\left[\left((1 + v_i + q_{\mathrm{Hc}i}) - (1 + v_i + q_{\mathrm{He}i})\right)^2\right] = \tfrac12\,\mathrm{E}\left[\left(z_{\mathrm{Hc}i} - z_{\mathrm{He}i}\right)^2\right] \\
&= \tfrac12\,\mathrm{E}\left[z_{\mathrm{Hc}i}^2 - 2\,z_{\mathrm{Hc}i}\,z_{\mathrm{He}i} + z_{\mathrm{He}i}^2\right] = \tfrac12\left(\mathrm{E}[z_{\mathrm{Hc}i}^2] + \mathrm{E}[z_{\mathrm{He}i}^2]\right) - \mathrm{E}[z_{\mathrm{Hc}i}\,z_{\mathrm{He}i}]
\end{aligned}
```

This leaves the intrinsic and extrinsic variances written as moments of the expression factors. To simplify, we multiply and divide each expression factor by its allele offset and per-genome rate, and replace each product with the conditional expected count. Because the two counts are independent given their expression factors, the product of the conditional expectations is the conditional expectation of the product. The law of total expectation then removes the conditioning.

```latex
\begin{aligned}
\mathrm{E}[z_{\mathrm{Hc}i}\,z_{\mathrm{He}i}] &= \mathrm{E}\left[\frac{O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}} \cdot \frac{O_{\mathrm{He}i}\,\mu_{\mathrm{He}}\,z_{\mathrm{He}i}}{O_{\mathrm{He}i}\,\mu_{\mathrm{He}}}\right] = \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i} \mid z_{\mathrm{Hc}i}]}{O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}} \cdot \frac{\mathrm{E}[Y_{\mathrm{He}i} \mid z_{\mathrm{He}i}]}{O_{\mathrm{He}i}\,\mu_{\mathrm{He}}}\right] \\
&= \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i} \mid z_{\mathrm{Hc}i}, z_{\mathrm{He}i}]}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] = \mathrm{E}\left[\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\,\middle|\, z_{\mathrm{Hc}i}, z_{\mathrm{He}i}\right]\right] = \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right]
\end{aligned}
```

For the square moment, we write z\_Hc in terms of the conditional Poisson mean of its count. A Poisson count has variance equal to its mean, so its conditional mean square is the mean plus the squared mean. Subtracting the conditional mean leaves the squared mean. The denominator is fixed, and the law of total expectation removes the conditioning, giving:

```latex
\begin{aligned}
\mathrm{E}[z_{\mathrm{Hc}i}^2] &= \mathrm{E}\left[\frac{\left(O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i}\right)^2}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\right] = \mathrm{E}\left[\frac{\left(O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i} + \left(O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i}\right)^2\right) - O_{\mathrm{Hc}i}\,\mu_{\mathrm{Hc}}\,z_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\right] = \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1) \mid z_{\mathrm{Hc}i}]}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\right] \\
&= \mathrm{E}\left[\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\,\middle|\, z_{\mathrm{Hc}i}\right]\right] = \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\right], \qquad \mathrm{E}[z_{\mathrm{He}i}^2] = \mathrm{E}\left[\frac{Y_{\mathrm{He}i}(Y_{\mathrm{He}i} - 1)}{O_{\mathrm{He}i}^2\,\mu_{\mathrm{He}}^2}\right]
\end{aligned}
```

Substituting these moments into the definitions gives σ²\_v and σ²\_q in terms of the counts, the allele offsets and the per-genome rates. We then replaced each expectation by the mean over cells and each expected count by its estimate. These estimators take the allele offsets as exact and set aside the error of the library split, which biases them, so we marked them with a star and corrected them.

```latex
\begin{aligned}
\sigma^2_v &= \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] - 1 \quad\Rightarrow\quad \hat{\sigma}^{2*}_v = \frac{1}{N_{\mathrm{H}}}\sum_i \frac{y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\hat{\mu}_{\mathrm{Hc}}\,\hat{\mu}_{\mathrm{He}}} - 1 \\
\sigma^2_q &= \tfrac12\left(\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{O_{\mathrm{Hc}i}^2\,\mu_{\mathrm{Hc}}^2}\right] + \mathrm{E}\left[\frac{Y_{\mathrm{He}i}(Y_{\mathrm{He}i} - 1)}{O_{\mathrm{He}i}^2\,\mu_{\mathrm{He}}^2}\right]\right) - \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] \\
\hat{\sigma}^{2*}_q &= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}(y_{\mathrm{Hc}i} - 1)}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} + \frac{y_{\mathrm{He}i}(y_{\mathrm{He}i} - 1)}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2}\right] - \frac{1}{N_{\mathrm{H}}}\sum_i \frac{y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\hat{\mu}_{\mathrm{Hc}}\,\hat{\mu}_{\mathrm{He}}}
\end{aligned}
```

Expanding the factorial moments and completing the square in the scaled allele counts rewrites the estimator of the intrinsic variance.

```latex
\begin{aligned}
\hat{\sigma}^{2*}_q 
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}^2}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} - \frac{y_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} + \frac{y_{\mathrm{He}i}^2}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2} - \frac{y_{\mathrm{He}i}}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2} - \frac{2\,y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\hat{\mu}_{\mathrm{Hc}}\,\hat{\mu}_{\mathrm{He}}}\right] \\
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}^2}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} - \frac{2\,y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_{\mathrm{Hc}i}\,O_{\mathrm{He}i}\,\hat{\mu}_{\mathrm{Hc}}\,\hat{\mu}_{\mathrm{He}}} + \frac{y_{\mathrm{He}i}^2}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2}\right] - \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} + \frac{y_{\mathrm{He}i}}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2}\right) \\
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}\,\hat{\mu}_{\mathrm{Hc}}} - \frac{y_{\mathrm{He}i}}{O_{\mathrm{He}i}\,\hat{\mu}_{\mathrm{He}}}\right)^2 - \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_{\mathrm{Hc}i}^2\,\hat{\mu}_{\mathrm{Hc}}^2} + \frac{y_{\mathrm{He}i}}{O_{\mathrm{He}i}^2\,\hat{\mu}_{\mathrm{He}}^2}\right)
\end{aligned}
```

The estimator of the intrinsic variance is therefore half the mean squared difference of the scaled allele counts minus their average Poisson noise. The estimator of the extrinsic variance is the mean product of the scaled allele counts minus one. Each scaled count has expectation one, so this estimates their covariance, and it needs no Poisson correction because the two counts are independent given their expression factors.

With the true mean rates, every term is an unbiased estimate of the corresponding moment of the expression factors, so both estimators are unbiased. As for π, the fitted mean rates are ratios of sums. Replacing each mean rate by its fitted value changes each term with a bias of order 1/N\_H, which we dropped.

```latex
\mathrm{E}\left[\hat{\sigma}^{2*}_v\right] - \sigma^2_v = \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right), \qquad \mathrm{E}\left[\hat{\sigma}^{2*}_q\right] - \sigma^2_q = \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right)
```

**Split correction.** The star estimators take the allele offsets as exact. Each is built from cell averages that estimate the second moments of the expression factors, α\_Hc ≡ E\[z²\_Hc\] = 1 + φ\_H, α\_He ≡ E\[z²\_He\] = 1 + φ\_H and α\_× ≡ E\[z\_Hc z\_He\] = 1 + σ²\_v, with σ̂²\*\_v = α̂\*\_× − 1 and σ̂²\*\_q = ½(α̂\*\_Hc + α̂\*\_He) − α̂\*\_×. Their limits under the split error are derived in *Library split and the estimated shift*. Dividing each cell average by its split factor and recombining gives the corrected estimators, which we truncated at zero. The bar marks the plain mean over hybrid cells:

```latex
\begin{aligned}
\hat{\alpha}_{\mathrm{Hc}} &= \frac{\hat{\alpha}^*_{\mathrm{Hc}}}{1 + 3\bar{\eta}_{\mathrm{Hc}}}, \qquad \hat{\alpha}_{\mathrm{He}} = \frac{\hat{\alpha}^*_{\mathrm{He}}}{1 + 3\bar{\eta}_{\mathrm{He}}}, \qquad \hat{\alpha}_{\times} = \frac{\hat{\alpha}^*_{\times}}{1 + \bar{\eta}_{\times}} \\
\hat{\sigma}^2_v &= \hat{\alpha}_{\times} - 1, \qquad \hat{\sigma}^2_q = \tfrac12\left(\hat{\alpha}_{\mathrm{Hc}} + \hat{\alpha}_{\mathrm{He}}\right) - \hat{\alpha}_{\times}
\end{aligned}
```

Plugging these estimates into the definition of Λ̂ gives the log₂ ratio of the estimated dispersion of a single genome to the estimated dispersion of the summed hybrid.

```latex
\begin{aligned}
\hat{\Lambda} &\equiv \log_2\frac{\hat{\sigma}^2_v + \hat{\sigma}^2_q}{\hat{\sigma}^2_v + \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q}
\end{aligned}
```

Both truncated variances are nonnegative, so Λ̂ lies between 0 and 1. A zero estimate of the intrinsic variance gives Λ̂ = 0. When both variance estimates are zero the ratio is undefined, and we set Λ̂ to zero. Truncation therefore affects Λ̂ only for genes with a variance estimate near zero.

The estimated raw shift Λ̂\_raw is the difference between the estimated dispersion of a single genome and that of the summed hybrid.

```latex
\begin{aligned}
\hat{\Lambda}_{\mathrm{raw}} &\equiv \left(\hat{\sigma}^2_v + \hat{\sigma}^2_q\right) - \left(\hat{\sigma}^2_v + \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q\right) \\
&= \hat{\sigma}^2_v - \hat{\sigma}^2_v + \hat{\sigma}^2_q - \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q = \left[1 - \hat{\pi}^2 - (1-\hat{\pi})^2\right]\hat{\sigma}^2_q
\end{aligned}
```

**Library split and the estimated shift.** The split error biases the estimators of σ²\_v and σ²\_q, and through them the estimated shifts. The weight estimate π̂ is a ratio of summed counts, which the library split does not enter, so π̂ → π.

Each star estimator divides the counts of a cell by its allele offsets, which are the true allele offsets times 1 + ξ (*Library split correction*). Given the expression factors, the split error is independent of the counts, so each cell average converges to a moment of the expression factors times a moment of the split errors. Expanding 1/(1 + ξ) to second order, with E\[ξ\] = 0, the He error ξ\_Hei = −s\_i ξ\_Hci/(1 − s\_i), because He takes the rest of the library, its variance η\_Hei = s\_i/\[(1 − s\_i)L\_Hi\], and the cross moment E\[ξ\_Hci ξ\_Hei\] = −s\_i η\_Hci/(1 − s\_i) = −1/L\_Hi, gives the split factors:

```latex
\begin{aligned}
\mathrm{E}\left[\frac{1}{\left(1+\xi_{\mathrm{Hc}i}\right)^2}\right] &\approx \mathrm{E}\left[1 - 2\xi_{\mathrm{Hc}i} + 3\xi_{\mathrm{Hc}i}^2\right] = 1 + 3\eta_{\mathrm{Hc}i} \\
\mathrm{E}\left[\frac{1}{\left(1+\xi_{\mathrm{Hc}i}\right)\left(1+\xi_{\mathrm{He}i}\right)}\right] &\approx \mathrm{E}\left[1 - \xi_{\mathrm{Hc}i} - \xi_{\mathrm{He}i} + \xi_{\mathrm{Hc}i}^2 + \xi_{\mathrm{Hc}i}\,\xi_{\mathrm{He}i} + \xi_{\mathrm{He}i}^2\right] = 1 + \eta_{\mathrm{Hc}i} + \eta_{\mathrm{He}i} - \frac{1}{L_{\mathrm{H}i}} \equiv 1 + \eta_{\times i}
\end{aligned}
```

The moments of the expression factors are the same in every cell, so averaging over hybrid cells gives the limits of the cell averages:

```latex
\hat{\alpha}^*_{\mathrm{Hc}} \to \alpha_{\mathrm{Hc}}\left(1 + 3\bar{\eta}_{\mathrm{Hc}}\right), \qquad \hat{\alpha}^*_{\mathrm{He}} \to \alpha_{\mathrm{He}}\left(1 + 3\bar{\eta}_{\mathrm{He}}\right), \qquad \hat{\alpha}^*_{\times} \to \alpha_{\times}\left(1 + \bar{\eta}_{\times}\right)
```

Solving these limits for the moments gives the corrected estimators of *Split correction*, with ŝ\_i in place of s\_i in each η. Substituting the limits into the star estimators shows that the split error raises both of them. The extrinsic estimator rises by α\_× η̄\_×, and the intrinsic estimator by 3/2 times the square-moment terms less α\_× η̄\_×, so the correction matters most for the intrinsic variance:

```latex
\hat{\sigma}^{2*}_v \to \sigma^2_v + \alpha_{\times}\,\bar{\eta}_{\times}, \qquad \hat{\sigma}^{2*}_q \to \sigma^2_q + \tfrac32\left(\alpha_{\mathrm{Hc}}\,\bar{\eta}_{\mathrm{Hc}} + \alpha_{\mathrm{He}}\,\bar{\eta}_{\mathrm{He}}\right) - \alpha_{\times}\,\bar{\eta}_{\times}
```

## Appendix E. State-conditioned estimates

This appendix gives an estimate of every contrast that is conditioned on one cell state. It flags the significant contrasts that the state, rather than the genotype, carries.

**Why the culture matters.** We grew the two parents together in one culture and the hybrid in a second culture. Every contrast between parental and hybrid cells, which includes the *trans* and dominance contrasts, therefore also carries any difference between the two cultures, and the design cannot separate that difference from the genotype. The total contrast compares the two parents within one culture, and the *cis* contrast compares two alleles within the same cells, so neither carries the culture difference. One culture difference is visible in the cells. A program of glucose derepression and of the induced stress response (Gasch et al. 2000) varies among cells. Its stress-response score, standardized with depth taken into account, exceeds 2.5 in 4.0% of *S. cerevisiae* parental cells and 5.9% of *S. eubayanus* parental cells, but in only 0.25% of hybrid cells. Its derepression score exceeds 4 in 3 of 834 and 6 of 1240 parental cells and in none of 2791 hybrid cells (check ST1, analysis/checks/state\_groups.R). In simulations of a state present in 0.1% to 0.3% of cells with a 50-fold expression change, the type I error of the dispersion test at the 5% level averaged 0.12 for total and 0.08 to 0.10 for *trans*, and reached 0.155, because a group that samples no state cells cannot show them in its bootstrap. The *cis* test averaged at most 0.06 across all settings and reached 0.09, because both alleles of a cell share its state (check M3, 200 datasets per setting). We therefore kept the unconditioned estimates and tests as the main analysis, refitted every contrast with the program in the mean, and flagged the calls that the program carries.

**Program score.** Twenty marker genes, the set G, define the program. Ten mark glucose derepression (ICL1, MLS1, CIT2, MDH2, ACS1, FOX2, ADH2, IDP2, JEN1, ALD4) and ten the induced stress response (HSP12, HSP26, HSP42, HSP104, CTT1, TPS1, TPS2, GSY2, GRE2, SSA4). We chose them before any noise estimate, as the first members of ordered candidate lists with a one-to-one ortholog in both species. For gene g in cell i, the score S\_i is the log₂ of the summed marker counts per unit offset, with gene g's own count left out so that no gene enters its own covariate. In a hybrid cell the sum runs over both alleles, so both alleles of a cell share one score. We centered the score at its median within each culture, separately for the two parents, because the marker genes differ in expression between the species:

```latex
S_i = \log_2\frac{\tfrac12 + \sum_{g' \in G,\, g' \neq g} y_{g'i}}{O_i} - \mathrm{median}_{i'}\left(\log_2\frac{\tfrac12 + \sum_{g' \in G,\, g' \neq g} y_{g'i'}}{O_{i'}}\right)
```

The median runs over the cells i′ of the same parent, or over all hybrid cells. The pseudocount ½ keeps the score finite in a cell without marker counts. The offset of a hybrid cell's score is that of the whole cell, L\_Hi/D, and any constant factor in the offset cancels in the centering. A score of zero marks the typical cell of its culture, and each unit of S\_i doubles the program's expression. We computed S\_i once from the observed counts and held it fixed.

**Conditioned mean.** The mean count of a cell now depends on its score through a log-linear term with slope τ\_r, the log₂ change in the gene's mean for each log₂ unit of the program:

```latex
\mathrm{E}\left[Y_{ri} \mid S_i\right] = M_i = O_i\,\mu_r\,2^{\tau_r S_i}, \qquad \mathrm{Var}\left[Y_{ri} \mid S_i\right] = M_i + \varphi_r\,M_i^2
```

Here μ\_r is the mean rate of a cell at S\_i = 0, the typical state of its culture, and φ\_r is the dispersion around the conditioned mean. We estimated μ\_r and τ\_r from the two Poisson score equations, which extend the estimating equation of μ̂\_r in *Estimating μ* by one equation for τ\_r:

```latex
\begin{aligned}
\sum_i \left(y_{ri} - \hat{M}_i\right) &= 0, \qquad \sum_i S_i\left(y_{ri} - \hat{M}_i\right) = 0, \qquad \hat{M}_i = O_i\,\hat{\mu}_r\,2^{\hat{\tau}_r S_i} \\
\hat{\tau}_r = 0 \quad &\Rightarrow\quad \sum_i \left(y_{ri} - O_i\,\hat{\mu}_r\right) = 0 \quad\Rightarrow\quad \hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i}
\end{aligned}
```

These equations assume only the form of the mean, so they are the quasi-Poisson estimating equations (Wedderburn 1974, McCullagh and Nelder 1989), and their solution contains no φ. Without the covariate they give the ratio estimate of the main text. We solved them by Newton steps from the ratio estimate and τ\_r = 0, halving each step that lowered the Poisson log-likelihood. The allele means and the per-genome hybrid mean then follow from μ̂\_r as in the main text. In 300 simulated groups of 500 cells with τ\_r = 0.8 and φ = 0.5, τ̂\_r averaged 0.798 (SD 0.080).

**Conditioned dispersion.** The mean now uses two parameters, so the Cox–Reid adjustment uses the information matrix of both (Cox and Reid 1987, McCarthy et al. 2012). We write ℐ\_0, ℐ\_1 and ℐ\_2 for the weighted sums of the zeroth, first and second powers of the score. Up to constant factors of ln 2 that do not move the maximum, they are the entries of the information matrix of log μ\_r and τ\_r, and its determinant is ℐ\_0ℐ\_2 − ℐ\_1²:

```latex
\begin{aligned}
\ell^{\mathrm{CR}}_r(\varphi) &\equiv \ell_r(\varphi) - \tfrac12\log\left[\mathcal{I}_0(\varphi)\,\mathcal{I}_2(\varphi) - \mathcal{I}_1(\varphi)^2\right], \qquad \mathcal{I}_0(\varphi) = \sum_i \frac{\hat{M}_i}{1 + \varphi\hat{M}_i}, \quad \mathcal{I}_1(\varphi) = \sum_i \frac{S_i\hat{M}_i}{1 + \varphi\hat{M}_i}, \quad \mathcal{I}_2(\varphi) = \sum_i \frac{S_i^2\hat{M}_i}{1 + \varphi\hat{M}_i} \\
\hat{\tau}_r = 0 \ \text{fixed} \quad &\Rightarrow\quad \ell^{\mathrm{CR}}_r(\varphi) = \ell_r(\varphi) - \tfrac12\log\mathcal{I}_0(\varphi) = \ell_r(\varphi) - \tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}
\end{aligned}
```

Here ℓ\_r(φ) is the NB log-likelihood with every cell's mean fixed at M̂\_i. With the slope fixed at zero, the adjustment returns the one-parameter form of *Estimating φ*. A score with no variance among the cells of a group carries no information, so we fitted that group without it.

Whether the maximum lies at zero again depends on the adjusted slope at zero. The slope of ℓ\_r at zero is the sum of the per-cell slopes of *Estimating φ*, and the first score equation sets Σ\_i y\_ri equal to Σ\_i M̂\_i. At zero, ℐ\_0, ℐ\_1 and ℐ\_2 equal Σ\_i M̂\_i, Σ\_i S\_i M̂\_i and Σ\_i S\_i² M̂\_i, and each falls with φ at the rate of the same sum with M̂\_i² in place of M̂\_i:

```latex
\begin{aligned}
\ell_r'(0) &= \tfrac12\left[\sum_i \left(y_{ri} - \hat{M}_i\right)^2 - \sum_i y_{ri}\right] = \tfrac12\left[\sum_i \left(y_{ri} - \hat{M}_i\right)^2 - \sum_i \hat{M}_i\right] \\
\mathcal{I}_0'(0) &= -\sum_i \hat{M}_i^2, \qquad \mathcal{I}_1'(0) = -\sum_i S_i\hat{M}_i^2, \qquad \mathcal{I}_2'(0) = -\sum_i S_i^2\hat{M}_i^2 \\
\left.\frac{\partial}{\partial \varphi}\left(-\tfrac12\log\left[\mathcal{I}_0\mathcal{I}_2 - \mathcal{I}_1^2\right]\right)\right|_{\varphi = 0} &= \frac{\mathcal{I}_2(0)\sum_i \hat{M}_i^2 + \mathcal{I}_0(0)\sum_i S_i^2\hat{M}_i^2 - 2\,\mathcal{I}_1(0)\sum_i S_i\hat{M}_i^2}{2\left[\mathcal{I}_0(0)\,\mathcal{I}_2(0) - \mathcal{I}_1(0)^2\right]}
\end{aligned}
```

The estimator therefore keeps the two cases of the main text, with the adjustment's slope in place of Σ\_i M̂\_i²/Σ\_i M̂\_i:

```latex
\hat{\varphi}_r = \begin{cases} 0 & \sum_i \left(y_{ri} - \hat{M}_i\right)^2 \le \sum_i \hat{M}_i - \dfrac{\mathcal{I}_2(0)\sum_i \hat{M}_i^2 + \mathcal{I}_0(0)\sum_i S_i^2\hat{M}_i^2 - 2\,\mathcal{I}_1(0)\sum_i S_i\hat{M}_i^2}{\mathcal{I}_0(0)\,\mathcal{I}_2(0) - \mathcal{I}_1(0)^2} \\ \arg\max_{\varphi > 0}\ \ell^{\mathrm{CR}}_r(\varphi) & \text{otherwise} \end{cases}
```

We applied the split correction of the hybrid alleles to the conditioned dispersions unchanged, which assumes that the transcriptome share of a cell is independent of its program score. In simulations with a gamma expression factor and a log-linear score effect (600 cells, φ = 0.3, 200 replicates), the conditioned estimate averaged 0.300, while the unconditioned estimate averaged 0.82 because it absorbed the program's variance. With Poisson counts and a score effect (300 replicates), 48% of conditioned estimates fell on the boundary, close to the 50% expected for a parameter on the boundary.

**Contrasts, bootstrap and flags.** We built every contrast of every quantity from the conditioned fits as in *Contrasts*, so the log₂ and raw dispersion, burst frequency, burst size and kinetic balance contrasts follow from μ̂\_r and φ̂\_r exactly as before. Each bootstrap replicate used the resamples of the reported bootstrap and refitted every group twice, with and without the score, with the split correction recomputed on the resampled hybrid cells. A resampled cell kept its score. This gave each contrast a conditioned estimate Δ^S with a 95% percentile interval, and the paired change Δ − Δ^S between the unconditioned and conditioned contrasts on the same cells, with its own percentile interval. We built no null for the conditioned contrasts, so they carry intervals but no p-values. We flagged a contrast with a significant unconditioned test (adjusted p-value below 0.05, Multiple testing) when its conditioned interval contains zero and the interval of its change excludes zero. For the *trans* and dominance contrasts a flag marks an effect confounded with the culture. For the total and *cis* contrasts, whose groups share a culture, it marks an effect of the program itself. The figures draw flagged genes as open triangles.

**What conditioning removes.** Conditioning removes the program's contribution to each group's dispersion, and it evaluates each group's mean at the typical state of its culture, so the share of cells in the program no longer moves the mean contrasts. It does not remove a difference between the cultures in their typical program level, which the centering assigns to each culture's own median. It also leaves every culture difference that the score does not measure. The score is measured with counting noise, so the conditioning is partial and the conditioned dispersion stays above the dispersion given the true state. The marker genes keep their tests, and their conditioned dispersion measures the noise beyond what the other markers explain.

## Appendix F. Correlation of mean and dispersion errors on shared cells

When the mean and dispersion of a group are fitted on the same cells, as in the kinetic balance and burst size partitions, their errors can be correlated, so we asked when. Within one group, let the mean be fitted on a set of cells T and the dispersion on a set U, which may share cells. The error of the fitted mean is a sum of count deviations over T. When the counts follow the NB, the fitted dispersion sets the slope of the log-likelihood in φ to zero (*Estimating φ*), so to first order its error is proportional to the sum over U of the per-cell slopes ℓ′\_i, evaluated at the true mean and dispersion. Writing out the covariance of the two errors with these definitions, the expected counts M\_i are constants, a count and a slope from different cells are independent, and only cells in both sets remain:

```latex
\begin{aligned}
\hat{\mu}_T - \mu &= \frac{\sum_{i \in T}\left(Y_i - M_i\right)}{\sum_{i \in T} O_i}, \qquad \hat{\varphi}_U - \varphi \propto \sum_{j \in U} \ell'_j \\
\mathrm{Cov}\left[\hat{\mu}_T - \mu,\ \hat{\varphi}_U - \varphi\right] &\propto \mathrm{Cov}\left[\sum_{i \in T}\left(Y_i - M_i\right),\ \sum_{j \in U} \ell'_j\right] = \sum_{i \in T}\sum_{j \in U}\mathrm{Cov}\left[Y_i, \ell'_j\right] \\
&= \sum_{i \in T \cap U}\mathrm{Cov}\left[Y_i, \ell'_i\right], \qquad \text{since } \mathrm{Cov}\left[Y_i, \ell'_j\right] = 0 \text{ for } i \ne j
\end{aligned}
```

Each remaining term is zero when the counts follow the NB. A slope has mean zero, so the covariance is the mean of the count times its slope. Differentiating the sum of y times its probability turns this mean into the derivative of the expected count with respect to φ, and the expected count does not depend on φ. The step sums over the NB probability, so it holds only when the counts follow the NB:

```latex
\mathrm{Cov}\left[Y_i, \ell'_i\right] = \mathrm{E}\left[Y_i\,\ell'_i\right] = \sum_y y\,\frac{\partial \Pr(y \mid M_i, \varphi)}{\partial \varphi} = \frac{\partial}{\partial \varphi}\sum_y y\,\Pr(y \mid M_i, \varphi) = \frac{\partial M_i}{\partial \varphi} = 0
```

The mean and dispersion errors of one group are therefore uncorrelated when the two sets share no cells, whatever the distribution of the counts. When the sets share cells, the errors are uncorrelated only if the counts follow the NB, which requires a gamma expression factor. The factor of a hybrid allele is not gamma in general, because the allele offset carries the split error and the expression factor has a shared and a private part, so shared cells can correlate the two errors.

To size this covariance, we used the moment estimator of Appendix C, whose first-order error is explicit, with equal offsets scaled to one so that every cell has the same weight and M\_i = μ. Let z\_i be the expression factor of cell i, with mean one and variance φ, so that given z\_i the count Y\_i is Poisson with mean μ z\_i. The estimator averages \[(Y\_j − μ̂\_U)² − μ̂\_U\]/μ̂\_U² over U. As a function of the fitted mean μ̂\_U, the expected value of this bracket is \[μ + μ²φ + (μ − μ̂\_U)² − μ̂\_U\]/μ̂\_U², and its derivative with respect to μ̂\_U at μ̂\_U = μ is −(1/μ² + 2φ/μ). Expanding the estimator around the true mean gives its error. A count in T and a term in U covary only through the cells in both sets, and the fitted mean of U covaries with the mean of T through the variance of one count over those cells. We write N\_T, N\_U and N\_{T∩U} for the numbers of cells in T, in U and in both. Because E\[Y\_i − μ\] = 0, the covariance of a count with the bracket is its third central moment divided by μ², and the variance of a count is μ(1 + μφ) for any mean-one factor with variance φ:

```latex
\begin{aligned}
\hat{\varphi}_U - \varphi &\approx \frac{1}{N_U}\sum_{j \in U}\left[\frac{(Y_j - \mu)^2 - \mu}{\mu^2} - \varphi\right] - \left(\frac{1}{\mu^2} + \frac{2\varphi}{\mu}\right)\left(\hat{\mu}_U - \mu\right) \\
\mathrm{Cov}\left[\hat{\mu}_T - \mu,\ \hat{\varphi}_U - \varphi\right] &\approx \frac{N_{T \cap U}}{N_T N_U}\left\{\frac{\mathrm{E}\left[(Y_i - \mu)^3\right]}{\mu^2} - \left(\frac{1}{\mu^2} + \frac{2\varphi}{\mu}\right)\mu\left(1 + \mu\varphi\right)\right\}
\end{aligned}
```

The third central moment follows from writing Y\_i − μ as the Poisson deviation Y\_i − μ z\_i plus the factor deviation μ(z\_i − 1). Given z\_i, the Poisson deviation has mean zero and its second and third moments equal μ z\_i. Expanding the cube leaves three terms, because the term with one Poisson deviation and the square of the factor deviation has conditional mean zero, and E\[z\_i²\] = 1 + φ gives E\[z\_i(z\_i − 1)\] = φ:

```latex
\begin{aligned}
\mathrm{E}\left[(Y_i - \mu)^3\right] &= \mathrm{E}\left[(Y_i - \mu z_i)^3\right] + 3\,\mathrm{E}\left[(Y_i - \mu z_i)^2\,\mu(z_i - 1)\right] + \mu^3\,\mathrm{E}\left[(z_i - 1)^3\right] \\
&= \mu + 3\mu^2\,\mathrm{E}\left[z_i(z_i - 1)\right] + \mu^3\,\mathrm{E}\left[(z_i - 1)^3\right] \\
&= \mu + 3\mu^2\varphi + \mu^3\,\mathrm{E}\left[(z_i - 1)^3\right]
\end{aligned}
```

Substituting this moment and expanding the product in the second term leaves a covariance that depends on the third central moment of the expression factor alone:

```latex
\begin{aligned}
\mathrm{Cov}\left[\hat{\mu}_T - \mu,\ \hat{\varphi}_U - \varphi\right] &\approx \frac{N_{T \cap U}}{N_T N_U}\left\{\frac{1}{\mu} + 3\varphi + \mu\,\mathrm{E}\left[(z_i - 1)^3\right] - \left(\frac{1}{\mu} + 3\varphi + 2\mu\varphi^2\right)\right\} \\
&= \frac{N_{T \cap U}}{N_T N_U}\,\mu\left(\mathrm{E}\left[(z_i - 1)^3\right] - 2\varphi^2\right)
\end{aligned}
```

A gamma factor with mean one and variance φ has third central moment 2φ², so the covariance vanishes, as the slope argument above showed. Any other factor leaves a covariance that shrinks as 1/N like the variances of the two estimates, so the correlation depends on the overlap N\_{T∩U}/√(N\_T N\_U) and does not fall as the sample grows. Kinetic balance and burst size contrasts fit the mean and dispersion on the same cells and carry this correlation. It adds a covariance term to the standard errors derived in *Burst kinetics*, and the bootstrap standard errors of both contrasts include it. We derived the covariance for the moment estimator with equal offsets. The NB maximum likelihood fit depends on the third moment in the same way, through the covariance of a count with its slope, which no longer equals ∂M\_i/∂φ, but we did not derive its coefficient.

## References

Andrews DWK (2000) Inconsistency of the bootstrap when a parameter is on the boundary of the parameter space. Econometrica 68:399–405. https://doi.org/10.1111/1468-0262.00114

Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Benjamini Y, Yekutieli D (2005) False discovery rate–adjusted multiple confidence intervals for selected parameters. Journal of the American Statistical Association 100:71–81. https://doi.org/10.1198/016214504000001907

Cox DR, Reid N (1987) Parameter orthogonality and approximate conditional inference. Journal of the Royal Statistical Society, Series B 49:1–39. https://doi.org/10.1111/j.2517-6161.1987.tb01422.x

Elowitz MB, Levine AJ, Siggia ED, Swain PS (2002) Stochastic gene expression in a single cell. Science 297:1183–1186. https://doi.org/10.1126/science.1070919

Fraser HB (2019) Improving estimates of compensatory *cis*-*trans* regulatory divergence. Trends in Genetics 35:3–5. https://doi.org/10.1016/j.tig.2018.09.003

Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown PO (2000) Genomic expression programs in the response of yeast cells to environmental changes. Molecular Biology of the Cell 11:4241–4257. https://doi.org/10.1091/mbc.11.12.4241

Greenwood M, Yule GU (1920) An inquiry into the nature of frequency distributions representative of multiple happenings with particular reference to the occurrence of multiple attacks of disease or of repeated accidents. Journal of the Royal Statistical Society 83:255–279. https://doi.org/10.2307/2341080

Lawless JF (1987) Negative binomial and mixed Poisson regression. Canadian Journal of Statistics 15:209–225. https://doi.org/10.2307/3314912

McCarthy DJ, Chen Y, Smyth GK (2012) Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40:4288–4297. https://doi.org/10.1093/nar/gks042

McCullagh P, Nelder JA (1989) Generalized Linear Models, 2nd edn. Chapman and Hall, London. ISBN 0-412-31760-5. https://doi.org/10.1007/978-1-4899-3242-6

Phipson B, Smyth GK (2010) Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology 9:39. https://doi.org/10.2202/1544-6115.1585

Raj A, Peskin CS, Tranchina D, Vargas DY, Tyagi S (2006) Stochastic mRNA synthesis in mammalian cells. PLoS Biology 4:e309. https://doi.org/10.1371/journal.pbio.0040309

Schäfer J, Strimmer K (2005) A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical Applications in Genetics and Molecular Biology 4:32. https://doi.org/10.2202/1544-6115.1175

Wedderburn RWM (1974) Quasi-likelihood functions, generalized linear models, and the Gauss-Newton method. Biometrika 61:439–447. https://doi.org/10.1093/biomet/61.3.439

Wittes J, Brittain E (1990) The role of internal pilot studies in increasing the efficiency of clinical trials. Statistics in Medicine 9:65–72. https://doi.org/10.1002/sim.4780090113

Wittkopp PJ, Haerum BK, Clark AG (2004) Evolutionary changes in *cis* and *trans* gene regulation. Nature 430:85–88. https://doi.org/10.1038/nature02698

## Notes for revision

Kept in the document (Notes for revision tab, https://claude.ai/artifact/76UUgMeP9KaXjjcBGwMQa7); docs/validation_checks.md records the code items and their status.
