<!-- Exported from the Claude Docs document "Statistical models for allele-specific expression, noise, and ploidy"
     (https://claude.ai/artifact/79JHPV3XCH2K8VS7RJZQHS), revision 17945, retrieved 2026-11-03. Display equations are LaTeX in $$ blocks.
     This file is a snapshot for the repository; the document remains the specification and is edited there. -->

# Statistical models for allele-specific expression, noise, and ploidy

2026-10-20 · Brian Metzger

We compared expression and noise between groups of cells through quantities fitted to each gene. To put cells sequenced to different depths on a common scale, we used a depth offset, and then estimated a mean rate and a dispersion for each gene with a negative binomial model of the counts. We summarized burst kinetics with kinetic balance rather than burst size. In the hybrid, we removed the effects of cell-to-cell variation in the transcriptome share of each allele on the mean and the dispersion. We then formed contrasts between groups on the log₂ and raw scales, partitioned the hybrid cells to equalize *cis* and* trans* precision, and tested each contrast against a null distribution generated under its null hypothesis.

## Notation

Cells are indexed by i. Cell groups (r) include the haploid *S. cerevisiae* parent (Pc), the haploid* S. eubayanus* parent (Pe), the* S. cerevisiae* allele in the hybrid (Hc), and the* S. eubayanus* allele in the hybrid (He). H is used for the hybrid as a whole. We omit the group index when a quantity is taken within a single group. We fit each gene separately and omit gene indexes (g) unless needed for clarity. We use lowercase letters for latent variables and realized observations. We use matching uppercase letters for the corresponding random variables, for fixed values calculated from the data, and for expectations built from parameters. Model parameters use Greek letters. A bar marks the mean of a variable over cells, and σ² with the variable as its subscript marks its variance. Δ marks contrasts between groups. Parameters of the summed allele counts of a hybrid cell carry the subscript Σ. A hat marks an estimate of a latent variable or parameter, and a star marks a quantity that still carries the transcriptome share distortion. A tilde marks an intermediate or comparison estimator.

| Symbol | Definition |
|---|---|
| y, Y | UMI count of one gene in one cell |
| L | Library size of a cell |
| D | Reference library size, the median library size |
| S | Depth factor of a cell, L/D |
| O | Offset of a cell, S for a parent and S/2 for each hybrid allele |
| N | Number of cells in a group |
| μ | Mean rate, the expression per genome, which is the expected count at the reference depth. |
| M | Expected count of a cell, O μ |
| φ | Dispersion, the variability beyond Poisson counting. |
| s | Slope in φ of the log probability of one cell |
| T | Boundary likelihood ratio statistic for φ |
| θ | Burst frequency, 1/φ |
| ζ | Burst size, μφ |
| κ | Kinetic balance, θ/ζ |
| Δ | Contrast between groups |
| f | Fraction of hybrid cells assigned to the *cis* estimate |
| σ² | Variance of the variable in its subscript. σ²_Δ is the variance of a contrast |
| γ | Ratio of the parental to the hybrid contrast variance |
| c | Transcriptome share, the probability that a UMI from a hybrid cell comes from the *S. cerevisiae* genome |
| t | Transcriptome share multiplier of an allele, its share relative to the mean share |
| η | Variance of a transcriptome share multiplier |
| a | Expression factor, a cell's departure from the group mean |
| d | Extrinsic deviation, shared by both alleles of a cell |
| q | Intrinsic deviation, private to one allele |
| ω | Intrinsic fraction of the dispersion |
| π | Count share of Hc in the summed allele counts |
| δ | Ploidy shift, the fall in dispersion from summing two alleles |
| λ | Weight of a group in a contrast |
| k | Common level of the hybrid alleles relative to the parents in the trans null |

## Sequencing depth offset

Sequencing depth and capture efficiency vary among cells and groups, which causes two problems. First, deeply sequenced cells yield more counts of every gene, and this cell-to-cell spread inflates the dispersion estimates, so genes appear noisier than they are. Second, groups differ in sequencing depth, which shifts their contrasts. We solved both problems with cell-specific offsets. The offset of each cell is its library size, L, divided by a reference depth, D. We set D to the median library size across all cells, calculated each offset once, and held it fixed.

$$
L_{ri} = \sum_g y_{gri}, \qquad L_{\mathrm{H}i} = L_{\mathrm{Hc}i} + L_{\mathrm{He}i} \\
D = \mathrm{median}\left(L_{ri}\right)_{}, \qquad S_{ri} = \frac{L_{ri}}{D}
$$

A parental cell took the offset of its own library. The two alleles of a hybrid cell share one library, so each allele received half of the hybrid offset, and summed hybrid counts used the full hybrid offset.

$$
O_{\mathrm{Pc}i} = S_{\mathrm{Pc}i}, \qquad O_{\mathrm{Pe}i} = S_{\mathrm{Pe}i}, \qquad O_{\mathrm{Hc}i} = O_{\mathrm{He}i} = \frac{S_{\mathrm{H}i}}{2}, \qquad O_{\mathrm{H}i} = S_{\mathrm{H}i}
$$

## Negative binomial model of single-cell counts

We modeled the counts of each gene in each group as negative binomial (NB). The NB distribution arises from a model of bursty transcription in which a gene switches between an inactive state and a short-lived active state that releases several transcripts at once. In the bursty limit of this two-state model, the steady-state mRNA number is NB, with the inverse of the dispersion, θ = 1/φ, equal to the burst frequency, and ζ = μφ the average burst size (Raj et al. 2006).

$$
Y_{ri} \sim \mathrm{NB}\left(O_i\,\mu_r,\ \varphi_r\right) \\

\mathrm{E}[Y_{ri}] = O_i\,\mu_r \equiv M_i, \qquad \mathrm{Var}[Y_{ri}] = M_i + \varphi_r\,M_i^2
$$

where μ_r is the mean rate of group r, the expected count of a cell sequenced at the reference depth D, and M_i is the expected count of cell i at its own depth. The dispersion sets how far the variance exceeds the Poisson variance. A small dispersion means little excess noise, and zero dispersion recovers the Poisson distribution. The squared coefficient of variation of a count separates the two sources of variance:

$$
\mathrm{CV}^2_{ri} = \frac{\mathrm{Var}[Y_{ri}]}{\mathrm{E}[Y_{ri}]^2} = \frac{O_i\,\mu_r + \varphi_r\,(O_i\,\mu_r)^2}{(O_i\,\mu_r)^2} = \frac{1}{O_i\,\mu_r} + \varphi_r
$$

The first term is the Poisson contribution. It falls as depth or expression rises because counting noise shrinks as more molecules are sampled. The second term does not depend on depth and is the excess variability that bursting adds.

**Estimating μ.** The NB likelihood of group r is the product over its cells of the NB probability of each count. Setting the derivative of its log with respect to μ to zero gives the maximum likelihood estimate μ̃_r (Lawless 1987). When the offsets differ among cells, μ̃_r has no closed form and is the sum of a term that depends on φ and a term that does not, so it is usually found iteratively and jointly with φ:

$$
\begin{aligned}
\Pr\left(y_{ri} \mid O_i\mu, \varphi\right) &= \frac{\Gamma\left(y_{ri} + \varphi^{-1}\right)}{\Gamma\left(\varphi^{-1}\right)\,y_{ri}!}\left(\frac{1}{1 + \varphi O_i\mu}\right)^{\varphi^{-1}}\left(\frac{\varphi O_i\mu}{1 + \varphi O_i\mu}\right)^{y_{ri}} \\
L_r(\mu, \varphi) &= \prod_i \Pr\left(y_{ri} \mid O_i\mu, \varphi\right) \\
\ell_r(\mu, \varphi) &= \sum_i \left[\log\Gamma\left(y_{ri} + \varphi^{-1}\right) - \log\Gamma\left(\varphi^{-1}\right) - \log\left(y_{ri}!\right) - \varphi^{-1}\log\left(1 + \varphi O_i\mu\right) + y_{ri}\log\frac{\varphi O_i\mu}{1 + \varphi O_i\mu}\right] \\
\frac{\partial \ell_r}{\partial \mu} &= \frac{1}{\mu}\sum_i \frac{y_{ri} - O_i\mu}{1 + \varphi O_i\mu}, \qquad \sum_i \frac{y_{ri} - O_i\tilde{\mu}_r}{1 + \varphi O_i\tilde{\mu}_r} = 0 \\
\tilde{\mu}_r &= \frac{\sum_i y_{ri}}{\sum_i O_i} - \frac{\varphi}{\sum_i O_i}\sum_i \frac{O_i\tilde{\mu}_r\left(y_{ri} - O_i\tilde{\mu}_r\right)}{1 + \varphi O_i\tilde{\mu}_r}
\end{aligned}
$$

At the true mean, the φ dependent term has expectation zero:

$$
\mathrm{E}\left[\sum_i \frac{M_i\left(Y_{ri} - M_i\right)}{1 + \varphi M_i}\right] = \sum_i \frac{M_i\left(\mathrm{E}[Y_{ri}] - M_i\right)}{1 + \varphi M_i} = 0
$$

The term therefore adds no systematic shift to μ̃_r. Instead, it changes the relative contribution of cells when estimating the mean. Dropping this term gives a closed form estimator, μ̂_r, that is independent of φ. 

$$
\hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i}
$$

This is the Poisson score with a log link and an offset. It assumes only that the mean count is O_i μ, which makes it the quasi-Poisson estimating equation (Wedderburn 1974, McCullagh and Nelder 1989). When φ = 0, every weight equals one and μ̃_r and μ̂_r are equal. Similarly, when every cell has the same offset, every weight is the same and cancels from the equation. Otherwise the two differ, and μ̂_r has three advantages. First, it keeps the mean separate from the dispersion. In a joint fit, μ̃_r depends on φ̂_r, so sampling error in φ̂_r reaches the mean and every mean contrast built from it. Orthogonality of μ and φ makes this dependence second order in large samples, but μ̂_r removes it at every sample size. Second, a joint fit gives equal weights for genes whose φ̂_r is zero and NB weights for the rest, so the mean estimator would differ among genes according to a test on their dispersion. μ̂_r treats every gene the same way. Third, μ̂_r is computed once from the counts and offsets, before any dispersion fit, so the dispersion fit runs with every cell mean fixed rather than alternating between μ and φ. This saving repeats for every gene in every observed, bootstrap, null and power-simulation fit.

The cost is precision. μ̂_r gives every count equal weight, whereas μ̃_r gives less weight to cells with large φM_i, whose counts vary the most. The variance of μ̂_r follows from the NB definition. The large-sample variance of μ̃_r is the inverse of the Fisher information, the variance of the score. The ratio of these variances gives:

$$
\begin{aligned}
\mathrm{Var}[\hat{\mu}_r] &= \frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i O_i\right)^2} = \mu_r^2\,\frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2} \\
\mathrm{Var}\left[\frac{\partial \ell_r}{\partial \mu}\right] &= \frac{1}{\mu_r^2}\sum_i \frac{\mathrm{Var}[Y_{ri}]}{\left(1 + \varphi\,M_i\right)^2} = \frac{1}{\mu_r^2}\sum_i \frac{M_i}{1 + \varphi\,M_i} \quad\Rightarrow\quad \mathrm{Var}[\tilde{\mu}_r] \approx \frac{\mu_r^2}{\sum_i M_i/\left(1 + \varphi\,M_i\right)} \\
\frac{\mathrm{Var}[\hat{\mu}_r]}{\mathrm{Var}[\tilde{\mu}_r]} &\approx \mu_r^2\,\frac{\sum_i M_i\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2} \cdot \frac{\sum_i M_i/\left(1 + \varphi\,M_i\right)}{\mu_r^2} = \frac{\sum_i M_i\left(1 + \varphi\,M_i\right)\,\sum_i M_i/\left(1 + \varphi\,M_i\right)}{\left(\sum_i M_i\right)^2}
\end{aligned}
$$

The Cauchy–Schwarz inequality gives (Σ_i M_i)² ≤ Σ_i M_i(1 + φM_i) · Σ_i M_i/(1 + φM_i), so the ratio of variances is at least one and μ̂_r is never more precise than μ̃_r. However, this ratio has an upper bound. Both the quotient M_i/(1 + φM_i) and the product M_i(1 + φM_i) increase with M_i, so across cells they are positively correlated. The product of their averages is then at most the average of their product:

$$
\begin{aligned}
\frac{\mathrm{Var}[\hat{\mu}_r]}{\mathrm{Var}[\tilde{\mu}_r]} &= \frac{\frac{1}{N_r}\sum_i M_i\left(1 + \varphi\,M_i\right) \cdot \frac{1}{N_r}\sum_i \frac{M_i}{1 + \varphi\,M_i}}{\left(\frac{1}{N_r}\sum_i M_i\right)^2} \le \frac{\frac{1}{N_r}\sum_i M_i\left(1 + \varphi\,M_i\right)\frac{M_i}{1 + \varphi\,M_i}}{\left(\frac{1}{N_r}\sum_i M_i\right)^2} \\
&= \frac{N_r\sum_i M_i^2}{\left(\sum_i M_i\right)^2} = \frac{N_r\sum_i O_i^2}{\left(\sum_i O_i\right)^2} = 1 + \mathrm{CV}^2(O)
\end{aligned}
$$

Here CV(O) is the coefficient of variation of the depth offsets across the N_r cells of the group. The cap depends only on the variability of the offsets, whereas the position of a gene between one and the cap depends on φM_i. When φM_i is small, as for lowly expressed genes or genes with φ near zero, the variances of the two estimators nearly coincide, and the ratio approaches the cap only when φM_i is large. 

**Estimating φ.** Because μ̂_r contains no φ, we fixed each cell's mean at its fitted value and estimated φ by maximizing the profile log-likelihood over φ ≥ 0. Under the NB model μ and φ are orthogonal, so fixing μ̂_r leaves only a bias of order 1/N_r, the cost of the degree of freedom spent on the mean. Relative to φ this bias is of order 1/(N_r M̄_r φ). We therefore used a Cox–Reid adjustment to remove the first-order part of this bias (Cox and Reid 1987, McCarthy et al. 2012). 

$$
\begin{aligned}
\ell_r(\varphi) &\equiv \ell_r\left(\hat{\mu}_r, \varphi\right), \qquad \hat{M}_i \equiv O_i\,\hat{\mu}_r \\
\ell^{\mathrm{CR}}_r(\varphi) &\equiv \ell_r(\varphi) - \tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}
\end{aligned}
$$

The maximum of ℓ^CR_r has no closed form, so we located it numerically over log φ. We write s_ri = ∂ log Pr(y_ri | M̂_i, φ)/∂φ for the slope of the log probability of cell i, so the slope of ℓ_r is Σ_i s_ri. Whether the maximum is at zero depends on the slope of ℓ^CR_r at zero. For ℓ_r, starting from the log probability of one cell, we write the ratio of gamma functions as a finite sum, collect the two terms in log(1 + φM), and expand to first order in φ:

$$
\begin{aligned}

\log\Pr(y \mid M, \varphi) &= \log\Gamma\left(y + \varphi^{-1}\right) - \log\Gamma\left(\varphi^{-1}\right) - \log y! - \varphi^{-1}\log\left(1 + \varphi M\right) + y\log\frac{\varphi M}{1 + \varphi M} \\
&= \sum_{j=0}^{y-1}\log\left(\varphi^{-1} + j\right) - \log y! - \varphi^{-1}\log\left(1 + \varphi M\right) + y\log\varphi + y\log M - y\log\left(1 + \varphi M\right) \\
&= -y\log\varphi + \sum_{j=0}^{y-1}\log\left(1 + j\varphi\right) - \log y! + y\log\varphi + y\log M - \left(\varphi^{-1} + y\right)\log\left(1 + \varphi M\right) \\
&= \sum_{j=0}^{y-1}\log\left(1 + j\varphi\right) + y\log M - \log y! - \left(\varphi^{-1} + y\right)\log\left(1 + \varphi M\right) \\
&= \tfrac{\varphi}{2}\,y(y-1) + y\log M - \log y! - \left(\varphi^{-1} + y\right)\left(\varphi M - \tfrac{\varphi^2 M^2}{2}\right) + \mathcal{O}(\varphi^2) \\
&= \tfrac{\varphi}{2}\,y(y-1) + y\log M - \log y! - M + \tfrac{\varphi}{2}\,M^2 - \varphi\,y\,M + \mathcal{O}(\varphi^2) \\
&= y\log M - M - \log y! + \tfrac{\varphi}{2}\left[y(y-1) + M^2 - 2yM\right] + \mathcal{O}(\varphi^2) \\
&= y\log M - M - \log y! + \tfrac{\varphi}{2}\left[(y - M)^2 - y\right] + \mathcal{O}(\varphi^2)
\end{aligned}
$$

Summing over cells then gives ℓ_r and its slope, written as the sum of the per-cell slopes s_ri:

$$
\begin{aligned}
\ell_r(\varphi) &= \sum_i \left[y_{ri}\log \hat{M}_i - \hat{M}_i - \log y_{ri}!\right] + \tfrac{\varphi}{2}\sum_i \left[(y_{ri} - \hat{M}_i)^2 - y_{ri}\right] + \mathcal{O}(\varphi^2) \\
\frac{\partial \ell_r}{\partial \varphi} &= \sum_i s_{ri}, \qquad s_{ri} = \tfrac12\left[(y_{ri} - \hat{M}_i)^2 - y_{ri}\right] + \mathcal{O}(\varphi) \\
\left.\frac{\partial \ell_r}{\partial \varphi}\right|_{\varphi = 0} &= \left.\sum_i s_{ri}\right|_{\varphi = 0} = \tfrac12\sum_i \left[(y_{ri} - \hat{M}_i)^2 - y_{ri}\right] = \tfrac12\left[\sum_i (y_{ri} - \hat{M}_i)^2 - \sum_i y_{ri}\right] = \tfrac12\left[\sum_i (y_{ri} - \hat{M}_i)^2 - \sum_i \hat{M}_i\right]
\end{aligned}
$$

The Cox-Reid adjustment adds its own slope at zero:

$$
\left.\frac{\partial}{\partial \varphi}\left(-\tfrac12\log\sum_i \frac{\hat{M}_i}{1 + \varphi\,\hat{M}_i}\right)\right|_{\varphi = 0} = \left.\tfrac12\,\frac{\sum_i \hat{M}_i^2/\left(1 + \varphi\,\hat{M}_i\right)^2}{\sum_i \hat{M}_i/\left(1 + \varphi\,\hat{M}_i\right)}\right|_{\varphi = 0} = \tfrac12\,\frac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i}
$$

An adjusted slope that is not positive makes φ = 0 a maximum at the boundary. We therefore define the estimator by these two cases:

$$
\hat{\varphi}_r = \begin{cases} 0 & \sum_i (y_{ri} - \hat{M}_i)^2 \le \sum_i \hat{M}_i - \dfrac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i} \\ \arg\max_{\varphi > 0}\ \ell^{\mathrm{CR}}_r(\varphi) & \text{otherwise} \end{cases}
$$

We used this estimator for the observed fits and for every bootstrap, null and power-simulation replicate. *Appendix A* gives a moment estimator for the dispersion that does not assume an NB distribution and uses only the first two moments of the counts. It is less precise, and we used it only for sensitivity analyses.

**Burst kinetics.** In the NB model, burst frequency and burst size are not estimated independently. Burst frequency is the inverse dispersion, θ = 1/φ, and burst size follows by division, ζ = μ/θ. Both estimates are linear in log₂ μ̂_r and log₂ φ̂_r:

$$
\log_2\hat{\theta}_r = -\log_2\hat{\varphi}_r, \qquad \log_2\hat{\zeta}_r = \log_2\hat{\mu}_r - \log_2\hat{\theta}_r = \log_2\hat{\mu}_r + \log_2\hat{\varphi}_r
$$

Because ζ̂_r is built from θ̂_r, any sampling error in θ̂_r enters ζ̂_r with the opposite sign. Expanding their covariance and using the orthogonality of the NB mean and dispersion gives:

$$
\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\zeta}_r\right] = \mathrm{Cov}\left[\log_2\hat{\theta}_r,\ \log_2\hat{\mu}_r - \log_2\hat{\theta}_r\right] = \mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] - \mathrm{Var}\left[\log_2\hat{\theta}_r\right] \approx -\mathrm{Var}\left[\log_2\hat{\theta}_r\right]
$$

This covariance is large when the dispersion is poorly estimated. Estimated burst frequency and burst size are therefore negatively correlated across genes even when the underlying kinetics are unrelated, and within a gene an imprecise dispersion estimate appears as burst frequency traded against burst size.

We therefore summarized burst kinetics with the kinetic balance κ = θ/ζ. The pair μ̂_r and κ̂_r carries the same information as θ̂_r and ζ̂_r:

$$
\log_2\hat{\kappa}_r = \log_2\hat{\theta}_r - \log_2\hat{\zeta}_r = 2\log_2\hat{\theta}_r - \log_2\hat{\mu}_r, \qquad \log_2\hat{\theta}_r = \tfrac12\left(\log_2\hat{\mu}_r + \log_2\hat{\kappa}_r\right), \qquad \log_2\hat{\zeta}_r = \tfrac12\left(\log_2\hat{\mu}_r - \log_2\hat{\kappa}_r\right)
$$

And kinetic balance is still coupled: 

$$
\mathrm{Cov}\left[\log_2\hat{\kappa}_r, \log_2\hat{\mu}_r\right] = 2\,\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] - \mathrm{Var}\left[\log_2\hat{\mu}_r\right] \approx -\mathrm{Var}\left[\log_2\hat{\mu}_r\right]
$$

Kinetic balance is still coupled to the mean, but only weakly, because the mean is estimated far more precisely than the dispersion. Converting each covariance to a correlation puts the two couplings on a common scale. With the mean and dispersion estimates uncorrelated, log₂ ζ̂_r has variance Var[log₂ μ̂_r] + Var[log₂ θ̂_r], and the two correlations are:

$$
\begin{aligned}
\rho\left(\log_2\hat{\theta}_r, \log_2\hat{\zeta}_r\right) &\approx \frac{-\mathrm{Var}[\log_2\hat{\theta}_r]}{\sqrt{\mathrm{Var}[\log_2\hat{\theta}_r]\left[\mathrm{Var}[\log_2\hat{\mu}_r] + \mathrm{Var}[\log_2\hat{\theta}_r]\right]}} = -\left(1 + \frac{\mathrm{Var}[\log_2\hat{\mu}_r]}{\mathrm{Var}[\log_2\hat{\theta}_r]}\right)^{-1/2} \\
\rho\left(\log_2\hat{\kappa}_r, \log_2\hat{\mu}_r\right) &\approx \frac{-\mathrm{Var}[\log_2\hat{\mu}_r]}{\sqrt{\mathrm{Var}[\log_2\hat{\kappa}_r]\,\mathrm{Var}[\log_2\hat{\mu}_r]}} = -\frac{\mathrm{SE}(\log_2\hat{\mu}_r)}{\mathrm{SE}(\log_2\hat{\kappa}_r)}
\end{aligned}
$$

When the dispersion is the less precise estimate, burst frequency and burst size are almost perfectly coupled, while kinetic balance couples to the mean only by the ratio of two standard errors, which is small. The cost is precision, because kinetic balance carries the dispersion estimate twice:

$$
\mathrm{Var}\left[\log_2\hat{\kappa}_r\right] = 4\,\mathrm{Var}\left[\log_2\hat{\theta}_r\right] + \mathrm{Var}\left[\log_2\hat{\mu}_r\right] - 4\,\mathrm{Cov}\left[\log_2\hat{\theta}_r, \log_2\hat{\mu}_r\right] \approx 4\,\mathrm{Var}\left[\log_2\hat{\theta}_r\right]
$$

Kinetic balance therefore carries about twice the standard error of burst frequency. It gives up precision to weaken its coupling to the mean.

UMI capture thins each count binomially, which leaves φ and θ unchanged but scales μ and ζ by the capture rate, so burst size is measured in captured molecules. A difference in capture between parental cells, or between parental and hybrid cells, enters the total, *trans* and dominance contrasts in the mean and burst size but not in burst frequency. For the hybrid alleles, the transcriptome share correction of the mean removes a genome-wide capture difference.

**Fano factor.** The Fano factor of a count, its variance divided by its mean, is a common summary of bursting. However, under the NB, the Fano factor is one plus the dispersion times the expected count, so it rises with the sequencing depth of each cell, and we did not use it. 

$$
\qquad \frac{\mathrm{Var}[Y_{ri}]}{\mathrm{E}[Y_{ri}]} = \frac{O_i\,\mu_r + \varphi_r\,(O_i\,\mu_r)^2}{O_i\,\mu_r} = 1 + \varphi_r\,O_i\,\mu_r
$$

## Transcriptome share corrections in the hybrid

A hybrid cell carries two genomes that share one library. Each captured UMI comes from one of the two genomes, and we call the probability that it comes from the *S. cerevisiae* genome the transcriptome share, c_i. It reflects the share of the cell's mRNA from each genome, weighted by how efficiently each genome's molecules are captured and assigned. It varies from cell to cell and is not observed. The observed fraction of UMIs, ĉ_i, estimates it with binomial sampling noise. The offsets do not capture the transcriptome share, so its variation affects the estimated means and dispersions of the hybrid alleles, and we corrected both.

**Transcriptome share estimate.** For each hybrid cell, we estimated its transcriptome share as the fraction of its UMIs assigned to the *S. cerevisiae* genome. The average of ĉ_i over the N_H hybrid cells, c̄̂, estimates the mean transcriptome share c̄.

$$
\hat{c}_i = \frac{L_{\mathrm{Hc}i}}{L_{\mathrm{Hc}i}+L_{\mathrm{He}i}} = \frac{L_{\mathrm{Hc}i}}{L_{\mathrm{H}i}}, \qquad \hat{\bar{c}} = \frac{1}{N_{\mathrm{H}}}\sum_i \hat{c}_i
$$

**Mean.** The mean estimator of the NB model is the ratio of summed counts to summed offsets, where each allele's offset is O_i = S_Hi/2:

$$
\hat{\mu}^*_{\mathrm{Hc}} = \frac{\sum_i y_{\mathrm{Hc}i}}{\sum_i S_{\mathrm{H}i}/2}, \qquad \hat{\mu}^*_{\mathrm{He}} = \frac{\sum_i y_{\mathrm{He}i}}{\sum_i S_{\mathrm{H}i}/2}
$$

This estimator is biased. Given its transcriptome share, the expected Hc count of cell i is its library size times its transcriptome share and its expression per genome. Since the offsets are fixed, the expectation of the ratio is the sum of the expected counts over the sum of the offsets. Averaging over the transcriptome share gives:

$$
\begin{aligned}
\mathrm{E}[Y_{\mathrm{Hc}i} \mid c_i] &= S_{\mathrm{H}i}\,c_i\,\mu_{\mathrm{Hc}}, \qquad \mathrm{E}[Y_{\mathrm{Hc}i}] = S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}} \\
\mathrm{E}\left[\hat{\mu}^*_{\mathrm{Hc}}\right] &= \frac{\sum_i S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}}{\sum_i S_{\mathrm{H}i}/2} = 2\bar{c}\,\mu_{\mathrm{Hc}}, \qquad \mathrm{E}\left[\hat{\mu}^*_{\mathrm{He}}\right] = 2(1-\bar{c})\,\mu_{\mathrm{He}}
\end{aligned}
$$

The fitted mean therefore equals μ_Hc only when c̄ is one half. We call 2c̄ and 2(1 − c̄) the abundance factors of the two alleles. Every group has 800 cells or more, so terms of 𝒪(1/N) are small, and we dropped them throughout. Dividing each fitted mean by its abundance factor gives estimators of the per-genome rates that are unbiased up to 𝒪(1/N_H):

$$
\hat{\mu}_{\mathrm{Hc}} = \frac{\hat{\mu}^*_{\mathrm{Hc}}}{2\hat{\bar{c}}}, \qquad \hat{\mu}_{\mathrm{He}} = \frac{\hat{\mu}^*_{\mathrm{He}}}{2(1-\hat{\bar{c}})}
$$

Because c̄ averages over all genes, this correction removes the genome-wide allele imbalance of the hybrid. A genome-wide difference in mRNA output between the two genomes and a difference in their capture or assignment efficiency are removed together. The *cis* mean contrasts therefore measure each gene relative to the hybrid's genome-wide allele balance, and a* cis* effect shared by every gene cannot be detected.

**Dispersion.** To find how the transcriptome share affects the dispersion estimate, we first need to know what an NB dispersion estimates. The NB arises when a Poisson count is scaled by an expression factor a that has mean one and a gamma distribution:

$$
Y_{ri} \mid a_{ri} \sim \mathrm{Poisson}\left(O_i\,\mu_r\,a_{ri}\right), \qquad a_{ri} \sim \mathrm{Gamma}\left(\text{shape} = 1/\varphi_r,\ \text{rate} = 1/\varphi_r\right) \\
\mathrm{E}[a_{ri}] = 1, \qquad \mathrm{Var}[a_{ri}] = \varphi_r
$$

Integrating the Poisson probability against the gamma density, with θ = 1/φ, gives the marginal distribution of the count:

$$
\begin{aligned}
\Pr(y) &= \int_0^\infty \underbrace{\frac{e^{-M_i a}(M_i a)^y}{y!}}_{\text{Poisson pmf}}\;\underbrace{\frac{\theta^\theta a^{\theta-1}e^{-\theta a}}{\Gamma(\theta)}}_{\text{gamma pdf}}\,da \\
&= \frac{M_i^y\,\theta^\theta}{y!\,\Gamma(\theta)}\int_0^\infty a^{y+\theta-1}\,e^{-(M_i+\theta)a}\,da \\
&= \frac{M_i^y\,\theta^\theta}{y!\,\Gamma(\theta)}\cdot\frac{\Gamma(y+\theta)}{(M_i+\theta)^{y+\theta}} \qquad \text{(using } \textstyle\int_0^\infty a^{s-1}e^{-ba}\,da = \Gamma(s)/b^s) \\
&= \frac{\Gamma(y+\theta)}{y!\,\Gamma(\theta)}\left(\frac{\theta}{M_i+\theta}\right)^{\theta}\left(\frac{M_i}{M_i+\theta}\right)^{y}
\end{aligned}
$$

This is the NB with mean M_i and dispersion φ. An NB fit therefore estimates the variance of the mean-one factor that scales the Poisson mean. For an allele, that factor is its expression factor, so the target is φ_Hc ≡ Var[a_Hci]. However, the fit absorbs every factor that scales the Poisson mean, including the transcriptome share in the hybrid. Writing the library size through the allele offset, S_Hi = 2O_i, and the transcriptome share relative to its mean, c_i = c̄ t_Hci, separates the expected Hc count into the offset, the fitted mean 2c̄μ_Hc, and two mean-one factors:

$$
\mathrm{E}[Y_{\mathrm{Hc}i} \mid c_i, a_{\mathrm{Hc}i}] = S_{\mathrm{H}i}\,c_i\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i} = 2O_i\left(\bar{c}\,\mu_{\mathrm{Hc}}\right)t_{\mathrm{Hc}i}\,a_{\mathrm{Hc}i}, \qquad t_{\mathrm{Hc}i} \equiv \frac{c_i}{\bar{c}}
$$

We assumed the transcriptome share multiplier t_Hci to be independent of the expression factor and of the depth factor. The fit then treats the product a*_Hci ≡ t_Hci a_Hci as its expression factor, marked with a star because it still carries the transcriptome share. The product has mean one, so φ̂*_Hc estimates its variance. We define η_Hc as the squared coefficient of variation of the transcriptome share, and have E[t²_Hc] = 1 + η_Hc:

$$
\begin{aligned}
\mathrm{E}\left[\hat{\varphi}^*_{\mathrm{Hc}}\right] \approx \mathrm{Var}[a^*_{\mathrm{Hc}i}] = \mathrm{Var}[t_{\mathrm{Hc}i}\,a_{\mathrm{Hc}i}] &= \mathrm{E}[t_{\mathrm{Hc}i}^2]\,\mathrm{E}[a_{\mathrm{Hc}i}^2] - 1 \qquad \eta_{\mathrm{Hc}} \equiv \mathrm{Var}[t_{\mathrm{Hc}i}] = \frac{\sigma^2_c}{\bar{c}^2}, \qquad \mathrm{E}[t_{\mathrm{Hc}i}^2] = \frac{\mathrm{E}[c_i^2]}{\bar{c}^2} = \frac{\sigma^2_c + \bar{c}^2}{\bar{c}^2} = 1 + \eta_{\mathrm{Hc}}\\
&= (1+\eta_{\mathrm{Hc}})(1+\varphi_{\mathrm{Hc}}) - 1 = \varphi_{\mathrm{Hc}} + (1+\varphi_{\mathrm{Hc}})\,\eta_{\mathrm{Hc}}
\end{aligned}
$$

The fitted dispersion therefore exceeds φ_Hc by (1 + φ_Hc)η_Hc, which has a fixed part η_Hc and a part proportional to the dispersion. Removing it requires η_Hc, and thus the variance of the true transcriptome share, σ²_c.

The sample variance of ĉ_i overstates σ²_c, because each ĉ_i also carries binomial sampling noise. The law of total variance splits Var[ĉ_i] into the variance of its conditional mean and the mean of its conditional variance. Given the library size and the true transcriptome share, the number of *S. cerevisiae* UMIs is binomial, so ĉ_i has conditional mean c_i and conditional variance c_i(1 − c_i)/L_Hi. Taking the library size to be independent of the true transcriptome share and substituting these moments gives the variance of the true transcriptome share plus the average sampling variance:

$$
\begin{gathered}
\mathrm{Var}[\hat{c}_i] = \mathrm{Var}\left[\mathrm{E}[\hat{c}_i \mid c_i, L_{\mathrm{H}i}]\right] + \mathrm{E}\left[\mathrm{Var}[\hat{c}_i \mid c_i, L_{\mathrm{H}i}]\right] \\
L_{\mathrm{Hc}i} \mid L_{\mathrm{H}i}, c_i \sim \mathrm{Binomial}\left(L_{\mathrm{H}i}, c_i\right) \\
\mathrm{E}[\hat{c}_i \mid c_i, L_{\mathrm{H}i}] = \frac{L_{\mathrm{H}i}\,c_i}{L_{\mathrm{H}i}} = c_i, \qquad \mathrm{Var}[\hat{c}_i \mid c_i, L_{\mathrm{H}i}] = \frac{L_{\mathrm{H}i}\,c_i(1-c_i)}{L_{\mathrm{H}i}^2} = \frac{c_i(1-c_i)}{L_{\mathrm{H}i}} \\
\mathrm{Var}[\hat{c}_i] = \mathrm{Var}[c_i] + \mathrm{E}\left[\frac{c_i(1-c_i)}{L_{\mathrm{H}i}}\right]
\end{gathered}
$$

This correction needs c_i, which is unknown, so we substituted ĉ_i, which adds a bias of its own. The library size is fixed under the conditioning, so the expectation of ĉ_i(1 − ĉ_i)/L_Hi is the mean of ĉ_i minus the mean of its square, divided by L_Hi, and the mean square is the binomial variance plus the squared mean. The plug-in term is therefore too small by the factor (L_Hi − 1)/L_Hi, and dividing by L_Hi − 1 instead of L_Hi removes it. Each cell then gives an unbiased estimate of its own sampling variance. Averaging these over cells and subtracting the average from the sample variance of ĉ_i gives the estimator of the variance of the transcriptome share.

$$
\begin{aligned}
\mathrm{E}\left[\frac{\hat{c}_i(1-\hat{c}_i)}{L_{\mathrm{H}i}} \,\middle|\, c_i, L_{\mathrm{H}i}\right] &= \frac{\mathrm{E}[\hat{c}_i \mid c_i, L_{\mathrm{H}i}] - \mathrm{E}[\hat{c}_i^2 \mid c_i, L_{\mathrm{H}i}]}{L_{\mathrm{H}i}} = \frac{1}{L_{\mathrm{H}i}}\left(c_i - \left[\frac{c_i(1-c_i)}{L_{\mathrm{H}i}} + c_i^2\right]\right) = \frac{c_i(1-c_i)}{L_{\mathrm{H}i}}\cdot\frac{L_{\mathrm{H}i}-1}{L_{\mathrm{H}i}} \\
\hat{\sigma}_c^2 &= \frac{1}{N_{\mathrm{H}} - 1}\sum_i \left(\hat{c}_i - \hat{\bar{c}}\right)^2 - \frac{1}{N_{\mathrm{H}}}\sum_i \frac{\hat{c}_i(1-\hat{c}_i)}{L_{\mathrm{H}i}-1}
\end{aligned}
$$

Solving the expectation of φ̂*_Hc for φ_Hc and replacing each quantity by its estimate gives the corrected allele dispersion:

$$
\hat{\varphi}_{\mathrm{Hc}} = \frac{\hat{\varphi}^*_{\mathrm{Hc}} - \hat{\eta}_{\mathrm{Hc}}}{1 + \hat{\eta}_{\mathrm{Hc}}}, \qquad \hat{\eta}_{\mathrm{Hc}} = \frac{\hat{\sigma}_c^2}{\hat{\bar{c}}^2}
$$

Both parts of η̂_Hc are unbiased, so η̂_Hc is unbiased up to a term of order 1/N_H. 

$$
\begin{aligned}
\mathrm{E}\left[\hat{\sigma}_c^2\right] &= \mathrm{Var}[\hat{c}_i] - \mathrm{E}\left[\frac{c_i(1-c_i)}{L_{\mathrm{H}i}}\right] = \sigma_c^2, \qquad \mathrm{E}\left[\hat{\bar{c}}\right] = \mathrm{E}[c_i] = \bar{c} \\
\hat{\eta}_{\mathrm{Hc}} &= \frac{\hat{\sigma}_c^2}{\hat{\bar{c}}^2} = \frac{\sigma_c^2}{\bar{c}^2} + \frac{\hat{\sigma}_c^2 - \sigma_c^2}{\bar{c}^2} - \frac{2\,\sigma_c^2\left(\hat{\bar{c}} - \bar{c}\right)}{\bar{c}^3} + \mathcal{O}_p\left(\frac{1}{N_{\mathrm{H}}}\right) \\
\mathrm{E}\left[\hat{\eta}_{\mathrm{Hc}}\right] &= \eta_{\mathrm{Hc}} + \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right)
\end{aligned}
$$

Analogous formulas apply for the He allele with 1 − c̄̂ in place of c̄̂. In our data the estimated mean transcriptome share c̄̂ is 0.5024, with η̂_Hc = 7.4 × 10⁻⁵ and η̂_He = 7.5 × 10⁻⁵. At the median dispersions of Hc (0.150) and He (0.152) the transcriptome share correction terms are 8.5 × 10⁻⁵ and 8.6 × 10⁻⁵, about 0.06% of each. Relative to the dispersion, the correction (1 + φ)η/φ grows as the dispersion falls, and it exceeds 1% of φ for dispersions below about 0.0075. We therefore corrected every hybrid allele dispersion and truncated negative values at zero.

## Per-genome hybrid estimates

The dominance contrasts compare the hybrid with its parents per genome, which requires the per-genome mean and dispersion of the hybrid.

**Mean.** We fitted the summed hybrid counts with the offset S_H, the full depth of the hybrid cell. When the two alleles contribute equally on average (c̄ = ½), the mean rate of the summed hybrid is the average of the two per-genome rates.

$$
\hat{\mu}^*_{\mathrm{H}} = \frac{\sum_i y_{\mathrm{H}i}}{\sum_i S_{\mathrm{H}i}} = \frac{\sum_i y_{\mathrm{Hc}i} + \sum_i y_{\mathrm{He}i}}{\sum_i S_{\mathrm{H}i}} = \tfrac12\left(\frac{\sum_i y_{\mathrm{Hc}i}}{\sum_i S_{\mathrm{H}i}/2} + \frac{\sum_i y_{\mathrm{He}i}}{\sum_i S_{\mathrm{H}i}/2}\right) = \tfrac12\left(\hat{\mu}^*_{\mathrm{Hc}} + \hat{\mu}^*_{\mathrm{He}}\right)
$$

When c̄ differs from one half, μ̂*_H converges instead to the transcriptome share-weighted average of the per-genome rates. We divided each allele fit by its abundance factor to remove this weighting:

$$
\hat{\mu}_{\mathrm{H}} = \tfrac12\left(\hat{\mu}_{\mathrm{Hc}} + \hat{\mu}_{\mathrm{He}}\right) = \tfrac12\left(\frac{\hat{\mu}^*_{\mathrm{Hc}}}{2\hat{\bar{c}}} + \frac{\hat{\mu}^*_{\mathrm{He}}}{2(1-\hat{\bar{c}})}\right)
$$

**Dispersion.** The per-genome dispersion of the hybrid is the dispersion of a single genome. We therefore averaged the two allele dispersions on the log₂ scale, which gives their geometric mean:

$$
\log_2\hat{\varphi}_{\mathrm{H}} = \tfrac12\left(\log_2\hat{\varphi}_{\mathrm{Hc}} + \log_2\hat{\varphi}_{\mathrm{He}}\right)
$$

Each allele carries the full extrinsic noise of its cell and its own intrinsic noise, so each corrected allele dispersion is a per-genome dispersion and their average is φ̂_H. Summing the two allele counts instead averages their intrinsic noise, which lowers the dispersion of the summed counts, φ_Σ, below φ_H by a ploidy shift δ. Appendix B derives δ, and we reported it for each gene as a measure of how much noise summing hides. 

Burst size and kinetic balance of the hybrid are built from the same per-genome quantities, μ̂_H and φ̂_H. The CV² is a raw-scale quantity, so it uses the arithmetic mean of the two allele dispersions:

$$
\log_2\hat{\zeta}_{\mathrm{H}} = \log_2\hat{\mu}_{\mathrm{H}} + \log_2\hat{\varphi}_{\mathrm{H}}, \qquad \log_2\hat{\kappa}_{\mathrm{H}} = -\log_2\hat{\mu}_{\mathrm{H}} - 2\log_2\hat{\varphi}_{\mathrm{H}}
$$

$$
\mathrm{CV}^2_{\mathrm{H}} = \frac{1}{\hat{\mu}_{\mathrm{H}}} + \tfrac12\left(\hat{\varphi}_{\mathrm{Hc}} + \hat{\varphi}_{\mathrm{He}}\right)
$$

On the raw scale, the burst size and kinetic balance of the hybrid are averages of the allele values rather than functions of μ̂_H and φ̂_H:

$$
\hat{\zeta}^{\mathrm{raw}}_{\mathrm{H}} = \tfrac12\left(\hat{\mu}_{\mathrm{Hc}}\,\hat{\varphi}_{\mathrm{Hc}} + \hat{\mu}_{\mathrm{He}}\,\hat{\varphi}_{\mathrm{He}}\right), \qquad \hat{\kappa}^{\mathrm{raw}}_{\mathrm{H}} = \tfrac12\left(\hat{\kappa}_{\mathrm{Hc}} + \hat{\kappa}_{\mathrm{He}}\right)
$$

The raw burst size and kinetic balance contrasts are the weighted sums of these per-group values, Σ_r λ_r ζ̂_r and Σ_r λ_r κ̂_r, with the group weights λ_r of the raw dispersion contrasts (*Null for the log₂ dispersion contrasts*).

Each bootstrap replicate refitted the allele means and dispersions, with their transcriptome share corrections, on its resampled hybrid cells, so the standard errors of every quantity built from them include the sampling error of the corrections.

## Contrasts

**Contrasts on two scales.** We compared a quantity x, such as the mean or the dispersion, between groups through contrasts on a log₂ scale for relative differences and a raw scale for absolute differences.

$$
\Delta(x) = \log_2 x_r - \log_2 x_{r'}, \qquad \Delta^{\mathrm{raw}}(x) = x_r - x_{r'}
$$

We used only the log₂ scale for mean contrasts and both scales for dispersion contrasts. On the log₂ scale, a doubling of φ counts the same whether the gene starts quiet or noisy, so log₂ dispersion contrasts compare genes on a common footing and can be compared directly with log₂ mean contrasts. When two groups differ by a constant factor in an underlying rate, their log₂ dispersions differ by a fixed amount regardless of the starting dispersion. The log₂ scale therefore suits questions about the mechanism behind a difference between groups.

On the raw scale, a contrast is the difference of two dispersions. Because the squared coefficient of variation of a count is 1/μ + φ, a raw difference in φ is the difference in cell-to-cell variability beyond Poisson sampling, so it measures how much noise a difference adds. The raw scale suits questions about the size of noise and its likely consequences for phenotype. It also handles genes at or near zero dispersion, which the log₂ contrast must exclude.

**Informative estimates.** Whether a dispersion estimate can support a contrast between groups depends on its precision. We therefore placed each estimate in one of three classes. First, an estimate is resolved above zero when a test against the Poisson model rejects at the 5% level. The test compares the adjusted log-likelihood at φ̂_r and at zero. Because zero lies on the boundary of the parameter space, the null distribution is an equal mixture of a point mass at zero and a χ² distribution with one degree of freedom:

$$
T_r = 2\left[\ell^{\mathrm{CR}}_r(\hat{\varphi}_r) - \ell^{\mathrm{CR}}_r(0)\right], \qquad T_r \sim \tfrac12\,\chi^2_0 + \tfrac12\,\chi^2_1 \ \text{under } \varphi_r = 0, \qquad \text{resolved above zero if } T_r > 2.71
$$

Second, an estimate is resolved near zero when it is not resolved above zero, but its upper 95% likelihood bound φ^up_r, scaled by the average fitted count per cell, falls below 0.1:

$$
\varphi^{\mathrm{up}}_r\,\bar{M}_r < 0.1, \qquad \bar{M}_r \equiv \frac{1}{N_r}\sum_i \hat{M}_i
$$

Since φM̄_r is the extra-Poisson variance as a fraction of the Poisson variance, these genes are effectively Poisson. We considered all other estimates uninformative.

The two scales need different precision from each dispersion estimate. A raw contrast can use estimates resolved above or near zero, while a log₂ contrast can use only estimates resolved above zero. A fit that returns zero can still be resolved near zero, so we recorded it as a valid estimate, separate from a failed fit. We formed a contrast only when every group it uses has an estimate in an allowed class:

| Class | Raw contrasts | Log₂ contrasts |
|---|---|---|
| Resolved above zero | yes | yes |
| Resolved near zero | yes | no |
| Uninformative | no | no |

**Contrasts between groups.** The total contrast compares the two parents, the *cis* contrast compares the two hybrid alleles, and the* trans* contrast is their difference. The midparent contrast compares the per-genome hybrid with the mean of the two parents, and the parental contrast compares it with one parent, P. For the mean:

$$
\begin{aligned}
\Delta_{\mathrm{total}}(\mu) &= \log_2\mu_{\mathrm{Pc}} - \log_2\mu_{\mathrm{Pe}}, \qquad \Delta_{\mathrm{cis}}(\mu) = \log_2\mu_{\mathrm{Hc}} - \log_2\mu_{\mathrm{He}}, \qquad \Delta_{\mathrm{trans}}(\mu) = \Delta_{\mathrm{total}}(\mu) - \Delta_{\mathrm{cis}}(\mu) \\
\Delta_{\mathrm{mid}}(\mu) &= \log_2\mu_{\mathrm{H}} - \log_2\mu_{MP}, \qquad \Delta_{\mathrm{par}}(\mu) = \log_2\mu_{\mathrm{H}} - \log_2\mu_{\mathrm{P}}, \qquad \mu_{MP} = \tfrac12\left(\mu_{\mathrm{Pc}} + \mu_{\mathrm{Pe}}\right)
\end{aligned}
$$

The dispersion contrasts take the same form, with geometric means for the per-genome hybrid and the midparent on the log₂ scale and arithmetic means on the raw scale:

$$
\begin{aligned}
\Delta_{\mathrm{total}}(\varphi) &= \log_2\varphi_{\mathrm{Pc}} - \log_2\varphi_{\mathrm{Pe}}, \qquad \Delta_{\mathrm{cis}}(\varphi) = \log_2\varphi_{\mathrm{Hc}} - \log_2\varphi_{\mathrm{He}}, \qquad \Delta_{\mathrm{trans}}(\varphi) = \Delta_{\mathrm{total}}(\varphi) - \Delta_{\mathrm{cis}}(\varphi) \\
\Delta_{\mathrm{mid}}(\varphi) &= \tfrac12\left(\log_2\varphi_{\mathrm{Hc}} + \log_2\varphi_{\mathrm{He}}\right) - \tfrac12\left(\log_2\varphi_{\mathrm{Pc}} + \log_2\varphi_{\mathrm{Pe}}\right), \qquad \Delta_{\mathrm{par}}(\varphi) = \tfrac12\left(\log_2\varphi_{\mathrm{Hc}} + \log_2\varphi_{\mathrm{He}}\right) - \log_2\varphi_{\mathrm{P}} \\
\Delta^{\mathrm{raw}}_{\mathrm{total}}(\varphi) &= \varphi_{\mathrm{Pc}} - \varphi_{\mathrm{Pe}}, \qquad \Delta^{\mathrm{raw}}_{\mathrm{cis}}(\varphi) = \varphi_{\mathrm{Hc}} - \varphi_{\mathrm{He}}, \qquad \Delta^{\mathrm{raw}}_{\mathrm{trans}}(\varphi) = \Delta^{\mathrm{raw}}_{\mathrm{total}}(\varphi) - \Delta^{\mathrm{raw}}_{\mathrm{cis}}(\varphi) \\
\Delta^{\mathrm{raw}}_{\mathrm{mid}}(\varphi) &= \tfrac12\left(\varphi_{\mathrm{Hc}} + \varphi_{\mathrm{He}}\right) - \tfrac12\left(\varphi_{\mathrm{Pc}} + \varphi_{\mathrm{Pe}}\right), \qquad \Delta^{\mathrm{raw}}_{\mathrm{par}}(\varphi) = \tfrac12\left(\varphi_{\mathrm{Hc}} + \varphi_{\mathrm{He}}\right) - \varphi_{\mathrm{P}}
\end{aligned}
$$

For each type of contrast, the burst frequency, burst size and kinetic balance contrasts follow from its mean and dispersion contrasts:

$$
\Delta(\theta) = -\Delta(\varphi), \qquad \Delta(\zeta) = \Delta(\mu) - \Delta(\theta), \qquad \Delta(\kappa) = \Delta(\theta) - \Delta(\zeta) = 2\,\Delta(\theta) - \Delta(\mu)
$$

Each estimated contrast replaces every quantity by its corrected estimate, so it converges to its true contrast within 𝒪(1/N) because each estimator does.

**Kinetic balance.** For a contrast between independent groups, the variances of the two groups add and the covariances between groups are zero, so the covariances, correlations and precision cost derived in *Burst kinetics* carry over, with each variance replaced by its sum over the two groups. For the* c**is* contrasts, whose two alleles share cells, we assumed that the covariance between the mean estimate of one allele and the dispersion estimate of the other is negligible, so the same results hold. We formed a kinetic balance contrast only when every group it uses has a dispersion estimate resolved above zero.

Together, the kinetic balance and mean contrasts separate two biological questions. Δ(μ) asks how much expression differs. Δ(κ) then asks how that difference divides between differences in burst frequency and burst size. Rewriting the relationship between μ, θ, and ζ, we have:

$$
\Delta(\mu) = \Delta(\theta) + \Delta(\zeta)
$$

Setting either Δ(θ) or Δ(ζ) to zero, or setting the two equal, and substituting produces:

$$
\begin{aligned}
\text{frequency only, } \Delta(\zeta) = 0 &: \quad \Delta(\theta) = \Delta(\mu) \ \Rightarrow\ \Delta(\kappa) = 2\,\Delta(\mu) - \Delta(\mu) = \Delta(\mu) \\
\text{size only, } \Delta(\theta) = 0 &: \quad \Delta(\kappa) = 0 - \Delta(\mu) = -\Delta(\mu) \\
\text{equal, } \Delta(\theta) = \Delta(\zeta) &: \quad \Delta(\theta) = \tfrac{1}{2}\Delta(\mu) \ \Rightarrow\ \Delta(\kappa) = \Delta(\mu) - \Delta(\mu) = 0
\end{aligned}
$$

A difference carried entirely by burst frequency therefore gives Δ(κ) = Δ(μ), one carried entirely by burst size gives Δ(κ) = −Δ(μ), and equal log₂ differences in both give Δ(κ) = 0.

## Statistical power of *cis* and* trans* estimates

In the hybrid, both alleles share one nucleus and one set of *trans*-acting factors, so their difference measures* cis* divergence alone (Wittkopp et al. 2004).* Trans* divergence is the remainder of the parental difference. The two estimates differ in statistical power because they draw on different data. Fluctuations shared by both alleles of a hybrid cell, such as cell size, sequencing depth and extrinsic fluctuations like cell-cycle state, largely cancel in the* cis* contrast. The* trans* contrast subtracts the hybrid ratio from the parental ratio, and the parental ratio comes from two independent sets of cells, so its variance adds to the* trans* contrast without cancellation. If all hybrid cells were used for both components,* trans* would always be the less precise estimate, with lower power than* cis* to detect differences, and its estimation error would correlate negatively with that of* cis* (Fraser 2019). We resolved both problems by using a fraction of the hybrid cells to estimate* cis* and the rest to estimate* trans*:

$$
N_{\mathrm{cis}} = f_x\,N_{\mathrm{H}}, \qquad N_{\mathrm{trans}} = (1-f_x)\,N_{\mathrm{H}}
$$

where N_H is the number of hybrid cells and f_x the fraction assigned to *cis*.

Since the two subsets share no hybrid cells, the *cis* and* trans* estimates are independent and estimation errors in one do not affect the other. We also chose the partition so that the two contrasts have equal sampling variance, and therefore equal power to detect an effect of a given size.

**Variance of the cis contrast.** Both hybrid alleles are fit on the same *cis* cells, so the variance of the difference of their log₂ estimates includes their covariance. In large samples each variance and the covariance are inversely proportional to the number of cells used, here f_x N_H. The whole sum is therefore a constant divided by f_x N_H. We call this constant σ²_Δ,H, the per-cell variance of the hybrid allele ratio.

$$
\begin{aligned}
\mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Hc}} - \log_2 \hat{x}_{\mathrm{He}}\right] &= \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Hc}}\right] + \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{He}}\right] - 2\,\mathrm{Cov}\left[\log_2 \hat{x}_{\mathrm{Hc}},\ \log_2 \hat{x}_{\mathrm{He}}\right] \\
&= \frac{\sigma^2_{\Delta,\mathrm{H}}}{f_x\,N_{\mathrm{H}}}
\end{aligned}
$$

Since the alleles share depth and extrinsic fluctuations, the covariance is positive. It therefore lowers the variance of the *cis* contrast below that of two alleles measured in separate cells.

**Variance of the trans contrast.** The *trans* contrast subtracts the hybrid allele ratio, measured on the remaining hybrid cells, from the parental ratio. The two parents are measured in different cells, so their estimates have no covariance, and the variance of the parental ratio is the sum of two terms.

$$
\sigma^2_{\Delta,P} = \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Pc}}\right] + \mathrm{Var}\left[\log_2 \hat{x}_{\mathrm{Pe}}\right]
$$

The parental cells and the *trans* hybrid cells are also disjoint, so the two parts of the* trans* contrast are independent and their variances add with no covariance. The fitted difference for the* trans* contrast is the parental log ratio minus the hybrid allele log ratio, with the two hybrid alleles fit only on the cells of the* trans* subset. Substituting these gives:

$$
\mathrm{Var}\left\{\left[\log_2 \hat{x}_{\mathrm{Pc}} - \log_2 \hat{x}_{\mathrm{Pe}}\right] - \left[\log_2 \hat{x}_{\mathrm{Hc}} - \log_2 \hat{x}_{\mathrm{He}}\right]\right\} = \sigma^2_{\Delta,P} + \frac{\sigma^2_{\Delta,\mathrm{H}}}{(1-f_x)\,N_{\mathrm{H}}}
$$

**Solving for the balanced statistical power.** To balance statistical power between *cis* and *trans* estimates, we set the *cis* variance equal to the* trans* variance and solved for f_x. Multiplying both sides by the number of hybrid cells and dividing by the per-cell variance gives:

$$
\begin{aligned}
\frac{\sigma^2_{\Delta,\mathrm{H}}}{f_x\,N_{\mathrm{H}}} &= \sigma^2_{\Delta,P} + \frac{\sigma^2_{\Delta,\mathrm{H}}}{(1-f_x)\,N_{\mathrm{H}}} \\
\frac{1}{f_x} &= \frac{N_{\mathrm{H}}\,\sigma^2_{\Delta,P}}{\sigma^2_{\Delta,\mathrm{H}}} + \frac{1}{1-f_x}
\end{aligned}
$$

The partition now depends on the data through a single ratio, γ. It divides the parental variance by the variance a hybrid allele ratio would have if all hybrid cells were used, so it measures the fixed parental variance in units of the best precision the hybrid can achieve.

$$
\gamma = \frac{N_{\mathrm{H}}\,\sigma^2_{\Delta,P}}{\sigma^2_{\Delta,\mathrm{H}}} \quad\Longrightarrow\quad \frac{1}{f_x} = \gamma + \frac{1}{1-f_x}
$$

Multiplying through by both denominators and collecting terms gives a quadratic in the *cis* fraction with two roots.

$$
\begin{aligned}
1 - f_x &= \gamma\,f_x\,(1-f_x) + f_x \quad\Longrightarrow\quad \gamma\,f_x^2 - (2 + \gamma)\,f_x + 1 = 0 \\
\Longrightarrow\quad f_x &= \frac{(2 + \gamma) \pm \sqrt{(2 + \gamma)^2 - 4\gamma}}{2\gamma} = \frac{(2 + \gamma) \pm \sqrt{\gamma^2 + 4}}{2\gamma}
\end{aligned}
$$

The root with the minus sign is the balanced partition, because the root with the plus sign exceeds one.

Two limits show how the balanced partition behaves. When the parental variance is negligible compared to the hybrid variance (γ → 0), the optimum returns to an even partition. When the parental variance is much larger than the hybrid variance, the partition falls as the inverse of the ratio. Power is therefore equalized by moving hybrid cells from estimating *cis* to estimating* trans* as the fixed parental variance grows.

$$
\begin{aligned}
\gamma \to 0&: \quad \sqrt{\gamma^2 + 4} \approx 2 + \frac{\gamma^2}{4} \;\Longrightarrow\; f_x \approx \frac12 - \frac{\gamma}{8} \\
\gamma \to \infty&: \quad f_x \approx \frac{1}{\gamma}
\end{aligned}
$$

**Separate partitions for each quantity.** The optimal partition depends on the parental variance relative to the hybrid variance, which differs among quantities. Estimates of φ carry larger and more skewed sampling error than estimates of μ, and kinetic balance carries about four times the dispersion variance, so we estimated γ separately for each. Mean contrasts used the partition at f_μ, dispersion contrasts the partition at f_φ, and kinetic balance contrasts the partition at f_κ, in which the mean and the dispersion are fitted on the same cells. Burst size contrasts combined Δ(μ) from the first partition with Δ(φ) from the second, and their standard errors treated the two as uncorrelated.

Sharing cells between the mean and dispersion partitions, or fitting both on the same cells as in the kinetic balance partition, does not correlate the mean and dispersion contrasts to first order. Within one group, let the mean be fitted on a set of cells A and the dispersion on a set B, which may share cells. The error of the fitted mean is a sum of count deviations over A. The fitted dispersion sets the slope of the log-likelihood in φ to zero (*Estimating φ*), so to first order its error is a constant c times the sum over B of the per-cell slopes s_i, evaluated at the true mean and dispersion. Writing out the covariance of the two errors with these definitions, the expected counts M_i are constants, a count and a slope from different cells are independent, and only cells in both sets remain:

$$
\begin{aligned}
\hat{\mu}_A - \mu &= \frac{\sum_{i \in A}\left(Y_i - M_i\right)}{\sum_{i \in A} O_i}, \qquad \hat{\varphi}_B - \varphi \approx c\sum_{j \in B} s_j \\
\mathrm{Cov}\left[\hat{\mu}_A - \mu,\ \hat{\varphi}_B - \varphi\right] &= \mathrm{Cov}\left[\frac{\sum_{i \in A}\left(Y_i - M_i\right)}{\sum_{i \in A} O_i},\ c\sum_{j \in B} s_j\right] = \frac{c}{\sum_{i \in A} O_i}\sum_{i \in A}\sum_{j \in B}\mathrm{Cov}\left[Y_i, s_j\right] \\
&= \frac{c}{\sum_{i \in A} O_i}\sum_{i \in A \cap B}\mathrm{Cov}\left[Y_i, s_i\right], \qquad \text{since } \mathrm{Cov}\left[Y_i, s_j\right] = 0 \text{ for } i \ne j
\end{aligned}
$$

Each remaining term is also zero. A slope has mean zero, so the covariance is the mean of the count times its slope, and the expected count does not depend on φ:

$$
\mathrm{Cov}\left[Y_i, s_i\right] = \mathrm{E}\left[Y_i\,s_i\right] = \sum_y y\,\frac{\partial \Pr(y \mid M_i, \varphi)}{\partial \varphi} = \frac{\partial}{\partial \varphi}\sum_y y\,\Pr(y \mid M_i, \varphi) = \frac{\partial M_i}{\partial \varphi} = 0
$$

The mean and dispersion errors of one group are therefore uncorrelated whether or not two partitions share its cells.

**Estimating partition fractions.** To estimate f_x for each quantity, we first estimated σ²_Δ,H and σ²_Δ,P for that quantity. We bootstrapped cells from the undivided hybrid and from both parents, refitting the hybrid alleles on the same resampled cells in each replicate to preserve their covariance. Because per-gene ratios of squared bootstrap standard errors are noisy and right-skewed, we pooled the variance components across genes by their medians. The pooled γ is the ratio of their medians across genes. For kinetic balance the medians ran over genes whose dispersion estimates are resolved above zero on both sides of the contrast, because κ is undefined otherwise. We reported γ̂_x and f_x for each quantity, with the resulting numbers of cis and trans cells, so that the partition sizes can be checked against the large-sample assumption.

$$
\begin{aligned}
\hat{\sigma}^2_{\Delta,\mathrm{H},g} &= N_{\mathrm{H}}\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Hc},g} - \log_2 \hat{x}_{\mathrm{He},g}\right], \qquad \hat{\sigma}^2_{\Delta,P,g} = \widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Pc},g} - \log_2 \hat{x}_{\mathrm{Pe},g}\right] \\
\hat{\gamma} &= \frac{N_{\mathrm{H}}\,\mathrm{median}_g\,\hat{\sigma}^2_{\Delta,P,g}}{\mathrm{median}_g\,\hat{\sigma}^2_{\Delta,\mathrm{H},g}} = \frac{\mathrm{median}_g\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Pc},g} - \log_2 \hat{x}_{\mathrm{Pe},g}\right]}{\mathrm{median}_g\,\widehat{\mathrm{Var}}\left[\log_2 \hat{x}_{\mathrm{Hc},g} - \log_2 \hat{x}_{\mathrm{He},g}\right]}
\end{aligned}
$$

These estimates use only variances, not *cis* or* trans* point estimates, so they do not bias the* cis* and* trans* tests (Wittes and Brittain 1990). Within each partition, the hybrid cells used for* cis* and* trans* span the full range of sequencing depth.

## Null distributions

To identify significant contrasts, we tested each contrast against a null distribution generated under its null hypothesis. Every contrast is a difference of per-group or per-allele estimates after the transcriptome share corrections, so we built one null for each of the three statistics, the log₂ mean contrast, the log₂ dispersion contrast, and the raw dispersion contrast, and applied it to every type of contrast. Each null draw was refitted with the same estimators and corrections as the observed data. The permuted and simulated draws keep the observed cells and their libraries, so the transcriptome share estimates c̄̂ and η̂_r are the observed values in every draw, and k̂ of the trans null is held fixed. When the sensitivity analysis used the weighted moment estimator (Appendix A), its null draws were refits of that estimator on the same permuted or simulated counts. With c the number of the B draws whose contrast is at least as extreme as the observed one, the p-value is (1 + c)/(1 + B), which is never zero (Phipson and Smyth 2010).

### Null for the log₂ mean contrasts

Under the null hypothesis the groups compared share one per-genome mean rate. The offsets put every cell on a common depth scale. For the *cis* and parental nulls, dividing each hybrid allele's counts by its abundance factor, 2c̄̂ for Hc and 2(1 − c̄̂) for He, puts them on the per-genome scale of a parent. The* trans* null places the same factor in the allele offsets instead. The cells of the groups compared then share one expected count per unit offset, and we permuted them. Exchangeability also requires the pooled cells to share a dispersion, so where the groups differ in dispersion the permutation null is approximate. We simulated the midparent mean null instead, because a pooled null would be miscalibrated.

**Cis.** The two alleles are measured in the same cells and share their extrinsic noise. We swapped the two rescaled allele counts within each cell at random, which keeps that shared noise and leaves the offsets unchanged.

**Trans.** Under the null the two ratios are equal, μ_Pc/μ_Pe = μ_Hc/μ_He. This allows the per-genome level of the hybrid to differ from the parents' by a common factor k, with μ_Hc = kμ_Pc and μ_He = kμ_Pe. We estimated k from the corrected allele means of the *trans* subset and the parental means, and gave each hybrid allele an offset that includes both its abundance factor and k, leaving its counts unchanged:

$$
\hat{k} = \sqrt{\frac{\hat{\mu}_{\mathrm{Hc}}\,\hat{\mu}_{\mathrm{He}}}{\hat{\mu}_{\mathrm{Pc}}\,\hat{\mu}_{\mathrm{Pe}}}}, \qquad O_{\mathrm{Hc}i} = \hat{k}\,\hat{\bar{c}}\,S_{\mathrm{H}i}, \qquad O_{\mathrm{He}i} = \hat{k}\,(1-\hat{\bar{c}})\,S_{\mathrm{H}i}
$$

The expected Hc count per unit offset is then S_Hi c̄ μ_Hc/(k c̄ S_Hi) = μ_Pc under the null, and likewise for He. We pooled the Pc cells with the Hc cells and split the pool at random, and did the same for Pe and He with an independent split. Each pool then shares one expected count per unit offset regardless of the dominance of the gene. The *trans* contrast does not change when both hybrid alleles are scaled by a common factor, so k̂ changes only the null distribution, not the observed contrast. We held k̂ fixed across permutations, since its sampling error is of order 1/N.

**Parental.** We pooled the hybrid cells with the cells of one parent. For each hybrid cell we summed its two rescaled allele counts. From the expected allele counts (*Transcriptome share corrections in the hybrid*), this sum has expectation S_Hi μ_H, so with offset S_Hi it shares the per-genome mean of a parental cell under the null:

$$
\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}}{2\hat{\bar{c}}} + \frac{Y_{\mathrm{He}i}}{2(1-\hat{\bar{c}})}\right] = \frac{S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}}{2\bar{c}} + \frac{S_{\mathrm{H}i}\,(1-\bar{c})\,\mu_{\mathrm{He}}}{2(1-\bar{c})} + \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right) = S_{\mathrm{H}i}\,\tfrac12\left(\mu_{\mathrm{Hc}} + \mu_{\mathrm{He}}\right) + \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right) = S_{\mathrm{H}i}\,\mu_{\mathrm{H}} + \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right)
$$

**Midparent.** A synthetic hybrid built from independent parental cells lacks the extrinsic noise that the two alleles of a hybrid cell share, so its noise differs from the hybrid's and a pooled permutation null would be miscalibrated. We therefore simulated the midparent mean null from the fitted model, as for the dispersion nulls. Under the null the per-genome hybrid mean equals the midparent, μ_H = μ_MP. We set the null allele means by scaling both corrected allele means by one common factor, which keeps their ratio and so the *cis* contrast:

$$
\mu^0_{\mathrm{Hc}} = \hat{\mu}_{\mathrm{Hc}}\,\frac{\hat{\mu}_{MP}}{\hat{\mu}_{\mathrm{H}}}, \qquad \mu^0_{\mathrm{He}} = \hat{\mu}_{\mathrm{He}}\,\frac{\hat{\mu}_{MP}}{\hat{\mu}_{\mathrm{H}}}, \qquad \tfrac12\left(\mu^0_{\mathrm{Hc}} + \mu^0_{\mathrm{He}}\right) = \hat{\mu}_{MP}
$$

Each draw simulated the two allele counts of every hybrid cell from these means with their abundance factors and offsets, with the dispersion each allele fit would show, φ̂_r(1 + η̂_r) + η̂_r, and with the shared extrinsic factor of the dispersion nulls. It simulated the parental counts from their fitted means and dispersions, refitted every group with the same estimators and corrections, and recomputed the midparent contrast.

### Null for the log₂ dispersion contrasts

A dispersion null cannot be built by pooling cells, because the groups compared may differ in their means under the null. Pooling cells with different means and refitting gives each pseudo-group a mixture of means, which inflates its fitted dispersion and widens the null. We therefore simulated each null draw from the fitted model with the contrast set to zero. Each group kept its own fitted expected counts M̂_i, and its dispersion was replaced by a null value. Writing a contrast as a weighted sum of log₂ dispersions, with weights λ_r of ±1 or ±½ for the groups r it uses, the null values are the fitted values moved the shortest distance that makes the contrast zero:

$$
\Delta(\hat{\varphi}) = \sum_r \lambda_r\log_2\hat{\varphi}_r, \qquad \log_2\varphi^0_r = \log_2\hat{\varphi}_r - \frac{\lambda_r\,\Delta(\hat{\varphi})}{\sum_{r'} \lambda_{r'}^2}, \qquad \sum_r \lambda_r\log_2\varphi^0_r = 0
$$

Parental counts were simulated as NB with mean M̂_i and dispersion φ⁰_r. A hybrid allele was simulated with the dispersion its uncorrected fit would show, φ⁰_r(1 + η̂_r) + η̂_r, so that the refit and its transcriptome share correction return φ⁰_r. The two alleles of a cell share extrinsic noise, which lowers the variance of the *cis* contrast. We kept it by drawing each allele's expression factor as the product of a factor shared by the cell, with variance σ̂²_d (Appendix B), and a private factor whose variance gives the allele its null dispersion. All factors have mean one and are gamma distributed:

$$
a_{ri} = d_i\,q_{ri}, \qquad \mathrm{Var}[d_i] = \hat{\sigma}^2_d, \qquad \mathrm{Var}[q_{ri}] = \frac{1 + \varphi^0_r(1 + \hat{\eta}_r) + \hat{\eta}_r}{1 + \hat{\sigma}^2_d} - 1, \qquad \mathrm{Cov}[a_{\mathrm{Hc}i}, a_{\mathrm{He}i}] = \hat{\sigma}^2_d
$$

We set a private variance below zero to zero. The burst size and kinetic balance contrasts, including the supplementary kinetic balance contrasts built from the f_μ and f_φ partitions, used the same construction with the fitted means held fixed and only the dispersions moved, until the contrast tested is zero:

$$
\log_2\varphi^0_r = \log_2\hat{\varphi}_r - \frac{\lambda_r\,\Delta(\hat{\zeta})}{\sum_{r'} \lambda_{r'}^2} \ \ \text{(burst size)}, \qquad \log_2\varphi^0_r = \log_2\hat{\varphi}_r + \frac{\lambda_r\,\Delta(\hat{\kappa})}{2\sum_{r'} \lambda_{r'}^2} \ \ \text{(kinetic balance)}
$$

Burst frequency needs no separate null, because Δ(θ) = −Δ(φ) and its test is the log₂ dispersion test. The CV² and the raw burst size and kinetic balance contrasts have no null, so we report their estimates and bootstrap intervals without p-values.

### Null for the raw dispersion contrasts

The raw null used the same simulation with the null values moved on the raw scale, φ⁰_r = φ̂_r − λ_r Δ^raw(φ̂)/Σ_r′ λ_r′², and null values below zero set to zero. When truncation moves a null value, the contrast is no longer exactly zero, so we repeated the projection over the groups that remain positive, holding the truncated groups at zero. For the two-group contrasts, total and cis, the raw and log₂ null hypotheses are the same, equal dispersions, so the raw tests reused the log₂ draws. Their common null value is the geometric mean of the two fitted dispersions rather than the arithmetic mean that the raw projection would give. For the *trans*, midparent and parental contrasts the two null hypotheses differ, because a ratio and a difference of dispersions are zero under different conditions, so the two nulls were generated separately.

## Model assumptions

| Assumption | Statement | Used in |
|---|---|---|
| Independent cells | Counts of different cells are independent given their offsets | Negative binomial model, Statistical power, Null distributions |
| Bursty limit | The NB arises from bursty transcription, so θ = 1/φ is the burst frequency and ζ = μφ the burst size | Negative binomial model, Burst kinetics |
| Fixed offsets | Offsets are library sizes divided by D and are held fixed, although a gene's own counts enter its library size | Sequencing depth offset, Transcriptome share corrections |
| Fitted dispersion | The NB fit converges to the variance of the expression factor only when that factor is gamma | Negative binomial model, Transcriptome share corrections |
| Large samples | Every group has about 800 cells or more, so terms of order 1/N are small and are dropped, except the plug-in bias of φ̂, which the Cox–Reid adjustment removes to first order. Partition sizes are reported in Statistical power | Estimating φ, Transcriptome share corrections, Contrasts |
| Capture efficiency | UMI capture thins counts binomially, which leaves φ and θ unchanged but scales μ and ζ. A capture difference between parental cells, or between parental and hybrid cells, enters the total, trans and dominance contrasts in μ and ζ. A genome-wide capture difference between the hybrid alleles is removed by the transcriptome share correction | Burst kinetics, Contrasts |
| Hybrid allele counts | Allele counts follow the parental NB form scaled by the transcriptome share multiplier c_i/c̄ for Hc and (1 − c_i)/(1 − c̄) for He | Transcriptome share corrections, Appendix B |
| Independent transcriptome share | c_i is independent of the offset and of every expression factor | Transcriptome share corrections, Appendix B |
| Binomial transcriptome share sampling | Genome-assigned UMIs are binomial given the hybrid library size and c_i, and library size is independent of c_i | Transcriptome share corrections |
| Transcriptome share correction | The transcriptome share corrections of the hybrid dispersions and variance estimates are estimated from the transcriptome share variance and subtracted | Transcriptome share corrections, Null distributions |
| Zero estimates | A zero dispersion estimate is a valid boundary value and is kept. Raw contrasts can use it when it is resolved near zero, and log₂ contrasts need estimates resolved above zero | Contrasts |
| Ploidy effects | Ploidy effects act equally on both alleles, so they cancel in the cis ratio and within each ratio of the trans contrast. A species-specific ploidy effect would remain in the trans divergence, and separating it requires diploid parental strains | Contrasts |
| Cross-allele covariance | In the cis kinetic contrasts, the covariance between the mean estimate of one allele and the dispersion estimate of the other is negligible | Contrasts, Kinetic balance |
| Common partition ratio | One ratio γ of parental to hybrid variance, pooled by medians across genes, sets the partition for every gene of a quantity | Statistical power |
| Exchangeable means | Under the mean null, pooled or swapped cells differ only in dispersion, so the permutation null is approximate where dispersions differ | Null distributions |
| Null simulation model | Dispersion nulls and the midparent mean null are simulated from gamma expression factors, with a shared extrinsic factor for the two alleles of a cell and null dispersions set by projection onto the null | Null distributions |
| Extrinsic and intrinsic parts | Allele expression factors split into a shared extrinsic part and equal private intrinsic parts | Null distributions, Appendix B |
| Working weights | The weighted moment estimator uses second moments alone and needs no gamma assumption. Its estimating equation is unbiased only with true means and fixed weights | Appendix A |
| Expression factors and offset | Expression factors are independent of the offset | Appendix B |

## Summary of estimators

The table gives where each estimator is defined and where it is used, and the display gives each estimator as it is computed. A corrected estimate below zero is set to zero after its correction.

| Quantity | Defined in | Used in |
|---|---|---|
| Depth offset | Sequencing depth offset | Every fit |
| Mean rate | Estimating μ | Every mean contrast, fitted means M̂_i |
| Dispersion before any correction | Estimating φ | Parental dispersions directly, and φ̂*_Hc and φ̂*_He for the hybrid alleles |
| Informativeness class | Contrasts, Informative estimates | Gating of raw, log₂ and kinetic balance contrasts |
| Transcriptome share and its variance | Transcriptome share corrections in the hybrid | η̂, corrected allele means |
| Multiplier moments | Transcriptome share corrections, Appendix B | Transcriptome share correction of the allele dispersions |
| Allele mean | Transcriptome share corrections in the hybrid | Cis, trans and dominance mean contrasts |
| Allele dispersion | Transcriptome share corrections in the hybrid | Cis and trans dispersion contrasts, per-genome hybrid dispersion, kinetic quantities |
| Per-genome hybrid mean | Per-genome hybrid estimates, Mean | Dominance mean contrasts, hybrid CV² and kinetic quantities |
| Per-genome hybrid dispersion | Per-genome hybrid estimates, Dispersion | Dominance dispersion contrasts, hybrid CV² and kinetic quantities |
| Partition fraction | Statistical power | Cis and trans partitions of the hybrid cells |
| Kinetic quantities | NB model, Burst kinetics | Burst frequency, burst size and kinetic balance contrasts |
| Ploidy shift | Appendix B | Reported only |

$$
\begin{aligned}
&\text{Depth offset} && O_i = L_{ri}/D, \qquad S_{\mathrm{H}i} = L_{\mathrm{H}i}/D, \qquad O_i = S_{\mathrm{H}i}/2 \ \text{for each allele} \\
&\text{Mean rate} && \hat{\mu}_r = \frac{\sum_i y_{ri}}{\sum_i O_i}, \qquad \hat{M}_i = O_i\,\hat{\mu}_r \\
&\text{Dispersion} && \hat{\varphi}_r = \arg\max_{\varphi > 0} \ell^{\mathrm{CR}}_r(\varphi), \qquad \hat{\varphi}_r = 0 \ \text{if} \ \sum_i (y_{ri} - \hat{M}_i)^2 \le \sum_i \hat{M}_i - \frac{\sum_i \hat{M}_i^2}{\sum_i \hat{M}_i} \\
&\text{Informativeness class} && T_r = 2\left[\ell^{\mathrm{CR}}_r(\hat{\varphi}_r) - \ell^{\mathrm{CR}}_r(0)\right] > 2.71, \qquad \varphi^{\mathrm{up}}_r\,\bar{M}_r < 0.1 \\
&\text{Transcriptome share} && \hat{c}_i = \frac{L_{\mathrm{Hc}i}}{L_{\mathrm{H}i}}, \qquad \hat{\bar{c}} = \frac{1}{N_{\mathrm{H}}}\sum_i \hat{c}_i, \qquad \hat{\sigma}^2_c = \widehat{\mathrm{Var}}[\hat{c}_i] - \frac{1}{N_{\mathrm{H}}}\sum_i \frac{\hat{c}_i(1-\hat{c}_i)}{L_{\mathrm{H}i} - 1} \\
&\text{Multiplier moments} && \hat{\eta}_{\mathrm{Hc}} = \frac{\hat{\sigma}^2_c}{\hat{\bar{c}}^2}, \qquad \hat{\eta}_{\mathrm{He}} = \frac{\hat{\sigma}^2_c}{(1-\hat{\bar{c}})^2} \\
&\text{Allele dispersion} && \hat{\varphi}_{\mathrm{Hc}} = \frac{\hat{\varphi}^*_{\mathrm{Hc}} - \hat{\eta}_{\mathrm{Hc}}}{1 + \hat{\eta}_{\mathrm{Hc}}}, \qquad \hat{\varphi}_{\mathrm{He}} = \frac{\hat{\varphi}^*_{\mathrm{He}} - \hat{\eta}_{\mathrm{He}}}{1 + \hat{\eta}_{\mathrm{He}}} \\
&\text{Allele mean} && \hat{\mu}_{\mathrm{Hc}} = \frac{\hat{\mu}^*_{\mathrm{Hc}}}{2\hat{\bar{c}}}, \qquad \hat{\mu}_{\mathrm{He}} = \frac{\hat{\mu}^*_{\mathrm{He}}}{2(1-\hat{\bar{c}})} \\
&\text{Per-genome hybrid mean} && \hat{\mu}_{\mathrm{H}} = \tfrac12\left(\hat{\mu}_{\mathrm{Hc}} + \hat{\mu}_{\mathrm{He}}\right) \\
&\text{Per-genome hybrid dispersion} && \log_2\hat{\varphi}_{\mathrm{H}} = \tfrac12\left(\log_2\hat{\varphi}_{\mathrm{Hc}} + \log_2\hat{\varphi}_{\mathrm{He}}\right) \\
&\text{Partition fraction} && f_x = \frac{(2 + \hat{\gamma}) - \sqrt{\hat{\gamma}^2 + 4}}{2\hat{\gamma}}, \qquad x = \mu, \varphi, \kappa \\
&\text{Kinetic quantities} && \hat{\theta} = 1/\hat{\varphi}, \qquad \hat{\zeta} = \hat{\mu}\,\hat{\varphi}, \qquad \hat{\kappa} = \hat{\theta}/\hat{\zeta} \\
&\text{Ploidy shift (reported)} && \hat{\delta} = -\log_2\left[1 - 2\hat{\pi}(1-\hat{\pi})\,\hat{\omega}\right], \qquad \hat{\delta}_{\mathrm{raw}} = 2\hat{\pi}(1-\hat{\pi})\,\hat{\sigma}^2_q \ \ \text{(Appendix B)}
\end{aligned}
$$

## Appendix A. Weighted method of moments for the dispersion

The Cox–Reid adjusted NB likelihood estimate is the primary estimator of φ. As a sensitivity analysis we also estimated φ from the first two moments of the counts, without requiring them to follow an NB distribution. For the count Y_i of cell i, the mean is M_i = O_i μ. The variance can exceed the Poisson variance M_i by a term proportional to M_i².

$$
\mathrm{E}[Y_i] = O_i\,\mu = M_i, \qquad \mathrm{Var}[Y_i] = \mathrm{E}\left[(Y_i - M_i)^2\right] = M_i + \varphi\,M_i^2
$$

Subtracting the Poisson variance M_i from the squared deviation leaves a random variable whose expectation is this excess variance:

$$
\mathrm{E}\left[(Y_i - M_i)^2 - M_i\right] = \mathrm{E}\left[(Y_i - M_i)^2\right] - M_i = \left(M_i + \varphi\,M_i^2\right) - M_i = \varphi\,M_i^2
$$

If the means were known, we could produce a method of moments estimate of φ by replacing the expectation with the observed values. Each cell provides one such moment, and we averaged them with non-negative weights w_i so that more precise cells count more:

$$
\mathrm{E}\left[\sum_i w_i\left((Y_i - M_i)^2 - M_i\right)\right] = \varphi\sum_i w_i\,M_i^2 \\
\sum_i w_i\left[(y_i - M_i)^2 - M_i\right] = \tilde{\varphi}\sum_i w_i\,M_i^2 \\
\tilde{\varphi} = \frac{\sum_i w_i\left[(y_i - M_i)^2 - M_i\right]}{\sum_i w_i\,M_i^2}
$$

The means are not known, so we replaced M_i by its fitted value M̂_i in φ̃ to give a naive estimator φ̃(M̂).

$$
\tilde{\varphi}(\hat{M}) = \frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{M}_i\right]}{\sum_i w_i\,\hat{M}_i^2}
$$

The offsets are fixed, so the expectation of the fitted mean is the offset times the sum of the expected counts over the sum of the offsets, which equals M_i. The residual therefore has expectation zero, and the expectation of its square is its variance.

$$
\begin{aligned}
\mathrm{E}[\hat{M}_i] &= O_i\,\mathrm{E}[\hat{\mu}] = O_i\,\frac{\sum_j \mathrm{E}[Y_j]}{\sum_j O_j} = O_i\,\frac{\sum_j O_j\,\mu}{\sum_j O_j} = O_i\,\mu = M_i \\
\mathrm{E}[Y_i - \hat{M}_i] &= \mathrm{E}[Y_i] - \mathrm{E}[\hat{M}_i] = M_i - M_i = 0 \\
\mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] &= \left(\mathrm{E}[Y_i - \hat{M}_i]\right)^2 + \mathrm{Var}[Y_i - \hat{M}_i] = \mathrm{Var}[Y_i - \hat{M}_i]
\end{aligned}
$$

To find this variance, we define p_i as the share of the total offset for cell i. 

$$
\hat{M}_i = O_i\,\hat{\mu} = O_i\,\frac{\sum_j y_j}{\sum_j O_j} = p_i\sum_j y_j, \qquad p_i \equiv \frac{O_i}{\sum_j O_j}
$$

Pulling out the part that depends on Y_i from the residual produces:

$$
\begin{aligned}
Y_i - \hat{M}_i &= Y_i - p_i\sum_j Y_j = Y_i - p_i\,Y_i - p_i\sum_{j\neq i} Y_j = (1-p_i)\,Y_i - p_i\sum_{j\neq i} Y_j \\
\mathrm{Var}[Y_i - \hat{M}_i] &= (1-p_i)^2\,\mathrm{Var}[Y_i] + p_i^2\sum_{j\neq i}\mathrm{Var}[Y_j] = (1-2p_i)\,\mathrm{Var}[Y_i] + p_i^2\sum_j \mathrm{Var}[Y_j]
\end{aligned}
$$

Substituting in the variance and using the definition of p (Σ_j M_j = M_i/p_i) separates a Poisson part from a part proportional to φ. 

$$
\begin{aligned}
\mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] &= (1-2p_i)\left(M_i + \varphi M_i^2\right) + p_i^2\left(\sum_j M_j + \varphi\sum_j M_j^2\right) \\
&= (1-2p_i)\left(M_i + \varphi M_i^2\right) + p_i^2\left(\frac{M_i}{p_i} + \varphi\sum_j M_j^2\right) \\
&= (1-2p_i)\,M_i + (1-2p_i)\,\varphi M_i^2 + p_i\,M_i + \varphi\,p_i^2\sum_j M_j^2 \\
&= M_i\,(1-p_i) + \varphi\left[M_i^2\,(1-2p_i) + p_i^2\sum_j M_j^2\right] = A_i + \varphi\,B_i \\
A_i &\equiv M_i\,(1-p_i), \qquad B_i \equiv M_i^2\,(1-2p_i) + p_i^2\sum_j M_j^2
\end{aligned}
$$

For the numerator of the naive estimator we then have: 

$$
\begin{aligned}
\mathrm{E}\left[(Y_i - \hat{M}_i)^2 - \hat{M}_i\right] &= \mathrm{E}\left[(Y_i - \hat{M}_i)^2\right] - \mathrm{E}[\hat{M}_i] = A_i + \varphi\,B_i - M_i \\
&= M_i\,(1 - p_i) + \varphi\,B_i - M_i = \varphi\,B_i - p_i\,M_i
\end{aligned}
$$

With known means the numerator would have expectation φM_i². With fitted means it does not, so this estimator is biased. The bias term p_i M_i does not scale with φ. We therefore defined a corrected numerator that adds its estimate, p_i M̂_i, whose expectation is p_i M_i because the fitted mean has expectation M_i. The expectation of the corrected numerator is then:

$$
\mathrm{E}\left[(Y_i - \hat{M}_i)^2 - \hat{M}_i + p_i\,\hat{M}_i\right] = \left(\varphi\,B_i - p_i\,M_i\right) + p_i\,M_i = \varphi\,B_i
$$

The weighted sum of these terms minus φ Σ w_i B_i therefore has expectation zero at the true φ. Solving this for φ gives an unbiased estimator (φ̌) if B_i is known. 

$$
\mathrm{E}\left[\sum_i w_i\left((Y_i - \hat{M}_i)^2 - \hat{M}_i + p_i\,\hat{M}_i\right) - \varphi\sum_i w_i\,B_i\right] = 0 \\
\check{\varphi} \equiv \frac{\sum_i w_i\left(\left[(y_i - \hat{M}_i)^2 - \hat{M}_i\right] + p_i\,\hat{M}_i\right)}{\sum_i w_i\,\hat{B}_i} = \frac{\sum_i w_i\left((y_i - \hat{M}_i)^2 - \hat{A}_i\right)}{\sum_i w_i\,\hat{B}_i}, \\ \hat{A}_i \equiv \hat{M}_i\,(1 - p_i), \quad \hat{B}_i \equiv \hat{M}_i^2\,(1 - 2p_i) + p_i^2\sum_j \hat{M}_j^2
$$

A check marks this second intermediate estimate. It uses the fitted means and the plug-in B̂_i, which keeps it distinct from φ̃, written with the true means, and from the final estimate φ̂, which also corrects the bias of B̂_i.

In practice we evaluated A_i, B_i and the weights at M̂_i. Â_i is linear in the means when p_i is known and is therefore unbiased for A_i. However, B̂_i is quadratic in the means, so it has to be corrected.

$$
\begin{aligned}
\mathrm{E}\left[\hat{B}_i\right] &= \mathrm{E}\left[(1 - 2p_i)\,\hat{M}_i^2 + p_i^2\sum_j \hat{M}_j^2\right] = (1 - 2p_i)\,\mathrm{E}\left[\hat{M}_i^2\right] + p_i^2\sum_j \mathrm{E}\left[\hat{M}_j^2\right] \\
&= (1 - 2p_i)\left(M_i^2 + \mathrm{Var}[\hat{M}_i]\right) + p_i^2\sum_j\left(M_j^2 + \mathrm{Var}[\hat{M}_j]\right) \\
&= (1 - 2p_i)\,M_i^2 + p_i^2\sum_j M_j^2 + (1 - 2p_i)\,\mathrm{Var}[\hat{M}_i] + p_i^2\sum_j \mathrm{Var}[\hat{M}_j] \\
&= B_i + (1 - 2p_i)\,\mathrm{Var}[\hat{M}_i] + p_i^2\sum_j \mathrm{Var}[\hat{M}_j]
\end{aligned}
$$

The variances of the fitted means remain. Since M̂_j = O_j μ̂, its variance is O_j² Var[μ̂]. Substituting it collects the terms into a multiple of B_i.

$$
\begin{aligned}
\mathrm{Var}[\hat{M}_j] &= \mathrm{Var}[O_j\,\hat{\mu}] = O_j^2\,\mathrm{Var}[\hat{\mu}] = \frac{M_j^2}{\mu^2}\,\mathrm{Var}[\hat{\mu}] = M_j^2\,\nu, \qquad \nu \equiv \frac{\mathrm{Var}[\hat{\mu}]}{\mu^2} \\
\mathrm{E}\left[\hat{B}_i\right] &= B_i + \nu\left[(1 - 2p_i)\,M_i^2 + p_i^2\sum_j M_j^2\right] = (1 + \nu)\,B_i
\end{aligned}
$$

Because μ̂ is the total count divided by the total offset, and the counts are independent across cells, its variance follows from the variance relationship for each count.

$$
\begin{aligned}
\nu &= \frac{\mathrm{Var}[\hat{\mu}]}{\mu^2} = \frac{\mathrm{Var}\left[\sum_j Y_j\right]}{\mu^2\left(\sum_j O_j\right)^2} = \frac{\sum_j\left(M_j + \varphi\,M_j^2\right)}{M_{\mathrm{tot}}^2} = \frac{1}{M_{\mathrm{tot}}} + \varphi\,h \\ 
M_{\mathrm{tot}} &\equiv \sum_j M_j = \mu\sum_j O_j, 
\qquad h \equiv \sum_j p_j^2 = \frac{\sum_j M_j^2}{M_{\mathrm{tot}}^2} \\
\mathrm{E}\left[\hat{B}_i\right] &= (1 + \nu)\,B_i = \left(1 + \frac{1}{M_{\mathrm{tot}}} + \varphi\,h\right) B_i
\end{aligned}
$$

Dividing B̂_i by the coefficient out front then gives an unbiased estimate of B_i, if the coefficient is known. 

$$
\tilde{B}_i \equiv \frac{\hat{B}_i}{1 + 1/\sum_j y_j + \hat{\varphi}\,h}
$$

Plugging B̃_i into the denominator of φ̌ gives φ̂. 

$$
\begin{aligned}
\hat{\varphi} &= \frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{A}_i\right]}{\sum_i w_i\,\tilde{B}_i} = \left(1 + \frac{1}{\sum_j y_j} + \hat{\varphi}\,h\right)\frac{\sum_i w_i\left[(y_i - \hat{M}_i)^2 - \hat{A}_i\right]}{\sum_i w_i\,\hat{B}_i} \\
&= \left(1 + \frac{1}{\sum_j y_j} + \hat{\varphi}\,h\right)\check{\varphi} \\
&= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} + h\,\check{\varphi}\,\hat{\varphi}
\end{aligned}
$$

Collecting the terms in φ̂ and solving gives the corrected estimator.

$$
\begin{aligned}
\hat{\varphi} - h\,\check{\varphi}\,\hat{\varphi} &= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} \\
\hat{\varphi}\left(1 - h\,\check{\varphi}\right) &= \left(1 + \frac{1}{\sum_j y_j}\right)\check{\varphi} \\
\hat{\varphi} &= \frac{\left(1 + 1/\sum_j y_j\right)\check{\varphi}}{1 - h\,\check{\varphi}}
\end{aligned}
$$

Finally, we chose the weights w_i. φ̂ is approximately unbiased for fixed weights, so the weights mainly affect the precision of φ̂ and only slightly its remaining bias. We minimized the variance of φ̃ instead, because with M_i and w_i known it can be written exactly.

$$
\begin{aligned}
\mathrm{Var}[\tilde{\varphi}] &= \mathrm{Var}\left[\frac{\sum_i w_i\left[(Y_i - M_i)^2 - M_i\right]}{\sum_i w_i\,M_i^2}\right] = \frac{\mathrm{Var}\left[\sum_i w_i\left[(Y_i - M_i)^2 - M_i\right]\right]}{\left(\sum_i w_i\,M_i^2\right)^2} \\
&= \frac{\sum_i \mathrm{Var}\left[w_i\left[(Y_i - M_i)^2 - M_i\right]\right]}{\left(\sum_i w_i\,M_i^2\right)^2} = \frac{\sum_i \mathrm{Var}\left[w_i\,(Y_i - M_i)^2 - w_i\,M_i\right]}{\left(\sum_i w_i\,M_i^2\right)^2} \\
&= \frac{\sum_i \mathrm{Var}\left[w_i\,(Y_i - M_i)^2\right]}{\left(\sum_i w_i\,M_i^2\right)^2} = \frac{\sum_i w_i^2\,\mathrm{Var}\left[(Y_i - M_i)^2\right]}{\left(\sum_i w_i\,M_i^2\right)^2}
\end{aligned}
$$

What remains unknown is V_i, the variance of the squared deviation of cell i. The variance of any random variable is the mean of its square minus the square of its mean. The deviation has expectation zero and we write E[(Y_i − M_i)⁴] as τ Var[Y_i]², where τ is the kurtosis of the counts of cell i. We assumed that τ is the same for every cell. 

$$
\begin{aligned}
V_i\equiv\mathrm{Var}\left[(Y_i - M_i)^2\right] &= \mathrm{E}\left[\left((Y_i - M_i)^2\right)^2\right] - \left(\mathrm{E}\left[(Y_i - M_i)^2\right]\right)^2 \\ 
&= \mathrm{E}\left[(Y_i - M_i)^4\right] - \mathrm{Var}[Y_i]^2  
= \tau\,\mathrm{Var}[Y_i]^2 - \mathrm{Var}[Y_i]^2 \\ 
&= (\tau - 1)\,\mathrm{Var}[Y_i]^2 = (\tau - 1)\,M_i^2(1 + \varphi M_i)^2
\\ \tau &\equiv \frac{\mathrm{E}\left[(Y_i - M_i)^4\right]}{\mathrm{Var}[Y_i]^2}
\end{aligned}
$$

To minimize the variance, we set the derivative with respect to the weights to zero and solved. 

$$
\begin{aligned}
\frac{\partial}{\partial w_i}\,\frac{\sum_j w_j^2\,V_j}{\left(\sum_j w_j\,M_j^2\right)^2} &= \frac{2\,w_i\,V_i}{\left(\sum_j w_j\,M_j^2\right)^2} - \frac{2\,M_i^2\sum_j w_j^2\,V_j}{\left(\sum_j w_j\,M_j^2\right)^3} = 0 \\
\Longrightarrow\quad w_i &= \frac{\sum_j w_j^2\,V_j}{\sum_j w_j\,M_j^2}\cdot\frac{M_i^2}{V_i}
\end{aligned}
$$

The first factor in the solution is a sum over all cells, so it is the same for every cell and only scales the weights. Similarly, τ − 1 is the same for every cell, so we dropped both factors.

$$
\begin{aligned}
w_i &= \frac{\sum_j w_j^2\,V_j}{\sum_j w_j\,M_j^2}\cdot\frac{M_i^2}{V_i} \propto \frac{M_i^2}{V_i} = \frac{M_i^2}{(\tau - 1)\,M_i^2\,(1 + \varphi M_i)^2} \\
&= \frac{1}{(\tau - 1)(1 + \varphi M_i)^2} \propto \frac{1}{(1 + \varphi M_i)^2}
\end{aligned}
$$

These weights are the optimum for φ̃ when the kurtosis is the same for every cell. We used them for φ̂, evaluated at M̂_i, which departs from the derivation in three ways. First, the fitted means lower the variance, because they remove the part of the numerator that moves with the total count. Second, the kurtosis does not have to be the same for every cell, for example for NB counts, where it falls as the mean rises. Third, evaluating the weights at fitted values makes them random, which adds bias.

These assumptions have little consequence at the cell numbers in this study. We tested this by simulation at N = 800 and N = 1600, with M of 0.5, 2 or 5, φ of 0.3 or 1, offset CV of 0.5 or 1, and either NB or Poisson-lognormal counts. We compared the weights from this derivation with more complete weights that use the true fourth moment of each cell and include the effect of the fitted means. Both used the true M, φ and count distribution. The weights from this derivation were nearly as precise. The median standard deviation of φ̂ was 0.5% higher at N = 800 and 0.3% higher at N = 1600. Evaluating them at the fitted M̂_i and an iterated φ̂ added at most 0.3% more at the median. The loss was concentrated when the true counts were NB with high mean and dispersion. Almost all of this loss comes from leaving out the effect of the fitted means. Bias differed from that of the more complete weights by at most 0.006, and it was at most 0.06 of the standard deviation for every choice.

The weights depend on φ, so we started from φ = 0, calculated the weights at M̂_i, updated φ̂ with the estimating equation, restricted it to nonnegative values, and repeated until convergence or a maximum number of iterations. We recorded the convergence status, whether the estimate sat at the boundary, and the untruncated estimate.

## Appendix B. The ploidy shift

Summing the two alleles of a hybrid cell averages their private noise, so the summed counts are less dispersed than a single genome. We derive this ploidy shift, δ, and its estimator, including the correction of the estimator for variation in the transcriptome share, and the transcriptome share term in the dispersion of the summed counts. We reported δ as a description of how much noise summing the two alleles hides.

**Summed counts.** Summing two alleles within a cell averages their private fluctuations, so the summed hybrid appears quieter than either allele alone. We derive the size of this effect from the decomposition of dispersion into intrinsic and extrinsic parts (Elowitz et al. 2002). Setting aside variation in the transcriptome share, the Hc count is Poisson around its library size times c̄, its expression per genome, and its expression factor a. The two allele counts of a cell are independent given these expression factors, and the sum of two independent Poisson counts is Poisson with the summed mean. The summed mean is therefore S_H times the sum of the two allele terms. 

$$
Y_{\mathrm{H}i} \mid a_{\mathrm{Hc}i}, a_{\mathrm{He}i} \sim \mathrm{Poisson}\left(S_{\mathrm{H}i}\left[\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i} + (1-\bar{c})\,\mu_{\mathrm{He}}\,a_{\mathrm{He}i}\right]\right)
$$

We write μ_Σ = c̄μ_Hc + (1 − c̄)μ_He for the mean rate of the summed counts and define π as the fraction of it contributed by the Hc allele. Rearranging gives a weighted sum for the expression factor of the summed hybrid, a_Hi. Each expression factor has mean one and the weights sum to one, so a_Hi has mean one.

$$
\begin{aligned}
Y_{\mathrm{H}i} \mid a_{\mathrm{Hc}i}, a_{\mathrm{He}i} &\sim \mathrm{Poisson}\left(S_{\mathrm{H}i}\,\mu_{\Sigma}\left[\pi\,a_{\mathrm{Hc}i} + (1-\pi)\,a_{\mathrm{He}i}\right]\right), \qquad \mu_{\Sigma} \equiv \bar{c}\,\mu_{\mathrm{Hc}} + (1-\bar{c})\,\mu_{\mathrm{He}}, \qquad \pi = \frac{\bar{c}\,\mu_{\mathrm{Hc}}}{\mu_{\Sigma}} \\
&= \mathrm{Poisson}\left(S_{\mathrm{H}i}\,\mu_{\Sigma}\,a_{\mathrm{H}i}\right), \qquad a_{\mathrm{H}i} \equiv \pi\,a_{\mathrm{Hc}i} + (1-\pi)\,a_{\mathrm{He}i} \\
\mathrm{E}[a_{\mathrm{H}i}] &= \pi\,\mathrm{E}[a_{\mathrm{Hc}i}] + (1-\pi)\,\mathrm{E}[a_{\mathrm{He}i}] = \pi + (1-\pi) = 1
\end{aligned}
$$

The NB fit of the summed counts treats a_Hi as its expression factor, and we call the variance of that factor φ_Σ. 

**Intrinsic and extrinsic noise.** In the hybrid, each allele has an intrinsic deviation (q) that arises from its own stochastic bursting, and the two alleles share an extrinsic deviation (d) due to shared cell size, cell-cycle state, and pool of transcription factors. The deviations are mutually independent and both deviations have a mean of zero. Since the dispersion is the variance of the expression factor, we have, assuming the two alleles have the same intrinsic variance:

$$
a_{\mathrm{Hc}i} = 1 + d_i + q_{\mathrm{Hc}i}, \qquad a_{\mathrm{He}i} = 1 + d_i + q_{\mathrm{He}i} \\
\mathrm{Var}[d_i] = \sigma^2_d, \qquad \mathrm{Var}[q_{\mathrm{Hc}i}] = \mathrm{Var}[q_{\mathrm{He}i}] = \sigma^2_q \\
\varphi_{\mathrm{Hc}} = \mathrm{Var}[a_{\mathrm{Hc}i}] = \mathrm{Var}[d_i] + \mathrm{Var}[q_{\mathrm{Hc}i}] = \sigma^2_d + \sigma^2_q = \varphi_{\mathrm{He}}
$$

The summed expression factor is a weighted average of the two allele expression factors. The weights sum to one, so the extrinsic deviation keeps a coefficient of one, and the independence of the deviations lets the variances add. This gives the dispersion of the summed hybrid.

$$
\begin{aligned}
\varphi_{\Sigma} &= \mathrm{Var}\left[\pi\,a_{\mathrm{Hc}i} + (1-\pi)\,a_{\mathrm{He}i}\right] = \mathrm{Var}\left[1 + d_i + \pi\,q_{\mathrm{Hc}i} + (1-\pi)\,q_{\mathrm{He}i}\right] \\
&= \sigma^2_d + \pi^2\sigma^2_q + (1-\pi)^2\sigma^2_q = \sigma^2_d + \left[\pi^2 + (1-\pi)^2\right]\sigma^2_q
\end{aligned}
$$

The extrinsic variance passes through in full, whereas the intrinsic variance enters with weight π² + (1 − π)². Expanding the square and completing it around π = 1/2 gives:

$$
\pi^2 + (1-\pi)^2 = 1 - 2\pi + 2\pi^2 = \tfrac12 + 2\left(\pi - \tfrac12\right)^2 \\
\tfrac12 \le \pi^2 + (1-\pi)^2 < 1 \qquad \text{for } 0 < \pi < 1
$$

Whenever both alleles contribute counts, the weight on the intrinsic variance is between one half and one. It equals one half at equal counts and approaches one as a single allele dominates. The summed hybrid therefore carries less intrinsic variance than a single allele, and this averaging can be estimated for each gene.

**Ploidy shift.** A parent carries one genome, so its expression factor has the same form as that of one allele, with the same extrinsic and intrinsic variances. We write φ_H for the dispersion of a single genome. The difference between φ_H and the dispersion of the summed hybrid is the averaging component. Dividing this difference by φ_H gives the fraction of dispersion lost in the summed hybrid. This fraction depends only on a single factor ω, the intrinsic variance share of φ_H.

$$
\begin{aligned}
\varphi_{\mathrm{H}} &= \mathrm{Var}\left[1 + d_i + q_{\mathrm{Hc}i}\right] = \sigma^2_d + \sigma^2_q \\
\varphi_{\mathrm{H}} - \varphi_{\Sigma} &= \left(\sigma^2_d + \sigma^2_q\right) - \left(\sigma^2_d + \left[\pi^2 + (1-\pi)^2\right]\sigma^2_q\right) = \left[1 - \pi^2 - (1-\pi)^2\right]\sigma^2_q \\
\frac{\varphi_{\Sigma}}{\varphi_{\mathrm{H}}} &= 1 - \frac{\varphi_{\mathrm{H}} - \varphi_{\Sigma}}{\varphi_{\mathrm{H}}} = 1 - \left[1 - \pi^2 - (1-\pi)^2\right]\omega, \qquad \omega \equiv \frac{\sigma^2_q}{\sigma^2_d + \sigma^2_q}
\end{aligned}
$$

The log₂ of the inverse of this ratio is the expected fall in log₂ dispersion from summing.

$$
\delta \equiv \log_2\frac{\varphi_{\mathrm{H}}}{\varphi_{\Sigma}} = -\log_2\frac{\varphi_{\Sigma}}{\varphi_{\mathrm{H}}} = -\log_2\left(1 - \left[1 - \pi^2 - (1-\pi)^2\right]\omega\right)
$$

The weight 1 − π² − (1 − π)² equals 2π(1 − π) and lies between zero and one half, and ω lies between zero and one, so the argument of the logarithm lies between one half and one and δ lies between 0 and 1. The shift equals 1 for two equally expressed alleles with purely intrinsic noise, the case in which averaging halves the noise. It falls to 0 when noise is entirely extrinsic, or when one allele dominates expression and leaves nothing to average.

On the raw scale, the same averaging lowers the dispersion of the summed hybrid below that of a single genome by the difference of the two dispersions. It is additive and proportional to the intrinsic variance, so it is larger in absolute terms for genes with more intrinsic noise.

$$
\delta_{\mathrm{raw}} \equiv \varphi_{\mathrm{H}} - \varphi_{\Sigma} = \left[1 - \pi^2 - (1-\pi)^2\right]\sigma^2_q
$$

**Estimating the ploidy shift.** We estimated the ploidy shift from estimators of π, σ²_d and σ²_q. For the weight π, we replaced each mean rate by its fitted value, the ratio of its summed counts to the summed offsets.

$$
\hat{\pi} = \frac{\hat{\mu}^*_{\mathrm{Hc}}}{\hat{\mu}^*_{\mathrm{Hc}} + \hat{\mu}^*_{\mathrm{He}}}
$$

This estimator is a ratio of sums, so it is consistent but biased. The bias comes from the joint fluctuation of the two sums around their means. Each sum runs over N_H cells, so the bias is 𝒪(1/N_H), and we dropped it.

To estimate the extrinsic variance, we note that it is also the covariance of the two alleles' expression factors, because the extrinsic deviation is the only part that they share. Each expression factor has mean one, so expanding the covariance gives:

$$
\begin{aligned}
\sigma^2_d &= \mathrm{Var}[d_i] = \mathrm{Cov}[a_{\mathrm{Hc}i}, a_{\mathrm{He}i}] = \mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}] - \mathrm{E}[a_{\mathrm{Hc}i}]\,\mathrm{E}[a_{\mathrm{He}i}] = \mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}] - 1
\end{aligned}
$$

The intrinsic variance is the same for both alleles and the expected deviations have mean zero, so expanding the variances gives: 

$$
\begin{aligned}
\sigma^2_q &= \mathrm{Var}[q_{\mathrm{Hc}i}] = \tfrac12\left[\mathrm{Var}[q_{\mathrm{Hc}i}] + \mathrm{Var}[q_{\mathrm{He}i}]\right] \\
&= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] - \mathrm{E}[q_{\mathrm{Hc}i}]^2 + \mathrm{E}[q_{\mathrm{He}i}^2] - \mathrm{E}[q_{\mathrm{He}i}]^2\right) 
= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] + \mathrm{E}[q_{\mathrm{He}i}^2]\right) 
\end{aligned}
$$

The two intrinsic deviations are independent, so E[q_Hc q_He] is zero, and subtracting twice this cross moment completes the square. Within the square, we add and subtract one plus the extrinsic deviation, and then substitute with the allele expression factors. Expanding and simplifying gives:

$$
\begin{aligned}
\sigma^2_q &= \tfrac12\left(\mathrm{E}[q_{\mathrm{Hc}i}^2] + \mathrm{E}[q_{\mathrm{He}i}^2] - 2\,\mathrm{E}[q_{\mathrm{Hc}i}\,q_{\mathrm{He}i}]\right) = \tfrac12\,\mathrm{E}\left[\left(q_{\mathrm{Hc}i} - q_{\mathrm{He}i}\right)^2\right]\\
&= \tfrac12\,\mathrm{E}\left[\left((1 + d_i + q_{\mathrm{Hc}i}) - (1 + d_i + q_{\mathrm{He}i})\right)^2\right] = \tfrac12\,\mathrm{E}\left[\left(a_{\mathrm{Hc}i} - a_{\mathrm{He}i}\right)^2\right] \\
&= \tfrac12\,\mathrm{E}\left[a_{\mathrm{Hc}i}^2 - 2\,a_{\mathrm{Hc}i}\,a_{\mathrm{He}i} + a_{\mathrm{He}i}^2\right] = \tfrac12\left(\mathrm{E}[a_{\mathrm{Hc}i}^2] + \mathrm{E}[a_{\mathrm{He}i}^2]\right) - \mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}]
\end{aligned}
$$

This leaves the intrinsic and extrinsic variances written as moments of the expression factors. To simplify, we multiply and divide each expression factor by its library size, transcriptome share, and per-genome rate, and replace each product with the conditional expected count. Because the two counts are independent given their expression factors, the product of the conditional expectations is the conditional expectation of the product. The law of total expectation then removes the conditioning.

$$
\begin{aligned}
\mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}] &= \mathrm{E}\left[\frac{S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i}}{S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}} \cdot \frac{S_{\mathrm{H}i}(1-\bar{c})\,\mu_{\mathrm{He}}\,a_{\mathrm{He}i}}{S_{\mathrm{H}i}(1-\bar{c})\,\mu_{\mathrm{He}}}\right] = \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i} \mid a_{\mathrm{Hc}i}]}{S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}} \cdot \frac{\mathrm{E}[Y_{\mathrm{He}i} \mid a_{\mathrm{He}i}]}{S_{\mathrm{H}i}(1-\bar{c})\,\mu_{\mathrm{He}}}\right] \\
&= \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i} \mid a_{\mathrm{Hc}i}, a_{\mathrm{He}i}]}{S_{\mathrm{H}i}^2\,\bar{c}(1-\bar{c})\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] = \mathrm{E}\left[\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{S_{\mathrm{H}i}^2\,\bar{c}(1-\bar{c})\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\,\middle|\, a_{\mathrm{Hc}i}, a_{\mathrm{He}i}\right]\right] = \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{S_{\mathrm{H}i}^2\,\bar{c}(1-\bar{c})\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right]
\end{aligned}
$$

For the square moment, we write a_Hc in terms of the conditional Poisson mean of its count. A Poisson count has variance equal to its mean, so its conditional mean square is the mean plus the squared mean. Subtracting the conditional mean leaves the squared mean. The denominator is fixed, and the law of total expectation removes the conditioning, giving:

$$
\begin{aligned}
\mathrm{E}[a_{\mathrm{Hc}i}^2] &= \mathrm{E}\left[\frac{\left(S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i}\right)^2}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\right] = \mathrm{E}\left[\frac{\left(S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i} + \left(S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i}\right)^2\right) - S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,a_{\mathrm{Hc}i}}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\right] = \mathrm{E}\left[\frac{\mathrm{E}[Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1) \mid a_{\mathrm{Hc}i}]}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\right] \\
&= \mathrm{E}\left[\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\,\middle|\, a_{\mathrm{Hc}i}\right]\right] = \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\right], \qquad \mathrm{E}[a_{\mathrm{He}i}^2] = \mathrm{E}\left[\frac{Y_{\mathrm{He}i}(Y_{\mathrm{He}i} - 1)}{S_{\mathrm{H}i}^2\,(1-\bar{c})^2\,\mu_{\mathrm{He}}^2}\right]
\end{aligned}
$$

Substituting these moments into the definitions gives σ²_d and σ²_q in terms of the counts, the library sizes, the transcriptome share and the per-genome rates. We then replaced each expectation by the mean over cells and each expected count by its estimate. These estimators set aside variation in the transcriptome share, which biases them, so we marked them with a star and corrected them.

$$
\begin{aligned}
\sigma^2_d &= \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{S_{\mathrm{H}i}^2\,\bar{c}(1-\bar{c})\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] - 1 \quad\Rightarrow\quad \hat{\sigma}^{2*}_d = \frac{1}{N_{\mathrm{H}}}\sum_i \frac{y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_i^2\,\hat{\mu}^*_{\mathrm{Hc}}\,\hat{\mu}^*_{\mathrm{He}}} - 1 \\
\sigma^2_q &= \tfrac12\left(\mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}(Y_{\mathrm{Hc}i} - 1)}{S_{\mathrm{H}i}^2\,\bar{c}^2\,\mu_{\mathrm{Hc}}^2}\right] + \mathrm{E}\left[\frac{Y_{\mathrm{He}i}(Y_{\mathrm{He}i} - 1)}{S_{\mathrm{H}i}^2\,(1-\bar{c})^2\,\mu_{\mathrm{He}}^2}\right]\right) - \mathrm{E}\left[\frac{Y_{\mathrm{Hc}i}\,Y_{\mathrm{He}i}}{S_{\mathrm{H}i}^2\,\bar{c}(1-\bar{c})\,\mu_{\mathrm{Hc}}\,\mu_{\mathrm{He}}}\right] \\
\hat{\sigma}^{2*}_q &= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}(y_{\mathrm{Hc}i} - 1)}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} + \frac{y_{\mathrm{He}i}(y_{\mathrm{He}i} - 1)}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2}\right] - \frac{1}{N_{\mathrm{H}}}\sum_i \frac{y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_i^2\,\hat{\mu}^*_{\mathrm{Hc}}\,\hat{\mu}^*_{\mathrm{He}}}
\end{aligned}
$$

Expanding the factorial moments and completing the square in the scaled allele counts rewrites the estimator of the intrinsic variance.

$$
\begin{aligned}
\hat{\sigma}^{2*}_q 
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}^2}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} - \frac{y_{\mathrm{Hc}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} + \frac{y_{\mathrm{He}i}^2}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2} - \frac{y_{\mathrm{He}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2} - \frac{2\,y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_i^2\,\hat{\mu}^*_{\mathrm{Hc}}\,\hat{\mu}^*_{\mathrm{He}}}\right] \\
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left[\frac{y_{\mathrm{Hc}i}^2}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} - \frac{2\,y_{\mathrm{Hc}i}\,y_{\mathrm{He}i}}{O_i^2\,\hat{\mu}^*_{\mathrm{Hc}}\,\hat{\mu}^*_{\mathrm{He}}} + \frac{y_{\mathrm{He}i}^2}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2}\right] - \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} + \frac{y_{\mathrm{He}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2}\right) \\
&= \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_i\,\hat{\mu}^*_{\mathrm{Hc}}} - \frac{y_{\mathrm{He}i}}{O_i\,\hat{\mu}^*_{\mathrm{He}}}\right)^2 - \frac{1}{2N_{\mathrm{H}}}\sum_i \left(\frac{y_{\mathrm{Hc}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{Hc}}}^2} + \frac{y_{\mathrm{He}i}}{O_i^2\,{\hat{\mu}^*_{\mathrm{He}}}^2}\right)
\end{aligned}
$$

The estimator of the intrinsic variance is therefore half the mean squared difference of the scaled allele counts minus their average Poisson noise. The estimator of the extrinsic variance is the mean product of the scaled allele counts minus one. Each scaled count has expectation one, so this estimates their covariance, and it needs no Poisson correction because the two counts are independent given their expression factors.

With the true mean rates, every term is an unbiased estimate of the corresponding moment of the expression factors, so both estimators are unbiased. As for π, the fitted mean rates are ratios of sums. Replacing each mean rate by its fitted value changes each term with a bias of order 1/N_H, which we dropped.

$$
\mathrm{E}\left[\hat{\sigma}^{2*}_d\right] - \sigma^2_d = \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right), \qquad \mathrm{E}\left[\hat{\sigma}^{2*}_q\right] - \sigma^2_q = \mathcal{O}\left(\frac{1}{N_{\mathrm{H}}}\right)
$$

**Transcriptome share correction.** The star estimators set aside variation in the transcriptome share. Their limits when the transcriptome share varies are derived in *Transcriptome share and the estimated shift*. Solving those limits for σ²_d and σ²_q gives the corrected estimators. They use η̂_r (*T**ranscriptome share corrections** in the hybrid*) and √(η̂_Hc η̂_He) (*Variance estimates*). We set any corrected variance below zero to zero.

$$
\begin{aligned}
\hat{\sigma}^2_d &= \frac{1 + \hat{\sigma}^{2*}_d}{1 - \sqrt{\hat{\eta}_{\mathrm{Hc}}\,\hat{\eta}_{\mathrm{He}}}} - 1, \qquad \sqrt{\hat{\eta}_{\mathrm{Hc}}\,\hat{\eta}_{\mathrm{He}}} = \frac{\hat{\sigma}^2_c}{\hat{\bar{c}}\left(1-\hat{\bar{c}}\right)} \\
\hat{\sigma}^2_q &= \frac{\hat{\sigma}^{2*}_q - \left[\tfrac12\left(\hat{\eta}_{\mathrm{Hc}} + \hat{\eta}_{\mathrm{He}}\right) + \sqrt{\hat{\eta}_{\mathrm{Hc}}\,\hat{\eta}_{\mathrm{He}}}\right]\left(1 + \hat{\sigma}^2_d\right)}{1 + \tfrac12\left(\hat{\eta}_{\mathrm{Hc}} + \hat{\eta}_{\mathrm{He}}\right)}
\end{aligned}
$$

Plugging these estimates into the definition of δ̂ gives the log₂ ratio of the estimated dispersion of a single genome to the estimated dispersion of the summed hybrid.

$$
\begin{aligned}
\hat{\delta} &\equiv \log_2\frac{\hat{\sigma}^2_d + \hat{\sigma}^2_q}{\hat{\sigma}^2_d + \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q}
\end{aligned}
$$

Both truncated variances are nonnegative, so δ̂ lies between 0 and 1. A zero estimate of the intrinsic variance gives δ̂ = 0. When both variance estimates are zero the ratio is undefined, and we set δ̂ to zero. Truncation therefore affects δ̂ only for genes with a variance estimate near zero.

The estimated raw shift δ̂_raw is the difference between the estimated dispersion of a single genome and that of the summed hybrid.

$$
\begin{aligned}
\hat{\delta}_{\mathrm{raw}} &\equiv \left(\hat{\sigma}^2_d + \hat{\sigma}^2_q\right) - \left(\hat{\sigma}^2_d + \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q\right) \\
&= \hat{\sigma}^2_d - \hat{\sigma}^2_d + \hat{\sigma}^2_q - \left[\hat{\pi}^2 + (1-\hat{\pi})^2\right]\hat{\sigma}^2_q = \left[1 - \hat{\pi}^2 - (1-\hat{\pi})^2\right]\hat{\sigma}^2_q
\end{aligned}
$$

**Transcriptome share and the estimated shift.** Variation in the transcriptome share biases the estimators of π, σ²_d and σ²_q, and through them the estimated shifts.

**Weight estimate.** Each allele count carries its transcriptome share multiplier, t_Hci = c_i/c̄ for Hc and t_Hei = (1 − c_i)/(1 − c̄) for He, which has mean one and is independent of depth and expression. The expected counts therefore keep the form used for the weight estimate, and the weight estimate is unaffected.

$$
\begin{aligned}
\mathrm{E}[Y_{\mathrm{Hc}i}] &= S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}\,\mathrm{E}[t_{\mathrm{Hc}i}] = S_{\mathrm{H}i}\,\bar{c}\,\mu_{\mathrm{Hc}}, \qquad \mathrm{E}[Y_{\mathrm{He}i}] = S_{\mathrm{H}i}\,(1-\bar{c})\,\mu_{\mathrm{He}}\,\mathrm{E}[t_{\mathrm{He}i}] = S_{\mathrm{H}i}\,(1-\bar{c})\,\mu_{\mathrm{He}} \\
\hat{\pi} &\to \pi
\end{aligned}
$$

**Variance estimates.** With transcriptome share variation, each allele count is scaled by its transcriptome share multiplier as well as its expression factor, so the count moments in the estimators pick up moments of the multipliers. The multipliers have mean one and are independent of the expression factors, so expectations of products factorize. Their second moments are E[t_r²] = 1 + η_r (*T**ranscriptome share corrections** in the hybrid*), and because the two transcriptome shares sum to one, E[t_Hc t_He] = E[c_i(1 − c_i)]/[c̄(1 − c̄)] = [c̄(1 − c̄) − σ²_c]/[c̄(1 − c̄)] = 1 − σ²_c/[c̄(1 − c̄)] = 1 − √(η_Hc η_He), because σ²_c/[c̄(1 − c̄)] is the geometric mean of η_Hc and η_He. The expression factors have E[a_r²] = 1 + φ_H and E[a_Hc a_He] = 1 + σ²_d. The estimators therefore converge to

$$
\begin{aligned}
\hat{\sigma}^{2*}_d &\to \mathrm{E}[t_{\mathrm{Hc}i}\,t_{\mathrm{He}i}]\,\mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}] - 1 = (1 - \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}})\left(1 + \sigma^2_d\right) - 1 = \sigma^2_d - \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\left(1 + \sigma^2_d\right) \\
\hat{\sigma}^{2*}_q &\to \tfrac12\left(\mathrm{E}[t_{\mathrm{Hc}i}^2]\,\mathrm{E}[a_{\mathrm{Hc}i}^2] + \mathrm{E}[t_{\mathrm{He}i}^2]\,\mathrm{E}[a_{\mathrm{He}i}^2]\right) - \mathrm{E}[t_{\mathrm{Hc}i}\,t_{\mathrm{He}i}]\,\mathrm{E}[a_{\mathrm{Hc}i}\,a_{\mathrm{He}i}] \\
&= \tfrac12\left(2 + \eta_{\mathrm{Hc}} + \eta_{\mathrm{He}}\right)\left(1 + \varphi_{\mathrm{H}}\right) - (1 - \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}})\left(1 + \sigma^2_d\right) = \sigma^2_q + \tfrac12\left(1 + \varphi_{\mathrm{H}}\right)\sigma^2_c\left[\frac{1}{\bar{c}^2} + \frac{1}{(1-\bar{c})^2}\right] + \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\left(1 + \sigma^2_d\right)
\end{aligned}
$$

**Bias of the shifts.** Since π̂ converges to π, substituting these limits into the raw shift built from the star estimators, δ̂*_raw, gives the bias that the correction removes. The log₂ shift inherits it through the ratio.

$$
\hat{\delta}^*_{\mathrm{raw}} \to \left[1 - \pi^2 - (1-\pi)^2\right]\left(\sigma^2_q + \tfrac12\left(1 + \varphi_{\mathrm{H}}\right)\sigma^2_c\left[\frac{1}{\bar{c}^2} + \frac{1}{(1-\bar{c})^2}\right] + \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\left(1 + \sigma^2_d\right)\right)
$$

Both biases raise the estimated intrinsic fraction and therefore the shift. Solving the two limits above for σ²_d and σ²_q gives the corrected estimators of this appendix, which remove both biases. A simulation of this model with Poisson counts, shifted exponential extrinsic and intrinsic deviations, 8 × 10⁶ cells, and the true transcriptome share variance and mean rates reproduced both limits to within 0.0005, and the corrected estimators recovered σ²_d and σ²_q to within 0.0005.

**Transcriptome share and the dispersion of the summed counts.** Variation in the transcriptome share also adds a term to the dispersion of the summed hybrid, and a correction removes it.

Each allele count is scaled by a transcriptome share multiplier t, defined in *T**ranscriptome share corrections** in the hybrid*, that has mean one and is independent of the expression factors. The sum of the two allele counts is again Poisson around the sum of the two allele means, so the rewriting of the summed counts at the start of this appendix applies with the multipliers included. The summed expression factor that the NB fit sees is then a*_Hi.

$$
\begin{aligned}
Y_{\mathrm{H}i} \mid c_i, a_{\mathrm{Hc}i}, a_{\mathrm{He}i} &\sim \mathrm{Poisson}\left(S_{\mathrm{H}i}\left[\bar{c}\,\mu_{\mathrm{Hc}}\,t_{\mathrm{Hc}i}\,a_{\mathrm{Hc}i} + (1-\bar{c})\,\mu_{\mathrm{He}}\,t_{\mathrm{He}i}\,a_{\mathrm{He}i}\right]\right) \\
&= \mathrm{Poisson}\left(S_{\mathrm{H}i}\,\mu_{\Sigma}\,a^*_{\mathrm{H}i}\right), \qquad a^*_{\mathrm{H}i} \equiv \pi\,t_{\mathrm{Hc}i}\,a_{\mathrm{Hc}i} + (1-\pi)\,t_{\mathrm{He}i}\,a_{\mathrm{He}i}
\end{aligned}
$$

Each multiplier and expression factor has mean one and they are independent, so a*_Hi has mean one. The fitted dispersion estimates its variance. The multipliers vary among cells only through c_i, so we split the variance into a part between transcriptome shares and a part within a transcriptome share.

$$
\begin{aligned}
\mathrm{E}[a^*_{\mathrm{H}i}] &= \pi\,\mathrm{E}[t_{\mathrm{Hc}i}]\,\mathrm{E}[a_{\mathrm{Hc}i}] + (1-\pi)\,\mathrm{E}[t_{\mathrm{He}i}]\,\mathrm{E}[a_{\mathrm{He}i}] = 1 \\
\mathrm{Var}[a^*_{\mathrm{H}i}] &= \underbrace{\mathrm{Var}\left[\mathrm{E}[a^*_{\mathrm{H}i} \mid c_i]\right]}_{\text{between transcriptome shares}} + \underbrace{\mathrm{E}\left[\mathrm{Var}[a^*_{\mathrm{H}i} \mid c_i]\right]}_{\text{within a transcriptome share}}
\end{aligned}
$$

For the between part, each expression factor has mean one and is independent of the transcriptome share, so conditioning on c_i leaves only the multipliers. Its variance follows from the moments of the multipliers, Var[t_r] = η_r and Cov[t_Hc, t_He] = −√(η_Hc η_He) (*Variance estimates*):

$$
\begin{aligned}
\mathrm{E}[a^*_{\mathrm{H}i} \mid c_i] &= \pi\,t_{\mathrm{Hc}i} + (1-\pi)\,t_{\mathrm{He}i} \\
\mathrm{Var}\left[\mathrm{E}[a^*_{\mathrm{H}i} \mid c_i]\right] &= \pi^2\,\eta_{\mathrm{Hc}} + (1-\pi)^2\,\eta_{\mathrm{He}} - 2\pi(1-\pi)\,\sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}
\end{aligned}
$$

For the within part, the multipliers are constants given the transcriptome share, so the conditional variance weights the allele variances and covariance by the squared multipliers. The intrinsic and extrinsic decomposition of this appendix gives the allele variances as σ²_d + σ²_q and the covariance as σ²_d.

$$
\mathrm{Var}[a^*_{\mathrm{H}i} \mid c_i] = \left[\pi^2\,t_{\mathrm{Hc}i}^2 + (1-\pi)^2\,t_{\mathrm{He}i}^2\right]\left(\sigma^2_d + \sigma^2_q\right) + 2\pi(1-\pi)\,t_{\mathrm{Hc}i}\,t_{\mathrm{He}i}\,\sigma^2_d
$$

Averaging over cells needs the second moments of the multipliers, given in *Variance estimates*:

$$
\mathrm{E}[t_{\mathrm{Hc}i}^2] = 1 + \eta_{\mathrm{Hc}}, \qquad \mathrm{E}[t_{\mathrm{He}i}^2] = 1 + \eta_{\mathrm{He}}, \qquad \mathrm{E}[t_{\mathrm{Hc}i}\,t_{\mathrm{He}i}] = 1 - \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}
$$

Substituting the second moments gives the averaged conditional variance.

$$
\mathrm{E}\left[\mathrm{Var}[a^*_{\mathrm{H}i} \mid c_i]\right] = \left[\pi^2\left(1 + \eta_{\mathrm{Hc}}\right) + (1-\pi)^2\left(1 + \eta_{\mathrm{He}}\right)\right]\left(\sigma^2_d + \sigma^2_q\right) + 2\pi(1-\pi)\left(1 - \sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\right)\sigma^2_d
$$

The terms without η or √(η_Hc η_He) reduce to φ_Σ, because π² + (1 − π)² + 2π(1 − π) equals one. Adding the between part and grouping the remaining terms gives the variance of the summed expression factor:

$$
\mathrm{Var}[a^*_{\mathrm{H}i}] = \underbrace{\sigma^2_d + \left[\pi^2 + (1-\pi)^2\right]\sigma^2_q}_{\varphi_{\Sigma}} + \left(1 + \sigma^2_d\right)\left[\pi^2\eta_{\mathrm{Hc}} + (1-\pi)^2\eta_{\mathrm{He}} - 2\pi(1-\pi)\,\sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\right] + \sigma^2_q\left[\pi^2\eta_{\mathrm{Hc}} + (1-\pi)^2\eta_{\mathrm{He}}\right]
$$

The transcriptome share correction is never negative, because its first bracket is the variance of the conditional mean and its second is a weighted sum of the η's. It is largest when one allele dominates, where it equals η_Hc(1 + φ_H) or η_He(1 + φ_H) for the dominant allele. The NB fit of the summed counts, which we write φ̂*_Σ, estimates this variance and so carries the transcriptome share correction, which we subtracted. A simulation of the summed expression factor with shifted exponential deviations and 8 × 10⁶ cells matched this variance to within 0.0002 in three settings.

$$
\hat{\varphi}^*_{\Sigma} \approx \mathrm{Var}[a^*_{\mathrm{H}i}] = \varphi_{\Sigma} + \left(1 + \sigma^2_d\right)\left[\pi^2\eta_{\mathrm{Hc}} + (1-\pi)^2\eta_{\mathrm{He}} - 2\pi(1-\pi)\,\sqrt{\eta_{\mathrm{Hc}}\,\eta_{\mathrm{He}}}\right] + \sigma^2_q\left[\pi^2\eta_{\mathrm{Hc}} + (1-\pi)^2\eta_{\mathrm{He}}\right]
$$

We subtracted the estimate of this term, built from the corrected variance estimates, the estimated weight π̂ and the estimates η̂_Hc, η̂_He and √(η̂_Hc η̂_He) of the multiplier moments, and set a corrected dispersion below zero to zero.

$$
\hat{\varphi}_{\Sigma} = \hat{\varphi}^*_{\Sigma} - \left(1 + \hat{\sigma}^2_d\right)\left[\hat{\pi}^2\hat{\eta}_{\mathrm{Hc}} + (1-\hat{\pi})^2\hat{\eta}_{\mathrm{He}} - 2\hat{\pi}(1-\hat{\pi})\,\sqrt{\hat{\eta}_{\mathrm{Hc}}\,\hat{\eta}_{\mathrm{He}}}\right] - \hat{\sigma}^2_q\left[\hat{\pi}^2\hat{\eta}_{\mathrm{Hc}} + (1-\hat{\pi})^2\hat{\eta}_{\mathrm{He}}\right]
$$

## References

Cox DR, Reid N (1987) Parameter orthogonality and approximate conditional inference. Journal of the Royal Statistical Society, Series B 49:1–39. https://doi.org/10.1111/j.2517-6161.1987.tb01422.x

Elowitz MB, Levine AJ, Siggia ED, Swain PS (2002) Stochastic gene expression in a single cell. Science 297:1183–1186. https://doi.org/10.1126/science.1070919

Fraser HB (2019) Improving estimates of compensatory *cis*-*trans* regulatory divergence. Trends in Genetics 35:3–5. https://doi.org/10.1016/j.tig.2018.09.003

Lawless JF (1987) Negative binomial and mixed Poisson regression. Canadian Journal of Statistics 15:209–225. https://doi.org/10.2307/3314912

McCarthy DJ, Chen Y, Smyth GK (2012) Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40:4288–4297. https://doi.org/10.1093/nar/gks042

McCullagh P, Nelder JA (1989) Generalized Linear Models, 2nd edn. Chapman and Hall, London. ISBN 0-412-31760-5. https://doi.org/10.1007/978-1-4899-3242-6

Phipson B, Smyth GK (2010) Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. Statistical Applications in Genetics and Molecular Biology 9:39. https://doi.org/10.2202/1544-6115.1585

Raj A, Peskin CS, Tranchina D, Vargas DY, Tyagi S (2006) Stochastic mRNA synthesis in mammalian cells. PLoS Biology 4:e309. https://doi.org/10.1371/journal.pbio.0040309

Wedderburn RWM (1974) Quasi-likelihood functions, generalized linear models, and the Gauss-Newton method. Biometrika 61:439–447. https://doi.org/10.1093/biomet/61.3.439

Wittes J, Brittain E (1990) The role of internal pilot studies in increasing the efficiency of clinical trials. Statistics in Medicine 9:65–72. https://doi.org/10.1002/sim.4780090113

Wittkopp PJ, Haerum BK, Clark AG (2004) Evolutionary changes in *cis* and* trans* gene regulation. Nature 430:85–88. https://doi.org/10.1038/nature02698

## Notes for revision

### Decisions

- **Informativeness thresholds.** The 5% test level and the near-zero threshold φM̄_r < 0.1 in Informative estimates are provisional. Report how class counts and results change across nearby thresholds.
- **Transcriptome share variance.** Whether to estimate σ²_c within replicates.
- **ML fit of non-gamma expression factors.** The NB fit converges to the variance of the expression factor only when that factor is gamma, and the factor seen by the allele fits, t·a, is not exactly gamma. A more systematic approach is needed to quantify the departure across the means, dispersions and allele weights in the data, and to decide whether the allele dispersions should use the moment estimator of Appendix A instead.
- **Coexpression.** Whether the coexpression results enter the paper. If they do, add the coexpression ploidy correction as a short paragraph in Per-genome hybrid estimates.

### Checks to run

- **Code.** In .fit_one(), replace the Pearson pre-check with the Cox–Reid adjusted fit and its slope rule, setting φ̂ = 0 when Σ(y − M̂_i)² ≤ ΣM̂_i − ΣM̂_i²/ΣM̂_i, and remove PHI.FLOOR. Keep the truncation at zero after the transcriptome share correction, and record fits whose maximum reaches either end of LOG.THETA.RANGE, which spans φ from about 3 × 10⁻⁷ to 55.
- **Code.** Classify every dispersion estimate as resolved above zero, resolved near zero or uninformative, using the boundary likelihood ratio test and the upper bound on φM̄_r in Informative estimates, and gate each raw and log₂ contrast and each kinetic balance contrast on the classes of the groups it uses. Classify after the transcriptome share correction. Report excluded gene counts per contrast, keep a gene set for each contrast and scale for the supplement, and flag a single matched set, class 1 in every group of every partition, for every main-text comparison. In bootstrap replicates, map a zero estimate to −∞ on the log₂ scale and take log₂ intervals and SEs from percentiles.
- **Code.** Estimate γ and f_κ for kinetic balance, add a third hybrid partition at f_κ, and fit the mean and dispersion on the same cells of that partition for the kinetic balance contrasts, pooling γ over genes in class 1 on both sides. Report a supplementary version from the f_μ and f_φ partitions.
- **Analysis.** Report γ̂_x, f_x and the cis and trans partition sizes for x = μ, φ and κ, and confirm that each partition keeps about 800 cells or more.
- **Analysis.** Flag genes with more than 2.5% of bootstrap replicates at zero, and report whether any of them has a significant contrast.
- **Code.** Normalize each allele by the ratio of summed counts to summed offsets rather than the average of the cell rates.
- **Code.** Reinstate the allele transcriptome share correction and add the variance corrections of Appendix B. Truncate negative variance estimates only after the correction, not at the star estimators.
- **Code.** Refit the transcriptome share corrections inside each bootstrap replicate, and apply the observed c̄̂ and η̂_r inside every refit of the null distributions, permuted or simulated.
- **Code.** Replace the permutation nulls for the dispersion contrasts with the simulated nulls of Null distributions: null dispersions by projection on each scale, inflated allele dispersions for the transcriptome share, and a shared extrinsic factor with variance σ̂²_d for the two alleles of a cell. Build the parental mean null from rescaled hybrid cells, the trans mean null with allele offsets that include the abundance factor and k̂, and the midparent mean null by simulation with null allele means scaled to the midparent, without δ̂ or δ̂_syn.
- **Code.** Confirm that c̄̂ is the unweighted mean of ĉ_i.
- **Code.** Rescale the allele fits by 2c̄̂ and 2(1 − c̄̂) before averaging in the midparent and parental mean contrasts.
- **Code.** Compute the dominance dispersion contrasts from the corrected allele dispersions, as a geometric mean on the log₂ scale and an arithmetic mean on the raw scale. Compute δ̂ and δ̂_raw for reporting only.
- **Code.** Rename the contrasts dom and dpar to mid and par in .contrast_value() and the output tables.
- **Simulation.** Simulate known depths, means, dispersions, mean transcriptome share and transcriptome share variance through the full pipeline to test the remaining bias of the dominance contrasts.
- **Simulation.** Add the Oct 23, 2026 simulation scripts for the NB fit of the summed counts, the mean and dispersion correlation, and Appendix B to the repository, and rerun them at the study's parameter estimates.
- **Simulation.** Check the calibration of the dispersion nulls. Simulate genes with known equal dispersions and different means, and confirm that the simulated nulls give uniform p-values on both scales, including the cis contrast with its shared extrinsic factor.
- **Analysis.** Test the observed correlation of Δ(κ) and Δ(μ) across genes against its structural value, as done for burst size.
- **Analysis.** Support the statement that φ is estimated imprecisely with the median bootstrap SE of log φ.
- **Analysis.** Confirm the calibration of the level-adjusted *trans* mean null and the simulated midparent mean null in the null calibration simulation, including genes whose hybrid level differs from the parents'.
- **Analysis.** Compare the parental NB fits with and without a dropout term by AIC, and report whether dropout can be ignored for UMI data (suggested by Giersdorf et al. 2026).
- **Analysis.** Fit cell-state clusters directly, with cluster-specific means and a shared noise parameter, to estimate how much of φ comes from cell-state variation, which bears on σ²_d (suggested by Giersdorf et al. 2026).
- **Analysis.** Report the correlation of φ̂ with μ̂ across genes as a check on the offset-NB fits.
- **Analysis.** To do after the splits are rerun. The offset caps are computed (CV and cap 1 + CV²(O) of 0.72 and 1.52 for Pc, 0.87 and 1.76 for Pe, 0.44 and 1.19 for H). Compute each gene's realized variance ratio with var_ratio() on the MIX.SC, MIX.SE and hybrid fits, replace the illustrative CV of 0.53 in Estimating μ with the per-group caps and realized ratios, and decide whether the parental means should use the NB estimate μ̃ if their realized cost is large.
- **Analysis.** Confirm the choice of the allele average for φ_H over the summed hybrid route with the Cox–Reid fits on a subset of genes. A moment-based comparison on all genes (compare_phiG_routes.R) gave a Spearman correlation of 0.97 and a median bootstrap SD ratio of 0.965, summed route over allele route.

### Reference checks

- Consider Shahrezaei V, Swain PS (2008) Analytical distributions for stochastic gene expression. PNAS 105:17256–17261, https://doi.org/10.1073/pnas.0803850105, as the primary citation for the analytical NB result, alongside Raj et al. 2006.
- Consider citing Larsson AJM et al. (2019) Genomic encoding of transcriptional burst kinetics. Nature 565:251–254, https://doi.org/10.1038/s41586-018-0836-1, where θ and ζ are introduced, as a precedent for allele-specific burst kinetics from single-cell data.
- DOIs for Lawless (1987), McCullagh and Nelder (1989) and Wittes and Brittain (1990) were added from publisher and Crossref records. For McCullagh and Nelder, the Springer DOI 10.1007/978-1-4899-3242-6 is the one Crossref uses in citing works; Routledge also lists the eBook as 10.1201/9780203753736.
- Add a full citation with DOI for Giersdorf et al. 2026 (bioRxiv) if its suggested analyses enter the paper.

