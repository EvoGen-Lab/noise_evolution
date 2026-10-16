# Yeast transcriptional noise pipeline

Reproducible R pipeline for genome-wide noise evolution in S. cerevisiae (YPS1000), S. eubayanus, and their F1 hybrid (10x scRNA-seq, allele-specific). Scientific context and settled decisions live in @docs/project_context.md. Read it before changing any analysis logic.

## Files
- `analysis/analysis.R` main analysis, organized in numbered sections (1, 2.1, 2.3, 4, 6.1 ... 8).
- `R/setup.R` session setup shared by `analysis.R`, `power_analysis.R` and `external_validation.R`: renv check, project paths and output folders, class names, color palettes, and `R/functions.R`. Each script starts with `library(here); source(here("R", "setup.R"))`. Defines no functions.
- `R/functions.R` helper functions sourced by the main script.
- `analysis/power_analysis.R` Section 10 (power analysis) with its functions in `R/functions_power.R`, also sourced by `power_grid.R`. Sections 10.1 to 10.8 simulate the total contrast; 10.9 to 10.14 simulate cis and trans on paired-allele data (`power_modes_grid.R`, `power_modes.sub`). Reads `PR`, `N.PERM` and the real design (hybrid and parent cell counts, `SPLIT.FRAC`) from the Section 2 checkpoint.
- `analysis/external_validation.R` Section 9 (comparison with published data) with its functions in `R/functions_external.R`. Reads `GENES` and `CONTRAST.FITS` from the Section 2 checkpoint.
- Cluster scripts, each a standalone R script plus a `.sub` SLURM file

## Environment
- Local: RStudio on Windows. Cluster: Negishi (RCAC), account `bphm`, partition `cpu`, module `r/4.5.2`.
- Both environments are in R 4.5.x with Seurat v5. Use one shared dated Posit Package Manager CRAN snapshot for both.
- Section 2 of `analysis/analysis.R` computes offsets as direct library-size factors (no scran): a parent's offset is its library over the median parental library, each hybrid allele is offset by half the cell's hybrid depth, and the summed hybrid by the full depth. Cell-to-cell variation of the Sc share is not an offset; `.fit_one(share_r)` removes its contribution from each hybrid allele's dispersion (`SHARE.R`, Section 2.0). `analysis/sim_inputs.R` and the estimator simulations still model the earlier scran haplotype offsets with `scran::computeSumFactors` (they need scran installed) and have not been updated to the library-size offsets.
- Claude Code cannot reach the cluster. Write the scripts and `.sub` files, never run `sbatch` or assume cluster output exists.

## Git workflow
- `main` holds the reviewed, runnable pipeline. Never commit or push to `main` directly.
- `claude/refactor` is the long-lived branch for large refactors. Commit and push here.
- Small fixes arrive on `feature/*` or `fix/*` branches. Merge them into `claude/refactor` when asked, resolve conflicts in favor of preserving both intents, and report each conflict resolved.
- Use small commits with messages that state the behavior change and the sections affected.
- Never force-push and never rewrite published history.
- Open a pull request from `claude/refactor` to `main` for review. Do not merge it.

## Rerun protocol (required in every change summary)
Cluster-dependent sections are 2.1 (pilot), 2.3 (bootstrap and permutation), 4 (coexpression), 6.1 (NuPoP), and 10 (power analysis).
Two further cluster jobs source `R/functions.R` and are resubmitted when the functions they call change: `cluster_stability.R` (7.3 resolution stability and the marker enrichment used in 7.4) and `go_enrich.R` (8.4 GO and KEGG enrichment).
Every summary of a change to `analysis/analysis.R` or `R/functions.R` must list
1. which sections and subsections need a rerun,
2. which cluster jobs need resubmission.
Prefer fixes that avoid cluster reruns. Offer a local recompute path (for example by editing checkpoints) whenever one exists.

## Architecture rules
- Cluster pattern: the main analysis packages inputs locally, saves an `.rda`, the cluster script computes, and the main analysis loads the output and verifies that keys match. Keep new cluster work in this pattern.
- NuPoP runs only on the Linux cluster. Its compiled Fortran aborts the hosting R process on Windows regardless of input.
- Bootstrap resampling draws all indices before any Seurat call (`make_boot_ky()`). `RunPCA()` reseeds the global RNG, so a single `set.seed()` before a loop collapses the draws. `collapse_score()` should stay near 0 (0.01 confirmed).
- Fix architectures, not individual crash triggers.
- Checkpoints follow `sectionN_checkpoint.rda`. If a refactor adds an object that downstream sections need, add it to the checkpoint

## Conventions
- R files use Windows CRLF line endings. `.sub` SLURM files use Unix LF line endings. `.gitattributes` enforces both. Keep it intact.
- Each section of an analysis script starts with `set.seed(SEED.SECTION + n)` (`SEED.SECTION` is defined in `R/setup.R`). A step that needs its own reproducible stream draws it through `with_local_seed(seed, expr)`, which restores the caller's random state.
- Check scalar conditions. Use `&&` and `||` only with length-one logical values. A parse-based scan of every `&&` and `||` in `R/`, `analysis/` and `cluster/scripts/` found only scalar operands.
- Write a function only when its code runs more than once: it is called from two or more places, or it is applied across many items (`lapply`, `sapply`, `vapply`, `Map`, `mapply`, `apply`, `parLapply`, `mclapply`). Code that runs once is written inline where it runs, so a reader follows the steps in order.
  - `analysis/analysis.R`, `analysis/power_analysis.R`, `analysis/external_validation.R`, `R/setup.R` and the cluster scripts define no named functions. A function that is applied across items, or called from two or more places, lives in `R/functions.R`, and is passed to the apply call with its other inputs as named arguments (`lapply(X, name, arg = value)`). A short one-line anonymous function in an apply call (`sapply(x, function(g) cor(a[g], b[g]))`) is fine inline and should not become a separate function. A multi-line body belongs in a named function in `R/functions.R`. Use `[[` for plain accessors.
  - A once-used step in `analysis.R` or a cluster script is written inline (a `local({ ... })` block when it needs private variables), with its concept comment above it.
  - A function that PSOCK workers need by name must be a top-level function in `R/functions.R` (`neg_binom_fit_offset`, `.fit_one`); `fit_counts_offset()` finds such dependencies from the code.
  - Before adding a function, check that it will run more than once; inline it if not.
- Noise axis: `DISP` is the dispersion phi (Var = mu + phi mu^2); theta = 1 / phi is the NB size and behaves as burst frequency, and `.fit_one()` returns both. Every noise contrast is reported on two scales from the same refits: a log2 contrast of phi (the `disp_*` columns of `BURST.CONTRASTS` and `PR`) and a raw difference of phi (the `rdisp_*` columns; defined when phi = 0). `fit_dispersion(method = "NB" or "weighted")` is the shared fitting interface (`.fit_one()` is its NB wrapper), and every dispersion contrast (total, cis, trans, dom, dpar) is reported for all four combinations of scale and estimator: log2 x NB (`disp_*`, the primary), raw x NB (`rdisp_*`), log2 x weighted method of moments (`wdisp_*`) and raw x weighted (`wrdisp_*`), each with a bootstrap estimate and SE in `BURST.CONTRASTS` and a permutation p-value and BH q-value in `PR` (the dpar ones also `_p_ploidy` and `_p_ind`). Means are log2 x NB only. Burst size and kinetic balance are log2 x NB (`bsize_*`, `kbal_*`) and log2 x weighted (`wbsize_*`, `wkbal_*`), with the same estimates, SEs, permutation p-values, ploidy-adjusted dpar p-values and resampling trans null; raw-scale versions are not defined yet. Bootstraps give confidence intervals and standard errors only; every p-value comes from the permutation job (`gene_perm`). Its trans test of each noise quantity is a resampling null (`permute_contrasts_one(resamples = NULL.DRAWS)`: the resample distribution of the trans contrast recentered on its observed value), because it tests that the trans contrast is zero on its own scale and works on the raw and log2 scales; the pooled-relabeling p-values are kept as `*_trans_p_pooled`. Burst size is b = mu phi, so `bsize` = mean + disp, and kinetic balance is kappa = theta / b, so `kbal` = -mean - 2 disp. phi = 0 (Pearson chi-square / df <= 1, or phi below `PHI.FLOOR`) is the Poisson limit and is excluded from every log2 contrast and from the current gene filter (`GENES`); `GENE.TABLE` keeps every gene with its fit status and per-comparison eligibility flags.
- Select Seurat columns by name, never by position (Seurat v5 changes `FindMarkers` output order).
- Validate numeric inputs at function entry (for example `BSIZE` in `SPLIT.FITS`).
- Do not commit data, `.rda` checkpoints, or cluster output. Keep them in `.gitignore`.

## Code annotation style
- Annotate each function and each analysis block with comments that explain the concept behind the step, such as why bootstrap draws must be independent.
- Write comments positively. State what the code does and why it works, not what used to break.
- Keep comments short and accurate. Update them whenever the code changes.

## Verification before finishing a task
- Run `Rscript -e "parse('R/functions.R'); parse('analysis/analysis.R')"` to confirm both files parse.
- Source `R/functions.R` and run any unit checks that exist for the functions you touched.
- Do not claim a result is unchanged unless a check confirms it. Expected result changes after the R and Seurat upgrade include `FindMarkers` fold changes, newer GO and KEGG annotations, and Monte Carlo seed behavior.

## Open work
1. Rerun the full pipeline after the R upgrade and review the expected result changes. Compare the cluster package versions with `renv.lock` through `check_pkg_versions()` (run by `R/setup.R` once the cluster scripts' `pkg_versions_*.csv` files are copied back).
2. Address the trans versus cis power asymmetry for noise (see @docs/project_context.md).
3. Add a filter for likely multiplets (Section 1 QC, before the offsets): UMI-count and gene-count outliers, a doublet score, and in any run that mixes the two species, cells with substantial reads from both genomes. Multiplets inflate library size and dispersion and bias hybrid-versus-parent noise contrasts. Every section reruns once it is in.
