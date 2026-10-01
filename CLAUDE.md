# Yeast transcriptional noise pipeline

Reproducible R pipeline for genome-wide noise evolution in S. cerevisiae (YPS1000), S. eubayanus, and their F1 hybrid (10x scRNA-seq, allele-specific). Scientific context and settled decisions live in @docs/project_context.md. Read it before changing any analysis logic.

## Files
- `analysis/analysis.R` main analysis, organized in numbered sections (1, 2.1, 2.3, 4, 6.1, 10 ...).
- `R/functions.R` helper functions sourced by the main script.
- Cluster scripts, each a standalone R script plus a `.sub` SLURM file

## Environment
- Local: RStudio on Windows. Cluster: Negishi (RCAC), account `bphm`, partition `cpu`, module `r/4.5.2`.
- Both environments are in R 4.5.x with Seurat v5. Use one shared dated Posit Package Manager CRAN snapshot for both.
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
- Check scalar conditions. Use `&&` and `||` only with length-one logical values. A scan of the full codebase for violations is still unfinished.
- Write a function only when its code runs more than once: it is called from two or more places, or it is applied across many items (`lapply`, `sapply`, `vapply`, `Map`, `mapply`, `apply`, `parLapply`, `mclapply`). Code that runs once is written inline where it runs, so a reader follows the steps in order.
  - `analysis/analysis.R` and the cluster scripts define no named functions. A function that is applied across items, or called from two or more places, lives in `R/functions.R`, and is passed to the apply call with its other inputs as named arguments (`lapply(X, name, arg = value)`). A short one-line anonymous function in an apply call (`sapply(x, function(g) cor(a[g], b[g]))`) is fine inline and should not become a separate function. A multi-line body belongs in a named function in `R/functions.R`. Use `[[` for plain accessors.
  - A once-used step in `analysis.R` or a cluster script is written inline (a `local({ ... })` block when it needs private variables), with its concept comment above it.
  - A function that PSOCK workers need by name must be a top-level function in `R/functions.R` (`neg_binom_fit_offset`, `.fit_one`); `fit_counts_offset()` finds such dependencies from the code.
  - Before adding a function, check that it will run more than once; inline it if not.
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
1. Finish the `&&` and `||` scalar condition scan.
2. Finalize `pkg_versions()` and `check_pkg_versions()` to record and compare package versions across local and cluster.
3. Rerun the full pipeline after the R upgrade and review the expected result changes.
4. Address the trans versus cis power asymmetry for noise (see @docs/project_context.md).