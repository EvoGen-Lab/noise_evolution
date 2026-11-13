# Section 7 revision note

Section 2.5 (`analysis/analysis.R`) now conditions every contrast on the derepression and stress program (Methods Appendix C draft,
`docs/methods_appendix_C_state_conditioning.md`). Several parts of Section 7 describe the same cell states with different definitions. This note lists what
to revise so that Section 7 describes the states that Section 2.5 conditions on, and states the two-culture design where it matters. Nothing below is
implemented yet.

## What changed in the premises

- The parents were grown in one tube and the hybrid in another. Every comparison of parental with hybrid cells in Section 7 (7.4 cluster consistency between
  hybrid alleles and parents, 7.8 `COMP.BOUND.TRANS`) compares cultures as well as genotypes. Comparisons within one culture (Pc against Pe, the two hybrid
  alleles, the clustering of each dataset) are unaffected.
- Check ST1 (`analysis/checks/state_groups.R`, depth-aware scores) found:
  - cell-cycle phase is the only state axis that explains clear transcriptome structure outside its marker genes (R² with PCs up to 0.52 in Pe and 0.58 in the
    hybrid), with phase proportions similar in all three groups;
  - derepression and stress are one program (score correlation 0.56; 10 of 12 derepressed Pc cells and 10 of 10 derepressed Pe cells are also stressed);
  - in Pc this program is the first PC (13.6% of variance, Spearman 0.83 with the stress score), a gradient across the culture with the extreme cells at its end;
  - the stressed state marks 4.0% of Pc, 5.9% of Pe and 0.25% of hybrid cells (z > 2.5);
  - PC1 of Pe and of the hybrid tracks the number of detected genes (Spearman 0.67), not a state.

## Proposed revisions

1. **7.6 module scoring.** Replace the NonFerm module (GO:0006097 and GO:0006635) with the 20-gene program of Section 2.5 (`STATE.MEMBERS`), so one definition
   serves the scoring, the conditioning and the figures. Revise the cell-cycle gene lists: SWI5 and ACE2 are transcribed in the CLB2 (G2/M) wave, and the M/G1
   list of three genes (SIC1, SWI5, ACE2) is too short to score; ST1 uses SIC1, PCL9, ASH1, CDC6, EGT2, CTS1, DSE1, DSE2, SCW11 and AMN1. Score with the
   depth-aware z of `module_z()` rather than `AddModuleScore()`, so groups sequenced to different depths share one scale.
2. **7.7 covariate-noise diagnostic.** Keep the cell-cycle axis. Report the program as the covariate that Section 2.5 conditions on, with its frequency per
   culture, and point to `STATE.FLAG.SUMMARY` for its effect on the contrasts instead of repeating a separate noise association.
3. **7.8 composition comparison.** State that the parent against hybrid comparisons (`COMP.BOUND.TRANS`) compare cultures. Report the program frequency
   difference between the parental and hybrid cultures as a culture-confounded phenotype, not as a hybrid effect.
4. **7.3 to 7.5 clustering.** Cross-tabulate the Seurat clusters of each dataset against the ST1 labels (cycle phase and program) to name the clusters by
   state. In Pc, check whether a cluster follows the program gradient of PC1. Add the two-culture caveat to the 7.4 consistency of hybrid allele clusters with
   parental clusters. Compare the 7.5 between-cluster share of variance with the change in dispersion under conditioning (`COND.CONTRASTS`, `d disp`), since
   both measure the variance that cell states explain.
5. **Methods main text.** Add one paragraph to the main Methods stating the two-culture design and which contrasts it affects, with the pointer sentence to
   Appendix C.

## Reruns that the revision would need

Section 7 only (local), and `cluster_stability.R` only if 7.3 changes its inputs. Checkpoint 7 would gain the ST1-style labels if 7.6 adopts them.
