# Notes for the Discussion

## The hybrid's genome-wide rise in dispersion (checks MID1 to MID4, Nov 20)

Each hybrid allele is noisier than its parent across the genome. At the same sequencing depth and the same number of cells, the median log2 dispersion is about 0.8 higher for the *S. cerevisiae* allele than for the *S. cerevisiae* parent (about 1.75 times) and about 1.4 higher for the *S. eubayanus* allele than for the *S. eubayanus* parent (about 2.6 times). The two hybrid alleles end at the same level (median log2 phi about -2.5 each), although the parents differ: at matched depth the *S. eubayanus* parent is the less noisy (about -4.05 against -3.32). So the rise is larger for the *S. eubayanus* allele because its parent starts lower, not because *S. eubayanus* expression is noisier in general. That the two alleles converge on one level points to a noise source that both alleles in a cell share.

What the checks exclude:
- A counting artefact in the hybrid: the variance term linear in the mean is smaller in the hybrid alleles than in the parents (MID2).
- Sequencing depth: thinning each parent to its allele's depth removes about 0.2 of the 1.4 log2 shift (MID3).
- Cell number: subsampling the hybrid to each parent's cell count changes the zero-fit shares but not the shift (MID4).
- The summing of two allele copies: parents are haploid and each hybrid allele is one copy, no contrast sums the alleles, and summing independent copies would lower the dispersion rather than raise it.

What the design cannot separate:
- Separate cultures. The parents grew in one tube and the hybrid in another, so any difference between the cultures acts on every hybrid gene.
- Ploidy. The hybrid is diploid and the parents haploid, so cell size, cell cycle and other properties of diploid cells differ between the groups.
- Hybridity. A hybrid cell runs two regulatory systems at once, and mismatches between them may make expression generally less precise.

Consequences for the results:
- The *cis* contrast is unaffected, because a rise shared by both alleles cancels in their ratio.
- The *trans* contrast is touched only through the parental side.
- The dominance contrasts of dispersion (midparent and parental) are mostly this shared rise. Calls on single genes should be read as departures from it.
