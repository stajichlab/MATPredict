# Companion-gene veto for the partial-coverage rule (2026-10-18)

Branch `companion-gene-veto` from main (5b9a9e0, includes PR #95). Curator's rule, 2026-10-18, from the flip review
(`origin/flip-review`: 46 flips, 41 right, 4 wrong: Monilinia vaccinii-corymbosi, Sclerotinia nivalis x2, S. sativa).

## The rule

When `--partial-core-policy coverage_profile` would flip a call away from idiomorph X, and a companion gene declared for X is intact in the
same locus (a modelled gene-evidence entry, not `lost_to`, coverage >= FULL = 90% of its reference), the call is `undetermined` instead of a flip, and
`idiomorph_vote_note` carries `decision: undetermined`, `veto: companion_intact`, `companion_gene`, `companion_coverage` and a reason ("companion gene
MAT1-1-5 intact ... the locus holds both idiomorphs"). The companion is read from the gene evidence, not the vote, so it works for `optional` genes
that `--idiomorph-policy core_first` keeps out of the vote (by construction; not run, because core_first is not on main).

Order of checks in `partial_core.evaluate`: trigger -> contig-edge guard (keeps the vote) -> companion veto -> profile confirmation. Thresholds,
the 25-bit margin (not changed, to be recalibrated), the default `off` and the narrow `medium` cap are unchanged. The both-strong veto was not implemented.

## Declaration

New optional locus key `idiomorph_companions: {idiomorph: [gene, ...]}` in `db/<Phylum>/order.yml`, in `db/_schema/order.schema.yaml` and as
`Family.idiomorph_companions`. A family with none behaves exactly as before. Declared so far: `Ascomycota:MAT` -> `{MAT1-1: [MAT1-1-5]}`. Nothing
for MAT1-2 (MAT1-2-10 was not reviewed) and nothing for other orders.

## Evaluation (flag on, with the veto vs without: the same code, `db/` with and without the key)

237 genomes x 2 variants, HPCC `cgv-run` array (preempt, 2 CPUs, 16 GB, 15 at once): the 205 truth-set runs (203 BioSample + T4 + B05.10), plus the
flip-review genomes not among them. Per-locus table of the 49 review loci: `results/2026-10-18_companion_gene_veto/flip_review_veto.tsv`.

| set | without veto | with veto |
|---|---|---|
| 49 flip-review loci: right flips (42 incl. unexamined) | 42 flipped | 42 flipped (none lost) |
| review "both idiomorphs present": S. nivalis x2, S. sativa | flipped to MAT1-2 | undetermined, `companion_intact` |
| review "both present": Monilinia vaccinii-corymbosi | MAT1-2 (no flip, see below) | same |
| Polytolypa, Calogaya (already undetermined); A. fumigatus edge-kept; S. spermophila ("unclear", MAT1-1-5 coverage blank) | unchanged | unchanged (spermophila still flips) |
| 34 Botrytis (correct/wrong) | 33 / 1 | 33 / 1, no change |
| 171 other truth genomes | 163 / 4 / 0 / 4 | identical |
| all 237 | | exactly 3 genomes change label |

The right flips kept: Botrytis x31 (all MAT1-2), Botryotinia calthae, Stromatinia, Ophiostoma x5, Leptographium, Intubia, Sporothrix. Gains and losses against
truth: 0 / 0 (none of the 3 changed genomes is in the truth set).

## Findings and caveats

- **Monilinia no longer flips on current main**, independent of the veto: after PR #94 (lost-to by location) its MAT1-1-1 is modelled as a 58% piece
  rather than a 97.7% span, so the trigger (other idiomorph full length) does not hold and the call is MAT1-2 (high), as the vote says. The veto was
  therefore NOT exercised on it, although MAT1-1-5 is intact there (97.6%) and would veto. The review verdict "both idiomorphs present" says MAT1-2 (high)
  is not clean either; that is not this change.
- All four wrong flips are one clade (Sclerotiniaceae, Leotiomycetes), and MAT1-1-5 is Leotiomycete-specific. The veto fires on only 3 loci here. It is a
  clade-limited fix from four loci; other fused/homothallic clades with different companions are not covered.
- The veto reads a companion's coverage as the same length ratio as the rule (model aa / reference aa, capped at 100). A reference that is a fragment
  would read as intact.
- Sclerotinia spermophila (review: unclear) still flips, its MAT1-1-5 is not intact.
- `core_first` composition not run (flag not on main); only the argument above. 
- Tests: `tests/detect/test_companion_gene_veto.py` (11: intact -> undetermined with the note; absent, partial or unknown companion -> flip; nothing declared, or
  declared for the other idiomorph -> unchanged; edge guard and non-trigger unaffected; Botrytis pattern still flips; roster key is read; the Ascomycota
  roster declares only MAT1-1-5), plus the existing partial-core tests. tests/detect + tests/db on the HPCC `test` env: 1116 passed, 4 skipped
  (test_draw/test_synteny excluded: no pygenomeviz; the 4 skips need `--idiomorph-policy`).
