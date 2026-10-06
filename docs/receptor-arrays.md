# Pheromone-receptor arrays and `receptor_array_support`

Report-only. Source: `analysis/2026-10-06_agaricomycetes-pr-arrays.md` (options 1 and 2).
Code: `src/MATPredict/detect/receptor_arrays.py`. Nothing here changes whether a locus is called,
its confidence, its verification label or any count.

## What is reported

For every family with a `pheromone_precursor_scan` (Basidiomycota `PR`):

1. **STE3-like loci**: the family's non-superseded receptor-gene search hits, merged per strand where
   they overlap or lie within 300 bp (intron-split HSPs, several references over one gene). Hits on
   opposite strands are never merged. A merged locus counts only when its hits cover at least 50% of one reference
   receptor (HSP lengths summed per reference, overlapping HSPs once) and it spans 8 kb or less, the
   study's miniprot filter; without it the raw tblastn fragments gave about ten times too many
   single-locus arrays (*Trametes versicolor*: 67 merged loci, 6 with 50% coverage).
2. **Arrays**: loci on one contig linked by single linkage, gap between neighbouring locus ends of
   50 kb or less (the study's rule; 80% of arrays are single loci at 50 kb, and the share barely moves
   from 10 kb to 200 kb). An array never crosses a contig.
3. **`receptor_array_support`**: `supported` when the array has 2 or more loci, or a pheromone-precursor homology
   hit (genes of class `pheromone_precursor`, other than the scan gene) within the scan window (10 kb)
   of a member, or 2 or more distinct strict-CAAX ORFs within that window of a member. Otherwise
   `unsupported`. Study chance rates: 2 or more CAAX ORFs 4% of arrays, one ORF 18%; size-2 arrays are
   weak evidence (about half are expected from contig fragmentation).

## Fields

On each PR call in `detection_report.yaml` (absent on other calls) and as columns of a per-locus
`loci.tsv` (empty for other calls; `receptor_array_members` and `receptor_array_support_reasons` are `|`-joined):

| field | meaning |
|---|---|
| `receptor_array_id` | `<phylum>:<locus>:<contig>:<start>-<end>` of the array; same id for every call in it |
| `receptor_array_size` | number of STE3-like loci in the array |
| `receptor_array_members` | `start-end:strand` of each locus |
| `receptor_array_support` | `supported` or `unsupported` |
| `receptor_array_support_reasons` | met criteria (`array_size>=2`, `precursor_homology`, `caax_orfs>=2`), or for `unsupported`: `single_locus`, `no_precursor_homology`, `strict_caax_orfs=N(<2)` |

A PR call that overlaps no array reports nulls. Top-level `receptor_arrays` lists every array of the
genome once (including arrays with no call), with `strict_caax_orfs`, `precursor_homology_hits` and
`calls`; its per-array keys use the same `receptor_array_*` names as the call fields; `receptor_arrays_note` carries the caveat below.

## Scope: PR family only, and the B locus

In tetrapolar Agaricomycetes the pheromone receptor genes are part of the **B mating-type locus**
(Casselton and Olesnicky 1998, *Microbiol Mol Biol Rev* 62:55). "PR" in this pipeline is the
pheromone-receptor family and may be one part of a B locus: in *Schizophyllum commune* the B locus has
two subloci (B-alpha and B-beta), each with one receptor and several pheromone genes; *Coprinopsis
cinerea* has a B locus of about 17 kb with three receptor-plus-two-pheromone cassettes.

Arrays are reported for the PR family only. B-alpha/B-beta receptors form arrays only when they are merged
into a PR call. How B-locus sublocus structure clusters is a planned separate assessment (analysis), not
done here.

## Parameters

Accepted and fixed: merge gap 300 bp; a locus needs at least 50% reference coverage and a span of at most
8 kb; precursor-homology and strict-CAAX window 10 kb around a member; array gap 50 kb. Each call
reports one array, the one with the most loci overlapping the call (ties: leftmost); an array is
identified by coordinates (`receptor_array_id` = `<phylum>:<locus>:<contig>:<start>-<end>`).

## Reading it

Arrays occur for both mating and non-mating receptors. In *Coprinopsis cinerea* and *Schizophyllum
commune* the mating receptors lie in an array that also holds paralogs, so **array membership does not
establish that a locus is a mating receptor**, and `receptor_array_support` says how much evidence sits behind
the array, not which copy mates. It is not used to lift or drop an `unverified` label.

The loci come from this run's reference-protein search, not from the study's miniprot scan of 1,094
STE3 queries, so array sizes can differ from the study's.
