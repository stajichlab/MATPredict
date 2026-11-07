# Idiomorph overlap resolver: bitscore rank basis (#275)

Status: opt-in and experimental, default unchanged. Written 2026-11-06.

## The option
`matpredict detect --idiomorph-rank-basis bitscore` makes the pre-polish overlap resolver (`resolve_idiomorph_overlaps`) compare the best
bitscore of each of two mutually exclusive idiomorph genes (MAT1-1-3 and MAT1-2-1, MATALPHA2 and MATA2, sexM and sexP) instead of the
identity of each gene's best hit (`identity`, the v0.8.x rule). The proteome-first standing is kept. A resolution made this way records
`basis: first_pass_bitscore`. It is in no preset.

## What it fixes and what it does not
- Fixes: a short HMG-box hit with a higher identity beating a longer, stronger hit of the other gene. *Microdochium nivale*: a 34.6% MAT1-2-1
  hit of 59.3 bits lost to a 36.6% MAT1-1-3 hit of 41.2 bits.
- Does not fix the *T. viride* case in the assembly-rescue pilot (#275): there no overlap resolution is made. The wrong call comes from the
  later idiomorph vote, where an unpolished 47-bit MAT1-2-1 fragment outscores a polished 39.7-bit MAT1-1-3 model. Same call under both bases.

## Measured (frozen v0.8.1 code, `--preset v0.8.1`, against the v0.8.1 reports)
| Panel | Identical | Coordinates or confidence only | Family, idiomorph or contig differs |
|---|---|---|---|
| Ascomycota, 300 stratified genomes | 255 | 32 | 13 |
| Basidiomycota, 120 | 120 | 0 | 0 |
| Mucoromycotina, 90 | 90 | 0 | 0 |

The 13 Ascomycota differences: 5 idiomorph flips in the MAT family (four MAT1-1 to MAT1-2, one MAT1-2 to MAT1-1), 2 genomes that gain an
undetermined MATsc locus, 6 yeast genomes (MATsc) that lose a duplicate or undetermined locus. Not yet examined: the 8 MATsc changes and the
32 coordinate or confidence changes.

The 5 flips against an independent `miniprot` model of the family's reference proteins on the locus plus 15 kb either side:

| Genome | v0.8.1 | Bitscore basis | Core models | Reading |
|---|---|---|---|---|
| GCA_964254555.1 (metabat2 MAG) | MAT1-1 low | MAT1-2 high | MAT1-2-1 (29% id, 188 aa) | new call supported |
| GCA_003073195.1 RamPxa02 | MAT1-1 medium | MAT1-2 high | MAT1-2-1 (26% id, 225 aa) | new call supported |
| GCA_022578595.1 Hyprub1 | MAT1-2 low | MAT1-1 low | none | neither supported |
| GCA_977970205.2 *Microdochium nivale* | MAT1-1 low | MAT1-2 low | none | neither supported |
| GCA_030077755.1 | MAT1-1 medium | MAT1-2 high | none | neither supported |

So 2 of 5 flips have a model behind the new call, 3 swap one weak call for another, and three flips raise the confidence to `high`, one of
them (GCA_030077755.1) with no core model at all.

## Why this is not a default
- Four of five flips go toward MAT1-2; a bitscore grows with protein length and the MAT1-2-1 references are longer than the MAT1-1-3 ones, so
  comparing the bitscores of two different genes may favour the longer reference. Not tested.
- It raises confidence where the evidence is one weak fragment. A tier cap for a call that rests on a single unpolished fragment under about
  50 bits is a separate scoring decision (curator ruling needed) and is not in this change.

## Open
The 8 MATsc yeast changes and the 32 coordinate or confidence changes; a length-normalised comparison; the tier cap; the *T. viride* vote.
