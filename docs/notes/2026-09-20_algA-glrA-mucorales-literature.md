# algA/algL and glrA at the Mucorales MAT locus — literature findings

Date: 2026-09-20. Prompted by a curation question: `algA` and `glrA` sat in
`db/Mucoromycota/order.yml`'s gene roster with **no reference protein anywhere
in `db/`**, so no genome could ever match them, and the roster spelling (`algA`)
disagreed with the spelling in the review that describes them (`algL`).

Sources were checked directly (PubMed/PMC full text, GenBank flat files via
E-utilities, Crossref). Claims below are marked with what actually supports
them. Several targets could not be retrieved; those are stated as not found
rather than guessed.

## Headline: algA and algL are one gene

Established by sequence and by the authors' own inconsistency, **not** by
naming authority:

- `AFA26122.1` (*Mucor mucedo* "putative alginate lyase", Wetzel et al. 2012)
  and `WZC35661.1` (AlgA, OR965930) are **both 384 aa and align at 62.8%
  identity** over 387 alignment positions. Unambiguous orthologs.
- Domain annotation agrees: `WZC35661` carries **Alginate_lyase pfam05426**
  (CDD:398861), matching the review's "predicted alginate lyase".
- **Idnurm (2025) writes `algL` in the paper text while that same paper's own
  GenBank deposits use `/gene="algA"`.** One author, one publication, two
  spellings.

NCBI confirms the split is purely prose-vs-deposit: `algL AND Mucorales[orgn]`
returns **0** nucleotide records; `algA AND Mucorales[orgn]` returns exactly the
three MAT-locus deposits.

**Decision for this project: use `algA`.** It is the spelling attached to
deposited sequence, and therefore the one a database keyed to accessions can
actually resolve. No rename, correction or erratum was found for either
spelling. A third spelling, `glrR`, appears in an unverified search snippet of
Schulz et al. 2016 — noted only so a future reader is not surprised by it.

## Who actually described these genes

**The Lee & Idnurm 2017 review is not a source.** Its sentence — "These include
the *algL* and *glrA* genes, which encode a predicted alginate lyase and
glutathione reductase, respectively" — **carries no citation marker**, and
neither does the sentence before it. Tracing "the reference the review cites"
was therefore a dead end; there is none.

| gene | first primary description | what it is called there |
|---|---|---|
| `glrA` | **Idnurm A (2011)** *Eukaryot Cell* 10(11):1485-1491, doi:[10.1128/EC.05149-11](https://doi.org/10.1128/EC.05149-11), PMC3209044 | glutathione **oxidoreductase** |
| alginate lyase | **Wetzel J, Burmester A, Kolbe M, Wöstemeyer J (2012)** *Microbiology* 158(4):1016-1023, doi:[10.1099/mic.0.054106-0](https://doi.org/10.1099/mic.0.054106-0) — deposits `JN587498`/`JN587499` | `/product="putative alginate lyase"`, **no gene symbol at all** |
| the symbol `algL` | **Schulz et al. (2016)** *Endocytobiosis and Cell Research* 27(4):39-57 | per Idnurm 2025's own citation — **could not be read** |

Note the product wording drift for `glrA`: "glutathione **oxidoreductase**"
(Idnurm 2011, and the GenBank records) vs "glutathione **reductase**" (2017
review). Same enzyme family; `WZC35665` carries `Pyr_redox_2` (CDD:476868),
which is the pyridine nucleotide-disulphide oxidoreductase family glutathione
reductase belongs to.

**Not retrievable:** Schulz et al. 2016 has no DOI and is absent from Crossref,
PubMed, PMC and Semantic Scholar; the publisher archive and ResearchGate both
block automated access. It is the one document that would settle both the
symbol's origin and the conservation question below, and it likely needs an
institutional library request.

Checked and confirmed **silent** on algL/algA/glrA: Gryganskyi et al. 2010
(*PLoS One* 5:e15273), Lee et al. 2010 (*MMBR* 74:298-340), Li et al. 2011
(*PLoS Pathog* 7:e1002086), Zhang et al. 2017 (*Sci Rep*, *M. irregularis*).
Idnurm et al. 2008 (*Nature* 451:193-196) could not be retrieved; its deposits
`EU009461`/`EU009462` carry only tpt/sex/rnhA features, so it is unlikely to be
the source — but that is an inference from the deposits, not a reading of the
paper.

## Gene order

From the GenBank flat files directly. `OR965930` is the only record carrying
all five:

```
algA   -   complement(984..2430)
tptA   -   complement(4734..6201)
sexP   +   8845..9759
rnhA   +   9891..14072
glrA   -   complement(14435..16043)
```

So **algA — tptA — [sex] — rnhA — glrA**: `algA` lies *distal to tptA*, outside
the core on the tptA side; `glrA` lies *distal to rnhA*, on the rnhA side.
Corroborated for `glrA` in *Rhizopus azygosporus* (`MG967659`) and *Syzygites
megalocarpus* (`JN112239`/`JN112240`, one copy a pseudogene).

**One documented exception.** In *Mucor mucedo* (`JN587498`/`JN587499`) the
alginate lyase sits ~259 bp from the sex locus with **no tptA between them** —
though neither 10.5 kb record contains tptA at all, so it may simply lie
outside the sequenced fragment. Reported as observed; the records offer no
explanation.

**Caution.** The 2017 review says *sexM* is "convergently transcribed with the
flanking *rnhA* gene". That does not match the coordinates: in `FJ009106/7`,
`JN587498/9` and `OR965930` the sex gene and `rnhA` are on the **same** strand
in tandem, not head-to-head. Trust the coordinates.

## How conserved are they, really?

Weaker than the review's phrasing suggests. No systematic survey across
Mucorales was found.

- Lee & Idnurm 2017 says only "conserved between **some** species" — no taxa
  named, no count.
- Idnurm 2011, for `glrA`: adjacent to the locus in *R. delemar* and
  *M. circinelloides*, and present in *S. megalocarpus* (one intact copy, one
  pseudogene).
- Idnurm 2025, for `algL`: "commonly observed adjacent to mating type loci in
  Mucorales species", citing Schulz et al. 2016 — the unreadable one.

Direct GenBank evidence at a MAT locus, verified: **`algA` in 4 species / 3
genera** (*Mooraboolomyces wintlei*, *Absidia urquhartii* ×2 strains, *Mucor
mucedo*); **`glrA` in ~5 species** (*Mooraboolomyces*, *Syzygites*, *Rhizopus
azygosporus*, plus the textual claim for *R. delemar* and *M. circinelloides*).

**A handful of genera, not "Mucorales-wide".** This matters for scoring: `glrA`
in particular now has exactly **one** curated reference protein, so a low hit
rate for it in a sweep is expected and is not evidence the detector is failing.

## No MAT-locus deposits exist for the genera we most need

Queried against nuccore with sexP/sexM/"mating type locus"/"sex locus"/rnhA/
tptA/HMG-domain, **zero hits** for: *Cunninghamella*, *Chaetocladium*,
*Actinomucor*, *Lichtheimia*, *Circinella*, *Backusella*, *Thamnidium*,
*Benjaminiella*, *Pilobolus*, *Gongronella*, *Zygorhynchus*, *Apophysomyces*,
*Saksenaea*, *Rhizomucor*.

This independently confirms the curation decision of 2026-09-20: the only way
to add *Cunninghamella* or *Chaetocladium* references is assembly-derived
tier-2 records, which the database design spec excludes in this phase. There
was no tier-1 alternative to choose.

## Naming collisions to guard against

1. **`alyA` is NOT an alginate lyase.** In `JN112240` (*S. megalocarpus*),
   `/gene="alyA"` carries `/note="putative alpha-arrestin"`. A prefix match on
   `aly*` would produce a false positive.
2. **`sagA` is ambiguous.** Idnurm 2011 defines it as an unknown-function gene
   with pfam04082; the *M. mucedo* deposits separately annotate a "putative
   serine-rich adhesion protein" in a comparable position. Whether these are
   the same gene could not be determined — **do not merge them**.
3. **`arbA`/`ArbB`**: `JN112239` uses `/gene="arbA"` with `/product="ArbB"`, an
   internal inconsistency in that record.

## What was done with this

Three tier-1 records ingested (commit `b6e920f`): `OR965930.1`, `PP971768.1`,
`PP971769.1`. `algA` and `glrA` have curated reference proteins for the first
time, and `order.yml` carries the naming resolution as a comment so it is not
re-derived.

## Further published Mucorales MAT-locus accessions

Recorded here as curation leads; none are ingested yet.

| genus | accessions | source |
|---|---|---|
| *Phycomyces* | `EU009461`, `EU009462`; `LN554891`-`LN554895` (*P. nitens*) | Idnurm et al. 2008 *Nature* 451:193-196; Camino et al. 2015 *Fungal Biol* 119:1007-1021 |
| *Mucor* | `FJ009106/7`, `HM565940/1`, `HM754261/2`, `JN587498/9` (*M. mucedo*), `KX966017` (*M. moelleri*), `KY434081`-`KY434097` (*M. irregularis*) | various; Wetzel et al. 2012; Schulz et al. 2016 |
| *Rhizopus* | `HQ450311`-`HQ450316`, `HQ435186`-`HQ435239`; `MG967658`-`MG967660` | Gryganskyi et al. 2010 *PLoS One* 5:e15273; Gryganskyi et al. 2018 *G3* 8:2007-2018 |
| *Syzygites* | `JN112226`-`JN112233`, `JN112239`, `JN112240` | Idnurm 2011 |
| *Parasitella* | `KY081664` | Schulz et al. 2016 |
| *Blakeslea* | `HG939557`, `HG939558` | — |

Ingesting the *Rhizopus azygosporus* (`MG967659`) and *Syzygites* (`JN112240`)
records would add second and third `glrA` references, which is the thinnest
part of the roster.
