# MAT locus visualization research — single-locus diagrams and group/synteny comparison

Research task, not an implementation plan. Grounds tool choice in MATPredict's
real current output formats before recommending anything.

## 1. What MATPredict already produces (the real starting point)

Read directly from source, 2026-09-18:

- **`src/MATPredict/db/gff_export.py`** (curated records): `write_gff3` emits
  flat `gene` features only (no CDS/mRNA/exon-level features) — 1-based
  fully-closed coordinates, `Name`/`role` attributes. `write_genbank` emits
  one `SeqRecord` per locus segment with the SAME flat `gene` features, BUT
  **the nucleotide sequence itself is a placeholder of all-`N` characters**
  (explicit in the docstring: "we do not have the actual nucleotide sequence
  for the segment... a placeholder until sub-project 2's tooling can fetch
  the real assembly sequence"). No CDS feature, no `translation` qualifier,
  no exon/intron structure anywhere in this writer.
- **`src/MATPredict/detect/report.py`**: `write_detection_gff3` emits one
  `MAT_locus` region feature per segment plus one `gene` feature per gene
  (found or explicitly `present=false`/`not_searchable=true`), same flat
  gene-feature model, Parent/child GFF3 relationships scoped per-contig to
  handle fragmented multi-contig loci. Also no CDS/exon features.
- **Domain fields worth visualizing** (`db/*/order.yml`, confirmed via
  `db/Ascomycota/order.yml`): `role` (`core_MAT`/`flanking_conserved`/
  `flanking_variable`, schema-enforced enum), `gene_class` (e.g.
  `alpha_box`/`HMG_box` — the two idiomorph-defining domain classes),
  `present_in_idiomorphs` (which idiomorph, e.g. `MAT1-1`/`MAT1-2`, a gene
  belongs to). A useful MAT diagram should color/label by at least `role`
  and ideally `gene_class`/idiomorph — none of this exists in the current
  GFF3/GenBank output today; it would need to be added as feature
  qualifiers/attributes for any downstream tool to use it directly.

**Consequence for every candidate tool below**: none of MATPredict's current
GenBank output carries real sequence or CDS/translation data. Any tool that
computes percent-identity links between homologous genes (clinker,
pyGenomeViz's alignment mode, MCScan-family tools) needs REAL protein or
nucleotide sequence to do that — this is not available today from
`write_genbank` as-is. It IS available elsewhere in the codebase (real
per-gene protein sequences already exist in each record's `proteins.faa`,
and — per this session's other work — `_independent_translation` can derive
one from curated coordinates when no accession exists), so the real gap is
in `gff_export.write_genbank`, not in the underlying data: it would need to
(a) fetch/embed the real segment nucleotide sequence instead of `N`s, and
(b) emit real `CDS` features (ideally with a `translation` qualifier from
the already-existing protein sequence) alongside the `gene` features.

## 2. cblaster + clinker

**cblaster** (Gilchrist & Chooi, *Bioinformatics* 2021; github.com/gamcil/cblaster,
995+ commits, GitHub Actions CI present): searches for CO-LOCATED homologous
gene clusters — you give it protein FASTA queries (or NCBI accessions), it
BLASTs them against a genome database (remote NCBI BLAST, or a local
DIAMOND database), applies identity/coverage/e-value filters, then reports
genomic regions where enough of your query proteins hit close together.
Installable via `pip install cblaster` (pip only per the fetched README;
did not confirm a bioconda recipe — see open question below).

**clinker** (Gilchrist & Chooi, same lab/paper family;
github.com/gamcil/clinker, 655 stars/80 forks, cited from 2020
*Bioinformatics*): takes a SET of GenBank or GFF3 files (GFF3 needs a
same-named FASTA alongside it) representing multiple gene clusters, aligns
them pairwise (BioPython's built-in aligner), and draws a `clustermap.js`
interactive SVG showing genes as arrows with percent-identity-colored links
between homologous genes across clusters. Supports custom coloring via a
`-cm/--colour_map` CSV (group→hex color) and custom function labels via
`-gf/--gene_functions` (gene ID→function name) — both are exactly the kind
of per-gene annotation MAT loci need (role/gene_class/idiomorph). Installable
via `pip install clinker` or `conda install -c conda-forge -c bioconda
clinker-py` — real bioconda availability, fits this project's existing
`pixi.toml` channel list (`conda-forge`, `bioconda`) directly.

**How they chain**: cblaster's own docs mention "both cblaster and clinker
can now be used without installation on the CAGECAT webserver" as a
complementary pairing, but the exact session-file handoff format between
cblaster's cluster-finding output and clinker's multi-GenBank input was not
confirmed from the fetched README content (see open questions). **Clinker
does NOT require cblaster** — it runs standalone on any folder of real
GenBank/GFF3+FASTA files, which is the more directly relevant path for
MATPredict, since this project already HAS curated loci and detection
results as discrete genomic regions; it doesn't need cblaster's own
cluster-DISCOVERY step (searching a genome database for co-located hits) —
that's essentially what `matpredict detect`'s own localize-then-polish
pipeline already does, via a different, already-built mechanism (tblastn
localization + exonerate/miniprot polishing). cblaster's value to this
project would be narrow: a from-scratch cluster search against genomes NOT
already processed by `matpredict detect`, which is redundant with existing
pipeline capability unless there's a reason to use cblaster's specific
remote-NCBI-BLAST convenience.

**Unresolved from documentation alone**: whether clinker's alignment/identity
computation uses protein or nucleotide sequence, and whether it strictly
requires a `CDS` feature with a `translation` qualifier or can work from
`gene`-only GenBank records with an externally-supplied translation. The
fetched README states it "automatically extract[s] protein translations"
and performs "global alignments... using the aligner built into BioPython"
but doesn't specify feature-type requirements precisely enough to commit
code against without a hands-on trial (see Section 6).

## 3. Single-locus gene-structure diagrams — compared

| Tool | Input | Bioconda? | Scriptable | Custom per-gene color | Maintenance |
|---|---|---|---|---|---|
| **DNA Features Viewer** (Python) | GenBank/GFF or BioPython `SeqRecord` objects directly | Yes (`bioconda::dna_features_viewer`) | Yes — pure Python API, matplotlib-based | Yes, per-feature color is a plain constructor kwarg on each `GraphicFeature` | Actively maintained (Edinburgh Genome Foundry); MIT-style license, real bioRxiv preprint (2020) |
| **pyGenomeViz** (Python) | GenBank AND GFF3 | Yes (`conda-forge::pygenomeviz` — note: conda-forge, not bioconda, but same channel list this project already uses) | Yes — both a Python API and a CLI, explicitly documented for "genome analysis scripts/workflow" | Yes — `fc=` per-feature and a `--feature_type2color` CLI flag for category-based coloring | Actively maintained; v1.0.0 released 2024-05, MIT license |
| **Biopython `GenomeDiagram`** | BioPython `SeqRecord` (so GenBank/GFF-parsed-into-SeqRecord) | Ships with Biopython itself (already likely a transitive dependency here via `Bio.Seq`/`Bio.SeqIO` already used in `gff_export.py`) | Yes, pure Python | Yes, per-feature color arguments | Maintained as part of core Biopython, but this specific submodule is older/less actively enhanced than the two above |
| **gggenomes** (R) | GFF3/GenBank/custom tibbles | Not applicable (CRAN/R, not conda-Python) | Yes — it's a ggplot2 extension, fully scriptable in an R script; general web commentary characterizes it as more "interactive/exploratory" in practice, but this is a soft framing, not a technical scriptability limitation | Yes, via ggplot2's own aesthetic mapping (very flexible) | Actively developed (Thomas Lin Pedersen / ggplot2-adjacent ecosystem) |
| **genoPlotR** (R) | Custom R objects built from GenBank/GFF import functions | Not applicable | Yes, scriptable R | Yes, per-segment color lists | Older, lower recent development velocity than gggenomes |

**For MATPredict specifically**: DNA Features Viewer and pyGenomeViz are the
two realistic choices — both are Python, both bioconda/conda-forge-installable
(matching this project's existing `pixi.toml` channel setup with zero new
channel additions), both take GenBank/GFF3 input MATPredict already produces
(modulo the CDS/sequence gap noted in Section 1), and both support arbitrary
per-gene coloring, so `role`/`gene_class`/idiomorph can drive color directly.
The R options add a cross-language dependency this project doesn't currently
have (no R anywhere in `pixi.toml`) for no clear capability gain over the
Python options for the single-locus case.

## 4. Group/synteny comparison diagrams — compared

| Tool | Input | Bioconda? | Scriptable | Custom per-gene color | Fit for MAT-locus scale | Maintenance |
|---|---|---|---|---|---|---|
| **clinker** (standalone) | GenBank or GFF3+FASTA, multiple files | Yes (`bioconda::clinker-py`) | Yes — CLI, non-interactive | Yes (`--colour_map`, `--gene_functions`) | Good — purpose-built for exactly this scale (a handful of genes per cluster, a handful of clusters), unlike whole-genome collinearity tools | Active (2020 paper, ongoing issue tracker) |
| **pyGenomeViz (synteny mode)** | GenBank/GFF, multiple genomes, needs an alignment tool for the identity links (BLAST/MUMmer/MMseqs/progressiveMauve) | Yes (`conda-forge::pygenomeviz`) | Yes | Yes | Good, same scale fit as clinker; adds a dependency on one of the alignment backends (BLAST is already a `pixi.toml` dependency here) | Active |
| **MCScan/MCScanX + SynVisio** | MCScanX's own `.collinearity`/`.gff` block-synteny output; SynVisio is a separate interactive web/JS viewer for that output | MCScanX itself: available via bioconda; SynVisio is a JS web app, not a Python/conda package | MCScanX itself is scriptable (a C program); SynVisio is primarily an interactive browser tool, not obviously scriptable into a static-image pipeline step | Not really designed for per-gene custom attribute coloring — it's built for large-scale genomic collinearity block detection (whole chromosomes/genomes), not small few-gene-cluster comparison | **Poor fit** — this toolchain is designed for a different scale problem (large-scale synteny blocks across whole genome assemblies), not a handful of MAT genes in one small locus; using it here would be forcing a mismatched tool | Actively used in the plant/comparative-genomics community, but SynVisio specifically has fewer recent commits than clinker/pyGenomeViz |
| **Easyfig** | GenBank | No confirmed bioconda recipe found in this research pass | GUI-first (Python2-era tool); has a command-line mode but the ecosystem around it appears less actively maintained | Limited | Poor fit for automation — a legacy, largely GUI-oriented tool | Not confirmed as actively maintained in this research pass — flag as an open question rather than a confirmed claim |
| **AliTV** | Custom JSON (built from GenBank/other via a converter script) | Not conda-distributed (standalone Perl/JS tool) | Partially — has a command-line JSON-generation step, but the viewer itself is a browser tool | Some, via its JSON config | Designed more for whole-genome dot-plot-style comparisons than small gene clusters | Lower recent activity than clinker/pyGenomeViz in this research pass |

**For MATPredict specifically**: **clinker** and **pyGenomeViz's synteny mode**
are the two tools genuinely fit for this problem's actual scale (a MAT locus
is a handful of genes spanning a few kb to tens of kb — not a whole-genome
collinearity problem). MCScanX/SynVisio, Easyfig, and AliTV are all built
for a different scale of comparison and would be a poor structural fit even
where technically usable.

## 5. Recommendation

**Single-locus diagrams**: **pyGenomeViz**, over DNA Features Viewer, for one
practical reason specific to this project: it is the SAME tool recommended
for the group/synteny case below, so choosing it for both use cases means
one dependency, one API to learn, and visual consistency between a
single-locus figure and a multi-locus comparison figure (useful if a
curator wants to go from "here's this one locus" to "here's how it compares
to 3 others" without switching rendering engines). DNA Features Viewer
remains a reasonable fallback if pyGenomeViz's GFF3/GenBank parsing turns
out to have friction with MATPredict's specific GFF3 flavor (Section 1's
per-segment/per-contig fragmented-locus model is somewhat unusual and
untested against either library).

**Group/synteny diagrams**: **clinker**, as the primary recommendation, with
pyGenomeViz's synteny mode as a credible alternative if clinker's
GenBank-with-real-sequence requirement (Section 2's open question) turns out
to be a hard blocker. clinker is purpose-built for exactly MATPredict's
scale (small gene clusters, not whole-genome synteny), is bioconda-installable
with zero new channels, supports the exact per-gene custom-coloring and
function-labeling MATPredict's `role`/`gene_class` fields need, and produces
publication-quality output (an interactive SVG, plus the ability to export
static images) with essentially no manual layout work.

**cblaster**: NOT recommended as part of this pipeline. Its core capability
— searching a genome database for co-located homologs — duplicates what
`matpredict detect`'s existing tblastn-localize + exonerate/miniprot-polish
pipeline already does, via a mechanism already tuned to this project's own
curated reference set and family-routing logic. Adopting cblaster would add
a second, redundant cluster-discovery pathway for no clear capability gain.

**What adaptation MATPredict's writers would need** (both recommendations
depend on this): `gff_export.write_genbank` needs two real changes before
either curated-record output or detection output can feed clinker/pyGenomeViz
meaningfully:
1. Write the SEGMENT's real nucleotide sequence (currently all-`N`
   placeholder) — the real sequence is fetchable via the same
   `NcbiClient.fetch_nucleotide_sequence` machinery already used elsewhere
   in this codebase (`db/validate.py`'s `_independent_translation` already
   does exactly this fetch for a gene's own span; the segment-level fetch is
   the same mechanism at a wider scope).
2. Add a real `CDS` feature (not just `gene`) per present gene, ideally
   carrying a `translation` qualifier from the gene's own already-available
   protein sequence (from `proteins.faa`, or from `_independent_translation`
   when no accession exists) and `role`/`gene_class` as qualifiers so
   clinker's `--colour_map`/`--gene_functions` CSVs (or pyGenomeViz's
   `--feature_type2color`) can be generated directly from this project's own
   schema fields without a separate manual mapping step.
`src/MATPredict/detect/report.py`'s `write_detection_gff3` would need the
analogous treatment for detection-result visualization (real sequence +
CDS features), though it currently emits GFF3 rather than GenBank — clinker
accepts GFF3+FASTA directly, so this may need less rework than the curated
side's GenBank writer.

## 6. Open questions / needs a hands-on trial

- **Clinker's exact feature-type/sequence requirement** (protein vs.
  nucleotide identity computation, `CDS`-with-`translation` vs. plain `gene`
  tolerance) could not be confirmed from README-level documentation alone —
  needs either reading clinker's actual Python source (`clinker/parsers.py`
  or equivalent) or a small real trial run against a hand-built GenBank file
  with a real CDS/translation to see what it actually accepts and rejects.
- **cblaster's real bioconda availability** was not confirmed — the fetched
  README only showed `pip install cblaster`; worth a direct `bioconda`
  channel search before treating pip-only as final (moot given the "not
  recommended" call above, but worth knowing if cblaster is reconsidered
  later for a different reason).
- **The exact cblaster→clinker session-file handoff mechanism** (if cblaster
  were ever reconsidered) was referenced only obliquely by the CAGECAT
  webserver mention, not confirmed via a primary source describing the real
  command/file format.
- **Whether pyGenomeViz's GFF3 parser tolerates MATPredict's specific
  fragmented-multi-contig-locus GFF3 shape** (per-segment `MAT_locus` parent
  features, `locus_group` attribute linking segments across contigs) is
  unconfirmed — this is a real, project-specific format that generic GFF3
  parsers are not guaranteed to handle gracefully; needs a real trial against
  an actual `detected_loci.gff3` file from a completed detection run.
- **Easyfig's and AliTV's actual current (2026) maintenance status** could
  not be confirmed with a real recent-commit-date lookup in this research
  pass — the comparison table's characterization of them as lower-priority
  should be read as "not competitive with clinker/pyGenomeViz for this
  project's scale," not as a confirmed claim of abandonment; a direct
  GitHub commit-history check would firm this up if it matters for a final
  decision.

## 7. Task 3 follow-up: real install + hand-trial results (2026-09-18)

This section resolves this research doc's two open questions plus a third
one raised during Task 1/2 code review, against real, freshly-generated
Task-1 output. Nothing here is inferred from documentation alone — every
claim below was reproduced against real files in this repo, or confirmed by
reading the actually-installed library source.

### 7.1 Installation — real package names/versions

```toml
# pixi.toml [dependencies]
clinker-py = ">=0.0.32,<0.0.33"
pygenomeviz = ">=1.7.0,<2"
```

- **`clinker` (bare name) does NOT resolve** on conda-forge/bioconda for
  this project's channel/platform combination — conda-forge's own
  `clinker` package is an unrelated tool that pins `python 2.7.*`, which
  cannot be solved alongside this project's `python >=3.11`. The correct
  bioconda package name is **`clinker-py`** (matching the research doc's
  prediction), confirmed by `pixi add clinker-py` resolving cleanly to
  `clinker-py 0.0.32`. `pixi run clinker --version` reports `clinker
  v0.0.32`.
- **`pygenomeviz`** resolved as predicted, from conda-forge, to `1.7.0`
  (`pixi run python -c "import pygenomeviz; print(pygenomeviz.__version__)"`
  → `1.7.0`). Note this is the v1.x API (class-based `GenomeViz` /
  `FeatureTrack` / `FeatureSegment` / `parser.Genbank` / `parser.Gff`), not
  the older function-based API some older pyGenomeViz examples online show
  — anything written against this library in Tasks 4-5 must target the
  v1.x class API.
- `pixi run pytest -v` after the `pixi.toml` change: **254 passed**, 0
  failed — the dependency addition alone did not break anything.

### 7.2 Real trial data used

Three real curated records, regenerated fresh via
`matpredict curate-db build-gff` (which fetches real NCBI sequence and
writes real `CDS`/`translation` features per Task 1's merged work):

- `Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1` (Coccidioides
  posadasii MAT1-1, EF512013.1) — 4 genes, 3 of them genuinely multi-exon
  (COX13: 5 exons, APN2: 6 exons, MAT1-1-4: 6 exons; MAT1-1-1 is 2 exons).
- `Ascomycota/Onygenales/199306_silveira_MAT_MAT1-2` (Coccidioides
  posadasii MAT1-2) — 1 gene (MAT1-2-1), used only as a second file for the
  clinker multi-file trial.
- `Ascomycota/Teloschistales/2903220_liq80xsp_MAT_combined` (Xanthoria sp.,
  homothallic, both idiomorphs) — the real fragmented-across-2-scaffolds
  record (JALAIM010000088.1 carries APN2/MAT1-2-1/MAT1-1-1, and
  JALAIM010000140.1 carries SLA2 alone), with a real 2-exon MAT1-1-1
  re-derivation from this session's earlier exonerate work.

### 7.3 Open question 1 (from Task 1/2 review): does the single-genomic-span
CDS (introns included) vs. correctly-spliced `translation=` qualifier
mismatch actually break either tool?

**Confirmed real bug in the writer, real behavior in the tools — the two
tools respond very differently:**

`gff_export.write_genbank` (`src/MATPredict/db/gff_export.py:114-141`)
builds the `CDS` feature's `FeatureLocation` from `gene["start"]`/
`gene["end"]` (the gene's own outer bounds), not from `gene["exons"]` —
confirmed directly by reading the freshly-generated
`db/Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1/locus.gbk`:

```
CDS             complement(885..1580)
                /gene="COX13"
                /translation="MFLQRSVIRAAQGGASRLLYSRVPLGLQRRSMASESKLRNLPDIK..."
```

885..1580 is 696 nt (232 codons), but the real `translation=` qualifier is
only 140 aa — correctly spliced from COX13's 5 real exons, confirming the
location and the translation genuinely disagree, exactly as flagged.

**clinker's actual behavior (read from
`clinker/classes.py:648`, the installed source):**

```python
translation = qualifiers.pop("translation", sequence.translate() if sequence.defined else "")
```

Python evaluates a function's default argument *eagerly*, before `.pop()`
runs, so `sequence.translate()` — a naive translation of the raw
intron-included CDS span — executes on **every** gene regardless of
whether a real `translation` qualifier exists. This is exactly why the
real trial run printed:

```
BiopythonWarning: Partial codon, len(sequence) not a multiple of three. ...
```

for APN2 and MAT1-1-1 (both real multi-exon genes whose intron-included
span length isn't a multiple of 3). **But `.pop()`'s actual *return value*
is still the real qualifier when the key is present** — the eagerly-computed
naive translation is thrown away unless `translation` is genuinely absent.
Verified directly:

```
COX13     clinker_matches_qualifier=True  warned=False len=140
APN2      clinker_matches_qualifier=True  warned=True  len=606
MAT1-1-4  clinker_matches_qualifier=True  warned=False len=246
MAT1-1-1  clinker_matches_qualifier=True  warned=True  len=389
```

A self-vs-self clinker alignment run (same GenBank file copied and aligned
against itself) confirmed clinker's real percent-identity computation is
correct end-to-end for all 4 genes including the two multi-exon ones:

```
Query     Target    Identity  Similarity
COX13     COX13     1.00      1.00
APN2      APN2      1.00      1.00
MAT1-1-4  MAT1-1-4  1.00      1.00
MAT1-1-1  MAT1-1-1  1.00      1.00
```

**Verdict for clinker: works as-is for correctness.** The
`translation=` qualifier is trusted over the (also computed, but discarded)
naive span translation, so clinker's alignment/identity output is NOT
corrupted by the intron-included span. The only real side effect is a
spurious/confusing `BiopythonWarning` printed once per multi-exon gene on
every run (cosmetic noise, not a data bug) — worth suppressing or
documenting in Tasks 4-5 so a user doesn't mistake it for a real problem.

**pyGenomeViz's actual behavior — this IS the visual bug the reviewers
predicted.** `write_genbank`'s CDS feature is a plain `SimpleLocation`
(single span), not a `CompoundLocation` (multi-part, one part per exon).
Confirmed directly:

```
['APN2'] <class 'Bio.SeqFeature.SimpleLocation'> [SimpleLocation(ExactPosition(1918), ExactPosition(4008), strand=1)]
```

pyGenomeViz has a purpose-built `add_exon_features` method (distinct from
plain `add_features`) specifically meant to draw multi-part
`CompoundLocation` features as separate exon boxes joined by intron lines
— but it operates on whatever `feature.location.parts` actually contains,
and MATPredict's CDS location has exactly one part (the whole span). A real
render (`single_locus.png`, using `add_exon_features` on all 4 real CDS
features) draws every gene, including 6-exon APN2, as a single unbroken
arrow spanning the full genomic range — intronic sequence is drawn
identically to coding sequence, with no visual indication of intron/exon
structure anywhere. This is a real, visually-confirmed instance of exactly
the risk flagged: **a viewer reader would see one solid "gene" box and have
no way to tell it isn't a single continuous ORF.**

**Verdict for pyGenomeViz: needs adaptation.** `write_genbank` would need
to emit a real multi-part `CompoundLocation` for the `CDS` feature (built
from `gene["exons"]`, same conversion-to-segment-relative-coordinates logic
the function already does for the flat span) before pyGenomeViz's
`add_exon_features` can draw real intron/exon structure. Until that change
lands, any single-locus diagram built from `write_genbank`'s current output
will visually misrepresent every multi-exon gene as intron-free. This is a
real, concrete, non-theoretical finding — not a guess.

### 7.4 Open question 2: does pyGenomeViz's GFF3 parser tolerate
MATPredict's fragmented-multi-contig-locus GFF3 shape?

Tested directly against the real, freshly-regenerated
`db/Ascomycota/Teloschistales/2903220_liq80xsp_MAT_combined/locus.gff3`:

```
##gff-version 3
##sequence-region JALAIM010000088.1 44505 53097
##sequence-region JALAIM010000140.1 36837 40643
JALAIM010000088.1  MATPredict  gene  44505  46498  ...  Name=APN2
JALAIM010000088.1  MATPredict  gene  47471  48625  ...  Name=MAT1-2-1
JALAIM010000088.1  MATPredict  gene  51979  53088  ...  Name=MAT1-1-1
JALAIM010000140.1  MATPredict  gene  36837  40643  ...  Name=SLA2
```

Two independent, real, reproduced problems, both confirmed against this
exact file (not a synthetic one):

**(a) Silent multi-seqid truncation.** `pygenomeviz.parser.Gff`'s own
docstring says plainly: *"If `target_seqid` is specified when the Gff
instance initialized, then the features of the target seqid are extracted.
Otherwise, extract the features of the seqid in the first row."* Calling
`Gff(path)` with no `target_seqid` on this real file silently returns only
APN2/MAT1-2-1/MAT1-1-1 (scaffold `...088.1`) — SLA2, on the second
scaffold, is silently dropped, with **no error or warning of any kind**.
The only way to get SLA2 is a second, separate `Gff(path,
target_seqid="JALAIM010000140.1")` call. **This confirms the parser is
single-seqid-oriented by design**; any caller (Tasks 4-5) must iterate over
every distinct seqid the record's segments reference and instantiate `Gff`
once per seqid, never relying on a single call to pull in a whole
fragmented locus.

**(b) Hard failure from non-1-based `##sequence-region` windows (a real
bug to design around, found only by reading real error output, per the
project's own instruction not to guess).** Even after fixing (a) by
iterating per-seqid, adding the real extracted features to a
`FeatureTrack`/`FeatureSegment` sized from `Gff.get_seqid2size()` raises:

```
pygenomeviz.exception.FeatureRangeError: feature_location='[44504:46498](+)'
is invalid (segment_range='0 - 8593')
```

Root cause, confirmed from `Gff.get_seqid2size`'s own docstring ("size is
defined by `##sequence-region` pragma of target seqid") plus direct
testing: pyGenomeViz computes a segment's size as `end - start + 1` from
the `##sequence-region` pragma (here: `53097 - 44505 + 1 = 8593`), but it
does **not** rebase the extracted `SeqFeature` locations to be relative to
that pragma's `start` — the features keep their raw, absolute GFF3
coordinates (44504-based, not 0-based-within-the-8593-window). Since
`write_gff3`/`write_detection_gff3` both declare a genuinely-windowed
`##sequence-region` (the true absolute genomic start on the real scaffold,
e.g. `44505`, not `1`), this is a **hard, unconditional failure**, not an
edge case — reproduced identically on the real file, and the same failure
would occur for `write_detection_gff3`'s output too, since its docstring
states it explicitly follows `write_gff3`'s same absolute-coordinate,
per-contig `##sequence-region` convention (`src/MATPredict/detect/report.py:7-9,86-91`).

**Confirmed real workaround** (tested and reproduced working, not assumed):
rewriting the `##sequence-region` pragma to state the seqid's *actual full
length starting at 1* (e.g. `##sequence-region JALAIM010000088.1 1 53097`
instead of `44505 53097`) makes `Gff.get_seqid2size()` return the full
scaffold length, which matches the features' own absolute coordinates, and
the render succeeds with no error. This means Tasks 4-5 have two viable
paths: (1) have MATPredict's writers emit `##sequence-region <seqid> 1
<full_scaffold_length>` instead of the windowed absolute range (requires
knowing the real full contig length, which `NcbiClient`/the rollout genome
FASTA already has available), or (2) bypass the `Gff` pragma-driven sizing
entirely and construct `FeatureTrack`/`FeatureSegment` objects directly in
Python from parsed coordinates, never relying on pyGenomeViz's own
`##sequence-region`-based size inference for a windowed/sliced-locus GFF3.

**Positive finding, also confirmed for completeness:** custom, non-standard
attributes and feature types used by `write_detection_gff3`'s fragmented-
locus model — `MAT_locus` as a feature type, `locus_group=`, `Parent=`
pointing at a same-contig segment ID — all parse and round-trip through
`Gff.extract_features()` with zero issues once the coordinate-rebasing
problem above is worked around (tested with a synthetic file built in this
record's exact shape). The multi-contig-fragmentation problem is entirely a
coordinate-rebasing issue (7.4b above and the single-seqid default in
7.4a), not an attribute/feature-type compatibility issue.

**Verdict for pyGenomeViz + fragmented loci: needs adaptation** — both (a)
per-seqid iteration in the calling code and (b) either a writer-side
`##sequence-region` fix or bypassing `Gff`'s pragma-based sizing — before
Tasks 4-5 can hand a fragmented-locus GFF3 to pyGenomeViz directly.

### 7.5 Summary for Tasks 4-5

| | clinker | pyGenomeViz |
|---|---|---|
| Accepts current `write_genbank` output at all | Yes, cleanly | Yes, cleanly |
| Multi-exon CDS-span-vs-translation mismatch | Harmless — real `translation=` is trusted for alignment; only a spurious per-gene `BiopythonWarning` (cosmetic) | **Real visual bug** — draws intron-included span as one solid box; needs a `CompoundLocation` CDS built from `gene["exons"]` before real exon/intron diagrams are trustworthy |
| Fragmented multi-contig locus GFF3 | Not tested here (clinker's own GFF3 mode needs a same-named companion FASTA per-cluster, not exercised in this trial — clinker was trialed via GenBank, which has no multi-contig ambiguity since each segment is its own SeqRecord) | **Two real, confirmed problems**: silently drops all but the first seqid unless the caller iterates `target_seqid` explicitly; hard `FeatureRangeError` from windowed (non-1-based) `##sequence-region` pragmas unless the writer emits full-scaffold-length pragmas or the caller bypasses `Gff`'s own sizing |

Net effect: clinker can be adopted essentially as designed once fed real
GenBank files (which Task 1 now produces). pyGenomeViz needs two concrete,
named adaptations before it is safe to build application code against: (1)
a real multi-part `CompoundLocation` CDS feature in `write_genbank` (and
the GFF3 writers, for the single-locus GFF3+FASTA path) so exon/intron
structure renders correctly, and (2) either a `##sequence-region`
full-scaffold-length fix in the GFF3 writers, or Tasks 4-5's own calling
code constructing `pygenomeviz` `FeatureTrack`/`FeatureSegment` objects
directly rather than relying on `Gff`'s pragma-driven auto-sizing for a
windowed/fragmented locus.
