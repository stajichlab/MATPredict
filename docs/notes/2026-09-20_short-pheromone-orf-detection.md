# Detecting short pheromone-precursor genes that tblastn cannot localise

Research note, not an implementation plan. Investigated 2026-09-20.

Every quantitative claim below is tagged as one of:

- **[SOURCE]** — stated by a cited primary source.
- **[MEASURED]** — computed here, in this repo, against real curated data or a
  real genome. Commands and exact numbers are in the appendix.
- **[INFERRED]** — my reasoning on top of the two above. Stated as inference.

---

## 1. The problem

MATPredict Stage 1 localisation is a genome-wide `tblastn` of curated MAT
proteins. A gene whose best curated reference protein is shorter than 60 aa is
flagged `not_searchable=true` (`_short_orf_genes`,
`src/MATPredict/detect/pipeline.py:249`; floor `short_orf_aa_floor=60` at
`pipeline.py:878`). Ten of 75 genes in the curated DB fall below that floor,
and eight are Basidiomycete pheromone precursors. `pheromone_precursor` is the
largest `gene_class` in the database.

**[MEASURED]** Across the whole curated DB there are 181 proteins; 20 of them
carry `gene_class: pheromone_precursor` (17 distinct sequences — some are exact
paralogous duplicates). Their lengths:

| length (aa) | 38 | 40 | 42 | 58 | 61 | 63 | 64 | 65 | 69 | 83 | 85 | 92 | 93 | 96 | 123 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| count | 3 | 1 | 3 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 2 | 1 | 1 | 1 | 1 |

So the problem is real but narrower than "all pheromone precursors are
untblastnable": the median is 65 aa. Only 8 of 20 are actually below the 60 aa
floor. The rest are short-but-searchable, and are a useful bridge.

**[MEASURED] The companion receptor is not short.** All 9 curated
`pheromone_receptor`/`STE3` proteins are 357–629 aa (357, 380, 385, 397, 417,
484, 576, 626, 629). This matters enormously for the recommendation in §7: the
7TM receptor sitting next to the precursor is comfortably tblastn-searchable.

---

## 2. Headline empirical result: the CAAX check

This was the falsifiable check the task asked for, and it passes decisively.

**[SOURCE]** Spellig, Bölker, Lottspeich, Frank & Kahmann (1994), *EMBO J*
13:1620-7, purified both *Ustilago maydis* a1 and a2 pheromones biochemically
and determined their structures. Verbatim from the abstract: *"Both pheromones
are post-translationally modified by farnesylation and carboxyl methyl
esterification of the C-terminal cysteine."* Also: *"The structure of the
secreted pheromones was determined to be 13 amino acids for a1 and nine amino
acids for a2"*, and *"An unmodified a1 peptide exhibits dramatically reduced
activity."*
[DOI](https://doi.org/10.1002/j.1460-2075.1994.tb06425.x) (retrieved from
PubMed; abstract and PMC landing page read — see §9 on full-text access).

This is primary biochemical evidence, not a review assertion: the modification
was determined on purified natural pheromone, and the unmodified peptide was
shown to be far less active. The C-terminal cysteine is the CAAX cysteine.

**[SOURCE]** Vaillancourt, Raudaskoski, Specht & Raper (1997), *Genetics*
146:541-51, describe the *Schizophyllum commune* Bβ1 locus as containing *"at
least three pheromone genes and one pheromone receptor gene"*, and the Bα1
locus as encoding *"three lipopeptide pheromones and a pheromone receptor with
a seven-transmembrane domain."*
[DOI](https://doi.org/10.1093/genetics/146.2.541) (PubMed).

### 2.1 What the real curated sequences actually show

**[MEASURED]** All 20 curated `pheromone_precursor` proteins, aligned on their
stop codon. `CAAX strict` = `C-[AVLIM]-[AVLIM]-X` at the exact C-terminus.

| gene | len | last 6 | C at −4 | CAAX strict |
|---|---|---|---|---|
| MFa1 / MFa2 / MFa3 | 42 | YSCVIA | yes | yes |
| MFalpha1 / MFalpha2 / MFalpha3 | 38 | TLCVIA | yes | yes |
| mfa1 | 40 | SSCVVA | yes | yes |
| bap3-1 | 58 | AFCVVA | yes | yes |
| bap3-3 | 63 | YFCVVA | yes | yes |
| bbp2-1 | 85 | GWCVVA | yes | yes |
| bbp2-2 | 65 | GWCIIA | yes | yes |
| bbp2-5 | 92 | GYCVVM | yes | yes |
| bbp2-6 | 85 | YFCEVM | yes | **no** (E at −3) |
| bbp2-8 | 93 | GYCVVA | yes | yes |
| bbp2_a | 83 | SYCTIA | yes | **no** (T at −3) |
| bbp2_b | 96 | GYCVVG | yes | yes |
| pheromone_B43 | 64 | AFCIIA | yes | yes |
| pheromone_B44 | 61 | AFCVIV | yes | yes |
| fungal_mating_type_pheromone | 69 | WFCVIA | yes | yes |
| fungal_mating_type_pheromone | 123 | SSIPSL | **no** | **no** |

**Answer to the question asked: yes, overwhelmingly.**

- **Cysteine at exactly position −4: 16 of 17 distinct sequences (94%).**
- **Strict CAAX at the exact C-terminus: 17 of 20 proteins / 14 of 17 distinct
  sequences (82%).**
- The two "near misses" (bbp2-6 `CEVM`, bbp2_a `CTIA`) still have the CAAX
  cysteine at −4 and differ only in the aliphatic-`a` positions. Under a
  relaxed `a` alphabet they pass. **[INFERRED]** These are real CAAX boxes with
  a non-canonical `a` residue, not absences — which matches the prior art's
  decision to relax the `a` alphabet (§4).
- **The single genuine exception is the 123 aa
  `fungal_mating_type_pheromone`** in
  `db/Basidiomycota/Agaricales/5346_a43-b43-okayama-7_PR_B43`. It has no
  C-terminal cysteine at all and ends `...PFVISILFMSALFIVSSIPSL`.

> **Curation flag, offered as a finding rather than a change (this note touched
> nothing under `db/`).** That 123 aa sequence's C-terminus is a ~20-residue
> uninterrupted hydrophobic stretch, which reads as a transmembrane helix, not
> as a farnesylation substrate. It is also nearly twice the length of every
> other precursor in the DB. **[INFERRED]** It is more likely a fragment of a
> 7TM receptor, or a mis-assigned ORF, than a pheromone precursor. It is worth
> a curator's eye. Note it is the *only* curated precursor that breaks an
> otherwise 16/17 rule, and it measurably degrades every motif learned from
> this set (§6.1).

### 2.2 There is more conserved signal than just the four CAAX residues

**[MEASURED]** Positional residue counts over the last 16 aa of the 17 distinct
precursors:

| position | −5 | −4 | −3 | −2 | −1 |
|---|---|---|---|---|---|
| top residues | F(6) Y(4) S(3) | **C(16)** | V(12) I(2) | V(8) I(8) | A(12) M(2) |

and positions −12…−6 are glycine/serine-rich (G at −12 in 7/17, G at −10 in
7/17, G at −8 in 5/17, G at −6 in 5/17).

**[MEASURED]** MEME 5.5.5 run on the 17 distinct precursors (`-protein -mod
zoops -nmotifs 3 -minw 5 -maxw 15`) independently recovers this:

| motif | sites | E-value |
|---|---|---|
| `GSGNGTGFCVVA` | 16/17 | 2.3e-40 |
| `MDSFTTL` | 13/17 | 3.0e-12 |
| `EAPRNQE` | 5/17 | 8.5e-03 |

So there are **three** real motifs, not one: an extended ~12-residue C-terminal
motif (the CAAX box plus a G-rich linker), a strong **N-terminal `MD[SA]F`
motif**, and a weaker internal one. §6.1 shows that the extended version is,
counterintuitively, *not* the most useful for detection.

---

## 3. Tooling: what is actually available in this environment

Checked by `command -v` on the bare login PATH, inside `pixi run -e default`,
and against the HPCC environment-modules tree. Nothing was installed.

| tool | bare PATH | pixi `default` env | HPCC module | verdict |
|---|---|---|---|---|
| `tblastn` (blast) | missing | **present** (`.pixi/envs/default/bin/tblastn`) | — | already a project dependency |
| `diamond`, `miniprot`, `exonerate` | missing | **present** | — | already project dependencies |
| `phmmer` / `hmmbuild` / `hmmsearch` | missing | **missing** | **`hmmer/3.4`** (also 3.3.2, 2.3.2, 1.8.5, `-mpi` variants) | usable today via module; would need adding to `pixi.toml` for pipeline use |
| `getorf` / `seqret` / `transeq` (EMBOSS) | missing | missing | **`emboss/6.6.0`** | usable today via module |
| `meme` / `mast` / `fimo` | missing | missing | **`meme/5.5.5`** (also 5.4.1) | usable today via module |
| `prodigal` | missing | missing | **`prodigal/2.6.3`** | available, but see §5 — wrong tool |
| `seqkit` | missing | missing | **`seqkit/2.4.0`** | available |
| `orfipy` | missing | missing | **no module** | would need installing; not evaluated hands-on |
| `torch` / `transformers` / `esm` | — | **all missing** | not checked | ESM2 not runnable today (§6.5) |

**Environment gotcha worth recording.** `module load X | tail -2` silently
fails to change `PATH`: `module` is a shell function, and piping it runs it in
a subshell that discards the environment change. The first EMBOSS trial here
failed with `/bin/bash: line 11: seqret: command not found` for exactly this
reason, while `module list` reported EMBOSS as loaded. Never pipe a `module
load`.

**`pixi.toml` note.** Everything the recommendation in §7 needs beyond the
current dependency set is `hmmer` and `emboss` (or a Python ORF finder), both
of which are on conda-forge/bioconda and would slot into the existing
`[dependencies]` table. Nothing was added in this investigation.

---

## 4. Prior art: this exact pipeline is published and validated

This is the single most important literature finding, and it was not in the
original hypothesis list.

**[SOURCE]** Srikant, Gaudet & Murray (2023), *Current Biology*
33:4098-4110.e3, "Extending the reach of homology by using successive
computational filters to find yeast pheromone genes."
[DOI](https://doi.org/10.1016/j.cub.2023.08.039) (PubMed; full text read via
PMC10592104).

They faced precisely MATPredict's problem, in the Ascomycete yeasts, and state
it plainly in the abstract: *"Peptide pheromones have been found by genetics or
biochemistry in a small number of fungi, but their short sequences and modest
conservation make it impossible to detect homologous sequences in most
species."*

Their four-step pipeline, verbatim from the abstract: *"we require that
candidate genes have a C-terminal prenylation motif, are shorter than 100 amino
acids long, and contain a proteolytic-processing motif upstream of the
potential mature pheromone sequence and that closely related species contain
highly conserved homologs of the potential mature pheromone sequence."*

Details quoted from the PMC full text:

- They did **not** use a strict regex: *"we relaxed this CAAX detection
  dictionary based on experimental studies of the combinatorial signatures of
  farnesylation, proteolysis and carboxymethylation."*
- ORF enumeration is **stop-anchored, not start-anchored**: candidates are
  those *"that have an in-frame methionine within 6 to 100 residues upstream of
  the stop codon."*
- On the multiple-start problem: *"Collapsing sets of candidate sequences that
  result from multiple Start codons upstream of a single CAAX to one sequence
  from the most upstream Start codon results in 87,326 unique farnesylated
  loci."*
- The proteolytic filter is *"an in-frame proteolytic-motif asparagine (N)
  important for the final step of maturation"*, with *"the C terminus of the
  mature pheromone ... between 5–20 residues downstream from the N-terminal
  residue X."*
- Funnel: 284,073 CAAX-stop candidates → 125,711 after the asparagine filter →
  87,326 unique loci after collapsing starts → *"500–2000 candidates per
  genome after initial filters"* → 812 manually curated candidates across 241
  species, from 332 genomes.
- Validated experimentally in *Yarrowia lipolytica*: deleting all four
  predicted genes in the a-mating type *"prevents mating."*

**Why this matters for MATPredict.** Three things transfer directly:

1. The overall architecture (ORF-first, CAAX-anchored, then filter) is
   published and experimentally validated. This is not a speculative design.
2. **Their per-genome burden after motif filtering was 500–2000 candidates.**
   Motif alone is nowhere near sufficient. They needed cross-species
   conservation *and* manual curation to get to 812 final calls across 241
   species. **[MEASURED]** my own independent number on a real 31.5 Mb fungal
   genome is 1,689 CAAX-stop ORFs — the same order of magnitude, which
   cross-validates both.
3. Their stop-anchored enumeration is the correct answer to a trap I fell into
   and measured (§5.2).

**Scope caveat, stated as a caveat:** Srikant et al. worked in Saccharomycotina
(Ascomycetes), not Basidiomycetes. The a-factor chemistry is shared, but their
specific relaxed dictionary and the asparagine filter were tuned on yeast
a-factors. **[MEASURED]** §6.1 shows the asparagine-equivalent processing-site
filters perform poorly on this repo's Basidiomycete set.

---

## 5. The Prodigal correction, and what to use instead

### 5.1 Prodigal is a prokaryotic gene finder — confirmed from its own docs

**[SOURCE]** The Prodigal README (github.com/hyattpd/Prodigal, fetched
2026-09-20) describes it verbatim as *"Fast, reliable protein-coding gene
prediction for prokaryotic genomes."* The documentation makes no mention of
eukaryotes, fungi, or introns.

**[INFERRED]** Using Prodigal on a fungal genome would be a category error for
two compounding reasons: its gene model has no intron concept, and — more
subtly — its scoring is trained on prokaryotic coding statistics (GC-frame
bias, Shine-Dalgarno ribosome binding site motifs) that simply do not apply to
a fungal genome. It would also *reject* candidates, which is the opposite of
what is wanted here: for this task we do not want a gene *finder* that makes
accept/reject decisions at all. We want an exhaustive, decision-free ORF
*enumerator* that hands every possibility downstream. That is a different
category of tool.

`prodigal/2.6.3` is available as an HPCC module. It should not be used for
this.

### 5.2 EMBOSS `getorf` works, is fast, and has one dangerous flag

**[MEASURED]** `getorf` (EMBOSS 6.6.0) on a real 31.5 Mb, 712-contig fungal
genome (`testset/Zygo/query/Absidia_blakesleeana_NRRL_1300.gbk`, which is
actually *Lichtheimia blakesleeana*, Mucoromycota):

| invocation | ORFs | wall time |
|---|---|---|
| `-find 1 -minsize 90 -maxsize 360` | 362,408 | 2.3 s |
| `-find 1 -minsize 90` (no cap) | 384,442 | 2.1 s |

Two seconds for a whole fungal genome. Performance is a non-issue. `seqret`
converted the GenBank file to FASTA in 1.3 s.

**The `-maxsize` trap — measured, not theorised.** The curator's hypothesis
suggested restricting to ORFs < 100 aa. **This actively destroys true
positives.** All three *Cryptococcus neoformans* `MFalpha1/2/3` genes are 38 aa
and intronless (gene span 117 bp = 38 codons + stop), yet they were **absent**
from the `-maxsize 360` ORF set:

| getorf options | ORFs from the jec21 locus | MFalpha recovered |
|---|---|---|
| `-find 1 -minsize 90 -maxsize 360` | 1,560 | **0 / 3** |
| `-find 1 -minsize 90` | 1,692 | **3 / 3**, inside a 139 aa ORF |
| `-find 0 -minsize 90` | 2,758 | 3 / 3, inside a 151 aa ORF |

The reason: `getorf -find 1` starts each ORF at the *first* ATG after the
preceding stop. For MFalpha that upstream in-frame ATG sits ~100 codons
further back, so the reported ORF is 139 aa — and a 120 aa cap discards it
*entirely*, taking the real 38 aa pheromone with it. The peptide is present in
a plain 6-frame `transeq` translation of the same locus (3 copies, one per
paralog), confirming the sequence is there and it is `getorf`'s reporting
that dropped it.

**[INFERRED]** This is exactly the artefact Srikant et al. avoid by anchoring
on the stop codon and collapsing starts, rather than filtering on ORF length.
The lesson is concrete and prescriptive: **never cap ORF length. Apply the
length constraint as the distance from a candidate in-frame Met to the stop
codon, not as a property of the ORF record.**

Conveniently, `getorf -find 1` already *is* the "most upstream start codon"
collapse that Srikant et al. describe, so `-find 1 -minsize 90` with no cap
gives the right enumeration directly.

### 5.3 Alternatives considered

- **`orfipy`** — no HPCC module, not installed, **not trialled**. Its
  advantage over `getorf` is claimed speed and native BED output. Given
  `getorf` processes a whole genome in 2.1 s, **[INFERRED]** speed is not a
  reason to prefer it here. No recommendation either way without a trial.
- **Biopython** (`Seq.translate` over 6 frames, split on `*`) — no new
  dependency (Biopython is already pulled in via the existing stack), full
  control over the stop-anchored enumeration, and trivially able to record
  exact genomic coordinates and strand for each ORF. **[INFERRED]** For a
  pipeline that must emit GFF3 coordinates anyway, this is likely preferable
  to shelling out to `getorf` and parsing its `[start - end]` headers — and it
  avoids adding EMBOSS as a dependency for one function. The measured result
  above is what validates the *approach*; the implementation does not have to
  be `getorf`.

---

## 6. The five hypotheses, assessed against measurement

The evaluation below uses a **leave-one-record-out (LOO)** protocol wherever it
reports sensitivity: for each of the 6 Basidiomycete records that contain
precursors, the method is built from the *other* records only, then tested on
the held-out record's own locus ORFs. This is the honest measure of "can this
find a pheromone in a species we have not curated?" A hit counts only if an ORF
carrying the true protein's C-terminal 15 aa is retrieved.

Two enabling measurements first.

**[MEASURED] ORF extraction recovers every curated precursor.** Using
`getorf -find 1 -minsize 90` (no cap) on each record's own `locus.gbk`:
**20/20 precursors have their C-terminal 15 aa present in the ORF set.** This
includes the 4 precursors whose gene span exceeds their CDS length (i.e. that
contain an intron or annotated UTR: `bbp2-1`, `pheromone_B44`, and both
`fungal_mating_type_pheromone` entries) — because the intron is not in the
terminal exon, the C-terminus survives intact in a single reading frame. Only
8/20 match the *full* protein exactly, for the start-codon reason in §5.2.

**[INFERRED]** This is the load-bearing result for the whole approach: it says
a naive, annotation-free ORF scan does not lose the signal, *provided* the
method keys on the C-terminus rather than on the full-length protein. Any
design that requires matching the complete precursor will lose 60% of them.

### 6.1 Hypothesis 2 — motif-based. **Best of the methods tested.**

**[MEASURED] Specificity**, on the 384,442-ORF *Lichtheimia* genome (a
Mucoromycota genome, where Basidiomycete-style B-locus pheromones are not
expected, so essentially all hits are false positives):

| filter | genome hits | rate |
|---|---|---|
| `C-[AVLIM]-[AVLIM]-X` (textbook CAAX) | 1,689 | 0.439% |
| `C-[VI]-[IV]-[AVMG]` (tightened on the curated data) | **72** | 0.019% |
| tightened + Met 6–100 residues upstream of stop | 70 | 0.018% |
| tightened + R/K in the −20…−5 window | 56 | 0.015% |
| tightened + `[ED]R` in the −20…−5 window | 8 | 0.002% |

Two things fall out. **Tightening the two `a` positions cuts false positives
23-fold at no cost in recall** — both the loose and tight forms retrieve 17/20
in-sample. And **the Srikant-style Met-window filter is redundant here**
(1,689→1,676; 72→70), because `getorf -find 1 -minsize 90` has already imposed
a Met start. Do not implement it; it buys nothing.

**[MEASURED] Sensitivity, leave-one-record-out.** Learning the four
C-terminal residue classes from the training records only (union of observed
residues), then testing on the held-out record: **14/20 (70%)**.

The LOO motifs learned, and their genome-wide false-positive burden, expose the
outlier problem concretely:

| held-out record | motif learned from the rest | recall | genome FPs |
|---|---|---|---|
| jec20 / jec21 / a1 / Balpha_3 | `[CI][EIPTV][ISV][AGLMV]` | 9/9 | 853 |
| Bbeta_2 | `[CI][IPV][ISV][ALV]` | 3/7 | 447 |
| PR_B43 | `[C][EITV][IV][AGM]` | 2/4 | **78** |

The `I` at position −4 and the sloppy `a` alphabets come entirely from the 123
aa outlier flagged in §2.1. **[INFERRED]** Dropping that one record from the
training set would tighten every learned motif and cut the false-positive
burden by roughly an order of magnitude, at no recall cost — which is a second,
independent reason to have a curator look at it.

**Verdict: promising, and the single most sensitive method tested.** But 72
candidates per genome is far too many to report as findings. It is a *recall*
filter, not a detector.

### 6.2 Hypothesis 1 — ORF-first, then profile search. **Half-right.**

The ORF-first half is correct and essential (§6, 20/20 recovery). The
profile-search half **underperforms the plain regex**, which was not the
expected result.

**[MEASURED] `phmmer` (HMMER 3.4), curated precursors as queries.**

- *Specificity is excellent.* Against all 384,442 genome ORFs: **0 hits at
  E<0.01**, 7 at E<1, 18 at E<10. Runtime 3.4 s.
- *Sensitivity is mediocre.* Leave-one-record-out at E<1: **10/20 (50%)**.

**[MEASURED] A C-terminal-anchored HMM.** Since all precursors align trivially
on their stop codon, I built a 20-column `hmmbuild` model from the right-
anchored last 20 aa and searched it against the last 20 aa of every ORF. This
is the natural way to exploit the extended MEME motif from §2.2.

- *Specificity is superb:* **0 hits at E<0.1** in 384,442 ORFs; 2 at E<1.
- *Sensitivity is poor:* leave-one-record-out, **5/20 (25%)** — the worst of
  the three methods.

**[INFERRED]** The extended C-terminal motif is genuinely conserved *within* a
genus and largely absent *between* genera: the G-rich −12…−6 linker that makes
the MEME motif look so strong (E = 2.3e-40) is fitted across all 17 sequences
at once, but its actual residues differ between *Cryptococcus*, *Coprinopsis*,
*Schizophyllum* and *Ustilago*. A profile model spends its information budget
on those divergent columns and dilutes the only signal that transfers — the
four CAAX residues. This is the mechanistic reason the 4-residue regex beats
the 20-column HMM, and it is consistent with Srikant et al.'s framing that
homology detection fails for these peptides beyond ~100 My of divergence.

**[MEASURED] Union of motif OR phmmer, leave-one-record-out: 15/20 (75%)** —
better than either alone (14 and 10). `phmmer` uniquely rescues `bbp2-6`, the
`CEVM` non-canonical CAAX. So profile search earns its place as a *parallel*
recall channel, not as the primary detector.

Per-gene LOO detail:

| record | gene | motif | phmmer | union |
|---|---|---|---|---|
| jec20 | MFa1/2/3 | yes | yes | yes |
| jec21 | MFalpha1/2/3 | yes | yes | yes |
| a1 | mfa1 | yes | no | yes |
| Balpha_3 | bap3-1, bap3-3 | yes | yes | yes |
| Bbeta_2 | bbp2-1, bbp2-2, bbp2-8 | yes | no | yes |
| Bbeta_2 | bbp2-6 | no | **yes** | yes |
| Bbeta_2 | bbp2_a, bbp2_b, bbp2-5 | no | no | **no** |
| PR_B43 | pheromone_B43, fmtp-69 | yes | partly | yes |
| PR_B43 | pheromone_B44, fmtp-123 | no | no | **no** |

The 5 LOO failures cluster in *Schizophyllum* Bβ2 (the most paralog-rich, most
divergent locus) and in the two PR_B43 entries, one of which is the suspect
outlier.

### 6.3 Hypothesis 3 — genomic context / proximity to other MAT genes. **The decisive filter.**

This is the hypothesis that converts an unusable method into a usable one, and
the numbers make the case cleanly.

**[MEASURED]** The motif's problem is burden, not recall: 72 candidates across
a whole 31.5 Mb genome. But across all 10 curated Basidiomycete loci —
4,252 ORFs in total, i.e. the ORFs that actually sit inside a MAT locus — the
tightened motif yields only **20 candidates**, and `phmmer` at E<1 only **14**.
Per locus that is roughly **2–6 candidates**, against 1–7 true precursors.

**[MEASURED]** And the anchor to define those windows already exists and is
already searchable: all 9 curated pheromone receptors are 357–629 aa (§1), far
above the 60 aa floor, so Stage 1 `tblastn` localises them today with no
changes.

**[SOURCE]** The biology justifies the linkage as a rule, not a coincidence.
Raudaskoski & Kothe (2010), *Eukaryotic Cell* 9:847-59: in agaricomycetes *"two
mating type loci, A, coding for homeodomain type transcription factors, and B,
encoding a pheromone/receptor system, regulate the four typical mating
interactions of tetrapolar species."*
[DOI](https://doi.org/10.1128/EC.00319-09) (PubMed). Vaillancourt et al. 1997
report the Bβ1 locus carrying three pheromone genes *and* a receptor gene, and
Bα1 likewise — pheromones and their receptor are clustered in the same locus.
[DOI](https://doi.org/10.1093/genetics/146.2.541) Fowler, Mitton, Vaillancourt
& Raper (2001), *Genetics* 158:1491-503, describe the Bβ2 locus as *"a large
cluster of genes encoding a single pheromone receptor and eight different
pheromones."* [DOI](https://doi.org/10.1093/genetics/158.4.1491)

That last figure — one receptor, eight pheromones — is worth dwelling on. It is
the curated `5334_h4-8_Bbeta_2` record in this repo (7 precursors), it is the
record where the motif does worst (3/7 LOO), and it tells you the expected
answer per locus is "several", not "one".

**Verdict: adopt.** **[INFERRED]** Receptor-anchored windowing is what makes
the whole thing work. It is also cheap: it reuses Stage 1 output, adds no
dependency, and the biology says the precursor is essentially always in the
same locus as the receptor.

One honest limitation: **[INFERRED]** this cannot find an orphan pheromone gene
far from any receptor, and the C. cinerea record's own `definition_note` warns
that this species has ~14 receptors and ~29 pheromones scattered genome-wide,
most of them paralogous clusters outside the mating-type locus. Receptor
anchoring will therefore surface *paralogous* pheromone/receptor clusters too.
That is a precision problem for locus assignment, not a recall problem.

### 6.4 Hypothesis 4 — promoter motifs. **Unpromising. Do not pursue now.**

**Not tested here**, and I am recommending against testing it, with reasons
rather than a shrug.

**[INFERRED]**, with the supporting source noted:

1. *It solves the wrong problem.* A pheromone-response element is a target of
   pheromone signalling, so it occurs upstream of the whole regulon, not
   specifically upstream of the precursor genes. As a genome-wide filter it
   would *add* candidates.
2. *A short degenerate DNA motif cannot be specific in a 31.5 Mb genome.* A
   6–10 bp site occurs by chance thousands of times. The measured CAAX filter
   already gives 0.019%; a promoter motif cannot improve on that and would
   dilute it.
3. *The regulatory logic is documented to be indirect and condition-dependent.*
   **[SOURCE]** Hartmann, Krüger, Lottspeich & Kahmann (1999), *Plant Cell*
   11:1293-306, show that in *U. maydis* the pheromone regulator Prf1 is itself
   controlled by carbon source and cAMP, and that *"pheromone and cAMP
   signaling regulate prf1 post-transcriptionally"* — i.e. the promoter layer
   is several steps removed from "is this ORF a pheromone gene."
   [DOI](https://doi.org/10.1105/tpc.11.7.1293)
4. *There is no training data.* n=17, spanning four genera with unaligned
   upstream regions, and this repo's `locus.gff3` files carry only `gene`
   features with genome-absolute coordinates, so extracting correctly-oriented
   upstream windows is itself non-trivial work before any motif search starts.

**[MEASURED]** supporting point 4: the `locus.gff3` files contain `gene`
features only — no `CDS`, `mRNA` or `exon` rows — so promoter extraction has no
TSS to anchor on.

### 6.5 Hypothesis 5 — ESM2 embeddings. **Not viable yet; revisit later.**

**[MEASURED]** Not runnable in this environment today: `torch`,
`transformers`, `fair_esm` and `esm` are all absent from the pixi env.

**[INFERRED]**, three reasons to defer rather than to try:

1. *n = 17 distinct sequences, 4 genera.* Any embedding-space decision boundary
   fitted on 17 points, then applied to 384,442 ORFs per genome, is a
   4-orders-of-magnitude extrapolation. Even a 0.1% error rate is 384 false
   positives — worse than the 72 the regex already achieves.
2. *These proteins are largely disordered.* The precursors are Ser/Thr/Pro/Ala-
   rich with no globular domain; the functional payload is a 9–13 aa
   farnesylated peptide (**[SOURCE]** Spellig et al. 1994,
   [DOI](https://doi.org/10.1002/j.1460-2075.1994.tb06425.x)). pLM embeddings
   encode structural and evolutionary context, which is precisely what these
   sequences lack. **[INFERRED]** the signal ESM2 would have to find is the
   CAAX box — which a four-character regex already finds for free.
3. *The measured ranking argues against it.* Going from a 4-residue regex
   (14/20 LOO) to a 20-column HMM (5/20 LOO) made things *worse* by adding
   model capacity over divergent columns. A 650M-parameter model has vastly
   more capacity to overfit the same 17 examples.

**[INFERRED]** The precondition for revisiting is data, not compute: if curation
grows the precursor set to, say, 100+ sequences across 15+ genera with held-out
genera available for validation, an ESM2 check becomes a reasonable experiment.
The MEMORY note on ESM2 as a future direction stands; this is just not the
problem where it pays off first.

---

## 7. Recommended approach

A four-stage funnel. Stages 1–2 are recall; stages 3–4 are precision.

**Stage A — anchor on the receptor (reuses existing Stage 1).**
Take `tblastn` hits for curated `pheromone_receptor`/`STE3` proteins (357–629
aa, already searchable). Define a window of ±20–30 kb around each. **[MEASURED]**
justification: curated PR/B loci span 3.9–20.7 kb, and the two *Cryptococcus*
MAT loci span 133–148 kb, so the window must be configurable, and for
*Cryptococcus*-like large loci it should be driven by the localised locus
bounds rather than a fixed offset.

**Stage B — enumerate ORFs in the window, with no upper length cap.**
Stop-anchored, all 6 frames, `-minsize 90` (30 aa) and **no `-maxsize`**
(§5.2). Either `getorf -find 1 -minsize 90` or ~30 lines of Biopython; the
latter is probably better since it keeps coordinates and strand for GFF3 output
without header parsing. **[MEASURED]** cost: ~2 s for a whole genome, so a
windowed run is free. Expected yield ~250 ORFs per 20 kb locus.

**Stage C — recall filter: motif OR profile, unioned.**
- C-terminal regex `C-[VI]-[IV]-[AVMG]` on the ORF's last 4 residues.
- **OR** `phmmer` of all curated precursors against the window's ORFs, E<1.

**[MEASURED]** LOO recall of the union is **15/20 (75%)**; the motif alone is
14/20, `phmmer` alone 10/20. **[MEASURED]** expected candidates surviving per
real curated locus: **2–6**.

**Stage D — rank and report, do not silently accept.**
Score each survivor by, in rough order of measured value: CAAX match quality
(exact `C-[VI]-[IV]-[AVMG]` > relaxed `a` positions); `phmmer` E-value;
distance to the anchoring receptor; presence of the N-terminal `MD[SA]F` motif
(**[MEASURED]** 13/17 curated, MEME E = 3.0e-12); and **near-identity to another
candidate in the same window** — **[MEASURED]** the `jec20` locus has 3
*identical* MFa copies and `jec21` has 2 identical MFalpha copies, so
within-locus duplication is a genuine positive signal, and it is the same
observation Srikant et al. exploited (*"many species carry more than one
a-factor gene, encoding identical or nearly identical pheromones"*).

Replace `not_searchable=true` for these genes with a ranked candidate list at
an explicitly lower confidence tier. **[INFERRED]** Given 75% LOO recall and
2–6 candidates per locus, presenting these as *curator-facing candidates* is
honest and useful; presenting them as confident gene calls would not be.

### Expected failure modes (each tied to a measurement)

| failure mode | evidence |
|---|---|
| **~25% of true precursors missed outright.** LOO recall 15/20. Failures concentrate in the most paralog-rich, most divergent loci (*Schizophyllum* Bβ2: 4/7 missed). | [MEASURED] §6.2 |
| **Intron in the terminal exon breaks it.** 20/20 recovery held only because no curated precursor has an intron interrupting its C-terminus. 4/20 do have introns elsewhere. A precursor spliced within its last ~15 codons is invisible to this method. | [MEASURED] §6 |
| **Paralogous non-MAT clusters will be surfaced.** *C. cinerea* has ~14 receptors and ~29 pheromones genome-wide; receptor anchoring cannot tell a real B-locus from a dispersed paralogous cluster on its own. | [SOURCE] the curated record's own `definition_note`; §6.3 |
| **Orphan pheromones far from any receptor are unreachable by design.** Stage A cannot window what it cannot anchor. | [INFERRED] §6.3 |
| **Un-anchored (whole-genome) fallback is unusable as a reporting mode.** 72 motif candidates and ~1,689 loose-CAAX candidates per genome; Srikant et al. independently report 500–2000 per genome. Only ever run this mode to produce curation candidates. | [MEASURED] §6.1; [SOURCE] §4 |
| **Non-canonical CAAX `a` residues are missed by the tight regex.** `bbp2-6` (`CEVM`) and `bbp2_a` (`CTIA`) both fail it; only `phmmer` rescued `bbp2-6`. This is why the union in Stage C is not optional. | [MEASURED] §2.1, §6.2 |
| **The curated set is small and taxonomically narrow.** 17 distinct sequences, 4 genera, 3 orders. Every motif and every LOO number here carries that sampling. | [MEASURED] §1 |

---

## 8. What remains unknown

1. **True recall on a real Basidiomycete genome is not measured.** Every
   sensitivity number here comes from running against curated *locus* sequence
   (3.9–148 kb), not a whole Basidiomycete genome. The specificity numbers come
   from a Mucoromycota genome used as a negative control. **No single test here
   measured recall and precision simultaneously on one genome.** The right next
   experiment is to download the *C. cinerea* Okayama-7 and *C. neoformans*
   JEC21 assemblies and run the full Stage A–D funnel end to end.
2. **The false-positive rate is estimated from one genome of the wrong phylum.**
   *Lichtheimia* (Mucoromycota) is a reasonable negative control but its amino
   acid composition and GC content differ from a Basidiomycete's; 72 could
   plausibly be 30 or 200 elsewhere. Unmeasured.
3. **Whether the 123 aa `fungal_mating_type_pheromone` is correctly classified.**
   Flagged in §2.1. It is the sole counterexample to a 16/17 rule and it
   measurably degrades every learned motif. A curator call, not mine.
4. **Whether the Srikant asparagine/proteolytic filter has a Basidiomycete
   equivalent.** The nearest analogues tested here (`R/K` or `[ED]R` upstream)
   cost far more recall than they bought: `[ED]R` cut genome hits from 72 to 8
   but recall from 14/17 to 5/17. A proper Kex2-site analysis on a larger
   curated set might do better. Untested.
5. **`orfipy` was never run.** No module, not installed. Its BED output might be
   more convenient than `getorf` headers, but given 2.1 s/genome the question
   is ergonomics, not performance.
6. **Whether the Ascomycete short genes generalise.** `cha1` (22 aa) and
   `mat1-Mi` (42 aa) were explicitly out of scope. `cha1` at 22 aa is short
   enough that even ORF enumeration at `-minsize 90` (30 aa) would miss it —
   the floor would have to drop, and the false-positive count would rise
   sharply. Unexamined.
7. **Full text of two key primary sources could not be retrieved.** See §9.

---

## 9. Source-access caveats

Stated explicitly, per the rule that a paywalled paper must never be summarised
from its abstract as though read.

- **Vaillancourt et al. 1997** *(Genetics 146:541-51,
  [DOI](https://doi.org/10.1093/genetics/146.2.541))* and **Fowler et al. 2001**
  *(Genetics 158:1491-503,
  [DOI](https://doi.org/10.1093/genetics/158.4.1491))*: the PubMed MCP
  `get_full_text_article` tool returned `"full_text": ""` for both PMC records
  (PMC1207996, PMC1461750). A direct fetch of the PMC landing page for
  PMC1207996 returned only the abstract, with the note that the body *"is only
  available as a PDF download."* **Everything cited from these two papers above
  is from their abstracts, and is quoted as such.** I did not read their
  methods, and in particular I did **not** verify from their full text how they
  originally identified the small pheromone ORFs.
- **Srikant et al. 2023** *(Curr Biol,
  [DOI](https://doi.org/10.1016/j.cub.2023.08.039))*: full text **was**
  successfully retrieved via PMC10592104, and the methods quotes in §4 are from
  that full text.
- **Spellig et al. 1994** *(EMBO J,
  [DOI](https://doi.org/10.1002/j.1460-2075.1994.tb06425.x))*: retrieved via
  PMC394992. The quoted sentences on farnesylation, carboxyl methyl
  esterification and peptide lengths are verbatim. I did **not** obtain a
  verbatim in-text statement of the literal string "CAAX" from this paper; the
  claim I draw from it is the biochemical one it does state — farnesylation and
  carboxyl methylation of the C-terminal cysteine — which is the CAAX
  modification by definition.
- The bioRxiv preprint of Srikant et al. (2021.09.28.462209) returned **HTTP 429
  Too Many Requests** and was not read; the published PMC version was used
  instead.
- Article metadata and abstracts in this note were retrieved **from PubMed**.

---

## Appendix A — trials run, with exact commands and results

All work was done in `$SCRATCH/caax_trial` (node-local scratch). Nothing was
installed; nothing under `src/`, `tests/` or `db/` was modified; no tool output
was committed.

### A.1 Tool availability

```
$ for t in tblastn phmmer hmmsearch hmmbuild meme getorf orfipy prodigal seqkit; do
    printf "%-10s %s\n" "$t" "$(command -v $t || echo MISSING)"; done
# bare PATH: ALL MISSING
# inside `pixi run -e default`: only tblastn present
$ module avail hmmer   -> hmmer/1.8.5 2.3.2 3.3.2 3.4(default) 3.3.2-mpi 3.4-mpi
$ module avail emboss  -> emboss/6.6.0
$ module avail meme    -> meme/5.4.1 meme/5.5.5
$ module avail prodigal-> prodigal/2.6.3
$ module avail seqkit  -> seqkit/2.4.0
$ module avail orfipy  -> (no output — not available)
```

### A.2 The `module load` piping failure (verbatim)

```
$ module load emboss/6.6.0 2>&1 | tail -2 ; which seqret
/usr/bin/which: no seqret in (...)
/bin/bash: line 11: seqret: command not found
```
`module` is a shell function; the pipe ran it in a subshell and the `PATH`
change was discarded. Without the pipe:
```
$ module load emboss/6.6.0 ; command -v seqret
/opt/linux/rocky/8.x/x86_64/pkgs/emboss/6.6.0/bin/seqret
```

### A.3 Genome conversion and ORF extraction

```
$ seqret -sequence testset/Zygo/query/Absidia_blakesleeana_NRRL_1300.gbk \
         -sformat genbank -outseq $W/genome.fa -osformat fasta
real 0m1.283s   -> 712 contigs, 31,516,467 bytes

$ getorf -sequence $W/genome.fa -outseq $W/orfs_30_120.faa \
         -find 1 -minsize 90 -maxsize 360 -table 1
real 0m2.283s   -> 362,408 ORFs

$ getorf -sequence $W/genome.fa -outseq $W/orfs_nomax.faa \
         -find 1 -minsize 90 -table 1
real 0m2.080s   -> 384,442 ORFs
```

### A.4 The `-maxsize` trap, traced to the sequence

```
$ transeq -sequence $W/loci/40410_jec21_MAT_alpha.fa -outseq $W/jec21.6f.faa -frame 6 -clean
$ grep -c 'GGMTLCVIA' $W/jec21.6f.faa        -> 3   (the 3 MFalpha paralogs, present)
$ grep -c 'GGMTLCVIA' <getorf -maxsize 360 output>  -> 0   (all three lost)
```
| getorf options | ORFs | MFalpha-containing ORFs |
|---|---|---|
| `-find 1 -minsize 90 -maxsize 360` | 1560 | 0 |
| `-find 1 -minsize 90` | 1692 | 3 (139 aa each) |
| `-find 0 -minsize 90` | 2758 | 3 (151 aa each) |

### A.5 CAAX check against curated data

Gene classes were read from `db/*/order.yml`; sequences from
`db/*/*/*/proteins.faa` (59 files, 181 proteins, 20 with
`gene_class: pheromone_precursor`).

```
CAAX strict (C-[AVLIM]-[AVLIM]-X) at exact C-terminus: 17/20
Cys at exactly position -4 (distinct sequences):       16/17
has Cys anywhere in last 6 aa:                         19/20
```

### A.6 Recovery from locus ORF sets

```
total precursors: 20
gene span > 3*(len+1)  (intron or UTR present):   4/20
C-terminal 15 aa present in ORF set (no maxsize): 20/20
full protein present exactly:                      8/20
```

### A.7 phmmer negative control

```
$ phmmer --cpu 4 -E 1000 --tblout g.tbl -o /dev/null pher_queries.faa orfs_nomax.faa
real 0m3.442s
hits E<10: 18 | E<1: 7 | E<0.01: 0        (of 384,442 ORFs)
```

### A.8 C-terminal-anchored HMM

`hmmbuild --amino` on a Stockholm alignment of the right-anchored last 20 aa of
all 17 distinct precursors; searched against the last 20 aa of every ORF.

```
genome negative control (384,442 ORFs): E<0.1 -> 0 hits ; E<1 -> 2 ; E<10 -> 4
leave-one-record-out recall:            E<0.01 -> 5/20 ; E<1 -> 5/20
```

### A.9 MEME

```
$ meme pher_queries.faa -protein -oc meme_out -nmotifs 3 -minw 5 -maxw 15 \
       -mod zoops -markov_order 0
MOTIF GSGNGTGFCVVA  width=12  sites=16  E-value=2.3e-040
MOTIF MDSFTTL       width= 7  sites=13  E-value=3.0e-012
MOTIF EAPRNQE       width= 7  sites= 5  E-value=8.5e-003
```

### A.10 Leave-one-record-out summary (the headline comparison)

| method | LOO recall (n=20) | genome hits (of 384,442) |
|---|---|---|
| C-terminal regex, classes learned per-fold | **14/20** | 78–853 (fold-dependent) |
| tightened regex `C-[VI]-[IV]-[AVMG]` (all-data) | 17/20 in-sample | **72** |
| `phmmer`, E<1 | 10/20 | 7 |
| C-terminal HMM, E<1 | 5/20 | 2 |
| **union: regex OR phmmer** | **15/20** | — |

### A.11 Within-locus paralogy

```
40410_jec20_MAT_a              n=3  distinct=1  identical_pairs=3
40410_jec21_MAT_alpha          n=3  distinct=2  identical_pairs=1
5334_h4-8_Bbeta_2              n=7  distinct=7  identical_pairs=0
5346_a43-b43-okayama-7_PR_B43  n=4  distinct=4  identical_pairs=0
```

### A.12 ESM2 feasibility

```
$ pixi run -e default python -c "import torch"          -> MISSING
$ ... transformers / fair_esm / esm                     -> MISSING
```
