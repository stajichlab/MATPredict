<!--
  Template for the About page. Edit this file and rebuild or redeploy the site; about.html renders it in the browser.
  Supported Markdown: # ## ### headings, paragraphs, - and 1. lists, > quotes, ``` code blocks, **bold**, *italic*,
  `code`, [text](url), bare https:// links, --- rules. HTML comments like this one are not shown.
  Replace anything marked [UPDATE] and delete this comment block when you are done.
-->

# About MATPredict

MATPredict finds and types the mating-type (MAT) locus in fungal genome assemblies and, where a good reference panel
exists, in raw reads. This site summarises how it works, how it performs, and which MAT alleles it finds across the
public [BFD](https://github.com/stajichlab/Fungi_BFD/) fungal genomes and the ZyGoLife Mucoromycotina genomes.

## How MATPredict works

MATPredict searches a genome assembly with a curated database of reference MAT loci, then reports where each locus is,
which genes it carries and which idiomorph (mating type) it is, with the evidence for each call.

1. **Route.** The genome's NCBI taxonomy picks the curated families to search (`db/<Phylum>/order.yml`). A genome outside
   every curated family falls back to its phylum, or is reported as not searched.
2. **Search.** The reference proteins (core MAT genes and the flanking genes) are searched against the genome with
   translated BLAST (tblastn), or against a proteome with DIAMOND.
3. **Cluster.** Hits on a contig are grouped into candidate loci using gap rules for each lineage.
4. **Polish.** Each candidate locus is modelled with exonerate and miniprot, on a window sliced around it.
5. **Verify.** A candidate must pass checks on its modelled genes, the fraction of the reference it covers, the presence of a
   MAT gene, flanking-gene evidence and paralog class (and the CAAX motif for pheromone receptors).
6. **Type.** The idiomorph comes from gene content, or from a profile-HMM classifier (Mucoromycota Plus/Minus), and each call
   gets a confidence tier (high, medium, low) with the reason it was not higher.

Every genome gets a report, including a reason when no locus was called. A genome with two idiomorphs is flagged for
review, not resolved. A polished gene model with a frameshift or an internal stop codon is marked defective; this is
reported only and does not change the confidence.

**Curation.** Each reference record is one strain and one idiomorph, and traces to a publication or a validated
GenBank deposit. Tier 1 is an experimentally validated, published locus. Tier 2 is a genome-derived locus admitted with
extra checks (domain evidence, alignment to curated orthologs or a gene tree). A candidate record is proposed with its
accession or PMID and the exact source text, checked automatically (the accession resolves, the translation matches the
cited protein), signed off by the curator, and then run through a regression panel before it enters the database. The
database is in the [db folder](https://github.com/stajichlab/MATPredict/tree/main/db); the rules are in the
[README](https://github.com/stajichlab/MATPredict#curation-process).

## How to cite

### This resource

> Stajich JE. Kingdom-wide Mating Type Locus Exploration. https://github.com/stajichlab/MATPredict

Please also give the date you accessed the site, and the version below if you used results from a particular release.

### The software

The command-line software is cited in the project README, under
[Citation](https://github.com/stajichlab/MATPredict#citation). The citation metadata is in
[CITATION.cff](https://github.com/stajichlab/MATPredict/blob/main/CITATION.cff), which GitHub shows under "Cite this repository".

> Stajich JE. MATPredict: MAT locus identification in Fungi. Version 0.8.0. https://github.com/stajichlab/MATPredict

```bibtex
@software{matpredict,
  author  = {Stajich, Jason E.},
  title   = {MATPredict: MAT locus identification in Fungi},
  version = {0.8.0},
  url     = {https://github.com/stajichlab/MATPredict},
  year    = {2026}
}
```

<!-- [UPDATE] When a paper or DOI exists, add it here and in CITATION.cff and the README; change the "No paper" wording. -->
No paper describes MATPredict yet. Cite the software and the release you used.

### The reference records

Every reference record in the database traces to a publication or a validated deposit. Please also cite the publications
behind the records you rely on; each record's `metadata.yaml` in the
[db folder](https://github.com/stajichlab/MATPredict/tree/main/db) lists its source (PMID, DOI or accession).

### Idiomorph names and labels

The database, the reports and the calls carry the names as the references use them (MATa, MAT1-2, A, P and so on). Names such as `Sacch__SacCer_alpha` are *labels*: they are computed on output from the family that made the call and the raw name, using a versioned registry, and they say which system's reference panel matched, not that two loci are homologous. The same raw name can mean different systems in different families, which is why the labels exist. A label is never stored in a record. To cite a call, give the raw name together with the registry version (shown in each download and in the explorer). A † beside a label means some loci behind it may change: they were called by two families, have an undetermined idiomorph, rest on a placeholder system with no curated reference, or on a family call that could not be confirmed (`label_status`). The tables list the reason for each locus.

## Data used on this site

- **[BFD](https://github.com/stajichlab/Fungi_BFD/)**: the Ascomycota, Basidiomycota and Mucoromycotina genomes of this fungal
  genome collection, as frozen in September 2026 (see also the
  [2026 run documentation](https://github.com/1KFG/Fungi_BFD_run2026/)). The v0.7.0 run searched 19,415 Ascomycota,
  3,275 Basidiomycota and 974 Mucoromycotina genomes in October 2026, on release commit `cleanrun-2026-10c` with the
  `--preset v0.7` policy set.
- **ZyGoLife**: the LCG and Jena Mucoromycotina genome sets, held out from training and curation.
- Numbers on the landing page and the explorer come from committed result tables in the
  [results folder](https://github.com/stajichlab/MATPredict/tree/main/results); they are detection counts, not ground truth.

<!-- [UPDATE] Add the ZyGoLife data citation or URL, and any BFD citation you want shown, one bullet each. -->

## What works, and what does not yet

The calls on this site are the output of one release, with its gaps stated here so they are not read as findings.

- **Coverage.** Only three phyla have curated families (Ascomycota, Basidiomycota, Mucoromycotina). An uncalled genome is
  weak evidence that it has no MAT locus: a genome far from every reference often fails the evidence bar, and a
  fragmented assembly can split a locus or lose its flanking genes.
- **Confidence tiers are conservative and rule-based.** Several known rules lower or blur a call and are reported, not
  changed, in this release: the idiomorph margin for loci with a single candidate (#137, #145), the vote on short fragments
  and tiny margins (#138), and a noise third gene breaking a tie (#129). The tables mark the affected loci (`tier_note`,
  `label_unsupported`). A better hit can also lower a confidence because a 39-bit floor for unpolished genes is a cliff
  (#149); this mostly affects *MTLA1* in CTG-clade yeasts.
- **Missed loci.** Known losses include some Hypocreales *MAT1-1-1* calls, a Helotiales idiomorph vote that blocks
  *MAT1-1-1* deposits, *Hyphopichia* and several yeast genera (*Metschnikowia*, *Sungouiella*, *Kurtzmaniella*,
  *Candida*) that lack a full-length *MTLA1* reference, and Basidiomycota pheromone-receptor genes needing introns longer
  than 10 kb (#154; one genome, GCA_900009495.1, stays lost).
- **Genomes not searched.** Seven rust genomes that v0.6.0 called (four *Austropuccinia psidii*, three *Phakopsora pachyrhizi*) and two other
  large rust genomes were skipped in this release because the BFD collection's suppress list excludes them for annotation-pipeline reasons, not
  because they lack a MAT locus. They have no report here. They are planned for the next release (#165).
- **Calls that may be wrong, accepted for this release.** 102 pheromone-receptor loci (87 in Agaricales) that v0.6.0 withheld are now called at
  medium confidence on pheromone-protein hits from one reference record (*Coprinopsis cinerea* B43); some may count a pheromone gene and the
  receptor's CAAX ORF as two genes. Six Basidiomycota *bE*/*bW* models in five genomes have introns of 5 to 10 kb at low coverage and may be
  chimeric (the 10 kb intron cap, #154). Neither has been checked against the genome sequence.
- **Candidate references not used.** Four reference records that would fix some of these (*Aciculoconidium aculeatum* and
  three *Lodderomyces* / *Candida parapsilosis*-group records) are held as candidates and are not used by detection until
  they pass a full regression.
- **Defective genes are not proven pseudogenes.** The flag means a model does not translate cleanly in this assembly. It
  fires on 7.5% of polished Ascomycota genes in a pilot, mostly weak models of distant flanking genes where an alignment
  artefact is likelier than a pseudogene. A real pseudogene shows in most isolates of a species (for example
  *Lodderomyces parapsilosis* *MTLA1*, 97 of 97 loci); an assembly error shows in few.
- **A two-idiomorph result is a flag for review.** It can be homothallism, a mixed culture or a duplication; the
  `homothallic_candidate` class is a candidate, not a statement of homothallism.
- **Policies are a preset.** The v0.7.0 results use `--preset v0.7`; the code defaults are unchanged, so a run without
  the preset gives different calls.

Held-out items and what each needs for the next release are tracked in
[issue #155](https://github.com/stajichlab/MATPredict/issues/155).

## Assessments and further reading

The v0.7.0 results were judged campaign by campaign before release. These notes say what changed against the previous
release, what was checked and what was not, and which decisions are still open. Read them for what the calls can and cannot support.

- **Ascomycota:** [judgement](https://github.com/stajichlab/MATPredict/blob/main/results/2026-10_cleanrun_preparation/v070_ascomycota_judgement.md) (changes against v0.6.0, the truth checks and the reason for each confidence change).
- **Basidiomycota:** [independent judgement](https://github.com/stajichlab/MATPredict/blob/main/results/2026-10_cleanrun_preparation/v070_basidiomycota_judgement.md) (lost and new calls, the long-intron models, the defective-gene flag, the nine open curator decisions).
- **Mucoromycotina:** [judgement](https://github.com/stajichlab/MATPredict/blob/main/results/2026-10_cleanrun_preparation/v070_mucoromycotina_judgement.md).
- **Before the run:** the [regression gate](https://github.com/stajichlab/MATPredict/blob/main/results/2026-10_cleanrun_preparation/regression_judgement.md) and the [pilot](https://github.com/stajichlab/MATPredict/blob/main/results/2026-10_cleanrun_preparation/pilot_judgement.md), and the [clean-run plan](https://github.com/stajichlab/MATPredict/blob/main/docs/notes/2026-10-22_clean-run-plan.md).
- **Benchmarks and held-out sets:** the [held-out benchmark method](https://github.com/stajichlab/MATPredict/blob/main/docs/holdout-benchmark.md), the [held-out sets](https://github.com/stajichlab/MATPredict/blob/main/analysis/2026-09-28_heldout-sets.md) and the [index of analyses](https://github.com/stajichlab/MATPredict/blob/main/analysis/INDEX.md).
- **Validation and known limits:** the project [README](https://github.com/stajichlab/MATPredict#validation-and-known-limits) and the [changelog](https://github.com/stajichlab/MATPredict/blob/main/CHANGELOG.md) (see "Known issues" under 0.7.0).
- **What is planned next:** [issue #155](https://github.com/stajichlab/MATPredict/issues/155) lists everything held out of this release and what each item needs.

## Source, issues and licence

- Source code and documentation: https://github.com/stajichlab/MATPredict
- Questions and bug reports: [GitHub issues](https://github.com/stajichlab/MATPredict/issues)
- Licence: MIT

## Author

Jason E. Stajich, University of California, Riverside.
[ORCID 0000-0002-7591-0020](https://orcid.org/0000-0002-7591-0020)

<!-- [UPDATE] Add funding and acknowledgements here, for example a "## Acknowledgements" section. -->
