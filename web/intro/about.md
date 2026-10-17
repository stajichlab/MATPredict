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

## How to cite

### This resource

> Stajich JE. Kingdom-wide Mating Type Locus Exploration. https://github.com/stajichlab/MATPredict

Please also give the date you accessed the site, and the version below if you used results from a particular release.

### The software

The command-line software is cited in the project README, under
[Citation](https://github.com/stajichlab/MATPredict#citation). The citation metadata is in
[CITATION.cff](https://github.com/stajichlab/MATPredict/blob/main/CITATION.cff), which GitHub shows under "Cite this repository".

> Stajich JE. MATPredict: MAT locus identification in Fungi. Version 0.6.3. https://github.com/stajichlab/MATPredict

```bibtex
@software{matpredict,
  author  = {Stajich, Jason E.},
  title   = {MATPredict: MAT locus identification in Fungi},
  version = {0.6.3},
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

## Data used on this site

- **[BFD](https://github.com/stajichlab/Fungi_BFD/)**: the Ascomycota, Basidiomycota and Mucoromycotina genomes of this fungal genome collection (see also the [2026 run documentation](https://github.com/1KFG/Fungi_BFD_run2026/)), run with MATPredict v0.6.0 in October 2026.
- **ZyGoLife**: the LCG and Jena Mucoromycotina genome sets, held out from training and curation.
- Numbers on the landing page and the explorer come from committed result tables in the
  [results folder](https://github.com/stajichlab/MATPredict/tree/main/results); they are detection counts, not ground truth.

<!-- [UPDATE] Add the ZyGoLife data citation or URL, and any BFD citation you want shown, one bullet each. -->

## Source, issues and licence

- Source code and documentation: https://github.com/stajichlab/MATPredict
- Questions and bug reports: [GitHub issues](https://github.com/stajichlab/MATPredict/issues)
- Licence: MIT

## Author

Jason E. Stajich, University of California, Riverside.
[ORCID 0000-0002-7591-0020](https://orcid.org/0000-0002-7591-0020)

<!-- [UPDATE] Add funding and acknowledgements here, for example a "## Acknowledgements" section. -->
