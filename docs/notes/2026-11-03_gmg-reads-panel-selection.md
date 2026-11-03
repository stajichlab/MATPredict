# Design: GmG and the reads path (`matpredict reads-type`)

Status: proposal, nothing built. 2026-11-03. Follows `2026-11-02_gmg-phylum-routing.md`, which covers the assembly path (`detect --phylum auto`).

## What the reads path is, and what GmG already does there

`reads-type` types the idiomorph of a strain from raw reads by counting idiomorph-unique k-mers of a **panel** (`--idiomorph NAME=FASTA`, one FASTA per
idiomorph; the panels are species-specific files under `panels/`). It does not search families and does not route by phylum, so `--phylum auto` has no meaning here.
The user (or the survey list) chooses the panel by naming the species.

GmG is already used on this path as a **screen** (#241, `datasets/reads/sra_survey/pilot/screen_reads.py`, calibrated in `results/2026-11-02_contam_spikein/NOTE.md`):
a read set that is another organism or a mixture is still typed by a panel's k-mers, so GmG's record is turned into flags (`genus_mismatch`, `mixed`,
`borderline`, `no_call`) that travel with the result and never stop a run. That covers "is this the species I think it is".

## The remaining gap: choosing the panel when the species is not known

If the reads come with no trustworthy species name (a project with mixed or mislabelled samples), the user must find the species first. GmG can say it, for reads as
for assemblies, from the same record: `guess.calls.genus` and `species` with their taxids, `qc.coverage_est`, `qc_flags`. Proposal, smallest first:

1. **A panel registry** `panels/index.yml`: one entry per panel set: the panel directory, the idiomorph names and files, and the taxon it applies to (an NCBI species or
   genus taxid, as `taxonomic_scope` does in `db/`). The existing panels (*A. fumigatus*, *Coccidioides*, *Cryptococcus* species, *Neurospora*, *Fusarium* species) are the first entries.
2. **`reads-type --panel auto`** (instead of `--idiomorph`): run GmG on the reads (the `gmg` binary, `--max-reads`, `--min-abund auto` as GmG's reads guidance says), take the species
   call, else the genus call, look it up in the registry by lineage, and use that panel. With no registry entry, or no species/genus call, the sample is **not typed** and the
   row says why (`no_panel_for <taxon>` or `gmg_no_call`), the same cautious default as `--phylum auto`. The GmG call and the screen flags go into the output row.
3. **Guard rails from what GmG measured for reads:** reads use a species ANI cutoff of 0.99 (errors lower the ANI), so a species call needs coverage and is often only a genus
   call; at about 0.5x coverage a genus call is reliable, below that `no_call`. A genus-level panel entry is therefore the common case. `weak_phylum` is never set for reads, so
   the assembly `usable` rule does not carry over: for reads the rule is "genus or species called, no `mixed`".

## Not proposed

- Routing reads through `detect`: that is a different method (assembly, tblastn), not this path.
- Typing when GmG says `mixed`: flag it and type nothing for that sample, or type and flag; the screen's current behaviour (flag, never stop) stays.

## Decisions needed (only if this is wanted)

1. Is there a use case for an unknown-species read set, or does the survey always know the species? If the latter, stop here: the screen is enough.
2. If wanted: registry keyed by species taxid, by genus taxid, or both (a genus entry with several species panels needs a rule for which panel).
3. Where the registry lives: `panels/index.yml`, with each panel's README (as `panels/coccidioides/README.md`) linked from its entry.

## Tests (if built)

The same shape as `tests/detect/test_gmg_guess.py`: the decision table with fake GmG records, the registry lookup by lineage (species, then genus, then none), and the
output row for each outcome; one real-reads check on the HPCC with the survey's pilot runs.
