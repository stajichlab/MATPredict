# Publication-notable findings

A curated, citeable list of findings from the MATPredict curation and
detection work that may merit a place in a publication: first detections,
homothallism candidates, hybrids, assembly or annotation artefacts, method
results and biology.

Started 2026-09-27 at the curator's request. It supersedes
`docs/notes/publication-highlights.md` on branch `curation-puccinio` once the
branches merge; move any new entry from that file here.

## Rules

1. **Every number cites its source**: a results path, a record ID, an
   accession with coordinates, or a commit. A number that cannot be traced
   to a file is marked `unverified`.
2. **Unpublished manuscript content is never copied here.** The group's
   Rhodotorula manuscript (local copy `resource/MBE_202608/`, git-excluded)
   is cited only through its public preprint, bioRxiv
   doi:10.1101/2025.09.11.675505, and only for facts stated there. Per-strain
   assignments from its tables stay out of the repository.
3. Plain language, short sentences. State limits. Do not inflate.
4. Results paths are relative to the repository root. Some result folders
   live only in the main checkout (`results/` is partly untracked); the
   entry says so when a path is not committed.

## Status vocabulary

| status | meaning |
|---|---|
| candidate | the data support it, but an independent check is still open |
| verified | checked against an independent source (reads, literature, known answers, synteny) |
| artefact-explained | an apparent signal traced to an assembly, annotation or code cause |

## Categories

curation-first, homothallism candidate, hybrid, assembly-or-annotation
artefact, method, biology.

## Adding an entry

1. Copy the template below to `NNN-short-slug.md` (next free number).
2. Fill every section; write `none` rather than leaving one out.
3. Add a row to `INDEX.md`.

## Entry template

```markdown
# NNN. Title

- **Category:** ...
- **Status:** candidate | verified | artefact-explained
- **Lineage:** ...

## Summary
Two or three plain sentences.

## Evidence
- Each bullet: accession/coordinates/record ID/results path, with gene,
  identity and coverage where relevant.

## Method that found it
Tool, rule, parameters, commit hash.

## Verification done / still open

## Limits

## Related
Records, notes, other entries.
```
