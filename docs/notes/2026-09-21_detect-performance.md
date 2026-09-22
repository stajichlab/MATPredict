# detect performance — where the time actually goes (2026-09-21)

Measured on the live repository, not estimated. Two caches landed; one
suspected problem turned out not to exist.

## What was NOT the problem

**Order-level scoping does not cost more, it costs much less.** The worry was
that scoping all six Agaricales loci to one order taxid makes every Agaricales
genome query six families instead of one or two. Measured query-set sizes from
`build_reference_fasta` against the live DB:

| routing | proteins | total aa |
|---|---|---|
| exhaustive (all phyla) | 268 | 101,138 |
| Ascomycota phylum fallback | 122 | 39,773 |
| Basidiomycota phylum fallback | 105 | 39,598 |
| Mucoromycota phylum fallback | 41 | 21,767 |
| **Agaricales (6 loci)** | **25** | 8,743 |
| **Tremellales (MAT)** | **10** | 2,109 |
| **Ustilaginales + Malasseziales** | **70** | 28,746 |

Order-level routing is 4.2x smaller than the phylum fallback it replaces for
Agaricales, and 10.5x smaller for Tremellales. The alternative was never
"one family" — it was the 105-protein phylum fallback, because no curated
family's scope covered those genomes at all.

**Sequence-level dedup is not worth doing.** Exact duplicate sequences across
the whole 268-protein set: 9. Within any single routed set: at most 3
(Tremellales). Collapsing identical queries saves nothing measurable.

## Problem 1 — the genome was re-indexed on every polish call

`search._extract_window` ran `SeqIO.index` over the WHOLE genome FASTA each
time it was called. It is called twice per gene polished (exonerate and
miniprot are deliberately given the same window so their models are
comparable), and again for every further gene and every further admitted
family in a cluster.

Measured on a synthetic 46 MB, 300-contig assembly — the shape of an ordinary
Agaricomycete:

* `SeqIO.index` + one fetch: **565 ms per call**
* 20 window extractions, re-indexing each time: **11.25 s**
* the same 20 with the index cached: **0.58 s** — **19.4x**

A cluster with 10 genes admitted to 3 families makes ~60 such calls, so this
was roughly 34 s per cluster spent re-reading one file before any aligner
started.

**Fixed** by `search._genome_index`, a bounded cache keyed on
(path, size, mtime_ns), evicting by closing. Deliberately not thread-safe:
`SeqIO.index` objects are not either, and nothing in this package polishes
concurrently within a process.

## Problem 2 — the curated DB was re-parsed once per genome

`run_pipeline` calls `load_record_families` for every genome, and that parses
every curated `metadata.yaml`.

* parse: **466 ms** (101 records)
* stat-only staleness stamp over the same files: **2.3 ms** — **201x cheaper
  to check than to redo**

| genomes | before | after |
|---|---|---|
| 23 (Zygo ground truth) | 10.7 s | 0.5 s |
| 283 (BFD Mucoromycota) | 131.9 s | 1.1 s |
| 3,174 (BFD Basidiomycota) | **1,479.9 s** | **7.8 s** |

**Fixed** by `family_registry._db_stamp` + `_RECORD_FAMILIES_CACHE`, keyed on
(file count, newest mtime_ns, summed size). All three components are needed:
count alone misses an edit, mtime alone misses a file swapped in with an older
timestamp, size catches a same-mtime rewrite of different length. It is a
cache-invalidation heuristic, not a content hash;
`clear_record_families_cache()` is the escape hatch. The cached dict is copied
on return so a caller mutating the result cannot poison it.

## What this note does NOT establish

* **No end-to-end genome run was made.** The previous validation runs lived in
  a node-local `$SCRATCH` that is gone with its job, and no genome was
  re-acquired. Every number above is a component measurement or a synthetic
  genome. The real per-genome wall time after these two caches is UNMEASURED.
* **The 6.6 min/genome Basidiomycota mean is not explained.** These two fixes
  remove a large constant, but nothing here profiles exonerate or miniprot
  themselves, which is where the remaining time most likely sits.
* **Nothing was done about per-(cluster, family, gene) polish fan-out.** That
  is the real combinatorial term and it is untouched.
