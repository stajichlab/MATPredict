# Design: guess the phylum with GmG, then route detection

Status: proposal, nothing built. Written 2026-11-02 after reading both repos (MATPredict v0.8.0, GmG v0.2.1).

## Goal

`matpredict detect` should be able to find out for itself which phylum a genome belongs to, so that it searches only that
phylum's curated families (the same narrowing `--phylum` gives today), without the user knowing the phylum or a taxid.
Step 1 is the CLI. Step 2 is a "Guess" button in the browser page that fills the existing Phylum pulldown.

GmG ("Guess My Genome", <https://github.com/stajichlab/GmG>) compares a FracMinHash sketch of the genome with a reference
index of 23,669 BFD genomes and returns a call at each rank (phylum to species), the NCBI taxid of each call, and a status.

## What already exists

MATPredict `route()` (`detect/family_registry.py`) tries, in order: explicit `--phylum`; direct taxid; taxid lineage;
phylum fallback; `not_searched` (or `--exhaustive`). `--phylum` costs no taxonomy lookup. Its choices are the directories in
`db/` (Ascomycota, Basidiomycota, Mucoromycota). The browser page already passes its `#phylum` select to `--phylum`
(`web/spike/page/app.mjs`).

GmG output used here (`results.jsonl` / `gmg guess` JSON; columns in GmG `docs/results-schema.md`):
`status` (`ok`, `low_confidence`, `no_match`, `error`), `call_phylum`, `conf_phylum`, `taxid_phylum`, `taxid_order`,
`taxid_family`, `top_hit_taxid`, `withheld`, `qc_flags` (for example `not_fungal`).

So the whole feature is: run GmG, translate its answer into the arguments `route()` already understands, and record what
was done. No new routing rule in `route()` is needed.

## Decision rules

Given a GmG record:

| GmG result | MATPredict action | `routing_mode` |
|---|---|---|
| `call_phylum` is a phylum in `db/` and passes the confidence rule | route with `phylum=` (and `taxid=` when usable, see below) | `explicit_phylum` or `lineage` (as `route()` decides), plus `phylum_source: gmg` |
| `call_phylum` is a fungal phylum with no curated families (for example Chytridiomycota) | not searched; reason names the GmG call | `not_searched` |
| no phylum call, `low_confidence`, `no_match`, `error` or QC flag `not_fungal` | not searched; reason says why. With `--exhaustive`, search everything, as today | `not_searched` / `exhaustive` |
| `--taxid` or `--phylum` given by the user | the user's value wins; GmG is not run (so no index is needed) | unchanged |

Confidence rule: use GmG's own `status`/`call_phylum` (a call is only present when it passed GmG's calibrated threshold),
and add no second threshold in MATPredict. A wrong phylum is worse than none, because the right family is never searched,
so the default is the cautious `not_searched`.

Taxid use. GmG's `taxid_order` or `taxid_family` lets `route()` take the `lineage` path, which also gives the
subphylum guard, the genetic code and scoped roster genes. A bare phylum loses all three. Proposal: pass the deepest
called taxid when the call exists, and only the phylum otherwise. Open question 3 covers how to test this.

## CLI design

```
matpredict detect --genome G.fa --out-dir out/ --phylum auto [--gmg-index PATH|URL]
```

- `--phylum auto` is added to the existing option's choices (`_phylum_choices()` adds the literal `auto`; the check in
  `_cmd_detect` that rejects a phylum with no families is skipped for `auto` until after the guess).
- New module `detect/gmg_guess.py`, one function `guess_routing(genome, index) -> GmgGuess` (dataclass: phylum, taxid,
  status, confidence, withheld, gmg_version, index_id). It returns data and never raises on a bad guess; an unreadable
  index is a usage error, an unreadable genome is `status=error`.
- `_cmd_detect` calls it before `route()`, then passes `phylum=guess.phylum` and `taxid=guess.taxid`.
- `detection_report.yaml` gets a `phylum_guess:` block (GmG version, index file and sha256, status, call, confidence,
  taxids, what MATPredict did with it). The HTML report shows one line, so a misrouted genome can be traced.
- `matpredict detect-batch` and `run_detection_batch.py` take the same option. A batch should call GmG once per genome and
  group by phylum, as `run_batch(phylum=...)` already expects one phylum per batch.

How MATPredict calls GmG (in order of preference):

1. **Python bindings**: `gmg.Index(path).classify(genome)` returns the same record as the CLI. Lazy import, so `matpredict --help`
   and other commands do not need it. **The name `gmg` on PyPI is an unrelated project (Green Mountain Grill), so there is no
   `pip install gmg` and no `MATPredict[gmg]` extra yet**; GmG's wheel will be published as **`guess-my-genome`** (import name `gmg`; GmG PR #136), and
   until then is built from the GmG repository (`maturin`). The import is accepted only if the module has `Index`.
2. **`gmg` binary on PATH** (`gmg run GENOME --index IDX -o DIR --no-report`, reading `results.jsonl`). GmG releases ship
   linux-64 and osx-arm64 binaries, so this is the practical route today, including the HPCC.
3. Neither found: clear error that says how to install, never a silent fallback to a different route.

## The index: fetch it from GmG, cache it

Yes, MATPredict can use the index GmG publishes. Checked on 2026-11-02:

- GmG releases `index-YYYY.M.N` carry `bfd-reps-<ver>-s20000.gmgi` (67 MB, 7,479 genomes) and `bfd-<ver>-s10000.gmgi`
  (213 MB, all 23,677), with `SHA256SUMS` and `MANIFEST.txt`.
- The GmG web site serves the same files from `https://gmg.fungalgenomes.org/index/` (`manifest.json` lists them) with
  `access-control-allow-origin: *`. The GitHub release URL redirects to `release-assets.githubusercontent.com` and showed
  no CORS header in my check, so a browser page must use the site URL, and the CLI can use either.

Proposal for the CLI:

- Default index: the representatives index (67 MB, enough for a phylum call; the full index only adds species resolution).
- `--gmg-index` accepts a path, a release file name (pinned), `latest`, or a URL; see Decisions 1 and 2. `latest` reads `manifest.json` (or the newest `index-*` release
  with `gh`/the GitHub API when offline-from-the-site), downloads once into `~/.cache/matpredict/gmg/` and verifies the sha256 against the release's `SHA256SUMS`.
- The index version in use is written to the report. `latest` can change the answer between runs, so batch and panel runs
  should pin it (`--gmg-index bfd-reps-2026.10.5-s20000.gmgi` or an env var) and say so in the run record.
- No download inside tests or CI by default; tests use a tiny index built with `gmg build` from the repo's fixtures.
- HPCC: compute nodes may not reach the internet, so document `matpredict fetch-gmg-index` (login node) and the shared
  path, the way the taxonomy table is handled today.

## Browser (WASM) step, after the CLI works

- The GmG WASM `Index` (`crates/gmg-wasm`) takes the index bytes and a genome as bytes, and returns the JSON record.
  Index loading, progress and Cache API caching already exist in GmG's `web/worker.js`; reuse that rather than rewrite.
- The MATPredict page runs Python in Pyodide on the main thread. The guess does not need Python: a small JS module
  loads the GmG wasm package in a Web Worker, sketches the user's genome, and returns the phylum. It then sets
  `#phylum` to the matching option, or leaves it and shows "no confident call".
- Same decision table as the CLI, implemented once in JS (`web/spike/page/gmg_route.mjs`) and unit-tested against the same
  cases as the Python function, so the two cannot drift. The mapping from GmG phylum to the page's options is a tiny table
  read from the page's own `<option>` values.
- Index: fetched from `https://gmg.fungalgenomes.org/index/` (CORS is open), default the 67 MB representatives index,
  cached by the browser. The page says up front that it is a one-time 67 MB download and offers "Guess" only on request,
  so nobody downloads it by accident.
- Pin the GmG wasm package and index versions in the page build, as the Pyodide version is pinned today.
- Memory: GmG measured about 180 MB in the browser for this index. The Pyodide genome and the MATPredict tools are also
  in memory; the combined peak is not measured yet. Run the guess first, then free the GmG index before starting detect.

## Genetic code

The translation table is a property of the taxon and MATPredict reads it from the taxonomy for `--taxid` (CTG-clade yeasts use
table 12). With `--phylum auto`, the table is read the same way for the family or order taxid GmG called (`default_genetic_code`),
unless `--genetic-code` or `--genetic-code-map` gives one. A failed lookup is recorded in `parameters.phylum_guess.genetic_code_error`
and the run uses table 1. A phylum-only route (no taxid) uses table 1, as `--phylum NAME` does today.

## Tests

- Unit: the decision table, with a fake GmG record per row (`ok`, `low_confidence`, `no_match`, `not_fungal`, a fungal
  phylum with no families, an `error`), using the injected-resolver style the `route()` tests already use.
- Contract: for every directory in `db/`, the NCBI name GmG returns for a genome of that phylum equals the directory name
  (guards against a naming mismatch, for example after an NCBI rename).
- Integration (marked slow, needs an index): the MATPredict test genomes (`testset/`, the 23 Zygo genomes) give the
  phylum curated for them, and a non-fungal genome is not searched.
- Report: `phylum_guess` is present, and a user-given `--phylum` is not overwritten.

## Decisions (2026-11-02)

1. **Cache location:** `~/.cache/matpredict/gmg/` (honouring `XDG_CACHE_HOME`), not the repo's `.matpredict_cache`.
   `MATPREDICT_CACHE_DIR` still overrides it.
2. **Index version:** production and batch runs pin an index. The pin is `--gmg-index <file>`, `$MATPREDICT_GMG_INDEX`, or a
   default written into MATPredict's release. When a newer `index-*` release exists, print one warning line (never change the
   answer silently, never fail). `--gmg-index latest` stays available for exploration and is recorded as such in the report.
3. **GmG routing hint:** GmG adds a stable field to its JSON record (and the `results.tsv` / Python / WASM outputs, which
   share the code in `gmg-core`): `routing_hint` = `{phylum, phylum_taxid, deepest_rank, deepest_taxid, usable, reason}`.
   `usable` is true only when the phylum is called with GmG's own thresholds and the input is not flagged `not_fungal`.
   MATPredict reads that field and keeps only the table "GmG phylum -> curated family directory" itself. This is a
   schema bump in GmG (documented in `docs/results-schema.md`, covered by the parity tests) and goes into a GmG release first.

## Measured: GmG phylum calls on the regression panel (2026-11-02)

165 of the 166 genomes of `testset/regression_panel.tsv` (29 Ascomycota, 33 Basidiomycota, 80 Mucoromycota and 23 `record_sources` genomes of
other lineages; one more Ascomycota genome is not in BFD). GmG 2026.10.5 full index, binary ca8e3b9, HPCC job 31607271 (`GmG_runs/matpredict_phylum/`). Truth is the NCBI phylum of the
genome's own BFD entry. `usable` is the proposed `routing_hint` rule (phylum called, not `not_fungal`, not `weak_phylum`). All of these genomes
are in BFD, so each was run with its own entry, its species, or its genus hidden from the reference.

| reference | usable | right and usable | wrong and usable | called but blocked (right / wrong) | no call |
|---|---|---|---|---|---|
| genome in the reference (upper bound) | 165 | 165 | 0 | 0 | 0 |
| own entry hidden | 134 | 134 | 0 | 21 / 8 | 2 |
| species hidden | 106 | 106 | 0 | 41 / 13 (in a curated phylum: 41 / 5) | 5 |
| genus hidden (genome of a new genus) | 40 | 36 | 4 | 81 in a curated phylum: 70 right, 11 wrong; 29 other, all wrong | 15 |

What this says:

- When the genome has relatives in the reference (own entry or species hidden), no usable call was wrong (0 of 240), and 64-81% of the
  genomes get one.
- A genome of a genus that is not in the reference is hard: only 40 of 165 get a usable call, and 4 of those are wrong. All 4 were called
  Neocallimastigomycota (3 *Cunninghamella*, flagged `mixture_suspected`; 1 *Pneumocystis* with no flag, best-hit ANI 0.852, just
  above the 0.85 weak cutoff). Neocallimastigomycota has no curated family, so the consequence is `not_searched` with a reason, never a search of
  the wrong family. No usable call named a wrong curated phylum.
- Most blocked calls in a curated phylum are right (species hidden 41 of 46, genus hidden 70 of 81). The `weak_phylum` block is therefore
  conservative for MATPredict: blocking costs a missed locus, where a wrong call costs a wasted search. This is a policy choice for
  you, see open question 5.
- Mucoromycota is the weak spot (few relatives in BFD): with the genus hidden 27 of 80 genomes get a wrong phylum (24 of them blocked) and 7 get none.
- The deepest taxid is not safe to pass blindly. With the species hidden, 12 of 12 species-level calls carry a taxid that is not the
  genome's own, and with the genus hidden every genus-level call (16) does, because the right genus is absent. Family and order taxids are
  right in 52 of 60 usable calls across the two runs (37 of 38 with the species hidden, 15 of 22 with the genus hidden). So pass the **order or family** taxid, not the genus or species one, to `route()`.
- Caveat: hiding a whole genus is harsher than real use, where a new MATPredict genome usually has a congener among the 23,669 BFD genomes.
  The "own entry hidden" and "species hidden" rows are closer to the typical case; the genus-hidden row is the worst case.

## Measured: `detect --phylum auto` end to end on the HPCC (2026-11-03)

30 genomes of the regression panel (10 each of Ascomycota, Basidiomycota and Mucoromycota), each run three ways with `--preset v0.8` on MATPredict
`main` (08022a2, with #233): `--phylum auto` (GmG v0.3.1, the representatives index 2026.10.5), `--taxid` (the panel taxid, the published v0.8.0
routing) and `--phylum <true phylum>`. The genomes are BFD genomes but are not in the representatives index. All 90 runs finished without error. Taxonomy:
NCBI snapshot 2026-11-01. Raw results: `GmG_runs/matpredict_phylum/` on the HPCC (`summary.tsv`, `loci.tsv`, `detect_out/`).

| | `auto` | `--taxid` | `--phylum` |
|---|---|---|---|
| loci found | 22 | 27 | 27 |
| families searched per genome (mean) | 1.4 | 2.0 | 6.0 |
| seconds per genome (mean) | 72 | 80 | 163 |

- **`auto` finds what `--taxid` finds, except where GmG declines.** On 27 of 30 genomes `auto` and `--taxid` agree on every locus (25 routed by lineage, 2 by the phylum
  fallback, as the taxid route also did). `auto` found nothing the taxid route did not. Its 5 missed loci come from 3 genomes it did not search: *Rhizoctonia solani*
  (weak call, 2 loci) and `Zycmex1` (weak call, 1 locus), both with the right phylum, and `CC3` (no phylum call, 2 loci).
- **`--phylum auto-weak` recovers the two weak ones.** Re-run on those two genomes, it found the same loci as `--taxid` (same coordinates, confidence and idiomorph;
  3 of 3). `CC3` has no phylum call and stays not searched. So the weak path is safe on this sample, and it stays opt-in (decision 5 below).
- **Lineage routing is narrower than phylum routing, and that is not new.** `--taxid` and `--phylum` differ by 5 loci in each direction. Plain phylum routing found two
  medium-confidence loci with an undetermined idiomorph that both lineage routings missed (a *Cryptococcus* PR locus and a MATsc locus in `GCA_002924885.1`); elsewhere the same
  locus carries a different family label (`MAT` vs `MATtub`, `redPR`/`redHD` vs `MAT`/`Aalpha`). `auto` inherits exactly the `--taxid` behaviour; the cases belong with #228.
- **Genetic code:** all 27 searched genomes got table 1 (none is a CUG-clade yeast; table 12 was seen only for *C. albicans* in the earlier routing probe).
- **Limits:** 30 genomes, the representatives index only; the hard case (a genus absent from the index) is covered by the panel table above, not by this run.

## Open questions

3. (Answered above for the regression panel: pass the phylum plus an order or family taxid. Not yet run: the 23 Zygo genomes and
   the chytrid negative control, which are outside this panel.)
4. Reads input: GmG classifies reads too, and MATPredict has a reads path (`src/MATPredict/reads`). Out of scope here, but
   the same function could serve it.

5. Policy for a called-but-blocked phylum (weak call) in a curated phylum. Options: (a) `not_searched` (as designed, safe, misses
   loci); (b) search that phylum and mark the run `phylum_source: gmg_weak` (about 86% right on the genus-hidden panel); (c) (b)
   only behind an explicit `--phylum auto-weak`. **Decided 2026-11-03: (c)**, opt-in; the HPCC run above supports it (3 of 3 weak-call loci recovered).
6. Should `mixture_suspected` also make `usable` false? Three of the four wrong usable calls had it, but they are one genus
   (*Cunninghamella*), so the evidence is thin; check it on a wider set before changing the rule.

## Order of work

1. This note, reviewed.
2. GmG: add `routing_hint` (decision 3) and release it (v0.3.0, since the schema changes).
3. Run GmG on MATPredict's test genomes; record phylum, taxid and status (answers open question 3).
4. `detect/gmg_guess.py` + `--phylum auto` + report block + tests, with the fetch and cache of the index.
5. Batch support and docs; release MATPredict.
6. WASM: JS route module, worker, `#phylum` Guess button, page tests (headless, as GmG does).
