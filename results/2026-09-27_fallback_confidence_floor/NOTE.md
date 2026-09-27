# Confidence floor for phylum_fallback calls: does core identity separate good calls?

2026-09-27. Read-only analysis of existing reports. Scripts: `collect.py`
(-> `calls.tsv`, 8,429 calls), `tally.py` (-> `tally.txt`), `adjacency.py`
(-> `adjacency_<window>.tsv`, `adjacency_summary_<window>.txt`),
`find_hd_anchors.slurm` (job 29118608, MIP1/beta-fg top tblastn hits for the
144 genomes with a fallback HD call, `hd_hits/`). Breakdowns:
`hd_by_class.txt`, `asco_by_family.txt`, `affected_high_below40.tsv`.

## 1. How many calls a floor would touch (best core identity, any status)

phylum_fallback HIGH calls below each floor:

| run | high n | <30 | <35 | <40 |
|---|---:|---:|---:|---:|
| basidio_full | 268 | 6 | 37 | 100 |
| basidio_capoff | 2 | 1 | 1 | 1 |
| pilots_0924 | 24 | 0 | 0 | 2 |
| polishcap_cap6 | 24 | 0 | 0 | 3 |
| early_diverging | 2 | 0 | 1 | 1 |

Modelled-only identity moves the Ascomycota counts to 4 (<35) and 7 (<40).
The 40% floor in basidio_full hits HD Hymenochaetales 30, bLocus Tilletiales
22, HD Sporidiobolales 14, HD Corticiales 9, others <=5 (`tally.txt`).
explicit_phylum (Mortierellomycota/Kickxellomycota) has 0 high calls.

## 2-3. Synteny by identity band (anchor top hit within 20 kb of the best core gene)

Measured from the best-identity core gene, not the call span (the span
contains SLA2/APN2 when they are roster flanks). All testable fallback calls
were tested, not a sample. Stable at 10 and 50 kb.

| group | <30 | 30-35 | 35-40 | 40-50 | >=50 |
|---|---|---|---|---|---|
| Ascomycota (SLA2/APN2) | 32/32 | 43/45 | 53/75 | 101/129 | 67/81 |
| Agaricomycetes HD (MIP1/beta-fg) | – | 10/17 | 32/38 | 37/50 | 24/38 |
| Microbotryomycetes HD (MIP1) | 0/13 | 4/9 | 0/11 | 7/8 | 0/1 |

- Ascomycota: low-identity calls are adjacent MORE often. The failures sit at
  35-50% and above, concentrated in MATsc in Phaffomycetales, Pichiales and
  Ascoideales, and PM (0/10) - lineages where the 2026-09-24 test already
  showed SLA2 linkage is weak (Pichiomycetes 61%).
- Agaricomycetes HD: support is 59-84% in every band; no drop below 40%.
  MIP1 adjacency held in 12/18 Agaricomycotina in the anchor pilot, so ~70% is
  near the test's ceiling.
- Microbotryomycetes: MIP1 is not an anchor in Pucciniomycotina (0/7 in the
  anchor pilot), so these are not a valid test.
- Not testable (no anchor): bLocus (139), Tremellales-type MAT (21), MTL (52;
  flanks are roster genes), Mucoromycota (15), aLocus (2).

Of the 107 high fallback calls below 40%: 37 supported, 25 not supported (14 of
them Microbotryomycetes, invalid test), 45 not testable.

## Reading

Core identity does not separate supported from unsupported fallback calls
between 30% and 50%. A 40% floor would demote ~100 Basidiomycota high calls
whose synteny support (Agaricomycetes 42/55 below 40%) is at least as good as
that of the calls it keeps (61/88 at >=40%), plus 3-7 supported Ascomycota
calls. The only band with no support anywhere is <30% in Microbotryomycetes,
where the test itself is invalid. A 30% floor touches 7 high calls
(Sporidiobolales bLocus 5, Kriegeriales HD 1, Microbotryum bLocus 1), none
testable by synteny.
