# Handoff: Xylariales MAT models, NC1011 interval, synteny (2026-10-05)

For agents joining this work. Written mid-run; jobs below were still running.
Read `docs/superpowers/specs/2026-10-05-xylariales-models-design.md` and the Xylariales
section of `docs/notes/2026-09-24_mat-reference-gap-literature.md` first (both merged, PR #30).

## Where things are
- Worktree (on /bigdata, visible to compute nodes): `/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval`
  Branch `xylariales-interval`, created from origin/main after PRs #28-#30 merged. Nothing committed yet.
- All new work: `results/2026-10-05_xylariales_nc1011_interval/` (untracked). Scripts `01`-`05*` there, shared paths in `common.sh`.
- Do not edit scripts in the worktree while a job that uses them runs (the project lesson from 2026-09-21).
- Pixi env for python/miniprot/diamond: `/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin`.
  Slurm is not on PATH in plain ssh: `export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin:$PATH`.

## Question being worked
*Xylaria flabelliformis* NC1011 (GCA_022453505.1, = JGI Xylcub1) gets one MAT1-2/medium call on a 734-aa HMG
protein (KAI0192626.1, contig JAJLYR010000012.1:276,138-285,941) that is probably not a MAT gene. Its conserved
flank block is intact but has no MAT gene. The curator's hypothesis: an inversion or translocation changed the
gene order and the MAT gene is missing. Robinson & Natvig 2019 (PMID 30557613) found no canonical MAT in 35 Xylariales.

## Results so far
1. **Interval check (done, job 29402338, `work/`).** NC1011 contig JAJLYR010000004.1: SLA2 KAI0195601.1 (243,984-247,733, +, partial)
   - COX13 KAI0195602.1 (248,209-249,501, -) - APN2 KAI0195603.1 (249,620-252,216, +). Gaps SLA2-COX13 476 bp, COX13-APN2 119 bp,
   APN2-next gene 784 bp, previous gene-SLA2 330 bp. None of the 14 annotated proteins in 225-272 kb has an HMG-box, MATalpha_HMGbox
   or MAT1-1-2 domain (Pfam PF00505, PF09011, PF04769, PF17043; hmmsearch E 1e-3). None of 139 six-frame ORFs (>=150 nt) hits them
   (E 1e-2). tBLASTn of 105 database MAT proteins: best 28 bits (E 0.005), 20-67 aa fragments (real Xylariales MAT hits were 39-47 bits).
2. **Order of APC5, SLA2, COX13, APN2 (job 29402406, `synteny_order.tsv`, annotated genomes only, 4 query proteins).**
   Normalised so APC5 is +. Canonical (N. crassa, Podospora, Chaetomium, Colletotrichum): APC5(+) COX13(-) APN2(+) SLA2(-).
   Class A: APC5(+) SLA2(+) APN2(-) COX13(+) (M. bolleyi, Hypoxylon CI_4A, Daldinia EC12). Class B: APC5(+) SLA2(+) COX13(-) APN2(+)
   (NC1011, G536 = same species, Rosellinia necatrix, Eutypa lata).
   Reading (an inference, not tested): A = canonical block inverted (MAT site stays between SLA2 and APN2; M. bolleyi has an HMG gene there
   per the paper). B = A plus a second inversion of [MAT-site, APN2, COX13], which would move the MAT site to the far side of APN2
   (the 784-bp gap in NC1011). In the paper's table, Rosellinia and E. lata have no linked HMG gene. 14 of 31 genomes had no NCBI annotation
   (including our M. paspali, D. rigidum, XT01) and are unclassified. P. fici (145 kb span) and Hypoxylon EC38 look like paralog hits.
3. **Paper supplement parsed** (PR #30): `docs/notes/2026-10-05_robinson-natvig-2019-xylariales.tsv`. 16 of 35 genomes have a MATA_HMG gene
   linked to SLA2/APN2 (7 between, 9 adjacent); nearly all are closest to N. crassa NCU03481, a non-MAT regulator.

## Jobs (Slurm, partition batch; exfab had a ~2-day queue estimate on 2026-10-05, so avoid it unless the curator insists)
| job | name | what | output |
|---|---|---|---|
| 29402370 | xyl-rnaseq | map NC1011 RNA-seq SRR8861595 (hisat2), depth over 225-272 kb and spliced reads in 247-250.5 kb | `work/rna/` (BAM, `*_region_depth.tsv`, `*_spliced_reads_*.tsv`, `hisat2_summary.txt`) |
| 29402372 | xyl-hmgtree | HMG-box tree: 19 of 31 proteomes downloaded (12 had no NCBI proteins, see `tree/proteome_download.log`), hmmsearch, mafft L-INS-i, trimal -gt 0.3, RAxML-NG LG+G4, 200 bootstraps | `tree/rx.raxml.support` when done; domains in `tree/hmg_domains.tsv` |
| 29402670 | xyl-mpsyn | miniprot (pixi build 0.18-r281) of the 14 NC1011 neighbourhood proteins (KAI0195597.1-KAI0195610.1) against 414 genomes: all 257 Xylariales, 98 Amphisphaeriales, up to 10 genera in each of 7 outgroup orders; 150-kb window with most genes; also an HMG-like search in the SLA2-APN2 vicinity | `genomes_miniprot.tsv`, `genes_miniprot.tsv`, `adjacency_miniprot.tsv` |
Check: `squeue -u $USER | grep xyl`; logs in `/bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_*`.
Slurm script gotchas already fixed: source `common.sh` by absolute path (Slurm runs a spool copy); `source /etc/profile` returns non-zero so guard it with `set +eu ... set -eu`.

## Update (later 2026-10-05)
- xyl-mpsyn (29402670) and the breakpoint alignments (29403041) are DONE. Results are summarised in `analysis/2026-10-05_xylariales-synteny.md`
  (read that first; it supersedes result 2 above, which was a 4-gene, annotated-only first pass and mis-described state A as an inversion).
  Tables: `genomes_miniprot.tsv`, `adjacency_miniprot.tsv`, `genes_miniprot.tsv`, `breakpoints/summary.md` (in the results directory).
- xyl-rnaseq (29402370) is DONE: summary in the analysis note and `work/rna/coverage_summary.tsv` (script `07_rna_summary.py`). xyl-hmgtree (29402372) is DONE: placement in `tree/hmg_placement.tsv` (script `08_hmg_placement.py`), summarised in section 6 of the analysis note. KAI0192626.1 is nearest NCU03481 (1.02 vs 2.22 for the nearest MAT reference) at 14% bootstrap; the tree is under-powered (9 MAT references, not monophyletic).
- Items 1 and 2 below are done at gene level; nucleotide-level breakpoints are NOT resolved (minimap2 found almost no alignable blocks).

## Still to do (in order)
1. When xyl-mpsyn finishes: classify each genome (canonical / A / B / other / incomplete), cross-tabulate against order, genus and the HMG-in-vicinity column,
   and against MATPredict v0.6.0 calls (`results/2026-10-03_ascomycota_v060/loci.tsv`, `genomes.tsv`).
2. Breakpoints: (a) gene level, from `adjacency_miniprot.tsv` (which NC1011 adjacencies are conserved in outgroups, class A, class B; report the intergenic
   spans in NC1011 coordinates); (b) nucleotide level, by aligning (minimap2 2.30 module, or mummer 4) the windows of congeneric genomes that differ in class.
3. Place KAI0192626.1 and the other Xylariales HMG proteins in the finished tree (clades: MAT1-2-1/MAT1-1-3 references from the database, NCU03481, fmf-1);
   write the placement table. Note the tree has only 19 proteomes; NCU03481/fmf-1 identity is by header text lost in the sed rename, so re-identify them by reciprocal hit.
4. RNA-seq: summarise coverage per gene and per gap (247,733-248,209, 249,501-249,620, 252,216-253,000) and any spliced reads.
5. Deliverable requested by the curator: a markdown summary of the tables and conclusions in `analysis/` (suggested `analysis/2026-10-05_xylariales-synteny.md`).
6. Commit on branch `xylariales-interval` (small files only: scripts, summary tables, markdown; not BAMs, PAFs, proteomes, `tree/prot`, `miniprot/` raw output,
   `data/ncbi_dataset`, `synteny/*.zip|gff|faa`). Ask before pushing or opening a PR.

## Decisions waiting for the curator (from the design spec, PR #30)
- Withhold or downgrade HMG-only Xylariales calls without clade support now (NC1011 is the case)?
- Label for hits that place with NCU03481/fmf-1: `hmg_regulator_like` or silent?
- *Pestalotiopsis* (Amphisphaeriales in our taxonomy) in the Xylariales seed set or held out?
- "miniprot2": the curator wrote this; the pixi build (0.18-r281) was used and there is also a `miniprot/0.2` module. Confirm which was meant.

## Repository and HPCC state (other agents should not trip over these)
- Merged today: #23-#30 (frameshift-aware, Mycotypha, spec-reads-population/decisions, curate-foxysporum #29 with its regression, xylariales-models #30).
  F. oxysporum regression: no call changes in 212 genomes.
- The top-level checkout `/bigdata/stajichlab/jstajich/projects/MATPredict` is 275+ commits behind origin/main and has untracked `biomni/ leuco/ logs/ results/ testset/Zygo/`.
  `git pull` fails because untracked results files collide with tracked ones. The curator started `git stash -u` on gpu12; no stash entry had appeared when last
  checked. Do not run git commands that touch its index until that finishes; do not pull there.
- Worktree cleanup done: 40 merged, clean worktrees removed; kept: `curate-foxysporum`, `curate-mycotypha`, `curation-mucor-dothideo`, `polish-scope-cuts`,
  `mucoromycota-scale-testing` (untracked logs), `spec-*`, measurement branches (`run-algcand`, `run-merge-measure*`, `run-nf-measure`, `run-umb-noguard`),
  and about 65 detached `run-*` frozen snapshots that results notes refer to. Do not delete the `run-*` ones without asking.
- Unmerged work still on branches: `curate-mycotypha` (sexM alignment findings, plan for alignment refinement and placement tests).
- The running BFD annotation pipeline (`nf-*` jobs, `do_annotation_wave1`) belongs to a different project; leave it alone.

## Conventions
- Commit trailer: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`; PR bodies end with the Claude Code line.
- Results directories track README, jobs.tsv and summary markdown, not large per-locus TSVs (see `results/2026-10-04_regression_frameshift`).
- Pushing and opening PRs have each been done only after the curator said so; keep asking.
- `/scratch` is node-local; keep anything a Slurm job reads on /bigdata.
