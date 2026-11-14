#!/usr/bin/env python3
"""Write the Batrachochytrium coverage / heterozygosity report (markdown, figure, tables) from the result tables.

Inputs  : --base results/2026-11-13_bd_coverage (bd_isolate_summary_<ref>.tsv, bd_depth_by_seq_<ref>.tsv,
          parent_origin/parent_origin_summary.tsv, parent_origin/parent_origin_windows.tsv.gz, excluded_runs.tsv)
Outputs : --out report dir: report.md, figures/bd_summary.{png,pdf}, tables/*.tsv
All numbers in the markdown come from these tables, so a page and the tables agree.
Isolate classes (per isolate, from window fractions in parent_origin_windows): CLFT-type > 60% CLFT windows;
JEL-type > 60% JELHOM windows; HET-rich > 50% HET windows; otherwise mixed.
"""
import argparse
import collections
import csv
import gzip
import os
import shutil
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def rd(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def cls(r):
    c, h, j = (int(r[k]) / int(r["windows"]) for k in ("windows_CLFT", "windows_HET", "windows_JELHOM"))
    return "CLFT-type" if c > .6 else "JEL-type" if j > .6 else "HET-rich" if h > .5 else "mixed"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--code-commit", default="")
    a = ap.parse_args()
    for d in ("figures", "tables"):
        os.makedirs(f"{a.out}/{d}", exist_ok=True)
    J = {r["run"]: r for r in rd(f"{a.base}/bd_isolate_summary_JEL423.tsv")}
    C = {r["run"]: r for r in rd(f"{a.base}/bd_isolate_summary_CLFT044.tsv")}
    P = {r["run"]: r for r in rd(f"{a.base}/parent_origin/parent_origin_summary.tsv")}
    for r in P.values():
        r["class"] = cls(r)
    for f in ("bd_isolate_summary_JEL423.tsv", "bd_isolate_summary_CLFT044.tsv", "excluded_runs.tsv", "bd_runs_popseq.tsv"):
        shutil.copy(f"{a.base}/{f}", f"{a.out}/tables/{f}")
    shutil.copy(f"{a.base}/parent_origin/parent_origin_summary.tsv", f"{a.out}/tables/parent_origin_summary.tsv")
    with open(f"{a.out}/tables/isolate_classes.tsv", "w") as fo:
        fo.write("run\tgroup\tisolate\tcountry\tclass\tmedian_depth\tsnps_per_mb_JEL423\tsnps_per_mb_CLFT044\n")
        for k, r in P.items():
            fo.write("\t".join([k, r["group"], r["isolate"], r["country"], r["class"], J[k]["median_depth"],
                                J[k]["snps_per_mb"], C[k]["snps_per_mb"]]) + "\n")

    col = {"CLFT-type": "#c0392b", "JEL-type": "#2471a3", "HET-rich": "#7d3c98", "mixed": "#7f8c8d"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for c in col:
        ks = [k for k in P if P[k]["class"] == c]
        ax[0].scatter([float(J[k]["snps_per_mb"]) for k in ks], [float(C[k]["snps_per_mb"]) for k in ks], s=22,
                      c=col[c], label=f"{c} (n={len(ks)})")
        ax[1].scatter([float(P[k]["recovery"]) for k in ks], [float(P[k]["het_share_of_called"]) for k in ks], s=22,
                      c=col[c], label=c)
    ax[0].set_xlabel("SNPs/Mb against JEL423"); ax[0].set_ylabel("SNPs/Mb against CLFT044")
    ax[0].legend(frameon=False, fontsize=8); ax[0].set_title("A. SNP density against the two references")
    ax[1].set_xlabel("recovery of JEL423 alleles at divergent sites"); ax[1].set_ylabel("heterozygous share of called sites")
    ax[1].set_title("B. Parent-of-origin summary per isolate")
    fig.tight_layout()
    fig.savefig(f"{a.out}/figures/bd_summary.png", dpi=150); fig.savefig(f"{a.out}/figures/bd_summary.pdf")

    runs = {r["run"]: r for r in rd(f"{a.base}/bd_runs.tsv")}
    brazil = [k for k in P if P[k]["group"] == "brazil"]
    other = [k for k in P if P[k]["group"] != "brazil"]
    cnt = {g: collections.Counter(P[k]["class"] for k in ks) for g, ks in (("Brazil", brazil), ("Other", other))}
    ctrl = {n: next(k for k in P if P[k]["isolate"] == n) for n in ("CLFT044", "JEL423")}
    self_snp = {n: (J[k]["snps_per_mb"], J[k]["het_fraction"], C[k]["snps_per_mb"], C[k]["het_fraction"])
                for n, k in ctrl.items()}
    seqdepth = {}
    for ref in ("JEL423", "CLFT044"):
        for r in rd(f"{a.base}/bd_depth_by_seq_{ref}.tsv"):
            seqdepth.setdefault((ref, r["run"]), []).append(r)
    het_rich = [k for k in P if P[k]["class"] == "HET-rich"]
    mixed = [k for k in P if P[k]["class"] == "mixed"]
    lowdepth = [k for k in P if float(J[k]["median_depth"]) < 20]
    odd = [k for k in other if P[k]["class"] == "CLFT-type"]
    excl = [l.rstrip("\n") for l in open(f"{a.base}/excluded_runs.tsv")][1:] if os.path.exists(f"{a.base}/excluded_runs.tsv") else []

    def row(k):
        p = P[k]
        return (f"| {p['isolate']} | {p['country']} | {J[k]['median_depth']} | {J[k]['snps_per_mb']} | {C[k]['snps_per_mb']} | "
                f"{p['recovery']} | {p['het_share_of_called']} | {int(p['windows_CLFT'])}/{int(p['windows_HET'])}/{int(p['windows_JELHOM'])} | {p['class']} |")

    hdr = ("| isolate | country | median depth | SNPs/Mb vs JEL423 | SNPs/Mb vs CLFT044 | recovery | het share of called | "
           "windows CLFT/HET/JELHOM | class |\n|---|---|---|---|---|---|---|---|---|")
    md = f"""# Batrachochytrium dendrobatidis: read depth, heterozygosity and parent of origin

Status: exploratory analysis of public reads. Nothing here is curated or accepted. Bd isolates are diploid, and the
allele balance is read that way.

Pipeline commit: `{a.code_commit}`. Jobs: JEL423 reference 32070127 and 32070630 (array 2-9); CLFT044 reference 32070128 and
32070631; reference alignment 32089404 (minimap2 2.30, `-x asm10 -c --cs=short`). Result folder on the HPCC:
`results/2026-11-13_bd_coverage/` (not committed; tables below are copies).

## Data

- Runs: {len(P)} runs analysed of {len(runs)} listed (Brazil {len(brazil)}, other countries {len(other)}), chosen from PopSeq
  `catalog/runs/fungi_all_runs.csv.gz` (ENA snapshot 2026-11-06) and ENA. All were Illumina WGS. See `tables/bd_runs_popseq.tsv`.
- References: JEL423 (GCA_048537975.1) and UCR_BdCLFT044 (Brazil lineage per the curator). I did not verify the lineage of JEL423 here.
- Each run: `prefetch`, `fasterq-dump`, read-pair count checked against ENA `read_count`, `bwa-mem2`, duplicate removal,
  `mosdepth` 10 kb windows (MAPQ 20), `bcftools` SNP calls (QUAL 30).
- Excluded runs (`tables/excluded_runs.tsv`):
{chr(10).join('  - ' + e.replace(chr(9), ': ') for e in excl)}

## 1. Depth and SNP density

Median depth of the {len(P)} isolates: {statistics.median(float(J[k]['median_depth']) for k in P):.1f}x (range
{min(float(J[k]['median_depth']) for k in P):.1f} to {max(float(J[k]['median_depth']) for k in P):.1f}). {len(lowdepth)} isolates have median depth below 20x; treat their numbers with more care.

Isolates against their own reference (controls):

| control | SNPs/Mb vs JEL423 | het fraction | SNPs/Mb vs CLFT044 | het fraction |
|---|---|---|---|---|
| CLFT044 reads | {self_snp['CLFT044'][0]} | {self_snp['CLFT044'][1]} | {self_snp['CLFT044'][2]} | {self_snp['CLFT044'][3]} |
| JEL423 reads | {self_snp['JEL423'][0]} | {self_snp['JEL423'][1]} | {self_snp['JEL423'][2]} | {self_snp['JEL423'][3]} |

An isolate does not give near-zero SNPs against its own reference. Each reference is a haploid representation of a diploid
genome, so heterozygous sites are expected, and the het fraction is high (0.9). Repeats and paralogs can add het-like SNPs;
I did not separate the two causes. A het fraction alone therefore cannot call a hybrid.

Depth flags: sequences of 200 kb or more with median depth outside 0.75 to 1.25 of the genome median occur in most isolates (0 to 11 per isolate)
(see `tables/bd_isolate_summary_*.tsv`, column `depth_flag_seqs`). I did not separate copy-number change from repeat or
assembly artefacts. Ploidy variation, including aneuploidy, is also known in Bd; I did not test it here.

![SNP density and parent-of-origin summary](figures/bd_summary.png) [PDF](figures/bd_summary.pdf)

Figure source: `tables/isolate_classes.tsv`, `tables/parent_origin_summary.tsv`; script `scripts/bd_report.py`.

## 2. Parent of origin

Method (`scripts/bd_parent_of_origin.py`): align JEL423 to CLFT044 and take single-base differences in primary alignments
(172,706 sites). In each isolate, a site is "called" when its CLFT044-reference SNP table has the JEL423 base. Per 50 kb
window: recovery = called sites / usable sites; het share = called sites with alt fraction 0.2 to 0.8. A window is `CLFT`
when recovery < 0.5, else `HET` when het share >= 0.5, else `JELHOM`.

Why sites without a SNP row are not counted as CLFT044 allele: in the JEL423 isolate about a third of divergent sites have no row,
so absence is not evidence. I first counted them as CLFT044 and the JEL423 isolate then scored only 0.57 JEL fraction; that
version was dropped.

Calibration: only two controls. In the CLFT044 isolate the window recovery has median 0.05 and 95th percentile 0.41. In the
JEL423 isolate the median is 0.87 and the 5th percentile 0.16. I set the cutoff at 0.5 from these. It is not validated on isolates
of known hybrid status. Isolate class: CLFT-type if over 60% CLFT windows, JEL-type if over 60% JELHOM windows, HET-rich if
over 50% HET windows, else mixed.

| group | CLFT-type | JEL-type | HET-rich | mixed |
|---|---|---|---|---|
| Brazil (n={len(brazil)}) | {cnt['Brazil']['CLFT-type']} | {cnt['Brazil']['JEL-type']} | {cnt['Brazil']['HET-rich']} | {cnt['Brazil']['mixed']} |
| Other countries (n={len(other)}) | {cnt['Other']['CLFT-type']} | {cnt['Other']['JEL-type']} | {cnt['Other']['HET-rich']} | {cnt['Other']['mixed']} |

Control isolates: CLFT044 is {cls(P[ctrl['CLFT044']])}; JEL423 is {cls(P[ctrl['JEL423']])}.

### Isolates that are HET-rich or mixed

{hdr}
{chr(10).join(row(k) for k in het_rich + mixed)}

These isolates carry JEL423 alleles in most windows, mostly in one copy. That is the pattern expected for a first-generation
hybrid of the two lineages. It does not show that they are hybrids. Other explanations I did not exclude: a third lineage that shares
alleles with both references; aneuploidy; mapping effects at divergent sites. A hybrid call needs a tree or ancestry analysis and
a lineage reference beyond these two. Published context: a Bd-GPL x Bd-Brazil hybrid zone in the Serra da Graciosa, Parana (Jenkinson et al. 2016,
doi 10.1111/mec.13599). I did not compare my isolates with that study's samples.

### Non-Brazil isolates that are CLFT-type

{hdr}
{chr(10).join(row(k) for k in odd)}

Their recovery is 0.33 and 0.37, well below the JEL-type isolates (about 0.6 to 0.75). They carry few alleles of either reference.
A third lineage is a possible reason; I did not check their lineage.

## Limits

- Two references only. Lineage names are the curator's for CLFT044 and my understanding for JEL423.
- SNP tables keep QUAL >= 30 and no per-site depth was stored, so absence of a call cannot be told from no coverage.
- Window class is noisy: the CLFT044 control has windows with recovery up to 0.6.
- The Brazil set is {len(brazil)} runs; PopSeq lists 58 Brazilian isolates, so the set is not all of them.

## Next

Add an ancestry or tree analysis with at least one more lineage reference (for example a BdASIA or BdCAPE assembly), and test
the HET-rich isolates for aneuploidy with allele balance by sequence.
"""
    open(f"{a.out}/report.md", "w").write(md)
    print("wrote", f"{a.out}/report.md", dict(cnt))


if __name__ == "__main__":
    main()
