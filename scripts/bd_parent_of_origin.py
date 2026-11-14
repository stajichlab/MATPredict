#!/usr/bin/env python3
"""Parent-of-origin scan of Batrachochytrium isolates in CLFT044 coordinates.

Idea: two reference assemblies (CLFT044 = Brazil lineage, JEL423 = reference for the other lineage) differ at
single-base sites. Along each sequence, the share of those sites where an isolate carries the JEL423 allele, and
whether it carries it in one or both copies, shows which lineage(s) the isolate resembles.

Diploid reading: Bd isolates are diploid. At a divergent site a homozygous isolate has f near 0 (CLFT044 allele) or
1 (JEL423 allele); a heterozygous one has f near 0.5. An F1 hybrid should be heterozygous at divergent sites
genome-wide. Later-generation hybrids give a mosaic of CLFT, HET and JELHOM windows. Aneuploid or changed-ploidy
sequences give other fractions (for example 1/3 and 2/3), which also fall in the 0.2-0.8 het range; compare with
bd_depth_by_seq_<ref>.tsv before reading HET as hybrid. Within-lineage heterozygosity is not removed.

Inputs
  --paf     minimap2 -c --cs=short alignment, query JEL423, target CLFT044 (JEL423_vs_CLFT044.paf.gz)
  --base    results folder with per_run_CLFT044/<run>.snps.tsv.gz and <run>.regions.bed.gz, and bd_runs.tsv
            (the SNP files come from scripts/bd_coverage_isolates.slurm run against CLFT044)
Method
  divergent sites   '*xy' cs substitutions of primary alignments with block >= --min-block; a CLFT044 position
                    seen in more than one alignment is dropped.
  site call         a divergent site is "called" when the isolate SNP table has a biallelic row whose alt equals the
                    JEL423 base; its fraction f = alt/(ref+alt). Sites with no SNP row are NOT counted as CLFT044
                    allele: on the JEL423 isolate itself about a third of divergent sites have no row (no reads
                    or no call), so absence is not evidence. Sites in 10 kb windows with depth < --min-depth, and
                    rows with another alt base or several alts, are skipped.
  window            --window bp in CLFT044 coordinates; kept when >= --min-sites usable divergent sites.
                    recovery r = called sites / usable sites; het share = called sites with 0.2 <= f <= 0.8.
                    class: r < --min-recovery -> 'CLFT' (JEL allele rare); else het share >= --het-share -> 'HET'
                    (JEL and CLFT alleles both present); else 'JELHOM'.
                    The thresholds are set from the two control isolates (JEL423 and CLFT044) and are not
                    validated on isolates of known hybrid status.
Outputs (in --out): parent_origin_summary.tsv, parent_origin_windows.tsv.gz
Limits: the SNP table keeps QUAL >= 30 only; recovery on the JEL423 control is below 1, so r is a relative score.
"""
import argparse
import collections
import gzip
import os
import re

CS_TOKEN = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")


def divergent_sites(paf_path, min_block=5000):
    """Return {(target, pos1): query_base} for single-base substitutions in primary alignments."""
    cover = collections.Counter()
    sites = {}
    op = gzip.open if paf_path.endswith(".gz") else open
    with op(paf_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or int(f[10]) < min_block:
                continue
            tags = {t[:2]: t[5:] for t in f[12:] if len(t) > 5}
            if tags.get("tp") not in (None, "P"):
                continue
            cs = tags.get("cs")
            if not cs:
                continue
            tname, pos = f[5], int(f[7])  # 0-based target start
            for tok in CS_TOKEN.findall(cs):
                if tok[0] == ":":
                    pos += int(tok[1:])
                elif tok[0] == "*":
                    k = (tname, pos + 1)
                    cover[k] += 1
                    sites[k] = tok[2].upper()
                    pos += 1
                elif tok[0] == "-":
                    pos += len(tok) - 1
    return {k: v for k, v in sites.items() if cover[k] == 1}


def read_depth_windows(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, v = line.rstrip("\n").split("\t")[:4]
            d[(c, int(s) // 10000 * 10000)] = float(v)
    return d


def read_snps(path):
    """{(chrom,pos): (alt, ref_count, alt_count)} for biallelic SNP rows."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or "," in f[3]:
                continue
            try:
                ad = [int(x) for x in f[6].split(",")]
            except ValueError:
                continue
            if len(ad) == 2:
                out[(f[0], int(f[1]))] = (f[3], ad[0], ad[1])
    return out


def site_calls(sites, snps, depth, min_depth):
    """Yield (chrom, pos, f) per usable divergent site; f is None when the JEL base was not called."""
    for (c, p), jb in sites.items():
        if depth.get((c, (p - 1) // 10000 * 10000), 0) < min_depth:
            continue
        row = snps.get((c, p))
        if row is None:
            yield c, p, None
        elif row[0] == jb and (row[1] + row[2]) > 0:
            yield c, p, row[2] / (row[1] + row[2])
        # another alt base: skipped


def classify(rec, het, min_recovery, het_share):
    return "CLFT" if rec < min_recovery else "HET" if het >= het_share else "JELHOM"


def window_table(calls, window, min_sites, min_recovery=0.2, het_share=0.5):
    """calls: iterable of (chrom,pos,f|None) -> sorted list of
    (chrom,start,usable,called,recovery,het_share_of_called,class)."""
    acc = collections.defaultdict(list)
    for c, p, f in calls:
        acc[(c, (p - 1) // window * window)].append(f)
    rows = []
    for (c, s), v in acc.items():
        if len(v) < min_sites:
            continue
        called = [x for x in v if x is not None]
        rec = len(called) / len(v)
        het = sum(0.2 <= x <= 0.8 for x in called) / len(called) if called else 0.0
        rows.append((c, s, len(v), len(called), round(rec, 4), round(het, 4),
                     classify(rec, het, min_recovery, het_share)))
    return sorted(rows)


def switches(rows):
    """Class changes between neighbouring kept windows on the same sequence."""
    n, prev = 0, None
    for r in rows:
        if prev and prev[0] == r[0] and prev[1] != r[6]:
            n += 1
        prev = (r[0], r[6])
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", default="CLFT044")
    ap.add_argument("--window", type=int, default=50000)
    ap.add_argument("--min-sites", type=int, default=20)
    ap.add_argument("--min-depth", type=float, default=10.0)
    ap.add_argument("--min-block", type=int, default=5000)
    ap.add_argument("--min-recovery", type=float, default=0.5)
    ap.add_argument("--het-share", type=float, default=0.5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sites = divergent_sites(a.paf, a.min_block)
    print(len(sites), "divergent single-base sites")
    meta = {}
    for l in list(open(f"{a.base}/bd_runs.tsv"))[1:]:
        f = l.rstrip("\n").split("\t")
        meta[f[1]] = (f[0], f[3], f[4])
    cols = ["run", "group", "isolate", "country", "usable_sites", "called_sites", "recovery", "het_share_of_called",
            "windows", "windows_CLFT", "windows_HET", "windows_JELHOM", "switches"]
    with open(f"{a.out}/parent_origin_summary.tsv", "w") as fs, \
            gzip.open(f"{a.out}/parent_origin_windows.tsv.gz", "wt") as fw:
        fs.write("\t".join(cols) + "\n")
        fw.write("run\tseq\tstart\tusable\tcalled\trecovery\thet_share\tclass\n")
        n = 0
        for run, (grp, iso, ctry) in meta.items():
            sp = f"{a.base}/per_run_{a.ref}/{run}.snps.tsv.gz"
            dp = f"{a.base}/per_run_{a.ref}/{run}.regions.bed.gz"
            if not (os.path.exists(sp) and os.path.exists(dp)):
                continue
            cl = list(site_calls(sites, read_snps(sp), read_depth_windows(dp), a.min_depth))
            if not cl:
                continue
            called = [x[2] for x in cl if x[2] is not None]
            rows = window_table(cl, a.window, a.min_sites, a.min_recovery, a.het_share)
            cnt = collections.Counter(r[6] for r in rows)
            het = sum(0.2 <= v <= 0.8 for v in called) / max(1, len(called))
            fs.write("\t".join(map(str, [run, grp, iso, ctry, len(cl), len(called), round(len(called) / len(cl), 4),
                                         round(het, 4), len(rows), cnt["CLFT"], cnt["HET"], cnt["JELHOM"],
                                         switches(rows)])) + "\n")
            for r in rows:
                fw.write(run + "\t" + "\t".join(map(str, r)) + "\n")
            n += 1
    print(n, "isolates")


if __name__ == "__main__":
    main()
