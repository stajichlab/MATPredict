#!/usr/bin/env python3
"""Lineage-diagnostic divergent sites from several Batrachochytrium reference assemblies.

Read-only analysis. The 172,706 sites where JEL423 and CLFT044 differ (ref/divergent_sites.tsv: contig, 1-based position in CLFT044,
JEL423 base) may hold changes private to the two assemblies. With more assemblies of each lineage, a site is kept as lineage-diagnostic
only when every Brazil-lineage assembly carries the CLFT044 base and every other-lineage assembly carries the JEL423 base.

Inputs
  --sites   divergent_sites.tsv (CLFT044 coordinates, JEL423 base)
  --paf     NAME=PATH, repeatable: minimap2 -c --cs=short alignment, query = that assembly, target = CLFT044
  --brazil  names (from --paf) of Brazil-lineage assemblies expected to carry the CLFT044 base
  --other   names of other-lineage assemblies expected to carry the JEL423 base (JEL423 itself defines the site list)
Per assembly and site the state is
  jel    primary block >= --min-block covers the site and the assembly base equals the JEL423 base
  clft   a block covers the site and the assembly has no substitution there (CLFT044 base)
  third  a block covers the site and the assembly carries a third base
  na     no block, or more than one block, covers the site
Outputs (--out folder): site_states.tsv.gz (contig, pos, jel_base, one state column per assembly, diagnostic 0/1),
  diagnostic_sites.tsv (contig, pos, jel_base; same format as the input), summary.tsv (assembly, n_sites, jel, clft, third, na).
Limit: an assembly that is lower quality or has a gap at a site is 'na' and does not veto the site; a site is diagnostic when no
  assembly contradicts it and at least one assembly of each group supports it. Lineage membership of each assembly comes from the
  user, not from this script.
"""
import argparse
import bisect
import collections
import gzip
import os
import re
import sys

CS_TOKEN = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+)")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_paf(path, min_block=5000):
    """-> (subs {(contig, pos1): base}, blocks {contig: [(start1, end1), ...]}) from primary alignments with block >= min_block."""
    subs = {}
    blocks = collections.defaultdict(list)
    with open_text(path) as fh:
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
            tname, start, end = f[5], int(f[7]), int(f[8])
            blocks[tname].append((start + 1, end))
            pos = start
            for tok in CS_TOKEN.findall(cs):
                if tok[0] == ":":
                    pos += int(tok[1:])
                elif tok[0] == "*":
                    subs.setdefault((tname, pos + 1), []).append(tok[2].upper())
                    pos += 1
                elif tok[0] == "-":
                    pos += len(tok) - 1
    for v in blocks.values():
        v.sort()
    return subs, blocks


def n_covering(blocks, contig, pos):
    """Number of blocks covering a 1-based position (blocks are sorted by start)."""
    v = blocks.get(contig)
    if not v:
        return 0
    i = bisect.bisect_right(v, (pos, float("inf")))
    return sum(1 for s, e in v[max(0, i - 50):i] if s <= pos <= e)


def site_state(subs, blocks, contig, pos, jel_base):
    """'jel', 'clft', 'third' or 'na' for one site in one assembly."""
    if n_covering(blocks, contig, pos) != 1:
        return "na"
    s = subs.get((contig, pos))
    if not s:
        return "clft"
    if len(s) > 1:
        return "na"
    return "jel" if s[0] == jel_base.upper() else "third"


def is_diagnostic(states, brazil, other):
    """No assembly contradicts the site and each group has at least one supporting assembly."""
    for n in brazil:
        if states[n] in ("jel", "third"):
            return False
    for n in other:
        if states[n] in ("clft", "third"):
            return False
    return any(states[n] == "clft" for n in brazil) and any(states[n] == "jel" for n in other)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--paf", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--brazil", nargs="+", required=True)
    ap.add_argument("--other", nargs="+", required=True)
    ap.add_argument("--min-block", type=int, default=5000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pafs = dict(p.split("=", 1) for p in a.paf)
    for n in a.brazil + a.other:
        if n not in pafs:
            sys.exit(f"no --paf for {n}")
    parsed = {n: parse_paf(p, a.min_block) for n, p in pafs.items()}
    names = list(pafs)
    tally = {n: collections.Counter() for n in names}
    n_sites = n_diag = 0
    with gzip.open(os.path.join(a.out, "site_states.tsv.gz"), "wt") as fs, \
            open(os.path.join(a.out, "diagnostic_sites.tsv"), "w") as fd:
        fs.write("contig\tpos\tjel_base\t" + "\t".join(names) + "\tdiagnostic\n")
        for line in open_text(a.sites):
            c, p, b = line.rstrip("\n").split("\t")[:3]
            p = int(p)
            st = {n: site_state(parsed[n][0], parsed[n][1], c, p, b) for n in names}
            for n in names:
                tally[n][st[n]] += 1
            d = is_diagnostic(st, a.brazil, a.other)
            n_sites += 1
            n_diag += d
            fs.write(f"{c}\t{p}\t{b}\t" + "\t".join(st[n] for n in names) + f"\t{int(d)}\n")
            if d:
                fd.write(f"{c}\t{p}\t{b}\n")
    with open(os.path.join(a.out, "summary.tsv"), "w") as fo:
        fo.write("assembly\tn_sites\tjel\tclft\tthird\tna\n")
        for n in names:
            t = tally[n]
            fo.write(f"{n}\t{n_sites}\t{t['jel']}\t{t['clft']}\t{t['third']}\t{t['na']}\n")
        fo.write(f"diagnostic\t{n_sites}\t{n_diag}\t\t\t\n")
    print(open(os.path.join(a.out, "summary.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
