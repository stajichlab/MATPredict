#!/usr/bin/env python3
"""Chance expectation for adjacent HMG-box / homeodomain gene pairs per lineage.

Read-only discovery tool. Input is stage-1 output of flank_first_scan.py (genes_all.tsv and candidates.tsv per genome).
For each lineage (a list of genome names) and for each pair type (HMG+Homeo, 2xHomeo, 2xHMG) it counts the genomes that
have a pair within --window genes on one contig, and compares that with random placement: in every genome the same
number of HMG genes and homeodomain genes is drawn at random from all genes (--perms draws), the pair is looked for
again, and the lineage count is summed. The reported p is the share of draws with a count at least as large as observed.
A homeodomain family that is tandem-duplicated gives an excess here without being a mating-type locus, so the table also
gives the share of observed pairs on the same strand.

Usage: adjacent_pair_null.py --stage1 DIR --out OUT.tsv LIST[:LABEL] [LIST[:LABEL] ...]
"""
import argparse
import collections
import os
import random


def table(path):
    return [l.rstrip("\n").split("\t") for l in open(path)][1:]


def has_pair(a_items, b_items, window, same_class):
    """a_items, b_items: lists of (contig, index). True when some a and b (different genes) lie within window."""
    for c1, i1 in a_items:
        for c2, i2 in b_items:
            if c1 == c2 and 0 < abs(i1 - i2) <= window:
                return True
    return False


def pair_strands(a_items, b_items, window):
    """Return list of 'same'/'opposite' for observed pairs (a, b with a index < b index)."""
    out = []
    for c1, i1, s1 in a_items:
        for c2, i2, s2 in b_items:
            if c1 == c2 and 0 < i2 - i1 <= window:
                out.append("same" if s1 == s2 else "opposite")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage1", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--perms", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("lists", nargs="+")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    D = a.stage1
    out = open(a.out, "w")
    out.write("lineage\tpair_type\tgenomes\tobserved_genomes\texpected_mean\texpected_max\tp\tshare_same_strand\n")
    for spec in a.lists:
        path, _, label = spec.partition(":")
        label = label or os.path.basename(path)
        gs = [x.strip() for x in open(path) if x.strip()]
        data = []
        for g in gs:
            gf, cf = f"{D}/{g}.genes_all.tsv", f"{D}/{g}.candidates.tsv"
            if not (os.path.exists(gf) and os.path.exists(cf)):
                continue
            genes = table(gf)
            st = {r[0]: (r[1], int(r[5]), r[4]) for r in genes}
            hmg = [st[r[0]] for r in table(cf) if r[3] in ("HMG_box", "HMG_box_2") and r[0] in st]
            hd = [st[r[0]] for r in table(cf) if r[3] == "Homeodomain" and r[0] in st]
            data.append((hmg, hd, [(x[0], x[1]) for x in st.values()]))
        for ptype, (ka, kb) in {"HMG+Homeo": ("hmg", "hd"), "2xHomeo": ("hd", "hd"), "2xHMG": ("hmg", "hmg")}.items():
            obs, strands, exp = 0, [], [0] * a.perms
            for hmg, hd, allg in data:
                A, B = {"hmg": hmg, "hd": hd}[ka], {"hmg": hmg, "hd": hd}[kb]
                if not A or not B or (ka == kb and len(A) < 2):
                    continue
                ai = [(x[0], x[1]) for x in A]
                bi = [(x[0], x[1]) for x in B]
                if has_pair(ai, bi, a.window, ka == kb):
                    obs += 1
                if ka == kb:
                    strands += pair_strands(A, A, a.window)
                else:
                    strands += pair_strands(A, B, a.window) + pair_strands(B, A, a.window)
                for p in range(a.perms):
                    if ka == kb:
                        s = rng.sample(allg, len(A))
                        if has_pair(s, s, a.window, True):
                            exp[p] += 1
                    else:
                        s = rng.sample(allg, len(A) + len(B))
                        if has_pair(s[:len(A)], s[len(A):], a.window, False):
                            exp[p] += 1
            p_val = (sum(e >= obs for e in exp) + 1) / (a.perms + 1)
            same = round(strands.count("same") / len(strands), 3) if strands else "NA"
            out.write(f"{label}\t{ptype}\t{len(data)}\t{obs}\t{sum(exp) / a.perms:.1f}\t{max(exp)}\t{p_val:.4f}\t{same}\n")
    out.close()
    print(open(a.out).read())


if __name__ == "__main__":
    main()
