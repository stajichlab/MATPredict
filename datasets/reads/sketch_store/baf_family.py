#!/usr/bin/env python3
"""The baf gene family (bafA, bafB, bafC, Pseudobaf) in a sketch store: build the screen targets, and compare the screen with assembly counts.

    baf_family.py targets  --hits baf-vs-Afum_genomes.FASTA.seqs_renamed --out baf_family.fa [--tsv baf_haplotypes.tsv]
    baf_family.py summarize --screen STORE/query/screen_baf_family/screen.tsv --calls starship_calls_wide.tsv
                            --counts strain_gene_counts.tab [strain_gene_counts_update.tab ...] --out baf_family_summary.json [--detail DIR]

targets: reads assembly hit sequences named gene|strain|locus (the pipeline's label bafD is bafA: 96% of bafA's 21-mers), collapses
exact-sequence haplotypes, and writes one record per gene, <gene>_all, with its haplotypes joined by 25 N (k-mers across a gap contain an N
and are skipped). One target per gene because the screen counts a hash as specific only if no other target has it: near-identical
haplotypes as separate targets would cancel each other. Run the screen with --max-solo 0 (allele-specific hashes appear without the
target's other hashes, which the solo-hash rule would drop):
    sketch_query.py --store STORE screen --name baf_family --targets baf_family.fa --max-solo 0

summarize: call a gene present when the screen's hit count >= CUT (bafA, bafC, Pseudobaf: 2; bafB has one specific hash: 1). Compares with
the assembly counts (strains matched by sample id local__<taxon> or normalised strain name; flagged samples left out), over a sweep of
cuts, against the earlier single-target bafA call, and gives prevalence over isolates and clone groups. The JSON holds counts only.
--detail DIR also writes concordance_<table>.tsv with strain names (local strains): keep it out of git."""
import argparse, collections, csv, json, os, re, sys

GENES = [("bafA", "bafD"), ("bafB", "bafB"), ("bafC", "bafC"), ("Pseudobaf", "bafPseudo")]   # (gene, column in the counts table)
CUT = {"bafA": 2, "bafB": 1, "bafC": 2, "Pseudobaf": 2}
LABEL_TO_GENE = {col: g for g, col in GENES}


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif name:
            buf.append(line)
    if name:
        yield name, "".join(buf).upper()


def cmd_targets(a):
    hap = collections.defaultdict(dict)                       # gene -> sequence -> set of strains (first-seen order kept)
    for name, seq in read_fasta(a.hits):
        label, strain = name.split("|")[:2]
        g = LABEL_TO_GENE.get(label)
        if g:
            hap[g].setdefault(seq, set()).add(strain)
    rows = []
    with open(a.out, "w") as fh:
        for g, _ in GENES:
            ordered = sorted(hap[g].items(), key=lambda x: -len(x[1]))
            fh.write(f">{g}_all\n{('N' * a.gap).join(s for s, _ in ordered)}\n")
            rows += [dict(gene=g, haplotype=f"{g}_hap{i}", strains=len(st), length=len(s)) for i, (s, st) in enumerate(ordered, 1)]
    if a.tsv:
        with open(a.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, ["gene", "haplotype", "strains", "length"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print({g: sum(1 for r in rows if r["gene"] == g) for g, _ in GENES}, "haplotypes per gene")


def rd(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(x):
    return re.sub(r"[^a-z0-9]", "", x.lower())


def cmd_summarize(a):
    scr, info = collections.defaultdict(dict), {}
    for r in rd(a.screen):
        t = r["target"]
        scr[r["sample_id"]][t] = (int(r["hits"]), float(r["rel_abund"])); info[t] = dict(hashes=int(r["target_hashes"]), specific=int(r["specific_hashes"]))
    calls = [r for r in rd(a.calls) if r["sample_id"] in scr]
    by_id = {r["sample_id"]: r for r in calls}
    by_strain = collections.defaultdict(list)
    for r in calls:
        if r["strain"]:
            by_strain[norm(r["strain"])].append(r)
    find = lambda taxon: [by_id["local__" + taxon]] if "local__" + taxon in by_id else by_strain.get(norm(taxon), [])
    hit = lambda s, g: scr[s["sample_id"]][f"{g}_all"][0]
    out = dict(screen=info, cuts=CUT, validation={}, population={})
    for path in a.counts:
        rows = rd(path)
        pairs = [(r, s) for r in rows for s in find(r["taxa"]) if not s["excluded"]]
        v = dict(taxa=len(rows), matched_taxa=len({r["taxa"] for r, _ in pairs}), pairs=len(pairs), genes={})
        for g, col in GENES:
            v["genes"][g] = {}
            for k in (1, 2, 3, 4):
                c = collections.Counter((int(r[col]) > 0, hit(s, g) >= k) for r, s in pairs)
                v["genes"][g][str(k)] = dict(tp=c[(True, True)], fn=c[(True, False)], fp=c[(False, True)], tn=c[(False, False)])
        c = collections.Counter((int(r["bafD"]) > 0, s["bafA_call"] == "P") for r, s in pairs)
        v["previous_single_target_bafA"] = dict(tp=c[(True, True)], fn=c[(True, False)], fp=c[(False, True)], tn=c[(False, False)])
        rel = collections.defaultdict(list)
        for r, s in pairs:
            if int(r["bafD"]) > 0 and hit(s, "bafA") >= CUT["bafA"]:
                rel[r["bafD"]].append(scr[s["sample_id"]]["bafA_all"][1])
        v["bafA_rel_abund_by_copy_number"] = {k: dict(n=len(x), median=round(sorted(x)[len(x) // 2], 2)) for k, x in sorted(rel.items())}
        out["validation"][os.path.basename(path)] = v
        if a.detail:
            os.makedirs(a.detail, exist_ok=True)
            with open(f"{a.detail}/concordance_{os.path.splitext(os.path.basename(path))[0]}.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["taxa", "sample_id"] + [col for _, col in GENES] + [f"{g}_hits" for g, _ in GENES] + ["previous_bafA_call"])
                for r, s in pairs:
                    w.writerow([r["taxa"], s["sample_id"]] + [r[col] for _, col in GENES] + [hit(s, g) for g, _ in GENES] + [s["bafA_call"]])
    kept = [r for r in calls if not r["excluded"]]
    pop = dict(kept=len(kept), genes={})
    for g, _ in GENES:
        f = lambda r, g=g: hit(r, g) >= CUT[g]
        groups = collections.defaultdict(list)
        for r in kept:
            groups[r["clone"] or "solo:" + r["sample_id"]].append(f(r))
        n = sum(f(r) for r in kept)
        pop["genes"][g] = dict(cut=CUT[g], isolates=n, clone_groups=sum(1 for x in groups.values() if sum(x) * 2 > len(x)), clone_groups_total=len(groups))
    old = lambda r: r["bafA_call"] == "P"; new = lambda r: hit(r, "bafA") >= CUT["bafA"]; neb = lambda r: r["Nebuchadnezzar-h1_call"] == "P"
    c = collections.Counter((old(r), new(r)) for r in kept)
    pop["old_vs_new_bafA"] = dict(both=c[(True, True)], new_only=c[(False, True)], old_only=c[(True, False)], neither=c[(False, False)])
    c = collections.Counter((new(r), neb(r)) for r in kept)
    pop["bafA_all_vs_neb_h1"] = dict(both=c[(True, True)], bafA_only=c[(True, False)], neb_only=c[(False, True)], neither=c[(False, False)])
    extra = [r for r in kept if new(r) and not old(r)]
    pop["newly_called"] = dict(isolates=len(extra), neb_h1_present=sum(neb(r) for r in extra))
    # two classes: cluster-linked = the AF293-type target (hrmA-cluster allele) is present; standalone divergent = only the all-alleles screen sees bafA
    classes = {}
    for name, sel in (("cluster_linked", old), ("standalone_divergent", lambda r: new(r) and not old(r))):
        rows = [r for r in kept if sel(r)]
        groups = collections.defaultdict(list)
        for r in kept:
            groups[r["clone"] or "solo:" + r["sample_id"]].append(sel(r))
        classes[name] = dict(isolates=len(rows), clone_groups=sum(1 for x in groups.values() if sum(x) * 2 > len(x)), clone_groups_total=len(groups),
                             neb_h1_present=sum(neb(r) for r in rows), neb_h1_absent=sum(not neb(r) for r in rows))
    pop["classes"] = classes
    if a.placement:
        pop["placement"] = json.load(open(a.placement))
    pop["bafA_hits_histogram"] = {str(k): n for k, n in sorted(collections.Counter(hit(r, "bafA") for r in kept if hit(r, "bafA") > 0).items())}
    pop["by_source"] = {s: dict(present=sum(new(r) for r in kept if r["source"] == s), total=sum(1 for r in kept if r["source"] == s)) for s in ("public", "local")}
    out["population"] = pop
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"kept isolates {pop['kept']}; " + "; ".join(f"{g} {x['isolates']} isolates / {x['clone_groups']} of {x['clone_groups_total']} clone groups" for g, x in pop["genes"].items()))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("targets"); t.add_argument("--hits", required=True); t.add_argument("--out", required=True); t.add_argument("--tsv"); t.add_argument("--gap", type=int, default=25)
    s = sub.add_parser("summarize"); s.add_argument("--screen", required=True); s.add_argument("--calls", required=True); s.add_argument("--counts", nargs="+", required=True)
    s.add_argument("--out", required=True); s.add_argument("--detail"); s.add_argument("--placement", help="optional counts-only JSON from the assembly check (blastn of bafA, hrmA, cgnA)")
    a = ap.parse_args(argv)
    return cmd_targets(a) if a.cmd == "targets" else cmd_summarize(a)


if __name__ == "__main__":
    sys.exit(main())
