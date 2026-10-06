"""Per flank-carried locus: strongest core tblastn hit within the locus +-PAD (bits, E, gene),
genome-wide rank of that gene's locus, genome length, independent labels. Adds the audit's cap6
Ascomycota calls (classified3.tsv; bits from the audit's own tblastn). Writes scored.tsv."""
import csv, collections, os, subprocess, yaml
from yaml import CSafeLoader as L
PAD = 20000
A = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_flank_rule_ascomycota/"
FL = {"PAP1","OBP1","PIK1"}; CORE_MTL = {"MTLA1","MTLA2","MTLalpha1","MTLalpha2"}
REAL = {"GCA_054906655.1_Umbelopsis_nana_v._1.0","GCA_964291815.1_UMBE_WA70503","GCA_053572175.1_ASM5357217v1",
        "GCA_977110975.1_gzUmbVina2","GCA_016758895.1_ASM1675889v1",          # Umbelopsis: Mucorales gene order (diagnosis)
        "GCA_000587855.1_B50","GCA_000697435.1_RhiVarB7584-1.0",              # Mucor irregularis (audit: best genome-wide sexM)
        "GCA_003707065.3"}                                                     # Trigonopsis variabilis (audit)
def modelled(e): return (e.get("status") or "").startswith("polished") or e.get("method") == "diamond_proteome"
def hits(g):
    p = f"hits/{g}.tsv.zst"
    if not os.path.exists(p): return None
    return [l.split("\t") for l in subprocess.run(["zstdcat", p], capture_output=True, text=True).stdout.splitlines()]
def loci_of(hs, gene):
    hh = sorted(((float(h[9]), float(h[8]), h[1], *sorted((int(h[6]), int(h[7]))), float(h[2])) for h in hs
                 if h[0].split("|")[-1] == gene), reverse=True)
    out = []
    for bs, ev, c, lo, hi, pid in hh:
        if any(c == l[2] and abs(lo - l[3]) < 5000 for l in out): continue
        out.append((bs, ev, c, lo, hi, pid))
    return out
w = csv.writer(open("scored.tsv", "w"), delimiter="\t")
w.writerow(["source","genome","phylum","order","species","family","locus","outcome","best_gene","bits","evalue","pid",
            "rank_genomewide","genome_len","label","label_basis"])
for r in csv.DictReader(open("loci.tsv"), delimiter="\t"):
    g = r["genome"]; hs = hits(g)
    if hs is None: continue
    c0, s0, e0 = r["contig"], int(r["start"]) - PAD, int(r["end"]) + PAD
    best = None
    for gene in sorted({h[0].split("|")[-1] for h in hs}):
        lo = loci_of(hs, gene)
        for i, l in enumerate(lo):
            if l[2] == c0 and s0 <= l[3] and l[4] <= e0:
                if best is None or l[0] > best[1]: best = (gene, l[0], l[1], l[5], i + 1)
                break
    L_ = int(open(f"hits/{g}.len").read().strip()) if os.path.exists(f"hits/{g}.len") else ""
    label, basis = "unknown", ""
    if any(g.startswith(x) for x in REAL): label, basis = "real", "curated diagnosis"
    if r["run"] == "serinales_882aa01":
        rep = yaml.load(open(f'{r["rundir"]}/detection_report.yaml'), Loader=L)
        det = rep.get("detected") or []
        has_core = any(any(e["gene"] in CORE_MTL and modelled(e) for e in x.get("gene_evidence") or []) for x in det)
        x = next(x for x in det if x["contig"] == r["contig"] and x["start"] == int(r["start"]))
        ev = x.get("gene_evidence") or []; f = [e for e in ev if e["gene"] in FL]
        inside = bool(f) and all(min(e["start"] for e in f) - 3000 <= e["start"] and e["end"] <= max(e["end"] for e in f) + 3000
                                 for e in ev if e["gene"] in CORE_MTL)
        if has_core and not inside: label, basis = "noise", "Serinales: second call outside PAP1-OBP1-PIK1 span in a genome with a core-modelled call"
        elif not has_core and inside: basis = "Serinales: only locus, core inside flank span"
    w.writerow([r["run"], g, r["phylum"], r["order"], r["species"], r["family"], f'{c0}:{r["start"]}-{r["end"]}',
                r["outcome"], *(best or ("none", 0, "", "", "absent")), L_, label, basis])
# audit cap6 Ascomycota calls (bits from the audit's tblastn; genome length not recorded there)
alen = dict(l.split() for l in open("audit_len.tsv"))
aud = collections.defaultdict(list)
for x in csv.DictReader(open(A + "classified3.tsv"), delimiter="\t"):
    if x["panel"] == "cap6": aud[(x["genome"], x["call"], x["mat_family"])].append(x)
for (g, call, fam), v in aud.items():
    v2 = [x for x in v if x["here_bits"]]
    b = max(v2, key=lambda x: float(x["here_bits"])) if v2 else None
    label = "real" if g.startswith("GCA_003707065.3") else "unknown"
    w.writerow(["audit_cap6", g, "Ascomycota", v[0]["order"], v[0]["species"], fam, call, v[0]["outcome"],
                b["gene"] if b else "none", b["here_bits"] if b else 0, b["here_evalue"] if b else "", b["here_pid"] if b else "",
                b["rank_here"] if b else "absent", alen.get(g, ""), label, "audit" if label == "real" else ""])
print("done")
