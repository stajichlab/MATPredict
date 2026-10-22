"""Pilot list for the clean full-BFD run (docs/notes/2026-10-22_clean-run-plan.md, phase 3). Deterministic (seed 20261022).
Run on the HPCC (needs the BFD tables):  python build_pilot_list.py V060_ASCO.tsv V060_BASIDIO.tsv MUCO_INPUTS.tsv MUCO_KEPT.tsv OUTDIR
Strata: (1) base: stratified by CLASS like the 246-genome v0.6.0 pilot (about 200 x class share, at least 3 per class, round-robin over orders);
(2) large genomes (>= 500 Mb) to exercise the 14400 s path; (3) known problem clades; (4) Mucoromycotina (the separate no-taxid path)."""
import csv, collections, os, random, sys, duckdb
SEED = 20261022
BFD = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
v_asco, v_bas, m_inputs, m_kept, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
rng = random.Random(SEED)
samples = {r["ASMID"]: r for r in csv.DictReader(open(BFD + "samples.csv"))}
size = dict(duckdb.sql(f"select ASMID, total_length_bp from read_parquet('{BFD}tables/asm_stats.parquet')").fetchall())
v = {}
for phylum, f in (("Ascomycota", v_asco), ("Basidiomycota", v_bas)):
    for r in csv.DictReader(open(f), delimiter="\t"):
        v[r["genome"]] = dict(phylum=phylum, v060_status=r["status"], v060_routing=r["routing"], v060_n_loci=r["n_loci"])
usable = {}
for g, info in v.items():
    s = samples.get(g)
    if not s or not (s.get("NCBI_TAXONID") or "").strip() or not os.path.exists(f"{LIB}{g}.fa.gz"):
        continue
    usable[g] = dict(info, taxid=s["NCBI_TAXONID"].strip(), cls=s["CLASS"], order=s["ORDER"], family=s["FAMILY"], genus=s["GENUS"],
                     species=s["SPECIES"], size_bp=size.get(g), size_bin="large" if (size.get(g) or 0) >= 500_000_000 else "small")
print("usable genomes:", len(usable), collections.Counter(u["phylum"] for u in usable.values()))
picked = collections.OrderedDict()   # genome -> list of strata
def add(g, stratum):
    picked.setdefault(g, []).append(stratum)

def stratified(pool, total, label):
    """about total x class share, at least 3 per class, round-robin over orders inside a class."""
    by_class = collections.defaultdict(list)
    for g in pool: by_class[usable[g]["cls"] or "NA"].append(g)
    n = len(pool); chosen = []
    for cls, gs in sorted(by_class.items()):
        want = max(3, round(total * len(gs) / n)); want = min(want, len(gs))
        by_order = collections.defaultdict(list)
        for g in sorted(gs): by_order[usable[g]["order"] or "NA"].append(g)
        orders = sorted(by_order); rng.shuffle(orders)
        for o in orders: rng.shuffle(by_order[o])
        got = []
        while len(got) < want and any(by_order[o] for o in orders):
            for o in orders:
                if by_order[o] and len(got) < want: got.append(by_order[o].pop())
        chosen += got
    for g in chosen: add(g, label)
    return chosen

def take(pool, k, label, prefer_uncalled=0):
    pool = sorted(pool); rng.shuffle(pool)
    unc = [g for g in pool if usable[g]["v060_status"] != "called"]; cal = [g for g in pool if usable[g]["v060_status"] == "called"]
    got = unc[:prefer_uncalled]; rest = [g for g in cal + unc[prefer_uncalled:] if g not in got]; rng.shuffle(rest)
    got += rest[: max(0, k - len(got))]
    for g in got: add(g, label)
    return got

small = [g for g, u in usable.items() if u["size_bin"] == "small"]
asco_small = [g for g in small if usable[g]["phylum"] == "Ascomycota"]; bas_small = [g for g in small if usable[g]["phylum"] == "Basidiomycota"]
stratified(asco_small, 180, "base_ascomycota_by_class"); stratified(bas_small, 60, "base_basidiomycota_by_class")
# large genomes: Basidiomycota (57 Pucciniales) and any Ascomycota >= 500 Mb
big = [g for g, u in usable.items() if u["size_bin"] == "large"]
take([g for g in big if usable[g]["phylum"] == "Basidiomycota"], 10, "large_genome_basidiomycota"); take([g for g in big if usable[g]["phylum"] == "Ascomycota"], 5, "large_genome_ascomycota")
G = lambda *names, **kw: [g for g, u in usable.items() if u["genus"] in names and u["phylum"] == kw.get("phylum", "Ascomycota")]
F = lambda *names: [g for g, u in usable.items() if u["family"] in names and u["phylum"] == "Ascomycota"]
O = lambda *names: [g for g, u in usable.items() if u["order"] in names and u["phylum"] == "Ascomycota"]
# problem clades (counts chosen to keep the pilot near 500 genomes; uncalled genomes are over-represented where recall is the question)
take(G("Bipolaris", "Curvularia", "Exserohilum", "Cochliobolus"), 40, "problem_cochliobolus_clade", prefer_uncalled=30)
take(G("Pyrenophora"), 25, "problem_pyrenophora_homothallic_candidates")
take(G("Aureobasidium"), 25, "problem_aureobasidium_fragment_models", prefer_uncalled=5)
take(G("Botrytis", "Botryotinia", "Sclerotinia", "Monilinia", "Stromatinia"), 30, "problem_sclerotiniaceae_remnant_fused_loci")
take(F("Cordycipitaceae", "Clavicipitaceae", "Ophiocordycipitaceae"), 25, "problem_hypocreales_three_exon_MAT1-1-1", prefer_uncalled=3)
take(O("Dipodascales", "Lipomycetales"), 25, "problem_dipodascales_lipomycetales_guard", prefer_uncalled=8)
take(G("Ceratocystis"), 10, "problem_ceratocystis", prefer_uncalled=4)
take(G("Komagataella"), 8, "problem_komagataella", prefer_uncalled=4); take(G("Ogataea"), 8, "problem_ogataea")
take(O("Xylariales"), 15, "problem_xylariales_absent_MAT"); take(G("Alternaria"), 10, "problem_alternaria"); take(G("Parastagonospora"), 5, "problem_parastagonospora"); take(G("Fusarium"), 10, "problem_fusarium")
take(G("Hemileia", phylum="Basidiomycota"), 5, "problem_hemileia_rustHD"); take(G("Puccinia", phylum="Basidiomycota"), 5, "problem_puccinia_rustHD")
# write outputs
with open(f"{out}/pilot_list.tsv", "w") as fl, open(f"{out}/pilot_meta.tsv", "w", newline="") as fm:
    w = csv.writer(fm, delimiter="\t"); w.writerow(["ASMID", "taxid", "phylum", "class", "order", "family", "genus", "species", "size_bp", "size_bin", "v060_status", "v060_routing", "strata"])
    for g, strata in picked.items():
        u = usable[g]; fl.write(f"{g}\t{u['taxid']}\n")
        w.writerow([g, u["taxid"], u["phylum"], u["cls"], u["order"], u["family"], u["genus"], u["species"], u["size_bp"], u["size_bin"], u["v060_status"], u["v060_routing"], ";".join(strata)])
# Mucoromycotina: the separate path (detect --phylum Mucoromycota --emit-cds-fasta, no taxid): 15 BFD + 5 LCG + 5 Jena
kept = {l.split("\t")[0] for l in open(m_kept)}
rows = list(csv.DictReader(open(m_inputs), delimiter="\t"))
bfd = [r for r in rows if r["source"] == "BFD" and r["id"] in kept]; lcg = [r for r in rows if r["source"] == "LCG"]; jena = [r for r in rows if r["source"] == "Jena"]
mrows = []
for pool, k in ((bfd, 15), (lcg, 5), (jena, 5)):
    pool = sorted(pool, key=lambda r: r["id"]); rng.shuffle(pool); mrows += pool[:k]
with open(f"{out}/pilot_mucoromycotina.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(mrows)
# summary
print("\npilot genomes (Ascomycota + Basidiomycota, detect-per-genome path):", len(picked), "| Mucoromycotina path:", len(mrows), "| total", len(picked) + len(mrows))
print("by phylum:", dict(collections.Counter(usable[g]["phylum"] for g in picked)), " by size bin:", dict(collections.Counter(usable[g]["size_bin"] for g in picked)))
print("by v0.6.0 status:", dict(collections.Counter(usable[g]["v060_status"] for g in picked)))
print("strata (a genome can belong to several):"); [print(f"  {k:52s} {v}") for k, v in collections.Counter(s for ss in picked.values() for s in ss).most_common()]
cls = collections.Counter(usable[g]["cls"] for g in picked if usable[g]["phylum"] == "Ascomycota"); print("Ascomycota classes covered:", len(cls), " of 30; Basidiomycota classes:", len({usable[g]['cls'] for g in picked if usable[g]['phylum']=='Basidiomycota'}))
print("multi-stratum genomes:", sum(len(s) > 1 for s in picked.values()))
