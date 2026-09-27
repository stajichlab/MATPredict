"""sexM vs sexP profile-HMM discrimination test.

Scores MODELLED proteins (no 6-frame ORFs), so this measures discrimination only.

Labelled set (truth): 15 curated references (9 sexM, 6 sexP) + the HMG-box
protein annotated inside each of the 23 Zygo truth loci (7 Minus = sexM,
16 Plus = sexP).  Leave-one-GENUS-out: for each genus, models are built from
the labelled set minus that genus and scored on that genus.

Arms
  full_bal   full-length HMMs, labelled set only
  dom_bal    HMG-box-only HMMs (PF00505 envelope +-5 aa), labelled set only
  full_aug   full-length, sexP set augmented with sexP-clade tree members
             (UFBoot 99 clade), minus the held-out genus and minus every
             disputed-set protein.  sexM is never augmented (clade UFBoot 57).
  dom_aug    same, HMG-box only
Baseline
  blast      blastp of the full protein vs the same training proteins as
             full_bal; best bitscore per type (approximates detect's
             best-bitscore rule on proteins, not its tblastn on genomes).
"""
import collections, csv, glob, os, re, subprocess, sys, tempfile
import yaml

os.environ["PATH"] = ("/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.4/bin:"
                      "/opt/linux/rocky/8.x/x86_64/pkgs/mafft/7.505/bin:"
                      "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:"
                      + os.environ["PATH"])
PH = "../2026-09-26_sexMP_phylogeny"
DB = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
PFAM = f"{PH}/PF00505.hmm"
os.environ["MAFFT_BINARIES"] = "/opt/linux/rocky/8.x/x86_64/pkgs/mafft/7.505/libexec/mafft"
TMP = tempfile.mkdtemp(prefix="hmmloo_", dir=os.environ.get("SCRATCH", "."))


def read_fa(path):
    seqs, cur = {}, None
    for l in open(path):
        if l.startswith(">"):
            cur = l[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(l.strip())
    return {k: "".join(v).replace("*", "") for k, v in seqs.items()}


def write_fa(path, d):
    with open(path, "w") as fo:
        for k, v in d.items():
            fo.write(f">{k}\n{v}\n")


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# ---------------------------------------------------------------- sequences
allp = read_fa(f"{PH}/all_proteins.faa")
zygo = read_fa("zygo_locus_proteins.faa")
# keep only the HMG-box protein per Zygo locus (exactly one each, verified)
zdom = set()
for l in open("zygo_pf.domtbl"):
    if not l.startswith("#"):
        zdom.add(l.split()[0])
zygo = {k: v for k, v in zygo.items() if k in zdom}
seqs = {**allp, **zygo}

# HMG-box envelopes for every protein
write_fa(f"{TMP}/all.faa", seqs)
run(f"hmmsearch --cpu 8 --domtblout {TMP}/pf.domtbl -E 1e-2 --domE 1e-2 {PFAM} {TMP}/all.faa")
env = {}
for l in open(f"{TMP}/pf.domtbl"):
    if l.startswith("#"):
        continue
    f = l.split()
    sid, s, e, ce = f[0], int(f[19]), int(f[20]), float(f[11])
    if sid not in env or ce < env[sid][2]:
        env[sid] = (s, e, ce)
dom = {k: seqs[k][max(0, s - 6):e + 5] for k, (s, e, _) in env.items()}

# ---------------------------------------------------------------- labels, genus
rec_sp = {}
for f in glob.glob(f"{DB}/*/*/*/metadata.yaml"):
    try:
        rec_sp[f.split("/")[-2]] = yaml.safe_load(open(f))["organism"]["species"]
    except Exception:
        pass
labelled = {}   # id -> (type, genus, source)
for k in allp:
    if k.startswith("REF|"):
        rec, typ = k.split("|")[1], k.split("|")[-1]
        labelled[k] = (typ, rec_sp.get(rec, rec).split()[0], "ref")
for k in zygo:
    org, idio = k.split("|")[1], k.split("|")[2]
    labelled[k] = ("sexP" if idio == "Plus" else "sexM", org.split("_")[0], "zygo")

samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
cand = {r["candidate_id"]: r for r in csv.DictReader(open(f"{PH}/candidates.tsv"), delimiter="\t")}
loci = {r["locus_id"]: r for r in csv.DictReader(open(f"{PH}/loci.tsv"), delimiter="\t")}


def genus_of(cid):
    r = cand[cid]
    sp = r["species"] or samp.get(r["genome"], {}).get("SPECIES", "") or r["genome"]
    return sp.split()[0]


# ---------------------------------------------------------------- disputed sets
disputed = collections.OrderedDict()
disputed["minus_in_sexP_clade"] = [c for c, r in cand.items() if r["status"] == "called"
                                   and r["idiomorph_label"] == "Minus" and r["tree_clade"] == "sexP"]
rhiz = []
for c, r in cand.items():
    L = loci.get(r["locus_or_copy"], {})
    g = set(L.get("genes_found", "").split(","))
    if (r["kind"] == "cand" and r["status"] == "called" and r["idiomorph_label"] == "Plus"
            and r["species"].startswith("Rhizopus") and "sexM" in g and "sexP" not in g
            and c in env):
        rhiz.append(c)
disputed["rhizopus_sexM_only_Plus"] = rhiz
disputed["lichth_syncephal_sexP_clade"] = [c for c, r in cand.items() if r["tree_clade"] == "sexP"
                                           and r["family"] in ("Lichtheimiaceae", "Syncephalastraceae")]
disputed_ids = {c for v in disputed.values() for c in v}
# dedup groups: never train on a protein identical to a disputed one
members = collections.defaultdict(set)
for r in csv.DictReader(open(f"{PH}/dedup_map.tsv"), delimiter="\t"):
    members[r["rep"]].add(r["member"])
grp_of = {m: rep for rep, ms in members.items() for m in ms}
disputed_grp = {grp_of.get(c, c) for c in disputed_ids}
aug_P = [c for c, r in cand.items() if r["tree_clade"] == "sexP" and c in env
         and c not in disputed_ids and grp_of.get(c, c) not in disputed_grp]


# ---------------------------------------------------------------- helpers
def build(name, ids, use_dom):
    src = dom if use_dom else seqs
    d = {i: src[i] for i in ids if i in src and len(src[i]) >= 30}
    fa, aln, hmm = f"{TMP}/{name}.faa", f"{TMP}/{name}.afa", f"{TMP}/{name}.hmm"
    write_fa(fa, d)
    run(f"mafft --auto --quiet --thread 8 {fa} > {aln}")
    run(f"hmmbuild --amino -n {name} {hmm} {aln}")
    return hmm, len(d)


def score(hmm, ids, use_dom, tag):
    src = dom if use_dom else seqs
    fa, tbl = f"{TMP}/q_{tag}.faa", f"{TMP}/q_{tag}.tbl"
    write_fa(fa, {i: src[i] for i in ids if i in src})
    run(f"hmmsearch -E 1000 --domE 1000 -Z 1 --tblout {tbl} {hmm} {fa}")
    out = {}
    for l in open(tbl):
        if not l.startswith("#"):
            f = l.split(); out[f[0]] = float(f[5])
    return {i: out.get(i, 0.0) for i in ids}


def blast_scores(train, test, tag):
    db, q, o = f"{TMP}/bdb_{tag}", f"{TMP}/bq_{tag}.faa", f"{TMP}/bo_{tag}.tsv"
    write_fa(f"{db}.faa", {i: seqs[i] for i in train})
    run(f"makeblastdb -in {db}.faa -dbtype prot -out {db}")
    write_fa(q, {i: seqs[i] for i in test})
    run(f"blastp -query {q} -db {db} -outfmt '6 qseqid sseqid bitscore' -evalue 10 -max_target_seqs 500 -out {o}")
    best = collections.defaultdict(lambda: {"sexM": 0.0, "sexP": 0.0})
    for l in open(o):
        qq, ss, b = l.split("\t"); t = labelled[ss][0]
        best[qq][t] = max(best[qq][t], float(b))
    return {i: best[i]["sexP"] - best[i]["sexM"] for i in test}


def identity_spread(ids):
    d = {i: seqs[i] for i in ids}
    fa, aln = f"{TMP}/div.faa", f"{TMP}/div.afa"
    write_fa(fa, d)
    run(f"mafft --auto --quiet --thread 8 {fa} > {aln}")
    a = read_fa(aln); ks = list(a); ids_ = []
    for x in range(len(ks)):
        for y in range(x + 1, len(ks)):
            p, q = a[ks[x]], a[ks[y]]
            cols = [(i, j) for i, j in zip(p, q) if i != "-" and j != "-"]
            if cols:
                ids_.append(100 * sum(i == j for i, j in cols) / len(cols))
    ids_.sort()
    uniq = len(set(d.values()))
    if not ids_:
        return uniq, None
    return uniq, (round(ids_[0], 1), round(ids_[len(ids_) // 2], 1), round(ids_[-1], 1))


# ---------------------------------------------------------------- diversity
lab_M = [k for k, v in labelled.items() if v[0] == "sexM"]
lab_P = [k for k, v in labelled.items() if v[0] == "sexP"]
rep = open("results.txt", "w")
def say(*a):
    print(*a); print(*a, file=rep)
say("== Training-set diversity (full-length; pairwise identity over co-aligned columns: min/median/max)")
for nm, ids in (("sexM labelled", lab_M), ("sexP labelled", lab_P), ("sexP clade augment", aug_P)):
    u, s = identity_spread(ids)
    say(f"  {nm:22s} n={len(ids):4d} unique={u:4d} identity {s}")
say(f"  disputed sets: " + ", ".join(f"{k}={len(v)}" for k, v in disputed.items()))

# ---------------------------------------------------------------- LOO
genera = sorted({v[1] for v in labelled.values()})
arms = {"full_bal": (False, False), "dom_bal": (True, False), "full_aug": (False, True), "dom_aug": (True, True)}
margins = {a: {} for a in list(arms) + ["blast"]}
for G in genera:
    test = [k for k, v in labelled.items() if v[1] == G]
    trM = [k for k in lab_M if labelled[k][1] != G]
    trP = [k for k in lab_P if labelled[k][1] != G]
    augP = [c for c in aug_P if genus_of(c) != G]
    for arm, (use_dom, aug) in arms.items():
        hM, _ = build(f"M_{arm}", trM, use_dom)
        hP, _ = build(f"P_{arm}", trP + (augP if aug else []), use_dom)
        sM, sP = score(hM, test, use_dom, f"M{arm}"), score(hP, test, use_dom, f"P{arm}")
        for t in test:
            margins[arm][t] = sP[t] - sM[t]
    margins["blast"].update(blast_scores(trM + trP, test, G))

say("\n== Leave-one-genus-out on 38 labelled proteins (%d genera)" % len(genera))
say("  arm        correct/38  sexM ok  sexP ok  ties  margin(correct side) median [min..max]")
for arm, m in margins.items():
    ok = {"sexM": 0, "sexP": 0}; ties = 0; signed = []
    for t, v in m.items():
        typ = labelled[t][0]
        if v == 0:
            ties += 1; continue
        call = "sexP" if v > 0 else "sexM"
        ok[typ] += call == typ
        signed.append(v if typ == "sexP" else -v)
    signed.sort()
    say(f"  {arm:9s}  {ok['sexM']+ok['sexP']:3d}/38      {ok['sexM']:2d}/{len(lab_M)}    {ok['sexP']:2d}/{len(lab_P)}    {ties:2d}"
        f"    {signed[len(signed)//2]:.1f} [{signed[0]:.1f}..{signed[-1]:.1f}]")
with open("loo_margins.tsv", "w") as fo:
    fo.write("id\ttruth\tgenus\tsource\t" + "\t".join(margins) + "\n")
    for t in labelled:
        fo.write(f"{t}\t{labelled[t][0]}\t{labelled[t][1]}\t{labelled[t][2]}\t" +
                 "\t".join(f"{margins[a].get(t, 0):.1f}" for a in margins) + "\n")

# ---------------------------------------------------------------- disputed sets
say("\n== Disputed sets (models built on ALL labelled proteins; aug excludes disputed proteins)")
final = {}
for arm, (use_dom, aug) in arms.items():
    hM, nM = build(f"FM_{arm}", lab_M, use_dom)
    hP, nP = build(f"FP_{arm}", lab_P + (aug_P if aug else []), use_dom)
    final[arm] = (hM, hP, use_dom)
    say(f"  {arm}: sexM HMM from {nM} seqs, sexP HMM from {nP} seqs")
rows = []
for setname, ids in disputed.items():
    ids = [i for i in ids if i in seqs]
    res = {}
    for arm, (hM, hP, use_dom) in final.items():
        sM, sP = score(hM, ids, use_dom, "dM"), score(hP, ids, use_dom, "dP")
        res[arm] = {i: sP[i] - sM[i] for i in ids}
    res["blast"] = blast_scores(lab_M + lab_P, ids, "disp")
    say(f"  {setname} (n={len(ids)}):")
    for arm in list(arms) + ["blast"]:
        v = list(res[arm].values())
        nP_ = sum(x > 0 for x in v); nM_ = sum(x < 0 for x in v)
        s = sorted(v)
        say(f"    {arm:9s} sexP {nP_:3d}  sexM {nM_:3d}  tie {len(v)-nP_-nM_:2d}   margin median {s[len(s)//2]:.1f} [{s[0]:.1f}..{s[-1]:.1f}]" if v else f"    {arm}: none")
    for i in ids:
        r = cand[i]
        rows.append([setname, i, r["species"], r["family"], r["status"], r["idiomorph_label"], r["tree_clade"]] +
                    [f"{res[a][i]:.1f}" for a in list(arms) + ["blast"]])
with open("disputed_scores.tsv", "w") as fo:
    fo.write("set\tid\tspecies\tfamily\tstatus\tlabel\ttree_clade\t" + "\t".join(list(arms) + ["blast"]) + "\n")
    for r in rows:
        fo.write("\t".join(r) + "\n")

# agreement of the HMM with tree placement on all tree-placed BFD candidates
say("\n== Agreement with tree placement (BFD candidates in sexM or sexP clade; not accuracy)")
tp = [c for c, r in cand.items() if r["tree_clade"] in ("sexM", "sexP") and c in seqs and c not in disputed_ids]
for arm in ("full_bal", "dom_bal"):
    hM, hP, use_dom = final[arm]
    sM, sP = score(hM, tp, use_dom, "tM"), score(hP, tp, use_dom, "tP")
    agree = collections.Counter()
    for c in tp:
        v = sP[c] - sM[c]; call = "sexP" if v > 0 else "sexM" if v < 0 else "tie"
        agree[(cand[c]["tree_clade"], call)] += 1
    say(f"  {arm}: " + ", ".join(f"tree {a} -> hmm {b}: {n}" for (a, b), n in sorted(agree.items())))
rep.close()
print("tmp:", TMP)
