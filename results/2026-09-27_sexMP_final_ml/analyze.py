"""Analyse the IQ-TREE and RAxML-NG trees of the trimmed Mucoromycotina set.

Writes tip_names.tsv (for draw_tree_pdf.py), analysis.txt, per_leaf.tsv.
Clade definitions follow read_tree.py: the sexP clade is the smallest side of a
bipartition holding all sexP references and no sexM reference / MAT1-2-1; the
same for sexM, or, if not monophyletic, the pure side holding the most sexM refs.
"""
import collections, csv, glob
import yaml
from Bio import Phylo

DB = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db"
samp = {r["ASMID"]: r for r in csv.DictReader(open(
    "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
ids = dict(l.rstrip("\n").split("\t") for l in list(open("taxon_ids.tsv"))[1:])
sel = {r["seqid"]: r for r in csv.DictReader(open("selection.tsv"), delimiter="\t")}
rec_sp = {}
for f in glob.glob(f"{DB}/*/*/*/metadata.yaml"):
    try:
        rec_sp[f.split("/")[-2]] = yaml.safe_load(open(f))["organism"]["species"]
    except Exception:
        pass

meta = {}
for tid, sid in ids.items():
    r = sel[sid]
    if sid.startswith(("REF|", "OUT|")):
        p = sid.split("|")
        rec = p[2] if sid.startswith("OUT|") else p[1]
        kind = "REF_MAT1-2-1" if sid.startswith("OUT|") else f"REF_{p[-1]}"
        meta[tid] = dict(tip=tid, kind=kind, species=rec_sp.get(rec, rec), family="", status="reference",
                         idiomorph_label="", best_ref_gene=kind[4:], tree_clade="", n_collapsed=1,
                         source=sid)
    else:
        status = {"called": "call"}.get(r["status"], r["status"])
        meta[tid] = dict(tip=tid, kind=r["kind"], species=r["species"], family=r["family"], status=status,
                         idiomorph_label=r["label"] or "-", best_ref_gene="", tree_clade="",
                         n_collapsed=r["n_members"], source=r.get("rep_member") or sid, why=r["why"])

refs = {k: {t for t, m in meta.items() if m["kind"] == f"REF_{k}"} for k in ("sexM", "sexP", "MAT1-2-1")}
ROOT = next(t for t, m in meta.items() if m["kind"] == "REF_MAT1-2-1" and "5518_3639" in m["source"])


def load(path):
    t = Phylo.read(path, "newick")
    t.root_with_outgroup(next(x for x in t.get_terminals() if x.name == ROOT))
    return t


def conf(cl):
    """(SH-aLRT, bootstrap) from an IQ-TREE 'SH/UF' label or a plain support value."""
    s = cl.name if cl.name else (str(cl.confidence) if cl.confidence is not None else "")
    p = s.split("/")
    try:
        return (float(p[0]), float(p[1])) if len(p) == 2 else (None, float(p[0]))
    except ValueError:
        return None, None


def clades(tree):
    leaves = {x.name for x in tree.get_terminals()}
    sides = []
    for cl in tree.find_clades():
        if cl.is_terminal() or cl == tree.root:
            continue
        S = {x.name for x in cl.get_terminals()}
        a, b = conf(cl)
        sides += [(S, a, b), (leaves - S, a, b)]
    out = {}
    for T, O in (("sexP", "sexM"), ("sexM", "sexP")):
        X = refs[O] | refs["MAT1-2-1"]
        ok = [s for s in sides if refs[T] <= s[0] and not (s[0] & X)]
        if ok:
            S, a, b = min(ok, key=lambda s: len(s[0])); mono = True
        else:
            pure = [s for s in sides if not (s[0] & X)]
            S, a, b = max(pure, key=lambda s: (len(s[0] & refs[T]), -len(s[0]))); mono = False
        out[T] = dict(S=S, sh=a, bs=b, mono=mono, nref=len(S & refs[T]))
    return out


def leaf_support(tree, tid, clade):
    """Support of the smallest clade that holds this leaf AND at least one curated
    reference of `clade` (sexP or sexM), with no reference of the other type or
    MAT1-2-1 -- i.e. how well the tree joins this leaf to that reference set."""
    if clade not in ("sexP", "sexM"):
        return (None, None), 0
    other = refs["sexM" if clade == "sexP" else "sexP"] | refs["MAT1-2-1"]
    path = tree.get_path(next(x for x in tree.get_terminals() if x.name == tid))
    for node in reversed(path[:-1]):
        S = {x.name for x in node.get_terminals()}
        if S & other:
            break
        if S & refs[clade]:
            return conf(node), len(S)
    return (None, None), 0


out = open("analysis.txt", "w")
def say(*a):
    print(*a); print(*a, file=out)

trees = {}
for name, path in (("IQ-TREE", "iq.treefile"), ("RAxML-NG", "rx.raxml.support")):
    try:
        trees[name] = load(path)
    except FileNotFoundError:
        say(f"{name}: {path} missing")
res = {n: clades(t) for n, t in trees.items()}
for n, r in res.items():
    for T in ("sexP", "sexM"):
        c = r[T]
        say(f"{n} {T}: {'monophyletic' if c['mono'] else 'NOT monophyletic'} "
            f"({c['nref']}/{len(refs[T])} refs), size {len(c['S'])}, SH-aLRT {c['sh']}, bootstrap {c['bs']}")

# label / clade concordance on called loci (IQ-TREE clades)
if "IQ-TREE" in res:
    r = res["IQ-TREE"]
    for t, m in meta.items():
        m["tree_clade"] = "sexP" if t in r["sexP"]["S"] else "sexM" if t in r["sexM"]["S"] else "other_HMG"
        if m["status"] == "reference":
            m["tree_clade"] = ""
    conc = collections.Counter()
    dis = []
    for t, m in meta.items():
        if m["status"] != "call" or m["tree_clade"] == "other_HMG":
            continue
        lab = m["idiomorph_label"]
        agree = (m["tree_clade"] == "sexP" and lab == "Plus") or (m["tree_clade"] == "sexM" and lab == "Minus")
        conc[(m["tree_clade"], lab)] += 1
        if not agree:
            dis.append((m["species"], m["source"].split("|")[0], lab, m["tree_clade"]))
    say("concordance (called loci placed in a clade), IQ-TREE:", dict(conc))
    say("disagreements:", len(dis))
    for d in dis:
        say("   ", d)

# per-leaf support for the curator's tier-2 question
keys = ("Syncephalastrum", "Rhizomucor", "Lichtheimia")
with open("per_leaf.tsv", "w") as fo:
    fo.write("tip\tspecies\tgenome\tstatus\tlabel\tclade_IQ\tIQ_SH_join_ref\tIQ_UFBoot_join_ref\tIQ_join_size\tclade_RX\tRX_bootstrap_join_ref\n")
    for t, m in meta.items():
        if not any(k in m["species"] for k in keys):
            continue
        row = [t, m["species"], m["source"].split("|")[0], m["status"], m["idiomorph_label"]]
        for n in ("IQ-TREE", "RAxML-NG"):
            if n not in trees:
                row += ["", "", ""] if n == "IQ-TREE" else ["", ""]; continue
            cl = res[n]
            c = "sexP" if t in cl["sexP"]["S"] else "sexM" if t in cl["sexM"]["S"] else "other_HMG"
            (sh, bs), size = leaf_support(trees[n], t, c)
            row += [c, sh, bs, size] if n == "IQ-TREE" else [c, bs]
        fo.write("\t".join(str(x) for x in row) + "\n")
say(open("per_leaf.tsv").read())

# topology agreement
if len(trees) == 2:
    def bip(t):
        L = frozenset(x.name for x in t.get_terminals())
        s = set()
        for cl in t.find_clades():
            if cl.is_terminal() or cl == t.root:
                continue
            S = frozenset(x.name for x in cl.get_terminals())
            if 1 < len(S) < len(L) - 1:
                s.add(min(S, L - S, key=lambda z: sorted(z)))
        return s
    a, b = bip(trees["IQ-TREE"]), bip(trees["RAxML-NG"])
    rf = len(a ^ b)
    say(f"Robinson-Foulds IQ-TREE vs RAxML-NG: {rf} of {len(a) + len(b)} (normalised {rf / (len(a) + len(b)):.3f}); "
        f"shared splits {len(a & b)}")
    for T in ("sexP", "sexM"):
        say(f"{T} clade identical in both: {res['IQ-TREE'][T]['S'] == res['RAxML-NG'][T]['S']} "
            f"(IQ {len(res['IQ-TREE'][T]['S'])} tips, RX {len(res['RAxML-NG'][T]['S'])} tips)")

with open("tip_names.tsv", "w", newline="") as fo:
    cols = ["name", "tip", "kind", "species", "family", "status", "idiomorph_label", "best_ref_gene",
            "tree_clade", "n_collapsed", "source"]
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for t, m in meta.items():
        w.writerow({"name": t, **m})
out.close()
