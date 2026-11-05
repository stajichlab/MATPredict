"""Gene models of the six candidate MATog windows: the miniprot projection of the two curated MATog records (DL-1 and CBS 4732) onto each
contig. Checks every model (translation, internal stops, start/stop, coverage of the reference protein, distance to the contig end) and writes
models.tsv. The projection is our own (tier 3): nothing here is an annotation of the assembly."""
import glob, re, sys
from Bio.Seq import Seq
OGREFS = ("1005962_dl-1_MATog", "460523_cbs-4732_MATog")
REFLEN = {}
for fa in ("refs_og.faa",):
    n = None
    for l in open(fa):
        if l[0] == ">": n = l[1:].strip(); REFLEN[n] = 0
        else: REFLEN[n] += len(l.strip())

def genome(fa):
    return "".join(l.strip() for l in open(fa) if l[0] != ">")

def models(slug_glob):
    """Best clean model per gene: the highest-scoring projection of a MATog reference protein, from the default and the relaxed miniprot
    runs, that has no frameshift, no internal stop and a CDS length divisible by 3. Rejected projections are listed in `rejected`."""
    f = glob.glob(f"data/{slug_glob}*.miniprot.gff")[0]
    seq = genome(f.replace(".miniprot.gff", ".fa")); out = {}
    mrna = {}
    lines = [(l, "default") for l in open(f)] + [(l, "relaxed") for l in open(f.replace(".miniprot.gff", ".relaxed.gff"))]
    for l, runname in lines:
        if l.startswith("#"): continue
        c = l.rstrip().split("\t")
        if len(c) < 9 or c[2] not in ("mRNA", "CDS"): continue
        a = dict(kv.split("=", 1) for kv in c[8].split(";"))
        tgt = a["Target"].split(" ")[0]
        if not tgt.startswith(OGREFS): continue
        mm = re.search(r"name=([^|]+)", tgt)
        name = mm.group(1) if mm else tgt.rsplit("_", 1)[-1]
        if c[2] == "mRNA":
            t0, t1 = [int(x) for x in a["Target"].split(" ")[1:3]]
            mrna[(runname, a["ID"])] = dict(run=runname, frameshift="Frameshift" in a, name=name, ref=tgt.split("|")[0] if mm else tgt.rsplit("_", 1)[0], start=int(c[3]), end=int(c[4]), score=int(c[5]), strand=c[6], ident=float(a["Identity"]), t0=t0, t1=t1, exons=[])
        else:
            mrna[(runname, a["Parent"])]["exons"].append((int(c[3]), int(c[4])))
    rejected = []
    for m in mrna.values():
        ex = sorted(m["exons"]); cds = "".join(seq[a-1:b] for a, b in ex)
        if m["strand"] == "-": cds = str(Seq(cds).reverse_complement())
        m["cds"] = cds; m["prot"] = str(Seq(cds[:len(cds)//3*3]).translate()); m["len"] = len(seq)
        full = [v for k, v in REFLEN.items() if k.startswith(m["ref"]) and (f"name={m['name']}|" in k + "|" or k.endswith(m["name"]))]
        m["reflen"] = full[0] if full else None
        m["clean"] = not m["frameshift"] and m["prot"].rstrip("*").count("*") == 0 and len(cds) % 3 == 0
        if not m["clean"]: rejected.append(m); continue
        k = m["name"]
        if k not in out or m["score"] > out[k]["score"] or (m["score"] == out[k]["score"] and m["ref"] < out[k]["ref"]): out[k] = m
    out["_rejected"] = rejected
    return out

if __name__ == "__main__":
    print("window\tgene\tstrand\tstart\tend\texons\taa\tinternal_stop\tM_start\tstop\tref\tref_range\tref_len\tident\tscore\tedge_bp")
    for f in sorted(glob.glob("data/*.miniprot.gff")):
        slug = f.split("/")[1].split(".miniprot")[0]
        for g, m in sorted(((g, m) for g, m in models(slug).items() if g != "_rejected"), key=lambda x: x[1]["start"]):
            p = m["prot"]; ed = min(m["start"] - 1, m["len"] - m["end"])
            print(slug.split("__")[0] + ":" + slug.split("__")[1].split("_")[0], g, m["strand"], m["start"], m["end"], len(m["exons"]), len(p.rstrip("*")), p.rstrip("*").count("*"), p[0] == "M", p.endswith("*"),
                  m["ref"][:8], f"{m['t0']}-{m['t1']}", m["reflen"], m["ident"], m["score"], ed, sep="\t")
