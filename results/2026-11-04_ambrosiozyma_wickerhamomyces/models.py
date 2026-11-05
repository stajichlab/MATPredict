"""Gene models of the two candidate genomes: the miniprot projection of a curated reference protein, chosen per gene from the hits listed by
summarize_miniprot.py (window coordinates, strand and reference prefix below), checked for frameshift, internal stops and CDS length.
Our own projection (tier 3): no RNA-seq and no annotation. Coordinates are written in contig coordinates (window start - 1 added)."""
import glob, re
from Bio.Seq import Seq
# window prefix -> contig, [(gene, role, strand, start, end, reference-id prefix)]
SEL = {
 "edaphicus__JBSXOT010000020.1": [("sla2", "flanking_conserved", "+", 25657, 28803, "460519_cbs-7435_MATkp_MATa"),
                                  ("MATalpha1", "core_MAT", "-", 29060, 29578, "460523_cbs-4732_MATog"),
                                  ("MATalpha2", "core_MAT", "-", 30671, 31584, "4932_s288c_MATsc_MATalpha"),
                                  ("dic1", "flanking_conserved", "+", 32772, 33608, "1005962_dl-1_MATog")],
 "edaphicus__JBSXOT010000274.1": [("sla2", "flanking_conserved", "+", 4, 537, "1173061_clib918_MATdc_A"),
                                  ("MATa1", "core_MAT", "-", 624, 840, "4932_s288c_MATsc_HMRa"),
                                  ("MATa2", "core_MAT", "+", 1620, 2444, "1173061_clib918_MATdc_A")],
 "philentoma__PPKZ03000004.1": [("dic1", "flanking_conserved", "-", 12586, 13407, "460523_cbs-4732_MATog"),
                                ("MATa2", "core_MAT", "-", 14686, 15171, "5476_sc5314_MTL_A_MTLA2"),
                                ("MATa1", "core_MAT", "+", 19030, 19203, "460523_cbs-4732_MATog"),
                                ("sla2", "flanking_conserved", "-", 19453, 22578, "MH260578_SLA2")],
 "philentoma__PPKZ03000012.1": [("MATalpha2", "core_MAT", "+", 14935, 15255, "460523_cbs-4732_MATog"),
                                ("MATalpha1", "core_MAT", "+", 19783, 20226, "1005962_dl-1_MATog")],
}
def genome(fa): return "".join(l.strip() for l in open(fa) if l[0] != ">")
def build(prefix):
    fa = glob.glob(f"data/{prefix}_*.fa")[0]; seq = genome(fa); off = int(fa.rsplit("_", 1)[1].split("-")[0]) - 1
    contig = fa.split("__")[1].rsplit("_", 1)[0]
    mr = {}
    for fn in (fa.replace(".fa", ".miniprot.gff"), fa.replace(".fa", ".relaxed.gff")):
        for l in open(fn):
            c = l.rstrip().split("\t")
            if len(c) < 9 or c[2] not in ("mRNA", "CDS"): continue
            a = dict(kv.split("=", 1) for kv in c[8].split(";"))
            key = (fn, a.get("ID") or a["Parent"])
            if c[2] == "mRNA": mr[key] = dict(s=int(c[3]), e=int(c[4]), st=c[6], tgt=a["Target"], fs="Frameshift" in a, ex=[], id=float(a["Identity"]), score=int(c[5]))
            else: mr[key]["ex"].append((int(c[3]), int(c[4])))
    out = []
    for gene, role, strand, s, e, ref in SEL[prefix]:
        m = [v for v in mr.values() if v["s"] == s and v["e"] == e and v["st"] == strand and v["tgt"].startswith(ref)]
        if not m: out.append((gene, None)); continue
        v = max(m, key=lambda x: x["score"]); ex = sorted(v["ex"]); cds = "".join(seq[a - 1:b] for a, b in ex)
        if strand == "-": cds = str(Seq(cds).reverse_complement())
        p = str(Seq(cds[:len(cds) // 3 * 3]).translate())
        out.append((gene, dict(role=role, strand=strand, exons=[(a + off, b + off) for a, b in ex], aa=len(p.rstrip("*")), stops=p.rstrip("*").count("*"), prot=p.rstrip("*"),
                                fs=v["fs"], mod3=len(cds) % 3, ident=v["id"], ref=ref, contig=contig, tgt=v["tgt"].split()[1:3])))
    return out
if __name__ == "__main__":
    for p in SEL:
        print("==", p)
        for g, m in build(p):
            print("  ", g, "MISSING" if m is None else (m["strand"], m["exons"][0][0], m["exons"][-1][1], len(m["exons"]), "aa", m["aa"], "stops", m["stops"], "fs", m["fs"], "mod3", m["mod3"], "id %.2f" % m["ident"], m["ref"][:22], m["tgt"]))
