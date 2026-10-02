"""Extract Rhodotorula MAT genes (STE3, STE20, HD1, HD2) from the group's
funannotate/miniprot GFFs, translated from the PUBLIC contig sequence."""
import json, re, sys
from Bio import SeqIO
from Bio.Seq import Seq

G = "/bigdata/stajichlab/ctsai085/projects/Rhodotorula_comparative_genomics/synteny_map/data/gff"
N = "/scratch/jstajich/28984668/claude-1181/-bigdata-stajichlab-jstajich-projects-MATPredict/48b08d1a-4f03-410b-bf7c-431cd980f588/scratchpad/rhod/ncbi"

STRAINS = {
    # key: (gff stem, {local contig -> public accession})
    "CBS14": ("Rhodotorula_toruloides_CBS14", {}),
    "NBRC_0880": ("Rhodotorula_toruloides_NBRC_0880", {}),
    "CBS_20": ("Rhodotorula_glutinis_CBS_20", {}),
    "JY1105": ("Rhodotorula_alborubescens_JY1105", {}),
    "Y-7192": ("Rhodotorula_sphaerocarpa_NRRL_Y-7192",
               {"scaffold_13": "JBRFVS010000013.1", "scaffold_5": "JBRFVS010000005.1"}),
    "LS11": ("Rhodotorula_kratochvilovae_LS11", {}),
    "JJ10.1": ("Rhodotorula_sp._clade_I_JJ10.1", {}),
}
WANT = {"PR": re.compile(r"^(STE3_A[12]|STE20)$"), "HD": re.compile(r"^(HD1|HD2)$")}


def contig_seq(acc):
    return str(next(SeqIO.parse(f"{N}/{acc}.fa", "fasta")).seq).upper()


def parse(path):
    genes, cds = {}, {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            genes[attrs["ID"]] = dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                                      name=attrs.get("Name", ""), attrs=attrs)
        elif f[2] == "CDS":
            parent = attrs["Parent"].rsplit("-T", 1)[0]
            cds.setdefault(parent, []).append((int(f[3]), int(f[4]), int(f[7]) if f[7] != "." else 0))
    return genes, cds


out = []
for key, (stem, cmap) in STRAINS.items():
    for locus in ("PR", "HD"):
        genes, cds = parse(f"{G}/{locus}/{stem}.gff")
        for gid, g in genes.items():
            if not WANT[locus].match(g["name"]):
                continue
            acc = cmap.get(g["contig"], g["contig"] if "." in g["contig"] else g["contig"] + ".1")
            seq = contig_seq(acc)
            ex = sorted(cds.get(gid, []))
            if not ex:
                continue
            nt = "".join(seq[a - 1:b] for a, b, _ in ex)
            if g["strand"] == "-":
                nt = str(Seq(nt).reverse_complement())
            first_phase = ex[-1][2] if g["strand"] == "-" else ex[0][2]
            prot = str(Seq(nt[first_phase:len(nt) - (len(nt) - first_phase) % 3]).translate())
            out.append(dict(strain=key, locus=locus, gene=g["name"], id=gid, contig=acc,
                            start=g["start"], end=g["end"], strand=g["strand"],
                            exons=[[a, b] for a, b, _ in ex], codon_start=first_phase + 1,
                            aa_len=len(prot.rstrip("*")), starts_M=prot.startswith("M"),
                            internal_stops=prot.rstrip("*").count("*"), ends_stop=prot.endswith("*"),
                            protein=prot.rstrip("*"), attrs=g["attrs"]))
json.dump(out, open(sys.argv[1], "w"), indent=1)
for o in out:
    print(o["strain"], o["locus"], o["gene"], o["id"], o["contig"], o["start"], o["end"], o["strand"],
          "exons", len(o["exons"]), "aa", o["aa_len"], "M", o["starts_M"], "stops", o["internal_stops"],
          "endstop", o["ends_stop"], "cs", o["codon_start"])
