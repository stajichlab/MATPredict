"""Does the open reading frame continue past a model's end? python orf_check.py RUNS LIB GENOME_LIST OUT_PREFIX
For each genome's MAT1-1-1 evidence: translate (a) the canonical model, (b) the model extended in frame to the next stop, (c) the alternate (miniprot) model;
write the loci window FASTA for tblastn."""
import gzip, sys
from pathlib import Path
import yaml

runs, lib, lst, outp = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
CODON = {}
bases = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(aas):
    CODON[bases[i // 16] + bases[(i // 4) % 4] + bases[i % 4]] = a
comp = str.maketrans("ACGTacgtN", "TGCAtgcaN")
rc = lambda s: s.translate(comp)[::-1]
tr = lambda s: "".join(CODON.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def contig_seq(path, name):
    seq, on = [], False
    with gzip.open(path, "rt") as f:
        for l in f:
            if l[0] == ">":
                if on:
                    break
                on = l[1:].split()[0] == name
            elif on:
                seq.append(l.strip())
    return "".join(seq)


def cds(seq, exons, strand):
    s = "".join(seq[a - 1:b] for a, b in sorted(exons))
    return s if strand == "+" else rc(s)


win = open(outp + "_windows.fa", "w")
for line in open(lst):
    grp, g = line.split()
    d = yaml.load(open(runs / g / "detection_report.yaml"), Loader=L)
    for loc in d["detected"]:
        for e in loc["gene_evidence"]:
            if e["gene"] != "MAT1-1-1" or e["role"] != "core_MAT" or not e.get("exons"):
                continue
            seq = contig_seq(lib / f"{g}.fa.gz", e["contig"])
            ex = [(a["start"], a["end"]) for a in e["exons"]]
            st = e["strand"]
            c = cds(seq, ex, st)
            p = tr(c)
            # extend 3' in frame
            if st == "+":
                tail = seq[max(b for _, b in ex):max(b for _, b in ex) + 2400]
            else:
                tail = rc(seq[max(0, min(a for a, _ in ex) - 1 - 2400):min(a for a, _ in ex) - 1])
            ext = tr(tail)
            stop = ext.find("*")
            alt = e.get("alternate_model") or {}
            ap = tr(cds(seq, [(a["start"], a["end"]) for a in alt["exons"]], alt["strand"])) if alt else ""
            print(f"== {grp} {g} {e['contig']}:{min(a for a,_ in ex)}-{max(b for _,b in ex)}{st}")
            print(f"  canonical model {len(p)} aa; ends with stop in model: {p.endswith('*')}; internal stops {p[:-1].count('*')}")
            print(f"  in-frame extension past the model end: {stop if stop >= 0 else '>800'} aa to the next stop ({'stop found' if stop >= 0 else 'no stop within 2.4 kb'})")
            print(f"  extension: {ext[:max(stop, 0) if stop >= 0 else 120]}")
            print(f"  miniprot alt: {len(ap)} aa, internal stops {ap[:-1].count('*')}, exons {len(alt.get('exons', []))}, identity {alt.get('identity')}")
            lo, hi = min(a for a, _ in ex), max(b for _, b in ex)
            w0, w1 = max(1, lo - 3000), min(len(seq), hi + 3000)
            win.write(f">{grp}|{g}|{e['contig']}:{w0}-{w1}\n{seq[w0 - 1:w1]}\n")
            open(outp + f"_{grp}_{g[:15]}_models.faa", "w").write(f">canon\n{p}\n>ext\n{p.rstrip('*')}{ext[:stop] if stop >= 0 else ''}\n>miniprot\n{ap}\n")
