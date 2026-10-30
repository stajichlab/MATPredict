"""Helpers for hand-curating gene models on a locus window: genome coordinates, splice-site and protein checks."""
from Bio import SeqIO
from Bio.Seq import Seq
OFFSET = {}
def load(fa):
    r = next(SeqIO.parse(fa, "fasta")); return str(r.seq).upper()
def sub(S, a, b, strand="+"):
    """1-based inclusive window coordinates -> sequence on the given strand."""
    x = S[a - 1:b]
    return x if strand == "+" else str(Seq(x).reverse_complement())
def cds(S, exons, strand):
    """exons: list of (a,b) window coords; returns the spliced CDS in transcript orientation."""
    ex = sorted(exons) if strand == "+" else sorted(exons, reverse=True)
    return "".join(sub(S, a, b, strand) for a, b in ex)
def protein(S, exons, strand):
    return str(Seq(cds(S, exons, strand)).translate())
def introns(exons, strand):
    ex = sorted(exons); return [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
def splice_sites(S, exons, strand):
    out = []
    for a, b in introns(exons, strand):
        t = sub(S, a, b, strand); out.append((a, b, b - a + 1, t[:2] + ".." + t[-2:], t[-45:-1] if False else ""))
    return out
def find_orfs_around(S, a, b, strand):
    return None
