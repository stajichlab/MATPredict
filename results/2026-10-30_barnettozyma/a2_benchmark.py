"""How strong is the MATa2 call? Score the Saccharomycotina MATa2 profile against (i) its own training proteins and the P. kudriavzevii MATa2,
(ii) every stop-to-stop 6-frame segment >= 80 aa of the whole PB4207 window (66 kb, ~ 20 genes), (iii) the B. californica MATa2 model."""
import pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq
from pyhmmer.easel import Alphabet, TextSequence, SequenceFile
def nm(x): return x.decode() if isinstance(x, bytes) else x
alpha = Alphabet.amino()
with pyhmmer.plan7.HMMFile("/Users/jstajich/projects/MATPredict/db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm") as hf: hmm = hf.read()
def score(seqs):
    return {nm(h.name): h.score for hits in pyhmmer.hmmsearch(hmm, seqs, E=1e3) for h in hits}
with SequenceFile("a2_refs.faa", digital=True, alphabet=alpha) as sf: refs = sf.read_block()
print("curated MATa2 proteins (training + Pk):", {k: round(v) for k, v in score(refs).items()})
S = str(next(SeqIO.parse("data/GCA_023627795.1_JAMAIE010000005.1_589225-655853.fa", "fasta")).seq).upper()
segs = []
for strand in "+-":
    s = S if strand == "+" else str(Seq(S).reverse_complement())
    for fr in range(3):
        aa = str(Seq(s[fr:len(s) - ((len(s) - fr) % 3)]).translate())
        for part in aa.split("*"):
            if len(part) >= 80: segs.append(part)
sc = score([TextSequence(name=str(i).encode(), sequence=p).digitize(alpha) for i, p in enumerate(segs)])
vals = sorted(sc.values(), reverse=True)
print(f"{len(segs)} segments >= 80 aa in the 66.6 kb window; best profile scores:", [round(v) for v in vals[:6]], "median", round(vals[len(vals)//2]) if vals else None)
seqs = [x for x in SeqIO.parse("models_proteins.faa", "fasta")]
m = [TextSequence(name=x.id.encode(), sequence=str(x.seq)).digitize(alpha) for x in seqs]
print("models:", {k: round(v) for k, v in score(m).items()})
