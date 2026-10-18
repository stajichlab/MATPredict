"""Extract the HMG-box clan (CL0114) and the MAT-specific Pfam models from Pfam-A (current = 38.2, 2026-01-27)."""
import sys
import pyhmmer
D = "/srv/projects/db/pfam/current/"
want = {"PF04769", "PF17043", "PF28944", "PF04611"}
acc = clan = None
for line in open(D + "Pfam-A.hmm.dat"):
    if line.startswith("#=GF AC"):
        acc = line.split()[2].split(".")[0]
    elif line.startswith("#=GF CL") and line.split()[2] == "CL0114":
        want.add(acc)
print("wanted", sorted(want), file=sys.stderr)
n = 0
with pyhmmer.plan7.HMMFile(D + "Pfam-A.hmm") as fh, open(sys.argv[1], "wb") as out:
    for hmm in fh:
        a = hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession
        if a.split(".")[0] in want:
            hmm.write(out); n += 1
print("wrote", n, "models", file=sys.stderr)
