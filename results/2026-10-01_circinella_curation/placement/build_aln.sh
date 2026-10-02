#!/bin/bash
# Tree set = label-tree tree_set.faa (105) + the two out-of-group gains (detect's
# scored sexP model). hmmalign to PF00505 (match columns only, 69) and MAFFT
# L-INS-i full length, columns with >= 50% gaps removed (same rules as the
# label tree). Tip names: '|' -> '__'.
set -euo pipefail
source /etc/profile.d/modules.sh
module load hmmer/3.4 mafft/7.505
LT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_circinella_label_tree
P=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
cat $LT/tree_set.faa gains.faa | sed 's/|/__/g' > tree_set_plus.faa
hmmalign --trim --outformat afa /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_sexMP_phylogeny/PF00505.hmm tree_set_plus.faa > aln_hmm_raw.afa
mafft --localpair --maxiterate 1000 --thread 4 --quiet tree_set_plus.faa > aln_mafft.afa
$P - <<'PY'
from Bio import AlignIO
r = AlignIO.read('aln_hmm_raw.afa', 'fasta')
k = [i for i in range(r.get_alignment_length()) if all(not (x.seq[i].islower() or x.seq[i] == '.') for x in r)]
with open('aln_hmm.afa', 'w') as o:
    for x in r: o.write(f'>{x.id}\n{"".join(x.seq[i] for i in k)}\n')
a = AlignIO.read('aln_mafft.afa', 'fasta'); n = len(a)
k = [i for i in range(a.get_alignment_length()) if sum(1 for x in a if x.seq[i] == '-') / n < 0.5]
with open('aln_mafft_g50.afa', 'w') as o:
    for x in a: o.write(f'>{x.id}\n{"".join(x.seq[i] for i in k).upper()}\n')
print('hmm', len(r), 'x', len(AlignIO.read("aln_hmm.afa","fasta")[0]), '; mafft', n, 'x', len(k))
PY
