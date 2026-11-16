#!/usr/bin/bash
#SBATCH --job-name=hd-tree
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/hdtree_%A.log
# Homeodomain tree. Env: SRC, BASE (results/2026-11-13_flank_first_scan), OUT.
set -uo pipefail
: "${SRC:?}" "${BASE:?}" "${OUT:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load hmmer/3.4
module load iqtree/3.0.1
mkdir -p "$OUT"
hmmfetch "$BASE/pfam/hmg_homeo.hmm" Homeodomain > "$OUT/Homeodomain.hmm"
python3 "$SRC/scripts/homeodomain_tree.py" --base "$BASE" --hmm "$OUT/Homeodomain.hmm" --out "$OUT" || exit 1
hmmalign --trim --amino -o "$OUT/hd.sto" "$OUT/Homeodomain.hmm" "$OUT/homeodomain_domains.faa"
esl-reformat afa "$OUT/hd.sto" > "$OUT/hd_aln.afa"
python3 - "$OUT" <<PY
import sys
o=sys.argv[1]; seqs={}; k=None
for l in open(o+"/hd_aln.afa"):
    if l[0]==">": k=l[1:].strip(); seqs[k]=[]
    else: seqs[k].append(l.strip())
seqs={k:"".join(v).upper().replace(".","-") for k,v in seqs.items()}
L=len(next(iter(seqs.values())))
keep=[i for i in range(L) if sum(s[i]!="-" for s in seqs.values())>=0.5*len(seqs)]
with open(o+"/hd_trim.afa","w") as f:
    for k,s in seqs.items():
        t="".join(s[i] for i in keep)
        if sum(c!="-" for c in t)>=0.6*len(keep): f.write(f">{k}\n{t}\n")
print("alignment columns",L,"kept",len(keep),"sequences",len(seqs))
PY
cd "$OUT"
iqtree3 -s hd_trim.afa -m MFP -B 1000 -T 8 --prefix hd_tree -redo > hd_tree.stdout 2>&1 || iqtree -s hd_trim.afa -m MFP -B 1000 -T 8 --prefix hd_tree -redo > hd_tree.stdout 2>&1
echo done
