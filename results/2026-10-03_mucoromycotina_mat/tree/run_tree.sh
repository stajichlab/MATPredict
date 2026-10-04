#!/usr/bin/bash -l
#SBATCH -p epyc,batch -c 16 --mem 32gb -t 12:00:00 -J mm_sexMP_tree
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree/logs/%x_%j.log
# sexM/sexP gene trees for the Mucoromycotina MAT campaign (2026-10-03).
# 1. cd-hit 0.98 within each idiomorph (redundancy), cd-hit 0.90 on the outgroup.
# 2. Two alignments: (a) HMG box only (hmmalign to PF00505, match columns);
#    (b) full length (MAFFT L-INS-i, ClipKIT kpic-smart-gap).
# 3. IQ-TREE 3: ModelFinder (protein), UFBoot 1000, SH-aLRT 1000, on each.
# 4. RAxML-NG check on the HMG-box alignment (best model from IQ-TREE, 200 bootstraps).
set -eo pipefail
T=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree
PF=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_sexMP_fasttree/PF00505.hmm
# module activation scripts reference unset variables: load before `set -u`
module load cd-hit/4.8.1 hmmer/3.4 mafft/7.505 clipkit/1.3.0 iqtree/3.0.1 raxml-ng/2.0.2
set -u
W=${SCRATCH:?}/tree; mkdir -p $W; cd $W
N=${SLURM_CPUS_PER_TASK:-8}
cd-hit -i $T/ingroup_sexP.faa -o sexP.nr.faa -c 0.98 -n 5 -T $N -M 0 -d 0 > cdhit_P.log
cd-hit -i $T/ingroup_sexM.faa -o sexM.nr.faa -c 0.98 -n 5 -T $N -M 0 -d 0 > cdhit_M.log
cd-hit -i $T/outgroup.faa -o out.nr.faa -c 0.90 -n 5 -T $N -M 0 -d 0 > cdhit_O.log
cat sexP.nr.faa sexM.nr.faa out.nr.faa > all.faa
cp sexP.nr.faa.clstr sexM.nr.faa.clstr out.nr.faa.clstr all.faa $T/
# (a) HMG box
hmmalign --trim --outformat afa $PF all.faa > hmg_raw.afa
# keep match columns only (upper case and '-'), drop insert columns (lower case and '.')
python3 - <<'EOF'
seqs, n = {}, None
for l in open("hmg_raw.afa"):
    l = l.rstrip()
    if l.startswith(">"): n = l[1:].split()[0]; seqs[n] = []
    else: seqs[n].append(l)
with open("hmg.afa", "w") as fo:
    for k, v in seqs.items():
        s = "".join(c for c in "".join(v) if not (c.islower() or c == "."))
        if sum(c != "-" for c in s) >= 40: fo.write(f">{k}\n{s}\n")
EOF
# (b) full length
mafft --localpair --maxiterate 1000 --thread $N all.faa > full_raw.afa 2> mafft.log
clipkit full_raw.afa -m kpic-smart-gap -o full.afa > clipkit.log
cp hmg_raw.afa hmg.afa full_raw.afa full.afa $T/
for a in hmg full; do
  iqtree3 -s $a.afa --prefix iq_$a -T $N -m MFP -mset LG,WAG,JTT,VT,Q.pfam -B 1000 -alrt 1000 -seed 20261003 > /dev/null
  cp iq_$a.* $T/
done
M=$(grep "Best-fit model" iq_hmg.iqtree | sed 's/.*: //; s/ .*//')
raxml-ng --all --msa hmg.afa --model "$M" --bs-trees 200 --threads $N --seed 20261003 --prefix rx_hmg > /dev/null || echo "raxml-ng failed (model $M)"
cp rx_hmg.* $T/ 2>/dev/null || true
echo done
