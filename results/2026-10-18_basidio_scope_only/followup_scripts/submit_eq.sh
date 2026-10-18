#!/usr/bin/bash -l
cd /bigdata/stajichlab/jstajich/projects/mat-fam18
P=lists
A1="5003_cbs-491_redPR_A1 5005_cbs-483_redPR_A1 5286_cbs-14_redPR_A1 5535_cbs-20_redPR_A1 5537_jy1105_redPR_A1 86836_pycc-4818_redPR_A1"
A2="217517_cbs-9654_redPR_A2 29898_jj10-1_redPR_A2 5003_cbs-499_redPR_A2 5005_cbs-490_redPR_A2 5286_nbrc-0880_redPR_A2 86836_cbs-4785_redPR_A2 86839_nrrl-y-7192_redPR_A2"
excl() { for r in $A1 $A2; do case " $* " in *" $r "*) ;; *) printf "%s," $r;; esac; done | sed 's/,$//'; }
# sample of 80 genomes (seed 1) into 10 chunks
cat lists/sporidio.* | sort > lists/sporidio_all.txt
python3 - <<'E'
import random
L = open("lists/sporidio_all.txt").read().split("\n"); L = [l for l in L if l]
random.Random(1).shuffle(L)
S = sorted(L[:80])
for i in range(10):
    open(f"lists/samp80.{i:02d}", "w").write("".join(x + "\n" for x in S[i::10]))
E
# all 251 into 10 chunks
python3 - <<'E'
L = sorted(l for l in open("lists/sporidio_all.txt").read().split("\n") if l)
for i in range(10):
    open(f"lists/all251.{i:02d}", "w").write("".join(x + "\n" for x in L[i::10]))
E
E1=$(excl 5286_cbs-14_redPR_A1 5286_nbrc-0880_redPR_A2 5005_cbs-483_redPR_A1 217517_cbs-9654_redPR_A2)
E3=$(excl $A1 217517_cbs-9654_redPR_A2 5286_nbrc-0880_redPR_A2)
E4=$(excl $A2 5286_cbs-14_redPR_A1 5005_cbs-483_redPR_A1)
echo "E1 $E1"; echo "E3 $E3"; echo "E4 $E4"
a=$(EXCLUDE="$E1" sbatch --parsable --export=ALL --array=0-9%10 -o logs/e1.%A_%a.log run_ex.slurm lists/all251 runs_eq/E1)
b=$(EXCLUDE="$E3" sbatch --parsable --export=ALL --array=0-9%10 -d afterany:$a -o logs/e3.%A_%a.log run_ex.slurm lists/samp80 runs_eq/E3)
c=$(EXCLUDE="$E4" sbatch --parsable --export=ALL --array=0-9%10 -d afterany:$b -o logs/e4.%A_%a.log run_ex.slurm lists/samp80 runs_eq/E4)
echo $a $b $c
