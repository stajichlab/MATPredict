#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; R=$P/results/sra_reads_typing; S=$P/results/2026-10-26_v070_submit
export PYTHONPATH=$P/.claude/worktrees/run-b317ece/src MATPREDICT_DB_ROOT=$P/.claude/worktrees/run-b317ece/db PATH=$P/.pixi/envs/default/bin:$PATH
python $S/assembly_panel_check.py --panel $R/waves/wave2/plan/botrytis_cinerea/panel --prefix botrytis --assembly $R/wave_inputs/genomes/GCF_000143535.2_ASM14353v4.fa.gz
echo "== all 3 B05.10 reads results vs other MAT1-1 runs of the species (breadth / depth / shared)"
python3 - <<PY
import csv,glob
W="$R/waves/wave2/out"
rows=[]
for f in glob.glob(W+"/*.botrytis.tsv"):
    r=list(csv.DictReader(open(f).read().replace("\r","").splitlines(),delimiter="\t"))[0]; rows.append((r["sample"],r["call"],r["MAT1-1_breadth"],r["MAT1-1_depth"],r["MAT1-2_breadth"],r["MAT1-2_depth"],r["shared_depth"]))
print("MAT1-1 calls: median MAT1-1 breadth %.2f depth %.1f"%(sorted(float(x[2]) for x in rows if x[1]=="MAT1-1")[len([x for x in rows if x[1]=="MAT1-1"])//2],sorted(float(x[3]) for x in rows if x[1]=="MAT1-1")[len([x for x in rows if x[1]=="MAT1-1"])//2]))
for x in rows:
    if x[0] in ("SRR29713114","SRR29713115","SRR29713116"): print(x)
PY
