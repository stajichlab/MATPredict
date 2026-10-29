#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; I=$P/results/sra_reads_typing/wave_inputs
O=$P/results/2026-10-28_taphrinomycotina_mat_candidates; E=$O/e2e; mkdir -p $E
SRC=$P/.claude/worktrees/run-b317ece
export PYTHONPATH=$SRC/src
export PATH=$P/.pixi/envs/default/bin:$PATH
PY=$P/.pixi/envs/default/bin/python
# 1. database copy with close representatives appended to the Taphrina PM record (matMc, matPi; Mi too short)
if [ ! -d $E/db_taph ]; then cp -r $SRC/db $E/db_taph; DOAPPEND=1; fi
[ -n "$DOAPPEND" ] && $PY - <<PY2
import csv,sys
sys.path.insert(0,"$P/results/2026-10-26_v070_submit")
import taph_mat_candidates as t
O="$O"; seqs={k.split()[0]:v for k,v in t.read_fasta(O+"/candidates.faa")}
GENE={"matPi":"matPi","mat1-Pi":"matPi","matMc":"matMc","mat1-Mc":"matMc"}
add=[];seen=set()
for r in csv.DictReader(open(O+"/candidates.tsv"),delimiter="\t"):
    g=GENE.get(r["ref_gene"]); sp=r["species"]
    if r["ref_class"]!="PM" or not g or sp.startswith(("Schizosaccharomyces pombe","Neolecta","Saitoella")): continue
    if float(r["identity"])<0.40: continue
    q=seqs[r["candidate"]].rstrip("*")
    if q in seen or len(q)<100: continue
    seen.add(q); add.append((g,sp,q))
f="$E/db_taph/Ascomycota/Taphrinales/5011_unknown-1_PM_combined/proteins.faa"
with open(f,"a") as fh:
    for i,(g,sp,q) in enumerate(add,10): fh.write(f">5011_unknown-1_PM_combined|gene_index={i}|name={g}|role=core_MAT\n{q}\n")
print("appended",len(add),"proteins:",{g:sum(1 for a in add if a[0]==g) for g in ("matMc","matPi")},sorted({a[1] for a in add}))
PY2
# 2. detect on target genomes with both databases
for g in GCA_001929475.1_NeoIirr1.0 GCA_051599425.1_Sai_col_YB2330_v1 GCF_001661265.1_Saico1; do
  gunzip -c $I/genomes/$g.fa.gz > $E/$g.fna
  for db in base taph; do
    D=$SRC/db; [ $db = taph ] && D=$E/db_taph
    MATPREDICT_DB_ROOT=$D $PY -m MATPredict detect --genome $E/$g.fna --phylum Ascomycota --out-dir $E/${g}_$db > $E/${g}_$db.log 2>&1
    echo "$g $db exit=$?"
  done
done
