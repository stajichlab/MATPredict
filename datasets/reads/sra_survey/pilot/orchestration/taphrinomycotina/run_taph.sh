#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; I=$P/results/sra_reads_typing/wave_inputs
O=$P/results/2026-10-28_taphrinomycotina_mat_candidates; mkdir -p $O
export PATH=$P/.pixi/envs/default/bin:$PATH
python3 - <<PY
import csv
G=list(csv.DictReader(open("$I/asco_genomes.tsv"),delimiter="\t"))
keep=[r for r in G if r["order"] in ("Taphrinales","Schizosaccharomycetales","Pneumocystidales","Neolectales","Archaeorhizomycetales","Protomycetales","Saitoellales") or r["class_"] in ("Taphrinomycetes","Schizosaccharomycetes","Pneumocystidomycetes","Neolectomycetes","Saitoellomycetes") or r["genus"] in ("Saitoella","Taphrina","Protomyces","Neolecta","Schizosaccharomyces","Pneumocystis","Archaeorhizomyces","Savitreella")]
w=csv.writer(open("$O/genomes.tsv","w"),delimiter="\t",lineterminator="\n"); w.writerow(["genome","species","class_","order","status"])
for r in keep: w.writerow([r["genome"],r["species"],r["class_"],r["order"],r["status"]])
print(len(keep),"genomes")
PY
python3 $P/results/2026-10-26_v070_submit/taph_mat_candidates.py --db $P/.claude/worktrees/run-b317ece/db --genomes $O/genomes.tsv --genome-dir $I/genomes --out $O --reports-glob "$P/results/2026-10-26_v070_ascomycota/wave_*/small_*/runs/GENOME/detection_report.yaml" --threads 4 > $O/run.log 2>&1
