#!/bin/bash
P=/bigdata/stajichlab/jstajich/projects/MATPredict
I=$P/results/sra_reads_typing/wave_inputs; W=$P/results/sra_reads_typing/waves/wave1/plan
S=$P/results/2026-10-26_v070_submit
export PYTHONPATH=$P/.claude/worktrees/run-b317ece/src MATPREDICT_DB_ROOT=$P/.claude/worktrees/run-b317ece/db
PY=$P/.pixi/envs/default/bin/python
while read slug pfx asm; do
  echo "=== $slug $asm"
  f=$(ls $I/genomes/$asm.fa.gz $I/genomes/$asm.fasta.gz 2>/dev/null | head -1); [ -z "$f" ] && f=$(ls $I/genomes/${asm}* | head -1)
  $PY $S/assembly_panel_check.py --panel $W/$slug/panel --prefix $pfx --assembly $f
done <<LIST
macrophomina_phaseolina macrophomina GCA_047665475.1_ASM4766547v1
fusarium_proliferatum fusarium GCA_058129285.1_ASM5812928v1
aureobasidium_melanogenum aureobasidium GCA_022814665.1_ASM2281466v1
ustilaginoidea_virens ustilaginoidea GCA_002939685.2_ASM293968v2
dothistroma_septosporum dothistroma GCA_002236565.2_GUA1.long
dothistroma_septosporum dothistroma GCA_002236745.2_GRE1.long
LIST
