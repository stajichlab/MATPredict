#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; R=$P/results/sra_reads_typing; O=$R/waves/wave2/botrytis/b05map
export PATH=$P/.pixi/envs/default/bin:$PATH
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/minimap2/2.30/bin:$PATH; which minimap2
cd $O; [ -s b05.fa ] || gunzip -c $R/wave_inputs/genomes/GCF_000143535.2_ASM14353v4.fa.gz > b05.fa
minimap2 -t 4 -ax sr b05.fa SRR29713116.fq 2>/dev/null | samtools sort -@2 -o wt.bam - && samtools index wt.bam
samtools depth -a -r NC_037310.1:800000-826000 wt.bam > locus.depth
samtools depth -a wt.bam | awk '{s+=$3;n++} END{printf "genome mean depth %.1f over %d positions\n",s/n,n}'
python3 - <<PY
import statistics
d=[(int(l.split()[1]),int(l.split()[2])) for l in open("locus.depth")]
def m(a,b): v=[x for p,x in d if a<=p<=b]; return round(statistics.mean(v),1) if v else None
print("left flank 806478-812477:",m(806478,812477)," MAT1-1 core 812478-817230:",m(812478,817230)," right flank 817231-821303:",m(817231,821303))
for a in range(808000,820000,1000): print(a,a+999,m(a,a+999))
PY
