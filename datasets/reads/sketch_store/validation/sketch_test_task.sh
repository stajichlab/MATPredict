#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 12G -t 2:00:00 -J afum_sketch -o logs/afum_sketch.%A_%a.log
# One strain of the Afumigatus_Global set: two independent ~30x subsamples (read pairs 1..N and N+1..2N) of the local FASTQs,
# each sketched with sourmash (k21,k31 scaled=100 abund; kept raw and filtered at >=2), and subsample A typed with MATPredict reads-type (L1 panel).
set -uo pipefail
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-07_sketch_pcoa_mat
IN=/bigdata/stajichlab/shared/projects/Population_Genomics/Afumigatus_Global/input
SM=/opt/linux/rocky/8.x/x86_64/pkgs/sourmash/4.9/bin/sourmash
P=/bigdata/stajichlab/jstajich/projects/MATPredict
PY=$P/.pixi/envs/default/bin/python
export PYTHONPATH=$P/.claude/worktrees/run-b317ece/src MATPREDICT_DB_ROOT=$P/.claude/worktrees/run-b317ece/db
PANEL=$P/results/sra_reads_typing/waves/wave1/plan/aspergillus_fumigatus/panel/L1
GENOME_KB=${GENOME_KB:-29400}; COV=${COV:-30}
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/strains.tsv)
S=$(cut -f1 <<< "$line"); PAT=$(cut -f2 <<< "$line")
OUT=$W/sigs; mkdir -p $OUT $W/mat $W/stats
[[ -s $OUT/$S.A.min2.sig.zip && -s $W/mat/$S.tsv ]] && { echo "skip $S"; exit 0; }
work=$(mktemp -d ${SCRATCH:-/tmp}/afsk.XXXXXX); trap 'rm -rf $work' EXIT
r1=(); r2=()
IFS=';' read -ra parts <<< "$PAT"
for p in "${parts[@]}"; do
  ex=( $IN/$p )      # the sample sheet uses shell patterns ([12], R?): expand, first file = read 1, second = read 2
  r1+=("${ex[0]}"); [[ ${#ex[@]} -ge 2 ]] && r2+=("${ex[1]}")
done
for f in "${r1[@]}" "${r2[@]}"; do [[ -s $f ]] || { echo "missing $f" >&2; exit 1; }; done
L=$(zcat -f "${r1[0]}" | head -2 | tail -1 | tr -d '\n' | wc -c)
if [[ ${#r2[@]} -gt 0 ]]; then N=$(( GENOME_KB * 1000 * COV / (2 * L) )) ; else N=$(( GENOME_KB * 1000 * COV / L )); fi
N=${N%.*}
t0=$(date +%s)
# subsample A = pairs 1..N, B = pairs N+1..2N; interleave pairs for reads-type
cut1() { zcat -f "$@" | awk -v n=$N -v a=$work/A.$k.fq -v b=$work/B.$k.fq 'NR<=4*n{print > a; next} NR<=8*n{print > b; next} {exit}'; }
k=1; cut1 "${r1[@]}"
if [[ ${#r2[@]} -gt 0 ]]; then k=2; cut1 "${r2[@]}"; fi
for h in A B; do
  [[ -s $work/$h.1.fq ]] || { echo "no reads left for subsample $h (library under 2x${COV}x)"; rm -f $OUT/$S.$h.min2.sig.zip; continue; }
  files=($work/$h.1.fq); [[ -s $work/$h.2.fq ]] && files+=($work/$h.2.fq)
  bases=$(awk 'NR%4==2{s+=length($0)}END{print s+0}' "${files[@]}")
  $SM sketch dna -p k=21,k=31,scaled=100,abund --name "$S.$h" -o $work/$S.$h.raw.sig.zip "${files[@]}" 2>/dev/null
  $SM sig filter -m 2 $work/$S.$h.raw.sig.zip -o $OUT/$S.$h.min2.sig.zip 2>/dev/null
  echo -e "$S\t$h\t$L\t$N\t$bases\t$(stat -c %s $work/$S.$h.raw.sig.zip)\t$(stat -c %s $OUT/$S.$h.min2.sig.zip)" > $W/stats/$S.$h.tsv
done
t1=$(date +%s)
# MAT typing of subsample A (interleaved) with the A. fumigatus lineage panel
if [[ -s $work/A.2.fq ]]; then paste <(paste - - - - < $work/A.1.fq) <(paste - - - - < $work/A.2.fq) | tr '\t' '\n' > $work/A.fq; else mv $work/A.1.fq $work/A.fq; fi
printf 'sample\treads\n%s\t%s\n' "$S" "$work/A.fq" > $work/samples.tsv
$PY -m MATPredict reads-type --idiomorph MAT1-1=$PANEL/aspergillus_MAT1-1.fasta --idiomorph MAT1-2=$PANEL/aspergillus_MAT1-2.fasta --samples $work/samples.tsv --out $W/mat/$S.tsv 2>$work/type.log || tail -3 $work/type.log
echo "$S done: sketch $((t1-t0)) s, type $(( $(date +%s) - t1 )) s"
