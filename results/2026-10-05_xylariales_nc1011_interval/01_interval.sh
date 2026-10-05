#!/bin/bash
#SBATCH -p exfab -c 8 --mem 16G -t 02:00:00 -J xyl-interval
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_interval_%j.log
# Annotation, HMM and tBLASTn survey of the NC1011 (GCA_022453505.1) SLA2-APN2 block.
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
set +eu; source /etc/profile 2>/dev/null; module load ncbi-blast/2.16.0+ hmmer/3.4 samtools/1.19.2 emboss/6.6.0; set -eu
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
cd "$D"; mkdir -p work
W=$D/work

# genome region
zcat "$GENOME_GZ" > "$SCRATCH/genome.fa"; samtools faidx "$SCRATCH/genome.fa"
samtools faidx "$SCRATCH/genome.fa" "$CONTIG:$REG_START-$REG_END" > "$W/region.fa"

# reference MAT proteins from the curated database (core_MAT genes, all Ascomycota records)
cat /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/db/Ascomycota/*/*/proteins.faa \
  | awk '/^>/{keep=($0 ~ /role=core_MAT/)} keep' > "$W/mat_refs.faa"
echo "mat refs: $(grep -c '>' "$W/mat_refs.faa")"

# HMG-related HMMs from Pfam: HMG_box PF00505, HMG_box_2 PF09011, MATalpha_HMGbox PF04769, PF17043
awk 'BEGIN{RS="//\n";ORS="//\n"} /ACC +PF(00505|09011|04769|17043)\./' "$PFAM" > "$D/data/hmg_related.hmm"
hmmpress -f "$D/data/hmg_related.hmm" >/dev/null
grep -E '^(NAME|ACC)' "$D/data/hmg_related.hmm" | paste - - > "$W/hmg_related_profiles.txt"; cat "$W/hmg_related_profiles.txt"

# gene table for the region, proteins of region genes
$PY - <<PYEOF
import re,gzip,sys
gff="$D/data/ncbi_dataset/data/GCA_022453505.1/genomic.gff"
faa="$D/data/ncbi_dataset/data/GCA_022453505.1/protein.faa"
c,s,e="$CONTIG",$REG_START,$REG_END
genes={};cds={}
for l in open(gff):
    if l.startswith('#'):continue
    f=l.rstrip('\n').split('\t')
    if len(f)<9 or f[0]!=c:continue
    a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
    if f[2]=='gene' and int(f[4])>=s and int(f[3])<=e:
        genes[a['ID']]=[int(f[3]),int(f[4]),f[6],a.get('locus_tag',''),'partial' if 'partial' in a else '']
    if f[2]=='CDS' and 'Parent' in a:
        cds.setdefault(a['Parent'],[]).append((int(f[3]),int(f[4]),a.get('Name',''),a.get('product','')))
# mRNA parent map
par={}
for l in open(gff):
    if l.startswith('#'):continue
    f=l.rstrip('\n').split('\t')
    if len(f)>=9 and f[0]==c and f[2]=='mRNA':
        a=dict(x.split('=',1) for x in f[8].split(';') if '=' in x)
        par[a['ID']]=a['Parent']
prot={}
for m,g in par.items():
    if g in genes and m in cds:
        prot[g]=(cds[m][0][2],cds[m][0][3])
seqs={};k=None
for l in open(faa):
    if l[0]=='>':k=l[1:].split()[0];seqs[k]=[]
    else:seqs[k].append(l.strip())
out=open("$W/genes.tsv",'w');pf=open("$W/region_proteins.faa",'w')
out.write("gene\tlocus_tag\tstart\tend\tstrand\tpartial\tprotein\tproduct\tlength_aa\n")
for g,(a,b,st,lt,p) in sorted(genes.items(),key=lambda x:x[1][0]):
    pa,pr=prot.get(g,('',''))
    L=len(''.join(seqs.get(pa,[])))
    out.write(f"{g}\t{lt}\t{a}\t{b}\t{st}\t{p}\t{pa}\t{pr}\t{L}\n")
    if pa in seqs:pf.write(f">{pa} {lt} {a}-{b}{st}\n{''.join(seqs[pa])}\n")
PYEOF
cat "$W/genes.tsv" | cut -c1-200

# HMG-family domains on annotated proteins in the region
hmmsearch --cpu 8 -E 1e-3 --domtblout "$W/region_proteins.hmg.domtbl" "$D/data/hmg_related.hmm" "$W/region_proteins.faa" > /dev/null

# tBLASTn of reference MAT proteins against the region (permissive)
tblastn -query "$W/mat_refs.faa" -subject "$W/region.fa" -evalue 10 -seg no -max_hsps 3 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore sstart send" > "$W/tblastn_refs_vs_region.tsv"

# six-frame ORFs (>= 150 nt, stop-to-stop) over the region and HMM scan
getorf -sequence "$W/region.fa" -minsize 150 -find 1 -outseq "$W/orfs.faa" 2>/dev/null
echo "orfs: $(grep -c '>' "$W/orfs.faa")"
hmmsearch --cpu 8 -E 1e-2 --domtblout "$W/orfs.hmg.domtbl" "$D/data/hmg_related.hmm" "$W/orfs.faa" > /dev/null
echo done
