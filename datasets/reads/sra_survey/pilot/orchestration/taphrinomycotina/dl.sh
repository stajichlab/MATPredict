#!/bin/bash
O=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/synteny; A=$O/annot; mkdir -p $A; cd $A
get() { # label ftpdir
  b=$(basename $2)
  for ext in genomic.gff protein.faa; do
    [ -s $1.$ext ] && continue
    curl -sS --fail -L -o $1.$ext.gz "$2/${b}_$ext.gz" 2>/dev/null && gunzip -f $1.$ext.gz && echo "ok $1 $ext $(grep -c . $1.$ext)" || echo "missing $1 $ext"
  done
}
G=https://ftp.ncbi.nlm.nih.gov/genomes/all
get Pjir_RU7 $G/GCF/001/477/535/GCF_001477535.1_Pneu_jiro_RU7_V2
get Pcar_B80 $G/GCF/001/477/545/GCF_001477545.1_Pneu_cari_B80_V3
get Scomp $G/GCF/001/661/265/GCF_001661265.1_Saico1
get Spombe $G/GCF/000/002/945/GCF_000002945.2_ASM294v3
get Tdef_A2 $G/GCA/005/281/805/GCA_005281805.1_ASM528180v1
get Tdef_JCM $G/GCA/000/836/155/GCA_000836155.1_Td_assembly_v1.1
get Tdef_PYCC $G/GCA/000/312/925/GCA_000312925.2_ASM31292v2
get Pinouyei $G/GCA/003/717/155/GCA_003717155.1_ASM371715v1
get Plactucae $G/GCA/002/105/105/GCA_002105105.1_Prola1
get Scoloradoensis $G/GCA/051/599/425/GCA_051599425.1_Sai_col_YB2330_v1
get Neolecta $G/GCA/001/929/475/GCA_001929475.1_NeoIirr1.0
get Pgravidus $G/GCA/003/717/255/GCA_003717255.1_ASM371725v1
get Parab $G/GCA/003/568/695/GCA_003568695.1_ASM356869v1
