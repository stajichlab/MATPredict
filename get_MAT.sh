#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 1 --mem 2gb --out get_MAT.%A.log
module load biopython
python get_MAT.py
