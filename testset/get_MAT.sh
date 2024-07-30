#!/usr/bin/bash -l
#SBATCH  -N 1 -c 1 --mem 2gb --time=1-00:00:00 --out logs/get_MAT.%A.log
module load biopython
python get_MAT.py