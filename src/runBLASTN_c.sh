#!/bin/bash
# author: ph-u
# script: runBLASTN_c.sh
# desc: sequence alignment via blastn
# in: bash runBLASTN_c.sh [dbNam]
# out: raw/p_aln/[query]-[dbNam].txt
# arg: 1
# date: 20230329

#SBATCH -J blastn
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p skylake-himem

##### blastn for SLURM #####
j2=$1
for j1 in `ls ../raw/p_paGenes/*.fa`;do
	j3=`echo -e "${j1}" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "."`
	blastn -query ${j1} -db ${j2} -out ../raw/p_aln/${j3}-${j2}.txt
done; rm ${j2}.*
exit
