#!/bin/env Rscript
# author: ph-u
# script: paACC.r
# desc: rearrange aligned accession numbers
# in: Rscript paACC.r
# out: raw/PAref_accessions.csv
# arg: 0
# date: 20230329

#SBATCH -J paACC
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p skylake-himem

##### env #####
options(warn=2)
a = read.csv("../raw/resBLASTN.csv",header=F, stringsAsFactors=F)

a1 = unique(a[,1]); a2 = unique(a[,2])
a0 = as.data.frame(matrix("", nr=length(a1), nc=length(a2)), stringsAsFactors=F)
colnames(a0) = a2

##### rearrange accessions #####
for(i in 1:nrow(a)){ cat(i,"/",nrow(a),"\r")
	a0[which(a1==a[i,1]),which(colnames(a0)==a[i,2])] = a[i,3]
};rm(i)
a0$sample = a1

write.csv(a0, "../raw/PAref_accessions.csv", quote=F, row.names=F)
