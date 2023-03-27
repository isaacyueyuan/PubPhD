#!/bin/env Rscript
# author: ph-u
# script: paREF.r
# desc: mark PA ref
# in: Rscript paREF.r
# out: ?
# arg: 0
# date: 20230323

##### env #####
options(warn=2); p9=2
source("src.r");library(ape);library(stringdist);library("xlsx2dfs")
rGene = as.character(read.FASTA(paste0(ptIN,"Pseudomonas_aeruginosa_PAO1_107.ffn")))
d = xlsx2dfs(paste0(ptIN,"Copy of Hypermutator_mutations_ASM-LB_dataframe.xlsx"), rowNames=F, detectDates = F)
cLinical = list.files(ptREF,"fna")

##### PA ref #####
a = unique(c(d[[1]]$LOCUS_TAG[which(d[[1]]$TYPE=="snp" & !is.na(d[[1]]$Strain))],d[[3]]$LOCUS_TAG[which(d[[3]]$TYPE=="snp" & !is.na(d[[3]]$Strain))]));
a = a[!is.na(a)]
a = a[-grep("a",a)] ## PA3574a has no confirmed reference sequence

if(p9 == 1){

## reference & standard of protein function
a = data.frame(locus=a[order(a)], func="", standard="")
for(i in 1:nrow(a)){a1 = strsplit(names(rGene)[grep(a[i,1],names(rGene))],";")
	if(length(grep("product",a1[[1]]))>0){a[i,2] = gsub(",",";",strsplit(a1[[1]][grep("product",a1[[1]])],"=")[[1]][2])}else{a[i,2] = ""}
};rm(i)
write.csv(a,paste0(ptIN,"PAref.csv"),quote=F,row.names=F)

}else if(p9 == 2){

cat("Identify PA gene reference",date(),"\n")
## identify useful gene equivalents by length (foundation assumption of dN/dS)
a0 = as.data.frame(matrix("",nr=length(a),nc=length(cLinical)))
a = data.frame(locus=a[order(a)], len=0); for(i in 1:nrow(a)){a$len[i] = length(rGene[[grep(paste0("=",a$locus[i],";"),names(rGene))]])};rm(i)
for(i in 1:ncol(a0)){
	cat("Processing:",cLinical[i],"(",date(),")\n")
	a1 = as.character(read.FASTA(paste0(ptREF,cLinical[i])))
	colnames(a0)[i] = strsplit(sub("_cds_",";",cLinical[i]),";")[[1]][1]
	a2 = rep(0,length(a1)); for(i0 in 1:length(a2)){a2[i0] = length(a1[[i0]])}
	for(i0 in 1:nrow(a0)){ ## identify (most?) possible gene identity
		a3 = names(a1)[which(a2==a[i0,2])]
		if(length(a3)==0){a4 = ""}else if(length(a3)>1){a4 = as.numeric(adist(names(rGene)[grep(a[i0,1],names(rGene))],a3)[1,])}else{a4 = 1}
		if(length(a3)==0){a0[i0,i] = ""}else{a0[i0,i] = strsplit(sub("[]] [[]pro","@",sub("_tag=","@",a3[which(a4==max(a4))])),"@")[[1]][2]}}
};rm(i,i0,a1,a2,a3,a4)
a0$sample = a$locus
write.csv(a0,paste0(ptIN,"PAref_accessions.csv"), quote=F, row.names=F)

}
