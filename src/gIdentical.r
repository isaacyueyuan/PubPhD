#!/bin/env Rscript
# author: ph-u
# script: gIdentical.r
# desc: identify identical gene
# in: Rscript gIdentical.r
# out: ?
# arg: 0
# date: 20230326

##### env #####
options(warn=2)
source("src.r");library(ape)
load(paste0(ptOT,"SNP_dNdS.rda"))
aCc = read.csv(paste0(ptIN,"PAref_accessions.csv"), header=T)
cLinical = list.files(ptREF,"fna")

##### identify fasta sources #####
i0 = c();for(i in 1:length(rEs)){if(length(grep("identical", rEs[[i]]))){i0 = c(i0,i)}};rm(i)
i1 = unique(substr(names(rEs)[i0],19,nchar(names(rEs)[i0])))
i2 = unique(substr(names(rEs)[i0],12,17))

##### gene name of identicals #####
i3 = c();for(i in 1:length(i1)){
	s = as.character(read.FASTA(paste0(ptREF,cLinical[grep(i1[i],cLinical)])))
	i3 = c(i3,strsplit(sub("[]] [[]loc","@",sub("gene=","@",names(s)[grep(aCc[which(aCc$sample==i2),which(colnames(aCc)==i1[i])],names(s))])),"@")[[1]][2])
};rm(i,s)
cat("Identical gene =",unique(i3),"\n")
#PAmut[i0,] # check identical entries
