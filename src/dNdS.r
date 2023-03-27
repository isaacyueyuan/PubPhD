#!/bin/env Rscript
# author: ph-u
# script: dNdS.r
# desc: calculate dN/dS ratio for gene
# in: Rscript dNdS.r
# out: ?
# arg: 0
# date: 20230317

##### env #####
options(warn=2)
source("src.r");library(ape);library("xlsx2dfs")
rGene = as.character(read.FASTA(paste0(ptIN,"Pseudomonas_aeruginosa_PAO1_107.ffn")))
d = xlsx2dfs(paste0(ptIN,"Copy of Hypermutator_mutations_ASM-LB_dataframe.xlsx"), rowNames=F, detectDates = F)
cLinical = list.files(ptREF,"fna")
aCc = read.csv(paste0(ptIN,"PAref_accessions.csv"), header=T)

##### Sequence data sources #####
# https://pseudomonas.com/strain/download
# https://www.ncbi.nlm.nih.gov/assembly?LinkName=bioproject_assembly_all&from_uid=325248
# download assemblies -> genomic FASTA from GenBank

## overall mutation data
cat("Preparing blank record sheet",date(),": ")
PAmut0 = unique(rbind(d[[1]],d[[3]])[,c("Strain","Medium","REGION","LOCUS_TAG")])
for(i in 1:ncol(PAmut0)){PAmut0 = PAmut0[!is.na(PAmut0[,i]),]};rm(i)
PAmut = PAmut0; for(i in 2:length(cLinical)){PAmut = rbind(PAmut,PAmut0)};rm(i)
PAmut$clinicalREF = rep(colnames(aCc)[-ncol(aCc)],each=nrow(PAmut0))
j0 = colnames(PAmut)[length(colnames(PAmut))]
cat(nrow(PAmut),"rows\n")

rEs = vector(mode="list", length=nrow(PAmut))
names(rEs) = apply(PAmut[,-3],1,paste0,collapse=".")
PAmut[,c("dNdS","pN","pS","Nd","Sd","N","S")] = NA

##### data record #####
cat("dNdS calculation",date(),"\n")
for(i in 1:nrow(PAmut)){ if(i%%1000==0){cat(i,"(",date(),")",round(i/nrow(PAmut),4)*100,"%\n")}
	j=0; if(i==1){j=1}else if(PAmut$clinicalREF[i]!=PAmut$clinicalREF[i-1]){j=1}
	if(j==1){kIn = as.character(read.FASTA(paste0(ptREF,cLinical[grep(PAmut$clinicalREF[i],cLinical)])))}
	sAm = rGene[[grep(paste0("=",PAmut$LOCUS_TAG[i],";"),names(rGene))]]
	k0 = aCc[grep(PAmut$LOCUS_TAG[i],aCc$sample),grep(PAmut$clinicalREF[i],colnames(aCc))]
	if(length(grep("_",k0))>0){
		refG = kIn[[grep(k0,names(kIn))]]
		refM = d[[ifelse(PAmut$Medium[i]=="ASM",1,3)]]
		refM = refM[which(refM$Strain==PAmut$Strain[i] & refM$LOCUS_TAG==PAmut$LOCUS_TAG[i] & !is.na(refM[,1])),]
		if(length(unique(refM$TYPE))==1){ if(unique(refM$TYPE)=="snp"){
			for(i0 in 1:nrow(refM)){sAm[as.numeric(strsplit(refM$NT_POS[i0],"/")[[1]][1])] = nFlip(refM$ALT[i0],refM$STRAND[i0])}
			if(length(sAm)>as.numeric(strsplit(refM$NT_POS[i0],"/")[[1]][1])){
				bIn = dNdS.rt(paste0(sAm, collapse=""), paste0(refG, collapse=""))
				}else{bIn=list(rep(NA,ncol(PAmut)-grep(j0,colnames(PAmut))),"Seq not match record")}
			}else{bIn=list(rep(NA,ncol(PAmut)-grep(j0,colnames(PAmut))),refM$TYPE)}
			}else{bIn=list(rep(NA,ncol(PAmut)-grep(j0,colnames(PAmut))),refM$TYPE)}
	}else{bIn=list(rep(NA,ncol(PAmut)-grep(j0,colnames(PAmut))),"Seq without identifiable clinical equivalent")}
	PAmut[i,-c(1:grep(j0,colnames(PAmut)))] = bIn[[1]]
	rEs[[i]] = bIn[[2]]
};rm(i,j,sAm,refG,refM,bIn,kIn);cat("Saving to RData",date(),"\n")
write.csv(PAmut,paste0(ptOT,"geneDNDS.csv"),row.names=F, quote=F)
capture.output(rEs,file=paste0(ptOT,"geneDNDS.txt"))
save(rEs, PAmut, file=paste0(ptOT,"SNP_dNdS.rda"))
