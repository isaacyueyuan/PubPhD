#!/bin/env Rscript
# author: ph-u
# script: IsaacPhD.r
# desc: log treatment data, pairwise t.test with False Discovery Rate (FDR) p-value adjustment
# in: Rscript IsaacPhD.r [data.csv] [refName] [threshold] [FDR]
# out: terminal output
# arg: 4
# date: 20211105, 20211106, 20211108, 20211110

argv = (commandArgs(T))
a = log(as.data.frame(t(read.csv(argv[1],header=F,row.names=1)))+1)
k = ncol(a)-1
cat(paste0("number of pairwise tests: ",k,"\n\n"))
pVal = rep(NA,ncol(a))
for(i in 1:ncol(a)){
	if(colnames(a)[i]!=argv[2]){
		cat(colnames(a)[i])
		s = t.test(a[,i],a[,argv[2]], paired=T)
		print(s)
		pVal[i] = s$p.value
#		j = s$p.value*k
#		cat(paste0("Bonferroni adjusted pairwise p-value (against threshold ",argv[3],"): ",round(j,4),"\nSignificant? "))
#		if(j < as.numeric(argv[3])){
#			cat("YES")
#		}else if (j < as.numeric(argv[3])*1.5){
#			cat("ALMOST")
#		}else{cat("NO")}
#		cat("\n\n%%%%%%%%%%%%\n\n")
	}
}

pCor = as.data.frame(matrix(NA,nr=length(pVal),nc=6))
colnames(pCor) = c("Paired_Treatment","raw_pVal","BH1995_criterion","BH_result","ST2003_result","Bonferroni")
pCor[,1] = colnames(a)[order(pVal)]
pCor[,2] = pVal[order(pVal)]
pCor = pCor[-nrow(pCor),]

##### Benjamini and Hochberg (1995) False Discovery Rate p-value correction #####
FDR = as.numeric(argv[4])
cat(paste0("##### Benjamini and Hochberg (1995)\nFalse Discovery Rate (FDR, arbitrary number ",FDR," by user)\n",FDR," = user believed that around ",FDR*100,"% of the multiple comparisons are false discovery (i.e. p-value significance by random chance)\ncriterion calculation formula: 1 / [#tests] * FDR * [#rank of accending-ordered p-values from the multiple comparisons]\n\n"))
pCor[,"BH1995_criterion"] = 1/k*FDR*as.numeric(row.names(pCor))
pCor[,"BH_result"] = ifelse(pCor[,"raw_pVal"]>pCor[,"BH1995_criterion"],"NO","Yes")

##### Storey and Tibshirani (2003) False Discovery Rate p-value correction #####
tHres = .5
dIs = sum(pCor[,2]>tHres)*2
cat(paste0("##### Storey and Tibshirani (2003)\nIf data is random, p-value should follow an uniform distribution around mean ",tHres,"\nany derivations can be categorized as non-random (hence reject respective H0)\ndiscard 2*[number of tests with raw p-val >",tHres,"] -- #tests to be discarded = ",dIs,"\n\n"))
pCor[,"ST2003_result"] = c(rep("Yes",nrow(pCor)-dIs),rep("NO",dIs))

##### Bonferroni p-value adjustment #####
lIm = as.numeric(argv[3])/k
pCor[,"Bonferroni"] = ifelse(pCor[,"raw_pVal"]>lIm,"NO","Yes")
cat(paste0("##### Bonferroni formula for significance: [raw p-value] < [threshold] / [# of repeated statistical tests]\n\n"))

print(pCor)

##### Justification #####
cat(paste0("The scale of the number of data points: ",nrow(a)," is differ from the number of pairwise test: ",k," by ratio of [(#data)/(#tests)] ",nrow(a)/k,"\n%%%%%%%%%%%%%%%%%%%%%\n\n"))
