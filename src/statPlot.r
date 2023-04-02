#!/bin/env Rscript
# author: ph-u
# script: statPlot.r
# desc: dN/dS statistics calculation + plotting
# in: Rscript statPlot.r
# out: ?
# arg: 0
# date: 20230325

##### env #####
options(warn=2);source("src.r")
load(paste0(ptOT,"SNP_dNdS.rda"))
#PAmut = read.csv(paste0(ptOT,"geneDNDS.csv"), header=T) # [visually identical]
cOl = c("#00000000","#00000033","#000000bb","#000000ff")

##### Analysis data extraction #####
a0 = PAmut[is.finite(PAmut$dNdS),]
a1 = list(log10(a0$dNdS[which(a0$Medium=="ASM")]),log10(a0$dNdS[which(a0$Medium=="LB")]))
a2 = list(summary(a1[[1]]),summary(a1[[2]]))
a3 = list(density(a1[[1]]),density(a1[[2]]))

##### Median difference test #####
sT = wilcox.test(a1[[1]],a1[[2]])

##### Plot #####
ySt = max(log10(PAmut$dNdS[which(is.finite(PAmut$dNdS))]))*.7 # stat result
tiff(paste0(ptOT,"ASMvsLB.tif"), width=700,height=1400, compression="lzw")
par(mar=c(4,3.7,.5,.5)+.1, mfrow=c(2,1), xpd=F, cex=3)
plot(as.factor(paste0(a0$Medium[which(is.finite(a0$dNdS) & a0$dNdS!=0)],"A")),log10(a0$dNdS[which(is.finite(a0$dNdS) & a0$dNdS!=0)]), pch=3, xlab="Medium", ylab="log10( dN/dS )", col=cOl[1], cex=.8)
#text(c(.7,2.3), c(a2[[1]][which(names(a2[[1]])=="Mean")],a2[[2]][which(names(a2[[2]])=="Mean")]), labels=paste0("Mean = ",round(c(a2[[1]][which(names(a2[[1]])=="Mean")],a2[[2]][which(names(a2[[2]])=="Mean")]),2)), cex=.63, srt=0)
text(c(.8,2.2), c(a2[[1]][which(names(a2[[1]])=="Median")],a2[[2]][which(names(a2[[2]])=="Median")])*.8, labels=paste0("Median = ",round(c(a2[[1]][which(names(a2[[1]])=="Median")],a2[[2]][which(names(a2[[2]])=="Median")]),2)), cex=.63, srt=0)
#text(rep(c(.7,2.3), each=length(a2[[1]])), c(a2[[1]],a2[[2]]), labels=paste(c(names(a2[[1]]),names(a2[[2]])),round(c(a2[[1]],a2[[2]]),2), sep="="), cex=.3, srt=45)
#segments(x=1,y=a2[[1]][which(names(a2[[1]])=="Mean")], x1=2,y1=a2[[2]][which(names(a2[[1]])=="Mean")], lwd=3)
#segments(x=1,y=a2[[1]][which(names(a2[[1]])=="Median")], x1=2,y1=a2[[2]][which(names(a2[[1]])=="Median")], lwd=3)
segments(x=1,y=ySt,x1=2,y1=ySt);text(1.5,ySt/7*9,labels=paste0("Wilcox\nW = ",sT$statistic,", p = ",round(sT$p.value,4)), cex=.55)

# https://stackoverflow.com/questions/38309547/r-density-plot-with-colors-by-group
plot(a3[[2]], col=cOl[2], xlab="log10( dN/dS )", main="", lwd=3)
lines(a3[[1]], col=cOl[3], lty=2, lwd=3)
legend("topleft", legend=c("ASMA","LBA"), lty=c(2,1), lwd=3, col=cOl[c(3,2)], bty="n")
abline(v=0, lwd=2)
invisible(dev.off())

##### Get density peaks & troughs #####
z = c(-.05,1.5,"#33ff00aa","#0033ff77")
z0 = c(-.51,-.815,-.66,-.93,-.995, -.46,-.57,-.763,-.955,-.86)
pdf(paste0(ptOT,"ASMvsLB_peaks.pdf"));par(mar=c(4,4,0,0)+.1)
plot(a3[[2]]$x,a3[[2]]$y, pch=20, cex=.1,col="#ff0000ff", xlab="log10( dN/dS )", ylab="Density")
points(a3[[1]]$x,a3[[1]]$y, pch=20, cex=.1)
abline(v=seq(-1.2,-.3,.1),h=seq(.3,1.5,.1), col="#0000ff44")
i1 = c();for(i in 1:length(z0)){
	rect(z0[i]-.01,as.numeric(z[1]),z0[i]+.01,as.numeric(z[2])-ifelse(i>5,.3,0),col=z[ifelse(i>5,4,3)], border=NA)
	i1 = c(i1,i0 <- unique(a0$LOCUS_TAG[which(log10(a0$dNdS)>(z0[i]-.01) & log10(a0$dNdS)<(z0[i]+.01) & a0$Medium==ifelse(i>5,"ASM","LB"))]))
	cat("Medium:",ifelse(i>5,"ASM","LB"),"; dNdS range:",z0[i]-.01,"-",z0[i]+.01,";Genes:",paste0(i0[order(i0)], collapse=","),"\n")
};rm(i,i0)
invisible(dev.off())

pdf(paste0(ptOT,"ASMvsLB_freq.pdf"), width=21);par(mar=c(4,4,0,0)+.1)
plot(table(i1), xlab="", ylab="Frequency",las=2)
invisible(dev.off())
