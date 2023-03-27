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
a1 = list(abs(log10(a0$dNdS[which(a0$Medium=="ASM")]+1)),abs(log10(a0$dNdS[which(a0$Medium=="LB")]+1)))
a2 = list(summary(a1[[1]]),summary(a1[[2]]))

##### Median difference test #####
sT = wilcox.test(a1[[1]],a1[[2]])

##### Plot #####
tiff(paste0(ptOT,"ASMvsLB.tif"), width=700,height=1400, compression="lzw")
par(mar=c(4,3.7,0,0)+.1, mfrow=c(2,1), xpd=F, cex=3)
plot(as.factor(a0$Medium),abs(log10(a0$dNdS+1)), pch=3, xlab="Medium", ylab="| log10( dN/dS + 1 ) |", col=cOl[1], cex=1)
text(c(.7,2.3), c(a2[[1]][which(names(a2[[1]])=="Mean")],a2[[2]][which(names(a2[[2]])=="Mean")]), labels=paste0("Mean = ",round(c(a2[[1]][which(names(a2[[1]])=="Mean")],a2[[2]][which(names(a2[[2]])=="Mean")]),2)), cex=.63, srt=0)
#text(rep(c(.7,2.3), each=length(a2[[1]])), c(a2[[1]],a2[[2]]), labels=paste(c(names(a2[[1]]),names(a2[[2]])),round(c(a2[[1]],a2[[2]]),2), sep="="), cex=.3, srt=45)
segments(x=1,y=a2[[1]][which(names(a2[[1]])=="Mean")], x1=2,y1=a2[[2]][which(names(a2[[1]])=="Mean")], lwd=3)
segments(x=1,y=1.3,x1=2,y1=1.3);text(1.5,1.4,labels=paste0("Wilcox\nW = ",sT$statistic,", p = ",round(sT$p.value,4)), cex=.55)

# https://stackoverflow.com/questions/38309547/r-density-plot-with-colors-by-group
plot(density(log10(a0$dNdS[which(a0$Medium!="ASM")]+1)), col=cOl[2], xlab="log10( dN/dS + 1 )", main="", lwd=3)
lines(density(log10(a0$dNdS[which(a0$Medium=="ASM")]+1)), col=cOl[3], lty=2, lwd=3)
legend("topright", legend=c("ASM","LB"), lty=c(2,1), lwd=3, col=cOl[c(3,2)], bty="n")
abline(v=0, lwd=2)
invisible(dev.off())
