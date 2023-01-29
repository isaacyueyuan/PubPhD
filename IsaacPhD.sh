#!/bin/bash
# author: ph-u
# script: IsaacPhD.sh
# desc: calling IsaacPhD.r for pairwise t.test with False Discovery Rate (FDR) p-value adjustment
# in: bash IsaacPhD.sh [directory/containing/data/csv] [threshold (decimal number 0~1, general .05)] [FDR (decimal number 0~1)]
# out: none [see respective script]
# arg: 3
# date: 20211105, 20211106, 20211108

[ -z $1 ] && grep "desc:\|in:" $0 | grep -v "grep" | cut -f 2 -d ":" && exit 0
pAth=`dirname $0`
cd $1
[ -f result.txt ] && rm result.txt
touch result.txt
for i in `ls *_*h.csv`;do
	nAm=`echo -e ${i} | cut -f 1 -d "."`
	rEf=`echo ${nAm} | cut -f 2 -d "_"`
	grep -e "${rEf}\|^${rEf%?}\+" ${i} > ${nAm}-d.csv #cut -f 2- -d "," ${i}
	echo -e "${i}" >> result.txt
	Rscript ${pAth}/IsaacPhD.r ${nAm}-d.csv ${rEf} $2 $3 >> result.txt
	rm ${nAm}-d.csv
done
echo -e "Small data with comparable amount of multiple comparisons using traditional p-correction can lead to false negative conclusions, so two False Discovery Rate (FDR) p-value correction methods are adopted;\nConclusions are drawn from the ensemble result of\n\n- 'Benjamini and Hochberg (1995)'[https://doi.org/10.1111/j.2517-6161.1995.tb02031.x]\n- 'Storey and Tibshirani (2003)'[https://doi.org/10.1214/aos/1074290335]\n\nRelated description article [https://doi.org/10.1016/b978-0-12-802101-9.00017-x]\nA useful simple lecture [https://evolution.gs.washington.edu/gs560/2011/lecture9.pdf]\nA useful online article [https://www.statisticshowto.com/benjamini-hochberg-procedure/]" >> result.txt
