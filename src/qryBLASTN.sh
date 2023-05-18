#!/bin/bash
# author: ph-u
# script: qryBLASTN.sh
# desc: section sequence for blastn local query
# in: bash qryBLASTN.sh
# out: ../raw/p_paGenes/*
# arg: 0
# date: 20230328

##### ref #####
# [blast exec] https://www.ncbi.nlm.nih.gov/books/NBK52640/
# [blast manual] https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
# [blast db] https://ftp.ncbi.nlm.nih.gov/blast/db/
# [blastdb build] https://ftp.ncbi.nlm.nih.gov/blast/documents/blastdb.html

##### env #####
i0="../raw/p_paGenes"
mkdir -p ${i0}
i1=`ls ../raw/Pse*.ffn`
[[ -f ${i0}/accession.txt ]]&& rm ${i0}/accession.txt

awk '{FS=","}{print $1}' < ../raw/p_PAref.csv | grep -v "a$" | tail -n +2 > ${i0}/paLocus.txt
grep -n ">" ${i1} > ${i0}/gLine.txt
#echo -e "`wc -l < ${i1}`: ${i1}" >> ${i0}/gLine.txt

##### blastn #####
#~/Desktop/ncbi_cli/ncbi-blast-2.13.0+/bin/blastn
while read -r L;do
#	grep -e "=${L};" ${i0}/gLine.txt | sed -e "s/replicon_accession=/@/" | sed -e "s/;product/@/" | cut -f 2 -d "@" >> ${i0}/accession.txt
	i2=`grep -e "=${L};" ${i0}/gLine.txt | cut -f 1 -d ":"`
	i3=`grep -n "=${L};" ${i0}/gLine.txt | cut -f 1 -d ":"`
	i4=$(( `tail -n +$(( ${i3}+1 )) ${i0}/gLine.txt | head -n 1 | cut -f 1 -d ":"` -1 ))
#	echo -e "${L},${i2},${i3},${i4}"
	head -n ${i4} ${i1} | tail -n +${i2} > ${i0}/${L}.fa
done < ${i0}/paLocus.txt
exit
