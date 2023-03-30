#!/bin/bash
# author: ph-u
# script: runBLASTN.sh
# desc: run blastn on each locus on each blastdb
# in: bash runBLASTN.sh
# out: ../raw/p_aln/*.txt
# arg: 0
# date: 20230329

##### env #####
#i1="$HOME/Desktop/ncbi_cli"; i2="blastDB"; i3="ncbi-blast-2.13.0+/bin"
i3="`pwd`/`dirname $0`../raw"; i1=`echo -e "${i3}" | sed -e "s/[.][.]/./g"`; i2="blastDB"
i4="../raw/ncbi-genomes-2023-03-23"
#echo -e "$PATH" > tmp.txt
#[[ `grep -e "${i1}" tmp.txt | wc -l` != 1 ]]&&export PATH="${i1}/${i2}:${i1}/${i3}:$PATH"
#rm tmp.txt

##### run blastn #####
#m1=0;m2=$(( `ls ${i1}/${i2}/*.ndb | wc -l` * `ls ../raw/p_paGenes/*.fa | wc -l` ))
for j0 in `ls ${i1}/${i2}/*.ndb`;do
	j2=`echo -e "${j0}" | rev | cut -f 1 -d "/" | rev | cut -f 1,2 -d "."`
	cp ${i1}/${i2}/${j2}.* .
	sbatch runBLASTN_c.sh ${j2} &
#	for j1 in `ls ../raw/p_paGenes/*.fa`;do
#		m1=$(( ${m1}+1 ))
#		j3=`echo -e "${j1}" | rev | cut -f 1 -d "/" | rev | cut -f 1 -d "."`
#		m3=$(( ${m1} / ${m2} *100 ))
#		printf "${j3}-${j2} (${m1}/${m2} - Percentage: ${m3})\r"
#		sbatch runBLASTN_c.sh ${j1} ${j2} ${j3} &
		#blastn -query ${j1} -db ${j2} -out ../raw/p_aln/${j3}-${j2}.txt
#	done; rm ${j2}.*
done;exit
