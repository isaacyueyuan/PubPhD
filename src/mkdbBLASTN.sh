#!/bin/bash
# author: ph-u
# script: mkdbBLASTN.sh
# desc: section sequence for blastn local query
# in: bash mkdbBLASTN.sh
# out: ~/Desktop/ncbi_cli/blastDB/*
# arg: 0
# date: 20230328

##### ref #####
# https://ncbi.github.io/magicblast/cook/blastdb.html

##### env #####
i0=`pwd`
i1="$HOME/Desktop/ncbi_cli"; i2="blastDB"; i3="ncbi-blast-2.13.0+/bin"
i4="../raw/ncbi-genomes-2023-03-23"
echo -e "$PATH" > tmp.txt
[[ `grep -e "${i1}/${i2}" tmp.txt | wc -l` != 1 ]]&&export PATH="${i1}/${i2}:$PATH"
[[ `grep -e "${i1}/${i3}" tmp.txt | wc -l` != 1 ]]&&export PATH="${i1}/${i3}:$PATH"
rm tmp.txt

##### make blastdb #####
for i in `ls ${i4}/*.fna`;do
	i5=`echo -e "${i}" | rev | cut -f 1 -d "/" | rev | cut -f 1,2,3 -d "_"`
	makeblastdb -in ${i} -dbtype nucl -parse_seqids -out ${i1}/${i2}/${i5} -title ${i5} 1> p_mkdbBLASTNmsg.log
done
exit
