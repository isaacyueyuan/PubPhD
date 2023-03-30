#!/bin/bash
# author: ph-u
# script: synBLASTN.sh
# desc: synthesize sequence alignment result
# in: bash synBLASTN.sh
# out: raw/resBLASTN.csv
# arg: 0
# date: 20230329

#SBATCH -J synBSTN
#SBATCH -A WELCH-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --requeue
#SBATCH -p skylake-himem

##### env #####
#i0="`pwd`/`dirname $0`../raw/p_aln"; i1=`echo -e "${i0}" | sed -e "s/[.][.]/./g"`
i1="`pwd`/../raw/p_aln"
[[ -f ${i1}/../resBLASTN.csv ]] && rm ${i1}/../resBLASTN.csv
[[ -f lSt.txt ]] || ls ${i1}/*.txt > lSt.txt

##### read & synthesize data #####
L=0
while read -r i;do
	L=$(( ${L} +1 ))
	i2=`echo -e "${i}" | rev | cut -f 1 -d "/" | rev | cut -f 1,2 -d "." | sed -e "s/-/,/"`
	printf "${i2} ${L}/`wc -l < lSt.txt`\n" #"\r"
	if [[ `grep -e ">" ${i} | wc -l` -gt 0 ]];then
		tail -n +`grep -n ">" ${i} | head -n 1 | cut -f 1 -d ":"` ${i} | head -n 12 > tmp.txt
		if [[ `grep -e "Gaps = 0/" tmp.txt | wc -l` -gt 0 ]];then
			i3=`grep -e "[[]locus_tag=" tmp.txt | sed -e "s/locus_tag=/@/" | sed -e "s/[]] [[]pro/@/" | cut -f 2 -d "@"`
			i4=`grep -e "Sbjct" tmp.txt | head -n 1 | sed -e "s/\t/ /g" | sed -e "s/  / /g" | cut -f 2 -d " "`
		else
			i3="";i4=""
		fi
	else
		i3="";i4=""
	fi
	echo -e "${i2},${i3},${i4}" >> ${i1}/../resBLASTN.csv
done < lSt.txt
rm lSt.txt; rm tmp.txt
printf "\nDone\n"
exit
