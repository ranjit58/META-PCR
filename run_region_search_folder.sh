#!/bin/sh
# Read in a top-level folder and all its subfolders
# Make a file for each of the subofolders
IN='/data/scratch/rkumar/work/METAG-PROJ/HMP/GATK-MULTISNP/'
OUT='/data/scratch/bdesilva/PCR_Primer_Discrimination/Output/'

# create the temp files in the OUT dir because assume has write permission already
>${OUT}subdirs.txt
# >${OUT}files.txt

find ${IN} -mindepth 1  -maxdepth 1 -type d | xargs -n 1 basename >> ${OUT}subdirs.txt;
while read sub; do 
#	find ${IN}${sub}/*.vcf | xargs -n 1 basename >> ${OUT}files.txt;
	python region_search.py -i ${IN}${sub}/${sub}.vcf -o ${OUT}${sub}.txt -w 300 -bp -1 -1 100
done < ${OUT}subdirs.txt

rm ${OUT}subdirs.txt
# rm ${OUT}files.txt
