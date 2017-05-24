#!/bin/sh
#SBATCH --partition=express
#SBATCH --job-name=VCFConfirmation
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --time=0-2:00:00
#SBATCH --mem=16000

#IN='/data/scratch/rkumar/work/METAG-PROJ/HMP/GATK-MULTISNP/'
IN='/data/scratch/bdesilva/PCR_Primer_Discrimination/RewrittenVCF/'
OUT='/data/scratch/bdesilva/PCR_Primer_Discrimination/RewrittenVCF2/'
SUBDIRS='subdirs.txt'
PARALLEL='/data/scratch/bdesilva/software/parallel-20160922/src/parallel' 
RW='/home/bdesilva/PCR_Primer_Discrimination/tools/rewrite_increment_by_one_vcf.py' # rewrite

# create out dir
mkdir -p ${OUT}

# create the temp files in the OUT dir because assume has write permission already
>${OUT}${SUBDIRS}
find ${IN} -mindepth 1  -maxdepth 1 -type d | xargs -n 1 basename >> ${OUT}${SUBDIRS};

while read directory; do
	mkdir ${OUT}${directory}
done < ${OUT}${SUBDIRS}

${PARALLEL} -a ${OUT}${SUBDIRS} python ${RW} -i ${IN}{}/{}.vcf -o ${OUT}{}/{}.vcf

rm ${OUT}${SUBDIRS}
