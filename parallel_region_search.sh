#!/bin/sh
#SBATCH --partition=express
#SBATCH --job-name=ParallelRegionSearch
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mem=1000

#IN='/data/scratch/rkumar/work/METAG-PROJ/HMP/GATK-MULTISNP/'
WINDOW_SIZE='200'
IN='/data/scratch/bdesilva/PCR_Primer_Discrimination/RewrittenVCF/'
OUT='/data/scratch/bdesilva/PCR_Primer_Discrimination/Analysis_Multiple_Step/'${WINDOW_SIZE}'/'
SUBDIRS='subdirs.txt'
PARALLEL='/data/scratch/bdesilva/software/parallel-20160922/src/parallel' 
RS='/home/bdesilva/PCR_Primer_Discrimination/tools/region_search.py' # region search

# create out dir
mkdir -p ${OUT}

# create the temp files in the OUT dir because assume has write permission already
>${OUT}${SUBDIRS}
find ${IN} -mindepth 1  -maxdepth 1 -type d | xargs -n 1 basename >> ${OUT}${SUBDIRS};

${PARALLEL} -a ${OUT}${SUBDIRS} python ${RS} -i ${IN}{}/{}.vcf -w ${WINDOW_SIZE} -bp -1 -1 25 -u . -l 50 -f 25 -p 0 -o ${OUT}{}.txt

rm ${OUT}${SUBDIRS}
