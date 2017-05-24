#!/bin/sh
#SBATCH --partition=express
#SBATCH --job-name=PCR_D_IND
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem 1000

#IN='/data/scratch/rkumar/work/METAG-PROJ/HMP/GATK-MULTISNP/Bacteroides_vulgatus/'
#IN='/data/scratch/bdesilva/PCR_Primer_Discrimination/RewrittenVCF/Bacteroides_vulgatus/'
IN='/home/bdesilva/PCR_Primer_Discrimination/tools/'
#finame='Bacteroides_vulgatus.vcf'
finame='B3.vcf'
#IN='/data/scratch/bdesilva/Bacteroides_vulgatus/'
#finame='Bacteroides_vulgatus.vcf'
FIN=${IN}${finame}
OUT='/data/scratch/bdesilva/PCR_Primer_Discrimination/Rewritten_VCF_Output/'
mkdir -p $OUT
foname='Bacteroides_vulgatus_new_output3.txt'
FOUT=${OUT}${foname}
delim='\t'
python -m cProfile -s time region_search.py -i ${FIN} -o ${FOUT} -w 3 -bp -1 -1 1 -u . -l 3 -f 1 -p 0
