#!/bin/sh
#SBATCH --partition=medium
#SBATCH --job-name=exp2igv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02-2:00:00
#SBATCH --mem=10000

IN='/data/scratch/bdesilva/PCR_Primer_Discrimination/'
FOLDER='Analysis_Multiple_Step/'
WINDOW_SIZE='550'

python export_to_igv.py -i ${IN}${FOLDER}${WINDOW_SIZE} -o ${IN}Export_IGV/${FOLDER}${WINDOW_SIZE} -e /200/ /500/test/ -d ${IN}RewrittenVCF/
