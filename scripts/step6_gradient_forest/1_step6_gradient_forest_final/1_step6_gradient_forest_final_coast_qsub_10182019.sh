#!/bin/bash
#$ -l highmem,highp,h_rt=03:00:00,h_data=96G
#$ -m bea
#$ -N gradient_forest_coast
#$ -o /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step4_gradient_forest/2_final/logs/coast_gf_10182019.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step4_gradient_forest/2_final/logs/coast_gf_10182019.err.txt
#$ -cwd #current working directory
#$ -M meixilin # username on hoffman

# this script takes in the path and name for the gradient forest object
# for setting up LOG


WORK="/u/project/rwayne/meixilin/abiotic_transect/"
# WORK="/Users/linmeixi/UCLA/Lab/abiotic_transect/"
RSCRIPT=${WORK}/r_codes/step4_gradient_forest/2_step4_gf_final/2_step4_gradient_forest_final_coast_10182019.R

source /u/local/Modules/default/init/bash
module load R/3.5.0

echo "job starts" `date` 
Rscript --vanilla ${RSCRIPT} 
echo "job done" `date`