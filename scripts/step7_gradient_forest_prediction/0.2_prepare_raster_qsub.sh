#!/bin/bash
#$ -l highmem,highp,h_rt=03:00:00,h_data=96G
#$ -m bea
#$ -N gradient_forest_prepare_gf
#$ -o /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step6_gf_prediction/logs/gradient_forest_prepare_gf_06062019.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step6_gf_prediction/logs/gradient_forest_prepare_gf_06062019.err.txt
#$ -cwd #current working directory
#$ -M meixilin # username on hoffman

# this script takes in the path and name for the gradient forest object
# for setting up LOG


WORK="/u/project/rwayne/meixilin/abiotic_transect"
# WORK="/Users/linmeixi/UCLA/Lab/abiotic_transect/"
RSCRIPT=${WORK}/r_codes/step6_gradient_forest_prediction/0.2_prepare_raster.R

source /u/local/Modules/default/init/bash
module load proj/5.2.0
module load gdal/2.3.2
module load geos/3.4.2
module load R/3.5.0

echo "job starts" `date` 
Rscript --vanilla ${RSCRIPT} 
echo "job done" `date`