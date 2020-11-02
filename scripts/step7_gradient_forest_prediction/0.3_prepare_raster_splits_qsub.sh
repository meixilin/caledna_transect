#!/bin/bash
#$ -l highp,highmem,h_rt=24:00:00,h_data=40G,h_vmem=96G
#$ -m bea
#$ -N prep_raster_splits
#$ -wd /u/project/rwayne/meixilin/caledna_transect 
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/0_gradient_forest_prep_raster_splits.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/0_gradient_forest_prep_raster_splits.err.txt
#$ -M meixilin # username on hoffman

# this script takes in the path and name for the gradient forest object
# prepare split raster files for gradient forest predictions 

WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/0.3_prepare_raster_splits.R

source /u/local/Modules/default/init/bash
module load proj/5.2.0
module load gdal/2.3.2
module load geos/3.4.2
module load R/3.5.0

echo "job starts" `date` 
Rscript --vanilla ${RSCRIPT} 
echo "job done" `date`