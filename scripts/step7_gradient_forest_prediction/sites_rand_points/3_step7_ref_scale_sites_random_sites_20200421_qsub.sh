#!/bin/bash
#$ -l highp,highmem,h_rt=05:00:00,h_data=40G,h_vmem=200G
#$ -m bea
#$ -N scale_sites_random_sites
#$ -wd /u/project/rwayne/meixilin/caledna_transect 
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/scale_sites_random_sites.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/scale_sites_random_sites.err.txt
#$ -M meixilin # username on hoffman

# this script is for reference transformation of the sites/random sites in california 

WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/sites_rand_points/3_step7_ref_scale_sites_random_sites_20200421.R

cd $WORK
source /u/local/Modules/default/init/bash

module load proj/5.2.0
module load gdal/2.3.2
module load R/3.5.0
module load gcc/6.3.0

echo "job starts" `date` 
Rscript --vanilla ${RSCRIPT} 
echo "job done" `date`