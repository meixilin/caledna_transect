#!/bin/bash
#$ -N x4.3_step7_store_ref_scale_grid
#$ -l highp,highmem,exclusive,h_rt=36:00:00,h_data=240G,h_vmem=INFINITY
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/4.3_step7_store_ref_scale_grid_20201019.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/4.3_step7_store_ref_scale_grid_20201019.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub 4.3_step7_store_ref_scale_grid_qsub.sh
# @description	store scaling of the CAgrid 
# @previous		step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_qsub.sh
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Oct 19 22:15:27 2020

###########################################################
## import packages 
source /u/local/Modules/default/init/bash
module load proj/5.2.0
module load gdal/2.3.2
module load R/3.5.0
module load gcc/6.3.0

set -euo pipefail

###########################################################
## def variables 
TODAY=$(date "+%Y%m%d")
WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/4.3_step7_store_ref_scale_grid_20201014.R
LOG="${WORK}/derive_data/step7_gradient_forest_prediction/logs/4.3_step7_store_ref_scale_grid_${TODAY}.log"

###########################################################
## main 
cd ${WORK}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" 
echo "Using R script:" ${RSCRIPT}
echo "Inputs: None" 
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" > ${LOG}
echo "Using R script:" ${RSCRIPT} >> ${LOG}
echo "Inputs: None" >> ${LOG}

Rscript --vanilla ${RSCRIPT} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
