#!/bin/bash
#$ -N util_compare_previous
#$ -l highp,h_rt=36:00:00,h_data=30G
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/util_compare_previous_20201020.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/util_compare_previous_20201020.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub util_compare_previous_qsub.sh
# @description	compare with the previous version 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Oct 20 16:23:21 2020

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
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/util_compare_previous.R
LOG="${WORK}/derive_data/step7_gradient_forest_prediction/logs/util_compare_previous_${TODAY}.log"

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