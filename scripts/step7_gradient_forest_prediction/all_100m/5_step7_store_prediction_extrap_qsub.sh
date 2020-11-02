#!/bin/bash
#$ -N x5_step7_store_prediction_extrap
#$ -l highp,highmem,exclusive,h_rt=36:00:00,h_data=240G,h_vmem=INFINITY
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/5_step7_store_prediction_extrap_20201019.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/5_step7_store_prediction_extrap_20201019.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub 5_step7_store_prediction_extrap_qsub.sh
# @description	store as raster format on gradient forest predictions (WITH EXTRAPOLATION)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Oct 18 20:19:11 2020

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
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_20201014.R
LOG="${WORK}/derive_data/step7_gradient_forest_prediction/logs/5_step7_store_prediction_extrap_${TODAY}.log"

# set parameters
extrap="TRUE"
gfname="gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
outdir="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/"

Trns_CA_filepath="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/predict_gf_Trns_CAgrid_extrap.RData"

###########################################################
## main 
cd ${WORK}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" 
echo "Using R script:" ${RSCRIPT}
echo "Inputs:" $extrap $Trns_CA_filepath $outdir
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" > ${LOG}
echo "Using R script:" ${RSCRIPT} >> ${LOG}
echo "Inputs:" $extrap $Trns_CA_filepath $outdir >> ${LOG}

Rscript --vanilla ${RSCRIPT} $extrap $Trns_CA_filepath $outdir &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}

