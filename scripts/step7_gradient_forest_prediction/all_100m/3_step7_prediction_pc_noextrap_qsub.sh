#!/bin/bash
#$ -N x3_step7_prediction_pc_noextrap
#$ -l highp,highmem,h_rt=23:00:00,h_data=90G,h_vmem=INFINITY
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/3_step7_prediction_pc_noextrap_20201018.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/3_step7_prediction_pc_noextrap_20201018.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub 3_step7_prediction_pc_noextrap_qsub.sh
# @description	perform principal component analyses on gradient forest predictions (WITHOUT EXTRAPOLATION)
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
WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/3_step7_prediction_pc_20201014.R

# set parameters
extrap="FALSE"
gfname="gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
outdir="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/"

Trns_site_filepath="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/predict_gf_Trns_site_noextrap.RData"
Trns_CA_filepath="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/predict_gf_Trns_CAgrid_noextrap.RData"

###########################################################
## main 
cd ${WORK}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" 
echo "Using R script:" ${RSCRIPT}
echo "Inputs:" $extrap $Trns_site_filepath $Trns_CA_filepath $outdir

Rscript --vanilla ${RSCRIPT} $extrap $Trns_site_filepath $Trns_CA_filepath $outdir
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"

