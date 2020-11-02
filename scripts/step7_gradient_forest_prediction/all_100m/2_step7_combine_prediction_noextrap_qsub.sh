#!/bin/bash
#$ -N x2_step7_combine_prediction_noextrap
#$ -l highp,highmem,h_rt=36:00:00,h_data=90G,h_vmem=INFINITY
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/2_step7_combine_prediction_noextrap_20201015.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/2_step7_combine_prediction_noextrap_20201015.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub 2_step7_combine_prediction_noextrap_qsub.sh
# @description	combine gradient forest predictions (WITHOUT EXTRAPOLATION)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Oct 15 16:32:55 2020

###########################################################
## import packages 
source /u/local/Modules/default/init/bash
module load R/3.5.0

set -euo pipefail
###########################################################
## def functions 

###########################################################
## def variables 
WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/2_step7_combine_prediction_20201014.R

# set parameters
extrap="FALSE"
gfname="gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
rdatadir="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/split_predictions/" 
outdir="${WORK}/derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/"
###########################################################
## main 
cd ${WORK}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}" 
echo "Using R script:" ${RSCRIPT}
echo "Inputs:" $extrap $rdatadir $outdir

Rscript --vanilla ${RSCRIPT} $extrap $rdatadir $outdir
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
