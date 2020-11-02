#! /bin/bash
#$ -N x1_step7_make_prediction_noextrap
#$ -l h_rt=05:00:00,h_data=15G,h_vmem=20G
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/1_step7_make_prediction_noextrap_20201013.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/1_step7_make_prediction_noextrap_20201013.err.txt
#$ -m abe

# @version 		v1
# @usage		qsub -t 1-51 1_step7_make_prediction_noextrap_qsub.sh
# @description	generate gradient forest predictions (WITHOUT EXTRAPOLATION)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Oct 14 13:22:38 2020

###########################################################
## import packages 
sleep $((RANDOM % 60)) 
source /u/local/Modules/default/init/bash
module load R/3.5.0

set -euo pipefail
###########################################################
## def functions 

###########################################################
## def variables 
# get input index (II)
IDX=$((SGE_TASK_ID-1))
II=$(printf %02d ${IDX})

WORK="/u/project/rwayne/meixilin/caledna_transect"
RSCRIPT=${WORK}/scripts/step7_gradient_forest_prediction/all_100m/1_step7_make_prediction_all_20201014.R
LOG=${WORK}/derive_data/step7_gradient_forest_prediction/logs/1_step7_make_prediction_noextrap_20201013_split${II}.log
# set parameters
extrap="FALSE"
gfdir="./derive_data/step6_gradient_forest/"
gfname="gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
outdir="./derive_data/step7_gradient_forest_prediction/$gfname/extrap_$extrap/split_predictions/" 

###########################################################
## main 
cd ${WORK}
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}" 
echo "Using R script:" ${RSCRIPT}
echo "Inputs:" $extrap $gfdir $gfname $outdir $II

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}" > ${LOG}
echo "Using R script:" ${RSCRIPT} >> ${LOG}
echo "Inputs:" $extrap $gfdir $gfname $outdir $II >> ${LOG}

Rscript --vanilla ${RSCRIPT} $extrap $gfdir $gfname $outdir $II &>> ${LOG}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

 
