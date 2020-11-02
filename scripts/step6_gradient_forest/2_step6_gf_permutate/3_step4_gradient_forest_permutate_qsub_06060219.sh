#!/bin/bash
#$ -m bea
#$ -cwd #current working directory
#$ -M meixilin # username on hoffman
#$ -t 1-25

echo "job start" `date`
echo ${SGE_TASK_ID}

WORK="/u/project/rwayne/meixilin/abiotic_transect/"
RSCRIPT=${WORK}/r_codes/step4_gradient_forest/3_step4_gf_permutate/3_step4_gradient_forest_permutate_06060219.R
GFSOURCE=$1 # please use absolute path
GFNAME=$2
NREP=$3

source /u/local/Modules/default/init/bash
module load R/3.5.0

Rscript --vanilla ${RSCRIPT} ${GFSOURCE} ${GFNAME} ${NREP} ${SGE_TASK_ID}

echo "job done" `date`