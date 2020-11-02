#!/usr/bin/env bash

# for finalizing paramter settings

WORK="/u/project/rwayne/meixilin/abiotic_transect/"
# WORK="/Users/linmeixi/UCLA/Lab/abiotic_transect/"
SCRIPT=${WORK}/r_codes/step4_gradient_forest/3_step4_gf_permutate/3_step4_gradient_forest_permutate_qsub_06060219.sh
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub
DATE=06062019
LOG=${WORK}/derive_data/step4_gradient_forest/3_permutate/logs
mkdir -p ${LOG}

# date: 06062019
# for this run on gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData
OUTNAME=gf_permutate_all
GFSOURCE=${WORK}/derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData # please use absolute path
GFNAME=gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17
NREP=4
SETTING="-l highp,h_rt=08:00:00,h_data=48G"

qsub -N ${OUTNAME}_${DATE} \
${SETTING} \
-o ${LOG}/${OUTNAME}_${DATE}.out.txt \
-e ${LOG}/${OUTNAME}_${DATE}.err.txt \
${SCRIPT} ${GFSOURCE} ${GFNAME} ${NREP}

# date: 10192019 ==================================
# for run without coast 
WORK=/u/project/rwayne/meixilin/abiotic_transect
# WORK="/Users/linmeixi/UCLA/Lab/abiotic_transect/"
SCRIPT=${WORK}/r_codes/step4_gradient_forest/3_step4_gf_permutate/3_step4_gradient_forest_permutate_qsub_06060219.sh
QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub
DATE=10192019
LOG=${WORK}/derive_data/step4_gradient_forest/3_permutate/logs

OUTNAME=gf_nocoast
GFSOURCE=${WORK}/derive_data/step4_gradient_forest/2_final/gf_nocoast_all_Family_Presence_2000_2_0.05_FALSE_17.RData # please use absolute path
GFNAME=gf_nocoast_all_Family_Presence_2000_2_0.05_FALSE_17
NREP=4
SETTING="-l highp,h_rt=08:00:00,h_data=48G"

qsub -N ${OUTNAME}_${DATE} \
${SETTING} \
-m bea \
-o ${LOG}/${OUTNAME}_${DATE}.out.txt \
-e ${LOG}/${OUTNAME}_${DATE}.err.txt \
${SCRIPT} ${GFSOURCE} ${GFNAME} ${NREP}

OUTNAME=gf_coast
GFSOURCE=${WORK}/derive_data/step4_gradient_forest/2_final/gf_coast_all_Family_Presence_2000_2_0.05_FALSE_17.RData # please use absolute path
GFNAME=gf_coast_all_Family_Presence_2000_2_0.05_FALSE_17
NREP=4
SETTING="-l highmem,h_rt=08:00:00,h_data=48G"

qsub -N ${OUTNAME}_${DATE} \
${SETTING} \
-m bea \
-o ${LOG}/${OUTNAME}_${DATE}.out.txt \
-e ${LOG}/${OUTNAME}_${DATE}.err.txt \
${SCRIPT} ${GFSOURCE} ${GFNAME} ${NREP}
