#!/bin/bash
#$ -l h_data=2G,h_rt=23:00:00
#$ -wd /u/project/rwayne/meixilin/caledna_transect
#$ -o /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/all_100m_get_md5sum_20201019.out.txt
#$ -e /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs/all_100m_get_md5sum_20201019.err.txt
#$ -m abe

# @version 		v0
# @usage		qsub -N <jobname> util_get_md5sum.sh <FILEDIR> <OUTPREFIX>
# @description	get md5sum of a designated directory
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Aug 10 16:05:43 2020

###########################################################
## import packages 

###########################################################
## def functions 
set -euo pipefail 
###########################################################
## def variables 
FILEDIR=${1} # The directory to put the output md5sum checkfile
OUTPREFIX=${2} # The prefix of the md5sum check_file

###########################################################
## main 
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; FILEDIR ${FILEDIR}; OUTPREFIX ${OUTPREFIX}" 

cd ${FILEDIR}
find -type f -exec md5sum "{}" + > ${OUTPREFIX}.md5sum

# remove the line that calculate the md5sum of md5sum (exact match of name)
sed -i "/${OUTPREFIX}.md5sum$/d" ${OUTPREFIX}.md5sum
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"  
