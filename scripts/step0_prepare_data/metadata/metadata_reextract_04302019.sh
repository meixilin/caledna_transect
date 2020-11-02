#!/bin/bash
#$ -l highp,h_rt=06:00:00,h_data=16G
#$ -cwd #current working directory
#$ -N metadata_reextract_04302019
#$ -M meixilin # username on hoffman
#$ -o /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/logs/metadata_reextract_04302019.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/logs/metadata_reextract_04302019.err.txt
#$ -m bea

# this script is used to 
# 1. extract values from 
echo "job start."
echo `date`

source /u/local/Modules/default/init/bash
module load proj/5.2.0
module load gdal/2.3.2
module load geos/3.4.2
module load R/3.5.0

RSCRIPT=/u/home/m/meixilin/project-rwayne/abiotic_transect/r_codes/step0_prepare_data/metadata/generate/2_extract_raster.R
WORK=/u/home/m/meixilin/project-rwayne/abiotic_transect

cd ${WORK}
Rscript --vanilla ${RSCRIPT}

echo "job done." 
echo `date`
