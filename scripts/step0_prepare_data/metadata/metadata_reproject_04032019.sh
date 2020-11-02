#!/bin/bash
#$ -l highmem,highp,h_rt=20:00:00,h_data=128G
#$ -cwd #current working directory
#$ -N metadata_reproject_04032019
#$ -M meixilin # username on hoffman
#$ -o /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/logs/metadata_reproject_04032019.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/logs/metadata_reproject_04032019.err.txt
#$ -m bea
#$ -t 1-59

# this script is used to 
# 1. aggregate files with higher resolution to 100m
# 2. reproject to align the grids, using Earth Engine derived raster as a reference 
# 3. report the final directory of raster
echo "job start."
echo `date`

source /u/local/Modules/default/init/bash
module load proj/5.2.0
module load gdal/2.3.2
module load R/3.5.0

RSCRIPT=/u/home/m/meixilin/project-rwayne/abiotic_transect/r_codes/step0_prepare_data/metadata/generate/1_reproject_raster.R
WORK=/u/home/m/meixilin/project-rwayne/abiotic_transect

cd ${WORK}
Rscript --vanilla ${RSCRIPT} ${SGE_TASK_ID}

echo "job done." 
echo `date`
