#!/bin/bash
#$ -l highp,h_data=8G,h_rt=36:00:00
#$ -o /u/home/m/meixilin/project-rwayne/bioclim/get_bioclim.out.txt
#$ -e /u/home/m/meixilin/project-rwayne/bioclim/get_bioclim.err.txt
#$ -m bea
#$ -M meixilin

# get to directory 

cd /u/home/m/meixilin/project-rwayne/bioclim

# get bioclim variables 

wget --continue http://biogeo.ucdavis.edu/data/worldclim/v2.0/tif/base/wc2.0_30s_bio.zip


