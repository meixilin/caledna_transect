#!/bin/bash
##$ -l highp,h_data=16G,h_rt=06:00:00
##$ -o /u/home/m/meixilin/project-rwayne/GeoTiff/rasters/soilgrid_x100m/get_soilgrid_100m.out.txt
##$ -e /u/home/m/meixilin/project-rwayne/GeoTiff/rasters/soilgrid_x100m/get_soilgrid_100m.err.txt
##$ -m bea
##$ -M meixilin

# this is executed interactively 

# bash get_soilgrid_100m.sh

# get to directory 

# cd /u/home/m/meixilin/project-rwayne/GeoTiff/rasters/soilgrid_x100m/

# get soil grid variables and rename them 
# soil organic carbon
wget https://scholarsphere.psu.edu/downloads/1zc77sn34p 
mv 1zc77sn34p soc_M_sl1_100m.tif

# sand percentage 
wget https://scholarsphere.psu.edu/downloads/4tm70ms91p
mv 4tm70ms91p sand_M_sl1_100m.tif

# pH in water 
wget https://scholarsphere.psu.edu/downloads/cnp1938923
mv cnp1938923 ph_h2o_M_sl1_100m.tif

# total nitrogen
wget https://scholarsphere.psu.edu/downloads/5qf85n939n
mv 5qf85n939n n_tot_M_sl1_100m.tif

# clay weight percentage 
wget https://scholarsphere.psu.edu/downloads/xwd375x33w
mv xwd375x33w clay_M_sl1_100m.tif

# bulk density 
wget https://scholarsphere.psu.edu/downloads/mkk91fm281
mv mkk91fm281 bd_M_sl1_100m.tif

# HERE STARTS PRELIMINARY UNVALIDATED maps 
# Mg 
wget https://scholarsphere.psu.edu/downloads/tx920fx36s
mv tx920fx36s mg_pct_M_sl1_100m.tif

# K 
wget https://scholarsphere.psu.edu/downloads/573666467v
mv 573666467v k_pct_M_sl1_100m.tif

# electrical conductivity
wget https://scholarsphere.psu.edu/downloads/5qf85n933s
mv 5qf85n933s ec_12pre_M_sl1_100m.tif

# the other two were: 
# TAXOUSDA
wget https://files.isric.org/soilgrids/data/recent/TAXOUSDA_250m_ll.tif

# CECSOL
wget https://files.isric.org/soilgrids/data/recent/CECSOL_M_sl1_250m_ll.tif


