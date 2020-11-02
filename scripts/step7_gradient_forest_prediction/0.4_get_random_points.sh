#! /bin/bash
# get random points for prediction plotting
# In UCLA Hoffman2 cluster 

source /u/local/Modules/default/init/bash
module load perl

echo "job starts" `date` 
# now go to the output section and sample 1% of the sites 
cd /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction
time perl -ne 'print if (rand() < .01)' myvar_ras_dataframe.csv > gf_rand_grids.csv
echo "job done" `date`
