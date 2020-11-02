###########################################################
cd /u/project/rwayne/meixilin/caledna_transect/scripts
###########################################################
# Wed Oct 14 13:34:55 2020
qsub -t 51 step7_gradient_forest_prediction/all_100m/1_step7_make_prediction_noextrap_qsub.sh
qsub -t 1-50 step7_gradient_forest_prediction/all_100m/1_step7_make_prediction_noextrap_qsub.sh
# check job completion 
cd /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs
for ii in {00..50}; do
grep '] JOB ID .* Done' 1_step7_make_prediction_noextrap_20201013_split${ii}.log
done 

###########################################################
# Wed Oct 14 16:34:56 2020
qsub -t 1-51 step7_gradient_forest_prediction/all_100m/1_step7_make_prediction_extrap_qsub.sh

# check job completion 
cd /u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/logs
for ii in {00..50}; do
grep '] JOB ID .* Done' 1_step7_make_prediction_extrap_20201013_split${ii}.log
done 

###########################################################
# Fri Oct 16 09:48:49 2020
qsub step7_gradient_forest_prediction/all_100m/2_step7_combine_prediction_noextrap_qsub.sh

###########################################################
# Fri Oct 16 15:48:49 2020
# Max vmem         = 30.326G
qsub step7_gradient_forest_prediction/all_100m/2_step7_combine_prediction_extrap_qsub.sh

###########################################################
# Sun Oct 18 23:35:34 2020
qsub step7_gradient_forest_prediction/all_100m/3_step7_prediction_pc_noextrap_qsub.sh
# JOB ID: 4798978; Max vmem         = 94.343G
qsub step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_noextrap_qsub.sh

###########################################################
# Sun Oct 18 23:35:42 2020
qsub step7_gradient_forest_prediction/all_100m/3_step7_prediction_pc_extrap_qsub.sh
# JOB ID: 4798996; Max vmem         = 92.993G
qsub step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_extrap_qsub.sh

###########################################################
# Mon Oct 19 01:17:27 2020
# JOB ID: 4798979 had a warning, modify the memory usage (commit 9d3259cab3c40e3dc19714fda0c6514e6ba6f0c2)
qsub step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_noextrap_qsub.sh

###########################################################
# Mon Oct 19 11:15:47 2020
# commit fa94543
qsub step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_qsub.sh

###########################################################
# Mon Oct 19 16:19:46 2020
# commit 7085f23
qsub step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_qsub.sh

###########################################################
# Mon Oct 19 16:41:31 2020
# commit d275502
qsub step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_extrap_qsub.sh

###########################################################
# Mon Oct 19 22:24:48 2020
# commit df67d15
qsub step7_gradient_forest_prediction/all_100m/4.3_step7_store_ref_scale_grid_qsub.sh

qsub step7_gradient_forest_prediction/all_100m/5_step7_store_prediction_extrap_qsub.sh

###########################################################
# Mon Oct 19 22:31:57 2020
FILEDIR="/u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/extrap_FALSE"
OUTPREFIX="predict_gf_Trns_CAgrid_noextrap"
qsub step7_gradient_forest_prediction/all_100m/util_get_md5sum.sh ${FILEDIR} ${OUTPREFIX}

###########################################################
# Mon Oct 19 23:59:21 2020
# commit f3214d1
qsub step7_gradient_forest_prediction/all_100m/4.2_step7_ref_scale_grid_plot_qsub.sh

###########################################################
# Tue Oct 20 00:30:36 2020
FILEDIR="/u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/extrap_TRUE"
OUTPREFIX="predict_gf_Trns_CAgrid_extrap"
qsub step7_gradient_forest_prediction/all_100m/util_get_md5sum.sh ${FILEDIR} ${OUTPREFIX}

FILEDIR="/u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/nobio_CAgrid"
OUTPREFIX="predict_nogf_Scld_CAgrid"
qsub step7_gradient_forest_prediction/all_100m/util_get_md5sum.sh ${FILEDIR} ${OUTPREFIX}

###########################################################
# Tue Oct 20 00:32:26 2020
# ALL DONE. 

###########################################################
# Tue Oct 20 16:26:00 2020
# commit 5f9ab68
qsub step7_gradient_forest_prediction/all_100m/util_compare_previous_qsub.sh
