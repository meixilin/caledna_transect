

# CALeDNA transect analyses

*This repository stores the scripts and important data used to generate the analyses in "Landscape Analyses Using eDNA Metabarcoding and Earth Observation Predict Community Biodiversity in California" by Lin et al. Submitted.*

*The complete repository including raw data, maps and rasters, derived data and intermediate plots can be found at zenodo archive [link pending]().*

*The sequencing results are deposited in NCBI SRA [link pending ]()*

Contact: Meixi Lin (meixilin[at]ucla[d0t]edu)

# File structure

```
.
├── README.md
├── raw_data
│   ├── TAXO
│   ├── metadata
│   ├── transect_ASV
│   └── ucnrs
├── final_data
│   ├── Final_metadata.RData
│   ├── deco_3
│   ├── metadata
│   └── ucnrs
├── maps_vectors_rasters
│   ├── rasters
│   └── vectors
├── derive_data
│   ├── phy_deco
│   ├── phy_rare
│   ├── step0_prepare_data
│   ├── step1_create_phyloseq
│   ├── step2_data_description
│   ├── step3_alpha_diver
│   ├── step4_beta_diver
│   ├── step6_gradient_forest
│   └── step7_gradient_forest_prediction
├── plots
│   ├── step0_prepare_data
│   ├── step1_create_phyloseq
│   ├── step2_data_description
│   ├── step3_alpha_diver
│   ├── step4_beta_diver
│   ├── step5_zeta_diver
│   ├── step6_gradient_forest
│   └── step7_gradient_forest_prediction
└── scripts
    ├── function_transect.R
    ├── install_packages.R
    ├── other_scripts
    ├── step0_prepare_data
    ├── step1_create_phyloseq
    ├── step2_data_description
    ├── step3_alpha_diver
    ├── step4_beta_diver
    ├── step5_zeta_diver
    ├── step6_gradient_forest
    └── step7_gradient_forest_prediction
```

# Scripts

Please adapt the scripts input and output directory accordingly to your working directories. Important derived data can be found in the zenodo archive as well. 

## Set up the environment 

1. Computational environment
   1. Most downstream analyses can be conducted on a PC/Mac. 
   2. Some analyses (e.g. sequence processing and gradient forest analyses) were performed on [UCLA Hoffman2 cluster program](https://www.hoffman2.idre.ucla.edu/)
   3. Some environmental layers were obtained using [Google Earth Engine](https://code.earthengine.google.com/)

2. Most analyses are performed using R. The required R packages are listed in [install_packages.R](scripts/install_packages.R) 

3. The utility script that records necessary variables and functions can be sourced at  [function_transect.R](scripts/function_transect.R) 

## Step0: prepare data

```
.
├── asv_tables
│   ├── 1_Transect_data_decontamination.R
│   ├── 2_remove_5_dup_sites_1028.R
│   └── 3_import_taxonomy_detailed_10172019.R
├── metadata
│   ├── 1_reproject_raster.R
│   ├── 2_extract_raster.R
│   ├── 3_calculate_geographical_dist_0827.R
│   ├── 3_combine_old_meta.R
│   ├── 4_specify_biom_RData_0412.R
│   ├── 5_update_coastal_0531.R
│   ├── 6.1_correlation_test_metadata_1026.R
│   ├── 6.2_correlation_grouping_metadata_0523.R
│   ├── 7.1_density_distribution.R
│   ├── function_raster.R
│   ├── get_bioclim.sh
│   ├── get_earth_engine.js
│   ├── get_land_cover.js
│   ├── get_soilgrid_100m.sh
│   ├── get_uncertainty_earth_engine.js
│   ├── get_uncertainty_percent_earth_engine.js
│   ├── metadata_reextract_04302019.sh
│   ├── metadata_reproject_04032019.sh
│   ├── x.2_get_other_uncertainty_20200930.R
│   └── x_get_soilgrid_uncertainty_20200904.R
└── sample_map
    ├── ca_base_map.R
    └── sample_map.R
```

### Usage: 

1. decontaminate asv tables and remove five replicate sites 
2. obtain and extract environmental layers

### Output:

1. final taxonomy tables and metadata 
2. environmental layers with Statewide coverage 

## Step1: create phyloseq

```
.
├── 1_step1_read_newdeco_data_0412.R
├── 2_step1_rarefaction_1107.R
└── rarefaction_evaluation
    ├── 1_rarefaction_plot_0506.R
    ├── 2_evaluate_rarefaction_0529.R
    └── 3_iNext_evaluation.R
```

### Usage: 

1. create `phyloseq` object from final dataset

2. perform rarefaction evaluation and select the appropriate rarefaction depth. 


### Output: 

1. decontaminated phyloseq objects
2. rarefied phyloseq objects

```bash
├── derive_data
│   ├── phy_deco # phyloseq objects of each decontaminated metabarcode dataset 
│   │   ├── phy_deco_16S.RData
│   │   ├── phy_deco_18S.RData
│   │   ├── phy_deco_CO1.RData
│   │   ├── phy_deco_FITS.RData
│   │   ├── phy_deco_PITS.RData
│   │   ├── phy_deco_all.RData
│   │   ├── phydeco.RData
│   │   └── phydeco_uc.RData
│   ├── phy_rare # phyloseq objects of each rarefied metabarcode dataset 
│   │   ├── phy_1000_CO1.RData
│   │   ├── phy_1000_PITS.RData
│   │   ├── phy_2000_16S.RData
│   │   ├── phy_4000_18S.RData
│   │   ├── phy_4000_FITS.RData
│   │   └── phyrare.RData
│   ├── step1_create_phyloseq # evaluation output 
│   │   ├── eval_rarefaction
│   │   └── physeq_rare_otu 
```

## Step2: data description

```
.
├── concordance
│   ├── 1_summary_5_dup_sites_20200107.R
│   ├── 2_describe_5_dup_sites_20200107.R
│   ├── 2_describe_5_dup_sites_full_20200107.R
│   ├── 3_summary_plots.R
│   ├── 3_summary_plots_cutoff.R
│   ├── 4_pcoa_dup_sites_supp_fig6.R
│   └── describe_5_dup_sites_call.sh
├── gbif_tos
│   └── 1_check_gbif_tos_20200928.R
├── summary_stats
│   ├── 1_decontam_seqs_read_depths.R
│   ├── 2.1_metadata_distributions.R
│   ├── 2.2_raster_data_distribution.R
│   ├── 3_category_correlation.R
│   ├── 4_step1_physeq_basic_stats_0830.R
│   ├── 4_step1_physeq_basic_stats_rarefied.R
│   ├── 5_step1_physeq_plot_stats_06122019.R
│   ├── 6_step1_physeq_plot_by_sample_stats_09022019.R
│   ├── 7_README_step1_physeq_phylo_tree_plot.md
│   ├── 7_step1_physeq_phylo_tree_plot.R
│   └── 7_step1_physeq_phylo_tree_plot.sh
└── ucnrs
    ├── 1_mk_phyloseq_ucnrs_sites.R
    ├── 2.1_clean_ucnrs_records_format_typo.R
    ├── 2.2_clean_ucnrs_taxlevel_genus.R
    ├── 2.3_clean_ucnrs_taxlevel_fillin.R
    ├── 3.1_compare_ucnrs_subset_sites.R
    ├── 3.2_compare_ucnrs_intersect.R
    ├── 3.3_compare_ucnrs_overall.R
    ├── README_ucnrs.md
    ├── README_ucnrs.pdf
    └── x.1_match_loc_ucnrs_name.R
```

### Usage: 

1. Perform evaluation of stability by comparing the concordance of the five replicated sites
2. Generate summary statistics such as metadata distribution, read depth, read count per taxa, etc. 
3. Perform cleaning of ucnrs records and comparison of ucnrs records with eDNA

### Output:

1. Concordance analyses phyloseq objects and results
2. Summary statistics: metadata distribution, sequencing read distribution and taxonomic coverage distribution
3. UCNRS analyses: UCNRS cleaned species list and comparison results 

## Step3: alpha diversity 

```
.
├── 1.2_step2_alpha_diversity_summary_20200113.R
├── 1_step2_alpha_diversity_calc_0530.R
├── 2.1_step2_alpha_diversity_kruskal_20200113.R
├── 2.2_step2_alpha_diversity_kruskal_eval_20200113.R
├── 3.1_step2_alpha_diversity_individual_lm_20200125.R
├── 3.2_step2_alpha_diversity_indi_lm_plot_FITS_20200504.R
├── 4.1_step2_alpha_diversity_reduce_pls_20200125.R
├── 5.1_step2_alpha_diversity_plot_fig2.R
├── 5.2_step2_alpha_diversity_plot_supp_fig7_loc.R
├── 5.3_step2_alpha_diversity_plot_supp_fig8.R
└── 6.1_step2_alpha_diversity_map.R
```

### Usage: 

1. Calculate `observed` and `Shannon index` for rarefied dataset. 
2. Perform Kruskal Walis testing on the categorical variables. 
3. Perform individual linear regression on the continuous variables. 
4. Perfrom partial least square analyses on the continuous variables. 

### Output:

1. Alpha diversity values.
2. Test results. 

## Step4: beta diversity 

```
.
├── 0_util_get_var_color.R
├── 0_util_query_var_color.R
├── 1_step4_beta_diversity_cal_diss_60152019.R
├── 2.1_step4_beta_diversity_lcbd_06242019.R
├── 3.1_step4_beta_diversity_adonis_0804.R
├── 3.2_step4_beta_diversity_out_adonis_0804.R
├── 4.1_step4_beta_diversity_pcoa_generate_plotobject.R
├── 5.1_step4_beta_diversity_pcoa_envfit.R
├── 5.2_step4_beta_diversity_pcoa_envfit_ordisurf_plot.R
├── 6.1_step4_beta_diversity_cap_generate_plotobject.R
├── 6.2_step4_beta_diversity_out_cap_0805.R
├── manuscript_plotting
│   ├── m_step4_beta_diversity_combined_lcbd_fig3.R
│   ├── m_step4_beta_diversity_combined_ordinations_fig4.R
│   ├── m_step4_beta_diversity_combined_ordinations_resub_fig2.R
│   ├── m_step4_beta_diversity_combined_ordinations_supp_figs.R
│   └── m_step4_beta_diversity_pcoa_envfit_ordisurf_plot_fig6.R
├── no_coast
│   ├── x.1_step3_beta_diversity_cal_diss_nocoast_60152019.R
│   ├── x.2_step3_beta_diversity_nocoast_adonis_0804.R
│   ├── x.3_step3_beta_diversity_nocoast_out_adonis_0804.R
│   ├── x.4_step4_beta_diversity_nocoast_pcoa_plot_0805.R
│   └── x.5_step3_beta_diversity_pcoa_combine_0805.R
└── zoom_in
    ├── z.1_step4_beta_diversity_adonis_minor_in_major_09302019.R
    └── z.2_step4_beta_diversity_pcoa_minor_in_major_09302019.R
```

### Usage: 

1. Calculate binary Jaccard dissimilarity 
2. Generate relative abundance plots
3. Perform beta dispersion and PERMANOVA analysis
4. Perform PCoA ordination and plotting 
5. Perform envfit on the PCoA ordination *post hoc*
6. Perform additional beta dispersion and PERMANOVA analysis on samples 1) excluded coastal sites 2) according to minor habitats within major habitats 
7. Perform cap scale and varpart analyses to remove location effects 

### Output:

1. Jaccard dissimilarity measures
2. Analyses output 

## Step5: zeta diver

Author: Ariel Levi Simons. His [github repository](https://github.com/levisimons/CALeDNA) for this section

```
.
├── CALeDNAZetaFactors.R
├── LICENSE
├── PresenceAbsence.R
├── README.md
├── Zeta4ClassMap.R
├── Zeta4ClassMapCluster.R
├── Zeta4FamilyMap.R
├── Zeta4FamilyMapCluster.R
├── Zeta4eDNA.R
├── ZetaClassCluster.R
├── ZetaFamilyCluster.R
├── ZetaMap.R
├── ZetaeDNA.R
├── step5_glom_to_family.R
├── step5_subset_phyloseq.R
└── zeta.sh
```
### Usage: 

1. Calculate $\zeta_4$ diversity at family level 
2. Perform variable importance testing 

## Step6: gradient forest modeling

```
.
├── 1_step6_gradient_forest_final
│   ├── 1_step6_gradient_forest_final_all_06060219.R
│   ├── 1_step6_gradient_forest_final_all_qsub_06060219.sh
│   ├── 1_step6_gradient_forest_final_coast_10182019.R
│   ├── 1_step6_gradient_forest_final_coast_qsub_10182019.sh
│   ├── 1_step6_gradient_forest_final_nocoast_10182019.R
│   ├── 1_step6_gradient_forest_final_nocoast_qsub_10182019.sh
│   ├── 2_step4_gradient_forest_final_run_stability_10152019.R
│   └── 2_step4_gradient_forest_final_run_stability_10152019.sh
├── 2_step6_gf_permutate
│   ├── 3.2_plot_permutation_result.R
│   ├── 3_step4_gradient_forest_permutate_06060219.R
│   ├── 3_step4_gradient_forest_permutate_qsub_06060219.sh
│   └── 3_step4_gradient_forest_permutate_qsub_wrapper_06060219.sh
├── 3_step6_gf_plot_object
│   ├── 4.1_step6_gradient_forest_summary_stability.R
│   ├── 4.1_step6_gradient_forest_summary_stability_output.R
│   ├── 4.2_plot_permutation_result.R
│   ├── 4.3_validate_reads_10172019.R
│   ├── 4.4_check_occurrences_predictability.R
│   ├── 4_step6_gradient_forest_plotting_for_publication.R
│   ├── 4_step6_gradient_forest_plotting_function_local.R
│   ├── 4_step6_gradient_forest_species_response_10162019.R
│   ├── other_save_gf_x_y_manuscript.R
│   └── species_response_legend_key.R
└── functions
    ├── change_density_plot.R
    ├── change_performance.plot.R
    ├── change_species_cumulative_plot.R
    └── function_gfprep.R
```

### Usage: 

1. Perform gradient forest data preparation 
2. Perform gradient forest analyses using classification tree setting
3. Perform gradient forest validation
3. Plot gradient forest output 

### Output:

1. gradient forest input table
2. gradient forest result 

## Step7: gradient forest predictions

```
.
├── 0.1_prepare_latlong_raster.R
├── 0.2_prepare_raster.R
├── 0.2_prepare_raster_qsub.sh
├── 0.3_prepare_raster_splits.R
├── 0.3_prepare_raster_splits_qsub.sh
├── 0.4_get_random_points.sh
├── all_100m
│   ├── 1_step7_make_prediction_all_20201014.R
│   ├── 1_step7_make_prediction_extrap_qsub.sh
│   ├── 1_step7_make_prediction_noextrap_qsub.sh
│   ├── 2_step7_combine_prediction_20201014.R
│   ├── 2_step7_combine_prediction_extrap_qsub.sh
│   ├── 2_step7_combine_prediction_noextrap_qsub.sh
│   ├── 3_step7_prediction_pc_20201014.R
│   ├── 3_step7_prediction_pc_extrap_qsub.sh
│   ├── 3_step7_prediction_pc_noextrap_qsub.sh
│   ├── 4.1_step7_ref_scale_grid_20201014.R
│   ├── 4.1_step7_ref_scale_grid_qsub.sh
│   ├── 4.2_step7_ref_scale_grid_plot_20201014.R
│   ├── 4.2_step7_ref_scale_grid_plot_qsub.sh
│   ├── 4.3_step7_store_ref_scale_grid_20201014.R
│   ├── 4.3_step7_store_ref_scale_grid_qsub.sh
│   ├── 5_step7_store_prediction_20201014.R
│   ├── 5_step7_store_prediction_extrap_qsub.sh
│   ├── 5_step7_store_prediction_noextrap_qsub.sh
│   ├── record_all_100m_predictions.sh
│   ├── util_compare_previous.R
│   ├── util_compare_previous_qsub.sh
│   └── util_get_md5sum.sh
├── manuscript_plotting
│   ├── 1_fig9_nobio_prediction_pc_plots.R
│   ├── 1_fig9_noextrap_prediction_pc_plots.R
│   └── record_manuscript_plotting_20201025.sh
└── sites_rand_points
    ├── 1_step7_make_prediction_gf_sites.R
    ├── 2_step7_make_prediction_random_sites.R
    ├── 3_step7_ref_scale_sites_random_sites_20200421.R
    ├── 3_step7_ref_scale_sites_random_sites_20200421_qsub.sh
    ├── 4_evaluate_prediction.R
    └── 5_evaluate_prediction_uncertainty.R
```

### Usage: 

1. Perform predictions based on the gradient forest output

### Output:

1. Stacked layers of 33 environmental variables in California at 100 m x 100 m resolution. 
2. Community turnover map predicted from the gradient forest model. 

