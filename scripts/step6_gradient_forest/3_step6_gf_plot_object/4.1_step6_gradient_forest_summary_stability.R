# Title: get a gf and summarize the outputs --------
# pair with stability test 
# Author: Meixi Lin
# Date: Fri Oct 18 21:19:38 2019
# Author:
# Date:
# Modification:

# call from hoffman2 ---------
# module load R/3.5.0 
# Rscript --vanilla /u/home/m/meixilin/project-rwayne/abiotic_transect/r_codes/step4_gradient_forest/4_step4_gf_plot_object/4.1_step4_gradient_forest_summary_function_10182019.R > /u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step4_gradient_forest/2_final/stability_test.log 2>&1

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

# setwd("~/Lab/caledna_transect/")
setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {
    source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))
    return(0)
    })

# define a function for getting gf needed --------
get_gf_summary <- function(gfpath) {
    gf <- loadRData(gfpath)
    most_imp <- importance(gf)
    # print(most_imp)
    rsq <- gf_res(gf)
    # print(rsq)
    imp_families <- gf$result
    gf_summary <- list(most_imp, rsq, imp_families)
    names(gf_summary) <- c("most_imp", "rsq", "imp_families")
    class(gf_summary) <- "gf_sum"
    return(gf_summary)
}

# define variables --------
gfpathes <- c("./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData", 
              paste0("./derive_data/step4_gradient_forest/2_final/gf_deco_stab__all_Family_Presence_2000_2_0.05_FALSE_", 1:19, ".RData"))

# main --------
gf_sumlist <- lapply(gfpathes, function(xx) {
    return(get_gf_summary(xx))
})

save(gf_sumlist, file = "./derive_data/step4_gradient_forest/2_final/stab_res_gf_deco_all_Family_Presence_2000_2_0.05_FALSE.RData")
