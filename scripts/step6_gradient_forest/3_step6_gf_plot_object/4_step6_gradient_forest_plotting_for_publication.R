# Title: Plot the result from "step4_gradient_forest_1018" --------
# Author: Meixi Lin
# Date: Mon Sep  3 11:55:04 2018
# Author: Meixi Lin
# Date:
# Modification: Pair with the step4_gradient_forest_1018.R

# source for finding local maxima 
# https://stats.stackexchange.com/questions/30750/finding-local-extrema-of-a-density-function-using-splines

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
require(pastecs)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {
    source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))
    return(0)
})

# get arguments --------
gf <- loadRData("./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData") # RENAME WHAT EVER GF to "gf"
filename <- "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17" # this should be the name of the "gf" object tested
plotdir <- "./plots_important/step4_gradient_forest/for_publication/"

# get basic parameters --------
print_gf_summary(gf)
dir.create(plotdir, recursive = T)

# plot the standard plots --------
# (accuracy importance for each parameters) =================================
# start an importance plot dataframe
overimp1 <- sort(importance(gf, type = "Weighted"), decreasing = F) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "var") %>%
    dplyr::mutate(var = factor(var, levels = var))
colnames(overimp1) <- c("var", "imp")
overimp1 <- overimp1
overimp2 <- sort(gf$result, decreasing = F) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var = "var") %>%
    dplyr::mutate(var = factor(var, levels = var)) %>%
    dplyr::top_n(30, var)
colnames(overimp2) <- c("species", "imp")

pp1 <- ggplot(overimp1, aes(x = var, y = imp)) +
    geom_bar(stat = "identity", colour = "black", fill = "lightgrey") +
    coord_flip() +
    theme_bw() + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Predictor", y = bquote(R^2~' weighted importance'))

pp2 <- ggplot(overimp2, aes(x = species, y = imp)) +
    geom_point(colour = "black", fill = "khaki", shape = 21, size = 2) +
    coord_flip() +
    # geom_line(linetype = 2) +
    theme_bw() + 
    # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = 'Family', y = 'Goodness-of-fit')

ggsave(filename = paste0(plotdir, filename, "_A1.pdf"), plot = pp1, height = 4, width = 3.5)
ggsave(filename = paste0(plotdir, filename, "_A2.pdf"), plot = pp2, height = 4, width = 4)

# (importance density for 8 most important parameters) =================================
most_important <- names(importance(gf))[1:3]
pdf(file = paste0(plotdir, filename, "_B", ".pdf"), height = 4, width = 8)
dStdNorm <- mysplit.density.plot(obj = gf, imp.vars = most_important, mfrow = c(1,3), nbin = 50, barwidth = 2, cex.legend = 1.2, imp.vars.names = c("Elevation", "Sand Percent", "NDVI32"))
# get local maximum 
mylocalmax <- lapply(dStdNorm, local_maxi)
print(mylocalmax)

dev.off()

# (predictor cumulative plot, not show speciess) =================================
# predictor cumulative plot (plot.type="C", show.species=F),
# which for each predictor shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species.
pdf(file = paste0(plotdir, filename, "_C", ".pdf"), height = 4, width = 8)
myspecies.cumulative.plot(obj = gf, imp.vars = most_important, imp.vars.names = c("Elevation", "Sand Percent", "NDVI32"), show.species = F, common.scale = T, mfrow = c(1,3))
dev.off()

# # output from local maxima at density ========
# $elev
# x            y Acut
# 1   10.97577 4.933017e-05    1
# 2  228.31989 7.929962e-06    0
# 3  341.85786 1.452442e-05    0
# 4  555.95804 1.256882e-05    0
# 5  740.86274 1.332491e-05    0
# 6 1039.30542 4.022896e-06    0
# 7 1220.96618 1.751285e-05    1
# 8 1431.82241 4.141573e-05    1
# 9 1645.92259 3.434640e-05    1
# 
# $sndppt
# x            y Acut
# 1  22.79557 0.0012790078    1
# 2  30.01128 0.0003292856    1
# 3  32.93940 0.0003749148    1
# 4  38.69105 0.0002966686    1
# 5  42.97865 0.0003990798    1
# 6  49.46233 0.0001989059    0
# 7  52.39045 0.0001908370    0
# 8  56.99177 0.0001720392    0
# 9  64.10291 0.0001268323    0
# 10 73.51471 0.0005205619    1
# 
# $NDVI32
# x          y Acut
# 1  -0.15640137 0.03636337    1
# 2  -0.03713826 0.01574690    1
# 3   0.05183581 0.03879711    1
# 4   0.16920586 0.01499688    0
# 5   0.27711058 0.01963544    1
# 6   0.38690837 0.01258408    0
# 7   0.48913390 0.01177827    0
# 8   0.58946636 0.01200945    0
# 9   0.62922073 0.01241226    0
# 10  0.70683641 0.01478348    0