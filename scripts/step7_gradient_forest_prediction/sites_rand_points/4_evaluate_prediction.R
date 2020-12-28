# Title: quantifying the performance of biologically informed mapping --------
# Author: Meixi Lin
# Date: Tue Jun 11 13:31:12 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
options(echo = T)
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# # if you need to set things locally 
setwd("~/Lab/caledna_transect/")
library(gradientForest)
library(vegan)
library(MASS)

# some functions 
source("./scripts/function_transect.R")
# source("./r_codes/step0_prepare_data/metadata/generate/function_raster.R")
# source("./r_codes/step6_gradient_forest_prediction/function_gfpredict.R")
# functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
# lapply(functions.file, function(xx) {
#     source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))
#     return(0)
# })

# define a function for calculating ED --------
calc_ED <- function(datadir) {
    load(file =  datadir)
    ED <- stats::dist(Trns_site, method = "euclidean")
    return(ED)
}

calc_ED_all <- function(datadir, gfname, primer, id) {
    load(file = paste0(datadir, gfname, "/predict_gf_output_", primer,".RData"))
    Trns_grid <- Trns_grid[id,]
    ED <- stats::dist(Trns_grid, method = "euclidean")
    return(ED)
}

calc_stress <- function(shep) {
    stress <- sum((shep$y - shep$yf)^2)/sum(shep$y^2)
    return(stress)
}

plot_shep <- function(shep, xxlab, yylab) {
    stress <- sum((shep$y - shep$yf)^2)/sum(shep$y^2)
    rstress <- 1 - stress
    ralscal <- cor(shep$y, shep$yf)^2
    plot(shep, pch = ".", col = "lightblue", xlab = xxlab, ylab = yylab)
    if (length(shep) == 3) {
        lines(shep$x, shep$yf, type = "S", col = "red", lwd = 1)
    }
    lab <- paste("Stress =", format(stress, digits = 5),
                 "\nNon-metric fit, R2 =", format(rstress, digits=3),
                 "\nLinear fit, R2 =", format(ralscal, digits=3))
    text(min(shep$x), 0.90*max(shep$y), lab, pos=4)
    return(stress)
}

# set variables -------
gfname <- "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
indir <- paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/")
gfdir <- "./derive_data/step6_gradient_forest/"

plotdir <- "./plots/step7_gradient_forest_prediction/"
extrap <- FALSE
dir.create(paste0(plotdir, gfname, "/", extrap), recursive = T)

# real --------
gf <- loadRData(file = paste0(gfdir, gfname, ".RData"))
bio.data <- apply(gf$Y, 2, function(xx) {yy = as.numeric(xx == "TRUE")})
rownames(bio.data) <- rownames(gf$Y)
bio.diss <- vegan::vegdist(bio.data, method = "jaccard")

# convert "ED" matrix of RF transformed sites --------
load(file =  "./derive_data/step6_gf_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/FALSE/predict_gf_output_all.RData")
# load(file =  "./derive_data/step6_gf_prediction/gf_nocoast_all_Family_Presence_2000_2_0.05_FALSE_17/FALSE/predict_gf_output_all.RData")
rm(list = c('Trns_grid', 'forplot'))
ED.bio.no <- stats::dist(Trns_site, method = "euclidean")
rm(Trns_site)

load(file =  "./derive_data/step6_gf_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/TRUE/predict_gf_output_all.RData")
rm(list = c('Trns_grid', 'forplot'))
ED.bio.ex <- stats::dist(Trns_site, method = "euclidean")

load(file =  "./derive_data/step6_gf_prediction/nobio/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/predict_gf_output_nobio.RData")
rm(list = c('Scld_grid', 'forplot'))
ED.nobio <- stats::dist(Scld_site, method = "euclidean")

# check if the rowname is the same --------
table(rownames(Trns_site) == rownames(bio.data))
table(rownames(Scld_site) == rownames(bio.data))

# perform stress analysis -------
shep1 <- isoreg(x = ED.bio.no, y = bio.diss)
shep2 <- isoreg(x = ED.bio.ex, y = bio.diss)
shep3 <- isoreg(x = ED.nobio, y = bio.diss)
shep4 <- isoreg(x = ED.nobio, y = ED.bio.no)

mantel(xdis = ED.bio.no, y = bio.diss) # r = 0.2747; p = 0.001
mantel(xdis = ED.bio.ex, y = bio.diss) # r = 0.2692; p = 0.001
mantel(xdis = ED.nobio, y = bio.diss) # r = 0.2601; p = 0.001 
mantel(xdis = ED.nobio, y = ED.bio.no)

# (0.2747-0.2601)/0.2601 # mantel 5.6% increase
# (0.023896-0.024233)/0.024233 # stress performance 1.4% decrease (improvement)

PCsites <- predict(PCs, Trns_site)
ED.bio.no.pc <- stats::dist(PCsites[,1:2], method = "euclidean")
mantel(xdis = ED.bio.no, y = ED.bio.no.pc)
mantel(xdis = ED.bio.no.pc, y = bio.diss)

# plot the result --------
pdf(paste0(plotdir, gfname, "/evaluate_performance_noext.pdf"), width = 10, height = 10)
pp <- par(mfrow= c(2,2))
plot_shep(shep = shep1, xxlab = "ED of GF transformed environment (extrap = F)", yylab = "Jaccard dissimilarity")
plot_shep(shep = shep2, xxlab = "ED of GF transformed environment (extrap = T)", yylab = "Jaccard dissimilarity")
plot_shep(shep = shep3, xxlab = "ED of uninformed environment", yylab = "Jaccard dissimilarity")
plot_shep(shep = shep4, xxlab = "ED of uninformed environment", yylab = "ED of GF transformed environment (extrap = F)")
par(pp)
dev.off()

# load the overall sites --------
# NOT RUN ########
pdf(paste0(plotdir, gfname, "/evaluate_performance_env_var.pdf"), width = 6, height = 6)
plotid <- sample(1:505012, 500)
ED_all.bio <- calc_ED_all(datadir, gfname, primer = "all", id = plotid)
ED_all.nobio <- calc_ED_all(datadir = "./derive_data/step6_gf_prediction/nobio/",
                    gfname = "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17", 
                    primer = "nobio",
                    id = plotid)
shep5 <- isoreg(x = ED_all.nobio, y = ED_all.bio)
plot_shep(shep = shep5, xxlab = "ED of raw grid", yylab = "ED of transformed grid")
dev.off()
