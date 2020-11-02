# Title: scale the variables as comparison for gradient forest predictions # Perform Principal Component Analyses on the scaled CA_grid (step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_20201014.R) --------
# Author: Meixi Lin
# Date: Thu Oct 31 14:04:57 2019
# Author: Meixi Lin
# Date: Mon Oct 19 23:13:12 2020
# Modification: check before publications

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/u/project/rwayne/meixilin/caledna_transect/")
library(data.table)
library(sp)
library(raster)
source("./scripts/function_transect.R")
print(sessionInfo())

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
# mark the transformation as "nobio_CAgrid"
nobio = "nobio_CAgrid"
# extent for plotting 
myextent <- raster::extent(CAlimit)
myseed = 17
today = format(Sys.Date(), "%Y%m%d")

# the gradient forest setting
gfdir = "./derive_data/step6_gradient_forest/"
gfname = "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"

# output directory 
outdir = paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/", nobio, "/")
# the Scld_CA filename (step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_20201014.R)
Scld_CA_filepath = paste0(outdir, "predict_nogf_Scld_CAgrid.RData")
# the scaled sample sites (step7_gradient_forest_prediction/sites_rand_points/3_step7_ref_scale_sites_random_sites_20200421.R)
Scld_site_filepath = paste0(outdir, "predict_nogf_Scld_site.RData")

# load data --------
# load the Scaled CAgrid
print(paste(date(), "Loading Scld_CA ...", Scld_CA_filepath))

load(file = Scld_CA_filepath) # should be named as "Scld_CA"
print(str(Scld_CA))
Scld_site=loadRData(Scld_site_filepath) # the scaled sample sites 
Scld_site=data.frame(Scld_site) # convert to data frame
print(str(Scld_site)) 

# get PC and convert to RGB --------
print(paste(date(), "Principal Component Analyses ..."))

# run PCA on all the variables available (the first two cols are 'lon' 'lat', after that the variables should have been ordered by relative importance in previous scale script)
PCs = prcomp(Scld_CA[,3:ncol(Scld_CA)])
print(summary(PCs))

print(paste(date(), "RGB band from PCA ..."))
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
rr <- a1+a2
gg <- -a2
bb <- a3+a2-a1
rr <- (rr-min(rr)) / (max(rr)-min(rr)) * 255
gg <- (gg-min(gg)) / (max(gg)-min(gg)) * 255
bb <- (bb-min(bb)) / (max(bb)-min(bb)) * 255
Scld_CArgb = cbind(Scld_CA[,c("lon", "lat")], rr, gg, bb)
print(str(Scld_CArgb))

# convert to raster --------
print(paste(date(), "Convert Scld_CArgb file to raster rsScld_CArgb ..."))
# Note that the function raster::rasterFromXYZ does not work as expected in this case. 
Scld_CArgb1 <- Scld_CArgb # make a copy of Scld_CArgb
Scld_CArgb1 <- data.table::setDF(Scld_CArgb1) # convert to data frame to work best with sp library
sp::coordinates(Scld_CArgb1) <- ~ lon + lat
sp::proj4string(Scld_CArgb1) <- CRS(mycrs)
# coerce to SpatialPixelsDataFrame
sp::gridded(Scld_CArgb1) <- TRUE
# coerce to raster
rsScld_CArgb <- raster::stack(Scld_CArgb1)
print(rsScld_CArgb)
rm(Scld_CArgb1) # release space

# Plot the output RGB band (for quality control) --------
Scld_CA_plotRGBname=paste0(outdir,"plotRGB_predict_nogf_Scld_CAgrid_", today,".png")
Scld_CA_plotRGBtitle=paste0("plotRGB_predict_nogf_Scld_CAgrid_", today)
print(paste(date(), "Plotting RGB bands ...", Scld_CA_plotRGBname))
# Create an RGB image from the raster stack (quality control)
png(filename = Scld_CA_plotRGBname, height = 7, width = 7, units = "in", res = 150)
raster::plotRGB(rsScld_CArgb, r = 1, g = 2, b = 3, 
                ext = myextent, axes = TRUE, main = Scld_CA_plotRGBtitle)
dev.off()

# Plot the PC rotations and loading (for quality control) --------
Scld_CA_plotLOADname=paste0(outdir,"plotLOAD_predict_nogf_Scld_CAgrid_", today,".pdf")
Scld_CA_plotLOADtitle=paste0("plotLOAD_predict_nogf_Scld_CAgrid_", today)
print(paste(date(), "Plotting PCs rotations and loading ...", Scld_CA_plotLOADname))

set.seed(myseed)
id <- sample(x = 1:nrow(Scld_CArgb), size = 50000)
Scld_CArgbid = data.table::setDF(Scld_CArgb[id,]) # convert the sampled Scld_CArgb to data.frame to work best with the rgb() function

nvs <- dim(PCs$rotation)[1]
vec <- colnames(Scld_CA)[3:(3+7)] 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <-0.025
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1

pdf(file = Scld_CA_plotLOADname, height = 7, width = 7)
# plot the back ground sample location 
plot((PCs$x[id, 1:2]), xlim = xrng, ylim = yrng,
     pch = ".", cex = 4, asp = 1,
     col = rgb(Scld_CArgbid[,'rr'], Scld_CArgbid[,'gg'], Scld_CArgbid[,'bb'], max = 255),
     xlab = paste0("PC1 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC1']), ")"),
     ylab = paste0("PC2 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC2']), ")"),
     main = Scld_CA_plotLOADtitle
)
# now plot the 278 sites we analyzed
PCsites <- predict(PCs, Scld_site)
# combine them on the plots 
points(PCsites[, 1:2], col = "grey")
# add arrow vector 
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec,1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec,2]), 
     labels = vec)
dev.off()

# save the files ---------
Scld_CA_PCdata= paste0(outdir, "predict_nogf_Scld_CAgrid_prcomp.RData")
print(paste(date(), "Exporting Scld_CA Principal Component object (PCs) ...", Scld_CA_PCdata))
save(PCs, file = Scld_CA_PCdata)

Scld_CA_PCrgbdata=paste0(outdir, "predict_nogf_Scld_CAgrid_prcomp_RGB.RData")
print(paste(date(), "Exporting Scld_CArgb ...", Scld_CA_PCrgbdata))
save(Scld_CArgb, file = Scld_CA_PCrgbdata)

Scld_CA_PCrgbtif=paste0(outdir, "predict_nogf_Scld_CAgrid_prcomp_RGB.tif")
print(paste(date(), "Exporting rsScld_CArgb ...", Scld_CA_PCrgbtif))
raster::writeRaster(rsScld_CArgb, filename = Scld_CA_PCrgbtif, format = "GTiff")

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))

