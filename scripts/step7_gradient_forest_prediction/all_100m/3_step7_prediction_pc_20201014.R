# Title: perform principal component (PC) calculation for the predict_gf_Trns_CAgrid --------
# Author: Meixi Lin
# Date: Sun Nov  3 16:38:59 2019
# Author: Meixi Lin
# Date: Thu Oct 15 23:28:18 2020
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

# get arguments ---------
args <-  commandArgs(trailingOnly=TRUE)

extrap <- as.logical(args[1]) # extrapolation 
Trns_site_filepath <- as.character(args[2]) # Path to file "predict_gf_Trns_site_extrap.RData/predict_gf_Trns_site_noextrap.RData" Transformed sample sites
Trns_CA_filepath <- as.character(args[3]) # Path to file "predict_gf_Trns_CAgrid_extrap.RData/predict_gf_Trns_CAgrid_noextrap.RData" Tranformed CAgrid in 100m x 100m
outdir <- as.character(args[4]) # Output path 

dir.create(outdir, recursive = FALSE, showWarnings = FALSE)

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84" 
# Extent of CA
myextent <- raster::extent(CAlimit)
myseed = 17

# define a function --------
getTrns_CA_names <- function(extrap, outdir) {
    today = format(Sys.Date(), "%Y%m%d")
    if (extrap == TRUE) {
        Trns_CA_plotRGBname = paste0(outdir, "plotRGB_predict_gf_Trns_CAgrid_extrap_", today,".png")
        Trns_CA_plotRGBtitle = paste0("plotRGB_predict_gf_Trns_CAgrid_extrap_", today)
        Trns_CA_plotLOADname = paste0(outdir, "plotLOAD_predict_gf_Trns_CAgrid_extrap_", today,".pdf")
        Trns_CA_plotLOADtitle = paste0("plotLOAD_predict_gf_Trns_CAgrid_extrap_", today)
        Trns_CA_PCdata = paste0(outdir, "predict_gf_Trns_CAgrid_extrap_prcomp.RData")
        Trns_CA_PCrgbdata = paste0(outdir, "predict_gf_Trns_CAgrid_extrap_prcomp_RGB.RData")
        Trns_CA_PCrgbtif = paste0(outdir, "predict_gf_Trns_CAgrid_extrap_prcomp_RGB.tif")
    } else { 
        if (extrap == FALSE) {
            Trns_CA_plotRGBname = paste0(outdir, "plotRGB_predict_gf_Trns_CAgrid_noextrap_", today,".png")
            Trns_CA_plotRGBtitle = paste0("plotRGB_predict_gf_Trns_CAgrid_noextrap_", today)
            Trns_CA_plotLOADname = paste0(outdir, "plotLOAD_predict_gf_Trns_CAgrid_noextrap_", today,".pdf")
            Trns_CA_plotLOADtitle = paste0("plotLOAD_predict_gf_Trns_CAgrid_noextrap_", today)
            Trns_CA_PCdata = paste0(outdir, "predict_gf_Trns_CAgrid_noextrap_prcomp.RData")
            Trns_CA_PCrgbdata = paste0(outdir, "predict_gf_Trns_CAgrid_noextrap_prcomp_RGB.RData")
            Trns_CA_PCrgbtif = paste0(outdir, "predict_gf_Trns_CAgrid_noextrap_prcomp_RGB.tif")
        } else {
            stop("Wrong extrap value!")
        }
    } 
    output = list(Trns_CA_plotRGBname,Trns_CA_plotRGBtitle,Trns_CA_plotLOADname,Trns_CA_plotLOADtitle,Trns_CA_PCdata,Trns_CA_PCrgbdata,Trns_CA_PCrgbtif)
    names(output) = c("Trns_CA_plotRGBname","Trns_CA_plotRGBtitle","Trns_CA_plotLOADname","Trns_CA_plotLOADtitle","Trns_CA_PCdata","Trns_CA_PCrgbdata","Trns_CA_PCrgbtif")
    return(output)
}

# load data --------
print(paste(date(), "Loading ...", Trns_CA_filepath))

load(Trns_CA_filepath) # should be named as Trns_CA
print(str(Trns_CA))
Trns_site=loadRData(Trns_site_filepath) 
print(str(Trns_site))

# get PC and convert to RGB --------
print(paste(date(), "Principal Component Analyses ..."))

# run PCA on all the variables available (the first two cols are 'lon' 'lat')
PCs = prcomp(Trns_CA[,3:ncol(Trns_CA)])
print(summary(PCs))

print(paste(date(), "RGB band from PCA ..."))
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]
rr <- a1+a2
gg <- -a2
bb <- a3+a2-a1
rr <- (rr-min(rr)) / (max(rr)-min(rr)) * 255 # maximize differences 
gg <- (gg-min(gg)) / (max(gg)-min(gg)) * 255
bb <- (bb-min(bb)) / (max(bb)-min(bb)) * 255
Trns_CArgb = cbind(Trns_CA[,c("lon", "lat")], rr, gg, bb)

# convert to raster --------
print(paste(date(), "Convert Trns_CArgb file to raster rsTrns_CArgb ..."))
# Note that the function raster::rasterFromXYZ does not work as expected in this case. 
Trns_CArgb1 <- Trns_CArgb # make a copy of Trns_CArgb
Trns_CArgb1 <- data.table::setDF(Trns_CArgb1) # convert to data frame to work best with sp library
sp::coordinates(Trns_CArgb1) <- ~ lon + lat
sp::proj4string(Trns_CArgb1) <- CRS(mycrs)
# coerce to SpatialPixelsDataFrame
sp::gridded(Trns_CArgb1) <- TRUE
# coerce to raster
rsTrns_CArgb <- raster::stack(Trns_CArgb1)
print(rsTrns_CArgb)
rm(Trns_CArgb1) # release space

# Plot the output RGB band (for quality control) --------
Trns_CA_plotRGBname=getTrns_CA_names(extrap, outdir)[['Trns_CA_plotRGBname']]
Trns_CA_plotRGBtitle=getTrns_CA_names(extrap, outdir)[['Trns_CA_plotRGBtitle']]
print(paste(date(), "Plotting RGB bands ...", Trns_CA_plotRGBname))
# Create an RGB image from the raster stack (quality control)
png(filename = Trns_CA_plotRGBname, height = 7, width = 7, units = "in", res = 150)
raster::plotRGB(rsTrns_CArgb, r = 1, g = 2, b = 3, 
                ext = myextent, axes = TRUE, main = Trns_CA_plotRGBtitle)
dev.off()

# Plot the PC rotations and loading (for quality control) --------
Trns_CA_plotLOADname=getTrns_CA_names(extrap, outdir)[['Trns_CA_plotLOADname']]
Trns_CA_plotLOADtitle=getTrns_CA_names(extrap, outdir)[['Trns_CA_plotLOADtitle']]
print(paste(date(), "Plotting PCs rotations and loading ...", Trns_CA_plotLOADname))

set.seed(myseed)
id <- sample(x = 1:nrow(Trns_CArgb), size = 50000)
Trns_CArgbid = data.table::setDF(Trns_CArgb[id,]) # convert the sampled Trns_CArgb to data.frame to work best with the rgb() function

nvs <- dim(PCs$rotation)[1]
vec <- colnames(Trns_CA)[3:(3+7)] 
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <-25
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1

pdf(file = Trns_CA_plotLOADname, height = 7, width = 7)
# plot the back ground sample location 
plot((PCs$x[id, 1:2]), xlim = xrng, ylim = yrng,
     pch = ".", cex = 4, asp = 1,
     col = rgb(Trns_CArgbid[,'rr'], Trns_CArgbid[,'gg'], Trns_CArgbid[,'bb'], max = 255),
     xlab = paste0("PC1 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC1']), ")"),
     ylab = paste0("PC2 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC2']), ")"),
     main = Trns_CA_plotLOADtitle
)
# now plot the 278 sites we analyzed
PCsites <- predict(PCs, Trns_site)
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
Trns_CA_PCdata=getTrns_CA_names(extrap, outdir)[['Trns_CA_PCdata']]
print(paste(date(), "Exporting Trns_CA Principal Component object (PCs) ...", Trns_CA_PCdata))
save(PCs, file = Trns_CA_PCdata)

Trns_CA_PCrgbdata=getTrns_CA_names(extrap, outdir)[['Trns_CA_PCrgbdata']]
print(paste(date(), "Exporting Trns_CArgb ...", Trns_CA_PCrgbdata))
save(Trns_CArgb, file = Trns_CA_PCrgbdata)

Trns_CA_PCrgbtif=getTrns_CA_names(extrap, outdir)[['Trns_CA_PCrgbtif']]
print(paste(date(), "Exporting rsTrns_CArgb ...", Trns_CA_PCrgbtif))
raster::writeRaster(rsTrns_CArgb, filename = Trns_CA_PCrgbtif, format = "GTiff")

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))
closeAllConnections()