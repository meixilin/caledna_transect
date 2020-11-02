# Title: calculate distance matrix between sites --------
# Author: Meixi Lin
# Date: Mon Aug 27 16:16:25 2018

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(geosphere)
library(sp)
library(ggmap)
library(rgeos)
date() # the execution date

# read data ---------
metadata <- read.csv("./final_data/metadata/Final_metadata_05312019.csv", stringsAsFactors = F)
latlong <- cbind(metadata$Longitude, metadata$Latitude)
rownames(latlong) <- metadata$MatchName
# define points
distmatrix <- matrix(data = NA, nrow = nrow(latlong), ncol = nrow(latlong), dimnames = list(row.names(latlong), row.names(latlong)))

# calculate distance --------
for (ii in 1:nrow(distmatrix)) {
    for (jj in 1:ncol(distmatrix)) {
        distmatrix[ii,jj] <- geosphere::distGeo(latlong[ii,], latlong[jj,])
    }
}

# here is the value of the distance matrix 
# image(t(distmatrix))

save(distmatrix, file = "./derive_data/step0_prepare_data/dist_caledna_metadata.RData")

# using hierachical clustering to cluster any localities < 1000 m apart --------
hc <- hclust(as.dist(distmatrix), method="complete")
# add to distance group 
dd <- 1000
clust <- as.data.frame(cutree(hc, h = dd)); colnames(clust) <- "clust"
latlong <- merge(latlong, clust, by = "row.names")
colnames(latlong)[2:4] <- c("Longitude", "Latitude", "cluster")
# Test if it was the same as the final metadata 
table(clust == metadata$clust)
latlong <- latlong[,-1]

# get each cluster's centroid --------
# get the centroid coords for each cluster
pts <- sp::SpatialPointsDataFrame(coords = latlong[,1:2], data = clust, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) # code for WGS1984
cent <- matrix(ncol=2, nrow=max(latlong$cluster))
for (ii in 1:max(latlong$cluster)) {
    # gCentroid from the rgeos package
    cent[ii,] <- gCentroid(subset(pts, clust == ii))@coords
}
colnames(cent) <- c("Longitude", "Latitude")

# plot the result on a basemap --------
cabasemap <- get_map(location = "california", zoom = 6)

clustmap <- ggmap(cabasemap) + 
    geom_point(aes(x = Longitude, y = Latitude, colour = as.factor(cluster)), data = latlong) +
    scale_color_hue() + 
    geom_point(aes(x = Longitude, y = Latitude), data = as.data.frame(cent), shape = 1, size = 5) 
                 
ggsave("clustmap.pdf", plot = clustmap, path = "./plots/", height = 20, width = 20, device = "pdf")

# the number of clusters 
max(clust)

# # write the cluster data out --------
# check if they are in the same order 
table(row.names(clust) == metadata$MatchName)
metadata <- cbind(metadata, clust)
write.csv(metadata, file = "./final_data/Final_metadata.csv", row.names = F, quote = F)
    
# code adapted from: 
# https://gis.stackexchange.com/questions/17638/how-to-cluster-spatial-data-in-r
    


