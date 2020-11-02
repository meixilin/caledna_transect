# Title: create a sample map --------
# Author: Meixi Lin
# Date: Tue Sep 10 10:57:04 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

library(dplyr)
library(ggplot2) # 
library(ggmap)
library(rgdal)
library(ggsn)

# souurce function
source("./scripts/function_transect.R")


# note that you need to enable billing account for it to work
# register a google api key and enable: 
# 1. Maps Static API
# 2. Geocoding API
ggmap::register_google(key = "YOUR_API_KEY")

date() # the execution date


# load data --------
ucnrs <- readOGR(dsn = "./maps_vectors_rasters_local/ucnrs/shapefiles/NRS_Boundaries.shp") 
ucnrs_df <- ggplot2::fortify(ucnrs)
raster::crs(ucnrs)
calim <- readOGR(dsn = "./maps_vectors_rasters_local/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
raster::crs(calim) 
calim_df <- ggplot2::fortify(calim)

# load data frame --------
# load ucnrs sites 
sites.nrs <- read.csv("./final_data/ucnrs/sites_ucnrs_final.csv", row.names = 1, stringsAsFactors = F)
sites.all <- loadRData(file = "./final_data/Final_metadata.RData")

# # convert to plotting frames
# forplot <- sites.all %>%
#     dplyr::select(MatchName, Longitude, Latitude, transect, clust) %>%
#     dplyr::mutate(innrs = MatchName %in% sites.nrs$MatchName) %>%
#     dplyr::group_by(innrs, clust) %>%
#     dplyr::summarise(long = mean(Longitude), 
#                      lat = mean(Latitude), 
#                      transect = first(transect),
#                      count = n())

forplot <- sites.all %>%
        dplyr::select(MatchName, Longitude, Latitude, transect, loc) %>%
        dplyr::mutate(innrs = loc %in% sites.nrs$loc) %>%
        dplyr::group_by(loc, transect, innrs) %>%
        dplyr::summarise(long = mean(Longitude),
                         lat = mean(Latitude),
                         count = n())
colnames(forplot)
forplot$transect <- factor(forplot$transect, labels = c("Coast", "Forest", "Shrub"))

# start plotting --------
fac.labels <- c(Coastal = "Coast", Forest = "Forest", ShrubScrub = "Shrub")
cabase <- get_map(location = "california", maptype = "roadmap", zoom = 6)
pp <- ggmap(cabase) + 
    geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "snow4", fill = 'snow4', alpha = 0.5) +
    geom_polygon(data = ucnrs_df, aes(x = long, y = lat, group = group), colour = "khaki", fill = 'khaki', alpha = 0.7, size = 1.3) +
    geom_point(data = forplot, aes(x = long, y = lat, colour = transect, fill = transect, size = count, shape = innrs), alpha = 0.8) +
    scale_shape_manual(name = "In UCNRS", values = c(24,21)) +
    scale_colour_discrete(name = "Transect", l = 20) +
    scale_size(name = "Sites count") +
    facet_grid(. ~ transect, labeller = labeller(transect = fac.labels)) +
    # scale_size_area(max_size = 5) +
    # ggsn::north(data = calim_df, location = "topleft") +
    ggsn::scalebar(data = calim_df,
                   dist = 150, dist_unit = "km", border.size = 0.5, st.size = 2.5, st.dist = 0.05,
                   transform = TRUE, model = "WGS84", location = "bottomleft") +
    labs(x = "Longitude", y = "Latitude") +
    guides(fill = "none") +
    theme_bw() 


ggsave(filename = "./plots_important/step0_prepare_data/sample_sites_leg.pdf", plot = pp, device = "pdf", width = 10, height = 4)

pp1 <- pp + 
    theme(legend.position = "none") 
ggsave(filename = "./plots_important/step0_prepare_data/sample_sites.pdf", plot = pp1, device = "pdf", width = 9, height = 4)

# to get a legend for background
dummy <- data.frame(matrix(nrow = 2, ncol = 2, c("CA Boundary", "UCNRS Boundary", 2, 1)))
pp2 <- ggplot(data = dummy, aes(x = X1, y = X2, fill = X1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("snow4", "khaki")) +
    theme_bw()

ggsave(filename = "./plots_important/step0_prepare_data/sample_sites_bb_leg.pdf", plot = pp2, device = "pdf", width = 4, height = 4)
