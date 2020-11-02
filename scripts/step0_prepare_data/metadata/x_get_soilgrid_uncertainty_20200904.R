# Title: Get soilgrid layers
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Sep  4 17:24:18 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(leaflet)
library(leaflet.extras)
library(mapview)
library(rgdal)

source("./scripts/function_transect.R")
today = format(Sys.Date(), "%Y%m%d")

# def functions --------
get_soilgrid_var <- function(var) {
    soilgrid_var <- paste0("https://maps.isric.org/mapserv?map=/map/",var, ".map")
    pp <- leaflet(ca) %>% 
        setView(lng = mean(CAlimit[,'Longitude']), lat = mean(CAlimit[,'Latitude']), zoom = 6) %>% 
        addPolylines() %>%
        addWMSTiles(
            soilgrid_var,
            layers = paste0(var,"_0-5cm_uncertainty"),
            options = WMSTileOptions(format = "image/png", transparent = TRUE)
        ) %>% 
        addWMSLegend(uri = paste0(soilgrid_var, 
                                  "&version=1.3.0&service=WMS&request=GetLegendGraphic&sld_version=1.1.0&layer=",var,"_0-5cm_uncertainty&format=image/png&STYLE=default&"))
    
    mapshot(pp, file = paste0("./plots/step0_prepare_data/",var,"_uncertainty_",today,".pdf"),remove_controls = c("zoomControl"))
    return(0)
}
# def variables --------
mycrs = "+proj=longlat +datum=WGS84"
varlist = c("cec", "bdod", "clay", "nitrogen", "sand", "silt", "soc")


# load data --------
# get ca boundary 
ca0 = readOGR("./maps_vectors_rasters/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
ca = spTransform(ca0, CRS(mycrs)); rm(ca0)

# main --------
lapply(varlist, get_soilgrid_var)

# cleanup --------
