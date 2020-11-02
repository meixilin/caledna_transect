# Title: function def for generate metadata --------
# Author: Meixi Lin
# Date: Thu Mar 14 19:51:55 2019
# Author:
# Date:
# Modification:

# define a function that return napt from a base raster --------
# napt: a spatial point object 
# traster.o: original raster to extract value from 
# mywidth: the size in degree to crop a raster to extract value 
na.near.neighbor <- function(napt, traster.o, mywidth) {
    nabuffer = rgeos::gBuffer(spgeom = napt, width = mywidth)
    traster = crop(traster.o, nabuffer)
    nadistance <- distanceFromPoints(object = traster, xy = napt) 
    names(nadistance) <- "dist_na_near"
    ss <- stack(nadistance, traster)
    ss <- mask(ss, traster)
    ssma <- as.matrix(ss)
    mindis <- which(ssma[,'dist_na_near'] 
          == min(ssma[,'dist_na_near'], na.rm = T))[1]
    # if mindis not working 
    if (length(mindis) == 0) {
        print("something wrong in get nearest na value, break out.")
        output <- rep(NA, times = 4)
    } else {
        # extract values
        value <- extract(x = ss, y = mindis)
        # extract the moved coordinates 
        coor <- xyFromCell(object = ss, cell = mindis)
        # start writing the navalues dataframe 
        # data frame for missing points: Long_corrected; Lat_corrected; matrix index; value; distance 
        output <- c(coor, value)
    }
    return(output)
}

# a function to locate na values redo extraction and return a data frame with all info --------
# myrr: my raster for extraction
# pts: the points used
# myex: the extraction result (possibly with NAs)
fix.napts <- function(myrr, pts, myex, mywidth, maxdist = NA) {
    naid = which(is.na(myex[,'raster_value']))
    napts <- pts[naid]
    navalues <- as.data.frame(matrix(nrow = length(napts), ncol = 7))
    colnames(navalues) <- c("MatchName", "Longitude_c", "Latitude_c", "raster_distance", "raw_raster_value", "naid", "raster_value")
    
    for (ii in 1:length(napts)) {
        xx = napts[ii]
        navalues[ii,1] <- row.names(xx@coords)
        navalues[ii,2:5] <- na.near.neighbor(xx, myrr, mywidth)
        navalues[ii,6] = naid[ii]
        # cut out the values with too much distance 
        if (!is.na(maxdist)) {
            if (navalues[ii, 'raster_distance'] > maxdist) {
                navalues[ii, 'raster_value'] <- NA 
            } else {
                navalues[ii, 'raster_value'] <- navalues[ii, 'raw_raster_value']
            }
        } else {
            # do nothing 
            navalues[ii, 'raster_value'] <- navalues[ii, 'raw_raster_value']
        }
    }

    output <- napts@coords %>% as.data.frame() %>% 
        tibble::rownames_to_column(var = "MatchName") %>%
        dplyr::inner_join(navalues, by = "MatchName") %>%
        dplyr::arrange(naid)
    return(output)
} 
