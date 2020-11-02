// This is a script for downloading land cover dataset 
// Author: Meixi Lin 

/* first import the geometry for downloading */
var geometry = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-126.22539062499999, 30.887798884911486],
          [-111.63554687499999, 30.887798884911486],
          [-112.07499999999999, 43.24896238414401],
          [-126.31328124999999, 43.24896238414401]]]);

/* define a palette for mapping tree cover */
var myParams = {min: 0, max: 100, palette: ['blue', 'white', 'green']};

/***************************************************
 * get ca boundary 
 ***************************************************/
var tiger = ee.FeatureCollection('TIGER/2016/States');
var cabound = tiger
            .filterMetadata('NAME', 'equals', 'California')
            .first()
            .geometry();
print(cabound);
print(cabound.projection());
Map.addLayer(cabound, {color:'ff0000'});

/***************************************************
 * import nlcd dataset 
 ***************************************************/
var dataset = ee.Image('USGS/NLCD/NLCD2011')
              .clip(cabound);

var dataset_wgs84 = dataset.reproject({crs: 'EPSG:4326', scale: 100});
print(dataset);
print(dataset.projection());
print(dataset_wgs84);
print(dataset_wgs84.projection());

// add layer for percent tree cover and error rate 
Map.addLayer(dataset.select('percent_tree_cover'), myParams, 'percent_tree');
Map.addLayer(dataset_wgs84.select('percent_tree_cover'), myParams, 'percent_tree_wgs84');

Export.image.toDrive({image: dataset_wgs84.select('landcover'), description: 'EE_NLCD_CA_wgs84_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: dataset_wgs84.select('impervious'), description: 'EE_imprv_CA_wgs84_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: dataset_wgs84.select('percent_tree_cover'), description: 'EE_ptrcv_CA_wgs84_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: dataset_wgs84.select('percent_tree_error'), description: 'EE_ptrer_CA_wgs84_100m', scale: 100, region: geometry, maxPixels: 1e10});

