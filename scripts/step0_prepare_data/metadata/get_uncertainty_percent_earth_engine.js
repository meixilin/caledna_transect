// This is a script for generating rasters for uncertainty layers in percentage standard deviation
// NOT USED
// Author: Meixi Lin
// Date: 2020/09/30

/***************************************************
 * preparations
 ***************************************************/
var s2Params = {min: -1, max: 1, palette: ['blue', 'white', 'green']};

/***************************************************
 * get ca boundary 
 ***************************************************/
var tiger = ee.FeatureCollection('TIGER/2016/States');
var cabound = tiger
            .filterMetadata('NAME', 'equals', 'California')
            .first()
            .geometry();
print(cabound);
Map.addLayer(cabound, {color:'ff0000'});

/* first import the geometry for downloading */
var exportgeometry = ee.Geometry.Polygon(
        [[[-126.22539062499999, 30.887798884911486],
          [-111.63554687499999, 30.887798884911486],
          [-112.07499999999999, 43.24896238414401],
          [-126.31328124999999, 43.24896238414401]]]);

/***************************************************
 * import sentinel dataset 
 ***************************************************/

/**
 * Function to mask clouds using the Sentinel-2 QA band
 * @param {ee.Image} image Sentinel-2 image
 * @return {ee.Image} cloud masked Sentinel-2 image
 */
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}

// Map the function over the time of sampling and take the median.
// Load Sentinel-2 TOA reflectance data.
var cloudMasked_sentinel = ee.ImageCollection('COPERNICUS/S2')
                  .filterDate('2017-01-01','2018-07-01')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds);
// print(cloudMasked_sentinel);

var image_sentinel_stddev = cloudMasked_sentinel
                            .reduce(ee.Reducer.stdDev())
                            .clip(cabound);
                            
var image_sentinel_stddev_wgs84 = image_sentinel_stddev.reproject({crs: 'EPSG:4326', scale: 100});
var image_sentinel_mean = cloudMasked_sentinel
                            .reduce(ee.Reducer.mean())
                            .clip(cabound);
var image_sentinel_mean_wgs84 = image_sentinel_mean.reproject({crs: 'EPSG:4326', scale: 100});

// Get the percent variation = STDDEV/MEAN 
var image_sentinel_perstddev_wgs84_B1 = image_sentinel_stddev_wgs84.select('B1').divide(landsat1999.select('B4')

print(image_sentinel_stddev); // already WGS84 
print(image_sentinel_stddev_wgs84);

/***************************************************
* Export the images to google drive

// Export.image.toDrive({image: images, description: 'EE_layer_CA_100m', scale: 100, maxPixels: 1e10,
//                       region: geometry});

Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B1_stdDev'), description: 'EE_B1_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B2_stdDev'), description: 'EE_B2_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B3_stdDev'), description: 'EE_B3_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B4_stdDev'), description: 'EE_B4_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B5_stdDev'), description: 'EE_B5_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B6_stdDev'), description: 'EE_B6_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B7_stdDev'), description: 'EE_B7_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B8_stdDev'), description: 'EE_B8_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B8A_stdDev'), description: 'EE_B8A_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B9_stdDev'), description: 'EE_B9_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B10_stdDev'), description: 'EE_B10_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B11_stdDev'), description: 'EE_B11_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('B12_stdDev'), description: 'EE_B12_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('QA10_stdDev'), description: 'EE_QA10_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('QA20_stdDev'), description: 'EE_QA20_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
Export.image.toDrive({image: image_sentinel_stddev_wgs84.select('QA60_stdDev'), description: 'EE_QA60_stdDev_CA_wgs84_100m', scale: 100, region: exportgeometry, maxPixels: 1e10});
***************************************************/