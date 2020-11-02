// This is a script for generating rasters for predictions 
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
var image_sentinel = cloudMasked_sentinel.median().clip(cabound);
print(image_sentinel);

/***************************************************
 * caculate NDVI value 
 ***************************************************/
var addNDVI = function(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands(ndvi);
};
 
var withNDVI = cloudMasked_sentinel.map(addNDVI);
 
// Make a "greenest" pixel composite.
var ndvis2 = withNDVI.qualityMosaic('NDVI').clip(cabound);
print(ndvis2);

Map.addLayer(ndvis2.select('NDVI'), 
            {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 
            'ndvi_sentinel');

/***************************************************
 * import LANDSAT derived values 
 ***************************************************/
// .filterDate('2017-01-01','2018-07-01')
// note that for all LANDSAT derived dataset, I used a shorter time frame to better capitulate the time frame where sampling occurred 
var ndvi32 = ee.ImageCollection('LANDSAT/LC08/C01/T1_32DAY_NDVI')
                  .filterDate('2017-04-01', '2017-06-15')
                  .median()
                  .clip(cabound);

Map.addLayer(ndvi32, 
            {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 
            'ndvi_landsat');

var evi = ee.ImageCollection("LANDSAT/LC08/C01/T1_32DAY_EVI")
                  .filterDate('2017-04-01', '2017-06-15')
                  .median()
                  .clip(cabound);
                  
Map.addLayer(evi, 
            {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 
            'evi_landsat');

var nbrt = ee.ImageCollection("LANDSAT/LC08/C01/T1_32DAY_NBRT")
                  .filterDate('2017-04-01', '2017-06-15')
                  .median()
                  .clip(cabound);

Map.addLayer(nbrt, 
            {min: -1, max: 1, palette: ['blue', 'white', 'green']}, 
            'nbrt_landsat');

var greenness = ee.Image('LANDSAT/LC8_L1T_ANNUAL_GREENEST_TOA/2017')
                .clip(cabound);
print(greenness);
Map.addLayer(greenness.select('greenness'));

/***************************************************
* Get all images 
***************************************************/
var images = image_sentinel
.addBands(ndvis2.select(['NDVI'],['NDVIS2']))
.addBands(ndvi32.select(['NDVI'],['NDVI32']))
.addBands(evi.select('EVI'))
.addBands(nbrt.select('NBRT'))
.addBands(greenness.select('greenness'))
.toFloat();

print(images);

/***************************************************
* Export the images to google drive
***************************************************/
// Export.image.toDrive({image: images, description: 'EE_layer_CA_100m', scale: 100, maxPixels: 1e10,
//                       region: geometry});

Export.image.toDrive({image: images.select('B1'), description: 'EE_B1_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B2'), description: 'EE_B2_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B3'), description: 'EE_B3_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B4'), description: 'EE_B4_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B5'), description: 'EE_B5_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B6'), description: 'EE_B6_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B7'), description: 'EE_B7_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B8'), description: 'EE_B8_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B8A'), description: 'EE_B8A_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B9'), description: 'EE_B9_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B10'), description: 'EE_B10_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B11'), description: 'EE_B11_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('B12'), description: 'EE_B12_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('QA10'), description: 'EE_QA10_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('QA20'), description: 'EE_QA20_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('QA60'), description: 'EE_QA60_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('NDVIS2'), description: 'EE_NDVIS2_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('NDVI32'), description: 'EE_NDVI32_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('greenness'), description: 'EE_greenness_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('EVI'), description: 'EE_EVI_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});
Export.image.toDrive({image: images.select('NBRT'), description: 'EE_NBRT_CA_100m', scale: 100, region: geometry, maxPixels: 1e10});














