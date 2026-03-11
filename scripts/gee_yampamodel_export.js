//YAMPA DAILY HYDROCLIMATE EXPORT FOR STREAMFLOW MODEL
//Use in Google Earth Engine code editor
//Exports a daily CSV with:
//  - date
//  - tmean_c
//  - prcp_mm
//  - srad_wm2
//  - swe_mm
//  - swe_b1_mm
//  - swe_b2_mm
//  - swe_b3_mm
//SWE bands are 3 equal-area elevation bands derived from DEM percentiles within the basin.
//Editable settings are at the top.
//speed is certainly not optimized fully, be patient with export. runs on my work internet take ~20 minutes


//SETTINGS

//Basin geometry: model is tuned for the Yampa at Deerlodge park, feel free to change at risk of not being applicable
var basinFC = ee.FeatureCollection('projects/jacksonparkerprofdev/assets/yampabasin');

//Start and end dates
var start  = ee.Date('2004-10-01');
var endReq = ee.Date('2024-09-30');   // These were dates used on original model tuning

//export name
var exportDescription = 'Yampa_GEE_model_inputs_2004-10-01_to_2024-09-30';

//Datasets
var DAYMET_ID = 'NASA/ORNL/DAYMET_V4';
var DEM_ID    = 'USGS/SRTMGL1_003';   // 30 m DEM



var basin = basinFC.union().geometry();

Map.centerObject(basin, 7);
Map.addLayer(basin, {}, 'Basin');


//setting end date to DAYMET availability

var dmAll = ee.ImageCollection(DAYMET_ID).filterBounds(basin);

var lastMillis    = ee.Number(dmAll.aggregate_max('system:time_start'));
var lastAvailDate = ee.Date(lastMillis);

print('Daymet last available date:', lastAvailDate);

var endMillis = ee.Number(endReq.millis()).min(lastMillis);
var endIncl   = ee.Date(endMillis);
var endExcl   = endIncl.advance(1, 'day');  

print('Requested start:', start);
print('Requested end (inclusive):', endReq);
print('Actual export end (inclusive):', endIncl);


//3 equal area elevation bands, calculated by pixel count

var dem = ee.Image(DEM_ID).select('elevation').clip(basin);

var elevPct = dem.reduceRegion({
  reducer: ee.Reducer.percentile([33.3333, 66.6667]),
  geometry: basin,
  scale: 30,
  maxPixels: 1e13,
  tileScale: 4,
  bestEffort: true
});

var p33 = ee.Number(elevPct.get('elevation_p33'));
var p67 = ee.Number(elevPct.get('elevation_p67'));

print('Elevation band break 1 (m):', p33);
print('Elevation band break 2 (m):', p67);

//masks for elevation bands
var band1Mask = dem.lte(p33);                          // low
var band2Mask = dem.gt(p33).and(dem.lte(p67));        // mid
var band3Mask = dem.gt(p67);                          // high

Map.addLayer(dem, {min: 1500, max: 4000}, 'DEM', false);
Map.addLayer(band1Mask.selfMask(), {palette: ['#2b83ba']}, 'Band 1 low', false);
Map.addLayer(band2Mask.selfMask(), {palette: ['#abdda4']}, 'Band 2 mid', false);
Map.addLayer(band3Mask.selfMask(), {palette: ['#d7191c']}, 'Band 3 high', false);



//Daymet collection
//  tmean_c   -> for PDD / rain partition in R
//  prcp_mm   -> precipitation
//  srad_wm2  -> shortwave radiation
//  swe_mm    -> basin mean SWE
//  swe_b1-3  -> SWE in 3 elevation bands

var daymet = ee.ImageCollection(DAYMET_ID)
  .filterBounds(basin)
  .filterDate(start, endExcl)
  .select(['tmax', 'tmin', 'prcp', 'swe', 'srad']);

print('Daymet image count after filters:', daymet.size());


//map to table

var dailyTable = daymet.map(function(img) {
  
  //Daymet units:
  //  tmax, tmin  -> °C
  //  prcp        -> mm/day
  //  swe         -> kg/m^2 ~ mm water equivalent
  //  srad        -> W/m^2
  
  var tmeanC = img.select('tmax')
    .add(img.select('tmin'))
    .divide(2)
    .rename('tmean_c');
  
  var prcpMM = img.select('prcp').rename('prcp_mm');
  var sweMM  = img.select('swe').rename('swe_mm');
  var srad   = img.select('srad').rename('srad_wm2');
  
  var stack = tmeanC.addBands([prcpMM, srad, sweMM]);
  var proj  = img.select('swe').projection();
  var scale = proj.nominalScale();
  
  //basin-wide means
  var basinStats = stack.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: basin,
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  });
  
  //SWE by elevation band
  var sweB1 = sweMM.updateMask(band1Mask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: basin,
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  }).get('swe_mm');
  
  var sweB2 = sweMM.updateMask(band2Mask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: basin,
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  }).get('swe_mm');
  
  var sweB3 = sweMM.updateMask(band3Mask).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: basin,
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4,
    bestEffort: true
  }).get('swe_mm');
  
  return ee.Feature(null, {
    date: img.date().format('YYYY-MM-dd'),
    tmean_c: basinStats.get('tmean_c'),
    prcp_mm: basinStats.get('prcp_mm'),
    srad_wm2: basinStats.get('srad_wm2'),
    swe_mm: basinStats.get('swe_mm'),
    swe_b1_mm: sweB1,
    swe_b2_mm: sweB2,
    swe_b3_mm: sweB3
  });
});

print('Preview:', dailyTable.limit(10));

//export as CSV

Export.table.toDrive({
  collection: ee.FeatureCollection(dailyTable),
  description: exportDescription,
  fileFormat: 'CSV',
  selectors: [
    'date',
    'tmean_c',
    'prcp_mm',
    'srad_wm2',
    'swe_mm',
    'swe_b1_mm',
    'swe_b2_mm',
    'swe_b3_mm'
  ]
});
