// MODULES DECLARATION -----------------------------------------------------------
// Total Precipitable Water 
var NCEP_TPW = require('users/rharod4/LST_for_users:NCEP_NEW.js')
//cloud mask
var cloudmask = require('users/rharod4/LST_for_users:Cloud_Mask.js')
//surface emissivity
var EM = require('users/rharod4/LST_for_users:Emissivity.js')
// land surface temperature
var LST = require('users/rharod4/LST_for_users:SMW_algorithm.js')
// --------------------------------------------------------------------------------


var COLLECTION = ee.Dictionary({
  'L4': {
    'TOA': ee.ImageCollection('LANDSAT/LT04/C01/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LT04/C01/T1_SR'),
    'TIR': ['B6',]
  },
  'L5': {
    'TOA': ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LT05/C01/T1_SR'),
    'TIR': ['B6',]
  },
  'L7': {
    'TOA': ee.ImageCollection('LANDSAT/LE07/C01/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LE07/C01/T1_SR'),
    'TIR': ['B6_VCID_1','B6_VCID_2'],
  },
  'L8': {
    'TOA': ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LC08/C01/T1_SR'),
    'TIR': ['B10','B11']
  }
});


exports.collection = function(landsat, date_start, date_end, geometry, USE_URBAN) {

  // load TOA Radiance/Reflectance
  var collection_dict = ee.Dictionary(COLLECTION.get(landsat));

  var landsatTOA = ee.ImageCollection(collection_dict.get('TOA'))
                .filter(ee.Filter.date(date_start, date_end))
                .filterBounds(geometry)
                .map(cloudmask.toa);
        
        
  var count = landsatTOA.size();
  var availability=ee.Algorithms.If(count.eq(0),'Image Not available in given date range','');
  print(availability);
  
  //LULC MAP
  var year=date_start.split('-')[0];
  var d1=year+'-01-01';
  var dataset = ee.ImageCollection('MODIS/006/MCD12Q1').filterDate(d1)
  var LULC = dataset.select('LC_Type1').first();
  
  // load Surface Reflectance collection for NDVI
  var landsatSR = ee.ImageCollection(collection_dict.get('SR'))
                .filter(ee.Filter.date(date_start, date_end))
                .filterBounds(geometry)
                .map(cloudmask.sr)
                .map(NCEP_TPW.addBand)
                .map(EM.emissivity(landsat,USE_URBAN,LULC));
  // select TIR bands
  var tir = ee.List(collection_dict.get('TIR'));
  var landsatALL = (landsatSR.combine(landsatTOA.select(tir), true));
  
  // compute the LST
  var landsatLST = landsatALL.map(LST.addBand(landsat));

  return landsatLST;
};