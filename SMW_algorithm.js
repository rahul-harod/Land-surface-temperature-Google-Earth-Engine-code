// coefficients for the Statistical Mono-Window Algorithm
var SMWcoef = require('users/rharod4/LST_for_users:SMW_coefficients.js');

// Function to create a lookup between two columns in a 
// feature collection
var get_lookup_table = function(fc, prop_1, prop_2) {
  var reducer = ee.Reducer.toList().repeat(2);
  var lookup = fc.reduceColumns(reducer, [prop_1, prop_2]);
  return ee.List(lookup.get('list'));
};


exports.addBand = function(landsat){
  
  var wrap = function(image){
  
    // Select algorithm coefficients
    var coeff_SMW = ee.FeatureCollection(ee.Algorithms.If(landsat==='L4',SMWcoef.coeff_SMW_L4,
                                        ee.Algorithms.If(landsat==='L5',SMWcoef.coeff_SMW_L5,
                                        ee.Algorithms.If(landsat==='L7',SMWcoef.coeff_SMW_L7,
                                        SMWcoef.coeff_SMW_L8))));
    
    // Create lookups for the algorithm coefficients
    var A_lookup = get_lookup_table(coeff_SMW, 'TPWpos', 'A');
    var B_lookup = get_lookup_table(coeff_SMW, 'TPWpos', 'B');
    var C_lookup = get_lookup_table(coeff_SMW, 'TPWpos', 'C');
  
    // Map coefficients to the image using the TPW bin position
    var A_img = image.remap(A_lookup.get(0), A_lookup.get(1),0.0,'TPWpos').resample('bilinear');
    var B_img = image.remap(B_lookup.get(0), B_lookup.get(1),0.0,'TPWpos').resample('bilinear');
    var C_img = image.remap(C_lookup.get(0), C_lookup.get(1),0.0,'TPWpos').resample('bilinear');
    
    // select TIR band
    var tir = ee.String(ee.Algorithms.If(landsat==='L8','B10',
                        ee.Algorithms.If(landsat==='L7','B6_VCID_1',
                        'B6')));
    // compute the LST
    var lst = image.expression(
      'A*Tb1/em1 + B/em1 + C',
         {'A': A_img,
          'B': B_img,
          'C': C_img,
          'em1': image.select('Emissivity'),
          'Tb1': image.select(tir)
         }).updateMask((image.select('TPW').lt(0).not()).and(image.select('Emissivity').lte(0).not())) ;
         
    
    return image.addBands(lst.rename('LST'))
  };
  return wrap
}