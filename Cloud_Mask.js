// cloudmask for TOA data
exports.toa = function(image) {
  var qa = image.select('BQA');
  var mask = qa.bitwiseAnd(1 << 4).eq(0);
  return image.updateMask(mask);
};

// cloudmask for SR data
exports.sr = function(image) {
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(1 << 3)
    .or(qa.bitwiseAnd(1 << 5))
  return image.updateMask(mask.not());
};