/**
 * 02b_export_lulc_maps.js
 *
 * Exports the classified maps and Random Forest class probabilities produced by
 * 02_classify_and_assess.js as one GeoTIFF per reserve, for the data archive
 * and for the R calibration of the error-propagation draws
 * (R/calibration/calibrate_correlation.R).
 *
 * Self-contained (no require), so it can be pasted straight into the Code Editor.
 *
 * Output  Drive  lulc_<reserve>.tif   local UTM, 10 m, Byte, nodata 255
 *   bands  class_2019, p_no_2019, p_low_2019, p_high_2019, n_obs_2019,
 *          class_2024, p_no_2024, p_low_2024, p_high_2024, n_obs_2024
 *          class: 0 no vegetation, 1 low vegetation, 2 high vegetation
 *          p_*:   share of RF trees voting for the class (%)
 *          n_obs: clear Sentinel-2 observations in the composite (capped at 254)
 */

var ROOT = 'projects/YOUR_CLOUD_PROJECT/assets/nigeria_fragmentation';   // same as ASSET_ROOT in config.js
var FOLDER = 'nigeria_fragmentation_gee';
var YEARS = [2019, 2024];
var RESERVES = [
  {id: 'Adoru',        utm: 'EPSG:32632'},
  {id: 'Akanto',       utm: 'EPSG:32632'},
  {id: 'Gambari',      utm: 'EPSG:32631'},
  {id: 'Maradun',      utm: 'EPSG:32632'},
  {id: 'Mayo_Ndaga',   utm: 'EPSG:32632'},
  {id: 'Stubbs_Creek', utm: 'EPSG:32632'},
  {id: 'Upper_Orashi', utm: 'EPSG:32632'},
  {id: 'Yo',           utm: 'EPSG:32633'}
];
var BANDS = ['class', 'p_no', 'p_low', 'p_high', 'n_obs'];

var reserves = ee.FeatureCollection(ROOT + '/reserves');

RESERVES.forEach(function (r) {
  var geom = reserves.filter(ee.Filter.eq('reserve', r.id)).geometry();
  var img = ee.Image.cat(YEARS.map(function (y) {
    var m = ee.Image(ROOT + '/lulc/' + r.id + '_' + y);
    return m.select(BANDS).round().min(254)
      .rename(BANDS.map(function (b) { return b + '_' + y; }));
  }));
  Export.image.toDrive({
    image: img.clip(geom).unmask(255, false).toByte(),  // clip first so outside pixels are nodata
    description: 'lulc_' + r.id,
    folder: FOLDER,
    fileNamePrefix: 'lulc_' + r.id,
    region: geom,
    crs: r.utm,
    crsTransform: [10, 0, 0, 0, -10, 0],
    maxPixels: 1e13,
    formatOptions: {cloudOptimized: true, noData: 255}
  });
});
