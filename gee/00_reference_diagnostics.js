/**
 * 00_reference_diagnostics.js
 *
 * Checks what independent global products say at each photo-interpreted
 * reference point. Used to judge whether the class definitions behind the
 * manual labels (especially High vegetation in the savanna reserves) match
 * what Dynamic World, ESA WorldCover and canopy-height products record.
 *
 * Output
 *   Drive  reference_diagnostics.csv   one row per reference point with
 *     dw_<class>            mean Dynamic World probability of each of the nine classes
 *     dw_p_no/low/high      the same grouped into the three study classes
 *     wc_code               ESA WorldCover code (2020 map for 2019, 2021 map for 2024)
 *     canopy_height_m       ETH global canopy height 2020 (10 m)
 *     treecover2000         Hansen tree canopy cover in 2000 (%)
 *     NDVI                  dry-season Sentinel-2 composite of the same epoch
 */

// ---------------------------------------------------------------------------
var REPO = 'users/YOUR_USERNAME/nigeria-forest-fragmentation';
var cfg = require(REPO + ':config');
var lib = require(REPO + ':lib/common');
// ---------------------------------------------------------------------------

var DW_BANDS = ['water', 'trees', 'grass', 'flooded_vegetation', 'crops',
                'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];
var canopyHeight = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1').rename('canopy_height_m');
var treecover = ee.Image(cfg.GFC_ASSET).select('treecover2000').rename('treecover2000');

var rows = [];
cfg.RESERVES.forEach(function (r) {
  var geom = lib.reserveGeometry(cfg, r);
  cfg.EPOCHS.forEach(function (e) {
    var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
      .filterBounds(geom).filterDate(e.start, e.end).select(DW_BANDS).mean();
    var grouped = ee.Image.cat([
      dw.select(['water', 'built', 'bare', 'snow_and_ice']).reduce(ee.Reducer.sum()).rename('dw_p_no'),
      dw.select(['grass', 'shrub_and_scrub', 'crops', 'flooded_vegetation']).reduce(ee.Reducer.sum()).rename('dw_p_low'),
      dw.select('trees').rename('dw_p_high')
    ]);
    var wc = ee.ImageCollection(e.worldcover).first().select('Map').rename('wc_code');
    var ndvi = lib.predictors(cfg, geom, e).stack.select('NDVI');
    var stack = dw.rename(DW_BANDS.map(function (b) { return 'dw_' + b; }))
      .addBands([grouped, wc, canopyHeight, treecover, ndvi]);
    rows.push(stack.sampleRegions({
      collection: lib.referencePoints(cfg, r, e.year),
      properties: ['ref_id', 'reserve', 'year', 'class', 'class_name'],
      scale: cfg.SCALE, projection: r.utm, tileScale: 4
    }));
  });
});
var table = ee.FeatureCollection(rows).flatten();

print('First rows', table.limit(5));

Export.table.toDrive({
  collection: table,
  description: 'reference_diagnostics_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'reference_diagnostics',
  fileFormat: 'CSV',
  selectors: ['ref_id', 'reserve', 'year', 'class', 'class_name',
              'dw_p_no', 'dw_p_low', 'dw_p_high'].concat(DW_BANDS.map(function (b) { return 'dw_' + b; }))
    .concat(['wc_code', 'canopy_height_m', 'treecover2000', 'NDVI'])
});
