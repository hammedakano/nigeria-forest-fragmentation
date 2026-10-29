/**
 * 03_change_analysis.js
 *
 * Post-classification comparison of the 2019 and 2024 maps produced in step 02.
 *
 * Outputs
 *   Drive  change_transitions.csv   area (km2), pixel count and row probability of
 *                                   every from-to class pair, per reserve
 *   Drive  obs_<reserve>.tif        GeoTIFF (local UTM, 10 m) with bands
 *                                   class_2019, class_2024, transition (3*from + to),
 *                                   p_high_2019, p_high_2024; nodata = 255.
 *                                   Input for R/02_landscape_metrics.R.
 */

// ---------------------------------------------------------------------------
var REPO = 'users/YOUR_USERNAME/nigeria-forest-fragmentation';
var cfg = require(REPO + ':config');
var lib = require(REPO + ':lib/common');
// ---------------------------------------------------------------------------

function lulc(reserve, year) {
  return ee.Image(cfg.LULC_FOLDER + '/' + reserve.id + '_' + year);
}

var reserves = lib.reservesToRun(cfg);
var rows = [];

reserves.forEach(function (r) {
  var geom = lib.reserveGeometry(cfg, r);
  var a = lulc(r, cfg.T0), b = lulc(r, cfg.T1);
  var c0 = a.select('class'), c1 = b.select('class');
  var trans = c0.multiply(3).add(c1).rename('transition').toByte();

  var grouped = ee.Image.pixelArea().divide(1e6).addBands(ee.Image(1).rename('n'))
    .addBands(trans)
    .reduceRegion({
      reducer: ee.Reducer.sum().repeat(2).group({groupField: 2, groupName: 'transition'}),
      geometry: geom, crs: r.utm, scale: cfg.SCALE, maxPixels: 1e13, tileScale: 4
    });
  var groups = ee.List(grouped.get('groups'));
  var fc = ee.FeatureCollection(groups.map(function (g) {
    g = ee.Dictionary(g);
    var t = ee.Number(g.get('transition'));
    var sums = ee.List(g.get('sum'));
    return ee.Feature(null, {
      reserve: r.id, biome: r.biome,
      from_class: t.divide(3).floor(), to_class: t.mod(3),
      area_km2: sums.get(0), n_pixels: sums.get(1)
    });
  }));
  // Row probabilities p(to | from)
  fc = fc.map(function (f) {
    var fromTotal = fc.filter(ee.Filter.eq('from_class', f.get('from_class')))
      .aggregate_sum('n_pixels');
    return f.set('p_row', ee.Number(f.get('n_pixels')).divide(fromTotal));
  });
  rows.push(fc);

  var img = ee.Image.cat([
    c0.rename('class_' + cfg.T0), c1.rename('class_' + cfg.T1), trans,
    a.select('p_high').rename('p_high_' + cfg.T0), b.select('p_high').rename('p_high_' + cfg.T1)
  ]).clip(geom).unmask(255, false).toByte();  // clip first so outside pixels are written as nodata, not 0

  Export.image.toDrive({
    image: img,
    description: 'obs_' + r.id,
    folder: cfg.DRIVE_FOLDER,
    fileNamePrefix: 'obs_' + r.id,
    region: geom,
    crs: r.utm,
    crsTransform: [cfg.SCALE, 0, 0, 0, -cfg.SCALE, 0],
    maxPixels: 1e13,
    formatOptions: {cloudOptimized: true, noData: 255}
  });
});

var table = ee.FeatureCollection(rows).flatten();
print('Transitions', table);

Export.table.toDrive({
  collection: table,
  description: 'change_transitions_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'change_transitions',
  fileFormat: 'CSV',
  selectors: ['reserve', 'biome', 'from_class', 'to_class', 'area_km2', 'n_pixels', 'p_row']
});

// Preview
var r0 = reserves[0];
var g0 = lib.reserveGeometry(cfg, r0);
Map.centerObject(g0, 12);
Map.addLayer(lulc(r0, cfg.T0).select('class'), {min: 0, max: 2, palette: cfg.CLASS_PALETTE}, r0.id + ' ' + cfg.T0);
Map.addLayer(lulc(r0, cfg.T1).select('class'), {min: 0, max: 2, palette: cfg.CLASS_PALETTE}, r0.id + ' ' + cfg.T1);
Map.addLayer(lulc(r0, cfg.T0).select('class').eq(2).and(lulc(r0, cfg.T1).select('class').lt(2)).selfMask(),
  {palette: ['e31a1c']}, 'High vegetation lost');
