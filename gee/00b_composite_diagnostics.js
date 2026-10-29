/**
 * 00b_composite_diagnostics.js
 *
 * Checks whether differences between the 2019 and 2024 dry-season composites
 * come from fire, image availability or phenology rather than land-cover change.
 * One row per reserve and epoch.
 *
 * Output  Drive  composite_diagnostics.csv
 *   reserve_km2          reserve area
 *   burned_km2, burned_pct   area burned during the composite window (MODIS MCD64A1, 500 m)
 *   n_scenes             Sentinel-2 scenes in the window (after the scene cloud filter)
 *   clear_obs_mean       mean number of clear (Cloud Score+) observations per pixel
 *   clear_obs_min_p10    10th percentile of clear observations per pixel
 *   ndvi_median_p50      median of the composite NDVI over the reserve
 *   ndvi_median_p25/p75  quartiles of the composite NDVI
 *   ndvi_max_p50         median of the per-pixel maximum NDVI (greenest observation)
 *   nbr_median_p50       median of the composite NBR (low values indicate burn scars)
 *   high/low/no_km2      class areas in the step-02 map (for side-by-side reading)
 */

var ROOT = 'projects/YOUR_CLOUD_PROJECT/assets/nigeria_fragmentation';   // same as ASSET_ROOT in config.js
var EPOCHS = [
  {year: 2019, start: '2018-11-01', end: '2019-04-01'},
  {year: 2024, start: '2023-11-01', end: '2024-04-01'}
];
var RESERVES = ['Adoru', 'Akanto', 'Gambari', 'Maradun', 'Mayo_Ndaga',
                'Stubbs_Creek', 'Upper_Orashi', 'Yo'];
var CS = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var reserves = ee.FeatureCollection(ROOT + '/reserves');

var rows = [];
RESERVES.forEach(function (rid) {
  var geom = reserves.filter(ee.Filter.eq('reserve', rid)).geometry();
  EPOCHS.forEach(function (e) {
    // Fire
    var burned = ee.ImageCollection('MODIS/061/MCD64A1').filterDate(e.start, e.end)
      .select('BurnDate').max().gt(0).unmask(0);
    var burnedKm2 = ee.Number(burned.multiply(ee.Image.pixelArea()).divide(1e6)
      .reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: 500, maxPixels: 1e10})
      .get('BurnDate'));
    var areaKm2 = geom.area(1).divide(1e6);

    // Sentinel-2 availability and greenness
    var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(geom)
      .filterDate(e.start, e.end).filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 60))
      .linkCollection(CS, ['cs_cdf'])
      .map(function (i) {
        var m = i.updateMask(i.select('cs_cdf').gte(0.6)).divide(10000);
        return m.normalizedDifference(['B8', 'B4']).rename('NDVI')
          .addBands(m.normalizedDifference(['B8', 'B12']).rename('NBR'));
      });
    var stats = ee.Image.cat([
      col.select('NDVI').count().rename('nobs'),
      col.select('NDVI').median().rename('ndvi_med'),
      col.select('NDVI').max().rename('ndvi_max'),
      col.select('NBR').median().rename('nbr_med')
    ]).reduceRegion({
      reducer: ee.Reducer.percentile([10, 25, 50, 75]).combine(ee.Reducer.mean(), '', true),
      geometry: geom, scale: 30, maxPixels: 1e10, tileScale: 4
    });

    // Class areas of the step-02 map
    var lulc = ee.Image(ROOT + '/lulc/' + rid + '_' + e.year).select('class');
    var cls = ee.Image.pixelArea().divide(1e6).addBands(lulc).reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
      geometry: geom, scale: 10, maxPixels: 1e13, tileScale: 4
    });
    var grp = ee.List(cls.get('groups'));
    var byClass = ee.Dictionary.fromLists(
      grp.map(function (g) { return ee.Number(ee.Dictionary(g).get('class')).int().format(); }),
      grp.map(function (g) { return ee.Dictionary(g).get('sum'); }));
    var areaOf = function (c) { return byClass.get(String(c), 0); };

    rows.push(ee.Feature(null, {
      reserve: rid, year: e.year,
      reserve_km2: areaKm2,
      burned_km2: burnedKm2,
      burned_pct: burnedKm2.divide(areaKm2).multiply(100),
      n_scenes: col.size(),
      clear_obs_mean: stats.get('nobs_mean'),
      clear_obs_min_p10: stats.get('nobs_p10'),
      ndvi_median_p25: stats.get('ndvi_med_p25'),
      ndvi_median_p50: stats.get('ndvi_med_p50'),
      ndvi_median_p75: stats.get('ndvi_med_p75'),
      ndvi_max_p50: stats.get('ndvi_max_p50'),
      nbr_median_p50: stats.get('nbr_med_p50'),
      no_km2: areaOf(0), low_km2: areaOf(1), high_km2: areaOf(2)
    }));
  });
});

var table = ee.FeatureCollection(rows);
print('Composite diagnostics', table);

Export.table.toDrive({
  collection: table,
  description: 'composite_diagnostics_csv',
  folder: 'nigeria_fragmentation_gee',
  fileNamePrefix: 'composite_diagnostics',
  fileFormat: 'CSV',
  selectors: ['reserve', 'year', 'reserve_km2', 'burned_km2', 'burned_pct', 'n_scenes',
              'clear_obs_mean', 'clear_obs_min_p10', 'ndvi_median_p25', 'ndvi_median_p50',
              'ndvi_median_p75', 'ndvi_max_p50', 'nbr_median_p50', 'no_km2', 'low_km2', 'high_km2']
});
