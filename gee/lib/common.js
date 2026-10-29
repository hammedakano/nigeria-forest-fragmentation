/**
 * lib/common.js
 * Shared functions: reserve geometry, Sentinel-2 composite, predictors,
 * reference points and accuracy summaries.
 *
 * Loaded with:
 *   var lib = require('users/<YOUR_USERNAME>/<YOUR_REPO>:lib/common');
 * Every function takes the config object as its first argument so that the
 * module itself holds no settings.
 */

// ---------------------------------------------------------------------------
// Reserves
// ---------------------------------------------------------------------------
exports.reservesToRun = function (cfg) {
  if (!cfg.RUN_ONLY || cfg.RUN_ONLY.length === 0) return cfg.RESERVES;
  return cfg.RESERVES.filter(function (r) { return cfg.RUN_ONLY.indexOf(r.id) >= 0; });
};

exports.reserveGeometry = function (cfg, reserve) {
  var fc;
  if (cfg.RESERVE_SOURCE === 'wdpa') {
    fc = ee.FeatureCollection('WCMC/WDPA/current/polygons')
      .filter(ee.Filter.eq('WDPAID', reserve.wdpaid));
  } else {
    fc = ee.FeatureCollection(cfg.RESERVES_ASSET)
      .filter(ee.Filter.eq('reserve', reserve.id));
  }
  return fc.geometry();
};

// 10 m processing grid in the reserve's UTM zone. Grid origin is the UTM
// origin, so every image produced for a reserve shares the same pixel edges.
exports.gridProjection = function (cfg, reserve) {
  return ee.Projection(reserve.utm, [cfg.SCALE, 0, 0, 0, -cfg.SCALE, 0]);
};

// ---------------------------------------------------------------------------
// Sentinel-2 composite
// ---------------------------------------------------------------------------
function maskQA60(img) {
  var qa = img.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  return img.updateMask(mask);
}

exports.s2Collection = function (cfg, geom, epoch) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(geom)
    .filterDate(epoch.start, epoch.end)
    .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cfg.MAX_SCENE_CLOUD));

  if (cfg.CLOUD_MASK === 'CSPLUS') {
    var cs = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
    s2 = s2.linkCollection(cs, ['cs_cdf']).map(function (img) {
      return img.updateMask(img.select('cs_cdf').gte(cfg.CS_THRESHOLD));
    });
  } else {
    s2 = s2.map(maskQA60);
  }
  return s2.map(function (img) {
    return img.select(cfg.BANDS).divide(10000)
      .copyProperties(img, ['system:time_start']);
  });
};

// Median composite plus the ten predictors of Table 1 and a count band that
// records how many clear observations went into each pixel.
exports.predictors = function (cfg, geom, epoch) {
  var col = exports.s2Collection(cfg, geom, epoch);
  var comp = col.median();
  var ndvi = comp.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = comp.normalizedDifference(['B3', 'B8']).rename('NDWI');
  var nbr  = comp.normalizedDifference(['B8', 'B12']).rename('NBR');
  var savi = comp.expression('((NIR - RED) / (NIR + RED + L)) * (1 + L)', {
    NIR: comp.select('B8'), RED: comp.select('B4'), L: 0.5
  }).rename('SAVI');
  var stack = comp.select(cfg.BANDS).addBands([ndvi, ndwi, nbr, savi]);
  if (cfg.SEASONAL_FEATURES && cfg.SEASONAL_FEATURES.length > 0) {
    // Per-observation NDVI and NBR, summarised by percentiles over the window.
    // High percentiles keep woody vegetation visible after grass senescence and
    // burning, so maps from years with different fire and drying timing remain
    // comparable. Low percentiles capture the driest or burnt state.
    var vi = col.map(function (img) {
      return img.normalizedDifference(['B8', 'B4']).rename('NDVI')
        .addBands(img.normalizedDifference(['B8', 'B12']).rename('NBR'));
    });
    var pct = vi.reduce(ee.Reducer.percentile([10, 90]));   // NDVI_p10, NDVI_p90, NBR_p10, NBR_p90
    stack = stack.addBands(pct.select(cfg.SEASONAL_FEATURES));
  }
  var nObs = col.select('B4').count().rename('n_obs');
  return {stack: stack.clip(geom), nObs: nObs.clip(geom), nScenes: col.size()};
};

exports.predictorNames = function (cfg) {
  return cfg.BANDS.concat(cfg.INDICES).concat(cfg.SEASONAL_FEATURES || []);
};

// ---------------------------------------------------------------------------
// Reference data
// ---------------------------------------------------------------------------
exports.referencePoints = function (cfg, reserve, year) {
  return ee.FeatureCollection(cfg.REFERENCE_ASSET)
    .filter(ee.Filter.eq('reserve', reserve.id))
    .filter(ee.Filter.eq('year', year))
    .map(function (f) { return f.set('class', ee.Number(f.get('class')).int()); });
};

// Consensus samples for a reserve-epoch, with verified labels applied if a
// verification table is configured.
exports.trainingSamples = function (cfg, reserve, year) {
  var s = ee.FeatureCollection(cfg.SAMPLES_ASSET)
    .filter(ee.Filter.eq('reserve', reserve.id))
    .filter(ee.Filter.eq('year', year));
  if (cfg.VERIFIED_ASSET && cfg.VERIFIED_ASSET !== '') {
    var ver = ee.FeatureCollection(cfg.VERIFIED_ASSET);
    var joined = ee.Join.saveFirst('ver').apply(s, ver,
      ee.Filter.equals({leftField: 'sample_id', rightField: 'sample_id'}));
    var notChecked = ee.Join.inverted().apply(s, ver,
      ee.Filter.equals({leftField: 'sample_id', rightField: 'sample_id'}));
    var checked = ee.FeatureCollection(joined).map(function (f) {
      var v = ee.Number(ee.Feature(f.get('ver')).get('verified_class')).int();
      return f.set('class_consensus', f.get('class'), 'class', v, 'verified', 1);
    }).filter(ee.Filter.gte('class', 0));
    s = notChecked.map(function (f) { return f.set('verified', 0); }).merge(checked);
  }
  return s;
};

// ---------------------------------------------------------------------------
// Accuracy summaries
// ---------------------------------------------------------------------------
// Converts a 3 x 3 error matrix (rows = reference, columns = map) into a flat
// feature: OA with 95 % CI, kappa, PA, UA, F1 and every matrix cell.
exports.accuracyFeature = function (cm, props) {
  var arr = ee.ConfusionMatrix(cm).array();
  var n = ee.Number(arr.reduce(ee.Reducer.sum(), [0, 1]).get([0, 0]));
  var diag = ee.Number(arr.matrixDiagonal().reduce(ee.Reducer.sum(), [0]).get([0, 0]));
  var oa = diag.divide(n.max(1));
  var ci = oa.multiply(ee.Number(1).subtract(oa)).divide(n.max(1)).sqrt().multiply(1.96);
  var rowSum = arr.reduce(ee.Reducer.sum(), [1]);   // reference totals (3 x 1)
  var colSum = arr.reduce(ee.Reducer.sum(), [0]);   // map totals (1 x 3)
  var pe = ee.Number(rowSum.matrixMultiply(colSum).matrixDiagonal()
    .reduce(ee.Reducer.sum(), [0]).get([0, 0])).divide(n.max(1).pow(2));
  var kappa = oa.subtract(pe).divide(ee.Number(1).subtract(pe).max(1e-9));
  var out = ee.Dictionary(props).combine({
    n: n, OA: oa, OA_ci95: ci, kappa: kappa
  });
  var names = ['no', 'low', 'high'];
  for (var i = 0; i < 3; i++) {
    var tp = ee.Number(arr.get([i, i]));
    var pa = tp.divide(ee.Number(rowSum.get([i, 0])).max(1e-9));
    var ua = tp.divide(ee.Number(colSum.get([0, i])).max(1e-9));
    var f1 = pa.multiply(ua).multiply(2).divide(pa.add(ua).max(1e-9));
    out = out.set('PA_' + names[i], pa).set('UA_' + names[i], ua).set('F1_' + names[i], f1)
             .set('n_ref_' + names[i], rowSum.get([i, 0]));
    for (var j = 0; j < 3; j++) {
      out = out.set('cm_r' + i + '_m' + j, arr.get([i, j]));
    }
  }
  return ee.Feature(null, out);
};

// Class areas (km2) inside a geometry for a single-band class image.
exports.classAreas = function (cfg, classImg, geom, reserve, props) {
  var grouped = ee.Image.pixelArea().divide(1e6).addBands(classImg.rename('class'))
    .reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
      geometry: geom, crs: reserve.utm, scale: cfg.SCALE, maxPixels: 1e13, tileScale: 4
    });
  var groups = ee.List(grouped.get('groups'));
  return ee.FeatureCollection(groups.map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, ee.Dictionary(props).combine({
      class: g.get('class'), area_km2: g.get('sum')
    }));
  }));
};
