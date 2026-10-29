/**
 * 01_build_training_samples.js
 *
 * Builds 500 stratified training samples per reserve and epoch (8 reserves x
 * 2 epochs = 8,000 samples) with labels assigned by multi-source consensus,
 * and splits them into training and test sets by spatial blocks.
 *
 * Label rule (all conditions must hold for a pixel to be eligible):
 *   1. Dynamic World V1 scenes from the same dry-season window are averaged and
 *      their nine class probabilities are summed into the three study classes
 *        High vegetation = trees
 *        Low vegetation  = grass + shrub_and_scrub + crops + flooded_vegetation
 *        No vegetation   = water + built + bare (+ snow_and_ice)
 *      The label is the class with the highest mean probability, which must be
 *      >= DW_MIN_PROB.
 *   2. At least DW_MIN_CONSISTENCY of the individual DW scenes agree with it,
 *      and at least DW_MIN_SCENES scenes are available.
 *   3. The label agrees with ESA WorldCover (2020 map for the 2019 epoch, 2021
 *      map for the 2024 epoch). Pixels with Hansen forest loss between the epoch
 *      and the WorldCover year are exempt in the direction the loss implies, so
 *      that recently cleared areas stay eligible. For the 2024 epoch, pixels
 *      labelled High vegetation must have no loss in 2021 to 2024.
 *   4. All eight neighbours carry the same DW label (homogeneous 3 x 3 window),
 *      which avoids mixed pixels at class edges.
 *   5. A valid Sentinel-2 composite value exists, and the pixel lies more than
 *      MIN_DIST_FROM_REF metres from every photo-interpreted reference point,
 *      so the original reference points remain fully independent for
 *      validation.
 *
 * Eligible pixels are then sampled with equal allocation per class
 * (CLASS_POINTS). Each sample inherits a random number from the 1 km block it
 * falls in; blocks with a value below 1 - TEST_FRACTION go to training and the
 * rest to testing, so train and test samples never share a block.
 *
 * Outputs
 *   Asset  SAMPLES_ASSET                       all samples with labels and split
 *   Drive  training_samples_consensus.csv      same, as CSV
 *   Drive  sample_counts.csv                   counts per reserve, epoch, class, split
 *   Drive  verification_blind.kml / .csv       ~10 % subset without labels, for
 *                                              checking in Google Earth Pro
 *   Drive  verification_key.csv                consensus labels of that subset
 */

// ---------------------------------------------------------------------------
var REPO = 'users/YOUR_USERNAME/nigeria-forest-fragmentation';
var cfg = require(REPO + ':config');
var lib = require(REPO + ':lib/common');
// ---------------------------------------------------------------------------

// Dynamic World label (0..8) to study class (0 no, 1 low, 2 high)
var DW_BANDS = ['water', 'trees', 'grass', 'flooded_vegetation', 'crops',
                'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];
var DW_LABELS   = [0, 1, 2, 3, 4, 5, 6, 7, 8];
var DW_TO_CLASS = [0, 2, 1, 1, 1, 1, 0, 0, 0];
// ESA WorldCover code to study class. Tree cover and mangroves are High,
// shrub, grass, cropland, herbaceous wetland and moss are Low, built-up,
// bare, snow and water are No vegetation.
var WC_FROM = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100];
var WC_TO   = [ 2,  1,  1,  1,  0,  0,  0,  0,  1,  2,   1];

// Agreement with WorldCover, with an exemption where Hansen recorded forest
// loss between the epoch and the WorldCover year.
function wcAgreement(label, wc, loss, epoch) {
  if (epoch.lossSide === 'before') {
    // Epoch precedes the WorldCover map: forest seen in the epoch may have been lost since.
    return wc.eq(label).or(loss.and(label.eq(2)));
  }
  // Epoch follows the WorldCover map: cleared land may still be forest in WorldCover,
  // and forest labels are not accepted where loss was recorded.
  return wc.eq(label).or(loss.and(label.lt(2))).and(label.neq(2).or(loss.not()));
}

// Random value per 1 km block, drawn separately for each class so that every
// class is split roughly 70/30 even when its samples cluster in few blocks.
function blockRandom(label, reserve) {
  var p = ee.Projection(reserve.utm).atScale(cfg.BLOCK_SIZE);
  var r = cfg.CLASS_VALUES.map(function (c) {
    return ee.Image.random(cfg.SEED + 1000 * (c + 1)).reproject(p);
  });
  return ee.Image(0).where(label.eq(0), r[0]).where(label.eq(1), r[1])
    .where(label.eq(2), r[2]).rename('block_rand');
}

// Returns one label image per tier (bands: class, tier, dw_prob, dw_consistency,
// dw_scenes, wc_class, block_rand, block_id). Tiers are mutually exclusive.
//   Tier 1 (strict)      grouped DW probability >= DW_MIN_PROB, consistency >=
//                        DW_MIN_CONSISTENCY, WorldCover agreement, homogeneous 3 x 3.
//   Tier 2 (relaxed)     label from the single most probable DW class (of nine),
//                        its probability >= TIER2_MIN_PROB, consistency >=
//                        TIER2_MIN_CONSISTENCY, WorldCover agreement; no window rule.
//   Tier 3 (WorldCover)  WorldCover label, with grouped DW probability for that
//                        class >= TIER3_MIN_DW_PROB (DW does not contradict it).
// Tiers 2 and 3 are only drawn on when a class falls short of its quota in tier 1.
function labelTiers(reserve, epoch, geom) {
  var proj = lib.gridProjection(cfg, reserve);
  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
    .filterBounds(geom).filterDate(epoch.start, epoch.end);
  var dwMean = dw.select(DW_BANDS).mean().reproject(proj);
  var probs = ee.Image.cat([
    dwMean.select(['water', 'built', 'bare', 'snow_and_ice']).reduce(ee.Reducer.sum()),
    dwMean.select(['grass', 'shrub_and_scrub', 'crops', 'flooded_vegetation']).reduce(ee.Reducer.sum()),
    dwMean.select('trees')
  ]).rename(['p_no', 'p_low', 'p_high']);
  var dwScenes = dw.select('label').count().rename('dw_scenes').reproject(proj);
  var consistency = function (label) {
    return dw.select('label').map(function (img) {
      return img.remap(DW_LABELS, DW_TO_CLASS).eq(label);
    }).mean().rename('dw_consistency').reproject(proj);
  };

  var wc = ee.ImageCollection(epoch.worldcover).first().select('Map')
    .remap(WC_FROM, WC_TO).rename('wc_class');
  var lossYear = ee.Image(cfg.GFC_ASSET).select('lossyear').unmask(0);
  var loss = lossYear.gte(epoch.lossYears[0]).and(lossYear.lte(epoch.lossYears[1]));

  var s2valid = lib.predictors(cfg, geom, epoch).stack.select('B4').mask().gt(0);
  var refPts = ee.FeatureCollection(cfg.REFERENCE_ASSET).filter(ee.Filter.eq('reserve', reserve.id));
  var farFromRef = refPts.distance(cfg.MIN_DIST_FROM_REF * 2)
    .gte(cfg.MIN_DIST_FROM_REF).unmask(1).reproject(proj);
  var base = s2valid.and(farFromRef).and(dwScenes.gte(cfg.DW_MIN_SCENES));

  var xy = ee.Image.pixelCoordinates(ee.Projection(reserve.utm));
  var blockId = xy.select('x').divide(cfg.BLOCK_SIZE).floor().multiply(10000)
    .add(xy.select('y').divide(cfg.BLOCK_SIZE).floor()).rename('block_id');

  // Tier 1
  var l1 = probs.toArray().arrayArgmax().arrayGet([0]).reproject(proj);
  var p1 = probs.reduce(ee.Reducer.max());
  var c1 = consistency(l1);
  var homogeneous = l1.reduceNeighborhood({
    reducer: ee.Reducer.countDistinctNonNull(), kernel: ee.Kernel.square(1)
  }).reproject(proj).eq(1);
  var e1 = base.and(p1.gte(cfg.DW_MIN_PROB)).and(c1.gte(cfg.DW_MIN_CONSISTENCY))
    .and(wcAgreement(l1, wc, loss, epoch)).and(homogeneous);

  // Tier 2
  var top = dwMean.toArray().arrayArgmax().arrayGet([0]).reproject(proj);
  var l2 = top.remap(DW_LABELS, DW_TO_CLASS);
  var p2 = dwMean.reduce(ee.Reducer.max());
  var c2 = consistency(l2);
  var e2 = base.and(e1.not()).and(p2.gte(cfg.TIER2_MIN_PROB)).and(c2.gte(cfg.TIER2_MIN_CONSISTENCY))
    .and(wcAgreement(l2, wc, loss, epoch));

  // Tier 3
  var l3 = wc;
  var p3 = probs.select('p_no').multiply(l3.eq(0))
    .add(probs.select('p_low').multiply(l3.eq(1)))
    .add(probs.select('p_high').multiply(l3.eq(2)));
  var c3 = consistency(l3);
  var e3 = base.and(e1.not()).and(e2.not()).and(p3.gte(cfg.TIER3_MIN_DW_PROB))
    .and(epoch.lossSide === 'after' ? l3.neq(2).or(loss.not()) : ee.Image(1));

  function pack(label, prob, cons, eligible, tier) {
    label = label.rename('class').int();
    return label.updateMask(eligible)
      .addBands([ee.Image(tier).rename('tier'), prob.rename('dw_prob'), cons, dwScenes, wc,
                 blockRandom(label, reserve), blockId])
      .clip(geom);
  }
  return [pack(l1, p1, c1, e1, 1), pack(l2, p2, c2, e2, 2), pack(l3, p3, c3, e3, 3)];
}

function samplesFor(reserve, epoch) {
  var geom = lib.reserveGeometry(cfg, reserve);
  var tiers = labelTiers(reserve, epoch, geom);
  var drawn = tiers.map(function (img, t) {
    return img.stratifiedSample({
      numPoints: 0, classBand: 'class', region: geom, scale: cfg.SCALE,
      projection: reserve.utm, seed: cfg.SEED + t, classValues: cfg.CLASS_VALUES,
      classPoints: cfg.CLASS_POINTS, dropNulls: true, tileScale: 8, geometries: true
    });
  });
  // Fill each class quota from tier 1 first, then tier 2, then tier 3
  var perClass = cfg.CLASS_VALUES.map(function (c, ci) {
    var quota = ee.Number(cfg.CLASS_POINTS[ci]);
    var taken = ee.FeatureCollection([]);
    drawn.forEach(function (fc) {
      var need = quota.subtract(taken.size()).max(0);
      // toList() rejects a count of 0, so ask for at least 1 and slice back to need
      var add = ee.FeatureCollection(fc.filter(ee.Filter.eq('class', c)).toList(need.max(1)).slice(0, need));
      taken = taken.merge(add);
    });
    return taken;
  });
  var s = ee.FeatureCollection(perClass).flatten();
  var names = ee.List(cfg.CLASS_NAMES);
  return s.map(function (f) {
    var ll = f.geometry().transform('EPSG:4326', 1).coordinates();
    var split = ee.Algorithms.If(
      ee.Number(f.get('block_rand')).lt(1 - cfg.TEST_FRACTION), 'train', 'test');
    return f.set({
      reserve: reserve.id,
      biome: reserve.biome,
      year: epoch.year,
      class_name: names.get(ee.Number(f.get('class')).int()),
      split: split,
      sample_id: ee.String('S_' + reserve.id + '_' + epoch.year + '_T')
        .cat(ee.Number(f.get('tier')).format('%d')).cat('_').cat(f.id()),
      lon: ll.get(0),
      lat: ll.get(1)
    });
  });
}

// ---------------------------------------------------------------------------
// Build all samples
// ---------------------------------------------------------------------------
var reserves = lib.reservesToRun(cfg);
var parts = [];
var counts = [];
var verify = [];

reserves.forEach(function (r) {
  cfg.EPOCHS.forEach(function (e) {
    var s = samplesFor(r, e);
    parts.push(s);
    cfg.CLASS_VALUES.forEach(function (c) {
      var sc = s.filter(ee.Filter.eq('class', c));
      counts.push(ee.Feature(null, {
        reserve: r.id, year: e.year, class: c, class_name: cfg.CLASS_NAMES[c],
        n_total: sc.size(),
        n_train: sc.filter(ee.Filter.eq('split', 'train')).size(),
        n_test: sc.filter(ee.Filter.eq('split', 'test')).size(),
        n_tier1: sc.filter(ee.Filter.eq('tier', 1)).size(),
        n_tier2: sc.filter(ee.Filter.eq('tier', 2)).size(),
        n_tier3: sc.filter(ee.Filter.eq('tier', 3)).size()
      }));
      verify.push(sc.randomColumn('vrand', cfg.SEED + 7).sort('vrand')
        .limit(cfg.N_VERIFY_PER_CLASS));
    });
  });
});

var samples = ee.FeatureCollection(parts).flatten();
var countTable = ee.FeatureCollection(counts);
var verifySet = ee.FeatureCollection(verify).flatten();

print('Sample counts per reserve, epoch and class', countTable);
print('Total samples', samples.size());

// ---------------------------------------------------------------------------
// Quick visual check for the first reserve and epoch
// ---------------------------------------------------------------------------
var r0 = reserves[0], e0 = cfg.EPOCHS[0];
var g0 = lib.reserveGeometry(cfg, r0);
Map.centerObject(g0, 12);
Map.addLayer(lib.predictors(cfg, g0, e0).stack, {bands: ['B8', 'B4', 'B3'], min: 0, max: 0.35},
  r0.id + ' ' + e0.year + ' S2 false colour');
labelTiers(r0, e0, g0).forEach(function (img, t) {
  Map.addLayer(img.select('class'), {min: 0, max: 2, palette: cfg.CLASS_PALETTE},
    'Eligible pixels, tier ' + (t + 1), t === 0);
});
Map.addLayer(samples.filter(ee.Filter.eq('reserve', r0.id)).filter(ee.Filter.eq('year', e0.year)),
  {color: 'ff00ff'}, 'Samples');

// ---------------------------------------------------------------------------
// Exports
// ---------------------------------------------------------------------------
Export.table.toAsset({
  collection: samples,
  description: 'training_samples_consensus_asset',
  assetId: cfg.SAMPLES_ASSET
});

Export.table.toDrive({
  collection: samples,
  description: 'training_samples_consensus_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'training_samples_consensus',
  fileFormat: 'CSV',
  selectors: ['sample_id', 'reserve', 'biome', 'year', 'class', 'class_name', 'tier', 'split',
              'block_id', 'block_rand', 'dw_prob', 'dw_consistency', 'dw_scenes',
              'wc_class', 'lon', 'lat']
});

Export.table.toDrive({
  collection: countTable,
  description: 'sample_counts_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'sample_counts',
  fileFormat: 'CSV',
  selectors: ['reserve', 'year', 'class', 'class_name', 'n_total', 'n_train', 'n_test',
              'n_tier1', 'n_tier2', 'n_tier3']
});

// Blind verification set: no labels, so the interpreter is not anchored by them.
var blind = verifySet.map(function (f) {
  return f.select(['sample_id', 'reserve', 'year', 'lon', 'lat'])
    .set('name', f.get('sample_id'), 'verified_class', '');
});

Export.table.toDrive({
  collection: blind,
  description: 'verification_blind_kml',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'verification_blind',
  fileFormat: 'KML'
});

Export.table.toDrive({
  collection: blind,
  description: 'verification_blind_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'verification_blind',
  fileFormat: 'CSV',
  selectors: ['sample_id', 'reserve', 'year', 'lon', 'lat', 'verified_class']
});

Export.table.toDrive({
  collection: verifySet,
  description: 'verification_key_csv',
  folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'verification_key',
  fileFormat: 'CSV',
  selectors: ['sample_id', 'reserve', 'year', 'class', 'class_name', 'tier', 'dw_prob', 'dw_consistency']
});
