/**
 * 02_classify_and_assess.js
 *
 * Year-specific Random Forest classification of each reserve (8 reserves x
 * 2 epochs) from the 500 consensus samples built in step 01 plus the
 * photo-interpreted reference points, followed by three accuracy assessments:
 *
 *   A. "blocked_holdout"  consensus samples in the held-out 1 km blocks
 *                         (about 30 % of the 500), spatially separated from
 *                         the training blocks.
 *   B. "independent_ref_cv" the photo-interpreted reference points (59 to 78
 *                         per reserve-epoch). They are also used for training,
 *                         so accuracy comes from N_FOLDS-fold cross-validation:
 *                         each point is predicted by a model that did not see it.
 *                         This is the headline accuracy reported in the paper.
 *   C. "consensus_only_ref" a model trained on consensus samples only, tested on
 *                         all reference points (for comparison).
 *
 * Outputs
 *   Asset  LULC_FOLDER/<reserve>_<year>   bands: class (0/1/2),
 *                                         p_no, p_low, p_high (RF vote share x 100),
 *                                         n_obs (clear S2 observations)
 *   Drive  rf_accuracy.csv                one row per reserve, epoch and assessment
 *   Drive  rf_feature_importance.csv      raw and normalised (%) Gini importance
 *   Drive  rf_class_areas.csv             class areas (km2) of the observed maps
 *   Drive  rf_training_counts.csv         samples actually used per class
 *
 * Create the LULC_FOLDER in the Assets tab before running the exports.
 */

// ---------------------------------------------------------------------------
var REPO = 'users/YOUR_USERNAME/nigeria-forest-fragmentation';
var cfg = require(REPO + ':config');
var lib = require(REPO + ':lib/common');
// ---------------------------------------------------------------------------

var PRED = lib.predictorNames(cfg);
var ORDER = cfg.CLASS_VALUES;

function runOne(reserve, epoch) {
  var geom = lib.reserveGeometry(cfg, reserve);
  var proj = lib.gridProjection(cfg, reserve);
  var pr = lib.predictors(cfg, geom, epoch);
  var stack = pr.stack;
  var meta = {reserve: reserve.id, biome: reserve.biome, year: epoch.year};

  // Samples with predictor values
  var samples = lib.trainingSamples(cfg, reserve, epoch.year);
  var sampleArgs = function (fc) {
    return stack.sampleRegions({
      collection: fc, properties: ['class', 'sample_id'],
      scale: cfg.SCALE, projection: reserve.utm, tileScale: 4
    }).filter(ee.Filter.notNull(PRED));
  };
  var train = sampleArgs(samples.filter(ee.Filter.eq('split', 'train')));
  var test  = sampleArgs(samples.filter(ee.Filter.eq('split', 'test')));
  var ref   = sampleArgs(lib.referencePoints(cfg, reserve, epoch.year));

  // Photo-interpreted reference points join the training so that the class
  // boundaries follow the study's own definitions. Each point is repeated
  // REF_WEIGHT times, which raises its influence relative to the consensus
  // samples (Random Forest in GEE has no case weights). Every point is assigned
  // to one of N_FOLDS folds for cross-validation.
  var refF = ref.randomColumn('fold_r', cfg.SEED + 11).map(function (f) {
    return f.set('fold', ee.Number(f.get('fold_r')).multiply(cfg.N_FOLDS).floor().int());
  });
  var repeatRef = function (fc) {
    var copies = [];
    for (var k = 0; k < cfg.REF_WEIGHT; k++) copies.push(fc);
    return ee.FeatureCollection(copies).flatten();
  };
  var trainRF = function (features) {
    return ee.Classifier.smileRandomForest({numberOfTrees: cfg.N_TREES, seed: cfg.SEED})
      .train({features: features, classProperty: 'class', inputProperties: PRED});
  };

  // Final model (all consensus training samples + all reference points) -> maps
  var rf = trainRF(train.merge(repeatRef(refF)));
  var rfProb = rf.setOutputMode('MULTIPROBABILITY');

  var cls = stack.classify(rf).rename('class').toByte();
  var prob = stack.classify(rfProb).arrayFlatten([['p_no', 'p_low', 'p_high']])
    .multiply(100).round().toByte();
  var out = cls.addBands(prob).addBands(pr.nObs.min(255).toByte())
    .clip(geom)
    .set({reserve: reserve.id, year: epoch.year, n_trees: cfg.N_TREES, seed: cfg.SEED,
          ref_weight: cfg.REF_WEIGHT, composite_start: epoch.start, composite_end: epoch.end,
          cloud_mask: cfg.CLOUD_MASK});

  // Accuracy
  // (1) independent_ref_cv: k-fold cross-validation. Each reference point is
  //     predicted by a model trained on the consensus samples and the reference
  //     points of the other folds only, so no point is ever predicted by a model
  //     that saw it. This is the headline accuracy.
  var cvPred = ee.FeatureCollection(ee.List.sequence(0, cfg.N_FOLDS - 1).map(function (k) {
    var inFold = refF.filter(ee.Filter.eq('fold', k));
    var model = trainRF(train.merge(repeatRef(refF.filter(ee.Filter.neq('fold', k)))));
    return inFold.classify(model);
  })).flatten();
  var cmCV = cvPred.errorMatrix('class', 'classification', ORDER);
  // (2) consensus_only_ref: model trained on consensus samples only, tested on all
  //     reference points (the earlier design, kept for comparison).
  var rfCons = trainRF(train);
  var cmCons = ref.classify(rfCons).errorMatrix('class', 'classification', ORDER);
  // (3) blocked_holdout: consensus samples in held-out spatial blocks.
  var cmTest = test.classify(rf).errorMatrix('class', 'classification', ORDER);
  var acc = ee.FeatureCollection([
    lib.accuracyFeature(cmCV,   ee.Dictionary(meta).set('assessment', 'independent_ref_cv')),
    lib.accuracyFeature(cmCons, ee.Dictionary(meta).set('assessment', 'consensus_only_ref')),
    lib.accuracyFeature(cmTest, ee.Dictionary(meta).set('assessment', 'blocked_holdout'))
  ]);

  // Feature importance, normalised to sum to 100 %
  var imp = ee.Dictionary(rf.explain().get('importance'));
  var total = ee.Number(imp.values().reduce(ee.Reducer.sum()));
  var impFc = ee.FeatureCollection(imp.keys().map(function (k) {
    var v = ee.Number(imp.get(k));
    return ee.Feature(null, ee.Dictionary(meta).combine({
      predictor: k, importance_raw: v, importance_pct: v.divide(total).multiply(100)
    }));
  }));

  // Areas and counts
  var areas = lib.classAreas(cfg, cls, geom, reserve, meta);
  var countsFc = ee.FeatureCollection(ORDER.map(function (c) {
    return ee.Feature(null, ee.Dictionary(meta).combine({
      class: c,
      n_train: train.filter(ee.Filter.eq('class', c)).size(),
      n_test: test.filter(ee.Filter.eq('class', c)).size(),
      n_ref: ref.filter(ee.Filter.eq('class', c)).size()
    }));
  }));
  var scenes = ee.Feature(null, ee.Dictionary(meta).set('n_scenes', pr.nScenes));

  return {image: out, proj: proj, geom: geom, acc: acc, imp: impFc, areas: areas,
          counts: countsFc, scenes: scenes};
}

// ---------------------------------------------------------------------------
var reserves = lib.reservesToRun(cfg);
var accAll = [], impAll = [], areaAll = [], countAll = [], sceneAll = [];

reserves.forEach(function (r) {
  cfg.EPOCHS.forEach(function (e) {
    var res = runOne(r, e);
    accAll.push(res.acc);
    impAll.push(res.imp);
    areaAll.push(res.areas);
    countAll.push(res.counts);
    sceneAll.push(res.scenes);

    Export.image.toAsset({
      image: res.image,
      description: 'lulc_' + r.id + '_' + e.year,
      assetId: cfg.LULC_FOLDER + '/' + r.id + '_' + e.year,
      region: res.geom,
      crs: r.utm,
      crsTransform: [cfg.SCALE, 0, 0, 0, -cfg.SCALE, 0],
      maxPixels: 1e13,
      pyramidingPolicy: {'class': 'mode', '.default': 'mean'}
    });
  });
});

var accFc   = ee.FeatureCollection(accAll).flatten();
var impFc   = ee.FeatureCollection(impAll).flatten();
var areaFc  = ee.FeatureCollection(areaAll).flatten();
var countFc = ee.FeatureCollection(countAll).flatten();
var sceneFc = ee.FeatureCollection(sceneAll);

print('Accuracy', accFc);
print('Scenes per composite', sceneFc);

var ACC_COLS = ['reserve', 'biome', 'year', 'assessment', 'n', 'OA', 'OA_ci95', 'kappa',
  'PA_high', 'PA_low', 'PA_no', 'UA_high', 'UA_low', 'UA_no', 'F1_high', 'F1_low', 'F1_no',
  'n_ref_high', 'n_ref_low', 'n_ref_no',
  'cm_r0_m0', 'cm_r0_m1', 'cm_r0_m2', 'cm_r1_m0', 'cm_r1_m1', 'cm_r1_m2',
  'cm_r2_m0', 'cm_r2_m1', 'cm_r2_m2'];

Export.table.toDrive({collection: accFc, description: 'rf_accuracy_csv',
  folder: cfg.DRIVE_FOLDER, fileNamePrefix: 'rf_accuracy', fileFormat: 'CSV', selectors: ACC_COLS});
Export.table.toDrive({collection: impFc, description: 'rf_feature_importance_csv',
  folder: cfg.DRIVE_FOLDER, fileNamePrefix: 'rf_feature_importance', fileFormat: 'CSV',
  selectors: ['reserve', 'biome', 'year', 'predictor', 'importance_raw', 'importance_pct']});
Export.table.toDrive({collection: areaFc, description: 'rf_class_areas_csv',
  folder: cfg.DRIVE_FOLDER, fileNamePrefix: 'rf_class_areas', fileFormat: 'CSV',
  selectors: ['reserve', 'biome', 'year', 'class', 'area_km2']});
Export.table.toDrive({collection: countFc.merge(sceneFc), description: 'rf_training_counts_csv',
  folder: cfg.DRIVE_FOLDER, fileNamePrefix: 'rf_training_counts', fileFormat: 'CSV',
  selectors: ['reserve', 'biome', 'year', 'class', 'n_train', 'n_test', 'n_ref', 'n_scenes']});

// ---------------------------------------------------------------------------
// Preview the first reserve
// ---------------------------------------------------------------------------
var r0 = reserves[0];
var p0 = runOne(r0, cfg.EPOCHS[1]);
Map.centerObject(p0.geom, 12);
Map.addLayer(p0.image.select('class'), {min: 0, max: 2, palette: cfg.CLASS_PALETTE},
  r0.id + ' ' + cfg.EPOCHS[1].year + ' RF classes');
Map.addLayer(p0.image.select('p_high'), {min: 0, max: 100, palette: ['ffffff', '1a9850']},
  'P(high vegetation) %', false);
