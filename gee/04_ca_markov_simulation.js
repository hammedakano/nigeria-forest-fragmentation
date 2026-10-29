/**
 * 04_ca_markov_simulation.js
 *
 * Stochastic CA-Markov simulation run entirely in Earth Engine, with a Monte
 * Carlo ensemble that yields per-pixel probability and agreement maps (the
 * spatial uncertainty layer requested for the 2033 projection).
 *
 * Markov component
 *   P   = row-normalised cross-tabulation of the 2019 and 2024 maps (5-year step).
 *   A   = annual transition matrix, the 5th matrix root of P computed with a
 *         binomial series. Negative entries are set to zero and rows are
 *         renormalised. If A^5 does not reproduce P within ROOT_TOLERANCE, the
 *         first-order approximation A = I + (P - I) / 5 is used instead. The
 *         method used is written to markov_matrices.csv.
 *
 * Cellular automata component (one iteration = one year)
 *   For a pixel of class i, the probability of moving to class j (j != i) is
 *     p_ij(x) = A_ij * (n_j(x)^k + e) / mean_i(n_j^k + e)
 *   where n_j(x) is the share of the 8 Moore neighbours in class j (equal
 *   weights), k = NEIGHBOUR_POWER sharpens the preference for patch edges and
 *   e = NUCLEATION lets new patches appear away from existing ones,
 *   and mean_i(.) is the mean over all class-i pixels of the current map
 *   (re-estimated every iteration), so that the expected annual quantity of
 *   each transition equals the Markov rate. If p_i1 + p_i2 exceeds 0.99 both
 *   are scaled down. One uniform random draw per pixel decides the outcome.
 *   The draws are spatially correlated (Gaussian-smoothed noise, CA_CORR_SIGMA
 *   pixels). With independent draws the simulation scatters single-pixel
 *   transitions and roughly doubles patch and edge density within five years;
 *   CA_CORR_SIGMA is calibrated so that the hindcast reproduces the observed
 *   2024 patch and edge density (R/calibration/). No suitability, constraint
 *   or driver layers are used.
 *
 * Runs
 *   Hindcast  2019 map -> 5 iterations -> simulated 2024, validated against
 *             (a) the independent 2024 photo-interpreted reference points and
 *             (b) the observed 2024 map (pixel agreement and Pontius figure of
 *             merit on change pixels). A persistence null model (2019 map used
 *             as the 2024 prediction) is evaluated the same way.
 *   Forward   2024 map -> 4 iterations (2028) and 9 iterations (2033).
 *   Each is repeated N_RUNS times with different random seeds, starting from
 *   the classified map (CA_START = 'classified').
 *   With PROPAGATE_CLASSIFICATION_ERROR = true, N_RUNS realisations of the
 *   observed 2019 and 2024 maps are also drawn from the RF class probabilities
 *   and exported. R uses them for the classification-uncertainty intervals of
 *   the observed landscape metrics. They are not used as CA starting maps by
 *   default: in the swamp reserves the realisations have several times the
 *   patch density of the classified map, and the CA neighbourhood rule would
 *   smooth this noise away and simulate a spurious decline in fragmentation.
 *
 * Outputs
 *   Drive  cam_<reserve>.tif        GeoTIFF, local UTM, 10 m, nodata 255. Bands are
 *                                   listed in cam_band_index.csv:
 *                                   <prefix>_rNN       class map of run NN
 *                                   <prefix>_p_no/low/high  share of runs in each class (%)
 *                                   <prefix>_mode      modal class
 *                                   <prefix>_agree     share of runs equal to the mode (%)
 *                                   prefixes: h2024 (hindcast), y2028, y2033
 *                                   i2019_rNN, i2024_rNN  starting-map realisations
 *                                   (only with PROPAGATE_CLASSIFICATION_ERROR)
 *   Drive  cam_hindcast_accuracy_<reserve>.csv one row per run (+ mode, + persistence)
 *   Drive  cam_markov_matrices.csv   P, annual A, root method and error per reserve
 *   Drive  cam_quantity_check_<reserve>.csv  Markov-expected vs simulated class areas
 *   Drive  cam_band_index.csv        band order of the GeoTIFFs
 */

// ---------------------------------------------------------------------------
var REPO = 'users/YOUR_USERNAME/nigeria-forest-fragmentation';
var cfg = require(REPO + ':config');
var lib = require(REPO + ':lib/common');
// ---------------------------------------------------------------------------

var ORDER = cfg.CLASS_VALUES;
var INTERVAL = cfg.T1 - cfg.T0;
var KERNEL = ee.Kernel.fixed(3, 3, [[1, 1, 1], [1, 0, 1], [1, 1, 1]], -1, -1, false);
var I3 = ee.Array.identity(3);

// Constant arrays built client-side, so no scalar broadcasting is needed
function full3(v) { return ee.Array([[v, v, v], [v, v, v], [v, v, v]]); }
function col3(v) { return ee.Array([[v], [v], [v]]); }

function pad2(n) { return (n < 10 ? '0' : '') + n; }
function lulcFull(r, y) { return ee.Image(cfg.LULC_FOLDER + '/' + r.id + '_' + y); }
function lulc(r, y) { return lulcFull(r, y).select('class'); }

// ---------------------------------------------------------------------------
// Markov component
// ---------------------------------------------------------------------------
function countMatrix(c0, c1, geom, r) {
  var t = c0.multiply(3).add(c1).toInt().rename('t');
  var hist = ee.Dictionary(t.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(), geometry: geom,
    crs: r.utm, scale: cfg.SCALE, maxPixels: 1e13, tileScale: 4
  }).get('t'));
  var counts = ee.List.sequence(0, 8).map(function (k) {
    return ee.Number(hist.get(ee.Number(k).int().format(), 0));
  });
  return ee.Array(counts).reshape([3, 3]);
}

function rowNormalise(M) {
  var rows = M.reduce(ee.Reducer.sum(), [1]).max(col3(1e-12));
  return M.divide(rows.repeat(1, 3));
}

function matPow(M, n) {
  var out = I3;
  for (var i = 0; i < n; i++) out = out.matrixMultiply(M);
  return out;
}

function annualMatrix(P) {
  var Q = P.subtract(I3);
  var a = 1 / INTERVAL;
  var A = I3, term = I3, coef = 1;
  for (var k = 1; k <= cfg.ROOT_SERIES_TERMS; k++) {
    coef = coef * (a - k + 1) / k;
    term = term.matrixMultiply(Q);
    A = A.add(term.multiply(full3(coef)));
  }
  A = rowNormalise(A.max(full3(0)));
  var err = ee.Number(matPow(A, INTERVAL).subtract(P).abs()
    .reduce(ee.Reducer.max(), [0, 1]).get([0, 0]));
  var linear = I3.add(Q.divide(full3(INTERVAL)));
  var ok = err.lt(cfg.ROOT_TOLERANCE);
  return {
    A: ee.Array(ee.Algorithms.If(ok, A, linear)),
    method: ee.Algorithms.If(ok, 'matrix_root', 'linear_approximation'),
    err: err
  };
}

// ---------------------------------------------------------------------------
// CA component
// ---------------------------------------------------------------------------
// CA parameters of a reserve: per-reserve values in cfg.CA_PARAMS (calibrated
// on the hindcast, R/calibration/) override the defaults.
function caParams(r) {
  var d = {sigma: cfg.CA_CORR_SIGMA, nucleation: cfg.NUCLEATION, power: cfg.NEIGHBOUR_POWER};
  var o = (cfg.CA_PARAMS || {})[r.id] || {};
  Object.keys(o).forEach(function (k) { d[k] = o[k]; });
  return d;
}

function neighbourShares(M, proj, par) {
  return ee.Image.cat(ORDER.map(function (j) {
    return M.eq(j).reduceNeighborhood({reducer: ee.Reducer.mean(), kernel: KERNEL})
      .unmask(0).pow(par.power).rename('n' + j);
  })).reproject(proj).add(par.nucleation);
}

// mean_i(n_j + e) for every from-class i in the current map. Estimated on a
// NORM_SCALE grid (a systematic subsample of the 10 m pixels) to keep each
// iteration cheap; recomputed every iteration so that simulated quantities
// track the Markov rates as the landscape changes.
function normalisers(M, N, geom, r) {
  var res = ee.Dictionary(N.addBands(M.rename('cls')).reduceRegion({
    reducer: ee.Reducer.mean().repeat(3).group({groupField: 3, groupName: 'cls'}),
    geometry: geom, crs: r.utm, scale: cfg.NORM_SCALE, maxPixels: 1e13, tileScale: 4
  }));
  var groups = ee.List(res.get('groups'));
  var byClass = ee.Dictionary.fromLists(
    groups.map(function (g) { return ee.Number(ee.Dictionary(g).get('cls')).int().format(); }),
    groups.map(function (g) { return ee.Dictionary(g).get('mean'); }));
  return ORDER.map(function (i) {
    var m = ee.List(byClass.get(String(i), [1, 1, 1]));
    return ORDER.map(function (j) { return ee.Number(m.get(j)); });
  });
}

function caStep(M, A, seed, geom, r, proj) {
  var par = caParams(r);
  var N = neighbourShares(M, proj, par);
  var norms = normalisers(M, N, geom, r);
  // Spatially correlated draws (CA_CORR_SIGMA, calibrated on the hindcast) so
  // that transitions occur as small coherent groups of pixels rather than as
  // isolated pixels, which would inflate patch and edge density.
  var u = correlatedUniform(seed, proj, par.sigma);
  var out = M;
  ORDER.forEach(function (i) {
    var js = ORDER.filter(function (j) { return j !== i; });
    var p = js.map(function (j) {
      var rate = ee.Number(A.get([i, j])).divide(norms[i][j]);
      return N.select('n' + j).multiply(rate);
    });
    var shrink = ee.Image(0.99).divide(p[0].add(p[1]).max(0.99));
    var p0 = p[0].multiply(shrink), p1 = p[1].multiply(shrink);
    var isI = M.eq(i);
    out = out.where(isI.and(u.lt(p0)), js[0])
             .where(isI.and(u.gte(p0)).and(u.lt(p0.add(p1))), js[1]);
  });
  return out.reproject(proj);
}

// ---------------------------------------------------------------------------
// Classification error propagation
// ---------------------------------------------------------------------------
// Realisations of the observed map, drawn from the Random Forest class
// probabilities (vote shares). The uniform draws are
// spatially correlated (Gaussian-smoothed white noise, standardised and passed
// through the normal CDF), so uncertain areas change as coherent groups of
// pixels rather than as salt-and-pepper noise that would inflate PD and ED.
// Pixels classified with near-unanimous votes almost never change.
// Spatially correlated uniform field: white noise smoothed with a Gaussian of
// sigma pixels, standardised with the exact SD of the smoothed noise and passed
// through the normal CDF, so values stay uniform on [0, 1] but neighbouring
// pixels receive similar draws. sigma = 0 gives independent draws.
function gaussKernel(s) {
  var rad = Math.ceil(3 * s), w = [], tot = 0, sq = 0;
  for (var i = -rad; i <= rad; i++) {
    for (var j = -rad; j <= rad; j++) {
      var v = Math.exp(-(i * i + j * j) / (2 * s * s)); w.push(v); tot += v;
    }
  }
  w.forEach(function (v) { sq += (v / tot) * (v / tot); });
  return {kernel: ee.Kernel.gaussian({radius: rad, sigma: s, units: 'pixels', normalize: true}),
          sd: Math.sqrt(sq / 12)};
}
var KERNELS = {};
function correlatedUniform(seed, proj, sigma) {
  var u = ee.Image.random(seed).reproject(proj);
  if (!sigma || sigma <= 0) return u;
  var k = KERNELS[sigma] || (KERNELS[sigma] = gaussKernel(sigma));
  var z = u.convolve(k.kernel).reproject(proj).subtract(0.5).divide(k.sd);
  return z.divide(Math.SQRT2).erf().add(1).multiply(0.5);
}

function drawStart(full, seed, proj) {
  var hard = full.select('class');
  if (!cfg.PROPAGATE_CLASSIFICATION_ERROR) return hard;
  var tot = full.select(['p_no', 'p_low', 'p_high']).reduce(ee.Reducer.sum()).max(1);
  var pNo = full.select('p_no').divide(tot);
  var pLow = full.select('p_low').divide(tot);
  var u = correlatedUniform(seed, proj, cfg.ERROR_CORR_SIGMA);
  return ee.Image(2).where(u.lt(pNo.add(pLow)), 1).where(u.lt(pNo), 0)
    .updateMask(hard.mask()).rename('class').toByte().reproject(proj);
}

// Returns {iteration: image} for the iterations listed in keep
function simulate(M0, A, nIter, keep, seedBase, geom, r, proj) {
  var M = M0, kept = {};
  for (var y = 1; y <= nIter; y++) {
    M = caStep(M, A, seedBase + y, geom, r, proj);
    if (keep.indexOf(y) >= 0) kept[y] = M;
  }
  return kept;
}

function ensembleBands(runs, prefix) {
  var col = ee.ImageCollection(runs);
  var names = runs.map(function (img, k) { return prefix + '_r' + pad2(k + 1); });
  var stack = ee.Image.cat(runs).rename(names);
  var probs = ee.Image.cat(ORDER.map(function (j) {
    return col.map(function (img) { return img.eq(j); }).mean().multiply(100).round()
      .rename(prefix + '_p_' + cfg.CLASS_SHORT[j]);
  }));
  var mode = col.mode().rename(prefix + '_mode');
  var agree = col.map(function (img) { return img.eq(mode); }).mean().multiply(100).round()
    .rename(prefix + '_agree');
  return {image: stack.addBands(probs).addBands(mode).addBands(agree), names: names, mode: mode};
}

// ---------------------------------------------------------------------------
// Hindcast evaluation
// ---------------------------------------------------------------------------
function evaluate(sim, base, obs, refFc, geom, r, label) {
  var sampled = sim.rename('sim').sampleRegions({
    collection: refFc, properties: ['class'], scale: cfg.SCALE,
    projection: r.utm, tileScale: 4
  });
  var cm = sampled.errorMatrix('class', 'sim', ORDER);
  var obsCh = obs.neq(base), simCh = sim.neq(base);
  var parts = ee.Image.cat([
    obsCh.and(simCh).and(sim.eq(obs)).rename('hits'),
    obsCh.and(simCh).and(sim.neq(obs)).rename('wrong_hits'),
    obsCh.and(simCh.not()).rename('misses'),
    obsCh.not().and(simCh).rename('false_alarms'),
    sim.eq(obs).rename('agree'),
    ee.Image(1).updateMask(sim.mask().and(obs.mask())).rename('n')
  ]);
  var s = ee.Dictionary(parts.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: geom, crs: r.utm, scale: cfg.SCALE,
    maxPixels: 1e13, tileScale: 4
  }));
  var h = ee.Number(s.get('hits')), w = ee.Number(s.get('wrong_hits'));
  var m = ee.Number(s.get('misses')), f = ee.Number(s.get('false_alarms'));
  var fom = h.divide(h.add(w).add(m).add(f).max(1));
  return lib.accuracyFeature(cm, {reserve: r.id, biome: r.biome, run: label})
    .set({OA_pixel: ee.Number(s.get('agree')).divide(ee.Number(s.get('n')).max(1)),
          FoM: fom, hits: h, wrong_hits: w, misses: m, false_alarms: f});
}

// ---------------------------------------------------------------------------
// Main loop
// ---------------------------------------------------------------------------
var reserves = lib.reservesToRun(cfg);
var matRows = [];
var ACC_COLS = ['reserve', 'biome', 'run', 'n', 'OA', 'OA_ci95', 'kappa',
  'PA_high', 'PA_low', 'PA_no', 'UA_high', 'UA_low', 'UA_no', 'F1_high', 'F1_low', 'F1_no',
  'OA_pixel', 'FoM', 'hits', 'wrong_hits', 'misses', 'false_alarms'];

var bandIndex = null;

var fwdIter = cfg.TARGET_YEARS.map(function (y) { return y - cfg.T1; });
var nFwd = Math.max.apply(null, fwdIter);

reserves.forEach(function (r, ri) {
  var accR = [], qtyR = [];
  var geom = lib.reserveGeometry(cfg, r);
  var proj = lib.gridProjection(cfg, r);
  var c0 = lulc(r, cfg.T0), c1 = lulc(r, cfg.T1);

  var counts = countMatrix(c0, c1, geom, r);
  var P = rowNormalise(counts);
  var am = annualMatrix(P);
  var A = am.A;

  var hindRuns = [], fwdRuns = {}, start0 = [], start1 = [];
  fwdIter.forEach(function (n) { fwdRuns[n] = []; });
  for (var k = 1; k <= cfg.N_RUNS; k++) {
    var seedH = cfg.SEED * 100000 + ri * 10000 + k * 100;
    var seedF = seedH + 5000;
    var s0 = drawStart(lulcFull(r, cfg.T0), seedH, proj);
    var s1 = drawStart(lulcFull(r, cfg.T1), seedF, proj);
    start0.push(s0); start1.push(s1);
    // CA_START 'classified': every run starts from the classified map, so the
    // ensemble spread is simulation uncertainty and the calibrated CA keeps the
    // observed landscape configuration. 'realisation': runs start from s0/s1.
    var fromReal = cfg.CA_START === 'realisation';
    var h = simulate(fromReal ? s0 : c0, A, INTERVAL, [INTERVAL], seedH, geom, r, proj);
    hindRuns.push(h[INTERVAL].toByte());
    var f = simulate(fromReal ? s1 : c1, A, nFwd, fwdIter, seedF, geom, r, proj);
    fwdIter.forEach(function (n) { fwdRuns[n].push(f[n].toByte()); });
  }

  // Ensemble images
  var hindEns = ensembleBands(hindRuns, 'h' + cfg.T1);
  var img = hindEns.image;
  var names = [].concat(hindEns.names, ['p_no', 'p_low', 'p_high', 'mode', 'agree']
    .map(function (s) { return 'h' + cfg.T1 + '_' + s; }));
  if (cfg.PROPAGATE_CLASSIFICATION_ERROR) {
    // Realisations of the observed maps, used in R for confidence intervals of
    // the observed 2019 and 2024 landscape metrics.
    [[start0, 'i' + cfg.T0], [start1, 'i' + cfg.T1]].forEach(function (x) {
      var nm = x[0].map(function (im, k) { return x[1] + '_r' + pad2(k + 1); });
      img = img.addBands(ee.Image.cat(x[0]).rename(nm));
      names = names.concat(nm);
    });
  }
  cfg.TARGET_YEARS.forEach(function (y, yi) {
    var e = ensembleBands(fwdRuns[fwdIter[yi]], 'y' + y);
    img = img.addBands(e.image);
    names = names.concat(e.names, ['p_no', 'p_low', 'p_high', 'mode', 'agree']
      .map(function (s) { return 'y' + y + '_' + s; }));

    // Quantity check: Markov-expected vs ensemble-mean class areas (km2)
    var startCounts = counts.reduce(ee.Reducer.sum(), [0]);              // 1 x 3, classes at T1
    var expected = startCounts.matrixMultiply(matPow(A, fwdIter[yi]));    // 1 x 3
    var sim = e.image.select(['y' + y + '_p_no', 'y' + y + '_p_low', 'y' + y + '_p_high'])
      .divide(100).reduceRegion({reducer: ee.Reducer.sum(), geometry: geom,
        crs: r.utm, scale: cfg.SCALE, maxPixels: 1e13, tileScale: 4});
    ORDER.forEach(function (j) {
      qtyR.push(ee.Feature(null, {
        reserve: r.id, biome: r.biome, year: y, class: j,
        markov_expected_km2: ee.Number(expected.get([0, j])).multiply(cfg.SCALE * cfg.SCALE / 1e6),
        simulated_mean_km2: ee.Number(sim.get('y' + y + '_p_' + cfg.CLASS_SHORT[j]))
          .multiply(cfg.SCALE * cfg.SCALE / 1e6)
      }));
    });
  });
  if (bandIndex === null) bandIndex = names;

  // Hindcast accuracy, per run, modal map and persistence null model
  var ref = lib.referencePoints(cfg, r, cfg.T1);
  hindRuns.forEach(function (sim, k) {
    accR.push(evaluate(sim, c0, c1, ref, geom, r, 'run_' + pad2(k + 1)));
  });
  accR.push(evaluate(hindEns.mode, c0, c1, ref, geom, r, 'ensemble_mode'));
  accR.push(evaluate(c0, c0, c1, ref, geom, r, 'persistence_null'));

  // Tables are exported per reserve: a single table for all reserves and runs
  // exceeds the request size limit.
  Export.table.toDrive({collection: ee.FeatureCollection(accR),
    description: 'cam_hindcast_accuracy_' + r.id, folder: cfg.DRIVE_FOLDER,
    fileNamePrefix: 'cam_hindcast_accuracy_' + r.id, fileFormat: 'CSV', selectors: ACC_COLS});
  Export.table.toDrive({collection: ee.FeatureCollection(qtyR),
    description: 'cam_quantity_check_' + r.id, folder: cfg.DRIVE_FOLDER,
    fileNamePrefix: 'cam_quantity_check_' + r.id, fileFormat: 'CSV',
    selectors: ['reserve', 'biome', 'year', 'class', 'markov_expected_km2', 'simulated_mean_km2']});

  // Matrices
  var mat = {reserve: r.id, biome: r.biome, root_method: am.method, root_error: am.err};
  ORDER.forEach(function (i) {
    ORDER.forEach(function (j) {
      mat['P_' + i + j] = P.get([i, j]);
      mat['A_' + i + j] = A.get([i, j]);
      mat['n_' + i + j] = counts.get([i, j]);
    });
  });
  matRows.push(ee.Feature(null, mat));

  if (cfg.EXPORT_CAM_IMAGES) Export.image.toDrive({
    image: img.clip(geom).unmask(255, false).toByte(),  // clip first so outside pixels are written as nodata, not 0
    description: 'cam_' + r.id,
    folder: cfg.DRIVE_FOLDER,
    fileNamePrefix: 'cam_' + r.id,
    region: geom,
    crs: r.utm,
    crsTransform: [cfg.SCALE, 0, 0, 0, -cfg.SCALE, 0],
    maxPixels: 1e13,
    formatOptions: {cloudOptimized: true, noData: 255}
  });

  if (ri === 0) {
    Map.centerObject(geom, 12);
    Map.addLayer(c1, {min: 0, max: 2, palette: cfg.CLASS_PALETTE}, r.id + ' observed ' + cfg.T1);
    var last = cfg.TARGET_YEARS[cfg.TARGET_YEARS.length - 1];
    Map.addLayer(img.select('y' + last + '_mode'), {min: 0, max: 2, palette: cfg.CLASS_PALETTE},
      r.id + ' ' + last + ' modal projection');
    Map.addLayer(img.select('y' + last + '_agree'), {min: 50, max: 100, palette: ['d7191c', 'ffffbf', '2c7bb6']},
      r.id + ' ' + last + ' run agreement (%)');
    print(r.id + ' P (5-year)', P, r.id + ' A (annual)', A, 'root method', am.method, 'root error', am.err);
  }
});

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------

var matCols = ['reserve', 'biome', 'root_method', 'root_error'];
['n', 'P', 'A'].forEach(function (p) {
  ORDER.forEach(function (i) { ORDER.forEach(function (j) { matCols.push(p + '_' + i + j); }); });
});
Export.table.toDrive({collection: ee.FeatureCollection(matRows), description: 'cam_markov_matrices_csv',
  folder: cfg.DRIVE_FOLDER, fileNamePrefix: 'cam_markov_matrices', fileFormat: 'CSV', selectors: matCols});


Export.table.toDrive({
  collection: ee.FeatureCollection(bandIndex.map(function (b, i) {
    return ee.Feature(null, {band_number: i + 1, band_name: b});
  })),
  description: 'cam_band_index_csv', folder: cfg.DRIVE_FOLDER,
  fileNamePrefix: 'cam_band_index', fileFormat: 'CSV', selectors: ['band_number', 'band_name']
});
