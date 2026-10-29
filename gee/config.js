/**
 * config.js
 * Single source of truth for every parameter used by the GEE scripts.
 * Edit ASSET_ROOT (and only that, in most cases) before running anything.
 *
 * Loaded by every script with:
 *   var cfg = require('users/<YOUR_USERNAME>/<YOUR_REPO>:config');
 */

// ---------------------------------------------------------------------------
// 1. Asset locations (EDIT THIS)
// ---------------------------------------------------------------------------
// Folder in your Cloud project that holds the uploaded inputs and the
// intermediate outputs. Create it once in the Assets tab.
var ASSET_ROOT = 'projects/YOUR_CLOUD_PROJECT/assets/nigeria_fragmentation';
exports.ASSET_ROOT = ASSET_ROOT;

// Uploaded inputs (see README, step 2)
exports.RESERVES_ASSET  = ASSET_ROOT + '/reserves';                    // zipped shapefile from gee_upload/
exports.REFERENCE_ASSET = ASSET_ROOT + '/original_reference_points';   // data/reference/original_reference_points.csv

// Produced by 01_build_training_samples.js
exports.SAMPLES_ASSET = ASSET_ROOT + '/training_samples_consensus';

// Optional. Blind verification results (CSV with sample_id, verified_class).
// Leave as '' until the verification step has been done. When set, verified
// labels replace consensus labels and samples marked -1 (unsure) are dropped.
exports.VERIFIED_ASSET = '';

// Produced by 02_classify_and_assess.js. One image per reserve and epoch:
//   ASSET_ROOT/lulc/<reserve>_<year>
exports.LULC_FOLDER = ASSET_ROOT + '/lulc';

// Google Drive folder for CSV and GeoTIFF exports
exports.DRIVE_FOLDER = 'nigeria_fragmentation_gee';

// Reserve boundaries source: 'asset' uses the uploaded shapefile (exact
// boundaries used in the paper); 'wdpa' pulls the polygons straight from the
// WDPA catalogue by WDPAID (current release, may differ slightly).
exports.RESERVE_SOURCE = 'asset';

// ---------------------------------------------------------------------------
// 2. Study reserves
// ---------------------------------------------------------------------------
// utm: local UTM zone used as the 10 m processing grid (Gambari lies in
// zone 31N, Yo in zone 33N, the rest in zone 32N).
exports.RESERVES = [
  {id: 'Adoru',        wdpaid: 36860,  biome: 'Savannah',       utm: 'EPSG:32632'},
  {id: 'Akanto',       wdpaid: 37026,  biome: 'Rainforest',     utm: 'EPSG:32632'},
  {id: 'Gambari',      wdpaid: 36817,  biome: 'Rainforest',     utm: 'EPSG:32631'},
  {id: 'Maradun',      wdpaid: 36323,  biome: 'Savannah',       utm: 'EPSG:32632'},
  {id: 'Mayo_Ndaga',   wdpaid: 36911,  biome: 'Savannah',       utm: 'EPSG:32632'},
  {id: 'Stubbs_Creek', wdpaid: 7851,   biome: 'Mangrove/swamp', utm: 'EPSG:32632'},
  {id: 'Upper_Orashi', wdpaid: 903112, biome: 'Mangrove/swamp', utm: 'EPSG:32632'},
  {id: 'Yo',           wdpaid: 36210,  biome: 'Savannah',       utm: 'EPSG:32633'}
];

// Run a subset while testing, e.g. ['Akanto']. Empty list = all reserves.
exports.RUN_ONLY = [];

// ---------------------------------------------------------------------------
// 3. Epochs (dry-season composites)
// ---------------------------------------------------------------------------
// end dates are exclusive. 2019 = Nov 2018 to Mar 2019; 2024 = Nov 2023 to Mar 2024.
// worldcover / lossYears / lossSide configure the label agreement rule used in
// 01_build_training_samples.js (see the comments there).
exports.EPOCHS = [
  {year: 2019, start: '2018-11-01', end: '2019-04-01',
   worldcover: 'ESA/WorldCover/v100', lossYears: [19, 20], lossSide: 'before'},
  {year: 2024, start: '2023-11-01', end: '2024-04-01',
   worldcover: 'ESA/WorldCover/v200', lossYears: [21, 24], lossSide: 'after'}
];

// ---------------------------------------------------------------------------
// 4. Classes
// ---------------------------------------------------------------------------
// High vegetation: canopy cover >= 60 %; Low vegetation: >= 20 % and < 60 %;
// No vegetation: < 20 % (bare soil, water, built-up).
exports.CLASS_VALUES  = [0, 1, 2];
exports.CLASS_NAMES   = ['No vegetation', 'Low vegetation', 'High vegetation'];
exports.CLASS_SHORT   = ['no', 'low', 'high'];
exports.CLASS_PALETTE = ['d9d9d9', 'a6d96a', '1a9850'];

// ---------------------------------------------------------------------------
// 5. Sentinel-2 compositing and predictors
// ---------------------------------------------------------------------------
exports.SCALE = 10;
// Cloud masking. 'CSPLUS' (Cloud Score+, recommended) or 'QA60'.
// QA60 is empty for most scenes processed between Jan 2022 and Feb 2024, so it
// does not mask clouds in the 2024 composite. Cloud Score+ is used instead.
exports.CLOUD_MASK = 'CSPLUS';
exports.CS_THRESHOLD = 0.60;          // keep pixels with cs_cdf >= threshold
exports.MAX_SCENE_CLOUD = 60;         // scene-level CLOUDY_PIXEL_PERCENTAGE filter
exports.BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
exports.INDICES = ['NDVI', 'NDWI', 'NBR', 'SAVI'];
// Within-season percentile features added to the ten median-composite predictors.
// In the savanna reserves the median composite depends on when grasses dried and
// burned in each year (2019 composite NDVI in Adoru was 0.44 against 0.56 in 2024,
// while the greenest observation was 0.67 and 0.68). The 90th percentiles keep
// woody cover visible through senescence and fire. Set [] to reproduce the
// ten-predictor model.
exports.SEASONAL_FEATURES = ['NDVI_p90', 'NBR_p90', 'NDVI_p10'];

// ---------------------------------------------------------------------------
// 6. Training samples (01_build_training_samples.js)
// ---------------------------------------------------------------------------
exports.N_SAMPLES = 500;              // per reserve and epoch
exports.CLASS_POINTS = [167, 167, 166]; // equal allocation (no, low, high)
exports.DW_MIN_PROB = 0.60;           // mean Dynamic World probability of the label class
exports.DW_MIN_CONSISTENCY = 0.60;    // share of DW scenes agreeing with the label
exports.DW_MIN_SCENES = 3;            // minimum number of DW scenes in the window
// Fallback tiers, used only for classes that fall short of their quota under the
// strict rule (mainly savanna reserves). Every sample records its tier.
exports.TIER2_MIN_PROB = 0.40;        // probability of the single most likely DW class
exports.TIER2_MIN_CONSISTENCY = 0.50;
exports.TIER3_MIN_DW_PROB = 0.25;     // DW grouped probability for the WorldCover class
exports.MIN_DIST_FROM_REF = 50;       // m; keeps training samples away from validation points
exports.BLOCK_SIZE = 1000;            // m; spatial block size for the train/test split
exports.TEST_FRACTION = 0.30;         // share of blocks held out
exports.N_VERIFY_PER_CLASS = 17;      // ~50 points per reserve-epoch for blind checking (10 %)
// Hansen Global Forest Change, used only to keep recently cleared pixels
// eligible for labelling (see 01_build_training_samples.js).
exports.GFC_ASSET = 'UMD/hansen/global_forest_change_2024_v1_12';

// ---------------------------------------------------------------------------
// 7. Random Forest
// ---------------------------------------------------------------------------
exports.N_TREES = 300;
exports.SEED = 101;
// Photo-interpreted reference points are added to the training data, each
// repeated REF_WEIGHT times; accuracy on them is estimated by N_FOLDS-fold
// cross-validation.
exports.REF_WEIGHT = 3;
exports.N_FOLDS = 5;

// ---------------------------------------------------------------------------
// 8. CA-Markov simulation (04_ca_markov_simulation.js)
// ---------------------------------------------------------------------------
exports.T0 = 2019;
exports.T1 = 2024;
exports.TARGET_YEARS = [2028, 2033];
// Set false to export only the CA-Markov tables (images already exported)
exports.EXPORT_CAM_IMAGES = true;
exports.N_RUNS = 20;                  // Monte Carlo ensemble size
exports.NORM_SCALE = 50;             // m; subsample grid for the per-iteration rate normalisation
// CA neighbourhood rule (see 04_ca_markov_simulation.js). Defaults, used for any
// reserve not listed in CA_PARAMS. With independent draws and NUCLEATION 0.05
// (the first version) the hindcast roughly doubled patch and edge density.
exports.NUCLEATION = 0.001;           // lets new patches appear away from existing ones
exports.NEIGHBOUR_POWER = 1;          // exponent on the neighbour share (2 = stronger edge preference)
// Classification error propagation. When true, N_RUNS realisations of the
// observed 2019 and 2024 maps are drawn from the RF class probabilities and
// exported; R uses them for uncertainty intervals of the observed landscape
// metrics.
exports.PROPAGATE_CLASSIFICATION_ERROR = true;
// Starting map of the CA runs: 'classified' (default) or 'realisation'.
// See the header of 04_ca_markov_simulation.js for why realisations are not used.
exports.CA_START = 'classified';
exports.ERROR_CORR_SIGMA = 8;         // pixels; spatial correlation of the probability draws
                                      // (calibrated in R/calibration/, calibration_error_sigma.csv)
// Spatial correlation of the CA transition draws (pixels). 0 = independent draws.
exports.CA_CORR_SIGMA = 4;
// Per-reserve values calibrated on the 2019->2024 hindcast so that the simulated
// 2024 map reproduces the patch and edge density of the classified 2024 map
// (R/calibration/calibrate_correlation.R, outputs/tables/calibration_ca_params.csv).
exports.CA_PARAMS = {
  Adoru:        {sigma: 4, nucleation: 0.001, power: 2},
  Akanto:       {sigma: 4, nucleation: 0.01,  power: 1},
  Gambari:      {sigma: 2, nucleation: 0.01,  power: 2},
  Maradun:      {sigma: 3, nucleation: 0.01,  power: 1},
  Mayo_Ndaga:   {sigma: 1, nucleation: 0.001, power: 2},
  Stubbs_Creek: {sigma: 3, nucleation: 0.001, power: 1},
  Upper_Orashi: {sigma: 4, nucleation: 0.001, power: 2},
  Yo:           {sigma: 3, nucleation: 0.001, power: 1}
};
exports.ROOT_SERIES_TERMS = 60;       // terms in the binomial series for the annual matrix root
exports.ROOT_TOLERANCE = 0.02;        // max abs error accepted for (annual matrix)^interval vs observed P
