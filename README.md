# Forest fragmentation in eight Nigerian forest reserves (Sentinel-2, Random Forest, CA–Markov)

Code and reference data for the article

> Akano, H., Igbaoreto, W., Odebiri, O., Young, M., Olamide, R., Adesina, A., Naicker, R., Sodiq, A., [year]. Forecasting forest fragmentation across Nigeria's major biomes using Sentinel-2, Random Forest, and CA–Markov. *Scientific African* [volume, article number]. https://doi.org/[article DOI]

The workflow maps high, low and no vegetation in eight forest reserves for the 2019 and 2024 dry seasons, measures change in landscape configuration, propagates classification uncertainty into the landscape metrics, and projects cover and configuration to 2028 and 2033 with a calibrated stochastic CA–Markov model.

The maps, tables and figures produced by this code are archived on Zenodo: https://doi.org/10.5281/zenodo.[DATA DOI]

## Classes

| Value | Class | Definition |
|---|---|---|
| 2 | High vegetation | canopy cover ≥ 60 % |
| 1 | Low vegetation | canopy cover ≥ 20 % and < 60 % |
| 0 | No vegetation | canopy cover < 20 % (bare soil, water, built-up) |

## Reserves

| Reserve | Biome group | State | WDPA ID |
|---|---|---|---|
| Stubbs Creek | Mangrove/swamp | Akwa Ibom | 7851 |
| Upper Orashi | Mangrove/swamp | Rivers | 903112 |
| Akanto | Rainforest | Ebonyi | 37026 |
| Gambari | Rainforest | Oyo | 36817 |
| Adoru * | Savannah | Kogi | 36860 |
| Maradun * | Savannah | Zamfara | 36323 |
| Mayo-Ndaga | Savannah | Taraba | 36911 |
| Yo | Savannah | Borno | 36210 |

\* Change between the two dry-season composites in these reserves is dominated by differences in burning and grass senescence. They are reported individually and excluded from aggregate results.

## Repository layout

```
gee/                              Google Earth Engine scripts (JavaScript)
  config.js                       all parameters (set ASSET_ROOT first)
  lib/common.js                   shared functions
  00_reference_diagnostics.js     checks of the reference points against global products
  00b_composite_diagnostics.js    burned area, observation counts and NDVI per composite
  01_build_training_samples.js    consensus training samples, spatial blocks
  02_classify_and_assess.js       Random Forest maps, class probabilities, accuracy, importance
  02b_export_lulc_maps.js         GeoTIFFs of the maps and class probabilities
  03_change_analysis.js           2019-2024 transitions and change GeoTIFFs
  04_ca_markov_simulation.js      CA-Markov ensemble, hindcast, uncertainty layers
R/
  00_setup.R                      packages, paths and helpers
  02_landscape_metrics.R          PD, ED, LPI and AI for every map, realisation and run
  03_tables_and_figures.R         Tables 2-4 and Figures 2, 5, 6, 7, S1
  04_map_plates.R                 Figures 3a-c and 4
  05_supplementary_tables.R       Tables S1-S8
  06_export_figure_maps.R         the rasters shown in Figures 3, 4 and 7 as styled GeoTIFFs
  07_package_maps.R               masked copies of the GEE exports for the data archive
  calibration/                    calibration of the CA rule and of the probability draws
scripts/
  prepare_inputs.py               builds the upload files from the original QGIS spectral libraries
  make_standalone.py              GEE scripts with config and lib pasted in (no require())
data/
  reserves.csv                    reserve, biome group and WDPA ID
  reference/                      1,088 photo-interpreted reference points
tests/make_synthetic_outputs.R    synthetic GEE outputs for testing the R scripts
```

## Running the workflow

### 1. Google Earth Engine

1. Upload the inputs to a folder in your Cloud project:
   - the reserve boundaries as `ASSET_ROOT/reserves` (shapefile with a `reserve` field matching `data/reserves.csv`; see Data below);
   - `data/reference/original_reference_points.csv` as `ASSET_ROOT/original_reference_points` (X column `longitude`, Y column `latitude`, EPSG:4326).
2. Create an empty image folder `ASSET_ROOT/lulc`.
3. Set `ASSET_ROOT` in `gee/config.js`.
4. Either copy `gee/` into a Code Editor script repository and set `REPO` at the top of each script, or build self-contained copies and paste each into a new script:
   ```bash
   python scripts/make_standalone.py --asset-root projects/<cloud-project>/assets/nigeria_fragmentation
   ```
5. Run the scripts in order (01, 02, 02b, 03, 04) and start all export tasks. Wait for the assets of 02 before running the later scripts. For a first test set `RUN_ONLY = ['Akanto']` in `config.js`.

### 2. R

Download the Drive folder `nigeria_fragmentation_gee` into `outputs/gee/` and run from the repository root:

```bash
Rscript R/02_landscape_metrics.R
Rscript R/03_tables_and_figures.R
Rscript R/04_map_plates.R
Rscript R/05_supplementary_tables.R
Rscript R/06_export_figure_maps.R
```

`02_landscape_metrics.R` caches each reserve and can run as parallel workers (see the script header). GeoTIFFs exported from Earth Engine code pixels outside each reserve as 0, which is also the No vegetation class, so the R scripts mask every layer with the reserve boundaries in `gee_upload/reserves.zip` when that file is present. The masked maps in the Zenodo archive can be used instead.

Tested with R 4.3, terra 1.7, landscapemetrics 2.2.1, sf 1.0, ggplot2 3.4 and patchwork. `Rscript tests/make_synthetic_outputs.R` writes random data into `outputs/gee/` for testing the R scripts without Earth Engine.

### 3. Calibration (optional)

The CA parameters (`CA_PARAMS`) and the correlation length of the probability draws (`ERROR_CORR_SIGMA`) in `gee/config.js` were set with `R/calibration/calibrate_correlation.R`. The script runs the same CA rule in R from the classified 2019 map and selects, for each reserve, the parameters whose simulated 2024 map best reproduces the observed 2024 patch and edge density.

## Main settings

| Step | Setting |
|---|---|
| Composites | Sentinel-2 SR, 1 Nov to 31 Mar (2018/19, 2023/24), Cloud Score+ cs_cdf ≥ 0.60, median |
| Predictors | B2, B3, B4, B8, B11, B12, NDVI, NDWI, NBR, SAVI, NDVI p90, NBR p90, NDVI p10 |
| Training data | 355–500 consensus samples per reserve and year (Dynamic World, ESA WorldCover, Hansen GFC) plus reference points (weight 3) |
| Validation | five-fold cross-validation on reference points; 1 km blocked hold-out; consensus-only model |
| Random Forest | 300 trees, seed 101 |
| Uncertainty | 20 realisations of each map drawn from the class probabilities (spatially correlated draws, σ = 8 pixels) |
| CA–Markov | annual matrix = 5th root of the 2019–2024 matrix (linear approximation for Yo); 3 × 3 neighbourhood; per-reserve σ, k and e calibrated on the hindcast; 20 runs from the classified 2024 map |
| Landscape metrics | PD, ED, LPI, AI; landscape level; eight-cell rule; local UTM zone; 10 m |

## Data

`data/reference/original_reference_points.csv` holds the photo-interpreted reference points. The high- and low-vegetation labels in the original source files were interchanged and are corrected by `scripts/prepare_inputs.py`; the original label is kept in the column `name_original`. See `data/README.md`.

Reserve boundaries come from the World Database on Protected Areas (UNEP-WCMC and IUCN) and are not redistributed. They can be obtained from https://www.protectedplanet.net using the WDPA IDs above, or loaded in Earth Engine from `WCMC/WDPA/current/polygons` by setting `RESERVE_SOURCE = 'wdpa'` in `gee/config.js`. The current WDPA release may differ slightly from the one used in the study.

## Licence

Code: MIT (see `LICENSE`). Data in `data/`: CC BY 4.0.

## Citation

Please cite the article and the archived release (see `CITATION.cff`).
