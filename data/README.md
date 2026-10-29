# Data

## reserves.csv

| Column | Description |
|---|---|
| reserve | identifier used throughout the workflow |
| biome | Rainforest, Mangrove/swamp or Savannah |
| wdpaid | World Database on Protected Areas identifier |
| wdpa_name | WDPA name |

Boundaries were taken from the WDPA (UNEP-WCMC and IUCN, December 2025 release). The WDPA terms of use do not allow redistribution, so the polygons are not included. They can be rebuilt from the WDPA IDs, or loaded in GEE directly from `WCMC/WDPA/current/polygons` by setting `RESERVE_SOURCE = 'wdpa'` in `gee/config.js` (the current release may differ slightly from the one used in the paper).

## reference/original_reference_points.csv

Photo-interpreted reference points. They are used for accuracy assessment by five-fold cross-validation and are also added to the training data (each point repeated three times). Labels were assigned by visual interpretation of high-resolution Google Earth imagery dated as close as possible to each Sentinel-2 composite. They are reference data, not field observations.

| Column | Description |
|---|---|
| ref_id | unique identifier, `REF_<reserve>_<year>_<nnn>` |
| reserve | reserve identifier |
| year | epoch (2019 = Nov 2018 to Mar 2019 composite; 2024 = Nov 2023 to Mar 2024) |
| class | 0 no vegetation, 1 low vegetation, 2 high vegetation |
| class_name | class label (corrected, see below) |
| name_original | label as stored in the source file |
| longitude, latitude | WGS 84 (EPSG:4326) |
| source_file | original QGIS spectral library file |

Cleaning applied by `scripts/prepare_inputs.py` to the original 1,092 points: three points without a label were removed, one exact duplicate was removed, and label spelling was harmonised (`High Vegetation` and `High vegetation` were both used in the source files). 1,088 points remain (59 to 78 per reserve and epoch).

The High and Low vegetation labels in the source files are reversed and are swapped back by `prepare_inputs.py` (`SWAP_HIGH_LOW = True`). Points named "Low vegetation" lie on closed canopy and points named "High vegetation" on open farmland and fallow mosaic. This was confirmed by visual checks over high-resolution imagery and by comparison with Dynamic World, ESA WorldCover, ETH canopy height and Sentinel-2 NDVI at all points. Agreement with Dynamic World was 20 % for the labels as named and 70 % after the swap (934 points with complete data). The original name is retained in `name_original`. `original_reference_points_summary.csv` gives the counts per class.

## Outputs archived on Zenodo (not in the GitHub repository)

The classified maps with class probabilities, change maps, CA-Markov ensembles, consensus training samples, all tables and the figures are archived on Zenodo (https://doi.org/10.5281/zenodo.[DATA DOI]) because of their size.

Licence for the data in this folder: CC BY 4.0.
