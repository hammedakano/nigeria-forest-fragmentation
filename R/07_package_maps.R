# 07_package_maps.R
# Writes masked, compressed copies of the GEE map exports for the data archive.
# Pixels outside each reserve are set to nodata (255) and every band keeps its
# name (GeoTIFF band description).
#
#   zenodo_maps/classification/lulc_<reserve>.tif   class and RF class probabilities, 2019 and 2024
#   zenodo_maps/change/obs_<reserve>.tif            classified maps, transition code, p_high
#   zenodo_maps/ca_markov/cam_<reserve>.tif         hindcast and projection ensembles (bands in cam_band_index.csv)
#
# Needs gee_upload/reserves.zip (scripts/prepare_inputs.py).

source("R/00_setup.R")
OUT <- "zenodo_maps"
stopifnot(file.exists(RESERVES_ZIP))
band_index <- read_csv(file.path(DIR_GEE, "cam_band_index.csv"), show_col_types = FALSE) |> arrange(band_number)
LULC_BANDS <- as.vector(outer(c("class", "p_no", "p_low", "p_high", "n_obs"), c(T0, T1), paste, sep = "_"))
OBS_BANDS <- c(paste0("class_", T0), paste0("class_", T1), "transition", paste0("p_high_", T0), paste0("p_high_", T1))
sets <- list(classification = list(prefix = "lulc_", bands = LULC_BANDS),
             change = list(prefix = "obs_", bands = OBS_BANDS),
             ca_markov = list(prefix = "cam_", bands = band_index$band_name))
for (res in RESERVES$reserve) for (s in names(sets)) {
  f <- file.path(DIR_GEE, paste0(sets[[s]]$prefix, res, ".tif"))
  r <- read_gee_tif(f, sets[[s]]$bands, res)
  dir.create(file.path(OUT, s), recursive = TRUE, showWarnings = FALSE)
  writeRaster(r, file.path(OUT, s, basename(f)), datatype = "INT1U", NAflag = NODATA, overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))
  message(s, " ", res)
}
file.copy(file.path(DIR_GEE, "cam_band_index.csv"), file.path(OUT, "ca_markov", "cam_band_index.csv"), overwrite = TRUE)
