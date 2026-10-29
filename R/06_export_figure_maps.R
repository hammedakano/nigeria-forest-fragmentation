# 06_export_figure_maps.R
# Writes the rasters shown in the map figures as separate single-band GeoTIFFs,
# masked to each reserve (nodata 255), with QGIS styles, so the figures can be
# checked or redrawn in a GIS.
#
#   figure_maps/
#     Figure3_maps/<Biome>/<Reserve>_<year>_classified.tif   (2019, 2024)
#     Figure3_maps/<Biome>/<Reserve>_<year>_projected.tif    (2028, 2033; modal class of 20 runs)
#     Figure4_hotspots/<Reserve>_<year>_window.tif           (1.5 km windows) + windows.gpkg
#     Figure7_uncertainty/<Reserve>_2033_run_agreement.tif   (% of runs = modal class)
#     styles/classes.qml, styles/agreement.qml
#
# Classes: 0 no vegetation, 1 low vegetation, 2 high vegetation.

source("R/00_setup.R")
suppressPackageStartupMessages(library(sf))

OUT <- "figure_maps"
OBS_BANDS <- c(paste0("class_", T0), paste0("class_", T1), "transition",
               paste0("p_high_", T0), paste0("p_high_", T1))
band_index <- read_csv(file.path(DIR_GEE, "cam_band_index.csv"), show_col_types = FALSE) |> arrange(band_number)
biome_dir <- c("Rainforest" = "a_Rainforest", "Mangrove/swamp" = "b_Mangrove_swamp", "Savannah" = "c_Savannah")
wopt <- c("COMPRESS=LZW", "TILED=YES")

save_tif <- function(r, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeRaster(r, path, datatype = "INT1U", NAflag = NODATA, overwrite = TRUE, gdal = wopt)
}

for (i in seq_len(nrow(RESERVES))) {
  res <- RESERVES$reserve[i]; biome <- RESERVES$biome[i]
  message("Exporting ", res)
  obs <- read_gee_tif(file.path(DIR_GEE, paste0("obs_", res, ".tif")), OBS_BANDS, res)
  cam <- read_gee_tif(file.path(DIR_GEE, paste0("cam_", res, ".tif")), band_index$band_name, res)
  d3 <- file.path(OUT, "Figure3_maps", biome_dir[[biome]])
  for (y in c(T0, T1)) save_tif(obs[[paste0("class_", y)]], file.path(d3, paste0(res, "_", y, "_classified.tif")))
  for (y in TARGET_YEARS) save_tif(cam[[paste0("y", y, "_mode")]], file.path(d3, paste0(res, "_", y, "_projected.tif")))
  save_tif(cam[[paste0("y", TARGET_YEARS[2], "_agree")]],
           file.path(OUT, "Figure7_uncertainty", paste0(res, "_", TARGET_YEARS[2], "_run_agreement.tif")))
}

# Figure 4 hotspot windows (same windows as R/04_map_plates.R)
win <- read_csv(file.path(DIR_TAB, "Figure4_hotspot_windows.csv"), show_col_types = FALSE)
polys <- list()
for (k in seq_len(nrow(win))) {
  res <- win$reserve[k]
  obs <- read_gee_tif(file.path(DIR_GEE, paste0("obs_", res, ".tif")), OBS_BANDS, res)
  cam <- read_gee_tif(file.path(DIR_GEE, paste0("cam_", res, ".tif")), band_index$band_name, res)
  ctr <- st_transform(st_sfc(st_point(c(win$centre_lon[k], win$centre_lat[k])), crs = 4326), crs(obs))
  xy <- st_coordinates(ctr); h <- win$window_m[k] / 2
  e <- ext(xy[1] - h, xy[1] + h, xy[2] - h, xy[2] + h)
  layers <- list(obs[[paste0("class_", T0)]], obs[[paste0("class_", T1)]], cam[[paste0("y", TARGET_YEARS[2], "_mode")]])
  for (j in seq_along(layers)) {
    y <- c(T0, T1, TARGET_YEARS[2])[j]
    save_tif(crop(layers[[j]], e, snap = "near"), file.path(OUT, "Figure4_hotspots", paste0(res, "_", y, "_window.tif")))
  }
  polys[[res]] <- st_sf(reserve = res, geometry = st_transform(st_as_sfc(st_bbox(c(xmin = xy[1] - h, ymin = xy[2] - h,
                         xmax = xy[1] + h, ymax = xy[2] + h), crs = st_crs(crs(obs)))), 4326))
}
st_write(do.call(rbind, polys), file.path(OUT, "Figure4_hotspots", "hotspot_windows.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# QGIS styles
dir.create(file.path(OUT, "styles"), showWarnings = FALSE)
writeLines(c(
'<!DOCTYPE qgis PUBLIC "http://mrcc.com/qgis.dtd" "SYSTEM">',
'<qgis version="3.28" styleCategories="Symbology">',
' <pipe><rasterrenderer type="paletted" band="1" opacity="1" nodataColor="">',
'  <colorPalette>',
'   <paletteEntry value="0" color="#d9d9d9" alpha="255" label="No vegetation"/>',
'   <paletteEntry value="1" color="#a6d96a" alpha="255" label="Low vegetation"/>',
'   <paletteEntry value="2" color="#1a9850" alpha="255" label="High vegetation"/>',
'  </colorPalette>',
' </rasterrenderer></pipe>',
'</qgis>'), file.path(OUT, "styles", "classes.qml"))
writeLines(c(
'<!DOCTYPE qgis PUBLIC "http://mrcc.com/qgis.dtd" "SYSTEM">',
'<qgis version="3.28" styleCategories="Symbology">',
' <pipe><rasterrenderer type="singlebandpseudocolor" band="1" opacity="1" classificationMin="30" classificationMax="100">',
'  <rastershader><colorrampshader colorRampType="INTERPOLATED" clip="0" minimumValue="30" maximumValue="100">',
'   <item value="30" color="#d7191c" alpha="255" label="30"/>',
'   <item value="47.5" color="#fdae61" alpha="255" label="47.5"/>',
'   <item value="65" color="#ffffbf" alpha="255" label="65"/>',
'   <item value="82.5" color="#abd9e9" alpha="255" label="82.5"/>',
'   <item value="100" color="#2c7bb6" alpha="255" label="100"/>',
'  </colorrampshader></rastershader>',
' </rasterrenderer></pipe>',
'</qgis>'), file.path(OUT, "styles", "agreement.qml"))

# Copy the classes style next to every class raster so QGIS applies it automatically
for (f in list.files(OUT, pattern = "(classified|projected|window)\\.tif$", recursive = TRUE, full.names = TRUE))
  file.copy(file.path(OUT, "styles", "classes.qml"), sub("\\.tif$", ".qml", f), overwrite = TRUE)
for (f in list.files(OUT, pattern = "agreement\\.tif$", recursive = TRUE, full.names = TRUE))
  file.copy(file.path(OUT, "styles", "agreement.qml"), sub("\\.tif$", ".qml", f), overwrite = TRUE)

writeLines(c(
"# Maps used in the figures",
"",
"Single-band GeoTIFFs, 10 m, local UTM zone of each reserve (Gambari EPSG:32631, Yo EPSG:32633, others EPSG:32632),",
"masked to the reserve boundary, nodata = 255. Each raster has a QGIS style (.qml) with the same name.",
"",
"- Figure3_maps/: classified maps (2019, 2024) and projected maps (2028, 2033), one folder per biome plate (a, b, c).",
"  Values: 0 no vegetation, 1 low vegetation, 2 high vegetation. Projected maps are the modal class of 20 CA-Markov runs.",
"- Figure4_hotspots/: the 1.5 km x 1.5 km windows in Stubbs Creek and Akanto (2019, 2024, 2033) and hotspot_windows.gpkg",
"  with the window outlines (WGS 84).",
"- Figure7_uncertainty/: percentage of the 20 runs that assign the modal class to each pixel in 2033 (0-100).",
"- styles/: the two QGIS styles.",
"",
"Made by R/06_export_figure_maps.R from the GEE exports in outputs/gee/."), file.path(OUT, "README.md"))
message("Done: ", OUT)
