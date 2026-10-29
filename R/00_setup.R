# 00_setup.R
# Shared settings for the R part of the workflow. Sourced by every other script.
#
# Expected layout (relative to the repository root):
#   outputs/gee/      everything exported by the GEE scripts to Google Drive
#                     (download the Drive folder 'nigeria_fragmentation_gee' here)
#   outputs/tables/   CSV tables written by the R scripts
#   outputs/figures/  figures written by the R scripts
#
# Run every script from the repository root, e.g.  Rscript R/02_landscape_metrics.R

required <- c("terra", "landscapemetrics", "dplyr", "tidyr", "readr", "ggplot2",
              "patchwork", "ggspatial")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(terra)
  library(landscapemetrics)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

DIR_GEE  <- file.path("outputs", "gee")
DIR_TAB  <- file.path("outputs", "tables")
DIR_FIG  <- file.path("outputs", "figures")
dir.create(DIR_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

RESERVES <- read_csv(file.path("data", "reserves.csv"), show_col_types = FALSE)
RESERVE_LABEL <- c(Adoru = "Adoru", Akanto = "Akanto", Gambari = "Gambari",
                   Maradun = "Maradun", Mayo_Ndaga = "Mayo-Ndaga",
                   Stubbs_Creek = "Stubbs Creek", Upper_Orashi = "Upper Orashi", Yo = "Yo")

CLASS_NAMES <- c("0" = "No vegetation", "1" = "Low vegetation", "2" = "High vegetation")
CLASS_COLS  <- c("No vegetation" = "#d9d9d9", "Low vegetation" = "#a6d96a",
                 "High vegetation" = "#1a9850")
BIOME_COLS  <- c("Rainforest" = "#1b9e77", "Mangrove/swamp" = "#7570b3", "Savannah" = "#d95f02")

T0 <- 2019; T1 <- 2024; TARGET_YEARS <- c(2028, 2033)

# Reserves whose observed 2019-2024 change is dominated by differences in
# dry-season state (senescence and burning) between the two composites rather
# than land-cover change (see gee/00b_composite_diagnostics.js). They are
# reported individually but excluded from biome and all-reserve totals.
FLAGGED <- c("Adoru", "Maradun")
NODATA <- 255

# Reserve boundaries, used to mask the GEE GeoTIFFs. The exports written before
# the fix in 03/04 (clip before unmask) code pixels outside the reserve as 0,
# which is also the No-vegetation class, so every layer is masked here.
# The shapefile is the one uploaded to GEE (gee_upload/reserves.zip, built by
# scripts/prepare_inputs.py); it is not redistributed (WDPA licence).
RESERVES_ZIP <- file.path("gee_upload", "reserves.zip")
reserve_polygon <- function(reserve) {
  if (!file.exists(RESERVES_ZIP)) return(NULL)
  v <- vect(paste0("/vsizip/", normalizePath(RESERVES_ZIP), "/reserves.shp"))
  v[v$reserve == reserve, ]
}
if (!file.exists(RESERVES_ZIP))
  warning("Reserve boundaries not found at ", RESERVES_ZIP, ". GeoTIFFs exported before the nodata fix ",
          "will not be masked; use the masked maps from the Zenodo archive or run scripts/prepare_inputs.py.",
          call. = FALSE)

# Read a GEE GeoTIFF, set nodata, attach band names and mask to the reserve
read_gee_tif <- function(path, band_names = NULL, reserve = NULL) {
  r <- rast(path)
  NAflag(r) <- NODATA
  if (!is.null(band_names)) {
    stopifnot(nlyr(r) == length(band_names))
    names(r) <- band_names
  }
  poly <- if (is.null(reserve)) NULL else reserve_polygon(reserve)
  if (!is.null(poly)) {
    inside <- rasterize(project(poly, crs(r)), r[[1]], field = 1)
    r <- mask(r, inside)
  }
  r
}

# Percentage change relative to a baseline
pct_change <- function(new, base) ifelse(base == 0, NA_real_, (new - base) / abs(base) * 100)

theme_paper <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.position = "bottom")
}
