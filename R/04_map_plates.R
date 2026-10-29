# 04_map_plates.R
# Figure 3: classified (2019, 2024) and projected (2028, 2033) maps, one plate
# per biome, with a coordinate grid, a scale bar and a north arrow in every
# panel. Figure 4: zoomed views of fragmentation hotspots in Stubbs Creek and
# Akanto at about 1:25,000.
#
# Projected maps show the modal class of the 20 CA-Markov runs at each pixel.
#
# Inputs   outputs/gee/obs_<reserve>.tif, cam_<reserve>.tif, cam_band_index.csv
# Outputs  outputs/figures/Figure3a_maps_rainforest.png
#          outputs/figures/Figure3b_maps_mangrove_swamp.png
#          outputs/figures/Figure3c_maps_savannah.png
#          outputs/figures/Figure4_hotspot_zooms.png
#          outputs/tables/Figure4_hotspot_windows.csv

source("R/00_setup.R")
suppressPackageStartupMessages({ library(sf); library(ggspatial) })

OBS_BANDS <- c(paste0("class_", T0), paste0("class_", T1), "transition",
               paste0("p_high_", T0), paste0("p_high_", T1))
band_index <- read_csv(file.path(DIR_GEE, "cam_band_index.csv"), show_col_types = FALSE) |>
  arrange(band_number)
YEARS <- c(T0, T1, TARGET_YEARS)
YEAR_LAB <- c(paste(T0, "(classified)"), paste(T1, "(classified)"),
              paste(TARGET_YEARS, "(projected)"))
MAX_CELLS <- 9e5   # display resolution per panel; larger maps are aggregated by mode

maps_for <- function(reserve) {
  obs <- read_gee_tif(file.path(DIR_GEE, paste0("obs_", reserve, ".tif")), OBS_BANDS, reserve)
  cam <- read_gee_tif(file.path(DIR_GEE, paste0("cam_", reserve, ".tif")), band_index$band_name, reserve)
  s <- c(obs[[paste0("class_", c(T0, T1))]], cam[[paste0("y", TARGET_YEARS, "_mode")]])
  names(s) <- as.character(YEARS)
  s
}

to_df <- function(r, fact = 1) {
  if (fact > 1) r <- aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  d <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  d |> pivot_longer(-c(x, y), names_to = "year", values_to = "class") |>
    filter(!is.na(class)) |>
    mutate(class = factor(CLASS_NAMES[as.character(class)], levels = CLASS_NAMES),
           year = factor(YEAR_LAB[match(year, YEARS)], levels = YEAR_LAB))
}

# Graticule breaks (degrees) inside the panel extent
geo_breaks <- function(d, crs_r, nx, ny) {
  bb <- st_bbox(c(xmin = min(d$x), ymin = min(d$y), xmax = max(d$x), ymax = max(d$y)), crs = crs_r)
  g <- st_bbox(st_transform(st_as_sfc(bb), 4326))
  inside <- function(lo, hi, n) { b <- pretty(c(lo, hi), n); b[b > lo & b < hi] }
  list(x = inside(g[["xmin"]], g[["xmax"]], nx), y = inside(g[["ymin"]], g[["ymax"]], ny))
}

map_panel <- function(d, crs_r, title = NULL, show_legend = FALSE, base = 6, xb = 3) {
  gb <- geo_breaks(d, crs_r, xb, 3)
  ggplot(d) +
    geom_raster(aes(x, y, fill = class), hjust = 0.5, vjust = 0.5) +
    scale_fill_manual(values = CLASS_COLS, name = NULL, drop = FALSE) +
    coord_sf(crs = crs_r, datum = st_crs(4326), expand = FALSE) +
    scale_x_continuous(breaks = gb$x) + scale_y_continuous(breaks = gb$y) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.45,
                     height = unit(0.8, "mm"), pad_x = unit(1, "mm"), pad_y = unit(1, "mm"),
                     bar_cols = c("black", "white"), line_width = 0.3) +
    annotation_north_arrow(location = "tr", height = unit(3.5, "mm"), width = unit(2.5, "mm"),
                           pad_x = unit(0.8, "mm"), pad_y = unit(0.8, "mm"),
                           style = north_arrow_orienteering(text_size = 3)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = base) +
    theme(axis.text = element_text(size = base - 2), axis.text.x = element_text(angle = 0),
          panel.grid = element_line(colour = "grey80", linewidth = 0.2),
          plot.title = element_text(size = base, face = "bold", hjust = 0.5),
          legend.position = if (show_legend) "bottom" else "none",
          legend.spacing.x = unit(3, "mm"),
          legend.text = element_text(margin = margin(r = 8)),
          plot.margin = margin(1, 2, 1, 2))
}

plate <- function(reserves, file, title) {
  rows <- lapply(reserves, function(res) {
    s <- maps_for(res)
    fact <- max(1, ceiling(sqrt(ncell(s) / MAX_CELLS)))
    d <- to_df(s, fact)
    crs_r <- st_crs(crs(s))
    lab <- paste0(unname(RESERVE_LABEL[res]), ifelse(res %in% FLAGGED, " *", ""))
    panels <- lapply(levels(d$year), function(yl) {
      map_panel(d |> filter(year == yl), crs_r,
                title = if (res == reserves[1]) yl else NULL)
    })
    panels[[1]] <- panels[[1]] + labs(y = lab) +
      theme(axis.title.y = element_text(size = 7, face = "bold", angle = 90))
    wrap_plots(panels, nrow = 1)
  })
  legend_src <- map_panel(to_df(maps_for(reserves[1]), 8), st_crs(crs(maps_for(reserves[1]))),
                          show_legend = TRUE) +
    guides(fill = guide_legend(override.aes = list(colour = "grey40", linewidth = 0.2)))
  leg <- cowplot::get_legend(legend_src + theme(legend.text = element_text(size = 7)))
  g <- wrap_plots(rows, ncol = 1) / wrap_elements(leg) +
    plot_layout(heights = c(rep(1, length(rows)), 0.08)) +
    plot_annotation(title = title, caption = if (any(reserves %in% FLAGGED))
      paste0("* ", paste(unname(RESERVE_LABEL[intersect(reserves, FLAGGED)]), collapse = " and "),
             ": change between composites dominated by dry-season state (Section 2.4).") else NULL,
      theme = theme(plot.title = element_text(size = 9, face = "bold"),
                    plot.caption = element_text(size = 6, hjust = 0)))
  asp <- sapply(reserves, function(r) { e <- ext(rast(file.path(DIR_GEE, paste0("obs_", r, ".tif"))))
                                         (e$ymax - e$ymin) / (e$xmax - e$xmin) })
  h <- 22 + sum(pmin(pmax(40 * asp, 22), 60) + 8)
  ggsave(file.path(DIR_FIG, file), g, width = 180, height = h, units = "mm", dpi = 600)
  message("Wrote ", file)
}

by_biome <- split(RESERVES$reserve, RESERVES$biome)
REDO <- "--redo" %in% commandArgs(trailingOnly = TRUE)   # re-render plates that already exist
plate0 <- plate
plate <- function(reserves, file, title) {
  if (!REDO && file.exists(file.path(DIR_FIG, file))) return(message("Exists, skipped: ", file))
  plate0(reserves, file, title)
}
plate(by_biome[["Rainforest"]], "Figure3a_maps_rainforest.png", "(a) Rainforest reserves")
plate(by_biome[["Mangrove/swamp"]], "Figure3b_maps_mangrove_swamp.png", "(b) Mangrove and freshwater swamp reserves")
plate(by_biome[["Savannah"]], "Figure3c_maps_savannah.png", "(c) Savannah reserves")

# ---------------------------------------------------------------------------
# Figure 4: hotspot zooms. The window is the 1.5 km x 1.5 km block with the
# largest loss of high vegetation between 2019 and 2024 (classified maps).
# ---------------------------------------------------------------------------
WIN <- 1500
hotspot <- function(reserve) {
  s <- maps_for(reserve)
  loss <- (s[[as.character(T0)]] == 2) & (s[[as.character(T1)]] != 2)
  agg <- aggregate(loss, fact = WIN / res(s)[1], fun = "sum", na.rm = TRUE)
  full <- aggregate(!is.na(s[[1]]), fact = WIN / res(s)[1], fun = "sum", na.rm = TRUE)
  agg[full < 0.9 * (WIN / res(s)[1])^2] <- NA          # keep windows fully inside the reserve
  i <- which.max(values(agg))
  xy <- xyFromCell(agg, i)
  e <- ext(xy[1] - WIN / 2, xy[1] + WIN / 2, xy[2] - WIN / 2, xy[2] + WIN / 2)
  list(stack = crop(s, e), ext = e, crs = crs(s), full = s,
       loss_km2 = values(agg)[i] * prod(res(s)) / 1e6)
}
zoom_rows <- list(); windows <- list()
for (res in c("Stubbs_Creek", "Akanto")) {
  h <- hotspot(res)
  d <- to_df(h$stack[[c(as.character(T0), as.character(T1), as.character(TARGET_YEARS[2]))]])
  crs_r <- st_crs(h$crs)
  ov <- to_df(h$full[[as.character(T1)]], max(1, ceiling(sqrt(ncell(h$full) / 4e5))))
  ev <- as.vector(h$ext)   # xmin, xmax, ymin, ymax
  box <- st_as_sf(st_as_sfc(st_bbox(c(xmin = ev[[1]], ymin = ev[[3]], xmax = ev[[2]], ymax = ev[[4]]),
                                    crs = crs_r)))
  overview <- map_panel(ov, crs_r, title = paste0(unname(RESERVE_LABEL[res]), ", ", T1, " (overview)")) +
    geom_sf(data = box, fill = NA, colour = "red", linewidth = 0.5, inherit.aes = FALSE) +
    coord_sf(crs = crs_r, datum = st_crs(4326), expand = FALSE)
  zp <- lapply(levels(droplevels(d$year)), function(yl)
    map_panel(d |> filter(year == yl), crs_r, title = yl, xb = 1))
  # Overview on its own row; the three zoom panels share the full page width so
  # that each 1.5 km window prints at about 50 mm (about 1:30,000).
  zoom_rows[[res]] <- (overview | plot_spacer()) / wrap_plots(zp, nrow = 1) +
    plot_layout(heights = c(0.8, 1))
  c84 <- st_coordinates(st_centroid(st_transform(box, 4326)))
  windows[[res]] <- tibble(reserve = res, centre_lon = c84[1], centre_lat = c84[2],
                           window_m = WIN, high_veg_loss_km2_2019_2024 = h$loss_km2)
}
leg <- cowplot::get_legend(map_panel(to_df(maps_for("Akanto"), 4), st_crs(crs(maps_for("Akanto"))),
                                     show_legend = TRUE) + theme(legend.text = element_text(size = 7)))
g5 <- wrap_plots(zoom_rows, ncol = 1) / wrap_elements(leg) +
  plot_layout(heights = c(1, 1, 0.08)) +
  plot_annotation(caption = paste(strwrap(paste0("Red box: 1.5 km x 1.5 km window with the largest loss of high vegetation between ",
                                   T0, " and ", T1, ". At the printed width of 180 mm each zoom panel is at about 1:30,000; ",
                                   T1 + 9, " shows the modal class of 20 CA-Markov runs."), 150), collapse = "\n"),
                  theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
ggsave(file.path(DIR_FIG, "Figure4_hotspot_zooms.png"), g5, width = 180, height = 230, units = "mm", dpi = 600)
write_csv(bind_rows(windows), file.path(DIR_TAB, "Figure4_hotspot_windows.csv"))
message("Wrote Figure4_hotspot_zooms.png")
