# 02_landscape_metrics.R
# Landscape metrics (FRAGSTATS definitions, landscapemetrics) for the observed
# maps, the starting-map realisations and every CA-Markov run, followed by
# summaries with uncertainty intervals and normalised percentage changes.
#
# Inputs (from Google Drive folder exported by GEE, copied to outputs/gee/)
#   obs_<reserve>.tif        03_change_analysis.js
#   cam_<reserve>.tif        04_ca_markov_simulation.js
#   cam_band_index.csv       04_ca_markov_simulation.js
#
# Outputs (outputs/tables/)
#   lsm_all_layers.csv       one row per reserve, layer, level, class and metric
#   class_area_all_layers.csv  class areas (km2) per layer
#   lsm_summary.csv          observed value, realisation interval and ensemble
#                            mean, SD and 2.5/97.5 % quantiles per reserve and year
#   lsm_pct_change.csv       percentage change 2019-2024 (classified maps, with the
#                            95 % interval across map realisations), 2024-2033
#                            (classified 2024 map to each run) and 2028-2033
#                            (paired by run), with 95 % intervals
#   class_area_summary.csv   class areas with ensemble mean, SD and intervals
#
# Metric definitions (eight-cell neighbourhood, landscape level unless stated)
#   PD  = N / A * 10,000 * 100        patches per 100 ha
#   ED  = E / A * 10,000              m per ha
#   LPI = max(a_ij) / A * 100         % of landscape
#   AI  = g_ii / max(g_ii) * 100      % (class level, then area-weighted at landscape level)

source("R/00_setup.R")

METRICS_L <- c("lsm_l_pd", "lsm_l_ed", "lsm_l_lpi", "lsm_l_ai")
METRICS_C <- c("lsm_c_pd", "lsm_c_ed", "lsm_c_lpi", "lsm_c_ai", "lsm_c_ca")
OBS_BANDS <- c(paste0("class_", T0), paste0("class_", T1), "transition",
               paste0("p_high_", T0), paste0("p_high_", T1))

band_index <- read_csv(file.path(DIR_GEE, "cam_band_index.csv"), show_col_types = FALSE) |>
  arrange(band_number)

# Describe a layer name: source, year and run
describe_layer <- function(nm) {
  tibble(layer_name = nm) |>
    mutate(
      source = case_when(
        grepl("^class_", nm) ~ "observed",
        grepl("^i\\d{4}_r", nm) ~ "observed_realisation",
        grepl("^h\\d{4}_r", nm) ~ "hindcast",
        grepl("^y\\d{4}_r", nm) ~ "projection"),
      year = as.integer(sub("^[a-z_]*?(\\d{4}).*$", "\\1", nm)),
      run = suppressWarnings(as.integer(sub("^.*_r(\\d+)$", "\\1", nm))))
}

metrics_for_reserve <- function(reserve) {
  message("Landscape metrics: ", reserve)
  obs <- read_gee_tif(file.path(DIR_GEE, paste0("obs_", reserve, ".tif")), OBS_BANDS, reserve)
  cam <- read_gee_tif(file.path(DIR_GEE, paste0("cam_", reserve, ".tif")), band_index$band_name, reserve)

  keep_cam <- grep("^(i\\d{4}|h\\d{4}|y\\d{4})_r\\d+$", names(cam), value = TRUE)
  stack <- c(obs[[c(paste0("class_", T0), paste0("class_", T1))]], cam[[keep_cam]])
  stopifnot(!is.lonlat(stack))

  lsm <- calculate_lsm(stack, what = c(METRICS_L, METRICS_C), directions = 8,
                       progress = FALSE)
  lsm <- lsm |>
    mutate(layer_name = names(stack)[layer]) |>
    filter(level == "landscape" | class == 2) |>
    select(layer_name, level, class, metric, value)

  areas <- lapply(names(stack), function(nm) {
    f <- freq(stack[[nm]])
    tibble(layer_name = nm, class = f$value, area_km2 = f$count * prod(res(stack)) / 1e6)
  }) |> bind_rows()

  list(lsm = lsm |> mutate(reserve = reserve) |> left_join(describe_layer(names(stack)), by = "layer_name"),
       area = areas |> mutate(reserve = reserve) |> left_join(describe_layer(names(stack)), by = "layer_name"))
}

# Per-reserve results are cached in outputs/tables/cache/, so the script can be
# resumed and several copies can run in parallel (each claims a reserve with a
# lock file). Run e.g. two copies at once:
#   Rscript R/02_landscape_metrics.R --worker & Rscript R/02_landscape_metrics.R --worker
# then once more without --worker to assemble the tables.
DIR_CACHE <- file.path(DIR_TAB, "cache")
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)
cache_file <- function(r) file.path(DIR_CACHE, paste0("lsm_", r, ".rds"))
worker_only <- "--worker" %in% commandArgs(trailingOnly = TRUE)

for (r in RESERVES$reserve) {
  if (file.exists(cache_file(r))) next
  lock <- file.path(DIR_CACHE, paste0("lsm_", r, ".lock"))
  if (!dir.create(lock, showWarnings = FALSE)) next   # another worker has it
  saveRDS(metrics_for_reserve(r), cache_file(r))
  unlink(lock, recursive = TRUE)
}
if (worker_only) quit(save = "no")
todo <- RESERVES$reserve[!file.exists(cache_file(RESERVES$reserve))]
if (length(todo) > 0) stop("Not yet computed (another worker may still be running): ",
                           paste(todo, collapse = ", "))
res <- lapply(RESERVES$reserve, function(r) readRDS(cache_file(r)))
lsm_all  <- bind_rows(lapply(res, `[[`, "lsm"))  |> left_join(RESERVES |> select(reserve, biome), by = "reserve")
area_all <- bind_rows(lapply(res, `[[`, "area")) |> left_join(RESERVES |> select(reserve, biome), by = "reserve")
write_csv(lsm_all, file.path(DIR_TAB, "lsm_all_layers.csv"))
write_csv(area_all, file.path(DIR_TAB, "class_area_all_layers.csv"))

# ---------------------------------------------------------------------------
# Summaries: observed hard map plus interval from realisations; ensemble stats
# ---------------------------------------------------------------------------
summ <- function(df, value_col) {
  df |>
    summarise(n_runs = sum(!is.na(.data[[value_col]])),
              mean = mean(.data[[value_col]], na.rm = TRUE),
              sd = sd(.data[[value_col]], na.rm = TRUE),
              q025 = quantile(.data[[value_col]], 0.025, na.rm = TRUE),
              q975 = quantile(.data[[value_col]], 0.975, na.rm = TRUE),
              .groups = "drop")
}

lsm_summary <- bind_rows(
  lsm_all |> filter(source == "observed") |>
    transmute(reserve, biome, source, year, level, class, metric, n_runs = 1L,
              mean = value, sd = NA_real_, q025 = NA_real_, q975 = NA_real_),
  lsm_all |> filter(source != "observed") |>
    group_by(reserve, biome, source, year, level, class, metric) |> summ("value")
) |> arrange(reserve, level, metric, year, source)
write_csv(lsm_summary, file.path(DIR_TAB, "lsm_summary.csv"))

area_summary <- bind_rows(
  area_all |> filter(source == "observed") |>
    transmute(reserve, biome, source, year, class, n_runs = 1L, mean = area_km2,
              sd = NA_real_, q025 = NA_real_, q975 = NA_real_),
  area_all |> filter(source != "observed") |>
    group_by(reserve, biome, source, year, class) |> summ("area_km2")
) |> mutate(class_name = CLASS_NAMES[as.character(class)])
write_csv(area_summary, file.path(DIR_TAB, "class_area_summary.csv"))

# ---------------------------------------------------------------------------
# Normalised percentage change (reviewer request: replace absolute deltas)
# ---------------------------------------------------------------------------
wide_runs <- function(src, yrs) {
  lsm_all |>
    filter(source == src, year %in% yrs) |>
    select(reserve, biome, level, class, metric, run, year, value) |>
    pivot_wider(names_from = year, values_from = value, names_prefix = "v")
}

# 2019 -> 2024 observed. Interval from paired realisations of both years.
obs_hard <- lsm_all |> filter(source == "observed") |>
  select(reserve, biome, level, class, metric, year, value) |>
  pivot_wider(names_from = year, values_from = value, names_prefix = "v") |>
  mutate(pct = pct_change(.data[[paste0("v", T1)]], .data[[paste0("v", T0)]]))
has_real <- any(lsm_all$source == "observed_realisation")
if (has_real) {
  obs_real <- wide_runs("observed_realisation", c(T0, T1)) |>
    mutate(pct = pct_change(.data[[paste0("v", T1)]], .data[[paste0("v", T0)]])) |>
    group_by(reserve, biome, level, class, metric) |> summ("pct")
} else {
  # PROPAGATE_CLASSIFICATION_ERROR was false in GEE: no intervals for observed change
  obs_real <- obs_hard |> distinct(reserve, biome, level, class, metric) |>
    mutate(n_runs = 0L, mean = NA_real_, sd = NA_real_, q025 = NA_real_, q975 = NA_real_)
}
chg_obs <- obs_hard |>
  select(reserve, biome, level, class, metric, pct_hard = pct) |>
  left_join(obs_real, by = c("reserve", "biome", "level", "class", "metric")) |>
  mutate(period = paste0(T0, "-", T1), basis = "observed")

# 2028 -> 2033 projected, paired within each run
y_a <- TARGET_YEARS[1]; y_b <- TARGET_YEARS[2]
chg_proj <- wide_runs("projection", TARGET_YEARS) |>
  mutate(pct = pct_change(.data[[paste0("v", y_b)]], .data[[paste0("v", y_a)]])) |>
  group_by(reserve, biome, level, class, metric) |> summ("pct") |>
  mutate(period = paste0(y_a, "-", y_b), basis = "projected", pct_hard = NA_real_)

# 2024 -> 2033: every run starts from the classified 2024 map (CA_START =
# 'classified' in gee/config.js), so each run is compared with that map
runs <- unique(na.omit(lsm_all$run[lsm_all$source == "projection"]))
base_t1 <- lsm_all |> filter(source == "observed", year == T1) |> select(-run) |>
  tidyr::crossing(run = runs)
chg_long <- bind_rows(
  base_t1,
  lsm_all |> filter(source == "projection", year == y_b)
) |>
  select(reserve, biome, level, class, metric, run, year, value) |>
  pivot_wider(names_from = year, values_from = value, names_prefix = "v") |>
  mutate(pct = pct_change(.data[[paste0("v", y_b)]], .data[[paste0("v", T1)]])) |>
  group_by(reserve, biome, level, class, metric) |> summ("pct") |>
  mutate(period = paste0(T1, "-", y_b), basis = "projected", pct_hard = NA_real_)

pct_tab <- bind_rows(chg_obs, chg_proj, chg_long) |>
  relocate(reserve, biome, period, basis, level, class, metric, pct_hard)
write_csv(pct_tab, file.path(DIR_TAB, "lsm_pct_change.csv"))

message("Done. Tables written to ", DIR_TAB)
