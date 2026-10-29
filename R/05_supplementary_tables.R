# 05_supplementary_tables.R
# Supplementary tables for the revised manuscript.
#
# Outputs (outputs/tables/)
#   TableS1_accuracy_all_assessments.csv   RF accuracy for the three assessments, with error matrices
#   TableS2_training_samples.csv           consensus training samples by class, label tier and split
#   TableS3_class_areas.csv                class areas 2019-2033 with 95 % intervals
#   TableS4_landscape_metrics.csv          PD, ED, LPI, AI (landscape level) with 95 % intervals
#   TableS5_transition_matrices.csv        transition areas 2019-2024, five-year and annual transition matrices
#   TableS6_ca_calibration.csv             calibrated CA parameters and hindcast fit
#   TableS7_composite_diagnostics.csv      burned area, observation counts and NDVI per composite
#   TableS8_realisation_calibration.csv    PD and ED of probability realisations relative to the classified maps

source("R/00_setup.R")
lab <- function(x) unname(RESERVE_LABEL[x])
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
iv <- function(m, lo, hi, d = 2) ifelse(is.na(lo), f2(m, d), paste0(f2(m, d), " (", f2(lo, d), "-", f2(hi, d), ")"))

# S1 -------------------------------------------------------------------------
acc <- read_csv(file.path(DIR_GEE, "rf_accuracy.csv"), show_col_types = FALSE)
s1 <- acc |>
  mutate(Reserve = lab(reserve),
         Assessment = recode(assessment,
                             independent_ref_cv = "Reference points, 5-fold CV",
                             consensus_only_ref = "Reference points, consensus-only model",
                             blocked_holdout = "Consensus labels, blocked hold-out")) |>
  transmute(Reserve, Biome = biome, Year = year, Assessment, n,
            OA = f2(OA), Kappa = f2(kappa),
            across(c(UA_high, UA_low, UA_no, PA_high, PA_low, PA_no, F1_high, F1_low, F1_no), f2),
            across(starts_with("cm_"), ~ as.integer(round(.x)))) |>
  arrange(Assessment, Reserve, Year)
names(s1) <- sub("^cm_r(\\d)_m(\\d)$", "Ref\\1_Map\\2", names(s1))
write_csv(s1, file.path(DIR_TAB, "TableS1_accuracy_all_assessments.csv"))

# S2 -------------------------------------------------------------------------
sc <- read_csv(file.path(DIR_GEE, "sample_counts.csv"), show_col_types = FALSE)
s2 <- sc |> transmute(Reserve = lab(reserve), Year = year, Class = class_name,
                      Total = n_total, Train = n_train, Test = n_test,
                      `Tier 1` = n_tier1, `Tier 2` = n_tier2, `Tier 3` = n_tier3)
write_csv(s2, file.path(DIR_TAB, "TableS2_training_samples.csv"))

# S3 -------------------------------------------------------------------------
ar <- read_csv(file.path(DIR_TAB, "class_area_summary.csv"), show_col_types = FALSE)
obs <- ar |> filter(source == "observed") |> select(reserve, biome, year, class, value = mean)
real <- ar |> filter(source == "observed_realisation") |> select(reserve, year, class, lo = q025, hi = q975)
proj <- ar |> filter(source == "projection") |> select(reserve, biome, year, class, value = mean, lo = q025, hi = q975)
s3 <- bind_rows(obs |> left_join(real, by = c("reserve", "year", "class")), proj) |>
  mutate(Reserve = paste0(lab(reserve), ifelse(reserve %in% FLAGGED, " *", "")),
         Class = CLASS_NAMES[as.character(class)], cell = iv(value, lo, hi)) |>
  select(Reserve, Biome = biome, Class, year, cell) |>
  pivot_wider(names_from = year, values_from = cell) |>
  arrange(Biome, Reserve, factor(Class, levels = CLASS_NAMES))
write_csv(s3, file.path(DIR_TAB, "TableS3_class_areas.csv"))

# S4 -------------------------------------------------------------------------
ls <- read_csv(file.path(DIR_TAB, "lsm_summary.csv"), show_col_types = FALSE) |> filter(level == "landscape")
lo <- ls |> filter(source == "observed") |> select(reserve, biome, year, metric, value = mean)
lr <- ls |> filter(source == "observed_realisation") |> select(reserve, year, metric, lo = q025, hi = q975)
lp <- ls |> filter(source == "projection") |> select(reserve, biome, year, metric, value = mean, lo = q025, hi = q975)
s4 <- bind_rows(lo |> left_join(lr, by = c("reserve", "year", "metric")), lp) |>
  mutate(Reserve = paste0(lab(reserve), ifelse(reserve %in% FLAGGED, " *", "")),
         Metric = recode(metric, pd = "PD (patches/100 ha)", ed = "ED (m/ha)", lpi = "LPI (%)", ai = "AI (%)"),
         cell = iv(value, lo, hi, 1)) |>
  select(Reserve, Biome = biome, Metric, year, cell) |>
  pivot_wider(names_from = year, values_from = cell) |>
  arrange(Biome, Reserve, Metric)
write_csv(s4, file.path(DIR_TAB, "TableS4_landscape_metrics.csv"))

# S5 -------------------------------------------------------------------------
mm <- read_csv(file.path(DIR_GEE, "cam_markov_matrices.csv"), show_col_types = FALSE)
tr <- read_csv(file.path(DIR_GEE, "change_transitions.csv"), show_col_types = FALSE) |>
  mutate(from_class = as.integer(from_class), to_class = as.integer(to_class))
cls <- c("No", "Low", "High")
s5 <- lapply(seq_len(nrow(mm)), function(k) {
  r <- mm[k, ]
  expand_grid(i = 0:2, j = 0:2) |>
    left_join(tr |> filter(reserve == r$reserve) |> select(i = from_class, j = to_class, area_km2),
              by = c("i", "j")) |>
    mutate(Reserve = lab(r$reserve), From = cls[i + 1], To = cls[j + 1],
           `Area 2019-2024 (km2)` = f2(area_km2),
           `P (2019-2024)` = f2(unlist(r[paste0("P_", i, j)]), 3),
           `Annual A` = f2(unlist(r[paste0("A_", i, j)]), 3),
           `Root method` = r$root_method, `Root error` = f2(r$root_error, 3))
}) |> bind_rows() |> select(-i, -j, -area_km2)
write_csv(s5, file.path(DIR_TAB, "TableS5_transition_matrices.csv"))

# S6 -------------------------------------------------------------------------
cp <- read_csv(file.path(DIR_TAB, "calibration_ca_params.csv"), show_col_types = FALSE)
s6 <- cp |> transmute(Reserve = paste0(lab(reserve), ifelse(reserve %in% FLAGGED, " *", "")),
                      `sigma (pixels)` = sigma, `Nucleation e` = nucleation, `Power k` = power,
                      `PD observed 2024` = f2(pd_obs, 1), `PD simulated 2024` = f2(pd, 1),
                      `ED observed 2024` = f2(ed_obs, 1), `ED simulated 2024` = f2(ed, 1),
                      `High-vegetation area error (%)` = f2(100 * rel_high, 1))
write_csv(s6, file.path(DIR_TAB, "TableS6_ca_calibration.csv"))

# S7 -------------------------------------------------------------------------
cd <- read_csv(file.path(DIR_GEE, "composite_diagnostics.csv"), show_col_types = FALSE)
s7 <- cd |> transmute(Reserve = lab(reserve), Year = year, `Scenes` = n_scenes,
                      `Clear observations per pixel (mean)` = f2(clear_obs_mean, 1),
                      `Burned area (%)` = f2(burned_pct, 1),
                      `Median-composite NDVI (median)` = f2(ndvi_median_p50),
                      `Greenest-observation NDVI (median)` = f2(ndvi_max_p50))
write_csv(s7, file.path(DIR_TAB, "TableS7_composite_diagnostics.csv"))
# S8 -------------------------------------------------------------------------
ce <- read_csv(file.path(DIR_TAB, "calibration_error_sigma.csv"), show_col_types = FALSE)
s8 <- ce |> group_by(reserve, year, sigma) |>
  summarise(pd = mean(pd), ed = mean(ed), pd_obs = first(pd_obs), ed_obs = first(ed_obs), .groups = "drop") |>
  mutate(Reserve = paste0(lab(reserve), ifelse(reserve %in% FLAGGED, " *", "")),
         cell = paste0(f2(100 * (pd / pd_obs - 1), 0), " / ", f2(100 * (ed / ed_obs - 1), 0))) |>
  select(Reserve, Year = year, sigma, cell) |>
  pivot_wider(names_from = sigma, values_from = cell, names_prefix = "sigma = ")
write_csv(s8, file.path(DIR_TAB, "TableS8_realisation_calibration.csv"))
message("Supplementary tables written.")
