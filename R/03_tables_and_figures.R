# 03_tables_and_figures.R
# Manuscript tables and figures rebuilt from the GEE and R outputs, following
# the reviewer requests (normalised importance, percentage-change metric plots on
# common axes, consolidated area table, uncertainty map).
#
# Outputs (outputs/tables/, outputs/figures/)
#   Table2_rf_accuracy.csv           RF accuracy (cross-validation, consensus-only model, blocked hold-out)
#   Table4_hindcast_accuracy.csv     CA-Markov hindcast, ensemble mean +/- SD, mode, persistence
#   Table3_high_vegetation_area.csv  high-vegetation area 2019-2033 with relative change
#   Figure2_feature_importance.png      importance normalised to 100 %
#   Figure5_pct_change_2019_2024.png    observed % change, common symmetric axis
#   Figure6_pct_change_2024_2033.png    projected % change from the classified 2024 map
#   FigureS1_pct_change_2028_2033.png   projected % change between the two target years
#   Figure7_uncertainty_2033.png        run agreement for the 2033 projection

source("R/00_setup.R")
suppressPackageStartupMessages({ library(sf); library(ggspatial) })

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
pm  <- function(m, s, d = 2) ifelse(is.na(s), fmt(m, d), paste0(fmt(m, d), " \u00b1 ", fmt(s, d)))
lab <- function(x) unname(RESERVE_LABEL[x])

# ---------------------------------------------------------------------------
# Table 2. RF accuracy
# ---------------------------------------------------------------------------
acc <- read_csv(file.path(DIR_GEE, "rf_accuracy.csv"), show_col_types = FALSE)
tab2 <- acc |>
  mutate(Reserve = lab(reserve),
         Assessment = factor(recode(assessment,
                             independent_ref_cv = "Reference points, 5-fold cross-validation",
                             consensus_only_ref = "Reference points, model trained on consensus labels only",
                             blocked_holdout = "Consensus labels, spatially blocked hold-out"),
                             c("Reference points, 5-fold cross-validation",
                               "Reference points, model trained on consensus labels only",
                               "Consensus labels, spatially blocked hold-out")),
         # Wilson score interval (stays within 0-1, unlike the Wald interval in OA_ci95)
         wz = 1.96, wc = (OA + wz^2 / (2 * n)) / (1 + wz^2 / n),
         wh = wz * sqrt(OA * (1 - OA) / n + wz^2 / (4 * n^2)) / (1 + wz^2 / n),
         OA = paste0(fmt(OA), " (", fmt(wc - wh), "-", fmt(wc + wh), ")"),
         Kappa = fmt(kappa)) |>
  select(Reserve, Biome = biome, Year = year, Assessment, n, OA, Kappa,
         UA_High = UA_high, UA_Low = UA_low, UA_No = UA_no,
         PA_High = PA_high, PA_Low = PA_low, PA_No = PA_no) |>
  mutate(across(starts_with(c("UA_", "PA_")), fmt)) |>
  arrange(Assessment, Reserve, Year)
write_csv(tab2, file.path(DIR_TAB, "Table2_rf_accuracy.csv"))

# ---------------------------------------------------------------------------
# Table 4. Hindcast accuracy
# ---------------------------------------------------------------------------
# One file per reserve (cam_hindcast_accuracy_<reserve>.csv)
hc <- list.files(DIR_GEE, pattern = "^cam_hindcast_accuracy.*\\.csv$", full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows() |> distinct()
runs <- hc |> filter(grepl("^run_", run)) |>
  group_by(reserve, biome) |>
  summarise(across(c(OA, kappa, OA_pixel, FoM, UA_high, UA_low, UA_no, PA_high, PA_low, PA_no),
                   list(m = mean, s = sd)), n = first(n), .groups = "drop")
null <- hc |> filter(run == "persistence_null") |>
  select(reserve, OA_null = OA, kappa_null = kappa, OA_pixel_null = OA_pixel)
mode_r <- hc |> filter(run == "ensemble_mode") |>
  select(reserve, OA_mode = OA, kappa_mode = kappa, OA_pixel_mode = OA_pixel, FoM_mode = FoM)
tab3 <- runs |> left_join(null, by = "reserve") |> left_join(mode_r, by = "reserve") |>
  transmute(Reserve = lab(reserve), Biome = biome, n_ref = n,
            OA = pm(OA_m, OA_s), Kappa = pm(kappa_m, kappa_s),
            `OA ensemble mode` = fmt(OA_mode), `OA persistence` = fmt(OA_null),
            `Kappa ensemble mode` = fmt(kappa_mode), `Kappa persistence` = fmt(kappa_null),
            `Pixel agreement` = pm(OA_pixel_m, OA_pixel_s),
            `Pixel agreement ensemble mode` = fmt(OA_pixel_mode),
            `Pixel agreement persistence` = fmt(OA_pixel_null),
            `Figure of merit` = pm(FoM_m, FoM_s, 3), `Figure of merit ensemble mode` = fmt(FoM_mode, 3),
            UA_High = pm(UA_high_m, UA_high_s), UA_Low = pm(UA_low_m, UA_low_s), UA_No = pm(UA_no_m, UA_no_s),
            PA_High = pm(PA_high_m, PA_high_s), PA_Low = pm(PA_low_m, PA_low_s), PA_No = pm(PA_no_m, PA_no_s)) |>
  arrange(Reserve)
write_csv(tab3, file.path(DIR_TAB, "Table4_hindcast_accuracy.csv"))

# ---------------------------------------------------------------------------
# Figure 2. Feature importance normalised to 100 %
# ---------------------------------------------------------------------------
imp <- read_csv(file.path(DIR_GEE, "rf_feature_importance.csv"), show_col_types = FALSE)
imp_s <- imp |> group_by(predictor) |>
  summarise(m = mean(importance_pct), s = sd(importance_pct), .groups = "drop") |>
  mutate(predictor = reorder(predictor, m))
g2 <- ggplot(imp_s, aes(m, predictor)) +
  geom_col(fill = "grey60", width = 0.7) +
  geom_errorbarh(aes(xmin = pmax(m - s, 0), xmax = m + s), height = 0.25) +
  labs(x = "Relative importance (%)", y = NULL,
       caption = "Mean \u00b1 SD across the 16 reserve-year models; each model sums to 100 %.") +
  theme_paper()
ggsave(file.path(DIR_FIG, "Figure2_feature_importance.png"), g2, width = 120, height = 90, units = "mm", dpi = 600)

if (!file.exists(file.path(DIR_TAB, "lsm_pct_change.csv"))) {
  message("Landscape-metric tables not found; run R/02_landscape_metrics.R first. ",
          "Only Tables 2 and 4 and Figure 2 were written.")
} else {
# ---------------------------------------------------------------------------
# Table 3. High-vegetation area dynamics
# ---------------------------------------------------------------------------
areas <- read_csv(file.path(DIR_TAB, "class_area_all_layers.csv"), show_col_types = FALSE) |>
  filter(class == 2)
obs_a <- areas |> filter(source == "observed") |> select(reserve, biome, year, area_km2)
proj_a <- areas |> filter(source == "projection") |> select(reserve, biome, year, run, area_km2)

# Per-run totals so that intervals for groups are computed correctly
group_runs <- function(df, grp) {
  df |> group_by(across(all_of(c(grp, "year", "run")))) |>
    summarise(area_km2 = sum(area_km2), .groups = "drop")
}
make_rows <- function(o, p, grp) {
  ob <- o |> group_by(across(all_of(c(grp, "year")))) |> summarise(a = sum(area_km2), .groups = "drop") |>
    pivot_wider(names_from = year, values_from = a, names_prefix = "obs_")
  pr <- group_runs(p, grp) |>
    left_join(ob, by = grp) |>
    group_by(across(all_of(c(grp, "year")))) |>
    summarise(m = mean(area_km2), lo = quantile(area_km2, .025), hi = quantile(area_km2, .975),
              chg_t1 = mean(pct_change(area_km2, .data[[paste0("obs_", T1)]])),
              .groups = "drop")
  pr_w <- pr |> pivot_wider(names_from = year, values_from = c(m, lo, hi, chg_t1))
  ob |> left_join(pr_w, by = grp)
}
fmt_row <- function(d) {
  b <- TARGET_YEARS[2]; a <- TARGET_YEARS[1]
  d |> transmute(
    `2019` = fmt(.data[[paste0("obs_", T0)]]),
    `2024` = fmt(.data[[paste0("obs_", T1)]]),
    `2028` = paste0(fmt(.data[[paste0("m_", a)]]), " (", fmt(.data[[paste0("lo_", a)]]), "-", fmt(.data[[paste0("hi_", a)]]), ")"),
    `2033` = paste0(fmt(.data[[paste0("m_", b)]]), " (", fmt(.data[[paste0("lo_", b)]]), "-", fmt(.data[[paste0("hi_", b)]]), ")"),
    `Change (%) 2019-2024` = fmt(pct_change(.data[[paste0("obs_", T1)]], .data[[paste0("obs_", T0)]]), 1),
    `Change (%) 2024-2033` = fmt(.data[[paste0("chg_t1_", b)]], 1),
    `Change (%) 2019-2033` = fmt(pct_change(.data[[paste0("m_", b)]], .data[[paste0("obs_", T0)]]), 1))
}
t_res <- make_rows(obs_a, proj_a, c("reserve", "biome"))
# Totals exclude the flagged reserves (see FLAGGED in 00_setup.R)
obs_k <- obs_a |> filter(!reserve %in% FLAGGED)
proj_k <- proj_a |> filter(!reserve %in% FLAGGED)
t_bio <- make_rows(obs_k, proj_k, "biome")
t_all <- make_rows(obs_k |> mutate(all = "All reserves"), proj_k |> mutate(all = "All reserves"), "all")
tab4 <- bind_rows(
  bind_cols(tibble(Reserve = paste0(lab(t_res$reserve), ifelse(t_res$reserve %in% FLAGGED, " *", "")),
                   Biome = t_res$biome), fmt_row(t_res)) |> arrange(Biome, Reserve),
  bind_cols(tibble(Reserve = paste("All", t_bio$biome, "reserves"), Biome = t_bio$biome), fmt_row(t_bio)),
  bind_cols(tibble(Reserve = "All reserves", Biome = "-"), fmt_row(t_all)))
attr(tab4, "note") <- paste0("* ", paste(lab(FLAGGED), collapse = " and "),
  ": change dominated by differences in dry-season state between composites; ",
  "excluded from biome and all-reserve totals.")
write_lines(attr(tab4, "note"), file.path(DIR_TAB, "Table3_note.txt"))
write_csv(tab4, file.path(DIR_TAB, "Table3_high_vegetation_area.csv"))

# ---------------------------------------------------------------------------
# Figures 5, 6 and S1. Percentage change of landscape metrics on a common axis
# ---------------------------------------------------------------------------
pct <- read_csv(file.path(DIR_TAB, "lsm_pct_change.csv"), show_col_types = FALSE) |>
  filter(level == "landscape") |>
  mutate(Metric = recode(metric, pd = "Patch density", ed = "Edge density",
                         lpi = "Largest patch index", ai = "Aggregation index"),
         Metric = factor(Metric, c("Patch density", "Edge density", "Largest patch index", "Aggregation index")),
         Reserve = paste0(lab(reserve), ifelse(reserve %in% FLAGGED, " *", "")),
         value = ifelse(basis == "observed", pct_hard, mean))
# Two-by-two panels, one per metric, reserves in alphabetical order, value
# labels on the bars. Each panel has its own x axis, symmetric about zero, so
# that small changes (e.g. aggregation index) remain visible.
# Observed change (dots = TRUE): the bar is the change between the classified
# maps; the dot and whisker are the mean and 95 % interval across the 20
# realisations drawn from the class probabilities. Realisations carry extra small patches, so
# their mean can lie away from the bar. Projected change: the bar is the
# ensemble mean and the whisker the 95 % interval across runs.
pct_plot <- function(d, dots = FALSE) {
  d <- d |> mutate(Reserve = factor(Reserve, levels = rev(sort(unique(Reserve)))))
  lims <- d |> group_by(Metric) |>
    summarise(lim = max(abs(c(value, q025, q975)), na.rm = TRUE) * 1.3, .groups = "drop")
  # The value label sits at the bar end. If the whisker would overlap it, the
  # label moves just beyond the whisker end on the same side.
  d <- d |> left_join(lims, by = "Metric") |>
    mutate(label = sub("-", "\u2212", formatC(value, format = "f", digits = 1)),
           # width of the label in data units (panel about 68 mm wide)
           zone = (nchar(label) * 1.4 + 3) * (2.1 * lim / 68),
           clash = !is.na(q025) & ifelse(value >= 0, q975 > value & q025 < value + zone,
                                                     q025 < value - 0 & q975 > value - zone),
           end = ifelse(clash, ifelse(value >= 0, q975, q025), value),
           lab_x = end + ifelse(value >= 0, 1, -1) * ifelse(clash & dots, 0.05, 0.03) * lim,
           hj = ifelse(value >= 0, 0, 1))
  blank <- bind_rows(lims |> mutate(x = -lim), lims |> mutate(x = lim)) |>
    mutate(Reserve = levels(d$Reserve)[1])
  inner_breaks <- function(l) {
    b <- scales::breaks_pretty(n = 6)(l)
    b[abs(b) <= 0.9 * max(abs(l))]   # no tick label at the panel edge
  }
  g <- ggplot(d, aes(value, Reserve)) +
    geom_blank(data = blank, aes(x = x, y = Reserve)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
    geom_col(aes(fill = "Classified maps"), width = 0.7) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.3,
                   linewidth = 0.35, colour = "black", na.rm = TRUE)
  if (dots) g <- g + geom_point(aes(x = mean, shape = "Realisations (mean and 95 % interval, n = 20)"),
                                size = 1.8, stroke = 0.5, fill = "white", colour = "black", na.rm = TRUE) +
    scale_shape_manual(values = 21, name = NULL)
  g + geom_text(aes(x = lab_x, label = label, hjust = hj), size = 2.6) +
    # free y as well, so that every panel shows the reserve names
    facet_wrap(~Metric, nrow = 2, scales = "free") +
    scale_fill_manual(values = "grey65", name = NULL, guide = if (dots) "legend" else "none") +
    scale_x_continuous(breaks = inner_breaks, labels = function(x) sub("-", "\u2212", x)) +
    labs(x = "Change relative to baseline (%)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey90", colour = "grey60"),
          strip.text = element_text(face = "bold", size = 9),
          axis.text.y = element_text(colour = "black", size = 8),
          panel.spacing = unit(6, "mm"),
          legend.position = "bottom", legend.text = element_text(size = 8),
          legend.key.size = unit(3.5, "mm"), legend.margin = margin(0, 0, 0, 0))
}
g6 <- pct_plot(pct |> filter(period == paste0(T0, "-", T1)), dots = TRUE)
g7 <- pct_plot(pct |> filter(period == paste0(T1, "-", TARGET_YEARS[2])))
g7s <- pct_plot(pct |> filter(period == paste(TARGET_YEARS, collapse = "-")))
ggsave(file.path(DIR_FIG, "Figure5_pct_change_2019_2024.png"), g6, width = 180, height = 130, units = "mm", dpi = 600)
ggsave(file.path(DIR_FIG, paste0("Figure6_pct_change_", T1, "_", TARGET_YEARS[2], ".png")), g7,
       width = 180, height = 130, units = "mm", dpi = 600)
ggsave(file.path(DIR_FIG, paste0("FigureS1_pct_change_", TARGET_YEARS[1], "_", TARGET_YEARS[2], ".png")), g7s,
       width = 180, height = 130, units = "mm", dpi = 600)

# ---------------------------------------------------------------------------
# Figure 7. Run agreement for the 2033 projection
# ---------------------------------------------------------------------------
band_index <- read_csv(file.path(DIR_GEE, "cam_band_index.csv"), show_col_types = FALSE) |> arrange(band_number)
last <- TARGET_YEARS[2]
panel <- function(reserve) {
  r <- read_gee_tif(file.path(DIR_GEE, paste0("cam_", reserve, ".tif")), band_index$band_name, reserve)
  a <- r[[paste0("y", last, "_agree")]]
  if (ncell(a) > 1.5e6) a <- aggregate(a, fact = 2, fun = "mean", na.rm = TRUE)
  df <- as.data.frame(a, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "agree"
  ggplot(df) +
    geom_tile(aes(x, y, fill = agree)) +
    scale_fill_gradientn(colours = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
                         limits = c(30, 100), oob = scales::squish, name = "Run agreement (%)") +
    coord_sf(crs = st_crs(crs(a)), datum = st_crs(4326), expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.5, height = unit(1, "mm")) +
    annotation_north_arrow(location = "tr", height = unit(5, "mm"), width = unit(4, "mm"),
                           style = north_arrow_orienteering(text_size = 4)) +
    labs(title = lab(reserve), x = NULL, y = NULL) +
    theme_paper(7) + theme(axis.text = element_text(size = 5))
}
g8 <- wrap_plots(lapply(RESERVES$reserve, panel), ncol = 2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(DIR_FIG, "Figure7_uncertainty_2033.png"), g8, width = 180, height = 240, units = "mm", dpi = 600)

}

message("Tables and figures written.")
