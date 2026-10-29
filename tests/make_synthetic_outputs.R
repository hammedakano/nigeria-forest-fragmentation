# tests/make_synthetic_outputs.R
# Writes small synthetic files with the same names, bands and columns as the
# GEE exports into outputs/gee/, so the R scripts can be tested end to end
# without running Earth Engine. Values are random and have no meaning.
#
#   Rscript tests/make_synthetic_outputs.R && Rscript R/02_landscape_metrics.R && Rscript R/03_tables_and_figures.R
#
# Run it in a copy of the repository: it overwrites outputs/gee/.

suppressPackageStartupMessages({ library(terra); library(readr); library(dplyr) })
set.seed(1)
dir.create("outputs/gee", recursive = TRUE, showWarnings = FALSE)
res <- read_csv("data/reserves.csv", show_col_types = FALSE)
utm <- c(Adoru = 32632, Akanto = 32632, Gambari = 32631, Maradun = 32632, Mayo_Ndaga = 32632,
         Stubbs_Creek = 32632, Upper_Orashi = 32632, Yo = 32633)
N_RUNS <- 3

patchy <- function(template, p_change = 0, base = NULL) {
  f <- focal(init(template, runif), w = 9, fun = mean)
  m <- classify(f, rbind(c(-Inf, 0.45, 0), c(0.45, 0.55, 1), c(0.55, Inf, 2)))
  if (!is.null(base)) m <- ifel(init(template, runif) < p_change, m, base)
  m
}

bands <- c()
for (pre in c("h2024", "y2028", "y2033")) {
  bands <- c(bands, sprintf("%s_r%02d", pre, seq_len(N_RUNS)),
             paste0(pre, c("_p_no", "_p_low", "_p_high", "_mode", "_agree")))
}
bands <- c(bands, sprintf("i2019_r%02d", seq_len(N_RUNS)), sprintf("i2024_r%02d", seq_len(N_RUNS)))
write_csv(tibble(band_number = seq_along(bands), band_name = bands), "outputs/gee/cam_band_index.csv")

for (r in res$reserve) {
  t <- rast(nrows = 120, ncols = 150, xmin = 500000, xmax = 501500, ymin = 800000, ymax = 801200,
            crs = paste0("EPSG:", utm[[r]]))
  c19 <- patchy(t); c24 <- patchy(t, 0.2, c19)
  obs <- c(c19, c24, c19 * 3 + c24, init(t, function(n) sample(0:100, n, TRUE)),
           init(t, function(n) sample(0:100, n, TRUE)))
  obs[1:5, ] <- NA
  writeRaster(obs, sprintf("outputs/gee/obs_%s.tif", r), datatype = "INT1U", NAflag = 255, overwrite = TRUE)
  layers <- lapply(bands, function(b) {
    if (grepl("_r\\d+$", b)) patchy(t, 0.3, c24)
    else if (grepl("_mode$", b)) c24
    else init(t, function(n) sample(30:100, n, TRUE))
  })
  cam <- rast(layers); cam[1:5, ] <- NA
  writeRaster(cam, sprintf("outputs/gee/cam_%s.tif", r), datatype = "INT1U", NAflag = 255, overwrite = TRUE)
}

acc_row <- function(reserve, biome, ...) {
  tibble(reserve, biome, ..., n = 70, OA = runif(1, .8, .95), OA_ci95 = .06, kappa = runif(1, .7, .9),
         PA_high = runif(1, .7, 1), PA_low = runif(1, .7, 1), PA_no = runif(1, .5, 1),
         UA_high = runif(1, .7, 1), UA_low = runif(1, .7, 1), UA_no = runif(1, .5, 1))
}
rf_acc <- bind_rows(lapply(seq_len(nrow(res)), function(i) bind_rows(
  lapply(c(2019, 2024), function(y) bind_rows(
    acc_row(res$reserve[i], res$biome[i], year = y, assessment = "independent_ref_cv"),
    acc_row(res$reserve[i], res$biome[i], year = y, assessment = "consensus_only_ref"),
    acc_row(res$reserve[i], res$biome[i], year = y, assessment = "blocked_holdout"))))))
write_csv(rf_acc, "outputs/gee/rf_accuracy.csv")

preds <- c("B2", "B3", "B4", "B8", "B11", "B12", "NDVI", "NDWI", "NBR", "SAVI", "NDVI_p90", "NBR_p90", "NDVI_p10")
imp <- expand.grid(reserve = res$reserve, year = c(2019, 2024), predictor = preds) |>
  as_tibble() |> mutate(importance_raw = runif(n(), 5, 30)) |>
  group_by(reserve, year) |> mutate(importance_pct = importance_raw / sum(importance_raw) * 100) |> ungroup()
write_csv(imp, "outputs/gee/rf_feature_importance.csv")

hc <- bind_rows(lapply(seq_len(nrow(res)), function(i) bind_rows(
  lapply(c(sprintf("run_%02d", seq_len(N_RUNS)), "ensemble_mode", "persistence_null"), function(run)
    acc_row(res$reserve[i], res$biome[i], run = run) |>
      mutate(OA_pixel = runif(1, .7, .9), FoM = runif(1, .05, .2))))))
write_csv(hc, "outputs/gee/cam_hindcast_accuracy.csv")

message("Synthetic outputs written to outputs/gee/")
