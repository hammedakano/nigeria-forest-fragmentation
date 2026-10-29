# calibrate_correlation.R
# Calibrates the two spatial-correlation lengths used in
# gee/04_ca_markov_simulation.js (config.js: CA_PARAMS and ERROR_CORR_SIGMA).
#
# Why: both the CA transitions and the classification-error realisations are
# decided by one uniform random draw per pixel. With independent draws, single
# pixels change on their own, which adds isolated patches and edges. Patch
# density (PD) and edge density (ED) then rise because of the random draws,
# not because of land-cover change. Smoothing the draws with a Gaussian of
# sigma pixels makes neighbouring pixels change together. Sigma is chosen so
# that simulated maps have the same configuration as the classified maps.
#
#   Part A  CA_PARAMS (sigma, nucleation, neighbour power). The CA (same rule
#           as in GEE) is run from the classified 2019 map for 5 annual steps
#           for every combination in GRID_CA, and PD and ED of the simulated
#           2024 map are compared with the classified 2024 map (hindcast
#           calibration). The best combination per reserve goes into CA_PARAMS.
#   Part B  ERROR_CORR_SIGMA. Realisations are drawn from the RF class
#           probabilities with several sigmas and compared with the classified
#           map of the same year. Needs lulc_<reserve>.tif from
#           gee/02b_export_lulc_maps.js; skipped if they are missing.
#
# Inputs   outputs/gee/obs_<reserve>.tif, cam_markov_matrices.csv, lulc_<reserve>.tif
# Outputs  outputs/tables/calibration_ca_sigma.csv
#          outputs/tables/calibration_error_sigma.csv
# Usage    Rscript R/calibration/calibrate_correlation.R [reserve ...]
#          CALIB_PARTS=B Rscript ...   runs Part B only (default: A and B)
#
# The sigma with the smallest mean absolute relative error in PD and ED
# across the non-flagged reserves is printed. Copy it into gee/config.js.

source("R/00_setup.R")
set.seed(2024)

GRID_CA <- expand.grid(sigma = 0:4, nucleation = c(0.01, 0.001), power = c(1, 2))
SIGMAS_ERR <- c(1.5, 3, 5, 8, 12, 20)
N_REP <- 2
args <- commandArgs(trailingOnly = TRUE)
todo <- if (length(args)) args else RESERVES$reserve
# Running a subset (e.g. one reserve per process, in parallel) writes files with
# a suffix; R/calibration/combine_calibration.R merges them.
tag <- if (length(args)) paste0("_", paste(todo, collapse = "_")) else ""
out_csv <- function(name) file.path(DIR_TAB, paste0(name, tag, ".csv"))
PARTS <- strsplit(Sys.getenv("CALIB_PARTS", "AB"), "")[[1]]

# --- helpers -----------------------------------------------------------------
neigh_sum <- function(X) {          # sum over the 8 Moore neighbours, NA -> 0
  nr <- nrow(X); nc <- ncol(X)
  X[is.na(X)] <- 0
  P <- matrix(0, nr + 2, nc + 2); P[2:(nr + 1), 2:(nc + 1)] <- X
  out <- matrix(0, nr, nc)
  for (di in -1:1) for (dj in -1:1) if (di != 0 || dj != 0)
    out <- out + P[(2 + di):(nr + 1 + di), (2 + dj):(nc + 1 + dj)]
  out
}

# Uniform field with Gaussian spatial correlation (as correlatedUniform in GEE)
corr_uniform <- function(nr, nc, sigma) {
  z <- matrix(runif(nr * nc), nr, nc)
  if (sigma <= 0) return(z)
  rad <- ceiling(3 * sigma)
  k <- exp(-(-rad:rad)^2 / (2 * sigma^2)); k <- k / sum(k)
  sm <- function(v) as.numeric(stats::filter(v, k, circular = TRUE))
  s <- t(apply(apply(z, 2, sm), 1, sm))
  sd0 <- sqrt(sum(outer(k, k)^2) / 12)
  pnorm((s - 0.5) / sd0)
}

# One annual CA step, same rule as caStep() in 04_ca_markov_simulation.js
ca_step <- function(M, A, sigma, nucleation, power) {
  valid <- !is.na(M); nb <- pmax(neigh_sum(valid * 1), 1)
  N <- lapply(0:2, function(j) (neigh_sum((M == j) * 1) / nb)^power + nucleation)
  u <- corr_uniform(nrow(M), ncol(M), sigma); out <- M
  for (i in 0:2) {
    isI <- valid & M == i
    if (!any(isI)) next
    js <- setdiff(0:2, i)
    p <- lapply(js, function(j) N[[j + 1]] * A[i + 1, j + 1] / mean(N[[j + 1]][isI]))
    sh <- 0.99 / pmax(p[[1]] + p[[2]], 0.99)
    p0 <- p[[1]] * sh; p1 <- p[[2]] * sh
    out[isI & u < p0] <- js[1]
    out[isI & u >= p0 & u < p0 + p1] <- js[2]
  }
  out
}

# Realisation drawn from class probabilities (as drawStart() in GEE)
draw_real <- function(pno, plow, phigh, sigma) {
  tot <- pmax(pno + plow + phigh, 1)
  u <- corr_uniform(nrow(pno), ncol(pno), sigma)
  out <- matrix(2L, nrow(pno), ncol(pno))
  out[u < (pno + plow) / tot] <- 1L
  out[u < pno / tot] <- 0L
  out[is.na(pno)] <- NA
  out
}

metrics <- function(M, template) {
  r <- rast(template); values(r) <- as.vector(t(M))
  m <- calculate_lsm(r, what = c("lsm_l_pd", "lsm_l_ed", "lsm_c_ca"), directions = 8, progress = FALSE)
  c(pd = m$value[m$metric == "pd"], ed = m$value[m$metric == "ed"],
    high_ha = sum(m$value[m$metric == "ca" & m$class == 2]))
}

pick <- function(d, by = "sigma", value_cols = c("rel_pd", "rel_ed")) {
  d |> filter(!reserve %in% FLAGGED) |>
    group_by(across(all_of(by))) |>
    summarise(mae = mean(abs(c(.data[[value_cols[1]]], .data[[value_cols[2]]]))), .groups = "drop") |>
    arrange(mae)
}

if ("A" %in% PARTS) {
  # --- Part A: CA_PARAMS --------------------------------------------------------
  mm <- read_csv(file.path(DIR_GEE, "cam_markov_matrices.csv"), show_col_types = FALSE)
  resA <- list()
  for (res in todo) {
    message("Part A (CA parameters): ", res)
    A <- matrix(unlist(mm[mm$reserve == res, paste0("A_", c("00", "01", "02", "10", "11", "12", "20", "21", "22"))]),
                3, byrow = TRUE)
    obs <- read_gee_tif(file.path(DIR_GEE, paste0("obs_", res, ".tif")),
                        c("c0", "c1", "tr", "ph0", "ph1"), res)
    M0 <- as.matrix(obs$c0, wide = TRUE)
    ref <- metrics(as.matrix(obs$c1, wide = TRUE), obs$c1)
    for (g in seq_len(nrow(GRID_CA))) for (k in seq_len(N_REP)) {
      M <- M0
      for (y in seq_len(T1 - T0)) M <- ca_step(M, A, GRID_CA$sigma[g], GRID_CA$nucleation[g], GRID_CA$power[g])
      m <- metrics(M, obs$c1)
      resA[[length(resA) + 1]] <- tibble(reserve = res, GRID_CA[g, ], rep = k,
        pd = m["pd"], ed = m["ed"], high_ha = m["high_ha"],
        pd_obs = ref["pd"], ed_obs = ref["ed"], high_ha_obs = ref["high_ha"])
    }
  }
  resA <- bind_rows(resA) |>
    mutate(rel_pd = pd / pd_obs - 1, rel_ed = ed / ed_obs - 1, rel_high = high_ha / high_ha_obs - 1)
  write_csv(resA, out_csv("calibration_ca_grid"))
  best <- resA |>
    group_by(reserve, sigma, nucleation, power) |>
    summarise(across(c(pd, ed, pd_obs, ed_obs, rel_pd, rel_ed, rel_high), mean), .groups = "drop") |>
    mutate(error = (abs(rel_pd) + abs(rel_ed)) / 2) |>
    group_by(reserve) |> slice_min(error, n = 1, with_ties = FALSE) |> ungroup()
  write_csv(best, out_csv("calibration_ca_params"))
  message("Best CA parameters per reserve (copy into CA_PARAMS in gee/config.js):")
  print(as.data.frame(best))
  message("Single parameter set for all non-flagged reserves:")
  print(head(as.data.frame(pick(resA, c("sigma", "nucleation", "power"))), 5))
}

# --- Part B: ERROR_CORR_SIGMA -------------------------------------------------
lulc_files <- file.path(DIR_GEE, paste0("lulc_", todo, ".tif"))
if (!"B" %in% PARTS) {
  message("Part B skipped (CALIB_PARTS).")
} else if (!all(file.exists(lulc_files))) {
  message("Part B skipped: lulc_<reserve>.tif not found (run gee/02b_export_lulc_maps.js).")
} else {
  resB <- list()
  for (res in todo) {
    message("Part B (error sigma): ", res)
    bands <- as.vector(outer(c("class", "p_no", "p_low", "p_high", "n_obs"), c(T0, T1), paste, sep = "_"))
    L <- read_gee_tif(file.path(DIR_GEE, paste0("lulc_", res, ".tif")), bands, res)
    for (yr in c(T0, T1)) {
      g <- function(b) as.matrix(L[[paste0(b, "_", yr)]], wide = TRUE)
      hard <- g("class"); ref <- metrics(hard, L[[1]])
      for (sg in SIGMAS_ERR) for (k in seq_len(N_REP)) {
        m <- metrics(draw_real(g("p_no"), g("p_low"), g("p_high"), sg), L[[1]])
        resB[[length(resB) + 1]] <- tibble(reserve = res, year = yr, sigma = sg, rep = k,
          pd = m["pd"], ed = m["ed"], pd_obs = ref["pd"], ed_obs = ref["ed"])
      }
    }
  }
  resB <- bind_rows(resB) |> mutate(rel_pd = pd / pd_obs - 1, rel_ed = ed / ed_obs - 1)
  write_csv(resB, out_csv("calibration_error_sigma"))
  message("ERROR_CORR_SIGMA, mean absolute relative error of PD and ED (non-flagged reserves):")
  print(as.data.frame(pick(resB)))
}
