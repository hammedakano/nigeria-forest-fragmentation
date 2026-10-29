# combine_calibration.R
# Merges per-reserve calibration files written by calibrate_correlation.R when it
# was run for one reserve at a time, e.g.
#   for r in Adoru Akanto ...; do Rscript R/calibration/calibrate_correlation.R $r; done
source("R/00_setup.R")
for (name in c("calibration_ca_grid", "calibration_ca_params", "calibration_error_sigma")) {
  f <- list.files(DIR_TAB, pattern = paste0("^", name, "_.+\\.csv$"), full.names = TRUE)
  if (!length(f)) next
  d <- bind_rows(lapply(f, read_csv, show_col_types = FALSE)) |> arrange(reserve)
  write_csv(d, file.path(DIR_TAB, paste0(name, ".csv")))
  file.remove(f)
  message("Wrote ", name, ".csv from ", length(f), " files")
}
