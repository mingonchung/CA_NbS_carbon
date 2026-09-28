### CA carbon revision pipeline
### Stage 3d. Time-invariant matching covariates
###
### Nine layers, no year dimension. Seven rasters extracted at the point, two
### polygon layers joined by point in polygon.
###
###   elevation, slope, aspect          SRTM 90 m
###   ppt_normal, tmean_normal          PRISM 1991-2020 normals
###   pop_density                       WorldPop, 2000 to 2002 mean
###   city_travel_time                  Travel time to cities, class 11
###   ecoregion_l3                      EPA Level III, exact matching stratum
###   huc8                              WBD HUC12 truncated, exact stratum
###
### The stage 3c and 3d boundary is dimensional rather than thematic. Stage 3c
### layers carry an event year and produce pixel-year records. Stage 3d layers
### produce one row per pixel, which is why the two zone joins sit here.
###
### Probe mode is a mode rather than a separate script, for the same reason it
### is in 3a and 3b. Keeping the registry and the probe in one file means they
### cannot drift.
###
### Run probe first. It reports the resolved files, CRS match, native
### resolution, declared NA flag, sampled value range, and the actual field
### names in the two shapefiles. Nothing should be written until the probe has
### been read, because the two zone field names are the only entries in
### ca_zone_layers that are assumed rather than verified.
###
### Usage
###   Rscript 3d_extract_covariates.R probe     inventory only
###   Rscript 3d_extract_covariates.R tasks     print the array size and exit
###   Rscript 3d_extract_covariates.R           run one array task

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

ca_require(c("terra", "arrow"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("probe", "tasks")) args[1] else "run"

# Existing output is skipped by default. Pass overwrite to force a rerun, which
# is needed whenever a derivation rule changes rather than a filename, since
# the file is present but its contents are stale.
#   sbatch --export=ALL,CA_OVERWRITE=1 3d_extract_covariates.sub
OVERWRITE <- "overwrite" %in% args ||
  nzchar(Sys.getenv("CA_OVERWRITE"))

# ---------------------------------------------------------------------------
# PROBE
# ---------------------------------------------------------------------------

if (mode == "probe") {

  ca_stamp("stage3d_probe")

  s <- ca_probe_static()
  print(s[, c("layer", "status", "n_files", "res_x", "same_crs",
              "grid_aligned", "datatype", "na_flag", "cfg_nodata",
              "smp_min", "smp_median", "smp_max")])

  z <- ca_probe_zones()
  print(z[, c("layer", "status", "n_features", "same_crs", "cfg_field",
              "actual_field", "n_distinct", "n_missing", "max_nchar",
              "example")])

  ca_log("Check three things before running the array.")
  ca_log("  1. actual_field against cfg_field for both zone layers.")
  ca_log("  2. smp_min for an undeclared sentinel, for instance a large ",
         "negative on aspect where flat cells are coded rather than masked.")
  ca_log("  3. same_crs. FALSE is expected and handled, TRUE means no ",
         "transform is needed.")

  quit(save = "no")
}

# ---------------------------------------------------------------------------
# TASKS
# ---------------------------------------------------------------------------

tasks <- ca_static_task_table()

if (mode == "tasks") {
  ca_log("Active raster layers: ",
         paste(ca_static_layers$layer[ca_static_layers$active], collapse = ", "))
  ca_log("Active zone layers:   ",
         paste(ca_zone_layers$layer[ca_zone_layers$active], collapse = ", "))
  ca_log("Array size: 1-", nrow(tasks))
  print(table(tasks$layer, tasks$type))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------

ca_stamp("stage3d")

task_id <- ca_task_id()
if (is.na(task_id)) {
  stop("No task id. Submit as an array job or pass an integer argument.")
}

ca_log("Source static.  Threads ", ca_threads(),
       ".  Array size ", nrow(tasks),
       if (OVERWRITE) ".  OVERWRITE on" else "")

ca_run_static_task(task_id, overwrite = OVERWRITE)

ca_log("Stage 3d task ", task_id, " complete.")
