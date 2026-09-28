### CA carbon revision pipeline
### Stage 4b. Matching covariate table
###
### The one materialised product of stage 4. Reads the nine stage 3d layers,
### derives northness and eastness from aspect, and writes one row per pixel.
### Everything else in stage 4 is a read-time operation in ca_clean.R.
###
### Output columns
###
###   pixel_id            terra cell index of the NLCD 2001 grid
###   lulc                Decid, Everg, Mixed. Exact matching stratum
###   elevation           m
###   slope               degrees
###   aspect              degrees from north, retained for the record only
###   northness           cos(aspect), 1 due north, 0 where flat
###   eastness            sin(aspect), 1 due east, 0 where flat
###   ppt_normal          mm/yr, PRISM 1991-2020
###   tmean_normal        degrees C, PRISM 1991-2020
###   pop_density         persons/km2, WorldPop 2000-2002 mean
###   city_travel_time    minutes, settlement class 11
###   ecoregion_l3        EPA Level III, character, width 2. Exact stratum
###   huc8                WBD HUC8, character, width 8. Exact stratum
###
### aspect is kept alongside its two components because it costs one column and
### it is the only way to audit the decomposition after the fact without
### rerunning stage 3d.
###
### Nothing is filtered here. A pixel with a missing covariate stays in the
### table and is dropped at stage 8, where the covariate set for the arm being
### run is known and the loss can be counted against that arm rather than
### against all of them.
###
### Usage
###   Rscript 4b_clean_covariates.R tasks     print the array size and exit
###   Rscript 4b_clean_covariates.R           run one array task
###   Rscript 4b_clean_covariates.R archive   copy completed output to /projects

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("tasks", "archive")) args[1] else "run"

OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

ca_covariate_path <- function(lulc) {
  dir <- ca_work("clean", "covariates")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("covariates_%s.parquet", lulc))
}

# ---------------------------------------------------------------------------

if (mode == "tasks") {
  cat("Array size: ", nrow(ca_lulc), "\n", sep = "")
  cat("Tasks: ", paste(ca_lulc$label, collapse = ", "), "\n", sep = "")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_work("clean", "covariates"), "covariates", subdir = "work")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id must be 1 to ", nrow(ca_lulc), ", got ", task)
}

lulc <- ca_lulc$label[task]
out_path <- ca_covariate_path(lulc)

if (file.exists(out_path) && !OVERWRITE) {
  ca_log("Present, skipping: ", basename(out_path))
  quit(save = "no")
}

ca_log("Stage 4b covariates, class ", lulc)

# The grid partition, not one of the extracts, defines the row set. Reading the
# row set from an extract would silently inherit whatever that extraction
# dropped, and a covariate table shorter than the grid is very hard to notice
# downstream.
grid <- ca_grid_points(lulc)
pid <- grid$pixel_id
rm(grid); gc()

n <- length(pid)
ca_log("Grid rows: ", format(n, big.mark = ","))

out <- data.frame(pixel_id = pid, lulc = lulc, stringsAsFactors = FALSE)

na_report <- data.frame(layer = character(0), n_na = integer(0),
                        pct_na = numeric(0), stringsAsFactors = FALSE)

record_na <- function(layer, v) {
  n_na <- sum(is.na(v))
  data.frame(layer = layer, n_na = n_na,
             pct_na = round(100 * n_na / length(v), 4),
             stringsAsFactors = FALSE)
}

# Continuous rasters -------------------------------------------------------

for (lyr in ca_static_layers$layer) {
  r <- ca_read_static(lyr, lulc, pixel_id = pid)
  out[[lyr]] <- r$value
  na_report <- rbind(na_report, record_na(lyr, r$value))
  ca_log("  ", lyr, "  na ",
         format(sum(is.na(r$value)), big.mark = ","))
  rm(r); gc()
}

# Aspect decomposition -----------------------------------------------------

d <- ca_derive_aspect(out$aspect)
out$northness <- d$northness
out$eastness  <- d$eastness

ca_log("  aspect decomposed, flat cells ", format(d$n_flat, big.mark = ","),
       " (", round(100 * d$n_flat / n, 3), " percent) set to 0, 0")

na_report <- rbind(na_report,
                   record_na("northness", out$northness),
                   record_na("eastness", out$eastness))

# Zone joins ---------------------------------------------------------------

for (lyr in ca_zone_layers$layer) {
  r <- ca_read_static(lyr, lulc, pixel_id = pid)
  out[[lyr]] <- r$value
  na_report <- rbind(na_report, record_na(lyr, r$value))
  ca_log("  ", lyr, "  na ", format(sum(is.na(r$value)), big.mark = ","),
         "  distinct ", length(unique(r$value[!is.na(r$value)])))
  rm(r); gc()
}

# Gates --------------------------------------------------------------------
# Three things are worth failing on rather than writing and discovering later.

if (nrow(out) != n) {
  stop("Row count drift: ", nrow(out), " against grid ", n)
}
if (anyDuplicated(out$pixel_id)) {
  stop("Duplicate pixel_id in the covariate table")
}

rng <- range(c(out$northness, out$eastness), na.rm = TRUE)
if (rng[1] < -1.0001 || rng[2] > 1.0001) {
  stop("northness or eastness outside [-1, 1]: ", paste(rng, collapse = " to "))
}

# A zone code that lost its leading zero is a wrong exact stratum, not a
# missing one, so width is checked rather than assumed.
for (k in seq_len(nrow(ca_zone_layers))) {
  lyr <- ca_zone_layers$layer[k]
  v <- out[[lyr]]
  v <- v[!is.na(v)]
  w <- unique(nchar(v))
  ca_log("  ", lyr, " code widths: ", paste(sort(w), collapse = ", "))
  if (length(w) > 1) {
    warning("Ragged code width in ", lyr,
            ". Expected fixed width from the stage 3d padding rule.")
  }
}

# Column order -------------------------------------------------------------

col_order <- c("pixel_id", "lulc",
               "elevation", "slope", "aspect", "northness", "eastness",
               "ppt_normal", "tmean_normal", "pop_density", "city_travel_time",
               "ecoregion_l3", "huc8")
missing_cols <- setdiff(col_order, names(out))
if (length(missing_cols)) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}
out <- out[, col_order]

# Write --------------------------------------------------------------------

arrow::write_parquet(arrow::as_arrow_table(out), out_path,
                     compression = ca_io$compression)

ca_log("Wrote ", basename(out_path), "  ",
       format(nrow(out), big.mark = ","), " rows  ",
       round(file.size(out_path) / 1e6, 1), " MB")

na_path <- ca_meta(sprintf("clean_4b_na_%s.csv", lulc))
utils::write.csv(na_report, na_path, row.names = FALSE)
ca_log("NA report: ", na_path)

# Complete-case rate under each matching arm, reported rather than applied.
for (arm in c("continuous_baseline", "continuous_augmented")) {
  vars <- intersect(ca_match[[arm]], names(out))
  if (!length(vars)) next
  ok <- stats::complete.cases(out[, c(vars, ca_match$exact)])
  ca_log("Complete cases, ", arm, " plus exact: ",
         format(sum(ok), big.mark = ","), " (",
         round(100 * mean(ok), 2), " percent)")
}

ca_stamp("4b_clean_covariates")
ca_log("Done, class ", lulc)
