### CA carbon revision pipeline
### Stage 3a. Carbon extraction, 2026 Wildland Almanac
###
### Extracts Carbon_AGB, Carbon_GPP, and Fire_LCP band 6 canopy height to the
### frozen forest point grid. Carbon_NPP, Carbon_NEP, and Carbon_NBP are
### present in the registry as inactive scaffolds. When a future Almanac release
### carries them, set active to TRUE in ca_config.R and rerun. Nothing else
### changes, because the array size is derived from the registry.
###
### Values are written in native units. Unit conversion and the 0.47 carbon
### fraction are applied at stage 4 so the conversion lives in one auditable
### place rather than being baked into the extraction.
###
### Usage
###   Rscript 3a_extract_carbon_2026.R probe     inventory only, no extraction
###   Rscript 3a_extract_carbon_2026.R tasks     print the array size and exit
###   Rscript 3a_extract_carbon_2026.R           run one array task

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
})

SOURCE <- "almanac2026"
REGISTRY <- ca_almanac_layers

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("probe", "tasks")) args[1] else "run"

# ---------------------------------------------------------------------------
# PROBE
# ---------------------------------------------------------------------------
# Run this before the first extraction array. It reports band count, datatype,
# NA flag, and any scale and offset carried in each file, which is what closes
# the outstanding Almanac metadata items in ca_pending(). Inactive layers are
# probed too, so a future release is detected as soon as the files land.

if (mode == "probe") {
  ca_stamp("stage3a_2026_probe")
  out <- ca_probe_layers(REGISTRY, SOURCE)
  print(out[, c("layer", "status", "n_files", "first_year", "last_year",
                "gaps", "n_bands", "datatype", "na_flag", "scale_file")])
  ca_log("Probe complete. Set ca_sources nodata and scale from this table.")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# TASK TABLE
# ---------------------------------------------------------------------------

tasks <- ca_task_table(REGISTRY, SOURCE)

if (mode == "tasks") {
  ca_log("Active layers: ",
         paste(REGISTRY$layer[REGISTRY$active], collapse = ", "))
  ca_log("Array size: 1-", nrow(tasks))
  print(table(tasks$layer))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# RUN ONE TASK
# ---------------------------------------------------------------------------

ca_stamp("stage3a_2026")

task_id <- ca_task_id()
if (is.na(task_id)) {
  stop("No task id. Submit as an array job or pass an integer argument.")
}

ca_log("Source ", SOURCE, ".  Threads ", ca_threads(),
       ".  Array size ", nrow(tasks))

ca_run_extract_task(REGISTRY, SOURCE, task_id)

ca_log("Stage 3a 2026 task ", task_id, " complete.")
