### CA carbon revision pipeline
### Stage 3b. Disturbance and vegetation screens, 2026 Wildland Almanac
###
### Structurally identical to 3a_extract_carbon_2026.R. Only the registry
### differs, because these layers are screens rather than outcomes.
###   Disturbance_TreeFrac  1986-2024  control verification, thinning detection
###   Disturbance_AGB       1986-2024  secondary screen
###   Veg_TreeFrac          1985-2025  continuous forest cover against NLCD
###
### These replace the record-based canopy loss assertion in the submitted
### version with a detection-based screen. Control pixels must show zero
### detected disturbance across all available years, applied at stage 6.
###
### MTBS runs through this same script under the mtbs selector. It is a
### year-keyed raster with the identical output schema, so it needs a registry
### row and a selector rather than a stage of its own. It is not part of 3c,
### which handles polygons with an event date.
###
### Usage
###   Rscript 3b_extract_screens_2026.R probe            screens, inventory
###   Rscript 3b_extract_screens_2026.R tasks            screens, array size
###   Rscript 3b_extract_screens_2026.R                  screens, one task
###   Rscript 3b_extract_screens_2026.R mtbs probe       MTBS, inventory
###   Rscript 3b_extract_screens_2026.R mtbs tasks       MTBS, array size
###   Rscript 3b_extract_screens_2026.R mtbs             MTBS, one task

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)

# The registry selector is read from any position, so mtbs and probe can be
# given in either order.
if ("mtbs" %in% args) {
  SOURCE   <- "mtbs"
  REGISTRY <- ca_mtbs_layers
} else {
  SOURCE   <- "almanac2026_screens"
  REGISTRY <- ca_screen_layers
}

mode <- if (any(args %in% c("probe", "tasks"))) {
  args[args %in% c("probe", "tasks")][1]
} else {
  "run"
}

if (mode == "probe") {
  ca_stamp(paste0("stage3b_probe_", SOURCE))
  out <- ca_probe_layers(REGISTRY, SOURCE)
  print(out[, c("layer", "status", "n_files", "first_year", "last_year",
                "gaps", "n_bands", "datatype", "na_flag", "scale_file")])
  if (SOURCE == "mtbs") {
    ca_log("Probe complete. Confirm that classes run 1 to 6 and that the ",
           "background is already NA, since ca_sources lists mtbs as ",
           "unverified.")
  } else {
    ca_log("Probe complete. Confirm whether the disturbance layers encode ",
           "magnitude or a binary flag before stage 6 sets the screen rule.")
  }
  quit(save = "no")
}

tasks <- ca_task_table(REGISTRY, SOURCE)

if (mode == "tasks") {
  ca_log("Active layers: ",
         paste(REGISTRY$layer[REGISTRY$active], collapse = ", "))
  ca_log("Array size: 1-", nrow(tasks))
  print(table(tasks$layer))
  quit(save = "no")
}

ca_stamp(paste0("stage3b_", SOURCE))

task_id <- ca_task_id()
if (is.na(task_id)) {
  stop("No task id. Submit as an array job or pass an integer argument.")
}

ca_log("Source ", SOURCE, ".  Threads ", ca_threads(),
       ".  Array size ", nrow(tasks))

ca_run_extract_task(REGISTRY, SOURCE, task_id)

ca_log("Stage 3b task ", task_id, " complete.")
