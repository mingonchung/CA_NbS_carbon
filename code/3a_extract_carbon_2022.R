### CA carbon revision pipeline
### Stage 3a. Carbon extraction, 2022 NCSDA plus eMapR and LEMMA
###
### Extracts the 2022 NCSDA flux stacks and the two legacy aboveground biomass
### products to the frozen forest point grid. All three are local downloads read
### the same way, which is why they share one script.
###
### NPP, NEP, and NBP are the operative outcomes here until an Almanac release
### carries equivalents. GPP is extracted as a cross-comparison against Almanac
### Carbon_GPP rather than as a separate outcome. CStocks_Live is an inactive
### scaffold, needed only if the aboveground versus total live divergence has to
### be quantified directly.
###
### The submitted pipeline inferred year from alphabetical file order and set
### every zero to NA. Both are retired. Years are parsed from filenames and
### zero is a legitimate value, which matters because NEP and NBP are signed.
###
### Usage
###   Rscript 3a_extract_carbon_2022.R probe     inventory only, no extraction
###   Rscript 3a_extract_carbon_2022.R tasks     print the array size and exit
###   Rscript 3a_extract_carbon_2022.R           run one array task

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
})

SOURCE <- "ncsda2022"
REGISTRY <- ca_ncsda_layers

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("probe", "tasks")) args[1] else "run"

if (mode == "probe") {
  ca_stamp("stage3a_2022_probe")
  out <- ca_probe_layers(REGISTRY, SOURCE)
  print(out[, c("layer", "status", "n_files", "first_year", "last_year",
                "gaps", "n_bands", "datatype", "na_flag", "crs_proj")])
  ca_log("Probe complete. eMapR and LEMMA CRS may differ from the grid. ",
         "Points are transformed into the raster CRS, rasters are never ",
         "reprojected, so a mismatch is handled rather than corrected.")
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

ca_stamp("stage3a_2022")

task_id <- ca_task_id()
if (is.na(task_id)) {
  stop("No task id. Submit as an array job or pass an integer argument.")
}

ca_log("Source ", SOURCE, ".  Threads ", ca_threads(),
       ".  Array size ", nrow(tasks))

ca_run_extract_task(REGISTRY, SOURCE, task_id)

ca_log("Stage 3a 2022 task ", task_id, " complete.")
