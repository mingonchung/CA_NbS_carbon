### CA carbon revision pipeline
### Stage 3c. Management, protection, and ownership polygons
###
### Eight layers by three forest classes, 24 tasks.
###
###   thin_nto         CAL FIRE notices of timber operations
###   thin_thp         CAL FIRE timber harvest plans, historical and current
###   thin_facts_th    FACTS timber harvest
###   thin_facts_hfr   FACTS hazardous fuel treatments
###   offset           CARB compliance offset project boundaries
###   tribal           federally recognised tribal lands, exclusion mask
###   cpad             protected areas, every managing agency level
###   ownership        public ownership, private is the complement
###
### Output is long. One row per pixel per intersecting polygon, so a pixel with
### three recorded entries yields three rows. This is required by
### ca_thin_rules$single_event_only, which cannot be evaluated from a first-hit
### join, and it is why stage 3c files are not comparable to stage 3d files.
###
### Two things this stage does not do. It does not assign treatment groups,
### which happens at stage 6 from record_class and intensity. It does not drop
### records for being non_disturbing or Variable, per ca_config section 0.6,
### which retains every archival record through extraction. The only filter is
### the archival status flag, which separates work done from work planned and
### is a property of the record rather than an analytical choice.
###
### Run probe first. It resolves every layer and estimates the join cardinality
### from a point sample, which is the number that decides whether the array
### fits in memory. Output length is not knowable in advance from the polygon
### count, because overlapping records multiply.
###
### Usage
###   Rscript 3c_extract_vectors.R probe     inventory and cardinality
###   Rscript 3c_extract_vectors.R tasks     print the array size and exit
###   Rscript 3c_extract_vectors.R           run one array task
###   Rscript 3c_extract_vectors.R overwrite rerun a task whose file exists

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

OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

# ---------------------------------------------------------------------------
# PROBE
# ---------------------------------------------------------------------------

if (mode == "probe") {

  ca_stamp("stage3c_probe")

  xw <- ca_intensity_table()
  ca_log("Crosswalk: ", nrow(xw), " entries, ",
         sum(xw$origin == "knight"), " from Knight, ",
         sum(xw$origin == "alias_alias"), " wording variants, ",
         sum(xw$origin == "alias_assigned"), " assigned by analogy")

  out <- ca_probe_vector_layers()

  print(out[, c("layer", "status", "n_read", "n_kept", "n_invalid_fixed",
                "n_status_dropped", "pct_year_na", "yr_min", "yr_max",
                "pct_unresolved", "pct_variant", "pct_assigned")])
  print(out[, c("layer", "smp_pts", "smp_rows", "smp_pct_hit",
                "rows_per_hit", "max_per_hit", "est_rows_class",
                "est_rows_everg")])

  ca_log("est_rows_everg is the number that sizes the job. Output is long, ",
         "so a pixel inside several overlapping records produces several ",
         "rows, and the row count is not knowable from the polygon count.")
  ca_log("The Decid sample understates the ownership and cpad layers for ",
         "Evergreen, because national forest land is predominantly ",
         "evergreen, so treat est_rows_everg for those two as a floor.")
  ca_log("pct_unresolved above zero means an activity string is present in ",
         "the shapefile and absent from both the Knight tables and ",
         "activity_alias.csv. Those records carry no record_class and would ",
         "be silently dropped at stage 6.")

  quit(save = "no")
}

# ---------------------------------------------------------------------------
# TASKS
# ---------------------------------------------------------------------------

tasks <- ca_vector_task_table()

if (mode == "tasks") {
  ca_log("Active layers: ",
         paste(ca_vector_layers$layer[ca_vector_layers$active],
               collapse = ", "))
  ca_log("Array size: 1-", nrow(tasks))
  print(table(tasks$layer))
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------

ca_stamp("stage3c")

task_id <- ca_task_id()
if (is.na(task_id)) {
  stop("No task id. Submit as an array job or pass an integer argument.")
}

ca_log("Source vector.  Threads ", ca_threads(),
       ".  Array size ", nrow(tasks),
       if (OVERWRITE) ".  OVERWRITE on" else "")

ca_run_vector_task(task_id, overwrite = OVERWRITE)

ca_log("Stage 3c task ", task_id, " complete.")
