### CA carbon revision pipeline
### Stage 5a. Treatment event table
###
### Normalises the eight stage 3c polygon layers and the MTBS rasters into one
### event schema. It classifies nothing. No pixel is excluded, no event is
### collapsed, no group code is assigned. Every analytical rule lives in
### ca_thin_rules, ca_fire_rules, ca_offset_rules, ca_cpad_agency, and
### ca_ownership_level, and every one of them is applied once, at stage 6.
###
### That division is the reason the submitted pipeline's exclusion counts could
### not be audited. Rules applied inside the recode leave no record of what
### they removed. Here the event table holds everything and stage 6 subtracts
### from it in the open.
###
### OUTPUT, one file per arm per class
###
###   events_thin_{lulc}    thin_nto, thin_thp, thin_facts_th, thin_facts_hfr
###   events_fire_{lulc}    MTBS, classes 1 to 5, class 6 dropped as the mask
###   events_offset_{lulc}  CARB compliance offset projects
###   events_pa_{lulc}      CPAD, all ten managing agency levels
###   events_own_{lulc}     Multi-Source Land Ownership, all seven levels
###   events_tribal_{lulc}  federally recognised tribal lands
###
### Per arm rather than one union, because ownership and MTBS alone would push
### the Evergreen table past a hundred million rows, and stage 6 reads arms
### selectively. Same reasoning as the shard-per-layer-year pattern in stage 3.
###
### SCHEMA
###
###   pixel_id       join key
###   arm            thin, fire, offset, pa, own, tribal
###   source_layer   the registry layer the row came from
###   poly_uid       provenance back to the source polygon, NA for MTBS
###   event_year     cohort year, NA where the source carries no usable date
###   class_value    the analytically meaningful class. Intensity label for
###                  thinning, severity code for fire, MNG_AG_LEV for PA,
###                  Own_Level for ownership, ARB id for offsets
###   key_value      raw activity or key field from 3c, provenance
###   key_value_2    second silvicultural prescription, thinning only
###   method_value   FACTS method field, thinning only
###   record_class   Knight classification, thinning only
###   xwalk_origin   knight, alias_alias, or alias_assigned
###   aux_year       offsets only, last_year from the issuance workbook, kept
###                  as provenance and for the SI sensitivity. NOT used to end
###                  treatment, see ca_offset_rules
###
### Usage
###   Rscript 5a_recode_treatment.R tasks
###   Rscript 5a_recode_treatment.R              run one array task
###   Rscript 5a_recode_treatment.R archive

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

ca_event_path <- function(arm, lulc) {
  dir <- ca_work("recode", "events")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("events_%s_%s.parquet", arm, lulc))
}

EVENT_COLS <- c("pixel_id", "arm", "source_layer", "poly_uid", "event_year",
                "class_value", "key_value", "key_value_2", "method_value",
                "record_class", "xwalk_origin", "aux_year")

# Every arm returns this shape, so the six files stack without a schema
# reconciliation step at stage 6.
ca_event_frame <- function(n) {
  data.frame(
    pixel_id = integer(n), arm = character(n), source_layer = character(n),
    poly_uid = rep(NA_integer_, n), event_year = rep(NA_integer_, n),
    class_value = rep(NA_character_, n), key_value = rep(NA_character_, n),
    key_value_2 = rep(NA_character_, n), method_value = rep(NA_character_, n),
    record_class = rep(NA_character_, n), xwalk_origin = rep(NA_character_, n),
    aux_year = rep(NA_integer_, n),
    stringsAsFactors = FALSE)
}

ca_write_events <- function(df, arm, lulc) {
  miss <- setdiff(EVENT_COLS, names(df))
  for (m in miss) df[[m]] <- NA
  df <- df[, EVENT_COLS]
  df$pixel_id <- as.integer(df$pixel_id)
  df$event_year <- as.integer(df$event_year)
  df$aux_year <- as.integer(df$aux_year)
  for (ch in c("arm", "source_layer", "class_value", "key_value",
               "key_value_2", "method_value", "record_class",
               "xwalk_origin")) {
    df[[ch]] <- as.character(df[[ch]])
  }
  p <- ca_event_path(arm, lulc)
  arrow::write_parquet(arrow::as_arrow_table(df), p,
                       compression = ca_io$compression)
  ca_log("  wrote ", basename(p), "  ", format(nrow(df), big.mark = ","),
         " rows  ", round(file.size(p) / 1e6, 1), " MB")
  invisible(p)
}

# ---------------------------------------------------------------------------

if (mode == "tasks") {
  cat("Array size: ", nrow(ca_lulc), "\n", sep = "")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_work("recode", "events"), "treat", subdir = "work")
  quit(save = "no")
}

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id must be 1 to ", nrow(ca_lulc), ", got ", task)
}
lulc <- ca_lulc$label[task]
ca_log("Stage 5a treatment events, class ", lulc)

read_vector <- function(layer) {
  p <- ca_vector_path(layer, lulc)
  if (!file.exists(p)) stop("Missing stage 3c layer: ", p)
  as.data.frame(arrow::read_parquet(p))
}

# ---------------------------------------------------------------------------
# THINNING
# ---------------------------------------------------------------------------
# class_value is the intensity label, verbatim. Low, Medium, High, Variable,
# Unknown, and the compound Variable labels all survive to stage 6, where
# ca_thin_is_eligible() decides which enter an arm. Collapsing them here to a
# numeric code is what made the submitted exclusions unauditable.

arm <- "thin"
if (!file.exists(ca_event_path(arm, lulc)) || OVERWRITE) {

  thin_layers <- ca_vector_layers$layer[grepl("^thin_", ca_vector_layers$layer)]
  parts <- list()

  for (lyr in thin_layers) {
    v <- read_vector(lyr)
    d <- ca_event_frame(nrow(v))
    d$pixel_id     <- v$pixel_id
    d$arm          <- arm
    d$source_layer <- lyr
    d$poly_uid     <- v$poly_uid
    d$event_year   <- v$event_year
    d$class_value  <- v$intensity
    d$key_value    <- v$key_value
    d$key_value_2  <- v$key_value_2
    d$method_value <- v$method_value
    d$record_class <- v$record_class
    d$xwalk_origin <- v$xwalk_origin
    parts[[lyr]] <- d
    ca_log("  ", lyr, ": ", format(nrow(v), big.mark = ","), " rows, ",
           format(length(unique(v$pixel_id)), big.mark = ","), " pixels, ",
           sum(is.na(v$event_year)), " undated")
    rm(v, d); ca_gc()
  }

  ca_write_events(do.call(rbind, parts), arm, lulc)
  rm(parts); ca_gc()
} else ca_log("  thin present, skipping")

# ---------------------------------------------------------------------------
# WILDFIRE
# ---------------------------------------------------------------------------
# MTBS is a categorical raster, so it is read through ca_read_categorical()
# rather than ca_read_layer(), which would demand a unit conversion that a
# class code does not have.
#
# Class 6 is the non-processing mask and produces no event. Classes 1 and 5 do
# produce events, because a pixel inside a perimeter is not an undisturbed
# control even when its severity is not analysable. ca_fire_rules decides what
# happens to them at stage 6.

arm <- "fire"
if (!file.exists(ca_event_path(arm, lulc)) || OVERWRITE) {

  yrs <- ca_mtbs_layers$first_year:ca_mtbs_layers$last_year
  parts <- list()

  for (y in yrs) {
    r <- ca_read_categorical("MTBS", y, lulc)
    keep <- which(!is.na(r$value) & r$value != ca_fire_nodata_code)
    if (!length(keep)) { rm(r); next }
    d <- ca_event_frame(length(keep))
    d$pixel_id    <- r$pixel_id[keep]
    d$arm         <- arm
    d$source_layer <- "MTBS"
    d$event_year  <- y
    d$class_value <- as.character(r$value[keep])
    parts[[as.character(y)]] <- d
    rm(r, d, keep); ca_gc()
  }

  fire <- do.call(rbind, parts)
  ca_log("  MTBS: ", format(nrow(fire), big.mark = ","), " events, ",
         format(length(unique(fire$pixel_id)), big.mark = ","), " pixels")
  print(table(fire$class_value, dnn = "severity class"))
  ca_write_events(fire, arm, lulc)
  rm(fire, parts); ca_gc()
} else ca_log("  fire present, skipping")

# ---------------------------------------------------------------------------
# OFFSETS
# ---------------------------------------------------------------------------
# The offset layer carries no date field in 3c, so the cohort year comes from
# the parsed issuance workbook, joined on ARB id. start_year is the Reporting
# Period 1 start date.
#
# last_year travels as aux_year and is not used to end treatment. Absence of
# issuance in a year is a reporting gap, not a terminated project, and CARB
# forest protocols carry hundred-year permanence obligations. See
# ca_offset_rules.

arm <- "offset"
if (!file.exists(ca_event_path(arm, lulc)) || OVERWRITE) {

  v <- read_vector("offset")
  ob <- utils::read.csv(ca_offset_years_csv, stringsAsFactors = FALSE)

  key_v <- toupper(trimws(as.character(v$key_value)))
  key_o <- toupper(trimws(as.character(ob[[ca_offset_rules$join_key]])))
  j <- match(key_v, key_o)

  n_unmatched <- sum(is.na(j))
  ca_log("  offset: ", format(nrow(v), big.mark = ","), " rows, ",
         length(unique(key_v)), " distinct ARB ids, ",
         sum(!is.na(unique(match(unique(key_v), key_o)))),
         " matched to the workbook's ", nrow(ob), " projects")
  if (n_unmatched) {
    warning("Offset rows with no workbook match: ", n_unmatched,
            ". Unmatched ids: ",
            paste(utils::head(unique(key_v[is.na(j)]), 10), collapse = ", "))
  }

  d <- ca_event_frame(nrow(v))
  d$pixel_id     <- v$pixel_id
  d$arm          <- arm
  d$source_layer <- "offset"
  d$poly_uid     <- v$poly_uid
  d$event_year   <- ob[[ca_offset_rules$cohort_field]][j]
  d$class_value  <- key_v
  d$key_value    <- v$key_value
  d$aux_year     <- ob$last_year[j]

  ca_write_events(d, arm, lulc)
  rm(v, ob, d); ca_gc()
} else ca_log("  offset present, skipping")

# ---------------------------------------------------------------------------
# PROTECTED AREAS, OWNERSHIP, TRIBAL
# ---------------------------------------------------------------------------
# All three are carried at full label resolution. A pixel intersecting several
# CPAD superunits produces several rows, which is what lets stage 6 apply a
# collapse rule rather than inherit whatever a first-hit join happened to pick.
#
# CPAD YR_EST of zero was already written as NA at stage 3c and stays NA here.
# The submitted pipeline's ifelse(YR_EST < 1986, 1985, YR_EST) asserted that
# every unknown-year protected area is old, which affects 73.5 percent of
# superunits. Assigning NA to the 1985 cohort is a stage 6 decision, and the
# sensitivity that drops those units depends on the NA surviving to that point.

for (spec in list(
  list(arm = "pa",     layer = "cpad",      class = "key_value"),
  list(arm = "own",    layer = "ownership", class = "key_value"),
  list(arm = "tribal", layer = "tribal",    class = NA))) {

  if (file.exists(ca_event_path(spec$arm, lulc)) && !OVERWRITE) {
    ca_log("  ", spec$arm, " present, skipping"); next
  }

  v <- read_vector(spec$layer)
  d <- ca_event_frame(nrow(v))
  d$pixel_id     <- v$pixel_id
  d$arm          <- spec$arm
  d$source_layer <- spec$layer
  d$poly_uid     <- v$poly_uid
  d$event_year   <- v$event_year
  d$key_value    <- v$key_value
  if (!is.na(spec$class)) d$class_value <- as.character(v[[spec$class]])

  n_pix <- length(unique(v$pixel_id))
  ca_log("  ", spec$layer, ": ", format(nrow(v), big.mark = ","), " rows, ",
         format(n_pix, big.mark = ","), " pixels, ",
         round(100 * mean(is.na(v$event_year)), 2), " percent undated")
  if (!is.na(spec$class)) {
    print(table(d$class_value, useNA = "ifany", dnn = spec$layer))
  }

  ca_write_events(d, spec$arm, lulc)
  rm(v, d); ca_gc()
}

ca_stamp("5a_recode_treatment")
ca_log("Done, class ", lulc)
