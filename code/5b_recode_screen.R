### CA carbon revision pipeline
### Stage 5b. Disturbance screen
###
### Emits one row per pixel per year in which the Almanac detected canopy loss
### above the derived threshold. Long rather than wide, because the measured
### detection rate is 0.08 to 0.73 percent per year, so Evergreen produces
### roughly ten million rows instead of the 3.7 billion a full pixel-year panel
### would carry. The absence of a row is the information.
###
### THRESHOLD
###
### Read from meta/screen_threshold.csv through ca_screen_threshold(), never
### hardcoded. The derived value is 0, meaning any detected loss disqualifies a
### control pixel.
###
### That is the submitted specification, and the diagnostic in section 4 of
### 4_verify_clean.R confirmed it against the alternative. 99.68 percent of
### reference pixel-years are exactly zero, retention at threshold 0 is 90.28
### percent, and raising the threshold to 0.10 buys 1.9 points. The objection to
### threshold 0, that a per-year false positive rate compounds over 39 years,
### assumed a noisy continuous field. Disturbance_TreeFrac is a sparse detection
### product with a 0.27 percent annual rate in reference pixels, so the
### compounding argument does not apply to it. The Methods text describing a 5
### percent canopy loss threshold does not describe this rule and needs
### correcting.
###
### SIGN
###
### Positive is loss. Negative values are canopy gain and never disqualify a
### pixel, which is why the test is value > threshold rather than abs(value).
###
### SECONDARY SCREEN
###
### Disturbance_AGB is read for the detected pixel-years only, as a second
### opinion rather than a second criterion. Stage 6 decides whether a detection
### in one layer and not the other counts. Reading it only where TreeFrac fired
### keeps the cost negligible.
###
### NO CUMULATIVE CRITERION
###
### Considered and rejected. Repeated entry is already caught by the stage 3c
### record and by single_event_only, and slow background decline from drought
### and beetle mortality is a common trend the DiD differences out while
### matching conditions on ecoregion, HUC8, and climate normals. See ca_screen.
###
### Usage
###   Rscript 5b_recode_screen.R tasks
###   Rscript 5b_recode_screen.R              run one array task
###   Rscript 5b_recode_screen.R archive

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

ca_screen_path <- function(lulc) {
  dir <- ca_work("recode", "screen")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("screen_%s.parquet", lulc))
}

if (mode == "tasks") {
  cat("Array size: ", nrow(ca_lulc), "\n", sep = "")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_work("recode", "screen"), "screen", subdir = "work")
  quit(save = "no")
}

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id must be 1 to ", nrow(ca_lulc), ", got ", task)
}
lulc <- ca_lulc$label[task]

out_path <- ca_screen_path(lulc)
if (file.exists(out_path) && !OVERWRITE) {
  ca_log("Present, skipping: ", basename(out_path))
  quit(save = "no")
}

# Errors rather than defaulting if the diagnostic has not run. A silent default
# is how an unresolved methodological question becomes a published number.
THR <- ca_screen_threshold("primary")
ca_log("Stage 5b screen, class ", lulc, ", threshold ", THR,
       " (", ca_unit_conversion$target[
         match(ca_screen$layer, ca_unit_conversion$layer)], ")")

yrs <- ca_screen$years
parts <- list()
tally <- list()

for (y in yrs) {

  # clean = TRUE applies ca_mask_sentinel first, so the -32767 pixel in 2024
  # becomes NA rather than entering as a -3.28 detection.
  r <- ca_read_layer("almanac2026", ca_screen$layer, y, lulc, clean = TRUE)

  hit <- which(!is.na(r$value) & r$value > THR)
  n_valid <- sum(!is.na(r$value))

  tally[[as.character(y)]] <- data.frame(
    year = y, n_valid = n_valid, n_detected = length(hit),
    pct_detected = round(100 * length(hit) / max(n_valid, 1), 5),
    n_gain = sum(!is.na(r$value) & r$value < 0),
    stringsAsFactors = FALSE)

  if (!length(hit)) { rm(r); ca_gc(); next }

  pid <- r$pixel_id[hit]
  d <- data.frame(pixel_id = as.integer(pid), year = as.integer(y),
                  treefrac_loss = as.numeric(r$value[hit]),
                  agb_loss = NA_real_, stringsAsFactors = FALSE)
  rm(r); ca_gc()

  # Secondary layer, read only where the primary fired. The pushdown is over a
  # few hundred thousand ids rather than the grid.
  ap <- ca_extract_path(ca_extract_store(ca_screen$secondary_layer),
                        ca_screen$secondary_layer, y, lulc)
  if (file.exists(ap)) {
    a <- ca_read_layer("almanac2026", ca_screen$secondary_layer, y, lulc,
                       pixel_id = pid, clean = TRUE)
    d$agb_loss <- a$value
    rm(a)
  }

  parts[[as.character(y)]] <- d
  ca_log("  ", y, "  detected ", format(length(hit), big.mark = ","),
         "  (", round(100 * length(hit) / max(n_valid, 1), 3), " percent)")
  rm(d, pid); ca_gc()
}

screen <- do.call(rbind, parts)
screen <- screen[order(screen$pixel_id, screen$year), ]

arrow::write_parquet(arrow::as_arrow_table(screen), out_path,
                     compression = ca_io$compression)

ca_log("Wrote ", basename(out_path), "  ",
       format(nrow(screen), big.mark = ","), " detections across ",
       format(length(unique(screen$pixel_id)), big.mark = ","), " pixels  ",
       round(file.size(out_path) / 1e6, 1), " MB")

n_grid <- ca_lulc_n_points(lulc)
ca_log("Control-eligible before any record screen: ",
       format(n_grid - length(unique(screen$pixel_id)), big.mark = ","),
       " of ", format(n_grid, big.mark = ","), " (",
       round(100 * (1 - length(unique(screen$pixel_id)) / n_grid), 2),
       " percent)")

tally_df <- do.call(rbind, tally)
tp <- ca_meta(sprintf("screen_detection_by_year_%s.csv", lulc))
utils::write.csv(tally_df, tp, row.names = FALSE)
ca_log("Detection rate by year: ", tp)

# Detections per pixel, which is the distribution stage 6 reads control
# eligibility off. Reported here so the shape is on the record before any rule
# consumes it.
per_pix <- table(table(screen$pixel_id))
ca_log("Detections per detected pixel:")
print(utils::head(per_pix, 10))

ca_stamp("5b_recode_screen")
ca_log("Done, class ", lulc)
