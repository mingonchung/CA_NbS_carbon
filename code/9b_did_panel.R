### CA carbon revision pipeline
### Stage 9b. DiD panel, one file per pairing and outcome
###
### WHAT THIS STAGE DOES
###
### Attaches carbon outcomes to the matched units written by 9a and writes the
### long panel stage 10 estimates on. One row per matched unit and year.
###
### THE ARRAY
###
### Pairing crossed with outcome, from ca_did_tasks(), filtered by
### ca_outcomes$scope so NBP appears for the two fire pairings only. 38 tasks.
###
### The submitted pipeline ran this as one pass per pairing over four carbon
### files (9_DiD_input_all.R L43). Seven outcomes and 36 years is 756 shard
### reads per pairing in one process, and the outcomes are independent, so the
### array is over the cross rather than over the pairing.
###
### THE LABEL COMES FROM THE REGISTRY, NOT FROM THE FILE ORDER
###
### 9_DiD_input_all.R L34 names outcomes by the position of their file in
### list.files(), so adding a product renames every column after it and nothing
### errors. Here the outcome is a task-table entry, its layer and source come
### from ca_outcome_layer(), its years from ca_did_years(), and the pairing,
### outcome, arm, and config are written to meta/did_panel_manifest.csv, which
### is what stage 10 reads.
###
### THE PANEL OPENS IN 1985, NOT 1990
###
### ca_did$panel_years is ca_years$study. The 1990 to 1994 cohorts have no
### pre-treatment period inside a panel that starts in 1990, so
### aggte(type = "dynamic") has nothing to test parallel trends against for
### them, and PA establishments dated 1986 to 1989 are dropped by att_gt() as
### already treated in the first period. Both are recovered by the five extra
### years. ca_years$analysis is untouched, since it governs which events may
### establish a cohort through ca_window_rules and that is settled at stages 5
### and 6.
###
### Per outcome the window is ca_did_years(), which is ca_layer_years()
### intersected with the panel. eMapR and LEMMA start in 1990 and NBP in 1986
### regardless, so the panel is ragged and the actual range written lands in
### the manifest rather than being assumed from the config.
###
### NO TRANSFORM
###
### The panel carries the converted value once. Stage 10 applies
### ca_log_transform() where ca_outcomes$log is TRUE. NPP, NEP, and NBP are
### signed and are non-log only, which is a property of the outcome and lives
### in the registry rather than in a column name.
###
### Usage
###   Rscript 9b_did_panel.R tasks     print the array and exit
###   Rscript 9b_did_panel.R           run one array task
###   Rscript 9b_did_panel.R archive   copy completed output to /projects
###
###   CA_OVERWRITE=1                   rebuild panels that already exist

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

ca_require(c("arrow", "dplyr"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("tasks", "archive")) args[1] else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

TASKS <- ca_did_tasks()

if (mode == "tasks") {
  TASKS$years <- vapply(TASKS$outcome, function(o) {
    y <- ca_did_years(o); sprintf("%d-%d", min(y), max(y))
  }, character(1))
  print(TASKS)
  cat("array size:", nrow(TASKS), "\n")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_did_dir(), "did_panel", subdir = "work")
  quit(save = "no")
}

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(TASKS)) {
  stop("Task id must be 1 to ", nrow(TASKS), ", got ", task)
}

PAIRING <- TASKS$pairing[task]
OUTCOME <- TASKS$outcome[task]
YEARS   <- ca_did_years(OUTCOME)
OL      <- ca_outcome_layer(OUTCOME)
ROW     <- ca_outcomes[match(OUTCOME, ca_outcomes$outcome), ]

ca_log("Stage 9b, task ", task, ".  pairing ", PAIRING, "  outcome ", OUTCOME)
ca_log("  source ", OL$source, "  layer ", OL$layer,
       "  years ", min(YEARS), " to ", max(YEARS), " (", length(YEARS), ")")

# SHARD PRE-FLIGHT.
#
# ca_read_layer() stops on a missing extract, and it stops inside the class
# loop after the reads that preceded it have already run. A year absent from
# disk is a stage 3 gap rather than a stage 9 fault, so it is detected up front
# and the window is trimmed. The trimmed range goes to the manifest, so a panel
# shorter than the registry claims is visible without reading the log.
have_year <- vapply(YEARS, function(y) {
  all(vapply(ca_lulc$label, function(cl) {
    file.exists(ca_extract_path(ca_extract_store(OL$layer), OL$layer, y, cl))
  }, logical(1)))
}, logical(1))

if (!all(have_year)) {
  ca_log("  WARNING extraction shards absent for: ",
         paste(YEARS[!have_year], collapse = ", "))
  YEARS <- YEARS[have_year]
  if (!length(YEARS)) {
    stop("No extraction shards on disk for layer ", OL$layer,
         ". Stage 3 has not run for this layer.")
  }
  ca_log("  window trimmed to ", min(YEARS), " to ", max(YEARS),
         " (", length(YEARS), " years)")
}

out_path <- ca_did_panel_path(PAIRING, OUTCOME)
if (file.exists(out_path) && !OVERWRITE) {
  ca_log("Present, not overwritten: ", basename(out_path))
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# UNITS
# ---------------------------------------------------------------------------

units_path <- ca_did_units_path(PAIRING)
if (!file.exists(units_path)) {
  stop("Missing units file: ", units_path,
       ". Run 9a_did_units.R before this stage.")
}

UNIT_COLS <- c("unit_id", "pair_id", "pixel_id", "lulc", "treat", "gvar",
               "cohort_year", "pa_year_known", "stratum", "analysis_group",
               "ecoregion_l3")

units <- as.data.frame(dplyr::collect(dplyr::select(
  arrow::open_dataset(units_path), dplyr::all_of(UNIT_COLS))))
units$lulc <- as.character(units$lulc)
units$pixel_id <- as.integer(units$pixel_id)

ca_log("  units ", format(nrow(units), big.mark = ","),
       ", pairs ", format(length(unique(units$pair_id)), big.mark = ","),
       ", distinct pixels ",
       format(length(unique(units$pixel_id)), big.mark = ","))


# ---------------------------------------------------------------------------
# BUILD
# ---------------------------------------------------------------------------
# Per class, because the extracts are sharded on class and the read has to name
# one. Values are read once per distinct pixel and then indexed out to the unit
# rows, so a control matched by three cohorts costs one read and three rows.

n_written <- 0L
n_cells   <- 0
n_missing <- 0
parts <- list()

for (cl in ca_lulc$label) {

  u <- units[units$lulc == cl, , drop = FALSE]
  if (!nrow(u)) {
    ca_log("  ", cl, ": no matched units")
    next
  }

  ids <- sort(unique(u$pixel_id))
  ca_log("  ", cl, ": ", format(nrow(u), big.mark = ","), " units, ",
         format(length(ids), big.mark = ","), " distinct pixels")

  m <- ca_read_outcome(OUTCOME, cl, pixel_id = ids, years = YEARS)

  # Expand to unit rows. match() rather than a join, since ids is the exact
  # row index of the matrix and the mapping is one-to-many by construction.
  rmap <- match(u$pixel_id, ids)
  v <- m[rmap, , drop = FALSE]
  rm(m); ca_gc()

  n_cells <- n_cells + length(v)
  keep <- which(!is.na(v))
  n_missing <- n_missing + (length(v) - length(keep))

  if (!length(keep)) {
    ca_log("    every value missing, nothing written for this class")
    rm(v, u, rmap); ca_gc()
    next
  }

  ri <- ((keep - 1L) %% nrow(v)) + 1L
  ci <- ((keep - 1L) %/% nrow(v)) + 1L

  d <- u[ri, , drop = FALSE]
  d$year  <- as.integer(YEARS[ci])
  d$value <- as.numeric(v[keep])
  rownames(d) <- NULL

  ca_log("    rows ", format(nrow(d), big.mark = ","), ", missing dropped ",
         format(length(v) - length(keep), big.mark = ","), " (",
         round(100 * (length(v) - length(keep)) / length(v), 2), " percent)")

  parts[[length(parts) + 1L]] <- d
  n_written <- n_written + nrow(d)
  rm(v, d, u, rmap, keep, ri, ci); ca_gc(paste0("after ", cl))
}

if (!n_written) {
  stop("No non-missing outcome values for ", PAIRING, " ", OUTCOME,
       ". Check that stage 3 wrote layer ", OL$layer, " for ",
       min(YEARS), " to ", max(YEARS), ".")
}

d <- do.call(rbind, parts)
rm(parts); ca_gc("panel assembled")

d <- d[order(d$unit_id, d$year, method = "radix"), , drop = FALSE]
rownames(d) <- NULL

for (v in c("lulc", "stratum", "analysis_group", "ecoregion_l3")) {
  d[[v]] <- factor(d[[v]])
}

arrow::write_parquet(arrow::as_arrow_table(d), out_path,
                     compression = ca_io$compression)
ca_log("Wrote ", basename(out_path), "  ",
       format(nrow(d), big.mark = ","), " rows  ",
       round(file.size(out_path) / 1e6, 1), " MB")


# ---------------------------------------------------------------------------
# MANIFEST
# ---------------------------------------------------------------------------
# One row per panel, appended under a lock. This is what stage 10 reads to know
# which file holds which outcome, replacing the positional csv.list[i] of the
# submitted pipeline.

acquire_lock <- function(lock_dir, timeout = 900, stale = 1800, interval = 0.5) {
  t0 <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) return(invisible(TRUE))
    mtime <- file.info(lock_dir)$mtime
    if (!is.na(mtime) &&
        difftime(Sys.time(), mtime, units = "secs") > stale) {
      warning("Removing stale lock: ", lock_dir)
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("Timed out acquiring lock: ", lock_dir)
    }
    Sys.sleep(interval + runif(1, 0, 0.3))
  }
}

sel <- ca_match_selected(PAIRING)
tr  <- d[d$treat == 1L, , drop = FALSE]

man <- data.frame(
  pairing        = PAIRING,
  outcome        = OUTCOME,
  label          = ROW$label,
  source         = OL$source,
  layer          = OL$layer,
  arm            = sel$arm,
  config         = sel$config,
  year_first     = min(d$year),
  year_last      = max(d$year),
  n_years        = length(YEARS),
  years_expected = length(ca_did_years(OUTCOME)),
  log_transform  = if (isTRUE(ROW$log)) ca_log_transform(OUTCOME) else "none",
  nonlog         = ROW$nonlog,
  signed         = ROW$signed,
  role           = ROW$role,
  n_rows         = nrow(d),
  n_units        = length(unique(d$unit_id)),
  n_pairs        = length(unique(d$pair_id)),
  n_treated_rows = nrow(tr),
  n_cohorts      = length(unique(tr$gvar)),
  pct_missing    = round(100 * n_missing / max(n_cells, 1), 3),
  value_min      = round(min(d$value), 4),
  value_median   = round(stats::median(d$value), 4),
  value_max      = round(max(d$value), 4),
  path           = out_path,
  written_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)

lock <- file.path(ca_did_dir(), ".lock_manifest")
acquire_lock(lock)
on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

mpath <- ca_did_manifest_path()
if (file.exists(mpath)) {
  old <- utils::read.csv(mpath, stringsAsFactors = FALSE)
  old <- old[!(old$pairing == PAIRING & old$outcome == OUTCOME), , drop = FALSE]
  # Column sets must agree, otherwise a manifest written by an older revision
  # of this script would silently rbind into a frame with shifted names.
  if (nrow(old) && !identical(sort(names(old)), sort(names(man)))) {
    stop("meta/", ca_did$manifest, " has a different column set from this ",
         "script. Delete it and rerun the array rather than appending.")
  }
  man <- rbind(old[, names(man), drop = FALSE], man)
}
man <- man[order(match(man$pairing, ca_pairings$pairing),
                 match(man$outcome, ca_outcomes$outcome)), , drop = FALSE]
utils::write.csv(man, mpath, row.names = FALSE)
unlink(lock, recursive = TRUE, force = TRUE)

ca_log("manifest updated, ", nrow(man), " of ", nrow(TASKS), " panels present")

ca_stamp(sprintf("9b_panel_%s_%s", PAIRING, OUTCOME))
ca_log("Stage 9b task ", task, " complete.  ", PAIRING, "  ", OUTCOME)
