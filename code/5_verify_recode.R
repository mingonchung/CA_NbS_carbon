### CA carbon revision pipeline
### Stage 5 verification
###
### Sections
###
###   events       per-arm row and pixel counts, cohort year distributions,
###                undated shares, class distributions
###   offsets      join coverage against the 173-project issuance workbook
###   screen       detection rate by year against the measured 0.08 to 0.73
###                percent range, detections per pixel
###   exclusion    the number that decides whether stage 6 has workable sample
###                sizes. single_event_only composes across arms, so a pixel
###                thinned once AND burned once is excluded from both
###   all          every section
###
### Usage
###   Rscript 5_verify_recode.R events     [class]
###   Rscript 5_verify_recode.R exclusion  [class]
###   Rscript 5_verify_recode.R all        [class]

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
section <- if (length(args)) args[1] else "all"
CLASS <- if (length(args) > 1) args[2] else "Everg"

if (!CLASS %in% ca_lulc$label) stop("Unknown class: ", CLASS)

STAMP <- format(Sys.time(), "%Y%m%d")
ARMS <- c("thin", "fire", "offset", "pa", "own", "tribal")

report <- function(df, name) {
  p <- ca_meta(sprintf("verify_stage5_%s_%s_%s.csv", name, CLASS, STAMP))
  utils::write.csv(df, p, row.names = FALSE)
  ca_log("Wrote ", basename(p))
  print(utils::head(df, 40))
  invisible(p)
}

ca_event_path <- function(arm, lulc) {
  file.path(ca_work("recode", "events"),
            sprintf("events_%s_%s.parquet", arm, lulc))
}
ca_screen_path <- function(lulc) {
  file.path(ca_work("recode", "screen"), sprintf("screen_%s.parquet", lulc))
}

read_arm <- function(arm, cols = NULL) {
  p <- ca_event_path(arm, CLASS)
  if (!file.exists(p)) { warning("Missing arm file: ", p); return(NULL) }
  ds <- arrow::open_dataset(p)
  if (is.null(cols)) return(as.data.frame(dplyr::collect(ds)))
  as.data.frame(dplyr::collect(dplyr::select(ds, dplyr::all_of(cols))))
}

# ---------------------------------------------------------------------------
# SECTION 1. EVENTS
# ---------------------------------------------------------------------------

verify_events <- function() {

  ca_log("Section 1. Events, class ", CLASS)
  n_grid <- ca_lulc_n_points(CLASS)

  out <- list(); cohort <- list(); classes <- list()

  for (arm in ARMS) {
    d <- read_arm(arm)
    if (is.null(d)) next

    out[[arm]] <- data.frame(
      arm = arm,
      n_rows = nrow(d),
      n_pixels = length(unique(d$pixel_id)),
      pct_of_class = round(100 * length(unique(d$pixel_id)) / n_grid, 3),
      n_layers = length(unique(d$source_layer)),
      n_undated = sum(is.na(d$event_year)),
      pct_undated = round(100 * mean(is.na(d$event_year)), 2),
      # A layer with no year field yields min() on an empty set, which
      # returns Inf rather than erroring. Reporting Inf as a year is worse
      # than reporting nothing.
      year_min = if (all(is.na(d$event_year))) NA_integer_ else
        as.integer(min(d$event_year, na.rm = TRUE)),
      year_max = if (all(is.na(d$event_year))) NA_integer_ else
        as.integer(max(d$event_year, na.rm = TRUE)),
      n_outside_window = sum(
        !is.na(d$event_year) &
          (d$event_year < ca_window_rules$first_year |
             d$event_year > ca_window_rules$last_year)),
      n_future = sum(!is.na(d$event_year) &
                       d$event_year > ca_window_rules$last_year),
      rows_per_pixel_max = max(table(d$pixel_id)),
      stringsAsFactors = FALSE)

    ct <- as.data.frame(table(d$event_year), stringsAsFactors = FALSE)
    if (nrow(ct)) {
      names(ct) <- c("event_year", "n_rows")
      ct$arm <- arm
      cohort[[arm]] <- ct[, c("arm", "event_year", "n_rows")]
    }

    cv <- as.data.frame(table(d$class_value, useNA = "ifany"),
                        stringsAsFactors = FALSE)
    if (nrow(cv)) {
      names(cv) <- c("class_value", "n_rows")
      cv$arm <- arm
      cv$n_pixels <- vapply(cv$class_value, function(k)
        length(unique(d$pixel_id[which(d$class_value == k)])), numeric(1))
      classes[[arm]] <- cv[, c("arm", "class_value", "n_rows", "n_pixels")]
    }

    ca_log("  ", arm, ": ", format(nrow(d), big.mark = ","), " rows, ",
           format(length(unique(d$pixel_id)), big.mark = ","), " pixels")
    rm(d); ca_gc()
  }

  report(do.call(rbind, out), "events_summary")
  report(do.call(rbind, cohort), "events_cohort")
  report(do.call(rbind, classes), "events_class")

  # Thinning intensity eligibility, which is what ca_thin_is_eligible() will
  # act on at stage 6. Reported before the rule runs so the exclusion is
  # visible rather than inferred from a shrunken sample.
  d <- read_arm("thin", c("pixel_id", "class_value", "record_class",
                          "event_year"))
  if (!is.null(d)) {
    elig <- ca_thin_is_eligible(d$class_value)
    ca_log("Thinning rows by intensity eligibility:")
    print(table(eligible = elig, useNA = "ifany"))
    ca_log("Thinning rows by record_class:")
    print(table(d$record_class, useNA = "ifany"))
    ca_log("Pixels whose ONLY thinning record_class is wildfire, which have ",
           "fire evidence but no severity and therefore no arm to enter:")
    wf <- unique(d$pixel_id[d$record_class == "wildfire"])
    oth <- unique(d$pixel_id[d$record_class != "wildfire"])
    ca_log("  ", format(length(setdiff(wf, oth)), big.mark = ","),
           " of ", format(length(wf), big.mark = ","),
           " pixels carrying a wildfire record")
    rm(wf, oth)
    ca_log("Excluded intensity labels present: ",
           paste(sort(unique(d$class_value[!elig & !is.na(d$class_value)])),
                 collapse = ", "))
    rm(d, elig); ca_gc()
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# SECTION 2. OFFSETS
# ---------------------------------------------------------------------------

verify_offsets <- function() {

  ca_log("Section 2. Offset join, class ", CLASS)

  d <- read_arm("offset")
  if (is.null(d)) return(invisible(NULL))
  ob <- utils::read.csv(ca_offset_years_csv, stringsAsFactors = FALSE)

  ids_pix <- unique(d$class_value)
  ids_wb  <- unique(toupper(trimws(ob[[ca_offset_rules$join_key]])))

  res <- data.frame(
    workbook_projects = nrow(ob),
    workbook_ids = length(ids_wb),
    ids_in_polygons = length(ids_pix),
    ids_matched = sum(ids_pix %in% ids_wb),
    ids_unmatched = sum(!ids_pix %in% ids_wb),
    rows_missing_cohort = sum(is.na(d$event_year)),
    pct_rows_missing_cohort = round(100 * mean(is.na(d$event_year)), 3),
    cohort_min = suppressWarnings(min(d$event_year, na.rm = TRUE)),
    cohort_max = suppressWarnings(max(d$event_year, na.rm = TRUE)),
    stringsAsFactors = FALSE)
  report(res, "offset_join")

  if (res$ids_unmatched > 0) {
    ca_log("Unmatched ARB ids in the polygon layer: ",
           paste(utils::head(setdiff(ids_pix, ids_wb), 20), collapse = ", "))
  }

  # aux_year is carried but not used. Reporting how much it would have changed
  # makes the absorbing-treatment decision auditable rather than assumed.
  has_end <- !is.na(d$aux_year) & !is.na(d$event_year)
  if (any(has_end)) {
    span <- d$aux_year[has_end] - d$event_year[has_end]
    ca_log("If treatment ended at aux_year, the mean treated span would be ",
           round(mean(span), 1), " years against ",
           round(mean(max(ca_years$analysis) - d$event_year[has_end]), 1),
           " under the absorbing rule now in force.")
  }
  rm(d, ob); ca_gc()
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 3. SCREEN
# ---------------------------------------------------------------------------

verify_screen <- function() {

  ca_log("Section 3. Screen, class ", CLASS)

  p <- ca_screen_path(CLASS)
  if (!file.exists(p)) { warning("Missing screen file, run 5b"); return(NULL) }
  sc <- as.data.frame(arrow::read_parquet(p))
  n_grid <- ca_lulc_n_points(CLASS)

  by_year <- as.data.frame(table(sc$year), stringsAsFactors = FALSE)
  names(by_year) <- c("year", "n_detected")
  by_year$pct <- round(100 * by_year$n_detected / n_grid, 5)
  by_year$agb_also <- vapply(by_year$year, function(y)
    sum(!is.na(sc$agb_loss[sc$year == as.integer(y)]) &
          sc$agb_loss[sc$year == as.integer(y)] > 0), numeric(1))
  report(by_year, "screen_by_year")

  n_det_pix <- length(unique(sc$pixel_id))
  res <- data.frame(
    threshold = ca_screen_threshold("primary"),
    n_detections = nrow(sc),
    n_detected_pixels = n_det_pix,
    n_grid = n_grid,
    pct_ever_detected = round(100 * n_det_pix / n_grid, 3),
    pct_control_eligible = round(100 * (1 - n_det_pix / n_grid), 3),
    median_loss = round(stats::median(sc$treefrac_loss), 5),
    p95_loss = round(stats::quantile(sc$treefrac_loss, 0.95), 5),
    agb_agreement_pct = round(100 * mean(!is.na(sc$agb_loss) &
                                           sc$agb_loss > 0), 2),
    stringsAsFactors = FALSE)
  report(res, "screen_summary")

  # The two screen layers disagreeing is expected, since one tracks canopy and
  # the other biomass. The size of the disagreement is what stage 6 needs in
  # order to decide whether the secondary layer is a criterion or a note.
  ca_log("Detections where Disturbance_AGB also fired: ",
         res$agb_agreement_pct, " percent")

  rm(sc); ca_gc()
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 4. EXCLUSION ACCOUNTING
# ---------------------------------------------------------------------------
# single_event_only now applies to thinning and to fire, and it composes across
# arms. A pixel thinned once and burned once has a contaminated durability
# window whichever arm it would have entered, so it leaves both.
#
# This section produces the number that decides whether the four pairings still
# have workable sample sizes. If the fire arms fall too far, the response is
# not to relax the rule but to state a shorter durability window, since a
# ten-year window excludes far fewer pixels than a thirty-year one.

verify_exclusion <- function() {

  ca_log("Section 4. Exclusion accounting, class ", CLASS)
  n_grid <- ca_lulc_n_points(CLASS)

  th <- read_arm("thin", c("pixel_id", "event_year", "class_value",
                           "record_class"))
  fi <- read_arm("fire", c("pixel_id", "event_year", "class_value"))
  if (is.null(th) || is.null(fi)) return(invisible(NULL))

  # Distinct events, not rows. Two polygons recording the same entry year on
  # the same pixel is one event, and counting rows would inflate every
  # multi-event figure.
  th_ev <- unique(th[, c("pixel_id", "event_year")])
  fi_ev <- unique(fi[, c("pixel_id", "event_year")])

  th_n <- table(th_ev$pixel_id)
  fi_n <- table(fi_ev$pixel_id)

  th_pix <- as.integer(names(th_n))
  fi_pix <- as.integer(names(fi_n))

  th_multi <- th_pix[th_n > 1]
  fi_multi <- fi_pix[fi_n > 1]
  both <- intersect(th_pix, fi_pix)

  # Undated and ineligible-intensity thinning, which disqualify a pixel as a
  # control under ca_thin_rules even though they cannot enter treatment.
  th_undated <- unique(th$pixel_id[is.na(th$event_year)])
  th_inelig <- unique(th$pixel_id[!ca_thin_is_eligible(th$class_value)])

  # Fire classes 1 and 5, in a perimeter but not an analysable severity.
  fi_unclass <- unique(fi$pixel_id[
    fi$class_value %in% as.character(ca_fire_disturb_codes)])

  res <- data.frame(
    quantity = c(
      "grid pixels",
      "pixels with any thinning record",
      "pixels with any fire record",
      "thinning pixels with more than one event",
      "fire pixels with more than one event (reburn)",
      "pixels with both a thinning and a fire event",
      "thinning pixels with an undated record",
      "thinning pixels with an ineligible intensity label",
      "fire pixels whose only severity is class 1 or 5",
      "union of all single_event_only exclusions"),
    n = c(
      n_grid,
      length(th_pix),
      length(fi_pix),
      length(th_multi),
      length(fi_multi),
      length(both),
      length(th_undated),
      length(th_inelig),
      length(setdiff(fi_unclass, fi_pix[fi_n > 1])),
      length(unique(c(th_multi, fi_multi, both)))),
    stringsAsFactors = FALSE)
  res$pct_of_grid <- round(100 * res$n / n_grid, 3)

  # Share of each treated population lost, which is the figure that matters for
  # sample size rather than the share of the grid.
  res$pct_of_arm <- NA_real_
  res$pct_of_arm[4] <- round(100 * length(th_multi) / length(th_pix), 2)
  res$pct_of_arm[5] <- round(100 * length(fi_multi) / length(fi_pix), 2)
  res$pct_of_arm[6] <- round(100 * length(both) /
                               length(union(th_pix, fi_pix)), 2)

  report(res, "exclusion")

  # Reburn detail, since the fire arms are the ones at risk and the severity of
  # the first fire determines which arm loses the pixel.
  if (length(fi_multi)) {
    fm <- fi[fi$pixel_id %in% fi_multi, ]
    first <- fm[order(fm$pixel_id, fm$event_year), ]
    first <- first[!duplicated(first$pixel_id), ]
    rb <- as.data.frame(table(first$class_value), stringsAsFactors = FALSE)
    names(rb) <- c("first_severity_class", "n_reburned_pixels")
    rb$label <- ca_fire_severity$label[match(as.integer(rb$first_severity),
                                             ca_fire_severity$code)]
    rb$group <- ca_fire_severity$group[match(as.integer(rb$first_severity),
                                             ca_fire_severity$code)]
    report(rb, "reburn_by_first_severity")
    rm(fm, first); ca_gc()
  }

  ca_log("Reburn share of fire pixels: ", res$pct_of_arm[5], " percent. ",
         "This is the scope condition for the durability claim and belongs in ",
         "the Discussion with this number attached.")

  rm(th, fi, th_ev, fi_ev); ca_gc()
  invisible(res)
}

# ---------------------------------------------------------------------------

ca_log("Stage 5 verification, section '", section, "', class ", CLASS)

if (section %in% c("events", "all"))    { verify_events();    ca_gc("events") }
if (section %in% c("offsets", "all"))   { verify_offsets();   ca_gc("offsets") }
if (section %in% c("screen", "all"))    { verify_screen();    ca_gc("screen") }
if (section %in% c("exclusion", "all")) { verify_exclusion(); ca_gc("excl") }

if (!section %in% c("events", "offsets", "screen", "exclusion", "all")) {
  stop("Unknown section '", section, "'")
}

ca_stamp("5_verify_recode")
ca_log("Done")
