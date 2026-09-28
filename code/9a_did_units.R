### CA carbon revision pipeline
### Stage 9a. Matched units, one row per matched pixel and subclass
###
### WHAT THIS STAGE DOES
###
### Reads the selected matched dataset for each pairing, restores the
### attributes the DiD needs from their authoritative source, and writes one
### units file per pairing. It reads no carbon. Outcomes are attached at 9b.
###
### It replaces the input half of 9_DiD_input_all.R, which merged the matched
### CSV against four carbon panels in one pass keyed on positional columns
### (L14, L57, L69). The split exists because the units are small and cheap and
### the outcomes are large and slow, so the decisions worth checking by eye
### should not be buried inside a six-hour read.
###
### FOUR THINGS ARE FIXED HERE
###
### 1. THE UNIT. A control pixel matched by two cohorts appears twice under two
###    subclasses. unit_id is one integer per matched row so att_gt() and
###    etwfe() do not see duplicate unit-year observations. pixel_id rides
###    along so the reuse is visible and clusterable.
###
### 2. THE COHORT YEAR. Stage 8 writes max(window) + 1. For thinning the window
###    ends at c-2 through ca_pretrt$event_drop, so that column is one year
###    early on TP_UP and TNP_UNP. cohort_year is re-read from the stage 6
###    groups file and the difference against stage 8 is reported per pairing.
###
### 3. THE STRATUM. Severity and intensity live on the treated side only. The
###    treated member's stratum is propagated to its control so a stage 10
###    subset keeps the pair rather than half of it.
###
### 4. gvar. Controls take 0. Treated take cohort_year, or pa_est_year for
###    UP_UNP, through ca_did_gvar(). Unknown CPAD years take the floor and
###    carry pa_year_known so the sensitivity is a filter, not a second panel.
###
### Usage
###   Rscript 9a_did_units.R tasks     print the pairings and exit
###   Rscript 9a_did_units.R           build every pairing
###   Rscript 9a_did_units.R archive   copy completed output to /projects
###
###   CA_PAIRINGS=UP_UNP,TP_UP         restrict to a subset
###   CA_OVERWRITE=1                   rebuild files that already exist

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("tasks", "archive")) args[1] else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

PAIRINGS <- ca_pairings$pairing
sel_env <- Sys.getenv("CA_PAIRINGS")
if (nzchar(sel_env)) {
  want <- trimws(strsplit(sel_env, "[,:]")[[1]])
  bad <- setdiff(want, PAIRINGS)
  if (length(bad)) stop("Unknown pairing in CA_PAIRINGS: ",
                        paste(bad, collapse = ", "))
  PAIRINGS <- want
}

if (mode == "tasks") {
  for (p in PAIRINGS) {
    s <- try(ca_match_selected(p), silent = TRUE)
    cat(sprintf("%-8s  %s\n", p,
                if (inherits(s, "try-error")) "NO SELECTION" else
                  paste(s$arm, s$config)))
  }
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_did_dir(), "did_panel", subdir = "work")
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# READS
# ---------------------------------------------------------------------------

# Stage 8 writes one file per config, and a level-subset run writes one file
# per range with the range in the name. Both are read by pattern rather than by
# a single constructed path, because a config that had to be split by CA_LEVELS
# is complete only when its parts are stacked.
matched_files <- function(pairing, arm, config) {
  pat <- sprintf("^matched_%s_%s_%s(_lev[0-9]+-[0-9]+)?\\.parquet$",
                 pairing, arm, config)
  f <- list.files(ca_match_dir(pairing), pattern = pat, full.names = TRUE)
  if (!length(f)) {
    stop("No matched file for ", pairing, " ", arm, " ", config, " under ",
         ca_match_dir(pairing), ". Stage 8 has not written this ",
         "specification, or ca_match$selected names a config that was never ",
         "run for this pairing.")
  }
  sort(f)
}

# Column set. distance is absent for the Mahalanobis configs and stratum is
# absent where stage 6 assigned none, so anything missing is filled rather than
# dropped, the same rule read_side() uses at stage 8.
MATCH_COLS <- c("pixel_id", "lulc", "treat", "subclass", "weights", "distance",
                "analysis_group", "stratum", "ecoregion_l3", "huc8",
                "cohort_year", "win_key")

read_matched <- function(files) {
  parts <- lapply(files, function(p) {
    ds <- arrow::open_dataset(p)
    have <- intersect(MATCH_COLS, names(ds))
    if (!all(c("pixel_id", "treat", "subclass") %in% have)) {
      stop("Matched file is missing a key column: ", basename(p))
    }
    d <- as.data.frame(dplyr::collect(
      dplyr::select(ds, dplyr::all_of(have))))
    for (v in setdiff(MATCH_COLS, have)) d[[v]] <- rep(NA, nrow(d))
    d[, MATCH_COLS, drop = FALSE]
  })
  d <- do.call(rbind, parts)
  for (v in c("lulc", "subclass", "analysis_group", "stratum",
              "ecoregion_l3", "huc8", "win_key")) {
    if (is.factor(d[[v]])) d[[v]] <- as.character(d[[v]])
  }
  d$pixel_id <- as.integer(d$pixel_id)
  d$treat    <- as.integer(d$treat)
  # Stage 8's column, kept under a name that cannot be mistaken for the event.
  names(d)[names(d) == "cohort_year"] <- "level_year_m8"
  d
}

# Stage 6 is the authority on cohort_year, pa_est_year, stratum, and
# analysis_group. One row per pixel across the three class partitions.
GROUP_COLS <- c("pixel_id", "analysis_group", "stratum", "eligibility",
                "cohort_year", "pa_est_year", "first_event_year", "n_events")

read_group_attrs <- function(ids) {
  ids <- as.integer(ids)
  parts <- lapply(ca_lulc$label, function(cl) {
    p <- ca_group_path(cl)
    if (!file.exists(p)) stop("Missing stage 6 groups file: ", p)
    ds <- arrow::open_dataset(p)
    have <- intersect(GROUP_COLS, names(ds))
    d <- as.data.frame(dplyr::collect(dplyr::select(
      dplyr::filter(ds, pixel_id %in% ids), dplyr::all_of(have))))
    for (v in setdiff(GROUP_COLS, have)) d[[v]] <- rep(NA, nrow(d))
    d[, GROUP_COLS, drop = FALSE]
  })
  d <- do.call(rbind, parts)
  for (v in c("analysis_group", "stratum", "eligibility")) {
    if (is.factor(d[[v]])) d[[v]] <- as.character(d[[v]])
  }
  d
}


# ---------------------------------------------------------------------------
# BUILD ONE PAIRING
# ---------------------------------------------------------------------------

build_units <- function(pairing) {

  sel <- ca_match_selected(pairing)
  out_path <- ca_did_units_path(pairing)

  if (file.exists(out_path) && !OVERWRITE) {
    ca_log("Present, not rebuilt: ", basename(out_path))
    return(NULL)
  }

  files <- matched_files(pairing, sel$arm, sel$config)
  ca_log(pairing, ", arm ", sel$arm, ", config ", sel$config, ", ",
         length(files), " matched file(s)")

  u <- read_matched(files)
  ca_log("  matched rows read: ", format(nrow(u), big.mark = ","))

  # A matched dataset is pairs. Anything else means two level ranges overlapped
  # or a file was written twice, and both are silent unless counted.
  tab <- table(u$subclass)
  if (any(tab != 2L)) {
    n_bad <- sum(tab != 2L)
    ca_log("  WARNING subclasses not holding exactly two rows: ",
           format(n_bad, big.mark = ","), " of ",
           format(length(tab), big.mark = ","))
    print(utils::head(sort(tab[tab != 2L], decreasing = TRUE), 10))
  }
  if (anyDuplicated(paste(u$subclass, u$treat))) {
    stop("A subclass holds two rows on the same side of the match in ",
         pairing, ". Two level ranges have been stacked twice.")
  }

  # Authoritative attributes.
  g <- read_group_attrs(unique(u$pixel_id))
  k <- match(u$pixel_id, g$pixel_id)
  if (anyNA(k)) {
    stop(sum(is.na(k)), " matched pixels absent from the stage 6 groups files")
  }
  u$cohort_year_g   <- as.integer(g$cohort_year[k])
  u$pa_est_year     <- as.integer(g$pa_est_year[k])
  u$stratum_pixel   <- g$stratum[k]
  u$analysis_group  <- g$analysis_group[k]
  u$eligibility     <- g$eligibility[k]
  rm(g, k)

  # THE OFF-BY-ONE AUDIT. Reported, never patched. Stage 8's level_year_m8 is
  # max(window) + 1, which equals the event for fire, offsets, and protection
  # dated 1986 or later, and the event minus one for thinning.
  tr_i <- which(u$treat == 1L)
  d8 <- u$cohort_year_g[tr_i] - u$level_year_m8[tr_i]
  d8_valid <- d8[!is.na(d8)]
  shift <- if (length(d8_valid) > 0L) {
    t8 <- table(d8_valid)
    as.integer(names(t8)[which.max(t8)])
  } else {
    NA_integer_
  }
  ca_log("  stage 6 cohort_year minus stage 8 level year, modal shift: ",
         shift, " on ", format(length(tr_i), big.mark = ","), " treated rows")
  if (length(d8) && length(unique(d8[!is.na(d8)])) > 1L) {
    print(table(d8, useNA = "ifany"))
  }

  # Pair-level attributes. The treated member defines the stratum, the cohort,
  # and the establishment year for its control.
  tr <- u[tr_i, , drop = FALSE]
  j <- match(u$subclass, tr$subclass)
  if (anyNA(j)) {
    n <- sum(is.na(j))
    ca_log("  WARNING subclasses with no treated member, rows dropped: ",
           format(n, big.mark = ","))
    u <- u[!is.na(j), , drop = FALSE]
    j <- j[!is.na(j)]
  }
  u$cohort_year <- tr$cohort_year_g[j]
  u$stratum     <- tr$stratum_pixel[j]
  u$pa_est_treat <- tr$pa_est_year[j]
  u$pa_year_known <- !is.na(tr$pa_est_year[j])
  rm(tr, j)

  # gvar, from config so the verifier tests the rule the script applied.
  u$gvar <- ca_did_gvar(pairing, u$treat, u$cohort_year, u$pa_est_treat)

  bad_g <- u$treat == 1L & (is.na(u$gvar) | u$gvar <= 0L)
  if (any(bad_g)) {
    ca_log("  treated rows with no usable gvar, pairs dropped: ",
           format(sum(bad_g), big.mark = ","))
    drop_sc <- unique(u$subclass[bad_g])
    u <- u[!u$subclass %in% drop_sc, , drop = FALSE]
  }

  rng <- ca_did_gvar_range(pairing)
  out_of_range <- u$treat == 1L & (u$gvar < rng[1] | u$gvar > rng[2])
  if (any(out_of_range)) {
    ca_log("  WARNING treated gvar outside [", rng[1], ", ", rng[2], "]: ",
           format(sum(out_of_range), big.mark = ","))
    print(table(u$gvar[out_of_range]))
  }

  # pair_id. subclass strings are level:index and run to twenty characters, and
  # the panel repeats every row once per year. The integer carries the same
  # grouping at a quarter of the width and is what stage 10 clusters on.
  u <- u[order(u$subclass, -u$treat, method = "radix"), , drop = FALSE]
  u$pair_id <- as.integer(factor(u$subclass, levels = unique(u$subclass)))
  u$unit_id <- seq_len(nrow(u))

  u$pairing <- pairing
  u$arm     <- sel$arm
  u$config  <- sel$config

  keep <- c("unit_id", "pair_id", "pixel_id", "lulc", "treat", "gvar",
            "cohort_year", "pa_est_year", "pa_year_known", "stratum",
            "analysis_group", "ecoregion_l3", "huc8", "weights", "distance",
            "subclass", "win_key", "level_year_m8", "eligibility",
            "pairing", "arm", "config")
  u <- u[, intersect(keep, names(u)), drop = FALSE]

  for (v in c("lulc", "stratum", "analysis_group", "ecoregion_l3", "huc8",
              "eligibility", "pairing", "arm", "config")) {
    if (v %in% names(u)) u[[v]] <- factor(u[[v]])
  }

  arrow::write_parquet(arrow::as_arrow_table(u), out_path,
                       compression = ca_io$compression)
  ca_log("  wrote ", basename(out_path), "  ",
         format(nrow(u), big.mark = ","), " rows  ",
         round(file.size(out_path) / 1e6, 1), " MB")

  # SUMMARY. Control reuse is the number ca_pending() asks for in the
  # cohort-stratified matching sentence, so it is computed once here rather
  # than recovered later from the panel.
  ct <- u[u$treat == 0L, , drop = FALSE]
  trr <- u[u$treat == 1L, , drop = FALSE]
  s <- data.frame(
    pairing            = pairing,
    arm                = sel$arm,
    config             = sel$config,
    n_units            = nrow(u),
    n_pairs            = length(unique(u$pair_id)),
    n_treated_pixels   = length(unique(trr$pixel_id)),
    n_control_rows     = nrow(ct),
    n_control_pixels   = length(unique(ct$pixel_id)),
    control_reuse      = round(nrow(ct) / max(length(unique(ct$pixel_id)), 1), 4),
    n_cohorts          = length(unique(trr$gvar)),
    gvar_min           = min(trr$gvar), gvar_max = max(trr$gvar),
    n_always_treated   = sum(trr$gvar <= min(ca_did$panel_years)),
    pct_always_treated = round(100 * mean(trr$gvar <= min(ca_did$panel_years)), 2),
    pa_year_unknown    = sum(!trr$pa_year_known),
    cohort_shift_m8    = shift,
    emfx_type          = ca_did_emfx_type(pairing),
    stringsAsFactors   = FALSE
  )
  print(s)

  ca_log("  treated pairs per stratum:")
  print(table(as.character(trr$stratum), useNA = "ifany"))
  ca_log("  treated pairs per cohort:")
  print(table(trr$gvar))

  s
}


# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------

ca_log("Stage 9a, pairings: ", paste(PAIRINGS, collapse = ", "))

summ <- list()
for (p in PAIRINGS) {
  s <- build_units(p)
  if (!is.null(s)) summ[[length(summ) + 1L]] <- s
  ca_gc(paste0("after ", p))
}

if (length(summ)) {
  s <- do.call(rbind, summ)
  path <- ca_meta("did_units_summary.csv")
  if (file.exists(path) && !OVERWRITE) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    s <- rbind(old[!old$pairing %in% s$pairing, , drop = FALSE], s)
    s <- s[order(match(s$pairing, ca_pairings$pairing)), , drop = FALSE]
  }
  utils::write.csv(s, path, row.names = FALSE)
  ca_log("wrote meta/did_units_summary.csv")
  print(s)
}

ca_stamp("9a_did_units")
ca_log("Stage 9a complete.")
