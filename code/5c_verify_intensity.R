### CA carbon revision pipeline
### Stage 5c. Thinning intensity dose-response
###
### Tests whether the Knight Low/Medium/High assignment is ordered in measured
### canopy loss, using the stage 5b screen as the validator. Knight et al.
### (2022) define intensity as basal area removed detectable from above and
### verify the categorisation against remotely sensed change in their section
### 3.4, so this is their own design run on this pixel set.
###
### TWO-PART STRUCTURE
###
### The screen records detections only, so an event with no matching screen row
### had loss at or below threshold. Those enter as zero rather than as missing,
### because conditioning on detection selects on the outcome.
###
### At the observed detection rate of roughly 7 percent the unconditional
### median, quartiles and p90 are all exactly zero and carry no information.
### Every quantity below is therefore reported twice. The extensive margin is
### the detection rate. The intensive margin is the loss distribution among
### detected events. A first run that reported only the mixture is what made
### the pooled result look like an inversion.
###
### STRATIFICATION IS NOT OPTIONAL
###
### Low events are almost entirely FACTS and Medium and High are dominated by
### CAL FIRE THP, and the two archives differ in detection rate by a factor of
### five. Pooling across source_layer produces a composition artefact rather
### than a dose-response. Section 2 is the result; section 1 is context.
###
### Usage
###   sbatch 5c_verify_intensity.sub
###   sbatch --array=2 5c_verify_intensity.sub          # Everg only
###   CA_SECTION=layer Rscript 5c_verify_intensity.R 2  # one section
###
### Sections: summary layer origin separation activity excluded lag all

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

task <- ca_task_id()
if (is.na(task)) task <- if (length(args)) as.integer(args[1]) else 2L
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id must be 1 to ", nrow(ca_lulc), ", got ", task)
}
CLASS <- ca_lulc$label[task]

SECT  <- Sys.getenv("CA_SECTION", "all")
STAMP <- format(Sys.time(), "%Y%m%d")
THR   <- ca_screen_threshold("primary")

LAGS      <- 0:1        # primary detection window
LAG_SENS  <- -6:6       # sensitivity sweep. The first run peaked at the -2 edge
                        # for thin_thp, so the window has to open further back
                        # before the peak can be located rather than truncated.
SUBSAMPLE <- 2e5        # per class, rank statistics only
MIN_N     <- 500L       # per-activity reporting floor

SECTIONS <- c("summary", "layer", "origin", "separation", "activity",
              "excluded", "lag")
if (!SECT %in% c(SECTIONS, "all")) stop("Unknown section: ", SECT)
run <- function(x) SECT == "all" || SECT == x

report <- function(df, name) {
  p <- ca_meta(sprintf("verify_stage5c_%s_%s_%s.csv", name, CLASS, STAMP))
  utils::write.csv(df, p, row.names = FALSE)
  ca_log("Wrote ", basename(p))
  print(utils::head(as.data.frame(df), 40))
  invisible(p)
}

# ---------------------------------------------------------------------------
# KEYS
# ---------------------------------------------------------------------------
# A numeric composite key rather than paste(). The character version built two
# vectors of seventeen million strings inside every lag iteration, which is what
# killed the first run. pixel_id tops out near 5e8 and the year term is four
# digits, so the product stays under 2^53 and match() stays exact.
#
# Adding L to the key advances the year by L, because year occupies the low
# order term with 1e4 spacing. No second key is built per lag.

ca_key <- function(pixel_id, year) as.double(pixel_id) * 1e4 + as.double(year)

attach_loss <- function(ev_key, lags, sc_key, sc_tf, sc_ab) {
  # Maximum inside the window, not the sum. Detections in consecutive years
  # describe one entry, and summing would reward the width of the window.
  tf <- numeric(length(ev_key))
  ab <- rep(NA_real_, length(ev_key))
  for (L in lags) {
    k <- match(ev_key + L, sc_key)
    hit <- which(!is.na(k))
    if (!length(hit)) { rm(k); next }
    kk <- k[hit]
    tf[hit] <- pmax(tf[hit], sc_tf[kk])
    ab[hit] <- pmax(ab[hit], sc_ab[kk], na.rm = TRUE)
    rm(k, hit, kk); ca_gc()
  }
  ab[!is.finite(ab)] <- NA_real_   # pmax(NA, NA, na.rm = TRUE) returns -Inf
  list(tf = tf, ab = ab)
}

# ---------------------------------------------------------------------------
# 1. EVENT TABLE
# ---------------------------------------------------------------------------
# One row per pixel-event. Two polygons recording the same entry year on one
# pixel is one event, and keeping both would weight that pixel twice in every
# distribution below.

want <- c("pixel_id", "event_year", "class_value", "record_class",
          "source_layer", "activity_key", "xwalk_origin")
have <- intersect(want, names(arrow::open_dataset(ca_event_path("thin", CLASS))))
if (!all(c("pixel_id", "event_year", "class_value", "record_class") %in% have)) {
  stop("Stage 5a thinning output is missing a required column. Present: ",
       paste(have, collapse = ", "))
}
for (m in setdiff(want, have)) {
  ca_log("NOTE column absent, the matching section will be skipped: ", m)
}

th <- as.data.frame(arrow::read_parquet(
  ca_event_path("thin", CLASS), col_select = dplyr::all_of(have)))

th <- th[th$record_class == "treatment_LMH" & !is.na(th$event_year), ]
th$record_class <- NULL
th <- th[!duplicated(th[, c("pixel_id", "event_year", "class_value")]), ]

th$intensity <- factor(th$class_value, levels = c("Low", "Medium", "High"))
th <- th[!is.na(th$intensity), ]
th$class_value <- NULL
ca_gc()

ca_log("Treatment events: ", format(nrow(th), big.mark = ","),
       " across ", format(length(unique(th$pixel_id)), big.mark = ","),
       " pixels, class ", CLASS)

sc <- as.data.frame(arrow::read_parquet(
  ca_screen_path(CLASS),
  col_select = dplyr::all_of(c("pixel_id", "year", "treefrac_loss",
                               "agb_loss"))))
sc_key <- ca_key(sc$pixel_id, sc$year)
sc_tf  <- sc$treefrac_loss
sc_ab  <- sc$agb_loss
rm(sc); ca_gc()

ev_key <- ca_key(th$pixel_id, th$event_year)
L <- attach_loss(ev_key, LAGS, sc_key, sc_tf, sc_ab)
th$treefrac_loss <- L$tf
th$agb_loss      <- L$ab
th$detected      <- L$tf > THR
rm(L); ca_gc()

ca_log("Detection window ", min(LAGS), " to ", max(LAGS),
       ", threshold ", THR, ", overall detection rate ",
       round(100 * mean(th$detected), 2), " percent")

# ---------------------------------------------------------------------------
# THE SUMMARY
# ---------------------------------------------------------------------------
# Extensive and intensive margins side by side. det_ columns are conditional on
# detection and are NA where no event in the cell was detected. agb_loss is only
# populated where TreeFrac fired, so det_agb_p50 is conditional by construction
# and is in gC per square metre, not a fraction.

two_part <- function(d, ...) {
  d %>%
    group_by(...) %>%
    summarise(
      n_events     = n(),
      n_detected   = sum(detected),
      detect_rate  = round(100 * mean(detected), 2),
      loss_mean    = round(mean(treefrac_loss), 5),
      loss_p95     = round(stats::quantile(treefrac_loss, 0.95), 4),
      loss_p99     = round(stats::quantile(treefrac_loss, 0.99), 4),
      det_loss_p25 = round(stats::quantile(treefrac_loss[detected], 0.25), 4),
      det_loss_p50 = round(stats::median(treefrac_loss[detected]), 4),
      det_loss_p75 = round(stats::quantile(treefrac_loss[detected], 0.75), 4),
      det_agb_p50  = round(stats::median(agb_loss[detected], na.rm = TRUE), 3),
      .groups = "drop")
}

if (run("summary")) report(two_part(th, intensity), "by_intensity")

# ---------------------------------------------------------------------------
# 2. BY SOURCE LAYER
# ---------------------------------------------------------------------------
# The decisive table. Knight states that group selection is Low in FACTS and
# Medium in CAL FIRE, so the two archives cannot be assumed comparable, and the
# observed detection rates differ far more than the intensity classes do.
# Ordering is judged inside a layer, never across.

if (run("layer") && "source_layer" %in% names(th)) {
  lay <- two_part(th, source_layer, intensity)
  report(lay, "by_layer_intensity")

  # Monotonicity inside each layer, stated rather than left to the reader.
  #
  # Every monotonicity flag is computed BEFORE the matching column is collapsed
  # to a string. The first version reused the name det_loss_p50 for the pasted
  # version and then called diff() on it. dplyr evaluates summarise expressions
  # in order, so diff() received a length-one character vector, returned
  # length zero, and all(logical(0)) is TRUE. Every loss_monotone value in that
  # run was vacuously TRUE and none of them were tested.
  mono <- lay %>%
    filter(n_events >= MIN_N) %>%
    group_by(source_layer) %>%
    arrange(intensity, .by_group = TRUE) %>%
    summarise(
      n_levels        = n(),
      detect_monotone = n() > 1 && all(diff(detect_rate) > 0),
      loss_monotone   = n() > 1 && all(diff(det_loss_p50) > 0),
      agb_monotone    = n() > 1 && all(diff(det_agb_p50) > 0),
      levels_present  = paste(intensity, collapse = " "),
      detect_rates    = paste(detect_rate, collapse = " "),
      loss_p50s       = paste(det_loss_p50, collapse = " "),
      agb_p50s        = paste(det_agb_p50, collapse = " "),
      .groups = "drop")
  report(mono, "layer_monotonicity")
}

if (run("origin") && "xwalk_origin" %in% names(th)) {
  report(two_part(th, xwalk_origin, intensity), "by_origin_intensity")
  # Share of each class that rests on the alias table rather than the published
  # crosswalk, which bounds how much of the ordering is judgement.
  sh <- th %>%
    group_by(intensity, xwalk_origin) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(pct_of_class = round(100 * n / sum(n), 2)) %>%
    ungroup()
  report(sh, "origin_share")
}

# ---------------------------------------------------------------------------
# 3. SEPARATION
# ---------------------------------------------------------------------------
# No p-values. At seventeen million events every contrast is significant and the
# number carries nothing. Cliff's delta is 2 * AUC - 1, so 0 is complete overlap
# and 1 is complete separation.
#
# The unconditional delta is reported with its tie fraction, because with 93
# percent of events at exactly zero the statistic is compressed toward zero for
# arithmetic reasons rather than substantive ones. The conditional delta is the
# one to read for magnitude, and the detection-rate difference is the one to
# read for the extensive margin.

cliffs <- function(a, b, n = SUBSAMPLE) {
  if (!length(a) || !length(b)) return(NA_real_)
  if (length(a) > n) a <- sample(a, n)
  if (length(b) > n) b <- sample(b, n)
  na <- as.double(length(a)); nb <- as.double(length(b))
  r <- rank(c(a, b))
  u <- sum(r[seq_len(length(a))]) - na * (na + 1) / 2
  round(2 * (u / (na * nb)) - 1, 3)
}

if (run("separation")) {
  strata <- if ("source_layer" %in% names(th)) {
    c(list(pooled = rep(TRUE, nrow(th))),
      lapply(split(seq_len(nrow(th)), th$source_layer),
             function(i) { z <- logical(nrow(th)); z[i] <- TRUE; z }))
  } else {
    list(pooled = rep(TRUE, nrow(th)))
  }

  sep <- do.call(rbind, lapply(names(strata), function(s) {
    z <- strata[[s]]
    pick <- function(lvl, det_only) {
      j <- z & th$intensity == lvl & (!det_only | th$detected)
      th$treefrac_loss[j]
    }
    do.call(rbind, lapply(list(c("Medium", "Low"), c("High", "Medium"),
                               c("High", "Low")), function(p) {
      hi_u <- pick(p[1], FALSE); lo_u <- pick(p[2], FALSE)
      hi_d <- pick(p[1], TRUE);  lo_d <- pick(p[2], TRUE)
      # n_det_hi and n_det_lo are the sample sizes delta_detected actually
      # rests on. Without them a delta of 0.45 computed on 631 detected events
      # reads the same as one computed on half a million.
      data.frame(
        stratum = s,
        contrast = paste(p[1], "vs", p[2]),
        n_hi = length(hi_u), n_lo = length(lo_u),
        n_det_hi = length(hi_d), n_det_lo = length(lo_d),
        detect_rate_diff_pp = round(
          100 * (mean(hi_u > THR) - mean(lo_u > THR)), 2),
        tie_frac = round(mean(c(hi_u, lo_u) == 0), 3),
        delta_uncond = cliffs(hi_u, lo_u),
        delta_detected = cliffs(hi_d, lo_d),
        stringsAsFactors = FALSE)
    }))
  }))
  report(sep, "separation")
  rm(strata); ca_gc()

  s <- th[sample(nrow(th), min(nrow(th), 3 * SUBSAMPLE)), ]
  ca_log("Spearman rho, intensity ordinal against canopy loss, pooled: ",
         round(stats::cor(as.integer(s$intensity), s$treefrac_loss,
                          method = "spearman"), 3))
  ca_log("  pooled rho is confounded by source_layer. Per-layer rho:")
  if ("source_layer" %in% names(s)) {
    for (lv in sort(unique(s$source_layer))) {
      j <- s$source_layer == lv
      if (sum(j) > 1000 && length(unique(s$intensity[j])) > 1) {
        ca_log("    ", lv, "  rho ",
               round(stats::cor(as.integer(s$intensity[j]),
                                s$treefrac_loss[j], method = "spearman"), 3),
               "  n ", format(sum(j), big.mark = ","))
      }
    }
  }
  rm(s); ca_gc()
}

# ---------------------------------------------------------------------------
# 4. PER-ACTIVITY
# ---------------------------------------------------------------------------
# The table that finds a single miscoded activity. An activity whose detected
# loss sits with a different class than the one it was assigned is the
# candidate, and the fix is the crosswalk row rather than the stratum.

if (run("activity") && "activity_key" %in% names(th)) {
  act <- two_part(th, source_layer, intensity, activity_key)
  act <- act[act$n_events >= MIN_N, ]
  cls <- two_part(th, intensity)
  act$class_det_p50 <- cls$det_loss_p50[match(act$intensity, cls$intensity)]
  act$class_detect_rate <- cls$detect_rate[match(act$intensity, cls$intensity)]
  act$ratio_to_class <- round(act$det_loss_p50 /
                                pmax(act$class_det_p50, 1e-6), 2)
  report(act[order(act$intensity, -act$n_events), ], "by_activity")

  off <- act[!is.na(act$ratio_to_class) &
               (act$ratio_to_class < 0.5 | act$ratio_to_class > 2), ]
  if (nrow(off)) {
    ca_log("Activities more than a factor of two from their class median:")
    print(as.data.frame(off[, c("source_layer", "activity_key", "intensity",
                                "n_events", "detect_rate", "det_loss_p50",
                                "class_det_p50", "ratio_to_class")]))
  }
  rm(act, cls, off); ca_gc()
}

# ---------------------------------------------------------------------------
# 5. EXCLUDED LABELS
# ---------------------------------------------------------------------------
# Where Variable and Unknown sit. If they overlay High, excluding them is
# conservative. If they overlay Medium, the exclusion is costing sample for
# nothing and that belongs in the Methods rather than left implicit.

if (run("excluded")) {
  tv <- as.data.frame(arrow::read_parquet(
    ca_event_path("thin", CLASS),
    col_select = dplyr::all_of(c("pixel_id", "event_year", "class_value",
                                 "record_class"))))
  tv <- tv[!is.na(tv$event_year) & !ca_thin_is_eligible(tv$class_value), ]
  tv <- tv[!duplicated(tv[, c("pixel_id", "event_year", "class_value")]), ]
  if (nrow(tv)) {
    Lx <- attach_loss(ca_key(tv$pixel_id, tv$event_year), LAGS,
                      sc_key, sc_tf, sc_ab)
    tv$treefrac_loss <- Lx$tf
    tv$agb_loss      <- Lx$ab
    tv$detected      <- Lx$tf > THR
    rm(Lx); ca_gc()
    report(two_part(tv, record_class, class_value) %>%
             arrange(desc(n_events)), "excluded_labels")
  }
  rm(tv); ca_gc()
}

# ---------------------------------------------------------------------------
# 6. LAG SENSITIVITY AND CUMULATIVE DETECTION
# ---------------------------------------------------------------------------
# FACTS completion dates and CAL FIRE plan dates need not align with the Almanac
# compositing year, and a THP is valid for several years, so an operation can
# post-date its plan year well outside a one-year window.
#
# Two questions. Does a single lag carry the detections, which would mean the
# window is misplaced. Or does the cumulative rate stay flat as the window
# widens, which would mean the detections are absent rather than late and the
# problem is polygon geometry rather than dates.

if (run("lag")) {
  by_lag <- do.call(rbind, lapply(LAG_SENS, function(l) {
    tf <- attach_loss(ev_key, l, sc_key, sc_tf, sc_ab)$tf
    d <- data.frame(source_layer = if ("source_layer" %in% names(th))
      th$source_layer else "all",
      intensity = th$intensity, det = tf > THR,
      stringsAsFactors = FALSE)
    out <- d %>%
      group_by(source_layer, intensity) %>%
      summarise(n_events = n(),
                detect_rate = round(100 * mean(det), 2), .groups = "drop")
    out$lag <- l
    rm(tf, d); ca_gc()
    out
  }))
  report(by_lag[, c("lag", "source_layer", "intensity", "n_events",
                    "detect_rate")], "lag_sensitivity")

  # PEAK LAG. The lag at which each series detects most. A peak away from zero
  # means the archive date is offset from the disturbance by that many years,
  # and the sign says which way. The first run put thin_thp at the -2 edge of
  # the window, which is why LAG_SENS now opens to -6.
  peak <- by_lag %>%
    group_by(source_layer, intensity) %>%
    arrange(desc(detect_rate), .by_group = TRUE) %>%
    summarise(n_events = first(n_events),
              peak_lag = first(lag),
              peak_rate = first(detect_rate),
              rate_at_zero = detect_rate[lag == 0],
              at_window_edge = first(lag) %in% range(LAG_SENS),
              .groups = "drop")
  report(peak, "peak_lag")
  if (any(peak$at_window_edge)) {
    ca_log("WARNING peak sits at the LAG_SENS edge, widen the window: ",
           paste(peak$source_layer[peak$at_window_edge],
                 peak$intensity[peak$at_window_edge], collapse = "; "))
  }

  # Cumulative over a SYMMETRIC expanding window centred on the event year. The
  # forward-only version missed every layer whose date runs late, which is the
  # failure mode this section exists to detect.
  cum <- NULL
  acc <- attach_loss(ev_key, 0, sc_key, sc_tf, sc_ab)$tf > THR
  base <- data.frame(source_layer = if ("source_layer" %in% names(th))
    th$source_layer else "all", intensity = th$intensity,
    stringsAsFactors = FALSE)
  for (w in 0:max(abs(range(LAG_SENS)))) {
    if (w > 0) {
      for (l in c(-w, w)) {
        if (l < min(LAG_SENS) || l > max(LAG_SENS)) next
        acc <- acc | (attach_loss(ev_key, l, sc_key, sc_tf, sc_ab)$tf > THR)
      }
    }
    o <- base %>%
      mutate(det = acc) %>%
      group_by(source_layer, intensity) %>%
      summarise(cum_detect_rate = round(100 * mean(det), 2), .groups = "drop")
    o$half_width <- w
    cum <- rbind(cum, o)
    rm(o); ca_gc()
  }
  report(cum[, c("half_width", "source_layer", "intensity",
                 "cum_detect_rate")], "cumulative_detection")
  rm(acc, base); ca_gc()
}

ca_stamp(paste0("5c_verify_intensity_", SECT))
ca_log("Done, class ", CLASS, ", section ", SECT)
