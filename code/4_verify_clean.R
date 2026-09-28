### CA carbon revision pipeline
### Stage 4 verification and the disturbance screen threshold diagnostic
###
### Sections
###
###   paths        every registry layer present in its expected store
###   units        cross-source agreement after ca_unit_conversion
###   sentinel     undeclared INT2S type limits passing through as data
###   rules        what ca_clean_rules changes, counted per layer and year
###   covariates   the stage 4b table, completeness and stratum integrity
###   degradation  late-record Carbon_GPP and Veg_TreeFrac failures, and
###                whether they sit in burn scars or scatter
###   threshold    derive the disturbance screen threshold from the data
###   all          paths, units, rules, covariates. Not threshold, expensive
###
### Usage
###   Rscript 4_verify_clean.R paths      [class]
###   Rscript 4_verify_clean.R units      [class]
###   Rscript 4_verify_clean.R sentinel   [class]
###   Rscript 4_verify_clean.R rules      [class]
###   Rscript 4_verify_clean.R covariates
###   Rscript 4_verify_clean.R degradation [class]
###   Rscript 4_verify_clean.R threshold  [class]
###   Rscript 4_verify_clean.R all
###
### class defaults to Everg, which is 95.6 of the 104.2 million pixels and is
### the partition every quoted figure comes from.
###
###
### WHY THE THRESHOLD IS DERIVED RATHER THAN ASSERTED
###
### The submitted code screens control pixels at thrd = 0 while Methods L225 to
### L226 state greater than 5 percent canopy loss. Neither survives scrutiny.
###
### At thrd = 0 a pixel is retained only by detecting exactly zero in all 39
### years, so retention of a truly undisturbed pixel is roughly (1 - p)^39 in
### the per-year false positive rate. The survivors are not a random subset.
### Change detection false positives concentrate on steep, heterogeneous, and
### high-biomass canopy, so the control pool ends up flatter, more homogeneous,
### and lower in biomass than the treated pixels it is the baseline for. That
### is a matching problem and an external validity problem at once.
###
### At 5 percent the error runs the other way. Five percent tree fraction is a
### real canopy opening in a closed conifer stand, and MTBS misses fires below
### roughly 400 ha in the West, so this screen is the only thing catching them.
### Admitting real disturbance to the business-as-usual baseline attenuates
### every treatment effect toward zero.
###
### The diagnostic below replaces both with a value read off this product's own
### error structure, using three independent lines of evidence.
###
###   1. Noise floor. Pixels with no MTBS record in any year and no stage 3c
###      management record of any kind. Their Disturbance_TreeFrac upper tail is
###      by construction mostly detection noise.
###   2. Separation. MTBS class 2, unburned to low, is the mildest real
###      disturbance in the data. The threshold maximising Youden's J against
###      the reference set is the value that best separates real low-severity
###      canopy loss from noise.
###   3. Year effects. Landsat 5 ends in 2011, Landsat 7 carries the SLC-off
###      gap from 2003, and Landsat 8 begins in 2013. If the reference noise
###      floor steps at those transitions, one absolute threshold is wrong for
###      the whole record and a year-specific percentile is used instead.
###
### Annual criterion only. A cumulative bound was considered and rejected.
### Repeated entry is already caught by the stage 3c record and by
### ca_thin_rules$multi_event_as_control = FALSE, and slow background decline
### from drought and beetle mortality is a common trend that the DiD
### differences out. See the note at ca_screen in ca_config.R.

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

if (!CLASS %in% ca_lulc$label) {
  stop("Unknown class '", CLASS, "'. One of: ",
       paste(ca_lulc$label, collapse = ", "))
}

STAMP <- format(Sys.time(), "%Y%m%d")
SEED  <- 20260801L

# Sample sizes. Tunable, and the only knobs in this script. The reference
# sample is the memory driver: n_ref by 39 years of doubles.
N_REF  <- 1e6L
N_FIRE <- 2e5L

report <- function(df, name) {
  path <- ca_meta(sprintf("verify_stage4_%s_%s_%s.csv", name, CLASS, STAMP))
  utils::write.csv(df, path, row.names = FALSE)
  ca_log("Wrote ", basename(path))
  print(utils::head(df, 40))
  invisible(path)
}

ca_covariate_path <- function(lulc) {
  file.path(ca_work("clean", "covariates"),
            sprintf("covariates_%s.parquet", lulc))
}

# ---------------------------------------------------------------------------
# SECTION 0. PATHS
# ---------------------------------------------------------------------------
# Every layer in every extraction registry, checked for the expected shard
# count in the store that ca_extract_store() resolves it to. This runs first
# and it runs fast, because an empty directory and a genuinely clean layer are
# indistinguishable in every other section, which is how a store name error
# passed as a clean result.

verify_paths <- function() {

  ca_log("Section 0. Paths, class ", CLASS)

  reg <- ca_layer_registry()

  out <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    lyr <- reg$layer[i]
    yrs <- reg$first_year[i]:reg$last_year[i]
    paths <- vapply(yrs, function(y)
      ca_extract_path(reg$store[i], lyr, y, CLASS), character(1))
    present <- file.exists(paths)
    act <- isTRUE(reg$active[i])
    data.frame(
      layer = lyr, source = reg$source[i], store = reg$store[i],
      active = act,
      first_year = reg$first_year[i], last_year = reg$last_year[i],
      n_expected = if (act) length(yrs) else 0L, n_present = sum(present),
      missing_years = if (!act || all(present)) "" else
        paste(yrs[!present], collapse = " "),
      stringsAsFactors = FALSE)
  }))

  # Static and vector stores have no year dimension and are checked separately.
  stat <- do.call(rbind, lapply(
    c(ca_static_layers$layer, ca_zone_layers$layer), function(lyr) {
      data.frame(layer = lyr, source = "static", store = "static",
                 active = TRUE,
                 first_year = NA_integer_, last_year = NA_integer_,
                 n_expected = 1L,
                 n_present = as.integer(file.exists(ca_static_path(lyr, CLASS))),
                 missing_years = "", stringsAsFactors = FALSE)
    }))

  vec <- do.call(rbind, lapply(ca_vector_layers$layer, function(lyr) {
    data.frame(layer = lyr, source = "vector", store = "vector",
               active = TRUE,
               first_year = NA_integer_, last_year = NA_integer_,
               n_expected = 1L,
               n_present = as.integer(file.exists(ca_vector_path(lyr, CLASS))),
               missing_years = "", stringsAsFactors = FALSE)
  }))

  res <- rbind(out, stat, vec)
  res$complete <- res$n_present >= res$n_expected
  report(res, "paths")

  inactive <- res[!res$active, ]
  if (nrow(inactive)) {
    ca_log("Inactive scaffold layers, absence expected: ",
           paste(inactive$layer, collapse = ", "))
  }

  bad <- res[res$active & !res$complete, ]
  if (nrow(bad)) {
    ca_log("INCOMPLETE stores, ", nrow(bad), " layers:")
    print(bad[, c("layer", "store", "n_expected", "n_present")])
    warning("Incomplete extract stores. Fix before running any other section.")
  } else {
    ca_log("All ", sum(res$active),
           " active layers complete in their expected stores.")
  }
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 1. UNITS
# ---------------------------------------------------------------------------
# The three biomass products converge on gC m-2 by different routes. If the
# conversion table is right they agree in level. If one factor is wrong the
# disagreement is a clean order of magnitude, which is why this is a level
# comparison rather than a correlation.

verify_units <- function() {

  ca_log("Section 1. Units, class ", CLASS)

  set.seed(SEED)
  grid <- ca_grid_points(CLASS)
  pid <- sort(sample(grid$pixel_id, min(2e5L, length(grid$pixel_id))))
  rm(grid); gc()

  # Years where all three biomass products exist.
  yrs_agb <- intersect(intersect(ca_outcome_years("agb_almanac"),
                                 ca_outcome_years("agb_emapr")),
                       ca_outcome_years("agb_lemma"))
  probe_agb <- yrs_agb[seq(1, length(yrs_agb), length.out = 5)]

  out <- list()

  for (y in round(probe_agb)) {
    for (oc in c("agb_almanac", "agb_emapr", "agb_lemma")) {
      ol <- ca_outcome_layer(oc)
      native <- ca_read_layer(ol$source, ol$layer, y, CLASS,
                              pixel_id = pid, clean = FALSE)
      raw <- ca_align(ca_read_shard(
        ca_extract_path(ca_extract_store(ol$layer), ol$layer, y, CLASS),
        pid), pid)$value
      out[[length(out) + 1]] <- data.frame(
        year = y, outcome = oc, layer = ol$layer,
        native_units = ca_unit_conversion$native[
          match(ol$layer, ca_unit_conversion$layer)],
        target_units = ca_target_units(ol$layer),
        n_valid = sum(!is.na(native$value)),
        raw_median = round(stats::median(raw, na.rm = TRUE), 3),
        conv_median = round(stats::median(native$value, na.rm = TRUE), 3),
        conv_p10 = round(stats::quantile(native$value, 0.10, na.rm = TRUE), 3),
        conv_p90 = round(stats::quantile(native$value, 0.90, na.rm = TRUE), 3),
        stringsAsFactors = FALSE)
      rm(native, raw); ca_gc()
    }
  }

  # GPP, two sources, same target units.
  yrs_gpp <- intersect(ca_outcome_years("gpp"),
                       ca_sources$first_year[ca_sources$source == "ncsda2022"]:
                         ca_sources$last_year[ca_sources$source == "ncsda2022"])
  for (y in round(yrs_gpp[seq(1, length(yrs_gpp), length.out = 3)])) {
    for (sl in list(c("almanac2026", "Carbon_GPP"), c("ncsda2022", "GPP"))) {
      v <- ca_read_layer(sl[1], sl[2], y, CLASS, pixel_id = pid,
                         clean = FALSE)$value
      out[[length(out) + 1]] <- data.frame(
        year = y, outcome = "gpp", layer = paste0(sl[1], ":", sl[2]),
        native_units = "gC/m2/yr", target_units = "gC/m2/yr",
        n_valid = sum(!is.na(v)),
        raw_median = round(stats::median(v, na.rm = TRUE), 3),
        conv_median = round(stats::median(v, na.rm = TRUE), 3),
        conv_p10 = round(stats::quantile(v, 0.10, na.rm = TRUE), 3),
        conv_p90 = round(stats::quantile(v, 0.90, na.rm = TRUE), 3),
        stringsAsFactors = FALSE)
      rm(v); gc()
    }
  }

  res <- do.call(rbind, out)
  report(res, "units")

  # One explicit gate. The three biomass medians should sit within a factor of
  # two of each other. eMapR runs about 20 percent below Carbon_AGB, which is a
  # product difference. A factor of 100 is a conversion error.
  agb <- res[res$outcome %in% c("agb_almanac", "agb_emapr", "agb_lemma"), ]
  by_oc <- tapply(agb$conv_median, agb$outcome, stats::median, na.rm = TRUE)
  ca_log("Biomass medians in gC/m2: ",
         paste(names(by_oc), round(by_oc, 1), sep = "=", collapse = "  "))
  if (max(by_oc, na.rm = TRUE) / min(by_oc, na.rm = TRUE) > 2.5) {
    warning("Biomass products differ by more than a factor of 2.5 after ",
            "conversion. Check ca_unit_conversion before proceeding.")
  }
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 2. CLEANING RULES
# ---------------------------------------------------------------------------
# ca_clean_rules counts in the config are from the Decid partition. This
# recomputes them on the requested class and on every year, so the numbers
# quoted in the Methods are the ones the pipeline actually applied.

verify_rules <- function() {

  ca_log("Section 2. Cleaning rules, class ", CLASS)

  out <- list()

  reg <- ca_layer_registry()

  for (lyr in unique(c(ca_clean_rules$layer, ca_sentinel_rules$layer))) {

    i <- match(lyr, reg$layer)
    if (is.na(i)) {
      warning("No registry entry for cleaning-rule layer ", lyr)
      next
    }
    years <- ca_layer_years(lyr)
    n_missing <- 0L

    for (y in years) {
      path <- ca_extract_path(ca_extract_store(lyr), lyr, y, CLASS)
      if (!file.exists(path)) { n_missing <- n_missing + 1L; next }
      r <- ca_read_layer(reg$source[i], lyr, y, CLASS, clean = TRUE)
      aff <- c(mask_sentinel = as.integer(r$n_sentinel), r$affected)
      aff <- aff[!is.na(aff)]
      if (!length(aff)) next
      out[[length(out) + 1]] <- data.frame(
        layer = lyr, year = y, rule = names(aff), n_affected = as.integer(aff),
        n_total = length(r$value),
        pct = round(100 * as.integer(aff) / length(r$value), 6),
        stringsAsFactors = FALSE)
      rm(r); ca_gc()
    }
    ca_gc(paste0("after ", lyr))
    # A layer that contributes no rows is either genuinely clean or absent from
    # disk, and those two cases must not look the same in the log. That is what
    # hid the almanac2026_screens store error.
    if (n_missing) {
      warning("Layer ", lyr, ": ", n_missing, " of ", length(years),
              " year shards missing from ", ca_extract_store(lyr))
    }
    ca_log("  ", lyr, " done, ", length(years) - n_missing, " of ",
           length(years), " year shards read")
  }

  res <- do.call(rbind, out)
  report(res, "rules")

  agg <- stats::aggregate(n_affected ~ layer + rule, data = res, FUN = sum)
  ca_log("Totals across all years:")
  print(agg)
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 3. COVARIATES
# ---------------------------------------------------------------------------

verify_covariates <- function() {

  ca_log("Section 3. Covariate table")

  out <- list()

  for (lulc in ca_lulc$label) {

    p <- ca_covariate_path(lulc)
    if (!file.exists(p)) {
      warning("Missing covariate table for ", lulc, ", run 4b first")
      next
    }
    cv <- as.data.frame(arrow::read_parquet(p))

    n_grid <- ca_lulc_n_points(lulc)

    for (v in setdiff(names(cv), c("pixel_id", "lulc"))) {
      x <- cv[[v]]
      is_chr <- is.character(x)
      out[[length(out) + 1]] <- data.frame(
        lulc = lulc, variable = v,
        n = nrow(cv),
        n_grid = n_grid,
        n_na = sum(is.na(x)),
        pct_na = round(100 * mean(is.na(x)), 4),
        n_distinct = if (is_chr) length(unique(x[!is.na(x)])) else NA_integer_,
        min = if (is_chr) NA_real_ else round(min(x, na.rm = TRUE), 4),
        median = if (is_chr) NA_real_ else
          round(stats::median(x, na.rm = TRUE), 4),
        max = if (is_chr) NA_real_ else round(max(x, na.rm = TRUE), 4),
        width = if (is_chr)
          paste(sort(unique(nchar(x[!is.na(x)]))), collapse = "/") else NA,
        stringsAsFactors = FALSE)
    }

    dup <- anyDuplicated(cv$pixel_id)
    ca_log("  ", lulc, "  rows ", format(nrow(cv), big.mark = ","),
           "  duplicate pixel_id ", dup)

    rm(cv); gc()
  }

  res <- do.call(rbind, out)
  report(res, "covariates")
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 1B. TYPE-LIMIT SENTINELS
# ---------------------------------------------------------------------------
# The degradation run showed Carbon_GPP reaching exactly -32768 in every year
# and Veg_TreeFrac reaching exactly 3.277 in every year. -32768 is the INT2S
# minimum. 3.277 is 32767 x 1e-4, the INT2S maximum after the registry scale.
# Both are type limits, not measurements.
#
# The registry declares nodata = -9999 and extraction removes it correctly, so
# these are a second, undeclared sentinel passing through as data. The current
# cleaning rules then convert them to real values: clip_lower turns a missing
# GPP into a GPP of zero, and clip_upper turns a saturated cover into a cover
# of 1.0. Both are fabricated numbers in the first case an outcome variable.
#
# What the rule should become depends on a split this measures rather than
# assumes. If every negative GPP is exactly the sentinel, the rule is negative
# to NA. If there is a population of small negatives between the sentinel and
# zero, those are genuine retrieval noise and the rule is sentinel to NA, small
# negative to zero. The same question applies at the top end of Veg_TreeFrac,
# and both limits are checked on every scaled layer rather than the two already
# known, because a sentinel that has not surfaced yet is still there.

# Both conventions are present. Disturbance_TreeFrac 2024 reaches -32767, not
# -32768, so a check looking only for the true type minimum misses it. All
# three candidates are tested and reported separately, since which one a layer
# uses is a fact about the product rather than something to assume.
ca_sentinel_codes <- c(-32768, -32767, 32767)

verify_sentinel <- function() {

  ca_log("Section 1b. Type-limit sentinels, class ", CLASS)

  reg <- ca_layer_registry()
  reg <- reg[reg$active %in% TRUE, , drop = FALSE]

  # Class layers carry no units and are not INT2S, so the limit check is
  # meaningless for them and ca_read_layer() correctly refuses to convert them.
  skip <- reg$layer[ca_is_categorical(reg$layer)]
  if (length(skip)) {
    ca_log("Skipping categorical layers, no units and not INT2S: ",
           paste(skip, collapse = ", "))
    reg <- reg[!ca_is_categorical(reg$layer), , drop = FALSE]
  }

  # Probe years spanning the record, not every year, since a sentinel present
  # at all is present throughout.
  out <- list()

  for (i in seq_len(nrow(reg))) {

    lyr <- reg$layer[i]
    yrs <- reg$first_year[i]:reg$last_year[i]
    probe <- unique(round(yrs[seq(1, length(yrs), length.out = 4)]))

    # The scale extraction already applied. Stored values carry it, so the
    # type limit is expressed in the same units. The unit conversion of stage 4
    # is NOT applied, because a limit lives in stored counts and comparing it
    # against a value already multiplied by 100 and 0.47 is a test that can
    # never fire. That invalidated Carbon_AGB, eMapR, and LEMMA on the first
    # run, all three of which reported zero because they could not report
    # anything else.
    sc <- c(ca_almanac_layers$scale[match(lyr, ca_almanac_layers$layer)],
            ca_ncsda_layers$scale[match(lyr, ca_ncsda_layers$layer)],
            ca_screen_layers$scale[match(lyr, ca_screen_layers$layer)],
            ca_mtbs_layers$scale[match(lyr, ca_mtbs_layers$layer)])
    sc <- sc[!is.na(sc)]
    sc <- if (length(sc)) sc[1] else 1

    codes <- ca_sentinel_codes * sc
    names(codes) <- c("m32768", "m32767", "p32767")

    for (y in probe) {

      v <- try(ca_read_layer(reg$source[i], lyr, y, CLASS,
                             clean = FALSE, convert = FALSE)$value,
               silent = TRUE)
      if (inherits(v, "try-error")) {
        warning("Read failed for ", lyr, " ", y, ": ", as.character(v))
        next
      }
      ok <- !is.na(v)

      hit <- vapply(codes, function(cd)
        sum(ok & abs(v - cd) <= abs(cd) * 1e-6 + 1e-9), numeric(1))

      n_neg <- sum(ok & v < 0)
      # Negatives that are not either low sentinel, which is the population
      # that decides whether the rule is invalid-to-NA or split.
      n_neg_real <- sum(ok & v < 0 &
                          abs(v - codes[["m32768"]]) > abs(codes[["m32768"]]) * 1e-6 &
                          abs(v - codes[["m32767"]]) > abs(codes[["m32767"]]) * 1e-6)

      out[[length(out) + 1]] <- data.frame(
        layer = lyr, year = y, scale_applied = sc,
        stored_units = ca_unit_conversion$native[
          match(lyr, ca_unit_conversion$layer)],
        n_valid = sum(ok),
        n_at_m32768 = as.integer(hit[["m32768"]]),
        n_at_m32767 = as.integer(hit[["m32767"]]),
        n_at_p32767 = as.integer(hit[["p32767"]]),
        n_negative = n_neg, n_negative_not_sentinel = n_neg_real,
        obs_min = round(min(v, na.rm = TRUE), 6),
        obs_max = round(max(v, na.rm = TRUE), 6),
        stringsAsFactors = FALSE)

      rm(v, ok); ca_gc()
    }
    ca_log("  ", lyr, " probed at ", paste(probe, collapse = ", "))
  }

  res <- do.call(rbind, out)
  res$n_sentinel <- res$n_at_m32768 + res$n_at_m32767 + res$n_at_p32767
  res$pct_sentinel <- round(100 * res$n_sentinel / res$n_valid, 8)
  report(res, "sentinel")

  hits <- res[res$n_sentinel > 0, ]
  if (nrow(hits)) {
    ca_log("Layers carrying an undeclared type-limit sentinel:")
    print(stats::aggregate(
      cbind(n_at_m32768, n_at_m32767, n_at_p32767, n_negative,
            n_negative_not_sentinel) ~ layer, data = hits, FUN = sum))
    ca_log("n_negative_not_sentinel decides the rule. Zero means every ",
           "negative is a sentinel and the rule is invalid to NA. Nonzero ",
           "means genuine values sit alongside it and the rule splits, ",
           "sentinel to NA and the remainder clipped.")
  } else {
    ca_log("No type-limit sentinel found on any active layer.")
  }
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 3B. LATE-RECORD DEGRADATION
# ---------------------------------------------------------------------------
# Two Almanac layers fail more often in recent years. Carbon_GPP returns
# negative values 19,343 times in 2021 and 17,331 in 2022 against a baseline in
# the hundreds, and Veg_TreeFrac exceeds 1.0 for 13,687 pixels in 2025 against
# a mid-record minimum near 3,000.
#
# Two explanations fit. Either the post-2020 fire seasons produce retrieval
# failures in burn scars, which is a property of the landscape and is confined
# to pixels that burned, or something in the Landsat stack degrades after 2020,
# which would affect every pixel and would make a robustness check run on the
# full record a different instrument at each end.
#
# The two are separable. If the bad pixels sit inside MTBS perimeters at a rate
# far above the base rate, it is the fires. If they scatter at roughly the base
# rate, it is the stack. This section measures that rather than arguing it.

verify_degradation <- function() {

  ca_log("Section 3b. Late-record degradation, class ", CLASS)

  yrs <- 2010:2025

  # Annual failure counts, raw values, so the clip is measured rather than
  # applied first and then counted.
  cnt <- do.call(rbind, lapply(yrs, function(y) {
    g <- ca_read_layer("almanac2026", "Carbon_GPP", y, CLASS, clean = FALSE)$value
    v <- ca_read_layer("almanac2026", "Veg_TreeFrac", y, CLASS,
                       clean = FALSE)$value
    r <- data.frame(
      year = y,
      n_valid_gpp = sum(!is.na(g)),
      gpp_negative = sum(!is.na(g) & g < 0),
      gpp_min = round(min(g, na.rm = TRUE), 3),
      n_valid_veg = sum(!is.na(v)),
      veg_above_1 = sum(!is.na(v) & v > 1),
      veg_max = round(max(v, na.rm = TRUE), 3),
      stringsAsFactors = FALSE)
    rm(g, v); gc()
    r
  }))
  cnt$gpp_neg_pct <- round(100 * cnt$gpp_negative / cnt$n_valid_gpp, 5)
  cnt$veg_hi_pct  <- round(100 * cnt$veg_above_1 / cnt$n_valid_veg, 5)
  report(cnt, "degradation_by_year")

  # Burn status against the failure, for the three worst years. A pixel counts
  # as burned if MTBS records any severity class in the year or the two before
  # it, which is the window a post-fire retrieval failure would sit in.
  grid <- ca_grid_points(CLASS)
  pid <- grid$pixel_id
  rm(grid); gc()

  cross <- list()

  for (y in c(2021, 2022, 2025)) {

    burned <- logical(length(pid))
    for (b in (y - 2):y) {
      p <- ca_extract_path(ca_extract_store("MTBS"), "MTBS", b, CLASS)
      if (!file.exists(p)) next
      mv <- ca_read_shard(p, NULL)
      j <- match(mv$pixel_id, pid)
      ok <- !is.na(j)
      hit <- !is.na(mv$value[ok]) & mv$value[ok] != ca_fire_nodata_code
      burned[j[ok][hit]] <- TRUE
      rm(mv); gc()
    }

    for (lyr in c("Carbon_GPP", "Veg_TreeFrac")) {
      x <- ca_read_layer("almanac2026", lyr, y, CLASS, clean = FALSE)$value
      bad <- if (lyr == "Carbon_GPP") !is.na(x) & x < 0 else !is.na(x) & x > 1
      base_rate <- mean(burned)
      in_burn <- if (sum(bad)) mean(burned[bad]) else NA_real_
      cross[[length(cross) + 1]] <- data.frame(
        year = y, layer = lyr,
        n_bad = sum(bad),
        pct_burned_base = round(100 * base_rate, 4),
        pct_burned_among_bad = round(100 * in_burn, 4),
        enrichment = round(in_burn / base_rate, 2),
        stringsAsFactors = FALSE)
      rm(x, bad); gc()
    }
    rm(burned); gc()
  }

  res <- do.call(rbind, cross)
  report(res, "degradation_burn")

  ca_log("Enrichment is the burn rate among failing pixels over the burn rate ",
         "in the class. Near 1 means the failures scatter and the cause is the ",
         "stack. Far above 1 means they sit in burn scars and the cause is the ",
         "fires.")
  invisible(res)
}

# ---------------------------------------------------------------------------
# SECTION 4. SCREEN THRESHOLD DIAGNOSTIC
# ---------------------------------------------------------------------------

# Pixels ever recorded as burned by MTBS. Class 6 is the non-processing mask
# and is not evidence of fire. Classes 1 and 5 are inside a perimeter, so they
# are evidence of fire even though they are not analysable severities, and a
# reference pixel must have none of them.
burned_ever <- function(pid) {

  ever <- logical(length(pid))
  fyear <- rep(NA_integer_, length(pid))   # earliest class 2 year
  yrs <- ca_mtbs_layers$first_year:ca_mtbs_layers$last_year

  for (y in yrs) {
    p <- ca_extract_path(ca_extract_store("MTBS"), "MTBS", y, CLASS)
    if (!file.exists(p)) next
    v <- ca_align(ca_read_shard(p, NULL), NULL)
    j <- match(v$pixel_id, pid)
    ok <- !is.na(j)
    val <- v$value[ok]; jj <- j[ok]

    inperim <- !is.na(val) & val != ca_fire_nodata_code
    ever[jj[inperim]] <- TRUE

    c2 <- !is.na(val) & val == 2
    hit <- jj[c2]
    upd <- is.na(fyear[hit])
    fyear[hit[upd]] <- y

    rm(v, val, jj); gc()
  }
  list(ever = ever, class2_year = fyear)
}

# Pixels with any stage 3c management record, of any activity class, dated or
# not. A reference pixel must have none.
recorded_ever <- function(pid) {
  ever <- logical(length(pid))
  layers <- ca_vector_layers$layer[grepl("^thin_", ca_vector_layers$layer)]
  for (lyr in layers) {
    p <- ca_vector_path(lyr, CLASS)
    if (!file.exists(p)) { warning("Missing 3c layer ", lyr); next }
    ids <- as.data.frame(
      dplyr::collect(dplyr::distinct(
        dplyr::select(arrow::open_dataset(p), pixel_id))))$pixel_id
    j <- match(ids, pid)
    ever[j[!is.na(j)]] <- TRUE
    ca_log("  ", lyr, " recorded pixels ", format(length(ids), big.mark = ","))
    rm(ids); gc()
  }
  ever
}

verify_threshold <- function() {

  ca_log("Section 4. Screen threshold diagnostic, class ", CLASS)

  grid <- ca_grid_points(CLASS)
  pid <- grid$pixel_id
  rm(grid); gc()

  ca_log("Building MTBS masks over ", format(length(pid), big.mark = ","),
         " pixels")
  mt <- burned_ever(pid)

  ca_log("Building stage 3c record mask")
  rec <- recorded_ever(pid)

  ref_ok <- !mt$ever & !rec
  ca_log("Reference set: ", format(sum(ref_ok), big.mark = ","),
         " pixels (", round(100 * mean(ref_ok), 2), " percent of class)")
  ca_log("MTBS class 2 pixels: ",
         format(sum(!is.na(mt$class2_year)), big.mark = ","))

  set.seed(SEED)
  ref_pid <- sort(sample(pid[ref_ok], min(N_REF, sum(ref_ok))))

  c2_idx <- which(!is.na(mt$class2_year))
  c2_take <- sort(sample(c2_idx, min(N_FIRE, length(c2_idx))))
  c2_pid <- pid[c2_take]
  c2_yr  <- mt$class2_year[c2_take]

  rm(mt, rec, ref_ok, pid); gc()

  # -- Reference panel ------------------------------------------------------
  yrs <- ca_screen$years
  ca_log("Reading Disturbance_TreeFrac for the reference sample, ",
         length(yrs), " years")

  ref_mat <- ca_read_years("almanac2026", ca_screen$layer, yrs, CLASS,
                           pixel_id = ref_pid, clean = TRUE)

  # -- 4a. Noise floor, pooled ----------------------------------------------
  pooled <- as.numeric(ref_mat)
  pooled <- pooled[!is.na(pooled)]

  qs <- c(0.50, 0.75, 0.90, 0.95, 0.975, 0.99, 0.999)
  noise <- data.frame(
    quantile = qs,
    value = round(as.numeric(stats::quantile(pooled, qs)), 6),
    stringsAsFactors = FALSE)
  ca_log("Reference pooled quantiles of Disturbance_TreeFrac:")
  print(noise)
  ca_log("Share of reference pixel-years at exactly zero: ",
         round(100 * mean(pooled == 0), 2), " percent")

  # -- 4b. Year effects -----------------------------------------------------
  by_year <- data.frame(
    year = yrs,
    n_valid = apply(ref_mat, 2, function(v) sum(!is.na(v))),
    zero_pct = round(100 * apply(ref_mat, 2,
                                 function(v) mean(v == 0, na.rm = TRUE)), 2),
    p95 = round(apply(ref_mat, 2, stats::quantile, 0.95, na.rm = TRUE), 6),
    p99 = round(apply(ref_mat, 2, stats::quantile, 0.99, na.rm = TRUE), 6),
    stringsAsFactors = FALSE)
  report(by_year, "threshold_by_year")

  # Undefined when the reference p95 is zero in every year, which is what a
  # sparse detection product gives. Reporting 0 there, as an earlier version
  # did, reads as perfect stability when it actually means the statistic does
  # not apply. Fall back to the share of nonzero pixel-years per year, which is
  # the quantity the ratio was standing in for.
  if (all(by_year$p95 == 0, na.rm = TRUE)) {
    nz <- 100 - by_year$zero_pct
    step <- NA_real_
    ca_log("Reference p95 is zero in every year, so the ratio is undefined. ",
           "Nonzero detection rate per year runs ",
           round(min(nz), 3), " to ", round(max(nz), 3),
           " percent, ratio ", round(max(nz) / max(min(nz), 1e-9), 2))
    ca_log("Disturbance_TreeFrac is a sparse detection product, not a noisy ",
           "continuous field. The compounding false-positive argument for ",
           "raising the threshold above zero does not apply to it.")
  } else {
    step <- max(by_year$p95, na.rm = TRUE) /
      max(min(by_year$p95, na.rm = TRUE), 1e-9)
    ca_log("Reference p95 ratio across years, max over min: ", round(step, 2))
  }
  if (!is.na(step) && step > 3) {
    ca_log("NOTE: the noise floor is not stable across the record. A single ",
           "absolute threshold is not defensible. Use the year-specific p99 ",
           "column of the by_year table as the screen instead.")
  }

  # -- 4c. Separation against MTBS class 2 ----------------------------------
  ca_log("Reading Disturbance_TreeFrac in the fire window for class 2 pixels")

  # The window comes from ca_screen$fire_detection_window rather than the
  # thinning one. Fire detection can lag differently from a harvest record, and
  # a sensitivity of 0.33 at plus or minus 1 year is low enough that the window
  # itself is a candidate explanation rather than a settled parameter.
  fw <- ca_screen$fire_detection_window
  c2_val <- rep(NA_real_, length(c2_pid))
  for (y in sort(unique(c2_yr))) {
    take <- which(c2_yr == y)
    win <- intersect((y - fw):(y + fw), yrs)
    if (!length(win)) next
    m <- ca_read_years("almanac2026", ca_screen$layer, win, CLASS,
                       pixel_id = c2_pid[take], clean = TRUE)
    c2_val[take] <- suppressWarnings(apply(m, 1, max, na.rm = TRUE))
    rm(m)
  }
  c2_val[!is.finite(c2_val)] <- NA_real_

  # One observation per reference pixel, year drawn at random, so the two
  # samples have the same design and the comparison is not inflated by the
  # reference set contributing 39 observations per pixel.
  set.seed(SEED + 1L)
  pick <- sample.int(ncol(ref_mat), nrow(ref_mat), replace = TRUE)
  ref_one <- ref_mat[cbind(seq_len(nrow(ref_mat)), pick)]
  ref_one <- ref_one[!is.na(ref_one)]
  c2_one <- c2_val[!is.na(c2_val)]

  ca_log("Separation samples: reference ", format(length(ref_one),
                                                  big.mark = ","),
         "  class 2 ", format(length(c2_one), big.mark = ","))

  cand <- seq(0, 0.10, by = 0.0025)
  roc <- data.frame(
    threshold = cand,
    sensitivity = vapply(cand, function(t) mean(c2_one > t), numeric(1)),
    specificity = vapply(cand, function(t) mean(ref_one <= t), numeric(1)),
    stringsAsFactors = FALSE)
  roc$youden <- roc$sensitivity + roc$specificity - 1
  report(roc, "threshold_roc")

  youden_t <- roc$threshold[which.max(roc$youden)]
  p99 <- noise$value[noise$quantile == 0.99]
  p95 <- noise$value[noise$quantile == 0.95]

  # The more conservative of the two lines of evidence becomes primary. A
  # smaller threshold keeps the control pool cleaner at the cost of pool size,
  # and pool size is the recoverable loss of the two.
  primary <- min(youden_t, p99)

  ca_log("Youden-optimal threshold: ", youden_t)
  ca_log("Reference p99: ", p99)
  ca_log("Primary threshold selected: ", primary)

  # -- 4d. Control pool retention under each candidate ----------------------
  annual_max <- suppressWarnings(apply(ref_mat, 1, max, na.rm = TRUE))
  annual_max[!is.finite(annual_max)] <- NA_real_

  grid_t <- sort(unique(c(ca_screen$threshold_grid, primary)))
  retain <- data.frame(
    threshold = grid_t,
    retention_pct = round(100 * vapply(
      grid_t, function(t) mean(annual_max <= t, na.rm = TRUE), numeric(1)), 3),
    stringsAsFactors = FALSE)
  retain$is_primary <- retain$threshold == primary
  report(retain, "threshold_retention")

  ca_log("Control pool retention among reference pixels:")
  print(retain)

  # -- Artifact -------------------------------------------------------------
  # The number alone is not the finding. Retention and sensitivity at the
  # chosen threshold are what make it defensible, so they travel with it.
  ret_primary <- retain$retention_pct[retain$threshold == primary][1]
  sens_primary <- roc$sensitivity[roc$threshold == primary][1]

  thr <- data.frame(
    name = c("primary", "youden", "noise_p95", "noise_p99"),
    value = c(primary, youden_t, p95, p99),
    basis = c("more conservative of youden and noise_p99",
              "max Youden J against MTBS class 2",
              "reference pooled 95th percentile",
              "reference pooled 99th percentile"),
    class = CLASS,
    n_reference = length(ref_pid),
    n_class2 = length(c2_one),
    year_stability_ratio = if (is.na(step)) NA_real_ else round(step, 3),
    retention_pct_at_primary = ret_primary,
    sensitivity_at_primary = round(sens_primary, 4),
    zero_share_reference = round(100 * mean(pooled == 0), 3),
    derived = as.character(Sys.Date()),
    stringsAsFactors = FALSE)

  utils::write.csv(thr, ca_screen$threshold_file, row.names = FALSE)
  ca_log("Wrote threshold artifact: ", ca_screen$threshold_file)
  print(thr)

  invisible(thr)
}

# ---------------------------------------------------------------------------

ca_log("Stage 4 verification, section '", section, "', class ", CLASS)

if (section %in% c("paths", "all")) { verify_paths(); ca_gc("after paths") }
if (section %in% c("units", "all")) { verify_units(); ca_gc("after units") }
if (section == "sentinel") { verify_sentinel(); ca_gc("after sentinel") }
if (section %in% c("rules", "all")) { verify_rules(); ca_gc("after rules") }
if (section %in% c("covariates", "all")) {
  verify_covariates(); ca_gc("after covariates")
}
if (section == "degradation") verify_degradation()
if (section == "threshold") verify_threshold()

if (!section %in% c("paths", "units", "sentinel", "rules", "covariates",
                    "degradation", "threshold", "all")) {
  stop("Unknown section '", section, "'")
}

ca_stamp("4_verify_clean")
ca_log("Done")
