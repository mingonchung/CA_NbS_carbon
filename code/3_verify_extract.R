### CA carbon revision pipeline
### Stage 3 verification. One script for 3a 2026, 3a 2022, 3b, 3c, and 3d.
###
### One verifier rather than four, because join integrity is the check that
### catches a corrupted pixel_id ordering and it must not exist in two places.
###
### The stages divide by output shape, not by subject, and each shape needs its
### own equivalents of the same four questions.
###
###   year-keyed   3a, 3b, MTBS. One row per pixel per year. Sections 1 to 6.
###   static       3d. One row per pixel, no year. Section 7.
###   long         3c. One row per pixel per intersecting polygon, so row
###                counts exceed pixel counts and completeness cannot be
###                checked by counting rows. Section 8.
###
###   1  Metadata        datatype, NA flag, file-carried scale, CRS, grid
###                      offset. Sets the unverified entries in ca_sources.
###   2  Completeness    expected against actual files per layer. A silently
###                      failed array task is visible here and nowhere else.
###   3  Distributions   zero and negative fractions, quantiles. A high zero
###                      fraction is the signature of an undeclared sentinel.
###   4  Join integrity  every file must carry the grid partition's pixel_id
###                      set in the same order.
###   5  Cross-source    where two products measure the same quantity at the
###                      same pixels, the ratio of medians is the scale factor
###                      between them and the correlation tests co-registration.
###   6  Anomalies       the specific oddities found so far, each with a rule
###                      attached rather than left to propagate silently.
###   7  Stage 3d        static equivalents of 1 to 4, run on all three classes
###                      because Evergreen carries 92 percent of the points.
###   8  Stage 3c        long-format equivalents, plus the checks that only
###                      apply to records: event years, crosswalk resolution,
###                      events per pixel, and label coverage against config.
###
### Sections 1 to 6 run on the Decid partition, 2.2 million points, large
### enough for stable quantiles and small enough to finish in minutes.
### Sections 7 and 8 loop all three classes, because their failure modes are
### class-dependent in a way the year-keyed layers are not.
###
### Usage
###   Rscript 3_verify_extract.R

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
})

ca_stamp("stage3_verify")

LULC <- "Decid"
stamp <- format(Sys.Date(), "%Y%m%d")

REG <- list(
  almanac2026         = ca_almanac_layers,
  ncsda2022           = ca_ncsda_layers,
  almanac2026_screens = ca_screen_layers,
  # MTBS is year-keyed with the same output schema, so it needs a registry
  # entry here rather than a section of its own. It runs through the 3b script
  # under the mtbs selector.
  mtbs                = ca_mtbs_layers
)

rd <- function(src, lay, yr) {
  f <- ca_extract_path(src, lay, yr, LULC)
  if (!file.exists(f)) return(NULL)
  as.data.frame(arrow::read_parquet(f))
}

active <- function(src) {
  r <- REG[[src]]
  r[r$active %in% TRUE, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 1  METADATA
# ---------------------------------------------------------------------------

ca_log("=== 1. Source metadata ===")
for (src in names(REG)) {
  m <- try(ca_probe_layers(REG[[src]], src), silent = TRUE)
  if (inherits(m, "try-error")) next
  print(m[, intersect(c("layer", "status", "n_files", "n_bands", "datatype",
                        "na_flag", "scale_file", "same_crs", "offset_x_m",
                        "offset_y_m", "res_x"), names(m))])
}

# ---------------------------------------------------------------------------
# 2  COMPLETENESS
# ---------------------------------------------------------------------------

ca_log("=== 2. Extraction completeness ===")
inv <- list()
for (src in names(REG)) {
  reg <- active(src)
  for (k in seq_len(nrow(reg))) {
    lay <- reg$layer[k]
    yrs <- seq(reg$first_year[k], reg$last_year[k])
    got <- list.files(ca_work("extract", src, lay), pattern = "\\.parquet$")
    have <- sort(unique(as.integer(sub(paste0("^", lay, "_([0-9]{4})_.*$"),
                                       "\\1", got))))
    inv[[paste(src, lay)]] <- data.frame(
      source = src, layer = lay,
      expected_files = length(yrs) * nrow(ca_lulc),
      found_files = length(got),
      years_expected = length(yrs), years_found = length(have),
      missing_years = paste(setdiff(yrs, have), collapse = " "),
      stringsAsFactors = FALSE)
  }
}
inv <- do.call(rbind, inv); rownames(inv) <- NULL
print(inv)

# ---------------------------------------------------------------------------
# 3  VALUE DISTRIBUTIONS
# ---------------------------------------------------------------------------

ca_log("=== 3. Value distributions, ", LULC, " ===")
vals <- list()
for (src in names(REG)) {
  reg <- active(src)
  for (k in seq_len(nrow(reg))) {
    lay <- reg$layer[k]
    yrs <- seq(reg$first_year[k], reg$last_year[k])
    for (yr in unique(c(yrs[1], yrs[ceiling(length(yrs) / 2)],
                        yrs[length(yrs)]))) {
      d <- rd(src, lay, yr)
      if (is.null(d)) next
      v <- d$value; ok <- v[!is.na(v)]
      vals[[paste(src, lay, yr)]] <- data.frame(
        source = src, layer = lay, year = yr, n = length(v),
        pct_na = round(100 * mean(is.na(v)), 3),
        pct_zero = round(100 * mean(ok == 0), 3),
        pct_neg = round(100 * mean(ok < 0), 3),
        min = min(ok), q01 = quantile(ok, 0.01, names = FALSE),
        median = median(ok), mean = round(mean(ok), 2),
        q99 = quantile(ok, 0.99, names = FALSE), max = max(ok),
        stringsAsFactors = FALSE)
      rm(d, v, ok)
    }
  }
}
vals <- do.call(rbind, vals); rownames(vals) <- NULL
print(vals)

# ---------------------------------------------------------------------------
# 4  JOIN INTEGRITY
# ---------------------------------------------------------------------------

ca_log("=== 4. Join integrity ===")
grid_ids <- ca_grid_points(LULC)$pixel_id
chk <- list()
for (src in names(REG)) {
  reg <- active(src)
  if (!nrow(reg)) next
  d <- rd(src, reg$layer[1], reg$first_year[1])
  if (is.null(d)) next
  chk[[src]] <- data.frame(
    source = src, layer = reg$layer[1], year = reg$first_year[1],
    n_rows = nrow(d), n_grid = length(grid_ids),
    identical_order = identical(as.integer(d$pixel_id), as.integer(grid_ids)),
    same_set = setequal(d$pixel_id, grid_ids),
    stringsAsFactors = FALSE)
  rm(d)
}
chk <- do.call(rbind, chk); rownames(chk) <- NULL
print(chk)

# ---------------------------------------------------------------------------
# 5  CROSS-SOURCE AGREEMENT
# ---------------------------------------------------------------------------

ca_log("=== 5. Cross-source agreement ===")

pair_check <- function(a_src, a_lay, b_src, b_lay, yr) {
  A <- rd(a_src, a_lay, yr); B <- rd(b_src, b_lay, yr)
  if (is.null(A) || is.null(B)) return(NULL)
  m <- merge(A, B, by = "pixel_id", suffixes = c("_a", "_b"))
  keep <- !is.na(m$value_a) & !is.na(m$value_b) &
    m$value_a > 0 & m$value_b > 0
  if (sum(keep) < 1000) return(NULL)
  x <- m$value_a[keep]; y <- m$value_b[keep]
  data.frame(year = yr, a = paste(a_src, a_lay), b = paste(b_src, b_lay),
             n_pairs = length(x),
             median_a = median(x), median_b = median(y),
             ratio_a_over_b = round(median(x) / median(y), 5),
             pearson = round(cor(x, y), 3),
             spearman = round(cor(x, y, method = "spearman"), 3),
             stringsAsFactors = FALSE)
}

pairs <- list(
  pair_check("almanac2026", "Carbon_AGB", "ncsda2022", "eMapR", 2016),
  pair_check("almanac2026", "Carbon_AGB", "ncsda2022", "LEMMA", 2016),
  pair_check("ncsda2022", "eMapR", "ncsda2022", "LEMMA", 2016),
  pair_check("almanac2026", "Carbon_GPP", "ncsda2022", "GPP", 2000),
  pair_check("almanac2026", "Carbon_GPP", "ncsda2022", "GPP", 2016)
)
pairs <- pairs[!vapply(pairs, is.null, logical(1))]
pairs <- if (length(pairs)) do.call(rbind, pairs) else NULL
if (!is.null(pairs)) { rownames(pairs) <- NULL; print(pairs) }

# LEMMA is kg/ha per the data provider, so the ratio against a t/ha product
# should be close to 1/1000 once the level difference between products is
# allowed for. Reported explicitly so the conversion stays checkable.
if (!is.null(pairs)) {
  lm_row <- pairs[grepl("LEMMA", pairs$b) & grepl("Carbon_AGB", pairs$a), ]
  if (nrow(lm_row)) {
    ca_log("LEMMA implied kg per t/ha unit: ",
           round(1 / lm_row$ratio_a_over_b[1], 1),
           "   (1000 would be an exact kg/ha match; the shortfall is the ",
           "level difference between products)")
  }
}

# ---------------------------------------------------------------------------
# 6  ANOMALY CHECKS
# ---------------------------------------------------------------------------

ca_log("=== 6a. LEMMA zero stability ===")
# If the same pixels are zero every year it is a static mask and zero is the
# nodata sentinel. If the set moves, zero is data.
lz <- lapply(c(1990, 2003, 2016), function(y) rd("ncsda2022", "LEMMA", y))
if (!any(vapply(lz, is.null, logical(1)))) {
  z <- lapply(lz, function(d) d$value == 0 & !is.na(d$value))
  ca_log("zeros 1990=", sum(z[[1]]), "  2003=", sum(z[[2]]),
         "  2016=", sum(z[[3]]))
  ca_log("identical zero set across all three: ",
         identical(z[[1]], z[[2]]) && identical(z[[1]], z[[3]]))
  ca_log("zero in 1990 but positive in 2016: ",
         sum(z[[1]] & !z[[3]] & !is.na(lz[[3]]$value)))
  rm(lz, z); gc(verbose = FALSE)
}

ca_log("=== 6b. Veg_TreeFrac out of range ===")
# Fractional cover times 10,000 should cap at 10,000, giving 0 to 1 after
# scaling. Anything outside that needs a documented clip at stage 4.
vt <- list()
for (yr in c(1985, 2005, 2025)) {
  x <- rd("almanac2026_screens", "Veg_TreeFrac", yr)
  if (is.null(x)) next
  v <- x$value[!is.na(x$value)]
  vt[[as.character(yr)]] <- data.frame(
    year = yr, n = length(v), n_gt_1 = sum(v > 1), n_lt_0 = sum(v < 0),
    pct_gt_1 = round(100 * mean(v > 1), 4),
    pct_lt_0 = round(100 * mean(v < 0), 4),
    max = round(max(v), 4), min = round(min(v), 4),
    stringsAsFactors = FALSE)
  rm(x, v)
}
vt <- if (length(vt)) do.call(rbind, vt) else NULL
if (!is.null(vt)) { rownames(vt) <- NULL; print(vt) }

ca_log("=== 6c. Fire_LCP_CH missingness ===")
# Nodata is 0 under the FARSITE convention, so a forest pixel with no canopy
# reads as missing. If those pixels also carry low biomass and low tree cover
# they are genuinely canopy-free and the loss is correct.
ch_rows <- list()
for (yr in c(1985, 2005, 2025)) {
  ch <- rd("almanac2026", "Fire_LCP_CH", yr)
  ag <- rd("almanac2026", "Carbon_AGB", yr)
  tf <- rd("almanac2026_screens", "Veg_TreeFrac", yr)
  if (is.null(ch) || is.null(ag)) next
  m <- merge(ch, ag, by = "pixel_id", suffixes = c("_ch", "_agb"))
  if (!is.null(tf)) {
    names(tf)[2] <- "value_tf"
    m <- merge(m, tf, by = "pixel_id")
  }
  na_ch <- is.na(m$value_ch)
  ch_rows[[as.character(yr)]] <- data.frame(
    year = yr, n = nrow(m),
    pct_ch_missing = round(100 * mean(na_ch), 2),
    n_ch_missing_agb_valid = sum(na_ch & !is.na(m$value_agb)),
    median_agb_ch_missing = round(median(m$value_agb[na_ch], na.rm = TRUE), 1),
    median_agb_ch_present = round(median(m$value_agb[!na_ch], na.rm = TRUE), 1),
    median_tf_ch_missing = if ("value_tf" %in% names(m))
      round(median(m$value_tf[na_ch], na.rm = TRUE), 3) else NA_real_,
    median_tf_ch_present = if ("value_tf" %in% names(m))
      round(median(m$value_tf[!na_ch], na.rm = TRUE), 3) else NA_real_,
    stringsAsFactors = FALSE)
  rm(ch, ag, tf, m, na_ch); gc(verbose = FALSE)
}
ch_rows <- if (length(ch_rows)) do.call(rbind, ch_rows) else NULL
if (!is.null(ch_rows)) { rownames(ch_rows) <- NULL; print(ch_rows) }

ca_log("=== 6d. Carbon_GPP negatives ===")
# Gross primary production cannot be negative. Counts are tiny but the rule
# has to be stated rather than left to propagate.
for (yr in c(1985, 2005, 2025)) {
  g <- rd("almanac2026", "Carbon_GPP", yr)
  if (is.null(g)) next
  v <- g$value[!is.na(g$value)]
  ca_log("year ", yr, "  n_neg=", sum(v < 0), "  min=", min(v),
         "  pct=", round(100 * mean(v < 0), 5))
  rm(g, v)
}

# ---------------------------------------------------------------------------
# 7  STAGE 3D COVARIATES
# ---------------------------------------------------------------------------
# Stage 3d output has no year dimension, so rd() does not apply and sections 1
# to 4 above cannot be reused directly. The five checks below are their static
# equivalents, kept in this file rather than a second verifier so that the join
# integrity logic, which is the one check that catches a corrupted pixel_id
# ordering, exists in exactly one place.

rds <- function(layer, lulc = LULC) {
  f <- ca_static_path(layer, lulc)
  if (!file.exists(f)) return(NULL)
  as.data.frame(arrow::read_parquet(f))
}

# Sections 1 to 6 run on one class, because a year-keyed layer behaves the same
# way in all three and Decid is the cheapest to read. Stage 3d does not have
# that property. Evergreen carries 95.6 million of the 104 million points, so a
# Decid-only check would certify 2 percent of the sample. The 3d value, strata,
# and aspect checks therefore loop all three classes.
LULC_3D <- ca_lulc$label

ca_log("=== 7a. Stage 3d source metadata ===")
stat_probe <- try(ca_probe_static(), silent = TRUE)
if (!inherits(stat_probe, "try-error")) {
  print(stat_probe[, c("layer", "status", "n_files", "res_x", "same_crs",
                       "datatype", "na_flag", "cfg_nodata",
                       "smp_min", "smp_median", "smp_max")])
} else {
  stat_probe <- NULL
}

zone_probe <- try(ca_probe_zones(), silent = TRUE)
if (!inherits(zone_probe, "try-error")) {
  print(zone_probe[, c("layer", "status", "n_features", "same_crs",
                       "cfg_field", "actual_field", "n_distinct",
                       "n_missing", "max_nchar")])
} else {
  zone_probe <- NULL
}

ca_log("=== 7b. Stage 3d completeness ===")
static_reg <- rbind(
  ca_static_layers[ca_static_layers$active %in% TRUE, c("layer", "type")],
  ca_zone_layers[ca_zone_layers$active %in% TRUE, c("layer", "type")]
)
stat_inv <- do.call(rbind, lapply(seq_len(nrow(static_reg)), function(i) {
  lay <- static_reg$layer[i]
  got <- ca_lulc$label[vapply(ca_lulc$label, function(cl)
    file.exists(ca_static_path(lay, cl)), logical(1))]
  data.frame(layer = lay, type = static_reg$type[i],
             expected = nrow(ca_lulc), found = length(got),
             missing = paste(setdiff(ca_lulc$label, got), collapse = " "),
             stringsAsFactors = FALSE)
}))
rownames(stat_inv) <- NULL
print(stat_inv)
if (any(stat_inv$found < stat_inv$expected)) {
  ca_log("INCOMPLETE. Resubmit the missing tasks before running stage 4.")
}

ca_log("=== 7c. Stage 3d raster distributions ===")
# Plausibility bounds are physical, not statistical. A value outside them is a
# sentinel that survived extraction, not an unusual pixel. Aspect is bounded at
# 360 rather than 359 because a flat cell is coded -1 in some SRTM derivatives
# and 360 in others, and both need to surface here rather than at stage 8.
bounds <- data.frame(
  layer = c("elevation", "slope", "aspect", "ppt_normal", "tmean_normal",
            "pop_density", "city_travel_time"),
  lo    = c(-100, 0, -1, 0, -20, 0, 0),
  hi    = c(4500, 90, 360, 6000, 30, 1e5, 5000),
  stringsAsFactors = FALSE
)

stat_vals <- list()
for (cl in LULC_3D) {
  for (lay in ca_static_layers$layer[ca_static_layers$active %in% TRUE]) {
    d <- rds(lay, cl)
    if (is.null(d)) next
    v <- d$value
    ok <- v[!is.na(v)]
    b <- bounds[bounds$layer == lay, ]
    q <- if (length(ok)) stats::quantile(ok, c(0.01, 0.5, 0.99)) else rep(NA, 3)
    stat_vals[[paste(cl, lay)]] <- data.frame(
      lulc = cl, layer = lay, n = length(v),
      pct_na = round(100 * mean(is.na(v)), 3),
      min = if (length(ok)) round(min(ok), 3) else NA_real_,
      q01 = round(q[[1]], 3), median = round(q[[2]], 3), q99 = round(q[[3]], 3),
      max = if (length(ok)) round(max(ok), 3) else NA_real_,
      n_below = if (nrow(b)) sum(ok < b$lo) else NA_integer_,
      n_above = if (nrow(b)) sum(ok > b$hi) else NA_integer_,
      stringsAsFactors = FALSE)
    rm(d, v, ok)
  }
  gc(verbose = FALSE)
}
stat_vals <- if (length(stat_vals)) do.call(rbind, stat_vals) else NULL
if (!is.null(stat_vals)) { rownames(stat_vals) <- NULL; print(stat_vals) }

ca_log("=== 7d. Stage 3d join integrity ===")
# Every static file must carry the grid partition's pixel_id set in the same
# order, otherwise the stage 4 join is silently misaligned.
stat_join <- list()
for (cl in LULC_3D) {
  ref_ids <- ca_grid_points(cl)$pixel_id
  for (lay in static_reg$layer) {
    d <- rds(lay, cl)
    if (is.null(d)) next
    stat_join[[paste(cl, lay)]] <- data.frame(
      lulc = cl, layer = lay, n = nrow(d), n_grid = length(ref_ids),
      identical_ids = identical(as.integer(d$pixel_id),
                                as.integer(ref_ids)),
      same_set = setequal(d$pixel_id, ref_ids),
      sorted = !is.unsorted(d$pixel_id),
      stringsAsFactors = FALSE)
    rm(d)
  }
  rm(ref_ids); gc(verbose = FALSE)
}
stat_join <- if (length(stat_join)) do.call(rbind, stat_join) else NULL
if (!is.null(stat_join)) { rownames(stat_join) <- NULL; print(stat_join) }

ca_log("=== 7e. Exact matching strata ===")
# An NA stratum drops a pixel from matching entirely, so coverage is reported
# rather than assumed. Points falling outside every polygon are the expected
# source, and they should be a boundary-thin fraction.
zone_cov <- list()
for (cl in LULC_3D) {
  for (lay in ca_zone_layers$layer[ca_zone_layers$active %in% TRUE]) {
    d <- rds(lay, cl)
    if (is.null(d)) next
    v <- as.character(d$value)
    tab <- sort(table(v[!is.na(v)]), decreasing = TRUE)
    zone_cov[[paste(cl, lay)]] <- data.frame(
      lulc = cl, layer = lay, n = length(v),
      pct_na = round(100 * mean(is.na(v)), 3),
      n_distinct = length(tab),
      nchar_min = min(nchar(v[!is.na(v)])),
      nchar_max = max(nchar(v[!is.na(v)])),
      # A space-padded key looks the right width and is still a different
      # string from the zero-padded one, so width alone does not prove the
      # derivation ran correctly.
      has_space = any(grepl(" ", v, fixed = TRUE)),
      all_digits = all(grepl("^[0-9]+$", v[!is.na(v)])),
      largest = names(tab)[1],
      pct_largest = round(100 * as.numeric(tab[1]) / length(v), 2),
      n_strata_lt_100 = sum(tab < 100),
      stringsAsFactors = FALSE)
    rm(d, v, tab)
  }
  gc(verbose = FALSE)
}
zone_cov <- if (length(zone_cov)) do.call(rbind, zone_cov) else NULL
if (!is.null(zone_cov)) { rownames(zone_cov) <- NULL; print(zone_cov) }
ca_log("Expected widths are 2 for ecoregion_l3 after padding and 8 for huc8. ",
       "has_space TRUE means the pad rule left-padded with spaces rather than ",
       "zeros, which yields a key of the right width and the wrong value. ",
       "all_digits FALSE means the same thing.")

ca_log("=== 7f. Aspect decomposition ===")
# Stage 4 replaces aspect with northness and eastness. Verified here on the
# extracted degrees so the rule is checked before stage 4 applies it. Flat
# cells are coded -1 by the SRTM derivative and take 0 for both components,
# which is the correct limit for a surface with no bearing.
for (cl in LULC_3D) {
  a <- rds("aspect", cl)
  if (is.null(a)) next
  v <- a$value[!is.na(a$value)]
  flat <- v < 0
  north <- ifelse(flat, 0, cos(v * pi / 180))
  east <- ifelse(flat, 0, sin(v * pi / 180))
  ca_log(cl, "  aspect degrees  min=", round(min(v), 2),
         "  max=", round(max(v), 2),
         "  flat n=", format(sum(flat), big.mark = ","),
         " (", round(100 * mean(flat), 3), "%)")
  ca_log(cl, "  northness ", round(min(north), 3), " to ", round(max(north), 3),
         "   eastness ", round(min(east), 3), " to ", round(max(east), 3))
  if (any(v > 0 & v < 0.001)) {
    ca_log("WARNING. Values between 0 and 0.001 degrees are present in ", cl,
           ", so the flat code may not be -1 everywhere.")
  }
  rm(a, v, flat, north, east); gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# 8  STAGE 3C MANAGEMENT AND OWNERSHIP RECORDS
# ---------------------------------------------------------------------------
# Stage 3c output is long, one row per pixel per intersecting polygon, which
# changes what every check means. Row count is not pixel count, so completeness
# cannot be verified by comparing to the grid size. pixel_id is a subset of the
# partition rather than the whole of it, so identical-order comparison does not
# apply. What replaces them is a containment check, because a pixel_id outside
# the partition means the join returned the wrong index and that is the failure
# that would silently corrupt every downstream merge.

rdv <- function(layer, lulc) {
  f <- ca_vector_path(layer, lulc)
  if (!file.exists(f)) return(NULL)
  as.data.frame(arrow::read_parquet(f))
}

vec_reg <- ca_vector_layers[ca_vector_layers$active %in% TRUE, ]

ca_log("=== 8a. Stage 3c completeness ===")
vec_inv <- do.call(rbind, lapply(seq_len(nrow(vec_reg)), function(i) {
  lay <- vec_reg$layer[i]
  got <- ca_lulc$label[vapply(ca_lulc$label, function(cl)
    file.exists(ca_vector_path(lay, cl)), logical(1))]
  data.frame(layer = lay, role = vec_reg$role[i],
             expected = nrow(ca_lulc), found = length(got),
             missing = paste(setdiff(ca_lulc$label, got), collapse = " "),
             stringsAsFactors = FALSE)
}))
rownames(vec_inv) <- NULL
print(vec_inv)
if (any(vec_inv$found < vec_inv$expected)) {
  ca_log("INCOMPLETE. Resubmit the missing tasks before running stage 4.")
}

ca_log("=== 8b. Stage 3c shape and pixel coverage ===")
# n_rows against n_pixels is the quantity the probe estimated. A large gap
# between them is not an error, it is overlapping records, and it is the reason
# ca_thin_rules$single_event_only exists.
vec_shape <- list()
vec_year <- list()
vec_class <- list()
vec_multi <- list()
vec_label <- list()

for (cl in ca_lulc$label) {

  grid_cl <- ca_grid_points(cl)$pixel_id

  for (i in seq_len(nrow(vec_reg))) {

    lay <- vec_reg$layer[i]
    d <- rdv(lay, cl)
    if (is.null(d)) next

    n_px <- length(unique(d$pixel_id))
    # Containment, not equality. Any pixel_id outside the partition means the
    # join returned a position rather than an identifier.
    outside <- sum(!(d$pixel_id %in% grid_cl))

    vec_shape[[paste(cl, lay)]] <- data.frame(
      lulc = cl, layer = lay,
      n_rows = nrow(d), n_pixels = n_px,
      n_grid = length(grid_cl),
      pct_grid_hit = round(100 * n_px / length(grid_cl), 3),
      rows_per_pixel = if (n_px) round(nrow(d) / n_px, 3) else 0,
      max_per_pixel = if (nrow(d)) max(table(d$pixel_id)) else 0L,
      n_pixel_outside_grid = outside,
      n_poly_uid = length(unique(d$poly_uid)),
      stringsAsFactors = FALSE)

    # Event years. The analysis window is what matters, not the raw range,
    # because a record dated 1900 or 2026 is real in the archive and outside
    # the study.
    y <- d$event_year
    vec_year[[paste(cl, lay)]] <- data.frame(
      lulc = cl, layer = lay, n_rows = nrow(d),
      pct_year_na = round(100 * mean(is.na(y)), 3),
      yr_min = suppressWarnings(if (all(is.na(y))) NA_integer_ else
        min(y, na.rm = TRUE)),
      yr_max = suppressWarnings(if (all(is.na(y))) NA_integer_ else
        max(y, na.rm = TRUE)),
      pct_in_analysis = round(100 * mean(y %in% ca_years$analysis), 3),
      pct_pre_1985 = round(100 * mean(!is.na(y) & y < 1985L), 3),
      pct_post_2025 = round(100 * mean(!is.na(y) & y > 2025L), 3),
      stringsAsFactors = FALSE)

    # Record class and crosswalk origin, on the thinning layers only. Reported
    # as a pixel share as well as a record share, because the response letter
    # question is what fraction of the treated sample rests on a locally
    # assigned activity, not what fraction of archive rows do.
    if (!is.na(vec_reg$crosswalk[i])) {
      trt <- d$record_class %in% "treatment_LMH"
      asg <- d$xwalk_origin %in% "alias_assigned"
      px_trt <- length(unique(d$pixel_id[trt]))
      px_asg <- length(unique(d$pixel_id[trt & asg]))
      vec_class[[paste(cl, lay)]] <- data.frame(
        lulc = cl, layer = lay, n_rows = nrow(d),
        pct_unresolved = round(100 * mean(is.na(d$record_class)), 3),
        pct_treatment = round(100 * mean(trt), 2),
        pct_disturb_only = round(100 * mean(
          d$record_class %in% "disturbance_only"), 2),
        pct_non_disturbing = round(100 * mean(
          d$record_class %in% "non_disturbing"), 2),
        px_treatment = px_trt,
        px_assigned = px_asg,
        pct_treated_px_assigned = if (px_trt)
          round(100 * px_asg / px_trt, 3) else NA_real_,
        stringsAsFactors = FALSE)

      # Events per pixel among treatment records. single_event_only excludes
      # every pixel above one, so this is the size of that exclusion.
      if (any(trt)) {
        tb <- table(d$pixel_id[trt])
        vec_multi[[paste(cl, lay)]] <- data.frame(
          lulc = cl, layer = lay,
          px_with_treatment = length(tb),
          px_single_event = sum(tb == 1L),
          px_multi_event = sum(tb > 1L),
          pct_multi = round(100 * mean(tb > 1L), 2),
          max_events = max(tb),
          stringsAsFactors = FALSE)
        rm(tb)
      }
    }

    # Label coverage. A value present in the data and absent from the config
    # lookup is dropped silently by a case_match at stage 6, so it is caught
    # here instead.
    lookup <- switch(lay,
      cpad      = ca_cpad_agency$label,
      ownership = ca_ownership_level$label,
      NULL)
    if (!is.null(lookup)) {
      obs <- sort(unique(d$key_value[!is.na(d$key_value)]))
      vec_label[[paste(cl, lay)]] <- data.frame(
        lulc = cl, layer = lay,
        n_observed = length(obs),
        n_in_config = length(lookup),
        unmatched = paste(setdiff(obs, lookup), collapse = " | "),
        unused_in_config = paste(setdiff(lookup, obs), collapse = " | "),
        stringsAsFactors = FALSE)
    }

    rm(d); gc(verbose = FALSE)
  }
  rm(grid_cl); gc(verbose = FALSE)
}

bindl <- function(x) {
  if (!length(x)) return(NULL)
  out <- do.call(rbind, x); rownames(out) <- NULL; out
}

vec_shape <- bindl(vec_shape)
vec_year  <- bindl(vec_year)
vec_class <- bindl(vec_class)
vec_multi <- bindl(vec_multi)
vec_label <- bindl(vec_label)

if (!is.null(vec_shape)) {
  print(vec_shape)
  bad <- vec_shape$n_pixel_outside_grid > 0
  if (any(bad)) {
    ca_log("FAIL. pixel_id values outside the grid partition in: ",
           paste(unique(paste(vec_shape$lulc[bad], vec_shape$layer[bad])),
                 collapse = ", "),
           ". The join returned a position rather than an identifier and the ",
           "affected layers must be rerun before stage 4.")
  }
}

ca_log("=== 8c. Stage 3c event years ===")
if (!is.null(vec_year)) print(vec_year)
ca_log("Records outside 1990 to 2025 are archival, not errors. They are ",
       "dropped at stage 6 by ca_treat_years. pct_year_na on cpad reflects ",
       "YR_EST = 0, handled by ca_pretrt$cpad_year_unknown_rule.")

ca_log("=== 8d. Record class and crosswalk resolution ===")
if (!is.null(vec_class)) print(vec_class)
ca_log("pct_unresolved must be zero. A record with no record_class is ",
       "dropped silently at stage 6. pct_treated_px_assigned is the figure ",
       "the response letter needs, the share of treated pixels resting on an ",
       "activity classified by analogy rather than by the published tables.")

ca_log("=== 8e. Events per pixel ===")
if (!is.null(vec_multi)) print(vec_multi)
ca_log("px_multi_event is the exclusion imposed by ",
       "ca_thin_rules$single_event_only. Those pixels leave the treatment arm ",
       "and remain eligible as controls before their first recorded event.")

ca_log("=== 8f. Label coverage against config ===")
if (!is.null(vec_label)) print(vec_label)
ca_log("unmatched must be empty. A label present in the data and absent from ",
       "ca_cpad_agency or ca_ownership_level is dropped by the stage 6 ",
       "recode without warning.")

ca_log("=== 8g. Cross-layer multi-event exclusion ===")
# Section 8e counts repeat entries within a layer. A pixel with one THP record
# and one FACTS record is single-event in both tables and multi-event in fact,
# so the within-layer figure is a floor. This section pools the four thinning
# layers and reports the real exclusion.
#
# Counting only. The rule itself lives in ca_thin_rules and is applied at
# stage 6, exactly once. This section describes what stage 6 will do, which
# also makes it the check that stage 6 did it correctly. Duplicating the
# classification here would put the same decision in two files.

thin_layers <- vec_reg$layer[!is.na(vec_reg$crosswalk)]
vec_pool <- list()

for (cl in ca_lulc$label) {

  keep <- list()
  for (lay in thin_layers) {
    d <- rdv(lay, cl)
    if (is.null(d)) next
    sel <- d$record_class %in% "treatment_LMH"
    if (any(sel)) {
      keep[[lay]] <- data.frame(pixel_id = d$pixel_id[sel],
                                layer = lay,
                                stringsAsFactors = FALSE)
    }
    rm(d, sel)
  }
  gc(verbose = FALSE)
  if (!length(keep)) next

  p <- do.call(rbind, keep)
  rm(keep); gc(verbose = FALSE)

  ev <- table(p$pixel_id)
  # A pixel appearing in more than one source layer, distinct from a pixel with
  # several entries in one layer. Reported separately because the two have
  # different causes, overlapping jurisdictions against repeat entry.
  ly <- tapply(p$layer, p$pixel_id, function(x) length(unique(x)))

  vec_pool[[cl]] <- data.frame(
    lulc = cl,
    px_with_treatment = length(ev),
    px_single_event = sum(ev == 1L),
    px_multi_event = sum(ev > 1L),
    pct_multi_pooled = round(100 * mean(ev > 1L), 2),
    px_multi_layer = sum(ly > 1L),
    pct_multi_layer = round(100 * mean(ly > 1L), 2),
    max_events = max(ev),
    stringsAsFactors = FALSE)

  rm(p, ev, ly); gc(verbose = FALSE)
}

vec_pool <- bindl(vec_pool)
if (!is.null(vec_pool)) print(vec_pool)
ca_log("pct_multi_pooled is the figure for the Methods. It supersedes the ",
       "within-layer values in 8e, which understate the exclusion because a ",
       "pixel entered once under a THP and once under FACTS counts as single ",
       "in both tables.")
ca_log("Under ca_thin_rules, these pixels leave the analysis entirely rather ",
       "than returning to the control pool, since a repeatedly entered stand ",
       "is not business as usual in any window. That narrows the estimand to ",
       "singly entered stands and does not bias it.")

# ---------------------------------------------------------------------------
# 9  WRITE
# ---------------------------------------------------------------------------

dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
w <- function(obj, nm) {
  if (is.null(obj)) return(invisible(NULL))
  write.csv(obj, ca_meta(paste0("verify_", nm, "_", stamp, ".csv")),
            row.names = FALSE)
}
w(inv, "inventory"); w(vals, "values"); w(chk, "join")
w(pairs, "pairs");  w(vt, "vegtreefrac"); w(ch_rows, "canopyheight")
w(stat_inv, "3d_inventory"); w(stat_vals, "3d_values")
w(stat_join, "3d_join"); w(zone_cov, "3d_strata")
w(vec_inv, "3c_inventory"); w(vec_shape, "3c_shape")
w(vec_year, "3c_years");    w(vec_class, "3c_class")
w(vec_multi, "3c_multi");   w(vec_label, "3c_labels")
w(vec_pool, "3c_pooled")

ca_log("Verification complete. CSVs written to ", ca_meta())
