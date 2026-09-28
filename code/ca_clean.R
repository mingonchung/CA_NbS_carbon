### CA carbon revision pipeline
### Stage 4 shared module. Unit conversion, cleaning rules, and reads.
###
### DESIGN
###
### Stage 4 does not materialise a converted copy of the carbon stacks. Trace
### what actually consumes carbon values and the reason is plain.
###
###   stage 5b   the two screen layers, full grid, all years
###   stage 6    nothing
###   stage 7    nothing
###   stage 8    gpp_pre5, agb_pre5, ch_pre5 for the sampled pool only
###   stage 9    all outcomes for matched pixels only
###
### Stages 5, 6, and 7 never read GPP, AGB, NPP, NEP, NBP, eMapR, or LEMMA. A
### full-grid converted panel would double roughly a hundred gigabytes for the
### sole purpose of being subset twice. Conversion is therefore a read-time
### operation living in one function, and the arithmetic runs more than once.
### Two float32 multiplies are free relative to the read.
###
### The rule that matters is that no caller applies a factor of its own. Every
### stage that touches a carbon value goes through ca_read_layer(), so a
### conversion cannot be applied twice or forgotten once.
###
### Paths go through ca_extract_store(layer). The extract directory is set by
### the SOURCE constant of the script that ran the registry, not by the per-row
### source column, and the two differ for eMapR, LEMMA, and all three screen
### layers. Nothing in this module builds a path any other way.
###
### WHAT IS APPLIED, AND IN WHICH ORDER
###
###   1. ca_almanac_layers$scale or ca_ncsda_layers$scale, at extraction.
###      Already in the stored value. Fire_LCP_CH is therefore metres on disk,
###      not decimetres, which is why ca_unit_conversion carries factor 1 for
###      it.
###   2. ca_unit_conversion$factor, here.
###   3. ca_unit_conversion$carbon_fraction, here, where TRUE. 0.47.
###   4. ca_clean_rules, here.
###
### Zero is never mapped to NA. The submitted C.sub[C.sub == 0] <- NA erased
### real values from NEP, NBP, the two Disturbance layers, and Fire_LCP_CH.
### LEMMA zero is a genuine sentinel and is removed at extraction through
### ca_sources, not here.

ca_require(c("arrow"))

# ---------------------------------------------------------------------------
# UNIT CONVERSION
# ---------------------------------------------------------------------------

ca_carbon_fraction <- 0.47

ca_convert <- function(x, layer) {

  i <- match(layer, ca_unit_conversion$layer)
  if (is.na(i)) {
    stop("No unit conversion entry for layer '", layer, "'. Every layer read ",
         "through this module must be declared in ca_unit_conversion, so a ",
         "missing entry is a config error rather than a reason to pass the ",
         "values through unchanged.")
  }

  f <- ca_unit_conversion$factor[i]
  if (!is.na(f) && f != 1) x <- x * f
  if (isTRUE(ca_unit_conversion$carbon_fraction[i])) x <- x * ca_carbon_fraction
  x
}

ca_target_units <- function(layer) {
  ca_unit_conversion$target[match(layer, ca_unit_conversion$layer)]
}

# ---------------------------------------------------------------------------
# CLEANING RULES
# ---------------------------------------------------------------------------
# ref is the reference vector a rule needs. Only mask_by_agb takes one, and it
# takes Carbon_AGB for the same year and class, aligned on pixel_id by the
# caller. Passing it as an argument rather than reading it here keeps the rule
# a pure function and keeps the read count visible at the call site.
#
# Returns the cleaned vector with an n_affected attribute per rule, so the
# verifier can report counts without recomputing them.

# Applied to STORED values, before conversion, because a type limit lives in
# stored counts. The extraction scale is already in the stored value, so the
# comparison target is value * scale from ca_sentinel_rules.
ca_mask_sentinel <- function(x, layer) {

  r <- ca_sentinel_rules[ca_sentinel_rules$layer == layer, , drop = FALSE]
  if (!nrow(r)) return(structure(x, n_sentinel = 0L))

  n <- 0L
  for (k in seq_len(nrow(r))) {
    target <- r$value[k] * r$scale[k]
    hit <- which(!is.na(x) & abs(x - target) <= abs(target) * 1e-6 + 1e-9)
    x[hit] <- NA_real_
    n <- n + length(hit)
  }
  structure(x, n_sentinel = n)
}

ca_clean <- function(x, layer, ref = NULL) {

  rules <- ca_clean_rules[ca_clean_rules$layer == layer, , drop = FALSE]
  if (!nrow(rules)) return(structure(x, affected = integer(0)))

  affected <- integer(nrow(rules))
  names(affected) <- rules$rule

  for (k in seq_len(nrow(rules))) {

    rule <- rules$rule[k]
    val  <- rules$value[k]

    if (rule == "clip_upper") {
      hit <- which(!is.na(x) & x > val)
      x[hit] <- val

    } else if (rule == "clip_lower") {
      hit <- which(!is.na(x) & x < val)
      x[hit] <- val

    } else if (rule == "mask_by_agb") {
      if (is.null(ref)) {
        stop("Rule mask_by_agb on layer '", layer, "' needs the Carbon_AGB ",
             "vector for the same year and class. Pass it as ref.")
      }
      if (length(ref) != length(x)) {
        stop("ref length ", length(ref), " does not match value length ",
             length(x), ". The two must be aligned on pixel_id.")
      }
      hit <- which(is.na(ref) & !is.na(x))
      x[hit] <- NA_real_

    } else {
      stop("Unknown cleaning rule '", rule, "' for layer '", layer, "'")
    }

    affected[k] <- length(hit)
  }

  structure(x, affected = affected)
}

# ---------------------------------------------------------------------------
# READS
# ---------------------------------------------------------------------------
# One shard is pixel_id int32 plus value float32, sorted on pixel_id. When
# pixel_id is supplied the filter is pushed into the Arrow scan, so only the
# requested rows are decompressed. That is what makes the lazy design cheap at
# stages 8 and 9, where the pixel set is a sample rather than the grid.
#
# The returned value is always aligned to the requested pixel_id, including
# positions the shard does not carry, which come back NA. Alignment is by
# match(), never by row order.

ca_align <- function(tbl, pixel_id) {

  if (is.null(pixel_id)) {
    return(list(pixel_id = tbl$pixel_id, value = as.numeric(tbl$value)))
  }

  # Fast path. Two shards written from the same sorted grid partition carry
  # identical pixel_id vectors, which is the common case at full-grid width.
  # identical() is a cheap comparison; match() on 95.6 million against 95.6
  # million builds a hash table of the same length and is the single largest
  # allocation in this module. Correctness is unchanged because the equality is
  # tested rather than assumed.
  if (identical(tbl$pixel_id, pixel_id)) {
    return(list(pixel_id = pixel_id, value = as.numeric(tbl$value)))
  }

  out <- rep(NA_real_, length(pixel_id))
  j <- match(tbl$pixel_id, pixel_id)
  ok <- !is.na(j)
  out[j[ok]] <- as.numeric(tbl$value)[ok]
  list(pixel_id = pixel_id, value = out)
}

ca_read_shard <- function(path, pixel_id = NULL) {

  if (!file.exists(path)) stop("Missing extract: ", path)

  if (is.null(pixel_id)) {
    tbl <- arrow::read_parquet(path, col_select = c("pixel_id", "value"),
                               as_data_frame = TRUE)
  } else {
    ids <- as.integer(pixel_id)
    ds <- arrow::open_dataset(path)
    tbl <- as.data.frame(
      dplyr::collect(dplyr::select(dplyr::filter(ds, pixel_id %in% ids),
                                   pixel_id, value))
    )
    rm(ds, ids)
  }
  tbl
}

# Arrow allocates outside the R heap, so gc() does not report or reclaim it and
# a long section can exhaust the node while R believes it has room. This is
# what an allocation failure on a 364 Mb vector at 128G actually means. Called
# between sections and inside the long loops.
ca_gc <- function(label = NULL) {
  gc(full = TRUE)
  if (!is.null(label) && requireNamespace("arrow", quietly = TRUE)) {
    pool <- try(arrow::default_memory_pool(), silent = TRUE)
    if (!inherits(pool, "try-error")) {
      ca_log(label, " | R ", round(sum(gc()[, 2]) / 1024, 2),
             " GB, arrow ", round(pool$bytes_allocated / 1e9, 2), " GB")
    }
  }
  invisible(TRUE)
}

# clean = FALSE returns converted but unrepaired values, which is what the
# verifier needs in order to count what the rules would have changed.
#
# source is retained in the signature because it names the product a caller
# believes it is reading, and a mismatch against the registry is worth catching
# rather than ignoring. The path itself comes from the layer.
# convert = FALSE returns stored values with no unit conversion at all, which
# is what a sentinel check needs. A type limit lives in stored counts, and
# comparing it against a value already multiplied by 100 and 0.47 is a test
# that can never fire.
ca_read_layer <- function(source, layer, year, lulc,
                          pixel_id = NULL, clean = TRUE, convert = TRUE) {

  reg_src <- ca_layer_registry()$source[
    match(layer, ca_layer_registry()$layer)]
  if (!is.na(reg_src) && !identical(reg_src, source)) {
    stop("Layer '", layer, "' belongs to source '", reg_src,
         "', not '", source, "'.")
  }

  path <- ca_extract_path(ca_extract_store(layer), layer, year, lulc)
  tbl  <- ca_read_shard(path, pixel_id)
  al   <- ca_align(tbl, pixel_id)

  # Order matters. Sentinel first, on stored values, then convert, then the
  # remaining rules on converted values. Reversing the first two would compare
  # a type limit against a number the data can no longer hold, and would leave
  # a missing-data marker to be clipped into a real value.
  v <- al$value
  n_sent <- 0L
  if (clean) {
    v <- ca_mask_sentinel(v, layer)
    n_sent <- attr(v, "n_sentinel")
    v <- as.numeric(v)
  }

  if (convert) v <- ca_convert(v, layer)

  if (!clean) {
    return(list(pixel_id = al$pixel_id, value = v,
                affected = integer(0), n_sentinel = 0L))
  }

  ref <- NULL
  if (any(ca_clean_rules$layer == layer &
          ca_clean_rules$rule == "mask_by_agb")) {
    ra <- ca_read_shard(
      ca_extract_path(ca_extract_store("Carbon_AGB"), "Carbon_AGB", year,
                      lulc), pixel_id)
    ref <- ca_align(ra, al$pixel_id)$value
  }

  v <- ca_clean(v, layer, ref = ref)
  list(pixel_id = al$pixel_id, value = as.numeric(v),
       affected = attr(v, "affected"), n_sentinel = n_sent)
}

# Class codes, read without conversion or cleaning. MTBS severity is the only
# member today. Values come back as integers rather than doubles, because a
# severity class compared with == against a float is a bug waiting for a
# floating point representation to change.
ca_read_categorical <- function(layer, year, lulc, pixel_id = NULL) {

  if (!ca_is_categorical(layer)) {
    stop("Layer '", layer, "' is not categorical. Read it with ",
         "ca_read_layer() so unit conversion and cleaning are applied.")
  }

  path <- ca_extract_path(ca_extract_store(layer), layer, year, lulc)
  al <- ca_align(ca_read_shard(path, pixel_id), pixel_id)
  list(pixel_id = al$pixel_id, value = as.integer(round(al$value)))
}

# Wide read. Rows are pixel_id, columns are years. Used by stage 8 for the
# pre-treatment windows and by stage 9 for the DiD panel, both against a small
# pixel set. Calling it on the full grid is possible and is not the intended
# use.
ca_read_years <- function(source, layer, years, lulc,
                          pixel_id, clean = TRUE) {

  if (is.null(pixel_id)) {
    stop("ca_read_years() requires an explicit pixel_id set. A full-grid wide ",
         "read is not the intended use of this module.")
  }

  m <- matrix(NA_real_, nrow = length(pixel_id), ncol = length(years),
              dimnames = list(NULL, as.character(years)))

  for (k in seq_along(years)) {
    m[, k] <- ca_read_layer(source, layer, years[k], lulc,
                            pixel_id = pixel_id, clean = clean)$value
  }
  m
}

ca_read_outcome <- function(outcome, lulc, pixel_id, years = NULL,
                            clean = TRUE) {
  ol <- ca_outcome_layer(outcome)
  if (is.null(years)) years <- ca_outcome_years(outcome)
  ca_read_years(ol$source, ol$layer, years, lulc, pixel_id, clean = clean)
}

# ---------------------------------------------------------------------------
# STATIC READS
# ---------------------------------------------------------------------------
# Stage 3d output. One row per pixel, no year. Zone layers carry character
# codes and must not be coerced to numeric, since that drops leading zeros and
# turns one exact-matching stratum into another.

ca_read_static <- function(layer, lulc, pixel_id = NULL) {

  path <- ca_static_path(layer, lulc)
  if (!file.exists(path)) stop("Missing static extract: ", path)

  if (is.null(pixel_id)) {
    tbl <- as.data.frame(arrow::read_parquet(path))
  } else {
    ids <- as.integer(pixel_id)
    tbl <- as.data.frame(
      dplyr::collect(dplyr::filter(arrow::open_dataset(path),
                                   pixel_id %in% ids))
    )
  }

  is_zone <- layer %in% ca_zone_layers$layer

  if (is.null(pixel_id)) {
    v <- if (is_zone) as.character(tbl$value) else as.numeric(tbl$value)
    return(list(pixel_id = tbl$pixel_id, value = v))
  }

  out <- if (is_zone) rep(NA_character_, length(pixel_id)) else
    rep(NA_real_, length(pixel_id))
  j <- match(tbl$pixel_id, pixel_id)
  ok <- !is.na(j)
  src <- if (is_zone) as.character(tbl$value) else as.numeric(tbl$value)
  out[j[ok]] <- src[ok]
  list(pixel_id = pixel_id, value = out)
}

# ---------------------------------------------------------------------------
# DERIVED COVARIATES
# ---------------------------------------------------------------------------
# Aspect in degrees is circular and cannot enter matching. 350 and 10 are 340
# apart numerically and 20 apart physically, which breaks both the standardised
# mean difference and the caliper at the wrap. Both components are returned in
# [-1, 1].
#
# Flat cells are coded -1 by the SRTM derivative. A flat cell has no bearing,
# so it is neither missing nor -1 degrees. Applying the formula would return
# northness 0.9998 and enter every flat pixel as due north. Both components are
# set to 0, which is the correct limit and the centroid of all bearings, and
# the pixel is retained.

ca_derive_aspect <- function(aspect) {

  flat_code <- ca_derive_rules$flat_code[1]
  flat_val  <- ca_derive_rules$flat_value[1]

  rad <- aspect * pi / 180
  northness <- cos(rad)
  eastness  <- sin(rad)

  flat <- !is.na(aspect) & aspect == flat_code
  northness[flat] <- flat_val
  eastness[flat]  <- flat_val

  list(northness = northness, eastness = eastness, n_flat = sum(flat))
}

invisible(TRUE)
