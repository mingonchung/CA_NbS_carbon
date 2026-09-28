### CA carbon revision pipeline
### Stage 3 shared extraction engine
###
### Sourced by 3a_extract_carbon_2026.R, 3a_extract_carbon_2022.R, and
### 3b_extract_screens_2026.R. Holds the logic those three share so that CRS
### handling, nodata handling, and output schema cannot drift between them.
###
### Rules enforced here.
###
### 1. Rasters are never reprojected. Point coordinates are transformed into the
###    raster's native CRS instead. In this pipeline no transform is actually
###    needed, because NLCD 2001, the Almanac, eMapR, and LEMMA all sit on NAD83
###    Conus Albers, but the check is enforced rather than assumed.
###
### 2. Zero is a valid value except where a layer declares nodata = 0. Only
###    Fire_LCP does, following the FARSITE convention. Everywhere else zero
###    must survive: NEP and NBP are signed and cross zero, and in the
###    disturbance layers zero means a valid pixel with no disturbance that
###    year. The submitted pipeline set every zero to NA, which destroyed both.
###
### 3. Year is parsed from the filename. The submitted pipeline inferred it from
###    alphabetical file order with `i + 1984`, which shifts every subsequent
###    year whenever a file is added, renamed, or missing.
###
### 4. Extraction runs in point chunks, so peak memory is set by chunk size
###    rather than class size. Evergreen alone carries 95.6 million points.

ca_require(c("terra", "arrow"))

# ---------------------------------------------------------------------------
# FILE DISCOVERY
# ---------------------------------------------------------------------------

ca_find_rasters <- function(dir, pattern = "\\.(tif|tiff|vrt)$") {
  if (!dir.exists(dir)) stop("Layer directory not found: ", dir)
  files <- list.files(dir, pattern = pattern, recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No rasters found under ", dir)
  files
}

# Extracts the first plausible four-digit year from each basename. Restricted to
# the study window, so a resolution or version number cannot be mistaken for a
# year. Handles mr200_1990_mn.tif and conus_biomass_yr1990_prj.tif alike.
ca_parse_year <- function(files, valid = ca_years$study) {
  base <- basename(files)
  hits <- regmatches(base, gregexpr("(19|20)[0-9]{2}", base))
  vapply(hits, function(h) {
    y <- as.integer(h)
    y <- y[y %in% valid]
    if (length(y) == 0L) NA_integer_ else y[1]
  }, integer(1))
}

# Returns file and year for one layer, sorted by year, after checking the
# sequence for duplicates and gaps.
ca_layer_index <- function(dir, first_year, last_year, pattern = NULL) {

  files <- if (is.null(pattern)) ca_find_rasters(dir) else
    ca_find_rasters(dir, pattern)

  years <- ca_parse_year(files)
  keep <- !is.na(years) & years >= first_year & years <= last_year

  idx <- data.frame(file = files[keep], year = years[keep],
                    stringsAsFactors = FALSE)

  if (!nrow(idx)) {
    stop("No rasters with a year in ", first_year, " to ", last_year,
         " under ", dir)
  }

  idx <- idx[order(idx$year), ]

  if (anyDuplicated(idx$year)) {
    stop("Duplicate years under ", dir, ": ",
         paste(unique(idx$year[duplicated(idx$year)]), collapse = ", "))
  }

  gap_years <- setdiff(seq(first_year, last_year), idx$year)
  if (length(gap_years)) {
    warning("Missing years under ", dir, ": ",
            paste(gap_years, collapse = ", "), call. = FALSE)
  }

  rownames(idx) <- NULL
  idx
}

# ---------------------------------------------------------------------------
# GRID ACCESS
# ---------------------------------------------------------------------------
# Reads one LULC partition of the frozen grid. Only pixel_id, x, and y are
# pulled, since lulc_code is implied by the partition. No geometry object is
# built at any point, because terra::extract and terra::project both take
# plain coordinate matrices.

ca_grid_points <- function(lulc) {

  part <- file.path(ca_grid$parquet_dir, paste0("lulc=", lulc))
  if (!dir.exists(part)) stop("Grid partition not found: ", part)

  ds <- arrow::open_dataset(part)

  tbl <- try(
    as.data.frame(
      arrow::Scanner$create(ds, projection = c("pixel_id", "x", "y"))$ToTable()
    ),
    silent = TRUE
  )
  if (inherits(tbl, "try-error")) {
    tbl <- as.data.frame(ds)[, c("pixel_id", "x", "y")]
  }

  # Stage 1 wrote blocks in ascending row order, so the partition is already
  # sorted. Sorting 95.6 million rows again would cost a needless copy.
  if (is.unsorted(tbl$pixel_id)) tbl <- tbl[order(tbl$pixel_id), ]

  list(pixel_id = tbl$pixel_id, xy = cbind(tbl$x, tbl$y))
}

# Transforms a coordinate matrix into a target CRS, returning it unchanged when
# the CRS already matches.
#
# Comparison uses terra::same.crs rather than string equality. Every product
# here sits on NAD83 Conus Albers but declares it differently: the Almanac as
# EPSG:5070, eMapR as Albers_Conic_Equal_Area with a Custom authority, LEMMA as
# National_Albers with a Custom authority, NLCD as a proj4 string. String
# comparison would see four different CRSs and run a pointless projection that
# introduces floating-point noise into coordinates that are already exact.
ca_xy_to <- function(xy, from_crs, to_crs) {
  if (terra::same.crs(from_crs, to_crs)) return(xy)
  terra::project(xy, from = from_crs, to = to_crs)
}

# ---------------------------------------------------------------------------
# EXTRACTION
# ---------------------------------------------------------------------------

ca_extract_points <- function(rast_path, xy, band = NA_integer_,
                              nodata = NA_real_, scale = NA_real_,
                              method = "simple",
                              chunk = ca_extract_chunk) {

  r <- terra::rast(rast_path)

  if (!is.na(band)) {
    if (band > terra::nlyr(r)) {
      stop("Band ", band, " requested but ", basename(rast_path),
           " has ", terra::nlyr(r), " layers")
    }
    r <- r[[band]]
  } else if (terra::nlyr(r) > 1L) {
    stop("Multi-band raster with no band specified: ", basename(rast_path))
  }

  if (!nzchar(terra::crs(r))) stop("Raster has no CRS: ", rast_path)

  # The sentinel is applied to the raster before extraction, not to the
  # extracted vector afterwards. Under bilinear interpolation a cell adjacent
  # to a sentinel is blended with it, producing a plausible-looking value that
  # never equals the sentinel and so survives any post-hoc equality test. A
  # travel time cell next to a 65535 sea cell would come back as roughly
  # 16,000 minutes and enter matching as real. Setting NAflag makes terra treat
  # those cells as missing during interpolation instead.
  if (!is.na(nodata)) terra::NAflag(r) <- nodata

  ref <- ca_ref_crs()
  transformed <- !terra::same.crs(ref, terra::crs(r))
  xy_r <- if (transformed) terra::project(xy, from = ref, to = terra::crs(r)) else xy

  n <- nrow(xy_r)
  out <- rep(NA_real_, n)

  for (s in seq(1L, n, by = chunk)) {
    e <- min(s + chunk - 1L, n)
    out[s:e] <- terra::extract(r, xy_r[s:e, , drop = FALSE],
                               method = method)[[1]]
  }
  rm(xy_r)

  # Backstop only. NAflag above has already removed the sentinel in almost
  # every case. This catches a file whose driver refuses the flag, and it is a
  # no-op otherwise. See rule 2 in the header: this is the only place zero may
  # become NA, and only when a layer explicitly declares nodata = 0.
  if (!is.na(nodata)) out[out == nodata] <- NA_real_

  if (!is.na(scale)) out <- out * scale

  attr(out, "transformed") <- transformed
  out
}

# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------
# One Parquet file per layer, year, and class. pixel_id is stored beside the
# value rather than relying on row order, so a downstream join cannot be
# silently misaligned. Sorted pixel_id compresses to near nothing under zstd.
# Values are float32, which is well inside the precision of every product here
# and halves the panel size.

ca_extract_path <- function(source, layer, year, lulc) {
  dir <- ca_work("extract", source, layer)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("%s_%d_%s.parquet", layer, year, lulc))
}

ca_write_extract <- function(pixel_id, value, source, layer, year, lulc) {
  path <- ca_extract_path(source, layer, year, lulc)
  tbl <- arrow::arrow_table(
    pixel_id = arrow::Array$create(as.integer(pixel_id), type = arrow::int32()),
    value    = arrow::Array$create(as.numeric(value), type = arrow::float32())
  )
  arrow::write_parquet(tbl, path, compression = ca_io$compression)
  path
}

# ---------------------------------------------------------------------------
# TASK TABLE
# ---------------------------------------------------------------------------
# Expands a layer registry into one array task per layer and year. Inactive
# layers are dropped before expansion, so deactivating a layer shrinks the array
# rather than leaving holes in it. Written to meta/ so the task-to-layer mapping
# is recoverable after the fact.

ca_task_table <- function(registry, source, write = TRUE) {

  reg <- registry[registry$active %in% TRUE, , drop = FALSE]
  if (!nrow(reg)) stop("No active layers in the registry for ", source)

  tasks <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
    data.frame(
      layer = reg$layer[i],
      band  = reg$band[i],
      year  = seq(reg$first_year[i], reg$last_year[i]),
      stringsAsFactors = FALSE
    )
  }))

  tasks <- tasks[order(tasks$layer, tasks$year), ]
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL

  if (write) {
    dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
    write.csv(tasks, ca_meta(paste0("task_table_", source, ".csv")),
              row.names = FALSE)
  }

  tasks
}

# ---------------------------------------------------------------------------
# PROBE
# ---------------------------------------------------------------------------
# Inventories a registry without extracting anything. Reports file count, year
# span, gaps, band count, datatype, NA flag, file-carried scale and offset, CRS
# match, and grid alignment. Run once per source before the first array.
#
# grid_aligned tests whether a product shares the reference resolution and sits
# on the reference grid origin. A product that shares the CRS but sits on a
# shifted origin was resampled at some point and no longer lines up pixel for
# pixel. That is the specific risk for eMapR, whose files carry a _prj suffix
# from an earlier projectRaster step that used bilinear interpolation.

ca_probe_layers <- function(registry, source) {

  ref <- if (file.exists(ca_grid$meta_rds)) readRDS(ca_grid$meta_rds) else NULL

  rows <- lapply(seq_len(nrow(registry)), function(i) {

    lay <- registry$layer[i]
    dir <- file.path(ca_input(registry$subdir[i]), registry$dir[i])

    blank <- function(msg) {
      data.frame(layer = lay, status = msg, dir = dir,
                 active = registry$active[i], stringsAsFactors = FALSE)
    }

    if (!dir.exists(dir)) return(blank("directory missing"))

    pat <- registry$file_pattern[i]
    idx <- try(ca_layer_index(dir, registry$first_year[i],
                              registry$last_year[i],
                              pattern = if (is.na(pat)) NULL else pat),
               silent = TRUE)
    if (inherits(idx, "try-error")) return(blank("no rasters found"))

    r <- terra::rast(idx$file[1])
    so <- terra::scoff(r)

    same_crs <- terra::same.crs(ca_ref_crs(), terra::crs(r))

    aligned <- NA
    if (!is.null(ref) && same_crs) {
      res_ok <- isTRUE(all.equal(terra::res(r), ref$resolution,
                                 tolerance = 1e-6))
      e <- as.vector(terra::ext(r))
      offx <- ((e[1] - ref$extent[1]) / terra::res(r)[1]) %% 1
      offy <- ((e[3] - ref$extent[3]) / terra::res(r)[2]) %% 1
      aligned <- res_ok &&
        min(offx, 1 - offx) < 1e-3 && min(offy, 1 - offy) < 1e-3
    }

    data.frame(
      layer        = lay,
      status       = "ok",
      dir          = dir,
      active       = registry$active[i],
      n_files      = nrow(idx),
      first_year   = min(idx$year),
      last_year    = max(idx$year),
      gaps         = paste(setdiff(seq(min(idx$year), max(idx$year)), idx$year),
                           collapse = " "),
      n_bands      = terra::nlyr(r),
      band_names   = paste(names(r), collapse = " "),
      datatype     = terra::datatype(r)[1],
      na_flag      = paste(terra::NAflag(r), collapse = " "),
      scale_file   = if (is.null(so)) NA_real_ else so[1, 1],
      offset_file  = if (is.null(so)) NA_real_ else so[1, 2],
      res_x        = terra::res(r)[1],
      same_crs     = same_crs,
      grid_aligned = aligned,
      cfg_nodata   = registry$nodata[i],
      cfg_scale    = registry$scale[i],
      crs_proj     = terra::crs(r, proj = TRUE),
      stringsAsFactors = FALSE
    )
  })

  cols <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))

  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  path <- ca_meta(paste0("probe_", source, "_",
                         format(Sys.Date(), "%Y%m%d"), ".csv"))
  write.csv(out, path, row.names = FALSE)
  ca_log("Probe written to ", path)

  out
}

# ---------------------------------------------------------------------------
# DRIVER
# ---------------------------------------------------------------------------
# Runs one array task, meaning one layer and one year across all three forest
# classes. Existing outputs are skipped, so a partially failed array can be
# resubmitted without recomputing completed tasks.

ca_run_extract_task <- function(registry, source, task_id, overwrite = FALSE) {

  tasks <- ca_task_table(registry, source, write = FALSE)
  if (is.na(task_id) || task_id < 1 || task_id > nrow(tasks)) {
    stop("task_id ", task_id, " outside 1 to ", nrow(tasks))
  }

  tk <- tasks[task_id, ]
  reg <- registry[registry$layer == tk$layer, ][1, ]

  dir <- file.path(ca_input(reg$subdir), reg$dir)
  idx <- ca_layer_index(dir, reg$first_year, reg$last_year,
                        pattern = if (is.na(reg$file_pattern)) NULL else
                          reg$file_pattern)

  hit <- idx[idx$year == tk$year, ]
  if (!nrow(hit)) stop("No raster for ", tk$layer, " ", tk$year)

  # Layer-level nodata and scale win. ca_sources supplies the fallback.
  src_row <- ca_sources[ca_sources$source == reg$source, ]
  if (!nrow(src_row)) stop("Unknown source in registry: ", reg$source)

  nodata <- if (!is.na(reg$nodata)) reg$nodata else src_row$nodata[1]
  scale  <- if (!is.na(reg$scale))  reg$scale  else src_row$scale[1]

  ca_log("Task ", task_id, "  ", tk$layer, " ", tk$year,
         "  file=", basename(hit$file),
         "  nodata=", nodata, "  scale=", scale)

  for (lulc in ca_lulc$label) {

    out_path <- ca_extract_path(source, tk$layer, tk$year, lulc)
    if (file.exists(out_path) && !overwrite) {
      ca_log("  ", lulc, "  already present, skipped")
      next
    }

    pts <- ca_grid_points(lulc)

    v <- ca_extract_points(hit$file, pts$xy, band = tk$band,
                           nodata = nodata, scale = scale)

    n_na <- sum(is.na(v))
    rng <- if (n_na == length(v)) c(NA, NA) else range(v, na.rm = TRUE)

    ca_log("  ", lulc,
           "  n=", format(length(v), big.mark = ","),
           "  NA=", format(n_na, big.mark = ","),
           " (", round(100 * n_na / length(v), 2), "%)",
           "  range=", paste(signif(rng, 5), collapse = " to "),
           "  crs_transform=", isTRUE(attr(v, "transformed")))

    if (n_na == length(v)) {
      stop("All values NA for ", tk$layer, " ", tk$year, " ", lulc,
           ". Check nodata and grid alignment before continuing.")
    }

    ca_write_extract(pts$pixel_id, v, source, tk$layer, tk$year, lulc)
    rm(pts, v)
    gc(verbose = FALSE)
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# STAGE 3D. STATIC COVARIATES AND ZONE JOINS
# ---------------------------------------------------------------------------
# Stage 3d output has no year dimension, so it cannot use ca_extract_path,
# ca_write_extract, ca_task_table, ca_probe_layers, or ca_run_extract_task,
# every one of which is keyed on year. The functions below mirror those five
# for one row per pixel, and share ca_grid_points and ca_extract_points so the
# CRS rule and the chunking behaviour cannot drift between stages.
#
# Zone values are character codes and are written as strings. Writing an EPA
# Level III code or a HUC8 as float32 would silently drop leading zeros, which
# is how an exact-matching stratum turns into a wrong one.

# Resolve the file or files backing one static layer. A layer is either a
# single named file or a pattern that matches several, in which case reduce
# states how they combine.
ca_static_files <- function(reg_row) {

  dir <- ca_input(reg_row$subdir)
  if (!dir.exists(dir)) stop("Input directory missing: ", dir)

  # A resolution failure is almost always a filename that differs from the one
  # recorded in the registry, so the actual directory contents are reported in
  # the error rather than leaving the user to go and look.
  present <- function() {
    f <- list.files(dir)
    if (!length(f)) "directory is empty" else paste(f, collapse = ", ")
  }

  if (!is.na(reg_row$file)) {
    f <- file.path(dir, reg_row$file)
    if (!file.exists(f)) {
      stop("File not found: ", f, "\n  present in ", dir, ": ", present())
    }
    return(f)
  }

  f <- list.files(dir, pattern = reg_row$file_pattern, full.names = TRUE)
  if (!length(f)) {
    stop("No file matched ", reg_row$file_pattern, " in ", dir,
         "\n  present: ", present())
  }
  sort(f)
}

ca_write_static <- function(pixel_id, value, layer, lulc) {

  path <- ca_static_path(layer, lulc)

  val <- if (is.character(value)) {
    arrow::Array$create(value, type = arrow::utf8())
  } else {
    arrow::Array$create(as.numeric(value), type = arrow::float32())
  }

  tbl <- arrow::arrow_table(
    pixel_id = arrow::Array$create(as.integer(pixel_id), type = arrow::int32()),
    value    = val
  )
  arrow::write_parquet(tbl, path, compression = ca_io$compression)
  path
}

# One task per layer and LULC class. Raster and vector layers share one array,
# because both write one file per layer per class and both are driven by the
# same skip-if-present rule.
ca_static_task_table <- function(write = TRUE) {

  reg <- rbind(
    ca_static_layers[ca_static_layers$active %in% TRUE,
                     c("layer", "type"), drop = FALSE],
    ca_zone_layers[ca_zone_layers$active %in% TRUE,
                   c("layer", "type"), drop = FALSE]
  )
  if (!nrow(reg)) stop("No active layers in the stage 3d registries")

  tasks <- expand.grid(lulc = ca_lulc$label, layer = reg$layer,
                       stringsAsFactors = FALSE)
  tasks$type <- reg$type[match(tasks$layer, reg$layer)]
  tasks <- tasks[order(tasks$layer, tasks$lulc), c("layer", "type", "lulc")]
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL

  if (write) {
    dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
    write.csv(tasks, ca_meta("task_table_static.csv"), row.names = FALSE)
  }

  tasks
}

# Inventory without extracting. Reports the resolved files, CRS match, grid
# alignment, native resolution, declared NA flag, and the value range read from
# a 100,000-cell sample. The sampled range is the check that matters, because
# several of these products ship no NA flag at all and an undeclared sentinel
# shows up only as an implausible minimum.
ca_probe_static <- function(sample_n = 1e5L) {

  ref <- if (file.exists(ca_grid$meta_rds)) readRDS(ca_grid$meta_rds) else NULL
  ref_crs <- ca_ref_crs()

  rows <- lapply(seq_len(nrow(ca_static_layers)), function(i) {

    reg <- ca_static_layers[i, ]
    blank <- function(msg) {
      data.frame(layer = reg$layer, status = msg, subdir = reg$subdir,
                 stringsAsFactors = FALSE)
    }

    f <- try(ca_static_files(reg), silent = TRUE)
    if (inherits(f, "try-error")) {
      return(blank(sub("\n.*$", "", conditionMessage(attr(f, "condition")))))
    }

    r <- terra::rast(f[1])
    so <- terra::scoff(r)
    same_crs <- terra::same.crs(ref_crs, terra::crs(r))

    aligned <- NA
    if (!is.null(ref) && same_crs) {
      res_ok <- isTRUE(all.equal(terra::res(r), ref$resolution,
                                 tolerance = 1e-6))
      e <- as.vector(terra::ext(r))
      offx <- ((e[1] - ref$extent[1]) / terra::res(r)[1]) %% 1
      offy <- ((e[3] - ref$extent[3]) / terra::res(r)[2]) %% 1
      aligned <- res_ok &&
        min(offx, 1 - offx) < 1e-3 && min(offy, 1 - offy) < 1e-3
    }

    s <- terra::spatSample(r, size = sample_n, method = "regular",
                           na.rm = TRUE, warn = FALSE)[[1]]

    data.frame(
      layer        = reg$layer,
      status       = "ok",
      subdir       = reg$subdir,
      n_files      = length(f),
      files        = paste(basename(f), collapse = " "),
      reduce       = reg$reduce,
      method       = reg$method,
      n_bands      = terra::nlyr(r),
      datatype     = terra::datatype(r)[1],
      na_flag      = paste(terra::NAflag(r), collapse = " "),
      scale_file   = if (is.null(so)) NA_real_ else so[1, 1],
      offset_file  = if (is.null(so)) NA_real_ else so[1, 2],
      res_x        = terra::res(r)[1],
      same_crs     = same_crs,
      grid_aligned = aligned,
      cfg_nodata   = reg$nodata,
      smp_min      = if (length(s)) min(s, na.rm = TRUE) else NA_real_,
      smp_median   = if (length(s)) stats::median(s, na.rm = TRUE) else NA_real_,
      smp_max      = if (length(s)) max(s, na.rm = TRUE) else NA_real_,
      units        = reg$units,
      crs_proj     = terra::crs(r, proj = TRUE),
      stringsAsFactors = FALSE
    )
  })

  cols <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
  # scoff() returns a named matrix, which rbind otherwise propagates into the
  # row names as scale, scale1, scale2.
  rownames(out) <- NULL

  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  path <- ca_meta(paste0("probe_3d_static_", format(Sys.Date(), "%Y%m%d"),
                         ".csv"))
  write.csv(out, path, row.names = FALSE)
  ca_log("Static probe written to ", path)
  # grid_aligned is expected FALSE for every stage 3d layer. All seven are
  # coarser than 30 m and none sits on the NLCD origin, which is why they are
  # sampled at the point rather than joined cell to cell. The column is kept
  # only so the same probe can be read against a future 30 m covariate.
  out
}

# Field discovery for the two zone layers. Confirms that the field named in
# ca_zone_layers exists and reports its distinct values, since a shapefile
# caps attribute names at ten characters and casing differs between vintages.
ca_probe_zones <- function() {

  ref_crs <- ca_ref_crs()

  rows <- lapply(seq_len(nrow(ca_zone_layers)), function(i) {

    reg <- ca_zone_layers[i, ]
    p <- ca_find_input(reg$subdir, reg$file)

    if (is.na(p)) {
      return(data.frame(layer = reg$layer, status = "file missing",
                        stringsAsFactors = FALSE))
    }

    v <- try(terra::vect(p), silent = TRUE)
    if (inherits(v, "try-error")) {
      return(data.frame(layer = reg$layer, status = "unreadable",
                        stringsAsFactors = FALSE))
    }

    nm <- names(v)
    hit <- nm[tolower(nm) == tolower(reg$field)]

    val <- if (length(hit)) as.character(v[[hit[1]]][, 1]) else character(0)

    data.frame(
      layer        = reg$layer,
      status       = if (length(hit)) "ok" else "field not found",
      file         = basename(p),
      location     = attr(p, "location"),
      n_features   = nrow(v),
      geom_type    = terra::geomtype(v),
      same_crs     = terra::same.crs(ref_crs, terra::crs(v)),
      cfg_field    = reg$field,
      actual_field = if (length(hit)) hit[1] else NA_character_,
      n_distinct   = length(unique(val)),
      n_missing    = sum(is.na(val) | !nzchar(val)),
      max_nchar    = if (length(val)) max(nchar(val), na.rm = TRUE) else NA_integer_,
      example      = paste(utils::head(sort(unique(val)), 5), collapse = " "),
      all_fields   = paste(nm, collapse = " "),
      stringsAsFactors = FALSE
    )
  })

  cols <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))

  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  path <- ca_meta(paste0("probe_3d_zones_", format(Sys.Date(), "%Y%m%d"),
                         ".csv"))
  write.csv(out, path, row.names = FALSE)
  ca_log("Zone probe written to ", path)
  out
}

# Point in polygon join. The polygon layer is transformed into the grid CRS
# rather than the points being transformed into the polygon CRS. This is not a
# departure from the raster rule. Reprojecting a vector layer resamples
# nothing and loses nothing, and transforming a few thousand polygons is far
# cheaper than transforming 95.6 million points.
#
# Chunked because each chunk builds a SpatVector, which raster extraction does
# not. Values are returned as character.
ca_join_zone <- function(poly_path, field, xy, derive = NA_character_,
                         chunk = ca_zone_chunk) {

  ref <- ca_ref_crs()
  v <- terra::vect(poly_path)

  if (!terra::same.crs(ref, terra::crs(v))) {
    v <- terra::project(v, ref)
    transformed <- TRUE
  } else {
    transformed <- FALSE
  }

  nm <- names(v)
  hit <- nm[tolower(nm) == tolower(field)]
  if (!length(hit)) {
    stop("Field ", field, " not present in ", basename(poly_path),
         ". Fields are: ", paste(nm, collapse = " "))
  }
  v <- v[, hit[1]]

  n <- nrow(xy)
  out <- rep(NA_character_, n)

  for (s in seq(1L, n, by = chunk)) {
    e <- min(s + chunk - 1L, n)
    pts <- terra::vect(xy[s:e, , drop = FALSE], type = "points", crs = ref)
    j <- terra::extract(v, pts)

    # Results are placed by the returned point index rather than by row order.
    # A point falling in no polygon may be returned as an NA row or omitted
    # entirely depending on the terra version, and assuming order would shift
    # every subsequent value in the chunk if it were omitted.
    idcol <- if ("id.y" %in% names(j)) "id.y" else names(j)[1]
    j <- j[!duplicated(j[[idcol]]), ]

    val <- rep(NA_character_, e - s + 1L)
    val[j[[idcol]]] <- as.character(j[[hit[1]]])
    out[s:e] <- val

    rm(pts, j, val)
  }

  if (!is.na(derive)) {
    out <- switch(derive,
      substr_8 = substr(out, 1, 8),
      # formatC(flag = "0") zero-pads numeric input but left-pads character
      # input with spaces, which produced " 6" rather than "06". The codes are
      # integers held as text, so they are converted before padding.
      pad_2    = ifelse(is.na(out), NA_character_,
                        sprintf("%02d", suppressWarnings(as.integer(out)))),
      stop("Unknown zone derivation rule: ", derive)
    )
    if (identical(derive, "pad_2") && any(grepl("NA", out, fixed = TRUE))) {
      stop("pad_2 received a non-integer code. Check the source field.")
    }
  }

  out[!is.na(out) & !nzchar(out)] <- NA_character_
  attr(out, "transformed") <- transformed
  out
}

# Runs one stage 3d array task, meaning one layer and one LULC class.
ca_run_static_task <- function(task_id, overwrite = FALSE) {

  tasks <- ca_static_task_table(write = FALSE)
  if (is.na(task_id) || task_id < 1 || task_id > nrow(tasks)) {
    stop("task_id ", task_id, " outside 1 to ", nrow(tasks))
  }

  tk <- tasks[task_id, ]
  out_path <- ca_static_path(tk$layer, tk$lulc)

  if (file.exists(out_path) && !overwrite) {
    ca_log("Task ", task_id, "  ", tk$layer, " ", tk$lulc,
           "  already present, skipped")
    return(invisible(TRUE))
  }

  pts <- ca_grid_points(tk$lulc)

  if (tk$type == "raster") {

    reg <- ca_static_layers[ca_static_layers$layer == tk$layer, ][1, ]
    files <- ca_static_files(reg)

    ca_log("Task ", task_id, "  ", tk$layer, " ", tk$lulc,
           "  files=", paste(basename(files), collapse = " "),
           "  method=", reg$method, "  nodata=", reg$nodata)

    vals <- lapply(files, function(f) {
      ca_extract_points(f, pts$xy, band = reg$band, nodata = reg$nodata,
                        scale = reg$scale, method = reg$method)
    })

    transformed <- isTRUE(attr(vals[[1]], "transformed"))

    v <- if (length(vals) == 1L) {
      vals[[1]]
    } else if (identical(reg$reduce, "mean")) {
      rowMeans(do.call(cbind, vals), na.rm = TRUE)
    } else {
      stop("Layer ", tk$layer, " resolves to ", length(vals),
           " files but reduce is ", reg$reduce)
    }
    # rowMeans returns NaN, not NA, where every input is NA.
    v[is.nan(v)] <- NA_real_

  } else {

    reg <- ca_zone_layers[ca_zone_layers$layer == tk$layer, ][1, ]
    p <- ca_find_input(reg$subdir, reg$file, require_scratch = TRUE)

    ca_log("Task ", task_id, "  ", tk$layer, " ", tk$lulc,
           "  file=", basename(p), "  field=", reg$field,
           "  derive=", reg$derive)

    v <- ca_join_zone(p, reg$field, pts$xy, derive = reg$derive)
    transformed <- isTRUE(attr(v, "transformed"))
  }

  n_na <- sum(is.na(v))

  if (is.character(v)) {
    ca_log("  n=", format(length(v), big.mark = ","),
           "  NA=", format(n_na, big.mark = ","),
           " (", round(100 * n_na / length(v), 2), "%)",
           "  distinct=", length(unique(v[!is.na(v)])),
           "  crs_transform=", transformed)
  } else {
    rng <- if (n_na == length(v)) c(NA, NA) else range(v, na.rm = TRUE)
    ca_log("  n=", format(length(v), big.mark = ","),
           "  NA=", format(n_na, big.mark = ","),
           " (", round(100 * n_na / length(v), 2), "%)",
           "  range=", paste(signif(rng, 5), collapse = " to "),
           "  crs_transform=", transformed)
  }

  if (n_na == length(v)) {
    stop("All values NA for ", tk$layer, " ", tk$lulc,
         ". Check CRS, extent, and nodata before continuing.")
  }

  ca_write_static(pts$pixel_id, v, tk$layer, tk$lulc)

  rm(pts, v)
  gc(verbose = FALSE)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# STAGE 3C. MANAGEMENT AND OWNERSHIP POLYGONS
# ---------------------------------------------------------------------------
# Stage 3c differs from 3d in one way that drives the whole design. A pixel can
# fall inside several polygons, and every one of them matters. A pixel with
# three recorded thinning entries is excluded from the treatment arm under
# ca_thin_rules$single_event_only, and a first-hit join cannot express that.
# Output is therefore long, one row per pixel per intersecting polygon.
#
# Per ca_config section 0.6, every archival record is retained through
# extraction. Nothing is dropped for being non_disturbing or Variable. The only
# filter applied here is the archival status flag, which distinguishes a record
# of work done from a record of work planned. That is not an analytical choice.

ca_vector_cols <- c("pixel_id", "poly_uid", "event_year", "key_value",
                    "key_value_2", "method_value", "intensity",
                    "record_class", "xwalk_origin")

# Knight crosswalk plus the locally resolved aliases, as one lookup keyed on
# table and normalised activity string.
ca_intensity_table <- function() {

  kn <- utils::read.csv(ca_knight_csv, stringsAsFactors = FALSE)
  al <- utils::read.csv(ca_alias_csv, stringsAsFactors = FALSE)

  kn <- data.frame(table = kn$table, key = kn$activity_key,
                   intensity = kn$intensity, method = kn$method,
                   record_class = kn$record_class, origin = "knight",
                   stringsAsFactors = FALSE)

  if (nrow(al)) {
    # The missing-value token must not be normalised. ca_norm_activity strips
    # the angle brackets and yields MISSING, while the layer side keeps
    # <MISSING>, and that mismatch left four THP polygons unresolved.
    is_miss <- is.na(al$observed_activity) |
      !nzchar(trimws(al$observed_activity)) |
      al$observed_activity == "<MISSING>"
    key <- ifelse(is_miss, "<MISSING>", ca_norm_activity(al$observed_activity))
    al <- data.frame(table = al$table, key = key,
                     intensity = al$intensity, method = al$method,
                     record_class = al$record_class,
                     origin = paste0("alias_", al$resolution),
                     stringsAsFactors = FALSE)
  } else {
    al <- kn[0, ]
  }

  out <- rbind(kn, al)
  out[!duplicated(out[, c("table", "key")]), ]
}

# Reads one registry row and returns geometry and attributes as two separate
# objects rather than one annotated layer.
#
# The first version pushed the resolved character columns through the geometry
# object and failed on coercion between sf and terra. Attributes never need to
# travel with the geometry here. The join only has to return which polygon a
# point fell in, so the geometry carries a single integer pid and everything
# else stays in a plain data frame indexed by that pid. This removes the
# conversion entirely, keeps sf out of the pipeline, and makes the join
# cheaper, because terra is not copying eight character columns per chunk.
ca_read_vector_layer <- function(reg, xwalk = NULL) {

  paths <- ca_find_input(reg$subdir, reg$file, require_scratch = TRUE)
  if (!is.na(reg$file_extra)) {
    paths <- c(paths,
               ca_find_input(reg$subdir, reg$file_extra, require_scratch = TRUE))
  }

  parts <- lapply(paths, function(p) terra::vect(as.character(p)))

  # The two CAL FIRE permit eras carry the same schema on every field read
  # here, so they are stacked on the intersection of their columns.
  if (length(parts) > 1L) {
    keep <- Reduce(intersect, lapply(parts, names))
    parts <- lapply(parts, function(x) x[, keep])
    v <- do.call(rbind, parts)
  } else {
    v <- parts[[1]]
  }
  rm(parts)

  # Winding-order and self-intersection problems are corrected here rather
  # than surfacing later during the join, where the failure is harder to trace.
  bad <- suppressWarnings(!terra::is.valid(v))
  n_invalid <- sum(bad, na.rm = TRUE)
  if (n_invalid) v <- terra::makeValid(v)

  att <- as.data.frame(v)
  n_read <- nrow(att)

  fld <- function(f) {
    if (is.na(f)) return(NULL)
    hit <- names(att)[tolower(names(att)) == tolower(f)]
    if (!length(hit)) stop("Field ", f, " not in ", reg$layer,
                           ". Fields are: ", paste(names(att), collapse = " "))
    att[[hit[1]]]
  }

  # 1. Archival status. A record of work planned is not a record of work done.
  n_status <- NA_integer_
  if (!is.na(reg$status_field)) {
    st <- as.character(fld(reg$status_field))
    keep <- !is.na(st) & trimws(st) == reg$status_keep
    n_status <- sum(!keep)
    v <- v[keep, ]
    att <- att[keep, , drop = FALSE]
  }
  if (!nrow(att)) stop("No records survived the status filter for ", reg$layer)

  # 2. Event year. The completion date is preferred over the cohort year,
  # because the cohort year records when a plan was filed and the completion
  # date records when the ground was worked. date_null values are placeholders
  # and must never be read as a year.
  yr <- rep(NA_integer_, nrow(att))

  if (!is.na(reg$date_field)) {
    d <- fld(reg$date_field)
    dc <- trimws(as.character(d))
    y <- suppressWarnings(as.integer(format(as.Date(d), "%Y")))
    null_tok <- reg$date_null
    if (!is.na(null_tok)) {
      null_yr <- suppressWarnings(as.integer(null_tok))
      is_null <- dc == null_tok
      if (!is.na(null_yr)) is_null <- is_null | (!is.na(y) & y == null_yr)
      y[is_null %in% TRUE] <- NA_integer_
    }
    y[!is.na(y) & (y < 1900L | y > 2100L)] <- NA_integer_
    yr <- y
  }

  if (!is.na(reg$year_field)) {
    y2 <- suppressWarnings(as.integer(as.character(fld(reg$year_field))))
    y2[!is.na(y2) & (y2 < 1800L | y2 > 2100L)] <- NA_integer_
    yr[is.na(yr)] <- y2[is.na(yr)]
  }

  # 3. Offsets carry no date field at all. The cohort year comes from the CARB
  # issuance lookup, joined on ARB_id and restricted to compliance projects.
  if (identical(reg$layer, "offset")) {
    ly <- utils::read.csv(ca_offset_years_csv, stringsAsFactors = FALSE)
    jf <- names(ly)[tolower(names(ly)) == tolower(ca_offset_join_field)][1]
    ly <- ly[toupper(trimws(ly$ea_cop)) == toupper(ca_offset_program), ]
    k <- trimws(as.character(fld(ca_offset_join_field)))
    yr <- suppressWarnings(
      as.integer(ly[[ca_offset_year_field]][match(k, trimws(ly[[jf]]))]))
  }

  chr <- function(f) {
    x <- fld(f)
    if (is.null(x)) rep(NA_character_, nrow(att)) else as.character(x)
  }

  # A polygon identifier for duplicate detection across stacked files.
  uid_candidates <- c("HD_NUM", "SUID", "ARB_id", "GEOID", "SUID_NMA",
                      "OBJECTID", "GlobalID")
  uid_hit <- uid_candidates[tolower(uid_candidates) %in% tolower(names(att))]
  poly_uid <- if (length(uid_hit)) chr(uid_hit[1]) else
    as.character(seq_len(nrow(att)))

  key  <- chr(reg$key_field)
  key2 <- chr(reg$key_field_2)
  meth <- chr(reg$method_field)

  # 4. Intensity. Only the thinning layers carry a crosswalk.
  intensity <- rep(NA_character_, nrow(att))
  rclass <- rep(NA_character_, nrow(att))
  origin <- rep(NA_character_, nrow(att))

  if (!is.na(reg$crosswalk)) {
    if (is.null(xwalk)) xwalk <- ca_intensity_table()
    tb <- xwalk[xwalk$table == reg$crosswalk, ]
    k <- ifelse(is.na(key) | !nzchar(trimws(key)), "<MISSING>",
                ca_norm_activity(key))
    i <- match(k, tb$key)
    intensity <- tb$intensity[i]
    rclass <- tb$record_class[i]
    origin <- tb$origin[i]
    n_unres <- sum(is.na(i))
    if (n_unres) {
      unres <- sort(unique(key[is.na(i)]))
      ca_log("  ", reg$layer, "  ", n_unres, " records unresolved by the ",
             reg$crosswalk, " crosswalk: ",
             paste(utils::head(unres, 10), collapse = " | "))
    }
  }

  # Geometry keeps one integer column and nothing else.
  terra::values(v) <- data.frame(pid = seq_len(nrow(att)))

  ref <- ca_ref_crs()
  transformed <- !terra::same.crs(ref, terra::crs(v))
  if (transformed) v <- terra::project(v, ref)

  out <- list(
    geom = v,
    attr = data.frame(
      pid          = seq_len(nrow(att)),
      poly_uid     = poly_uid,
      event_year   = as.integer(yr),
      key_value    = key,
      key_value_2  = key2,
      method_value = meth,
      intensity    = intensity,
      record_class = rclass,
      xwalk_origin = origin,
      stringsAsFactors = FALSE),
    n_read = n_read,
    n_invalid = n_invalid,
    n_status_dropped = n_status,
    n_dup_uid = sum(duplicated(poly_uid)),
    n_year_na = sum(is.na(yr)),
    transformed = transformed
  )
  out
}

# Point in polygon returning every intersection, not the first. Chunked over
# points, with the per-chunk result rbound. Empty chunks are skipped so a layer
# covering a small part of the state costs almost nothing.
ca_join_vector_all <- function(lay, pts, chunk = ca_vector_chunk) {

  ref <- ca_ref_crs()
  sv <- lay$geom
  att <- lay$attr
  xy <- pts$xy
  n <- nrow(xy)
  acc <- vector("list", ceiling(n / chunk))
  k <- 0L

  for (s in seq(1L, n, by = chunk)) {
    e <- min(s + chunk - 1L, n)
    p <- terra::vect(xy[s:e, , drop = FALSE], type = "points", crs = ref)
    j <- terra::extract(sv, p)
    idcol <- if ("id.y" %in% names(j)) "id.y" else names(j)[1]
    j <- j[!is.na(j$pid), , drop = FALSE]
    if (nrow(j)) {
      k <- k + 1L
      m <- match(j$pid, att$pid)
      acc[[k]] <- data.frame(
        pixel_id     = pts$pixel_id[s:e][j[[idcol]]],
        poly_uid     = att$poly_uid[m],
        event_year   = att$event_year[m],
        key_value    = att$key_value[m],
        key_value_2  = att$key_value_2[m],
        method_value = att$method_value[m],
        intensity    = att$intensity[m],
        record_class = att$record_class[m],
        # Carried through to the file rather than left in the probe. The
        # response letter needs the treated pixel share attributable to
        # locally assigned activities, and that cannot be recovered later
        # without rerunning the whole stage.
        xwalk_origin = att$xwalk_origin[m],
        stringsAsFactors = FALSE)
    }
    rm(p, j)
  }

  if (!k) {
    out <- data.frame(matrix(nrow = 0, ncol = length(ca_vector_cols),
                             dimnames = list(NULL, ca_vector_cols)))
    return(out)
  }
  do.call(rbind, acc[seq_len(k)])
}

ca_write_vector <- function(df, layer, lulc) {
  path <- ca_vector_path(layer, lulc)
  tbl <- arrow::arrow_table(
    pixel_id     = arrow::Array$create(as.integer(df$pixel_id),
                                       type = arrow::int32()),
    poly_uid     = arrow::Array$create(as.character(df$poly_uid),
                                       type = arrow::utf8()),
    event_year   = arrow::Array$create(as.integer(df$event_year),
                                       type = arrow::int32()),
    key_value    = arrow::Array$create(as.character(df$key_value),
                                       type = arrow::utf8()),
    key_value_2  = arrow::Array$create(as.character(df$key_value_2),
                                       type = arrow::utf8()),
    method_value = arrow::Array$create(as.character(df$method_value),
                                       type = arrow::utf8()),
    intensity    = arrow::Array$create(as.character(df$intensity),
                                       type = arrow::utf8()),
    record_class = arrow::Array$create(as.character(df$record_class),
                                       type = arrow::utf8()),
    xwalk_origin = arrow::Array$create(as.character(df$xwalk_origin),
                                       type = arrow::utf8())
  )
  arrow::write_parquet(tbl, path, compression = ca_io$compression)
  path
}

ca_vector_task_table <- function(write = TRUE) {
  reg <- ca_vector_layers[ca_vector_layers$active %in% TRUE, ]
  tasks <- expand.grid(lulc = ca_lulc$label, layer = reg$layer,
                       stringsAsFactors = FALSE)
  tasks <- tasks[order(tasks$layer, tasks$lulc), c("layer", "lulc")]
  tasks$task_id <- seq_len(nrow(tasks))
  rownames(tasks) <- NULL
  if (write) {
    dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
    write.csv(tasks, ca_meta("task_table_vector.csv"), row.names = FALSE)
  }
  tasks
}

# Inventory without joining. Reads and resolves every layer, then estimates the
# join cardinality from a sample of points. Cardinality is the number that
# decides whether the array is feasible, because output is long and a pixel
# inside many overlapping records produces many rows.
ca_probe_vector_layers <- function(sample_n = 2e5L, lulc = ca_lulc$label[1]) {

  xwalk <- ca_intensity_table()
  pts <- ca_grid_points(lulc)
  n_total <- length(pts$pixel_id)
  # Fixed seed. The sample drives numbers that are quoted rather than only
  # eyeballed, and three runs of this probe returned 1,763, 1,830, and 1,813
  # rows for the same layer purely from resampling.
  set.seed(20260731L)
  idx <- sort(sample.int(n_total, min(sample_n, n_total)))
  smp <- list(pixel_id = pts$pixel_id[idx],
              xy = pts$xy[idx, , drop = FALSE])
  rm(pts); gc(verbose = FALSE)

  # min and max over an all-NA year column return Inf and -Inf, which read as
  # data rather than as absence. Layers with no year field are the expected
  # case, so they are reported as NA.
  fin <- function(x) if (is.finite(x)) x else NA_integer_

  rows <- lapply(seq_len(nrow(ca_vector_layers)), function(i) {

    reg <- ca_vector_layers[i, ]
    lay <- try(ca_read_vector_layer(reg, xwalk), silent = TRUE)
    if (inherits(lay, "try-error")) {
      return(data.frame(layer = reg$layer, status = sub("\n.*$", "",
                        conditionMessage(attr(lay, "condition"))),
                        stringsAsFactors = FALSE))
    }

    j <- ca_join_vector_all(lay, smp)
    per_px <- if (nrow(j)) table(j$pixel_id) else integer(0)
    a <- lay$attr

    data.frame(
      layer            = reg$layer,
      status           = "ok",
      role             = reg$role,
      n_read           = lay$n_read,
      n_kept           = nrow(a),
      n_invalid_fixed  = lay$n_invalid,
      n_status_dropped = lay$n_status_dropped,
      n_dup_uid        = lay$n_dup_uid,
      pct_year_na      = round(100 * lay$n_year_na / nrow(a), 2),
      yr_min           = fin(suppressWarnings(min(a$event_year, na.rm = TRUE))),
      yr_max           = fin(suppressWarnings(max(a$event_year, na.rm = TRUE))),
      pct_unresolved   = if (is.na(reg$crosswalk)) NA_real_ else
        round(100 * mean(is.na(a$record_class)), 2),
      # Split, because the two resolution types carry very different weight.
      # A wording variant maps an observed string onto a published Knight row
      # with no judgement involved. An assignment classifies an activity that
      # the published tables do not contain, and that is the number a reviewer
      # would want stated.
      pct_variant      = if (is.na(reg$crosswalk)) NA_real_ else
        round(100 * mean(a$xwalk_origin %in% "alias_alias"), 2),
      pct_assigned     = if (is.na(reg$crosswalk)) NA_real_ else
        round(100 * mean(a$xwalk_origin %in% "alias_assigned"), 2),
      smp_pts          = length(smp$pixel_id),
      smp_rows         = nrow(j),
      smp_pct_hit      = round(100 * length(per_px) / length(smp$pixel_id), 2),
      rows_per_hit     = if (length(per_px)) round(mean(per_px), 2) else 0,
      max_per_hit      = if (length(per_px)) max(per_px) else 0L,
      # Linear extrapolation to the full class. Approximate, and the only
      # number that says whether the array fits in memory. Computed in double
      # precision, because the integer product overflows well before the row
      # counts involved here.
      est_rows_class   = if (length(smp$pixel_id))
        round(as.numeric(nrow(j)) * as.numeric(n_total) /
                as.numeric(length(smp$pixel_id))) else 0,
      # Evergreen carries roughly 43 times the Decid points and is the class
      # that sizes the job, so it is extrapolated explicitly rather than left
      # for the reader to scale.
      est_rows_everg   = if (length(smp$pixel_id))
        round(as.numeric(nrow(j)) * as.numeric(ca_lulc_n_points("Everg")) /
                as.numeric(length(smp$pixel_id))) else 0,
      stringsAsFactors = FALSE
    )
  })

  cols <- unique(unlist(lapply(rows, names)))
  out <- do.call(rbind, lapply(rows, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
  rownames(out) <- NULL

  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  path <- ca_meta(paste0("probe_3c_join_", format(Sys.Date(), "%Y%m%d"),
                         ".csv"))
  write.csv(out, path, row.names = FALSE)
  ca_log("Vector join probe written to ", path)
  out
}

ca_run_vector_task <- function(task_id, overwrite = FALSE) {

  tasks <- ca_vector_task_table(write = FALSE)
  if (is.na(task_id) || task_id < 1 || task_id > nrow(tasks)) {
    stop("task_id ", task_id, " outside 1 to ", nrow(tasks))
  }

  tk <- tasks[task_id, ]
  out_path <- ca_vector_path(tk$layer, tk$lulc)
  if (file.exists(out_path) && !overwrite) {
    ca_log("Task ", task_id, "  ", tk$layer, " ", tk$lulc,
           "  already present, skipped")
    return(invisible(TRUE))
  }

  reg <- ca_vector_layers[ca_vector_layers$layer == tk$layer, ][1, ]
  ca_log("Task ", task_id, "  ", tk$layer, " ", tk$lulc,
         "  file=", reg$file,
         if (!is.na(reg$file_extra)) paste0(" + ", reg$file_extra) else "")

  lay <- ca_read_vector_layer(reg)
  ca_log("  polygons read=", lay$n_read,
         "  kept=", nrow(lay$attr),
         "  invalid fixed=", lay$n_invalid,
         "  status dropped=", lay$n_status_dropped,
         "  duplicate uid=", lay$n_dup_uid,
         "  year NA=", lay$n_year_na,
         "  crs_transform=", lay$transformed)

  pts <- ca_grid_points(tk$lulc)
  j <- ca_join_vector_all(lay, pts)

  n_px <- length(unique(j$pixel_id))
  ca_log("  rows=", format(nrow(j), big.mark = ","),
         "  pixels hit=", format(n_px, big.mark = ","),
         " (", round(100 * n_px / length(pts$pixel_id), 2), "%)",
         "  rows per hit=", if (n_px) round(nrow(j) / n_px, 2) else 0)

  ca_write_vector(j, tk$layer, tk$lulc)

  rm(lay, pts, j); gc(verbose = FALSE)
  invisible(TRUE)
}
