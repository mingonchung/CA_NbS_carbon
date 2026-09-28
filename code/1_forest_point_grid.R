### CA carbon revision pipeline
### Stage 1. Forest point grid
###
### Converts NLCD 2001 forest classes to a frozen point grid that defines every
### pixel identifier used downstream. Runs once. Do not rerun after stage 3
### extraction begins, because every join key depends on this output.
###
### pixel_id is the terra cell index of the NLCD reference raster. It is stable,
### sortable, and invertible back to coordinates, so no lookup table is needed
### and the grid can be reconstructed from the raster alone.
###
### Reads the raster in row blocks rather than loading it whole and streams each
### block to Parquet on scratch. Peak memory is set by block size, not by forest
### extent. Part files are consolidated into one Parquet file per class on
### completion, because CURC filesystems perform poorly on many small files.
###
### Usage
###   Rscript 1_forest_point_grid.R            all forest classes, one pass
###   Rscript 1_forest_point_grid.R Everg      one class only, for reruns

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
ca_require(c("terra", "arrow", "sf"))

suppressPackageStartupMessages({
  library(terra)
  library(arrow)
  library(sf)
})

ca_stamp("stage1")
ca_log("Threads: ", ca_threads(), "   tempdir: ", ca_tmpdir())

# ---------------------------------------------------------------------------
# 1.1  ARGUMENTS AND TARGET CLASSES
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

lulc_tbl <- ca_lulc
if (length(args) && nzchar(args[1]) && args[1] %in% ca_lulc$label) {
  lulc_tbl <- ca_lulc[ca_lulc$label == args[1], , drop = FALSE]
}

ca_log("Target classes: ", paste(lulc_tbl$label, collapse = ", "))

# ---------------------------------------------------------------------------
# 1.2  STAGE INPUT TO SCRATCH
# ---------------------------------------------------------------------------
# The master raster lives on /projects, which must not carry job I/O. Copy it
# to scratch once, then read only from scratch.

if (!file.exists(ca_grid$nlcd_raster)) {
  ca_log("Staging NLCD raster from /projects to scratch")
  ca_stage(file.path("nlcd", basename(ca_grid$nlcd_master)))
}

# ---------------------------------------------------------------------------
# 1.3  REFERENCE RASTER
# ---------------------------------------------------------------------------

nlcd <- rast(ca_grid$nlcd_raster)

if (nlyr(nlcd) != 1L) {
  stop("Expected a single-layer raster, found ", nlyr(nlcd), " layers")
}

crs_wkt <- crs(nlcd)

if (!nzchar(crs_wkt)) {
  stop("Reference raster has no CRS. Downstream extraction cannot proceed.")
}
if (is.lonlat(nlcd)) {
  stop("Reference raster is geographic. A projected CRS is required so that ",
       "the 150 m sampling distance at stage 7 is metric.")
}
if (ncell(nlcd) > .Machine$integer.max) {
  stop("ncell exceeds 32-bit integer range. Store pixel_id as double.")
}

ca_log("CRS:        ", crs(nlcd, proj = TRUE))
ca_log("Resolution: ", paste(res(nlcd), collapse = " x "))
ca_log("Dimensions: ", nrow(nlcd), " rows x ", ncol(nlcd), " cols")

# ---------------------------------------------------------------------------
# 1.4  BLOCK-WISE CONVERSION TO POINTS
# ---------------------------------------------------------------------------

parquet_dir <- ca_grid$parquet_dir
part_root <- file.path(parquet_dir, "_parts")

# Clear any previous partial run so a rerun cannot mix vintages.
for (lab in lulc_tbl$label) {
  for (d in c(file.path(part_root, lab), file.path(parquet_dir, paste0("lulc=", lab)))) {
    if (dir.exists(d)) unlink(d, recursive = TRUE)
  }
  dir.create(file.path(part_root, lab), recursive = TRUE, showWarnings = FALSE)
}

bs <- if (is.null(ca_grid$chunk_rows)) {
  blocks(nlcd, n = 4)
} else {
  starts <- seq(1L, nrow(nlcd), by = ca_grid$chunk_rows)
  list(n = length(starts), row = starts,
       nrows = pmin(ca_grid$chunk_rows, nrow(nlcd) - starts + 1L))
}

ca_log("Blocks: ", bs$n)

n_col <- ncol(nlcd)
counts <- setNames(rep(0, nrow(lulc_tbl)), lulc_tbl$label)

readStart(nlcd)
on.exit(try(readStop(nlcd), silent = TRUE), add = TRUE)

for (b in seq_len(bs$n)) {

  vals <- readValues(nlcd, row = bs$row[b], nrows = bs$nrows[b])
  offset <- (as.numeric(bs$row[b]) - 1) * n_col

  for (k in seq_len(nrow(lulc_tbl))) {

    code <- lulc_tbl$code[k]
    lab <- lulc_tbl$label[k]

    idx <- which(vals == code)
    if (!length(idx)) next

    cells <- offset + idx
    xy <- xyFromCell(nlcd, cells)

    write_parquet(
      data.frame(
        pixel_id  = as.integer(cells),
        lulc_code = rep(as.integer(code), length(cells)),
        x         = xy[, 1],
        y         = xy[, 2]
      ),
      file.path(part_root, lab, sprintf("part-%05d.parquet", b)),
      compression = ca_io$compression
    )

    counts[lab] <- counts[lab] + length(cells)
    rm(cells, xy)
  }

  rm(vals)
  if (b %% 10 == 0 || b == bs$n) {
    ca_log("Block ", b, " of ", bs$n, ". Totals: ",
           paste(names(counts), format(counts, big.mark = ","),
                 sep = "=", collapse = "  "))
    gc(verbose = FALSE)
  }
}

readStop(nlcd)

# ---------------------------------------------------------------------------
# 1.5  CONSOLIDATE PART FILES
# ---------------------------------------------------------------------------
# One file per class replaces the per-block parts. Written through an Arrow
# dataset so the consolidation never loads a full class into memory.

for (k in seq_len(nrow(lulc_tbl))) {

  lab <- lulc_tbl$label[k]
  src <- file.path(part_root, lab)
  if (!length(list.files(src, pattern = "\\.parquet$"))) next

  dst <- file.path(parquet_dir, paste0("lulc=", lab))
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)

  write_dataset(
    open_dataset(src),
    path = dst,
    format = "parquet",
    basename_template = paste0("forest_grid_", lab, "-{i}.parquet"),
    max_rows_per_file = 0,
    compression = ca_io$compression
  )

  ca_log("Consolidated ", lab)
}

if (isTRUE(ca_io$consolidate)) unlink(part_root, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1.6  GEOPACKAGE EXPORT
# ---------------------------------------------------------------------------
# One layer per class, written from the Parquet store so the two representations
# cannot diverge. Downstream stages read Parquet, not the GeoPackage, so this
# step is skippable through ca_grid$write_gpkg.

if (isTRUE(ca_grid$write_gpkg)) {

  dir.create(ca_grid$gpkg_dir, recursive = TRUE, showWarnings = FALSE)

  for (k in seq_len(nrow(lulc_tbl))) {

    lab <- lulc_tbl$label[k]
    files <- sort(list.files(file.path(parquet_dir, paste0("lulc=", lab)),
                             pattern = "\\.parquet$", full.names = TRUE))
    if (!length(files)) next

    gpkg <- file.path(ca_grid$gpkg_dir, paste0("forest_grid_", lab, ".gpkg"))
    if (file.exists(gpkg)) unlink(gpkg)

    for (f in seq_along(files)) {
      tbl <- read_parquet(files[f])
      pts <- st_as_sf(tbl, coords = c("x", "y"), crs = crs_wkt, remove = FALSE)
      st_write(pts, gpkg, layer = paste0("forest_", lab),
               append = f > 1L, quiet = TRUE)
      rm(tbl, pts)
      gc(verbose = FALSE)
    }

    ca_log("GeoPackage written: ", basename(gpkg))
  }
}

# ---------------------------------------------------------------------------
# 1.7  FREEZE THE GRID REFERENCE
# ---------------------------------------------------------------------------
# Written to /projects, which is backed up. Every later script resolves the CRS
# and grid geometry through ca_ref_crs() and this object rather than by
# reopening the NLCD raster.

grid_ref <- list(
  source_raster = ca_grid$nlcd_master,
  crs_wkt       = crs_wkt,
  crs_proj      = crs(nlcd, proj = TRUE),
  extent        = as.vector(ext(nlcd)),
  resolution    = res(nlcd),
  dim           = c(nrow = nrow(nlcd), ncol = ncol(nlcd)),
  ncell         = ncell(nlcd),
  lulc          = lulc_tbl,
  counts        = counts,
  parquet_dir   = parquet_dir,
  built_on      = Sys.time(),
  r_version     = R.version.string,
  terra_version = as.character(packageVersion("terra"))
)

dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
saveRDS(grid_ref, ca_grid$meta_rds)

summary_tbl <- data.frame(
  lulc     = names(counts),
  code     = lulc_tbl$code[match(names(counts), lulc_tbl$label)],
  n_pixels = as.numeric(counts),
  area_km2 = as.numeric(counts) * prod(res(nlcd)) / 1e6,
  stringsAsFactors = FALSE
)

write.csv(summary_tbl, ca_grid$summary_csv, row.names = FALSE)

# The grid is expensive to rebuild and scratch purges after 90 days, so the
# consolidated Parquet store is copied to the backed-up archive.
ca_archive(parquet_dir, subdir = "work/grid")

ca_log("Grid frozen. Reference written to ", ca_grid$meta_rds)
print(summary_tbl)

# ---------------------------------------------------------------------------
# 1.8  VERIFICATION
# ---------------------------------------------------------------------------
# Confirms that pixel_id is unique, that reconstructed coordinates match the
# raster, and that the recorded class matches the raster value. Run on a sample
# so the check stays cheap.

verify_grid <- function(n_sample = 5e4) {

  report <- list()

  for (lab in lulc_tbl$label) {

    part <- file.path(parquet_dir, paste0("lulc=", lab))
    if (!dir.exists(part)) next

    ds <- open_dataset(part)
    smp <- as.data.frame(head(ds, n_sample))
    xy <- xyFromCell(nlcd, smp$pixel_id)

    report[[lab]] <- data.frame(
      lulc         = lab,
      n_rows       = ds$num_rows,
      n_checked    = nrow(smp),
      dup_id       = anyDuplicated(smp$pixel_id) > 0,
      coords_match = max(abs(xy[, 1] - smp$x)) < 1e-6 &&
                     max(abs(xy[, 2] - smp$y)) < 1e-6,
      class_match  = all(nlcd[smp$pixel_id][[1]] == smp$lulc_code),
      stringsAsFactors = FALSE
    )
  }

  report <- do.call(rbind, report)
  print(report)

  if (any(report$dup_id) || !all(report$coords_match) ||
      !all(report$class_match)) {
    stop("Grid verification failed.")
  }

  # Row counts written to Parquet must equal the counts tallied during the read
  # pass. A mismatch means a part file was lost or written twice.
  if (!isTRUE(all.equal(sort(report$n_rows),
                        sort(as.numeric(counts[report$lulc]))))) {
    stop("Parquet row counts disagree with the block-pass tally.")
  }

  ca_log("Total rows: ", format(sum(report$n_rows), big.mark = ","))
  invisible(report)
}

verify_grid()

ca_log("Stage 1 complete.")
