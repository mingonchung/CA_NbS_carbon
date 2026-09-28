### CA carbon revision pipeline
### Stage 3d inputs. Copy masters from /projects to /scratch/alpine.
###
### Same pattern and same reasons as 3c_stage_inputs.R. The 3d layers are small
### enough to keep a backed-up master on /projects, but CURC prohibits job I/O
### there and the 3d extraction runs as an array.
###
### Twelve files across six directories.
###   dem/        SRTM_90m_elevation_CA_prj.tif, _slope_, _aspect_
###   prism/      prism_ppt_us_30s_2020_avg_30y.tif, prism_tmean_...
###   worldpop/   CA_pd_2000_1km_UNadj_prj.tif, 2001, 2002
###   travel/     CA_travel_time_to_cities_4_prj.tif, 5, 11
###   ecoregion/  ca_eco_l3_prj.shp, as a file set
###   huc/        WBDHU12_CA_proj.shp, as a file set
###
### Travel time classes 4 and 5 are staged although only class 11 is active, so
### that a sensitivity run needs no second transfer. The boundary layer is not
### staged, since stage 3d does not use it. The grid is already California only.
###
### Run once, on a compile node or as a short job. Safe to rerun, since existing
### files are skipped unless OVERWRITE is TRUE.
###
### Usage
###   Rscript 3d_stage_inputs.R

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))

ca_stamp("stage3d_staging")

OVERWRITE <- FALSE

# ---------------------------------------------------------------------------
# RASTER COVARIATES
# ---------------------------------------------------------------------------
# Copied by directory rather than by registry row, because the registry names
# only the active files and the inactive travel time classes should travel with
# them.

raster_subdirs <- unique(ca_static_layers$subdir)

for (sub in raster_subdirs) {

  src <- ca_master(sub)
  dst <- ca_input(sub)

  if (!dir.exists(src)) {
    ca_log(sub, "  master directory not found: ", src)
    next
  }

  dir.create(dst, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(src, pattern = "\\.(tif|tiff|vrt)$",
                      full.names = TRUE, recursive = TRUE,
                      ignore.case = TRUE)

  n <- 0L
  for (f in files) {
    d <- file.path(dst, basename(f))
    if (file.exists(d) && !OVERWRITE) next
    if (file.copy(f, d, overwrite = TRUE)) n <- n + 1L
  }

  # Sidecars GDAL may rely on for GeoTIFFs.
  for (ext in c("tfw", "tif.aux.xml", "tif.ovr", "tif.vat.dbf")) {
    pat <- paste0("\\.", gsub("\\.", "\\\\.", ext), "$")
    aux <- list.files(src, pattern = pat, full.names = TRUE, recursive = TRUE)
    if (length(aux)) {
      file.copy(aux, file.path(dst, basename(aux)), overwrite = OVERWRITE)
    }
  }

  ca_log(sprintf("%-10s %d masters, %d copied, %d now on scratch",
                 sub, length(files), n,
                 length(list.files(dst, pattern = "\\.(tif|tiff|vrt)$",
                                   ignore.case = TRUE))))
}

# ---------------------------------------------------------------------------
# ZONE LAYERS
# ---------------------------------------------------------------------------
# Shapefiles are a file set. Copying only the .shp yields a layer that opens
# with no attributes or no CRS, which fails late and confusingly.

for (i in seq_len(nrow(ca_zone_layers))) {

  lay <- ca_zone_layers$layer[i]
  sub <- ca_zone_layers$subdir[i]
  fil <- ca_zone_layers$file[i]

  if (!file.exists(ca_master(sub, fil))) {
    ca_log(lay, "  master not found: ", ca_master(sub, fil))
    next
  }

  copied <- ca_stage_shapefile(sub, fil, overwrite = OVERWRITE)
  exts <- toupper(tools::file_ext(basename(copied)))

  ca_log(sprintf("%-14s %d files  %s", lay, length(copied),
                 paste(sort(unique(exts)), collapse = " ")))

  if (!"PRJ" %in% exts) {
    warning(lay, " has no .prj sidecar. CRS will be undefined on read.",
            call. = FALSE)
  }
  if (!"DBF" %in% exts) {
    warning(lay, " has no .dbf sidecar. Attributes will be missing.",
            call. = FALSE)
  }
}

# ---------------------------------------------------------------------------
# VERIFY
# ---------------------------------------------------------------------------
# Resolution goes through ca_find_input, which prefers scratch. Anything
# reported as archive would be read from /projects by an array task, which CURC
# forbids, so those must be resolved before submitting.

ca_log("--- verification, resolved from scratch first ---")

# A MISSING line here almost always means the filename on disk differs from the
# one recorded in ca_static_layers, not that the copy failed, since the copy
# counts above are reported per directory. The directory listing is printed for
# any layer that does not resolve so the mismatch is visible immediately.

miss <- character(0)

for (i in seq_len(nrow(ca_static_layers))) {
  reg <- ca_static_layers[i, ]
  f <- if (!is.na(reg$file)) {
    ca_find_input(reg$subdir, reg$file)
  } else {
    hit <- list.files(ca_input(reg$subdir), pattern = reg$file_pattern)
    if (length(hit)) structure(hit[1], location = "scratch") else
      structure(NA_character_, location = "missing")
  }
  ca_log(sprintf("%-18s %s", reg$layer,
                 if (is.na(f)) "MISSING" else attr(f, "location")))
  if (is.na(f)) miss <- c(miss, reg$subdir)
}

for (i in seq_len(nrow(ca_zone_layers))) {
  p <- ca_find_input(ca_zone_layers$subdir[i], ca_zone_layers$file[i])
  ca_log(sprintf("%-18s %s", ca_zone_layers$layer[i],
                 if (is.na(p)) "MISSING" else attr(p, "location")))
  if (is.na(p)) miss <- c(miss, ca_zone_layers$subdir[i])
}

for (sub in unique(miss)) {
  ca_log("--- contents of ", ca_input(sub), " ---")
  f <- list.files(ca_input(sub))
  if (!length(f)) ca_log("  directory is empty")
  for (x in f) ca_log("  ", x)
}

ca_log("Stage 3d staging complete.")
