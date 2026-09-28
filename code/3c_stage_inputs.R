### CA carbon revision pipeline
### Stage 3c inputs. Copy masters from /projects to /scratch/alpine.
###
### The 3c layers are small enough to keep a backed-up master on /projects, but
### CURC prohibits job I/O there and the 3c extraction runs as an array. This
### copies each layer to scratch once, before the array is submitted.
###
### Shapefiles are copied as a file set. Copying only the .shp yields a layer
### that opens with no attributes or no CRS, which fails late and confusingly.
###
### Run once, on a compile node or as a short job. Safe to rerun, since existing
### files are skipped unless overwrite = TRUE.
###
### Usage
###   Rscript 3c_stage_inputs.R

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))

ca_stamp("stage3c_staging")

OVERWRITE <- FALSE

# ---------------------------------------------------------------------------
# VECTOR LAYERS
# ---------------------------------------------------------------------------

for (i in seq_len(nrow(ca_vector_layers))) {

  lay <- ca_vector_layers$layer[i]
  sub <- ca_vector_layers$subdir[i]
  fil <- ca_vector_layers$file[i]

  if (!file.exists(ca_master(sub, fil))) {
    ca_log(lay, "  master not found: ", ca_master(sub, fil))
    next
  }

  copied <- ca_stage_shapefile(sub, fil, overwrite = OVERWRITE)

  exts <- toupper(tools::file_ext(basename(copied)))
  ca_log(lay, "  ", length(copied), " files  ",
         paste(sort(unique(exts)), collapse = " "))

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
# MTBS RASTER STACK
# ---------------------------------------------------------------------------

src <- ca_master(ca_mtbs_layers$subdir[1])
dst <- ca_input(ca_mtbs_layers$subdir[1])

if (dir.exists(src)) {

  dir.create(dst, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(src, pattern = ca_mtbs_layers$file_pattern[1],
                      full.names = TRUE, recursive = TRUE)

  n <- 0L
  for (f in files) {
    d <- file.path(dst, basename(f))
    if (file.exists(d) && !OVERWRITE) next
    if (file.copy(f, d, overwrite = TRUE)) n <- n + 1L
  }

  ca_log("MTBS  ", length(files), " masters, ", n, " copied, ",
         length(list.files(dst, pattern = ca_mtbs_layers$file_pattern[1])),
         " now on scratch")

  # Sidecars GDAL may rely on for GeoTIFFs.
  for (ext in c("tfw", "tif.aux.xml", "tif.ovr", "tif.vat.dbf")) {
    aux <- list.files(src, pattern = paste0("\\.", gsub("\\.", "\\\\.", ext), "$"),
                      full.names = TRUE, recursive = TRUE)
    if (length(aux)) {
      file.copy(aux, file.path(dst, basename(aux)), overwrite = OVERWRITE)
      ca_log("  ", length(aux), " ", ext, " sidecars copied")
    }
  }

} else {
  ca_log("MTBS master directory not found: ", src)
}

# ---------------------------------------------------------------------------
# VERIFY
# ---------------------------------------------------------------------------

ca_log("--- verification, resolved from scratch first ---")
for (i in seq_len(nrow(ca_vector_layers))) {
  p <- ca_find_input(ca_vector_layers$subdir[i], ca_vector_layers$file[i])
  ca_log(sprintf("%-16s %s", ca_vector_layers$layer[i],
                 if (is.na(p)) "MISSING" else attr(p, "location")))
}

ca_log("Stage 3c staging complete.")
