### CA carbon revision pipeline
### Stage 0. Global configuration
###
### Target system. CURC Alpine (University of Colorado Boulder).
### Sourced by every downstream script. Single authority for paths, years,
### LULC classes, source adapters, outcome registry, group codes, screens,
### matching specification, and shared helpers.
###
### FILESYSTEM POLICY. CURC enforces this and terminates violating jobs.
###   /home/$USER            2 GB, backed up. Scripts and job files only.
###   /projects/$USER        250 GB, backed up. Master data, archived results.
###   /scratch/alpine/$USER  10 TB, NOT backed up, purged 90 days after file
###                          creation. All compute I/O goes here.
###   $SLURM_SCRATCH         Node-local SSD, deleted at job end. Temp files.

# ---------------------------------------------------------------------------
# 0.1  USER AND ROOTS
# ---------------------------------------------------------------------------

ca_user  <- Sys.getenv("USER", unset = "mich9173")

ca_paths <- list(
  code    = file.path("/home", ca_user),
  log     = file.path("/home", ca_user, "log"),
  archive = file.path("/projects", ca_user, "CA_carbon"),
  scratch = file.path("/scratch/alpine", ca_user, "CA_carbon")
)

## --- host resolution: Alpine vs Alderaan --------------------------------
ca_host <- Sys.getenv("CA_HOST", unset = "")
if (!nzchar(ca_host)) {
  ca_host <- if (grepl("alderaan", Sys.info()[["nodename"]], fixed = TRUE))
    "alderaan" else "alpine"
}
if (identical(ca_host, "alderaan")) {
  ca_paths$archive <- "/data001/projects/chungmin/CA_carbon"
  ca_paths$scratch <- "/data002/scratch/chungmin/CA_carbon"
}

ca_paths$archive <- Sys.getenv("CA_ARCHIVE", unset = ca_paths$archive)
ca_paths$scratch <- Sys.getenv("CA_SCRATCH", unset = ca_paths$scratch)

ca_input <- function(...) file.path(ca_paths$scratch, "input", ...)
ca_work <- function(...) file.path(ca_paths$scratch, "work", ...)
ca_out <- function(...) file.path(ca_paths$scratch, "output", ...)
ca_master <- function(...) file.path(ca_paths$archive, "input", ...)
ca_meta <- function(...) file.path(ca_paths$archive, "meta", ...)
ca_archive_out <- function(...) file.path(ca_paths$archive, "output", ...)

ca_subdirs <- list(
  input = c(
    "nlcd", "almanac", "ncsda2022", "emapr", "lemma", "mtbs", "thinning",
    "cpad", "offsets", "ownership", "tribal", "prism", "dem", "worldpop",
    "travel", "ecoregion", "huc", "boundary"
  ),
  work = c(
    "grid", "extract", "clean", "panel", "pretrt", "groups", "sample",
    "match", "did_input"
  ),
  output = c("did", "balance", "figures", "tables")
)

ca_init_dirs <- function(verbose = TRUE) {
  paths <- c(
    ca_paths$log,
    ca_paths$archive, ca_master(), ca_meta(), ca_archive_out(),
    file.path(ca_paths$archive, "input", ca_subdirs$input),
    ca_paths$scratch, ca_input(), ca_work(), ca_out(),
    file.path(ca_paths$scratch, "tmp"),
    file.path(ca_paths$scratch, "input", ca_subdirs$input),
    file.path(ca_paths$scratch, "work", ca_subdirs$work),
    file.path(ca_paths$scratch, "output", ca_subdirs$output)
  )
  ok <- vapply(paths, dir.create, logical(1),
               recursive = TRUE, showWarnings = FALSE)
  if (verbose) message("Directories created: ", sum(ok), " of ", length(paths))
  invisible(paths)
}

# Node-local temp when available, scratch otherwise. Never /home, whose 2 GB
# quota would be exhausted by a single terra spill file.
ca_tmpdir <- function() {
  sl <- Sys.getenv("SLURM_SCRATCH", unset = "")
  if (nzchar(sl) && dir.exists(sl)) return(sl)
  d <- file.path(ca_paths$scratch, "tmp")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Stage a master file from /projects to /scratch/alpine. Call once before an
# array job, never inside array tasks.
ca_stage <- function(relative, overwrite = FALSE) {
  src <- ca_master(relative)
  dst <- ca_input(relative)
  if (!file.exists(src)) stop("Master file not found: ", src)
  if (file.exists(dst) && !overwrite) return(invisible(dst))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dst, overwrite = TRUE, recursive = dir.exists(src))
  invisible(dst)
}

# Shapefiles are a file set, not a file. Copying only the .shp yields a layer
# that opens with no attributes or no CRS, which fails late and confusingly.
ca_shp_sidecars <- c("shp", "shx", "dbf", "prj", "cpg", "sbn", "sbx",
                     "qix", "fix", "shp.xml", "qmd")

ca_stage_shapefile <- function(subdir, file, overwrite = FALSE) {
  stem <- sub("\\.shp$", "", file, ignore.case = TRUE)
  src_dir <- ca_master(subdir)
  dst_dir <- ca_input(subdir)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- character(0)
  for (ext in ca_shp_sidecars) {
    s <- file.path(src_dir, paste0(stem, ".", ext))
    if (!file.exists(s)) next
    d <- file.path(dst_dir, basename(s))
    if (file.exists(d) && !overwrite) { copied <- c(copied, d); next }
    if (file.copy(s, d, overwrite = TRUE)) copied <- c(copied, d)
  }
  invisible(copied)
}

# Resolve an input path. Scratch is preferred, because that is where array jobs
# must read from. The /projects master is the fallback, which is fine for a
# single-task probe but not for an array. Returns NA when neither exists, with
# the location recorded as an attribute either way.
ca_find_input <- function(subdir, file = NULL, require_scratch = FALSE) {
  s <- if (is.null(file)) ca_input(subdir) else ca_input(subdir, file)
  a <- if (is.null(file)) ca_master(subdir) else ca_master(subdir, file)
  if (file.exists(s)) return(structure(s, location = "scratch"))
  if (file.exists(a)) {
    if (require_scratch) {
      stop("Found only on /projects: ", a,
           "\nStage it to scratch before running an array job. ",
           "See ca_stage_shapefile() and 3c_stage_inputs.R.")
    }
    return(structure(a, location = "archive"))
  }
  structure(NA_character_, location = "missing")
}

# ---------------------------------------------------------------------------
# 0.2  REFERENCE GRID
# ---------------------------------------------------------------------------
# NLCD 2001 defines the analysis grid. pixel_id is the terra cell index of this
# raster, so it is stable, sortable, and invertible back to coordinates.

ca_grid <- list(
  nlcd_master = ca_master("nlcd", "Annual_NLCD_LndCov_2001_CU_C1V2_prj.tif"),
  nlcd_raster = ca_input("nlcd", "Annual_NLCD_LndCov_2001_CU_C1V2_prj.tif"),
  parquet_dir = ca_work("grid", "forest_grid_parquet"),
  gpkg_dir    = ca_work("grid", "gpkg"),
  meta_rds    = ca_meta("grid_reference.rds"),
  summary_csv = ca_meta("grid_summary.csv"),

  # Stage 3 extracts from the Parquet store, not from GeoPackage. terra::extract
  # takes a plain two-column coordinate matrix, and terra::project transforms
  # that matrix into each raster's native CRS, so no geometry object is needed
  # at any point in the pipeline. A full-grid GeoPackage of roughly 10^8 points
  # is therefore write-only. Set TRUE if a QGIS view of the whole grid is
  # wanted. The sampled subset is written to GeoPackage at stage 7 regardless.
  write_gpkg  = FALSE,

  chunk_rows  = NULL
)

ca_ref_crs <- function() {
  if (file.exists(ca_grid$meta_rds)) return(readRDS(ca_grid$meta_rds)$crs_wkt)
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("terra required to read the reference CRS")
  }
  src <- if (file.exists(ca_grid$nlcd_raster)) ca_grid$nlcd_raster else ca_grid$nlcd_master
  terra::crs(terra::rast(src))
}

# ---------------------------------------------------------------------------
# 0.3  LAND COVER CLASSES
# ---------------------------------------------------------------------------
# Forest only. NLCD 52 Shrub and 71 Herbaceous are dropped from the pipeline.

ca_lulc <- data.frame(
  code  = c(41L, 42L, 43L),
  label = c("Decid", "Everg", "Mixed"),
  name  = c("Deciduous forest", "Evergreen forest", "Mixed forest"),
  stringsAsFactors = FALSE
)

ca_lulc_codes <- ca_lulc$code
ca_lulc_labels <- ca_lulc$label

# Point count per class, read from the grid partition and cached for the
# session. Used by the stage 3c probe to extrapolate join cardinality from a
# Decid sample to Evergreen, which carries roughly 43 times the points and is
# the class that sizes the job.
.ca_lulc_n <- new.env(parent = emptyenv())

ca_lulc_n_points <- function(lulc) {
  key <- as.character(lulc)
  if (!is.null(.ca_lulc_n[[key]])) return(.ca_lulc_n[[key]])
  n <- length(ca_grid_points(key)$pixel_id)
  assign(key, n, envir = .ca_lulc_n)
  n
}

# ---------------------------------------------------------------------------
# 0.4  YEARS
# ---------------------------------------------------------------------------
# Baseline 1985 to 1989 is reserved for pre-treatment covariates on units
# treated before 1990. Outcomes are analyzed 1990 to 2025. Treatment cohorts
# run to 2025 because MTBS, CARB offsets, and the rebuilt FACTS and CAL FIRE
# thinning dataset all now extend through 2025.

ca_years <- list(
  study      = 1985:2025,   # full extraction window
  baseline   = 1985:1989,   # pre-treatment covariates for pre-1990 treatments
  analysis   = 1990:2025,   # DiD outcome panel
  treatment  = 1990:2025,   # admissible event years, default
  pretrt_n   = 5L,          # pre-treatment window length
  pretrt_min = 3L           # minimum years accepted where coverage is short
)

# Per-arm treatment bounds. Each treatment dataset has its own coverage, so the
# admissible cohort window is stated per arm rather than globally.
ca_treat_years <- data.frame(
  arm   = c("pa", "fire", "thin", "offset"),
  first = c(1990L, 1990L, 1990L, 2012L),
  last  = c(2025L, 2025L, 2025L, 2019L),
  source_note = c(
    "CPAD establishment dates",
    "MTBS through 2025",
    "FACTS and CAL FIRE rebuild, 1986 to 2025",
    "CARB reporting period 1 start dates, observed 2012 to 2019"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 0.5  DATA SOURCE ADAPTERS
# ---------------------------------------------------------------------------
# Nodata sentinels and scale factors differ between releases. Each adapter
# carries its own conversion so downstream stages are version-agnostic.
#
# Carbon fraction 0.47 is applied uniformly to all three AGB products.
# Carbon_AGB arrives in metric tons per hectare of total aboveground live
# biomass, not carbon, and requires x100 to reach g m-2 then x0.47 to reach
# gC m-2. Carbon_GPP is already gC m-2 yr-1 and takes no carbon fraction. The
# factor-of-2 correction present in the 2022 pipeline is removed.
#
# Almanac nodata is -9999 for every layer except Fire_LCP, which uses 0 under
# the FARSITE convention. That single exception is declared at the layer level
# in ca_almanac_layers. Everywhere else zero is a valid value, so the
# C.sub[C.sub == 0] <- NA line from the submitted pipeline is retired.

ca_carbon_fraction <- 0.47

# Source-level fallbacks. A layer row in the stage 3 registries overrides
# these. verified records whether the value comes from documentation rather
# than assumption.
# Values below are set from the stage 3 verification run, which recovered the
# scaling from the data rather than from documentation. Where two products
# measure the same quantity at the same pixels, the ratio of their medians is
# the scale factor between them.
#
# ncsda2022  Almanac GPP over NCSDA GPP is 1.02 in 2000 and 1.13 in 2016, with
#            correlation 0.79 to 0.80, so NCSDA fluxes are already gC m-2 yr-1
#            and take no scale factor. NPP sits at roughly half of GPP, the
#            expected carbon-use-efficiency ratio, which corroborates it.
# emapr      Almanac Carbon_AGB over eMapR is 0.93 with correlation 0.76, so
#            eMapR is already t/ha and takes no scale factor.
# lemma      k=7 GNN biomass rasters from the LEMMA California biomass
#            project, in kg/ha rather than the Mg/ha shown on the LEMMA
#            website, confirmed by the data provider. The kg/ha to t/ha step
#            lives in ca_unit_conversion, not here, so extraction writes native
#            values for every product and all conversion happens once at stage
#            4. Zero occupies 7.383 percent
#            of pixels in 1990, 2003, and 2016 alike, a static mask rather
#            than real zero biomass, so zero is the nodata sentinel.
#            Verified at stage 4 on the Evergreen partition, 2026-08-01. After
#            conversion the three biomass products sit at 8,084 gC/m2
#            (Carbon_AGB), 8,811 (eMapR), and 7,573 (LEMMA), a spread of 1.16
#            from highest to lowest. LEMMA runs 6 percent below Carbon_AGB and
#            eMapR 9 percent above, both product differences rather than unit
#            errors. The earlier note here recorded LEMMA at roughly 20 percent
#            below, which was an estimate made before conversion was applied.
ca_sources <- data.frame(
  source     = c("almanac2026", "ncsda2022", "emapr", "lemma", "mtbs"),
  subdir     = c("almanac", "ncsda2022", "emapr", "lemma", "mtbs"),
  first_year = c(1985L, 1985L, 1990L, 1990L, 1985L),
  last_year  = c(2025L, 2021L, 2017L, 2016L, 2025L),
  nodata     = c(-9999, NA_real_, NA_real_, 0, NA_real_),
  scale      = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
  verified   = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  note       = c("v2026.1 documentation section 2.1",
                 "gC/m2/yr, cross-checked against Almanac GPP",
                 "t/ha, cross-checked against Almanac Carbon_AGB",
                 "kg/ha per LEMMA, zero is a static mask",
                 "INT1U, classes 1 to 6, background already NA, confirmed at extraction 2026-07-31"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# EXTRACT STORE
# ---------------------------------------------------------------------------
# Where an extract physically sits. This is set by the SOURCE constant of the
# script that ran the registry, NOT by the per-row source column, and the two
# differ in two places. Deriving a path from registry$source produced a wrong
# directory twice, so the mapping is stated once, keyed on layer, and every
# path in stages 4 and later is built from it.
#
#   3a_extract_carbon_2026.R   SOURCE almanac2026          ca_almanac_layers
#   3a_extract_carbon_2022.R   SOURCE ncsda2022            ca_ncsda_layers
#                              including eMapR and LEMMA, whose source column
#                              reads emapr and lemma and is used for nodata and
#                              scale lookup only
#   3b_extract_screens_2026.R  SOURCE almanac2026_screens  ca_screen_layers
#                              including Veg_TreeFrac, whose source column
#                              reads almanac2026
#   3b_extract_screens_2026.R  SOURCE mtbs                 ca_mtbs_layers
#                              under the mtbs selector
#
# Built from the registries rather than typed out, so a layer added to a
# registry inherits the right store instead of silently missing from the table.

.ca_store_map <- function() {
  rbind(
    data.frame(layer = ca_almanac_layers$layer, store = "almanac2026",
               first_year = ca_almanac_layers$first_year,
               last_year  = ca_almanac_layers$last_year,
               source     = ca_almanac_layers$source,
               active     = ca_almanac_layers$active,
               stringsAsFactors = FALSE),
    data.frame(layer = ca_ncsda_layers$layer, store = "ncsda2022",
               first_year = ca_ncsda_layers$first_year,
               last_year  = ca_ncsda_layers$last_year,
               source     = ca_ncsda_layers$source,
               active     = ca_ncsda_layers$active,
               stringsAsFactors = FALSE),
    data.frame(layer = ca_screen_layers$layer, store = "almanac2026_screens",
               first_year = ca_screen_layers$first_year,
               last_year  = ca_screen_layers$last_year,
               source     = ca_screen_layers$source,
               active     = ca_screen_layers$active,
               stringsAsFactors = FALSE),
    data.frame(layer = ca_mtbs_layers$layer, store = "mtbs",
               first_year = ca_mtbs_layers$first_year,
               last_year  = ca_mtbs_layers$last_year,
               source     = ca_mtbs_layers$source,
               active     = ca_mtbs_layers$active,
               stringsAsFactors = FALSE)
  )
}

ca_layer_registry <- function() {
  m <- .ca_store_map()
  d <- m$layer[duplicated(m$layer)]
  if (length(d)) {
    stop("Layer names must be unique across registries. Duplicated: ",
         paste(unique(d), collapse = ", "))
  }
  m
}

# Inactive registry rows are scaffolds for a layer that does not exist in the
# current data release, Carbon_NPP, Carbon_NEP, and Carbon_NBP being the three.
# They are expected to be absent from disk, so the path check reports them as
# inactive rather than as incomplete. Flagging a deliberate absence as a fault
# trains the reader to ignore the warning, which is the opposite of the point.
ca_layer_active <- function(layer) {
  m <- ca_layer_registry()
  isTRUE(m$active[match(layer, m$layer)])
}

ca_extract_store <- function(layer) {
  m <- ca_layer_registry()
  i <- match(layer, m$layer)
  if (anyNA(i)) {
    stop("No registry entry for layer(s): ",
         paste(layer[is.na(i)], collapse = ", "),
         ". Every layer read downstream must appear in one of the four ",
         "extraction registries, because the store cannot be guessed.")
  }
  m$store[i]
}

ca_layer_years <- function(layer) {
  m <- ca_layer_registry()
  i <- match(layer, m$layer)
  if (is.na(i)) stop("No registry entry for layer: ", layer)
  m$first_year[i]:m$last_year[i]
}

# ---------------------------------------------------------------------------
# STAGE 3C REGISTRIES
# ---------------------------------------------------------------------------
# Management, disturbance, and ownership layers. These are small relative to the
# carbon stacks, so masters live on /projects where they are backed up, and are
# staged to scratch with ca_stage() before any array job. The projected
# shapefiles are work product that would be tedious to rebuild, which is why
# they are worth a backed-up copy.
#
# Field names below are the expected ones. 3c_probe_vectors.R reports what is
# actually present, and nothing downstream should be written until the probe
# output has been checked against this table.

# The two CAL FIRE THP files are complementary, not redundant. THPs_Historical
# carries plan years to roughly 2010 and THPs carries 2010 onward, so using
# either alone loses half the private harvest record. Both are registered and
# stage 5 splits them on THP_YEAR rather than deduplicating, since a plan year
# belongs to exactly one file.
ca_vector_layers <- data.frame(
  layer      = c("thin_nto", "thin_thp",
                 "thin_facts_th", "thin_facts_hfr",
                 "offset", "tribal", "cpad", "ownership"),
  subdir     = c("thinning", "thinning", "thinning", "thinning",
                 "offsets", "tribal", "cpad", "ownership"),
  file       = c("NTOs_prj.shp",
                 "THPs_Historical_prj.shp",
                 "S_USA.Actv_TimberHarvest_CA_prj.shp",
                 "S_USA.Actv_HazFuelTrt_PL_CA_prj.shp",
                 "CA_forest_offset_compliance_prj.shp",
                 "CA_Federally_Recognized_Tribal_Lands_prj_yr.shp",
                 "CPAD_2026a_SuperUnits_prj.shp",
                 "Land_Ownership_prj.shp"),
  # CAL FIRE splits the timber harvest plan archive by permit era, not by
  # content. THPs_Historical_prj.shp carries THP_YEAR up to 2010 and
  # THPs_prj.shp carries 2011 onward, with identical schemas on every field
  # this pipeline reads. They are one layer covering one programme, so they are
  # stacked at read time rather than carried as two registry rows. Keeping them
  # separate would produce two treatment layers that must never be double
  # counted, which is a bug waiting to happen at stage 6.
  file_extra = c(NA, "THPs_prj.shp", NA, NA, NA, NA, NA, NA),
  date_field = c("COMPLETED", "COMPLETED", "DATE_COMPL",
                 "DATE_COMPL", NA, NA, NA, NA),
  # CAL FIRE plans carry an explicit status. Unlogged and Approved plans have a
  # placeholder COMPLETED date of 1899-12-30, which is why date_missing reads
  # zero while 37 percent of THP falls outside the study window. Only Completed
  # plans are treatment.
  #
  # FACTS does carry an equivalent flag, contrary to the earlier note here. The
  # probe shows STAGE_DESC on the timber harvest layer with values Accomplished,
  # Layout, and NEPA, so Layout and NEPA records are planned rather than
  # executed and must not enter the treatment arm. The hazardous fuels layer
  # carries STAGE_VALU with Accomplished as its only value, so the filter is a
  # no-op there and is stated anyway so the rule is uniform.
  status_field = c("PLAN_STAT", "PLAN_STAT", "STAGE_DESC", "STAGE_VALU",
                   NA, NA, NA, NA),
  status_keep  = c("Completed", "Completed", "Accomplished", "Accomplished",
                   NA, NA, NA, NA),
  # Placeholder dates that must be read as missing, never as a year.
  date_null    = c("1899-12-30", "1899-12-30", "0", "0",
                   NA, NA, NA, NA),
  # Second silvicultural prescription. Not used for intensity, but a polygon
  # carrying one is disturbed and must be screened out of the control pool.
  key_field_2  = c("SILVI_2", "SILVI_2", NA, NA, NA, NA, NA, NA),
  # FACTS records the method separately from the activity name. Where the two
  # disagree, for instance an activity that reads mechanical against
  # METHOD = Prescribed Burn, the method field is the more reliable evidence.
  method_field = c(NA, NA, "METHOD_DES", "METHOD", NA, NA, NA, NA),
  # Establishment or cohort year.
  year_field   = c("NTO_YEAR", "THP_YEAR", NA, NA, NA, "Year",
                   "YR_EST", NA),
  # Confirmed against probe_3c_fields_20260731.csv. CPAD and ownership are
  # classified by a managing agency level held as a character label, not as the
  # numeric code the submitted pipeline recoded to. MNG_AG_LEV carries 10
  # distinct values and Own_Level carries 7, matching ca_cpad_agency and
  # ca_ownership_level respectively.
  key_field  = c("SILVI_1", "SILVI_1", "ACTIVITY_N", "ACTIVITY",
                 "ARB_id", NA, "MNG_AG_LEV", "Own_Level"),
  code_field = c(NA, NA, "ACTIVITY_2", "ACTIVITY_C", NA, NA, NA, NA),
  crosswalk  = c("S5", "S5", "S4", "S4", NA, NA, NA, NA),
  role       = c("treatment", "treatment", "treatment",
                 "treatment", "treatment", "exclusion", "treatment",
                 "stratifier"),
  active     = TRUE,
  purpose    = c("CAL FIRE notices of timber operations",
                 "CAL FIRE timber harvest plans, historical and current",
                 "FACTS timber harvest",
                 "FACTS hazardous fuel treatment polygons",
                 "CARB compliance offset project boundaries",
                 "federally recognised tribal lands, exclusion mask only",
                 "protected areas, all managing agency levels",
                 "public land ownership, private is the complement"),
  stringsAsFactors = FALSE
)

# MTBS is a raster stack and runs through the same engine as stage 3a and 3b.
# Severity classes 2, 3, and 4 are low, moderate, and high, per ca_fire_severity.
ca_mtbs_layers <- data.frame(
  layer        = "MTBS",
  subdir       = "mtbs",
  dir          = ".",
  file_pattern = "^mtbs_CA_[0-9]{4}_prj\\.tif$",
  source       = "mtbs",
  band         = NA_integer_,
  # 1984 is absent from the staged rasters. The MTBS thematic product begins in
  # 1984 nationally, but the California subset staged here starts in 1985, which
  # the probe confirms as 41 files with no interior gaps. Leaving 1984 in the
  # registry produced a 42-task array whose first task had no file to read.
  # 1985 also matches the Almanac and NCSDA start, so no arm loses a year.
  first_year   = 1985L,
  last_year    = 2025L,
  nodata       = NA_real_,
  scale        = NA_real_,
  active       = TRUE,
  units        = "MTBS severity class",
  purpose      = "wildfire severity and event year",
  stringsAsFactors = FALSE
)

# Knight et al. (2022, J. Environ. Manage.) treatment intensity crosswalk.
# Table S4 covers FACTS timber harvest and hazardous fuel treatments.
# Table S5 covers CAL FIRE timber harvest plans and NTMPs.
# Transcribed verbatim from the published supplement. activity_key is the
# normalised join key, see ca_norm_activity().
ca_knight_csv <- ca_meta("knight_intensity_crosswalk.csv")

# Activity strings observed in the shapefiles that are absent from the published
# tables. Two resolution types.
#   alias     a wording variant of a published activity, inheriting its class
#   assigned  not in the published tables, classified by a stated rationale
# Kept in a separate file so knight_intensity_crosswalk.csv remains a verbatim
# transcription of the published supplement. The response letter can then state
# exactly how many activities came from Knight and how many were resolved here.
ca_alias_csv <- ca_meta("activity_alias.csv")

# CARB offset project start years, derived from the ARB Offset Credit Issuance
# workbook. The offset shapefile carries no date field, so the cohort year for
# the ONP pairing comes from this lookup, joined on ARB_id.
#
# start_year is the calendar year of the first reporting period start date. For
# a compliance Improved Forest Management project that is the project
# commencement date, which is the point management changed, and therefore the
# correct treatment year.
#
# Two alternatives were rejected. Vintage disagrees with the start year for 116
# of 144 forest compliance projects, because a reporting period spanning a year
# boundary is credited to the later vintage. Issuance date lags the start by one
# to eight years, median two, since it records when CARB issued credits rather
# than when management changed.
#
# Early Action projects are included in the lookup because a compliance project
# can have an Early Action predecessor whose management change is years earlier.
# ea_cop distinguishes them.
# The offset shapefile carries compliance projects only, so the join uses the
# COP rows. Early Action rows are retained in the lookup as a diagnostic: a
# compliance project can have an Early Action predecessor under which
# management already changed years earlier, which would put treated years
# inside the 2007 to 2011 pre-treatment window.
ca_offset_years_csv <- ca_meta("carb_offset_start_years.csv")
ca_offset_join_field <- "ARB_id"
ca_offset_year_field <- "start_year"
ca_offset_program <- "COP"


# Activity strings differ between the published tables and the shapefile
# attributes in case, spacing, punctuation, and occasionally truncation. The
# join runs on a normalised key rather than the raw string. Normalisation is
# deliberately aggressive, since it only has to separate genuinely different
# activities, and there are fewer than 150 of them.
ca_norm_activity <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]+", " ", x)
  trimws(gsub("\\s+", " ", x))
}

# Record retention. Three tiers, not two. Every archival record is kept through
# extraction. What differs is the role it plays at stage 6.
#
#   treatment_LMH     Low, Medium, or High intensity AND a mechanical method.
#                     These are the thinning treatment pixels.
#
#   disturbance_only  Everything else that touched the ground: the Variable and
#                     Unknown intensity classes, prescribed fire, chemical,
#                     planting, and site preparation. These pixels are never
#                     treatment, and they are also disqualified from the
#                     undisturbed control pool. A pixel that was underburned is
#                     not an undisturbed pixel, even though its intensity
#                     cannot be resolved to Low, Medium, or High.
#
#   non_disturbing    Knight assigns intensity N/A: surveys, inventories,
#                     monitoring, marking, layout, and range or watershed
#                     administration. These record the presence of an observer,
#                     not a disturbance, so they disqualify nothing. Using
#                     Knight's own N/A designation rather than a method
#                     judgement keeps this defensible in the response letter.
#
#   wildfire          FACTS activities that record wildfire rather than
#                     treatment. Dropped entirely, because MTBS is the wildfire
#                     source and retaining these would double count.
#
# The record screen is secondary to the detection screen. Any pixel that
# actually lost canopy is caught by Almanac Disturbance_TreeFrac regardless of
# what the archives say. The record screen adds coverage for low-severity
# treatment that COLD may not detect, which is exactly the prescribed fire case.

ca_intensity_analysis <- c("Low", "Medium", "High")

ca_record_classes <- data.frame(
  record_class = c("treatment_LMH", "disturbance_only", "non_disturbing",
                   "wildfire"),
  treatment    = c(TRUE, FALSE, FALSE, FALSE),
  disqualifies_control = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Methods admitted to the thinning treatment arm. Prescribed fire is a
# disturbance and a control disqualifier, but it is not forest thinning, so it
# is excluded from the treatment arm even where Knight assigns it a clean
# intensity. Broadcast Burning at Medium is the case this rule exists for.
ca_treatment_methods <- "mechanical"

# Five activities carry an override in the crosswalk, moving them out of the
# treatment arm despite satisfying the intensity and method rule. The reason is
# recorded per row in the override_reason column so the deviation is auditable.
#   No Harvest Area          designated no-harvest area within a THP
#   Rearrangement of Fuels   redistributes fuels, removes no basal area
#   Chipping of Fuels        processes existing material, removes no basal area
#   Permanent Land Clearing  land-use change, no recovery trajectory
#   Conversion               land-use change, no recovery trajectory
# All five remain control disqualifiers. A cleared or chipped pixel is not an
# undisturbed pixel.
#
# Final counts: 40 activities in the treatment arm, 24 from Table S4 and 16
# from Table S5, of 141 published.

# ---------------------------------------------------------------------------
# STAGE 3 LAYER REGISTRIES
# ---------------------------------------------------------------------------
# One row per raster stack. Columns are consumed by ca_extract.R.
#   subdir        folder under input/
#   dir           folder under that subdir holding the annual rasters
#   file_pattern  regex narrowing the file list, NA to take every raster
#   band          band index for multi-band stacks, NA for single band
#   nodata        sentinel applied on top of the file's own NA flag
#   scale         multiplier applied after nodata, NA leaves native units
#   active        FALSE keeps the row as a scaffold without expanding it into
#                 the array
#
# All values below are from Wildland Almanac California v2026.1 documentation,
# sections 2 and 3, not inferred. Every file is a single-band COG in EPSG:5070
# at 30 m on the same continental grid as NLCD, so no coordinate transformation
# is needed for the Almanac. Nodata is -9999 everywhere except Fire_LCP, which
# follows the FARSITE convention of 0.

# 3a, 2026 Wildland Almanac.
# Carbon_AGB is total aboveground live biomass in metric tons per hectare, not
# carbon. Conversion to gC m-2 is x100 then x0.47, applied at stage 4.
# Carbon_GPP is already gC m-2 yr-1 and takes no carbon fraction.
# Canopy height is Fire_LCP band CH, delivered as its own single-band COG in
# decimeters. Scale 0.1 converts to metres.
#
# Fire_LCP nodata is documented as 0 under the FARSITE convention, but zero is
# not applied as a sentinel here. Verification showed that pixels with CH = 0
# carry a median 58 t/ha of biomass and 31 percent tree cover, so they are real
# forest, not masked ground. Two further facts settle it. The zero fraction
# rises monotonically from 13.81 percent in 1985 to 14.62 percent in 2025,
# tracking canopy loss, whereas a data mask would be static as LEMMA's is. And
# 304,040 of the 305,479 zero pixels carry a valid Carbon_AGB value, so they
# sit inside the Almanac wildland mask.
#
# CH = 0 therefore means no canopy stratum in the fuel model, which is
# information rather than absence, and it is carried through as 0 metres.
# Applying the documented sentinel would have dropped 14 percent of pixels from
# the augmented matching arm, and dropped them non-randomly, biasing the
# matched sample toward tall-canopy stands. Genuinely masked pixels are removed
# at stage 4 through ca_clean_rules instead, keyed on Carbon_AGB being NA.
ca_almanac_layers <- data.frame(
  layer        = c("Carbon_AGB", "Carbon_GPP", "Fire_LCP_CH",
                   "Carbon_NPP", "Carbon_NEP", "Carbon_NBP"),
  subdir       = "almanac",
  dir          = c("Carbon_AGB", "Carbon_GPP", "Fire_LCP",
                   "Carbon_NPP", "Carbon_NEP", "Carbon_NBP"),
  file_pattern = c(NA, NA, "_Fire_LCP_CH_[0-9]{4}\\.tif$", NA, NA, NA),
  source       = "almanac2026",
  band         = NA_integer_,
  first_year   = 1985L,
  last_year    = c(2025L, 2025L, 2025L, 2025L, 2025L, 2025L),
  nodata       = c(-9999, -9999, NA_real_, -9999, -9999, -9999),
  scale        = c(NA_real_, NA_real_, 0.1, NA_real_, NA_real_, NA_real_),
  active       = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  units        = c("t/ha", "gC/m2/yr", "m after scaling",
                   "reserved", "reserved", "reserved"),
  purpose      = c("carbon stock outcome and pre-treatment covariate",
                   "flux outcome and pre-treatment covariate",
                   "canopy height, pre-treatment covariate",
                   "not present in v2026.1, scaffold for a future release",
                   "not present in v2026.1, scaffold for a future release",
                   "not present in v2026.1, scaffold for a future release"),
  stringsAsFactors = FALSE
)

# 3a, 2022 NCSDA plus the two legacy aboveground products. Grouped because all
# three are local downloads read the same way, unlike the streamed Almanac.
# NPP, NEP, and NBP are the operative flux outcomes, since v2026.1 carries only
# AGB and GPP under its carbon theme. GPP is extracted for cross-comparison
# against Almanac Carbon_GPP.
ca_ncsda_layers <- data.frame(
  layer        = c("GPP", "NPP", "NEP", "NBP", "eMapR", "LEMMA"),
  subdir       = c("ncsda2022", "ncsda2022", "ncsda2022", "ncsda2022",
                   "emapr", "lemma"),
  dir          = c("CFlux_GPP", "CFlux_NPP", "CFlux_NEP", "CFlux_NBP",
                   "tif", "tif"),
  file_pattern = NA_character_,
  source       = c("ncsda2022", "ncsda2022", "ncsda2022", "ncsda2022",
                   "emapr", "lemma"),
  band         = NA_integer_,
  # NBP starts 1986. The 1985 source raster is truncated, 13 MB against
  # the 500 to 600 MB of a complete year, and cannot be repaired. Nothing
  # is lost: the analysis window opens in 1990 and NBP is not a matching
  # covariate, so 1985 was never used.
  first_year   = c(1985L, 1985L, 1985L, 1986L, 1990L, 1990L),
  last_year    = c(2021L, 2021L, 2021L, 2021L, 2017L, 2016L),
  nodata       = NA_real_,
  scale        = NA_real_,
  active       = TRUE,
  units        = c("probe", "probe", "probe", "probe", "probe", "probe"),
  purpose      = c("flux cross-comparison against Almanac GPP",
                   "NPP tree outcome",
                   "NEP outcome",
                   "NBP outcome, wildfire SI only",
                   "aboveground biomass, secondary product",
                   "aboveground biomass, secondary product"),
  stringsAsFactors = FALSE
)

# 3b, 2026 Almanac screens. Not outcomes.
# Disturbance layers report annual loss at COLD-detected pixels for water years
# 1986 to 2024. Boundary years 1985 and 2025 are not produced, because each
# estimate needs a bracketing pre and post year. Zero means a valid pixel with
# no disturbance that year and must never be recoded to NA.
# Disturbance timing carries a detection-lag correction and may be attributed
# to a slightly different year than first flagged, which is why the thinning
# detection window is event year plus or minus one.
ca_screen_layers <- data.frame(
  layer        = c("Disturbance_TreeFrac", "Disturbance_AGB", "Veg_TreeFrac"),
  subdir       = "almanac",
  dir          = c("Disturbance_TreeFrac", "Disturbance_AGB", "Veg_TreeFrac"),
  file_pattern = NA_character_,
  source       = "almanac2026",
  band         = NA_integer_,
  first_year   = c(1986L, 1986L, 1985L),
  last_year    = c(2024L, 2024L, 2025L),
  nodata       = -9999,
  scale        = c(1e-4, 0.1, 1e-4),
  active       = TRUE,
  units        = c("absolute delta tree fraction, positive is loss",
                   "delta t/ha, positive is loss",
                   "tree fractional cover, 0 to 1 after scaling"),
  purpose      = c("control verification and thinning detection screen",
                   "growth versus loss decomposition, secondary screen",
                   "NLCD robustness check, continuous forest cover"),
  stringsAsFactors = FALSE
)

# Fire_LCP also ships static Elevation, Slope, and Aspect bands with no year in
# the filename. These are candidates to replace the separately derived
# topographic covariates at stage 3d, which would put every matching covariate
# on one co-registered stack. Handled there, not here.

# CRS. Every raster product in this pipeline sits on NAD83 Conus Albers, but
# they declare it differently. The Almanac states EPSG:5070. eMapR carries
# Albers_Conic_Equal_Area with a Custom authority. LEMMA carries
# National_Albers with a Custom authority. NLCD 2001 carries the same
# parameters as a proj4 string. All four share lat_0 23, lon_0 -96, lat_1 29.5,
# lat_2 45.5, x_0 0, y_0 0, NAD83, metres. No coordinate transformation is
# required anywhere in stage 3. Comparison uses terra::same.crs rather than
# string equality, because the three declarations are textually different and
# numerically identical.

# Raw rasters live on scratch, not /projects. The full carbon input set is
# roughly 290 GB against a 250 GB /projects quota, and CURC forbids job I/O on
# /projects regardless. Scratch purges 90 days after file creation, so
# ca_input_manifest() writes a file list and sizes to the backed-up meta
# directory, making a re-download verifiable rather than guesswork.
ca_input_manifest <- function() {
  dirs <- file.path(ca_paths$scratch, "input", ca_subdirs$input)
  files <- unlist(lapply(dirs[dir.exists(dirs)], list.files,
                         recursive = TRUE, full.names = TRUE))
  if (!length(files)) {
    message("No input files found under ", ca_input())
    return(invisible(NULL))
  }
  out <- data.frame(
    path = sub(paste0("^", ca_input(), "/?"), "", files),
    size_mb = round(file.size(files) / 1e6, 1),
    mtime = format(file.mtime(files), "%Y-%m-%d"),
    stringsAsFactors = FALSE
  )
  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  path <- ca_meta("input_manifest.csv")
  write.csv(out, path, row.names = FALSE)
  message("Manifest written: ", path, "  files ", nrow(out),
          "  total ", round(sum(out$size_mb) / 1000, 1), " GB")
  invisible(out)
}

# ---------------------------------------------------------------------------
# STAGE 4 UNIT CONVERSION
# ---------------------------------------------------------------------------
# Extraction writes native values for every product. All conversion happens
# here, once, so a factor can be checked in one place instead of being chased
# through the extraction code. factor is applied first, then carbon_fraction
# where TRUE.
#
# The three biomass products converge on gC m-2:
#   Carbon_AGB  t/ha  x100 -> g m-2  x0.47 -> gC m-2
#   eMapR       t/ha  x100 -> g m-2  x0.47 -> gC m-2
#   LEMMA       kg/ha x0.1 -> g m-2  x0.47 -> gC m-2   (0.001 to t/ha, then 100)
# The flux products are already gC m-2 yr-1 and take neither.

# Layers carrying a class code rather than a measurement. They have no units,
# no conversion, and no carbon fraction, so they are deliberately absent from
# ca_unit_conversion and must never be read through ca_read_layer(), which
# errors on a missing conversion entry by design. Use ca_read_categorical().
#
# MTBS is INT1U with classes 1 to 6, so the INT2S limit check does not apply to
# it either. The sentinel section skips this list rather than special-casing
# one layer name, so a future class layer inherits the same treatment.
ca_categorical_layers <- c("MTBS")

ca_is_categorical <- function(layer) layer %in% ca_categorical_layers

ca_unit_conversion <- data.frame(
  layer  = c("Carbon_AGB", "eMapR", "LEMMA",
             "Carbon_GPP", "GPP", "NPP", "NEP", "NBP",
             "Disturbance_AGB", "Disturbance_TreeFrac", "Veg_TreeFrac",
             "Fire_LCP_CH"),
  native = c("t/ha", "t/ha", "kg/ha",
             "gC/m2/yr", "gC/m2/yr", "gC/m2/yr", "gC/m2/yr", "gC/m2/yr",
             "delta t/ha", "delta fraction", "fraction", "m"),
  factor = c(100, 100, 0.1,
             1, 1, 1, 1, 1,
             1, 1, 1, 1),
  carbon_fraction = c(TRUE, TRUE, TRUE,
                      FALSE, FALSE, FALSE, FALSE, FALSE,
                      FALSE, FALSE, FALSE, FALSE),
  target = c("gC/m2", "gC/m2", "gC/m2",
             "gC/m2/yr", "gC/m2/yr", "gC/m2/yr", "gC/m2/yr", "gC/m2/yr",
             "delta t/ha", "delta fraction", "fraction", "m"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# STAGE 4 CLEANING RULES
# ---------------------------------------------------------------------------
# Anomalies found by 3_verify_extract.R, each with a stated rule rather than
# left to propagate. All are applied at stage 4, after extraction, so extracted
# values stay in native form and every correction is auditable in one place.
# Counts are from the Decid partition, 2.2 million points.

# LEMMA zero is not listed here. It is a nodata sentinel, so it is applied at
# extraction through ca_sources rather than as a cleaning rule.
# USE RESTRICTIONS created by these rules.
#
#   Veg_TreeFrac  clipping 3.0 to 1.0 is defensible only because the layer is
#                 used for one binary question, whether an NLCD 2001 forest
#                 pixel carried trees. Under any threshold for that question a
#                 clipped and an unclipped pixel answer identically, so the
#                 clip changes no downstream decision. Promote this layer to a
#                 continuous covariate or an outcome and the clipped value
#                 becomes a fabricated number. It is verification only.
#
#   Carbon_GPP    clipping negatives to 0 breaks log(), which is why the GPP
#                 outcome uses log1p rather than log. At values near 1,000
#                 gC m-2 yr-1 the two differ by about 0.1 percent, and log1p is
#                 defined at zero. Recorded in ca_outcomes usage, applied at
#                 stage 9. The alternative, setting negatives NA, was rejected:
#                 it removes observations selected on the post-treatment
#                 outcome, in the treated arm, concentrated in the most
#                 severely affected units, which is outcome-dependent attrition
#                 inside a DiD and attenuates the ATT harder than the clip does
#                 while being far less visible.
# ORDER OF OPERATIONS. mask_sentinel runs on STORED values, before unit
# conversion. Every other rule runs on CONVERTED values, after. A type limit
# lives in stored counts, so masking it after multiplying by 100 and 0.47 would
# compare against a number the data can no longer hold. ca_read_layer()
# enforces this order.
#
# Verified on the Evergreen partition, 2026-08-01, four probe years per layer.
#
#   Carbon_GPP            5 pixels at -32768 in every year. Static, so a mask
#                         rather than a fire artifact. The other 81 to 1,104
#                         negatives are genuine small values and are clipped
#   NBP                   1 to 55 pixels at +32767. 32,767 gC m-2 yr-1 is 327
#                         tC/ha in one year and is not physical. The near-limit
#                         negatives, obs_min -32,750, are left alone because
#                         -160 tC/ha in a stand-replacing fire year is
#                         plausible in a high-biomass stand
#   Veg_TreeFrac          2 to 72 pixels at +32767. The other ~13,000 above 1.0
#                         are genuine overshoot and are clipped
#   Disturbance_TreeFrac  1 pixel at -32767 in 2024. Note the convention here
#                         is -32767, not -32768
#
# eMapR and LEMMA are deliberately absent. LEMMA stores to 3,488,625 kg/ha, so
# 32767 is an ordinary value it passes through and the 247 to 356 pixels
# carrying it are real. eMapR maxes at 707 t/ha and cannot reach the limit at
# all. Masking either would delete data.
#
# Carbon_AGB maxes at exactly 800 t/ha in all four probe years, which is a
# producer-applied cap rather than a sentinel. Nothing to mask, but the upper
# tail of the biomass distribution is truncated by the product and that belongs
# in the limitations.
ca_sentinel_rules <- data.frame(
  layer = c("Carbon_GPP", "NBP", "Veg_TreeFrac", "Disturbance_TreeFrac"),
  value = c(-32768, 32767, 32767, -32767),
  scale = c(1, 1, 1e-4, 1e-4),   # extraction scale, so value is in stored units
  n_per_year = c("5", "1 to 55", "2 to 72", "1 in 2024"),
  stringsAsFactors = FALSE
)

# Applied after conversion. LEMMA zero is not listed: it is a nodata sentinel
# handled at extraction through ca_sources.
#
# NCSDA GPP is added alongside Almanac Carbon_GPP. Gross primary production
# cannot be negative in either product, and the submitted pipeline screened
# neither. 1,955 to 29,927 pixels per year, reaching -166 gC m-2 yr-1.
#
# NPP is deliberately NOT clipped. Negative NPP is physically real where
# respiration exceeds photosynthesis, so a clip would erase signal rather than
# error. See ca_pending() for the consequence, which is that log() of a
# negative is NaN and npp_tree currently requests a log transform.
ca_clean_rules <- data.frame(
  layer = c("Veg_TreeFrac", "Veg_TreeFrac", "Carbon_GPP", "GPP",
            "Fire_LCP_CH"),
  rule  = c("clip_upper", "clip_lower", "clip_lower", "clip_lower",
            "mask_by_agb"),
  value = c(1, 0, 0, 0, NA),
  note = c(
    "fractional cover cannot exceed 1; genuine overshoot between 1 and 3.28",
    "fractional cover cannot be negative; 8 to 233 pixels per year",
    "gross primary production cannot be negative; 81 to 1,104 per year after the sentinel is masked",
    "gross primary production cannot be negative; 1,955 to 29,927 per year, NCSDA product",
    "Fire_LCP covers the full extent with no missing values at all, so CH is set NA where Carbon_AGB is NA, which is the Almanac wildland mask; CH = 0 elsewhere is no canopy"
  ),
  stringsAsFactors = FALSE
)

# Points per extraction chunk. Peak memory is set by this, not by class size.
ca_extract_chunk <- 1e7L

# Climate covariates are PRISM 1991 to 2020 30-year normals, replacing the
# 2000 to 2002 three-year window. Normals are geographic descriptors of pixel
# comparability, not longitudinal variables, so the pre-treatment timing rule
# does not apply to them. The same holds for elevation, slope, aspect,
# population density, and travel time.
ca_prism <- list(
  normals_period = "1991-2020",
  variables      = c("ppt", "tmean"),
  native_res_m   = 4000,
  resample       = "bilinear"
)

# ---------------------------------------------------------------------------
# 0.5d  STAGE 3D COVARIATE LAYERS
# ---------------------------------------------------------------------------
# Time-invariant matching covariates. The 3c and 3d boundary is dimensional,
# not thematic. Stage 3c layers carry an event year and produce pixel-year
# records. Stage 3d layers produce exactly one row per pixel with no year
# dimension, which is why the ecoregion and HUC joins sit here rather than with
# the other vector layers.
#
# These are geographic descriptors of pixel comparability, not longitudinal
# variables, so the pre-treatment timing rule of Ho et al. (2011) does not
# apply. Only gpp_pre5, agb_pre5, and ch_pre5 carry that requirement, and they
# are built at stage 8 from the 3a panels.
#
# Topography is SRTM 90 m, as submitted. The Almanac Fire_LCP static Elevation,
# Slope, and Aspect bands are deliberately not used, since carrying the same
# covariate from two sources buys nothing and doubles the reconciliation work.
#
# Every file carries a _prj or _proj suffix but no CRS is assumed. Rasters are
# never reprojected. Points are transformed into the raster CRS instead, and
# stage 3d is the first stage where that actually fires, because PRISM at 30
# arc-seconds is a geographic grid whatever the filename says.
#
# method is the terra::extract interpolation. Every 3d raster is coarser than
# the 30 m grid, so the choice matters. Slope and aspect take nearest because
# both are derived from the 90 m elevation neighbourhood, and interpolating a
# derivative produces values the source DEM does not support. Aspect is
# additionally circular, so bilinear across the 0 to 360 wrap would return
# south for a north-facing pixel.
#
# reduce applies where a layer resolves to more than one file. WorldPop is the
# 2000 to 2002 mean, as submitted. Travel time keeps settlement class 11 only,
# the aggregate 50,000 to 50,000,000 threshold. Classes 4 and 5 are staged but
# inactive.

ca_static_layers <- data.frame(
  layer   = c("elevation", "slope", "aspect",
              "ppt_normal", "tmean_normal",
              "pop_density", "city_travel_time"),
  subdir  = c("dem", "dem", "dem", "prism", "prism", "worldpop", "travel"),
  file    = c("SRTM_90m_elevation_CA_prj.tif",
              "SRTM_90m_slope_CA_prj.tif",
              "SRTM_90m_aspect_CA_prj.tif",
              "prism_ppt_us_30s_2020_avg_30y_prj.tif",
              "prism_tmean_us_30s_2020_avg_30y_prj.tif",
              NA_character_,
              NA_character_),
  # Confirmed against the scratch listing on 2026-07-31, not from the download
  # notes. The travel time files carry no _prj suffix, so unlike every other
  # input in this pipeline they may sit on a geographic grid. The probe reports
  # same_crs for that layer, and ca_extract_points transforms the points into
  # the raster CRS if it differs. No raster is reprojected either way.
  file_pattern = c(NA, NA, NA, NA, NA,
                   "^CA_pd_200[0-2]_1km_UNadj_prj\\.tif$",
                   "^CA_travel_time_to_cities_11\\.tif$"),
  reduce  = c(NA, NA, NA, NA, NA, "mean", NA),
  source  = "static",
  type    = "raster",
  band    = NA_integer_,
  # terra::extract accepts "simple" (nearest cell) or "bilinear" only.
  method  = c("bilinear", "simple", "simple", "bilinear", "bilinear",
              "bilinear", "bilinear"),
  # Confirmed at probe on 2026-07-31, not assumed. No file in this set declares
  # an NA flag, so every sentinel has to be stated here or it enters matching
  # as data.
  #   elevation  SRTM voids are -32768, the INT2S minimum. None appeared in the
  #              probe sample, but the guard costs nothing and the sampled
  #              minimum of -84 m is real, Death Valley reaching -86 m.
  #   aspect     flat cells are coded -1, which is not a sentinel and is not
  #              set to NA. See ca_derive_rules. Sampled range is -1 to 359.88,
  #              so the alternative 360 coding is not in use.
  #   travel     65535 is the documented sea and no-data sentinel.
  nodata  = c(-32768, NA_real_, NA_real_, NA_real_, NA_real_,
              NA_real_, 65535),
  scale   = NA_real_,
  active  = TRUE,
  units   = c("m", "degrees", "degrees from north",
              "mm/yr", "degrees C", "persons/km2", "minutes"),
  purpose = c("matching covariate",
              "matching covariate",
              "matching covariate, circular, decomposed at stage 4",
              "matching covariate, PRISM 1991-2020 normal",
              "matching covariate, PRISM 1991-2020 normal",
              "matching covariate, WorldPop 2000-2002 mean",
              "matching covariate, settlement class 11"),
  stringsAsFactors = FALSE
)

# Exact-matching strata. Joined by point in polygon, not extracted, and written
# as character rather than float32. Writing a HUC8 or an EPA code as a float
# would silently drop leading zeros, turning one exact stratum into another.
#
# Confirmed at probe on 2026-07-31. ca_eco_l3_prj.shp carries 13 features with
# US_L3CODE as unpadded strings of width 1 or 2, so codes read 1, 4, 5, 13, 14.
# They are padded to width 2 for a key of fixed width. Padding changes no
# grouping, since character comparison of unpadded codes is already
# unambiguous, but it prevents a width-dependent bug when the key is pasted
# into a compound identifier at stages 8 and 9.
#
# WBDHU12_CA_proj.shp carries 4,473 features with huc12 present at width 12
# throughout, so the truncation to 8 is exact.
ca_zone_layers <- data.frame(
  layer   = c("ecoregion_l3", "huc8"),
  subdir  = c("ecoregion", "huc"),
  file    = c("ca_eco_l3_prj.shp", "WBDHU12_CA_proj.shp"),
  field   = c("US_L3CODE", "huc12"),
  derive  = c("pad_2", "substr_8"),
  source  = "static",
  type    = "vector",
  active  = TRUE,
  units   = c("EPA Level III code", "WBD HUC8 code"),
  purpose = c("exact matching stratum, EPA Level III",
              "exact matching stratum, derived from HUC12"),
  stringsAsFactors = FALSE
)

# Static layers have no year dimension, so ca_extract_path() does not apply.
ca_static_path <- function(layer, lulc) {
  dir <- ca_work("extract", "static", layer)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("%s_%s.parquet", layer, lulc))
}

# Vector joins build a geometry object per chunk, which raster extraction does
# not, so they run at a tenth of the raster chunk size.
ca_zone_chunk <- 2e6L

# Stage 3c output. Long format, one row per pixel per intersecting polygon, so
# a pixel with three recorded entries yields three rows. This is required by
# ca_thin_rules$single_event_only, which cannot be evaluated from a
# first-hit join. Row counts therefore exceed the pixel count and the file is
# not comparable to a stage 3d file.
ca_vector_path <- function(layer, lulc) {
  dir <- ca_work("extract", "vector", layer)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("%s_%s.parquet", layer, lulc))
}

# Stage 3c holds a chunk of points and every polygon match for that chunk, so
# it runs smaller than the stage 3d zone join.
ca_vector_chunk <- 5e5L

# Aspect in degrees is a circular quantity and cannot enter matching directly.
# 350 and 10 degrees are 340 apart numerically and 20 apart physically, which
# breaks both the standardised mean difference and the caliper at the wrap.
# Stage 3d writes native degrees. Stage 4 replaces aspect with its two
# orthogonal components, which are bounded, continuous, and behave correctly
# under any distance metric.
#   northness = cos(aspect * pi / 180)
#   eastness  = sin(aspect * pi / 180)
#
# The SRTM derivative codes flat cells as -1, confirmed at probe. A flat cell
# has no bearing, so it is neither missing nor -1 degrees. Applying the formula
# to -1 would return northness 0.9998 and enter every flat pixel into matching
# as due north. Both components are set to 0 instead, which is the correct
# limit and the centroid of all bearings, and the pixel is retained. Treating
# -1 as a sentinel would instead drop flat terrain from the matched sample
# non-randomly, which is worse.
#
# This departs from the submitted specification and is stated as a
# methodological correction in the response letter.
ca_derive_rules <- data.frame(
  source_layer = c("aspect", "aspect"),
  derived      = c("northness", "eastness"),
  rule         = c("cos_deg", "sin_deg"),
  flat_code    = c(-1, -1),
  flat_value   = c(0, 0),
  note         = c("cosine of aspect, 1 is due north, 0 where flat",
                   "sine of aspect, 1 is due east, 0 where flat"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 0.6  OUTCOME REGISTRY
# ---------------------------------------------------------------------------
# Single authority for outcome names, sources, and admissible transforms.
# Stage 9 joins outcomes to labels through this table, not by file order.
# NEP and NBP are signed and are therefore non-log only.
#
# NPP tree, NEP, and NBP are extracted from the 2022 NCSDA pipeline regardless.
# Stage 3 checks whether v2026.1 carries equivalents. If it does, source is
# promoted to almanac2026 and the 2022 series becomes cross-comparison. If it
# does not, ncsda2022 is primary and the year range shortens to 2021.

# layer is the extraction layer name, which is what ca_extract_path() is keyed
# on. Without it every stage 4 and stage 9 caller has to hardcode the mapping
# from outcome to layer, which is exactly the positional coupling this rewrite
# removes. Added 2026-08-01.
#
# resolved is TRUE for all seven as of 2026-08-01. NPP, NEP, and NBP stay on
# ncsda2022 for this revision. v2026.1 carries no equivalent under its carbon
# theme and no confirmation has come back from the dataset producers, so the
# scaffold rows in ca_almanac_layers stay active = FALSE and the upgrade column
# records the path without activating it. Promoting them later costs one edit
# here plus a stage 3a rerun for three layers, because nothing downstream
# branches on availability.
ca_outcomes <- data.frame(
  outcome   = c("agb_almanac", "gpp", "npp_tree", "nep", "nbp",
                "agb_emapr", "agb_lemma"),
  label     = c("Aboveground biomass (Almanac)", "GPP", "NPP (tree)",
                "NEP", "NBP", "Aboveground biomass (eMapR)",
                "Aboveground biomass (LEMMA)"),
  source    = c("almanac2026", "almanac2026", "ncsda2022", "ncsda2022",
                "ncsda2022", "emapr", "lemma"),
  layer     = c("Carbon_AGB", "Carbon_GPP", "NPP", "NEP", "NBP",
                "eMapR", "LEMMA"),
  upgrade   = c(NA, NA, "almanac2026", "almanac2026", "almanac2026", NA, NA),
  carbon_fr = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  signed    = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
  # npp_tree changed from TRUE to FALSE, 2026-08-01. NPP is legitimately
  # negative where respiration exceeds photosynthesis, verified at 1,784 to
  # 29,508 pixels per year reaching -130 gC m-2 yr-1, so log() returns NaN and
  # R drops those observations without comment. Clipping would erase real
  # signal rather than error. The existing rule that signed variables are
  # non-log only therefore covers NPP, and it is now applied. This departs from
  # the submitted specification, which reported NPP on a log scale, and is
  # stated as a correction in the response letter.
  log       = c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
  nonlog    = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  scope     = c("all", "all", "all", "all", "fire", "all", "all"),
  role      = c("main", "main", "main", "main", "si", "main", "main"),
  resolved  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Resolve an outcome to the source and layer that back it. One place, so a
# source promotion is a data edit rather than a code search.
ca_outcome_layer <- function(outcome) {
  i <- match(outcome, ca_outcomes$outcome)
  if (is.na(i)) stop("Unknown outcome: ", outcome)
  list(source = ca_outcomes$source[i], layer = ca_outcomes$layer[i])
}

# Per-outcome year windows for the stage 9 panel. eMapR ends 2017 and LEMMA
# ends 2016 while Almanac runs to 2025, so the panel is ragged by construction.
ca_outcome_years <- function(outcome) {
  src <- ca_outcomes$source[match(outcome, ca_outcomes$outcome)]
  row <- ca_sources[ca_sources$source == src, ]
  seq(max(row$first_year, min(ca_years$analysis)), row$last_year)
}

# Variables entering matching in pre-treatment mean form. Canopy height comes
# from Almanac Fire_LCP band 6 and is a covariate only, never an outcome.
# Named by the matching covariate they produce, valued by what they are read
# from. Canopy height has no ca_outcomes row because it is never an outcome, so
# it is given its source and layer directly.
ca_pretrt_vars <- data.frame(
  covariate = c("gpp_pre5", "agb_pre5", "ch_pre5"),
  source    = c("almanac2026", "almanac2026", "almanac2026"),
  layer     = c("Carbon_GPP", "Carbon_AGB", "Fire_LCP_CH"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 0.7  GROUP CODES
# ---------------------------------------------------------------------------
# Codes follow the submitted pipeline. CPAD protection levels 1 to 6 collapse
# to a single P suffix at stage 6.

# CPAD managing agency level. This is where land.cpad 1 to 6 came from in the
# submitted pipeline, 5_data_recode_management.R L107. It is agency level, not
# GAP status. All six are treated as protected areas and collapse to a single P
# suffix at stage 6, as in the submitted version.
# CPAD agency levels serve two different purposes, and conflating them was the
# error in an earlier draft of this file.
#
#   pa_treatment      enters the protected-area treatment groups UP, F2P, F3P,
#                     F4P, T1P, T2P, T3P. Federal and State only, matching the
#                     submitted analysis: 6_3_group_selection.R L70 keeps only
#                     the P1 and P2 suffixes, and 8_1_matching_UP_UNP_NN.R L52
#                     subsets the treated set to c("UP1", "UP2"). Levels 3 to 6
#                     were categorised but never matched.
#
#   protected_other   protected, but not a public protected forest in the sense
#                     the paper uses. Never treatment. Still removed from the
#                     private non-protected control pool, because a county park
#                     or a land-trust preserve is not business as usual.
#
#   tribal            handled by the tribal exclusion mask.
#
# So the treatment definition is unchanged from the submitted analysis, while
# the control definition is tightened. Private and Home Owners Association were
# dropped entirely in the submitted pipeline, which left them eligible as
# controls. Joint is public but was absent from the submitted case_match, so it
# sits in protected_other rather than being added to the treatment group.
# Confirmed against probe_3c_fields_20260731.csv. The classifying field is
# MNG_AG_LEV, a character label with 10 distinct values matching the labels
# below. The numeric codes 1 to 7 carried here previously were an artefact of
# the submitted pipeline's recode and did not exist in the shapefile.
#
# Every level is extracted at stage 3c and the treatment subset is taken at
# stage 6. Extraction stores what the source says, and role is an analytical
# decision that belongs where the groups are formed, not where the data are
# read. That also means a change of PA definition costs a stage 6 rerun rather
# than a stage 3c rerun.
ca_cpad_agency <- data.frame(
  label = c("Federal", "State", "County", "City", "Non Profit",
            "Special District", "Joint", "Private",
            "Home Owners Association", "Tribal"),
  role  = c("pa_treatment", "pa_treatment",
            "protected_other", "protected_other", "protected_other",
            "protected_other", "protected_other", "protected_other",
            "protected_other", "tribal"),
  stringsAsFactors = FALSE
)

ca_cpad_field <- "MNG_AG_LEV"

# Treatment group membership.
ca_cpad_pa_treatment <- ca_cpad_agency$label[
  ca_cpad_agency$role == "pa_treatment"]

# Everything that disqualifies a pixel from the private control pool.
ca_cpad_protected_all <- ca_cpad_agency$label[
  ca_cpad_agency$role != "tribal"]

# Multi-Source Land Ownership. This layer defines what is not private, so the
# private non-protected control is the complement of it and CPAD together.
#   public        public ownership, removed from the private control pool
#   exclude_both  conservation non-profits. Privately held but not managed for
#                 production, so neither a public protected area nor a
#                 business-as-usual private control
#   tribal        handled by the tribal exclusion mask
#
# Neither role contributes a treatment group. This layer only defines what is
# not private, so public and exclude_both are both removed from the control
# pool and neither enters a P group. That is what the submitted pipeline did by
# excluding every polygon in this layer from the private class.
ca_ownership_level <- data.frame(
  label = c("Federal", "State", "County", "City", "Special District",
            "Non Profit", "Tribal"),
  role  = c("public", "public", "public", "public", "public",
            "exclude_both", "tribal"),
  stringsAsFactors = FALSE
)

ca_ownership_public <- ca_ownership_level$label[
  ca_ownership_level$role == "public"]

# Everything in this layer that disqualifies a pixel from the private control
# pool, which is every level except tribal, handled by its own mask.
ca_ownership_nonprivate <- ca_ownership_level$label[
  ca_ownership_level$role != "tribal"]

# Confirmed against probe_3c_fields_20260731.csv. Own_Level carries exactly the
# seven labels above with no missing values across 57,404 polygons.
ca_ownership_field <- "Own_Level"

# MTBS thematic burn severity. The rasters carry six classes, not three, and
# the same three-tier logic used for thinning applies.
#   2, 3, 4        severity classes, the wildfire treatment arms
#   1, 5           inside a fire perimeter but not an analysable severity.
#                  Class 1 is unburned to low, class 5 increased greenness.
#                  Never treatment, and never an undisturbed control either,
#                  because the pixel did burn
#   6              non-processing mask, read as NA
ca_fire_severity <- data.frame(
  code  = 1:6,
  label = c("Unburned to low", "Low", "Moderate", "High",
            "Increased greenness", "Non-processing mask"),
  group = c(NA, "F2", "F3", "F4", NA, NA),
  role  = c("disturbance_only", "treatment", "treatment", "treatment",
            "disturbance_only", "nodata"),
  stringsAsFactors = FALSE
)

ca_fire_treatment_codes <- ca_fire_severity$code[
  ca_fire_severity$role == "treatment"]
ca_fire_disturb_codes <- ca_fire_severity$code[
  ca_fire_severity$role == "disturbance_only"]
ca_fire_nodata_code <- ca_fire_severity$code[
  ca_fire_severity$role == "nodata"]

# Thinning intensity. T1 Low, T2 Medium, T3 High, confirmed against the
# case_match block in 5_data_recode_management.R L252.
ca_thin_intensity <- data.frame(
  code  = 1:3,
  label = c("Low", "Medium", "High"),
  group = c("T1", "T2", "T3"),
  stringsAsFactors = FALSE
)

# Variable and Unknown intensity labels are preserved verbatim and never
# collapsed to Low, Medium, or High. Pixels carrying them are excluded from
# both the treatment and the control arm at stage 6.
#
# Label-based as of 2026-08-01. The numeric codes c(12, 99, 912, 923, 9123)
# were the submitted pipeline's case_match targets, where 9123 meant a record
# carrying all of Variable, Low, Medium, and High. Matching on those numbers
# requires reconstructing the recode that produced them, which is the coupling
# this rewrite removes. Anything not in ca_thin_intensity$label is excluded,
# so the list is a statement of what has been seen rather than a filter that
# must be complete.
# Observed on the Evergreen partition, 2026-08-01, by section 1 of
# 5_verify_recode.R. The earlier list was a guess and contained two labels that
# do not exist, VariableLM and VariableMH, while missing three that do. Nothing
# leaked, because ca_thin_is_eligible() is a whitelist on Low, Medium, and
# High rather than a blacklist. This list is documentation of what has been
# seen, not the mechanism.
ca_thin_excluded_labels <- c("N/A", "Unknown", "Variable",
                             "VariableHM", "VariableHML", "VariableML")

ca_thin_is_eligible <- function(intensity) {
  intensity %in% ca_thin_intensity$label
}

# Wildfire counterpart to ca_thin_rules. Added 2026-08-01, because nothing
# stated the multi-event rule for fire and the omission would have let reburned
# pixels enter F2, F3, and F4 with a contaminated post-event trajectory.
#
# WHY SINGLE EVENT ONLY, SYMMETRICALLY WITH THINNING
#
# Durability is defined in the manuscript as the post-event carbon trajectory
# tracked over three decades. A pixel that burns again inside that window has
# no such trajectory to track. Censoring at the second fire keeps the sample
# but is worse than excluding: the long-horizon end of the durability curve
# would then be estimated on the subset that happened not to reburn, while the
# short-horizon end uses everything, so the comparison group composition drifts
# across exactly the axis the central claim runs along.
#
# The cost is a scope condition rather than a bias. California reburn over 1985
# to 2025 is common, so the arms lose sample and the estimand becomes durability
# in singly burned stands. That belongs in the Discussion as a stated scope, and
# 5_verify_recode.R reports the reburn share so it carries a number.
#
# WHY A BURNED PIXEL IS NEVER A CONTROL, IN ANY YEAR
#
# The same reasoning that fixed multi_event_as_control for thinning. A pixel
# selected for future treatment is not business as usual, and its pre-event
# years are not a clean baseline for someone else's cohort. The ETWFE
# literature admits not-yet-treated units as comparisons, and this design
# declines that for two reasons specific to it: matching fixes the control set
# before the DiD runs, and the durability window requires a control whose own
# trajectory stays clean for thirty years. A not-yet-treated control has a
# treatment landing inside that window by construction.
#
# A treated pixel's OWN pre-event years remain its pre-period. That is what
# identifies the effect and is not affected by any of the above.
# CARB compliance offsets. Treatment is absorbing from the Reporting Period 1
# start year to the end of the record.
#
# The submitted pipeline ended treatment at max(end.year) from the issuance
# workbook. That treats absence of issuance as termination, which it is not.
# Reporting lags, verification scheduling, and the shift from annual to less
# frequent reporting after the first period all produce gaps indistinguishable
# from a stopped project. CARB forest protocols carry hundred-year permanence
# obligations, so the land is under commitment whether or not credits issued in
# a given year.
#
# last_year and n_reporting_periods stay in the event table as provenance and
# as an SI sensitivity. They are not used in the main specification.
#
# The cost, stated in the Methods. A terminated project or an unreversed
# reversal would be misclassified as treated. CARB records those separately and
# they are rare, and the resulting bias runs toward zero, which works against
# the paper's own finding rather than for it.
ca_offset_rules <- list(
  absorbing = TRUE,
  cohort_field = "start_year",
  cohort_basis = "Reporting Period 1 start date",
  join_key = "arb_id",
  end_year_source = NA,          # deliberately unused, see above
  end_year_sensitivity = "last_year"
)

# Cohort years outside the analysis window. Observed ranges on Evergreen are
# 1921 to 2026 for thinning and 1858 to 2026 for CPAD, so both tails are real.
#
#   after last analysis year   a planned entry with no observable outcome. It
#                              cannot be a cohort. It is still evidence that
#                              work is scheduled, not that work occurred, so it
#                              does NOT disqualify the pixel as a control
#   before first analysis year an event that happened, with no pre-period
#                              inside the panel. It cannot be a cohort, and it
#                              DOES disqualify the pixel as a control, because
#                              the stand was entered
#
# The asymmetry is the point. A future record describes intent and a past
# record describes history, and only one of them changes what the pixel is.
# PROTECTED AREAS ARE A STATE, NOT AN EVENT
#
# The manuscript estimates PA effects in calendar year, not in normalised event
# time, so the PA arm does not need a cohort and the pre-1990 establishment
# dates are not a problem to solve. Treatment is "protected throughout the
# analysis window", and the estimate is a per-calendar-year difference against
# matched unprotected controls. The other three arms are event-time and do need
# cohorts. This asymmetry is deliberate and is stated in the Methods.
#
# Three consequences, from the Evergreen counts measured 2026-08-01.
#
#   dated before 1990    58,047,938 rows. Treated in every analysis year
#   dated 1990 to 2025    2,341,637 rows. Treated from the establishment year
#                        onward, untreated before it. Time-varying treatment
#                        status is the normal case under a calendar-year
#                        specification, not an ambiguity. These pixels are NOT
#                        eligible as controls in their pre-establishment years,
#                        for the same reason not-yet-treated pixels are barred
#                        in the thinning and fire arms
#   undated               1,526,559 rows. Protection status is known, only the
#                        date is missing. Excluded from the treatment arm
#                        because status by year is undefined without a start
#                        date, and excluded from the control pool because the
#                        pixel is protected
#
# The undated records are retained through stage 5 precisely so that stage 6
# can exclude those pixels from the control pool. Discarding the records would
# make protected pixels indistinguishable from unprotected ones and would leak
# them into the business-as-usual baseline, which is the opposite of the
# intended effect.
ca_pa_rules <- list(
  specification = "calendar year, not event time",

  # Treatment is time-varying: protected_it = 1 where est_year <= t. An earlier
  # draft required protection across the whole window and dropped the 2.34
  # million mid-window pixels. That was wrong. A PA established in 2004 is
  # cleanly untreated through 2003 and cleanly treated from 2004, which is a
  # staggered adoption rather than an ambiguous case, and a calendar-year
  # specification accommodates it directly.
  treatment_requires_full_window = FALSE,
  treatment_time_varying = TRUE,
  mid_window_as_treatment = TRUE,

  # Unchanged, and for a reason unrelated to the above. A pixel selected for
  # future protection is not business as usual in its pre-establishment years,
  # and matching fixes the control set before estimation runs. Same rule as
  # ca_thin_rules and ca_fire_rules.
  mid_window_excludes_control = TRUE,

  # THE PROTECTION LABEL MUST HOLD AT THE COHORT YEAR
  #
  # Follows directly from treatment_time_varying. If protection is a
  # time-varying state, then a pixel that burned in 1995 inside a PA
  # established in 2010 was not in a PA when it burned, and calling it F2P
  # describes its status today rather than at the event.
  #
  # Applies to the P groups of the fire and thinning arms only. UP_UNP has no
  # cohort year, so the rule cannot touch it, which is what separates this from
  # the establishment-window restriction rejected above. Those pixels leave the
  # analysis rather than moving to the NP arm, because a pixel now inside a PA
  # is not business-as-usual private forest either and its post-event recovery
  # is shaped by protection that arrived later.
  #
  # Roughly 1.8 percent of protected pixels were established in-window, so this
  # is a consistency rule rather than a material filter. The count lands in the
  # exclusion ledger under pa_after_cohort.
  #
  # The symmetric condition on the control side, that a UP control be protected
  # before its matched treated unit's cohort year, depends on the pair and can
  # only be applied at matching. It is in ca_pending() as a stage 8 item.
  label_at_cohort_year = TRUE,

  undated_as_treatment = FALSE,
  undated_excludes_control = TRUE,

  # The undated records are retained through stage 5 precisely so stage 6 can
  # exclude those pixels. Dropping the records at 5a would make 1.5 million
  # protected pixels indistinguishable from unprotected forest and leak them
  # into the business-as-usual baseline.
  undated_records_retained = TRUE
)

ca_window_rules <- list(
  first_year = min(ca_years$analysis),
  last_year  = max(ca_years$analysis),
  future_record_as_cohort = FALSE,
  future_record_excludes_control = FALSE,
  pre_window_record_as_cohort = FALSE,
  pre_window_record_excludes_control = TRUE
)

ca_fire_rules <- list(

  single_event_only = TRUE,
  multi_event_as_control = FALSE,

  # Classes 1 and 5 are inside a perimeter but are not analysable severities.
  # They cannot be treatment and they disqualify the pixel as a control, the
  # same shape as an undated thinning record.
  unclassified_severity_as_treatment = FALSE,
  unclassified_severity_excludes_control = TRUE,

  # Class 6 is the non-processing mask, not evidence of anything.
  nodata_excludes_control = FALSE,

  severity_crosswalk = "MTBS thematic classes, see ca_fire_severity"
)

ca_thin_rules <- list(

  # A pixel with one recorded entry of known intensity and no other disturbance
  # is treatment. A pixel with two or more entries is neither treatment nor
  # control. A pixel with no recorded entry of any kind is control.
  single_event_only = TRUE,

  # Corrected 2026-08-01. The earlier value was TRUE, on the reasoning that a
  # multi-event pixel is undisturbed until its first entry and can serve as a
  # control until then. That contradicts the manuscript's own control
  # definition, which is undisturbed non-protected private forest. A stand
  # harvested twice is not business as usual in any window, and admitting it to
  # the control pool would put real disturbance into the counterfactual and
  # understate the thinning effect. Multi-event pixels leave the analysis.
  #
  # This narrows the estimand rather than biasing it. The ATT is identified for
  # singly entered stands, not for all thinned forest, which is a scope
  # condition for the Discussion and not a robustness problem. The within-layer
  # exclusion is 24.7 percent of Evergreen treated pixels. The cross-layer
  # figure is larger and is computed in section 8g of 3_verify_extract.R.
  multi_event_as_control = FALSE,

  # A record with no usable completion date cannot be assigned a cohort, so it
  # cannot enter the treatment arm. It is still evidence that something
  # happened on that pixel, so it disqualifies the pixel as a control. Any
  # other handling either invents a cohort year or discards real evidence of
  # disturbance. Affects 23.9 percent of Evergreen thin_facts_hfr rows and 14.2
  # percent of thin_facts_th.
  undated_record_treatment = FALSE,
  undated_record_excludes_control = TRUE,

  # Intensity unknown is not impact absent. A Variable record cannot be
  # assigned to Low, Medium, or High, so it is ineligible as treatment, and it
  # disqualifies the pixel as a control for the same reason an undated record
  # does.
  variable_intensity_as_treatment = FALSE,
  variable_intensity_excludes_control = TRUE,

  # Records classified non_disturbing stay in the control pool by default, and
  # that default is tested rather than assumed. Three range activities are
  # classified non_disturbing but alter standing vegetation by definition:
  # Range Cover Type Conversion, Range Cover Manipulation, and Range Forage
  # Improvement. Three more are indirect markers of operations rather than
  # impacts: Range Piling Slash, Cover brush pile for burning, and Leave Trees
  # (Wildlife Reasons), since slash and retention designations only exist where
  # cutting occurred.
  #
  # The test, at stage 6. For pixels whose only stage 3c record is
  # non_disturbing, compare Almanac Disturbance_TreeFrac and Disturbance_AGB
  # flag rates in the window around the record year against pixels with no
  # record, matched on ecoregion and ownership. A category flagging at an
  # elevated rate is disturbing whatever the label says and moves to
  # disturbance_only.
  # RESULT, Evergreen, 2026-08-02. Only three categories were separately
  # testable, because the test requires pixels whose sole record is
  # non_disturbing and the remaining categories never occur alone. All three
  # were retained.
  #
  #   Administrative Changes      90 px, obs 13.33 pct, exp 0.77 pct, ratio 17.4
  #   Range Control Vegetation    50,738 px, obs 1.24 pct, exp 0.62 pct, ratio 2.00
  #   Range Cover Manipulation    44,492 px, obs 0.17 pct, exp 0.66 pct, ratio 0.25
  #
  # Two calls sit on the decision boundary and neither changes the analysis.
  # Administrative Changes clears the ratio bar by a wide margin and fails the
  # n_pixels floor at 90 against 100, and 90 pixels are immaterial. Range
  # Control Vegetation lands at ratio 1.9953 against a cutoff of 2.
  #
  # The reason neither matters is that the screen sits downstream of this
  # decision. A pixel classified non_disturbing enters the control pool and is
  # then removed by screen_detection if it carries any detection in any year,
  # so the elevated-rate pixels are already excluded on their own evidence.
  # The classification only governs pixels that never flagged, which are
  # indistinguishable from undisturbed forest by every available measure.
  #
  # Where it does bite is pixels carrying a non_disturbing record alongside a
  # thinning or fire record, roughly 100,000 rows on Evergreen. There Tn
  # against Td decides between a valid treatment cohort and a multi_event
  # exclusion, and those are exactly the pixels the test cannot isolate. That
  # is what the sensitivity below is for.
  non_disturbing_as_control = TRUE,
  non_disturbing_verified = TRUE,

  # SI sensitivity. Treats every non_disturbing record as disturbance_only,
  # which is the conservative bound on the untestable categories. Read by
  # 6a_group_codes.R.
  non_disturbing_as_disturbing = FALSE,

  # record_class == "wildfire". 7,656,433 rows on Evergreen, 15.5 percent of
  # the thinning event table. These are FACTS records whose Knight activity
  # resolves to wildfire, so they are fire evidence arriving through a
  # management layer.
  #
  # They cannot be thinning treatment. Admitting them would attribute fire
  # effects to thinning, which inverts the manuscript's central contrast
  # between an NbS and a disturbance.
  #
  # They cannot be control either. The pixel burned, and discarding that is the
  # same error as ignoring an undated record.
  #
  # They also cannot establish a fire cohort. FACTS carries no severity, and
  # F2, F3, and F4 are severity-stratified arms, so a fire event with no
  # severity has no arm to enter. MTBS remains the sole source of fire cohorts.
  # A pixel carrying a FACTS wildfire record but no MTBS perimeter is therefore
  # excluded from every arm, which is the correct treatment of an event that is
  # known to have happened and cannot be classified.
  wildfire_record_as_treatment = FALSE,
  wildfire_record_excludes_control = TRUE,
  wildfire_record_establishes_cohort = FALSE,

  # disturbance_only behaves the same way and always has. Stated here so the
  # four record_class values are each accounted for in one place.
  #   treatment_LMH     22,098,444 rows, eligible for the thinning arm
  #   disturbance_only  19,522,624 rows, never treatment, excludes control
  #   wildfire           7,656,433 rows, see above
  #   non_disturbing       195,276 rows, control by default, tested at stage 6
  disturbance_only_as_treatment = FALSE,
  disturbance_only_excludes_control = TRUE,

  intensity_crosswalk = "Knight et al. 2022, Tables S4 and S5"
)

# ---------------------------------------------------------------------------
# 0.7b  CLASS CODE TAXONOMY
# ---------------------------------------------------------------------------
# Two products come out of stage 6 and they answer different questions.
#
#   class_code      complete. Every pixel receives one, there are no missing
#                   values, and repeat entry, reburn, and thinning-then-fire
#                   histories all appear as ordered token chains. Nothing in
#                   this manuscript is estimated on them. They exist because
#                   the data already distinguish these histories and
#                   discarding the distinction would cost a full rerun to
#                   recover for a follow-up study
#
#   analysis_group  the 16 codes this manuscript estimates on, plus NA
#
# Encoding an exclusion as a missing group code, which an earlier draft did,
# throws away exactly the information the taxonomy is for.

# EVENT TOKENS. Chain order is by year, then thinning before fire within a
# year, which is arbitrary and stated rather than left implicit.
ca_class_tokens <- data.frame(
  token = c("U",
            "T1", "T2", "T3", "Tv", "Td", "Tn", "Tu",
            "F2", "F3", "F4", "Fu"),
  arm   = c("none",
            "thin", "thin", "thin", "thin", "thin", "thin", "thin",
            "fire", "fire", "fire", "fire"),
  label = c(
    "No recorded event",
    "Thinning, Low intensity",
    "Thinning, Medium intensity",
    "Thinning, High intensity",
    "Thinning, intensity not resolvable",
    "Entry recorded, not a thinning treatment",
    "Recorded activity classified non-disturbing",
    "Thinning record with no usable date",
    "Wildfire, Low severity",
    "Wildfire, Moderate severity",
    "Wildfire, High severity",
    "Fire evidence with no analysable severity"),
  source = c(
    "absence of any record",
    "record_class treatment_LMH, Knight intensity Low",
    "record_class treatment_LMH, Knight intensity Medium",
    "record_class treatment_LMH, Knight intensity High",
    "Variable, VariableHM, VariableHML, VariableML, Unknown, N/A, or two eligible records disagreeing in one pixel-year",
    "record_class disturbance_only",
    "record_class non_disturbing, subject to ca_nondisturbing_classes()",
    "event_year NA, unordered, written last in the chain",
    "MTBS thematic class 2",
    "MTBS thematic class 3",
    "MTBS thematic class 4",
    "MTBS class 1 or 5, or a FACTS record whose Knight activity resolves to wildfire"),
  # counts_as_event drives eligibility. Tu has no date so it cannot be ordered,
  # and it disqualifies both arms through its own rule rather than through the
  # event count. Tn does not count, subject to the probe.
  counts_as_event = c(FALSE,
                      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
                      TRUE, TRUE, TRUE, TRUE),
  treatment_eligible = c(FALSE,
                         TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                         TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

ca_tokens_treatment <- ca_class_tokens$token[ca_class_tokens$treatment_eligible]

# RECORDS ARE NOT EVENTS
#
# Both source layers are activity-level. FACTS writes one row per activity, so
# a single stand entry produces a harvest row, a site preparation row, and a
# slash row. MTBS writes one row per perimeter year, and a fire also recorded
# through a FACTS wildfire activity line appears twice. The first stage 6 run
# counted rows as events and lost 18,646,561 pixels to multi_event, 6,253,412
# of which carried at most one measured severity and a duplicate Fu.
#
# Records at the same pixel, in the same year, in the same arm are one event.
# Records in different years are different events, because a stand entered in
# 2005 and again in 2012 was entered twice. Same-year collapse is the
# defensible floor and needs no tuning parameter. Whether it should widen is a
# question about how far apart duplicate rows actually sit, so 6a_group_codes.R
# prints the gap distribution between consecutive same-arm entries rather than
# guessing at a window.
#
# ca_token_info ranks tokens by how much they say about magnitude, lowest wins.
# A measured severity beats Fu, which carries none. A measured intensity beats
# Td, which records that the stand was entered without saying how hard. Tv only
# survives where no measured intensity is present, so two eligible records
# disagreeing in one year resolve to the higher intensity rather than to Tv.
#
# Fu arrives through the thinning layer but belongs to the fire arm, so it does
# not compete with thinning tokens. A thinning record and a wildfire record in
# one year are one entry and one fire, which is two events and is meant to be.
ca_token_info <- c(F4 = 1L, F3 = 2L, F2 = 3L, Fu = 4L,
                   T3 = 1L, T2 = 2L, T1 = 3L, Tv = 4L, Td = 5L, Tn = 6L)

# Tokens carrying no magnitude. Fu says a fire touched the pixel without a
# usable severity, Td says the stand was entered without saying how hard. Only
# these merge into a measured neighbour one year away, which is what the gap
# diagnostic supports. Two measured tokens never merge, at any gap.
ca_tokens_no_magnitude <- c("Fu", "Td")

# STATUS. One per pixel, evaluated in this order, first match wins.
ca_status_levels <- data.frame(
  status = c("TR", "P", "PO", "PUB", "NGO", "UN", "NP"),
  label  = c("Tribal land",
             "CPAD Federal or State",
             "CPAD protected, other agency level",
             "Public ownership outside CPAD",
             "Conservation non-profit outside CPAD",
             "Record present, label not declared in config",
             "Private non-protected"),
  role   = c("excluded", "treatment", "excluded", "excluded", "excluded",
             "excluded", "control"),
  source = c("tribal layer, or Own_Level role tribal",
             "MNG_AG_LEV in ca_cpad_pa_treatment",
             "MNG_AG_LEV, any other non-tribal level",
             "Own_Level role public",
             "Own_Level role exclude_both",
             "MNG_AG_LEV or Own_Level value absent from the tables above",
             "complement of everything above"),
  stringsAsFactors = FALSE
)

# UN exists because a silent default is worse than an excluded pixel. A value
# absent from ca_cpad_agency or ca_ownership_level makes match() return NA,
# every role test then fails, and the pixel falls through to NP and into the
# control pool without being examined. The stage 6 verifier caught 5,807 UNP
# pixels carrying a CPAD record that way, and a further 36,618 tribal-owned
# pixels reached the control pool because only the tribal layer was consulted
# and not the ownership layer's own tribal role. 6a_group_codes.R logs every
# undeclared label with its record count, so UN should be empty once the tables
# are complete. It is not a category to live with.

# ANALYSIS GROUPS. 16 codes. OP is carried because the Methods restrict offsets
# to non-protected private land and the group is the audit of that claim. If it
# is non-empty the count belongs in the response letter.
ca_analysis_groups <- data.frame(
  group   = c("UP", "UNP", "OP", "ONP",
              "F2P", "F3P", "F4P", "F2NP", "F3NP", "F4NP",
              "T1P", "T2P", "T3P", "T1NP", "T2NP", "T3NP"),
  arm     = c("pa", "control", "offset", "offset",
              "fire", "fire", "fire", "fire", "fire", "fire",
              "thin", "thin", "thin", "thin", "thin", "thin"),
  role    = c("treatment", "control", "treatment", "treatment",
              "treatment", "treatment", "treatment",
              "control", "control", "control",
              "treatment", "treatment", "treatment",
              "control", "control", "control"),
  status  = c("P", "NP", "P", "NP",
              "P", "P", "P", "NP", "NP", "NP",
              "P", "P", "P", "NP", "NP", "NP"),
  token   = c("U", "U", "U", "U",
              "F2", "F3", "F4", "F2", "F3", "F4",
              "T1", "T2", "T3", "T1", "T2", "T3"),
  stratum = c(NA, NA, NA, NA,
              "Low severity", "Moderate severity", "High severity",
              "Low severity", "Moderate severity", "High severity",
              "Low intensity", "Medium intensity", "High intensity",
              "Low intensity", "Medium intensity", "High intensity"),
  analysed = c(TRUE, TRUE, FALSE, TRUE,
               rep(TRUE, 12)),
  stringsAsFactors = FALSE
)

# EXCLUSION RULES, applied in this order. The first rule that fires owns the
# pixel, so the ledger reads as a sequence of disjoint removals rather than
# overlapping counts.
ca_exclusion_order <- c(
  "tribal",
  "pa_undated",
  "status_not_analysed",
  "undated_thin_record",
  "multi_event",
  "unclassifiable_event",
  "pre_window_event",
  "pa_after_cohort",
  "offset_event_conflict",
  "offset_out_of_window",
  "no_detection_at_event",
  "screen_detection"
)

# Tribal lands are an exclusion mask, not a treatment group.
ca_tribal <- list(exclude = TRUE, group_code = "TrNP")

# ---------------------------------------------------------------------------
# 0.8  DISTURBANCE SCREEN
# ---------------------------------------------------------------------------
# The record-based canopy loss assertion is replaced by a detection-based
# screen on Almanac Disturbance_TreeFrac. Control pixels must show zero
# detected disturbance across all available years. The submitted code used
# thrd <- 0 while Methods L225 to L226 state greater than 5 percent canopy
# loss. The discrepancy is unresolved and is listed in ca_pending().

# Resolved 2026-08-01. Neither of the two candidate values survives scrutiny.
#
#   thrd = 0, the submitted code. Any non-zero detection disqualifies a pixel,
#   so a pixel survives only by detecting exactly zero in all 39 years.
#   Retention of a truly undisturbed pixel is roughly (1 - p)^39 in the
#   per-year false positive rate p, which is 46 percent at p = 0.02 and 14
#   percent at p = 0.05. The loss is not random. Change detection false-positive
#   rates rise on steep, heterogeneous, and high-biomass canopy, so the
#   surviving control pool is systematically flatter, more homogeneous, and
#   lower in biomass than the treated pixels it is the baseline for.
#
#   5 percent, the Methods text. A round number with no relationship to this
#   product's error structure. In a closed conifer stand 5 percent tree fraction
#   is a real canopy opening. MTBS misses fires below roughly 400 ha in the
#   West, so this screen is the only thing catching them, and a generous
#   threshold admits real disturbance to the business-as-usual baseline and
#   attenuates every treatment effect toward zero.
#
# The threshold is therefore derived from the data rather than asserted, by
# section 4 of 4_verify_clean.R, and written to meta as screen_threshold.csv.
# ca_screen_threshold() reads that artifact. The threshold is data, not code,
# so it cannot be changed by editing a script.
#
# Annual criterion only. A cumulative bound across the window was considered
# and rejected. Repeated entry is already caught twice, by the stage 3c record
# and by ca_thin_rules$multi_event_as_control = FALSE, and slow background
# decline from drought and beetle mortality is a common trend that the DiD
# differences out and that matching on ecoregion, HUC8, and climate normals
# already conditions on. A cumulative bound would shrink the control pool for
# a signal handled elsewhere. Stated in the Methods so it is a decision rather
# than an omission.

ca_screen <- list(
  layer            = "Disturbance_TreeFrac",
  years            = 1986:2024,
  control_rule     = "annual detected loss at or below the derived threshold in every year",
  threshold_file   = ca_meta("screen_threshold.csv"),
  threshold_grid   = c(0, 0.01, 0.05, 0.10),   # SI sensitivity, 0.02 dropped

  # Disturbance_AGB is extracted and stored at the detected pixel-years by
  # stage 5b, and it is not a criterion anywhere. Two criteria would force a
  # tie-break rule for the disagreement cases and that rule would be arbitrary.
  # Its one use is an SI sentence giving the agreement rate between the two
  # products at detected pixel-years.
  secondary_layer  = "Disturbance_AGB",
  secondary_role   = "reporting only, never a criterion",

  cumulative_rule  = NULL,              # rejected, see the note in 0.8

  # TREATED PIXELS MUST SHOW A DETECTION
  #
  # A polygon record asserts that work was scheduled or completed inside a
  # boundary, not that it touched every 30 m pixel inside that boundary.
  # Unit-level records over-assign treatment to pixels, and requiring an
  # observed canopy change at the pixel is the correction. This is the
  # submitted specification, 6_2_categorization.R L58 required
  # pct.chng.tree < thrd for every treated group.
  #
  # The cost, stated in the Methods. Detection sensitivity falls with event
  # magnitude, so the requirement removes proportionally more of the low
  # severity and low intensity arms and selects the survivors toward the
  # stronger end within their own class. That compresses the contrast between
  # strata rather than inflating it, so it works against the paper's own
  # gradient finding rather than for it. Retention per group is reported by
  # 6_verify_groups.R and quoted in the Methods.
  require_treatment_detection = TRUE,

  # Archival completion dates should sit close to the detection. Fire
  # attribution can fall a year or more away, and measured sensitivity against
  # MTBS class 2 is 0.33 at plus or minus one, so the fire window is the
  # parameter most worth testing.
  thin_detection_window = 1L,
  fire_detection_window = 1L,
  detection_window_grid = c(0L, 1L, 2L, 3L)    # SI sensitivity
)

# The derived annual threshold, as a fraction of tree cover. Errors rather than
# falling back to a default, because a silent default is how an unresolved
# methodological question becomes a published number.
ca_screen_threshold <- function(which = "primary") {
  f <- ca_screen$threshold_file
  if (!file.exists(f)) {
    stop("Screen threshold not derived yet. Run:\n",
         "  Rscript 4_verify_clean.R threshold\n",
         "Expected artifact: ", f)
  }
  tb <- utils::read.csv(f, stringsAsFactors = FALSE)
  i <- match(which, tb$name)
  if (is.na(i)) stop("No threshold named '", which, "' in ", f)
  tb$value[i]
}

ca_nondist_file <- ca_meta("nondisturbing_decision.csv")

ca_nondisturbing_classes <- function() {
  if (!file.exists(ca_nondist_file)) {
    stop("Non-disturbing decision not derived yet. Run:\n",
         "  Rscript 6a_group_codes.R probe\n",
         "Expected artifact: ", ca_nondist_file)
  }
  tb <- utils::read.csv(ca_nondist_file, stringsAsFactors = FALSE)
  tb$activity[tb$decision %in% "disturbing"]
}

# ---------------------------------------------------------------------------
# 0.9  PRE-TREATMENT WINDOWS
# ---------------------------------------------------------------------------
# Computed per cohort at stage 8 from the full-year panels stored at stage 4.5.
# The window must itself be disturbance-free. Where a detected disturbance
# falls inside it, the window is truncated, subject to the 3-year minimum.
#
# Protected areas use three tiers because many predate the carbon record.
# Offsets are all post-2012 and use a fixed 2007 to 2011 window.
# Controls use the same window as their matched treatment pixel.

ca_pretrt <- list(
  # WINDOW LENGTH AND THE HALO DROP.
  #
  # For arms whose recorded event year lags the canopy signal in Landsat, year
  # c-1 is contaminated and the window drops it. The anchor stays at c-5 rather
  # than sliding to c-6, so the window is four years (c-5 to c-2) instead of
  # five, and it never reaches earlier than the submitted rule or runs off the
  # 1985 start of the record.
  #
  # EVIDENCE. The stage 8c diagnostic counted every COLD/CCDC detection on the
  # 60,357 capped TNP_UNP treated pixels by its gap from the FACTS cohort year.
  #
  #   gap   +6    +5    +4    +3    +2    +1     0    -1    -2
  #   pct  0.4   0.5   0.8   1.4   1.1  68.9  20.6   6.1   0.2
  #
  # 95.6 percent of detections sit at +1, 0, and -1, and the largest bin is +1,
  # the year BEFORE the recorded year, at 68.9 percent. From +2 outward the
  # distribution is flat and +3 carries more than +2, so there is no decaying
  # tail. This is a three-year detection halo around the treatment. It appears
  # at 65 to 70 percent in all three intensity groups.
  #
  # Year c-1 therefore contains the treatment rather than preceding it. A
  # covariate computed over a window that includes c-1 is contaminated by the
  # thing it is meant to condition on, and matching on it would remove the
  # variation the DiD is trying to detect. Dropping c-1 clears it.
  #
  # PER ARM. MTBS records the actual fire year and CARB records a reporting
  # period start, neither an administrative completion date. The treated-side
  # attrition under the unlagged rule confirms it without a second diagnostic
  # run. Thinning lost 74 percent of its treated pixels to in-window detections
  # while fire lost 5.2 percent, offsets 0.05 percent, and protection 0.6
  # percent. Had fire carried a +1 halo like thinning, its loss would have been
  # of the same order. It is not, so fire detections sit at gap 0 and no drop
  # is warranted. CPAD establishment is a legal date on a parcel with no canopy
  # detection involved.
  #
  # event_drop = 1L means "drop one year before the event", giving a window of
  # window_n minus event_drop = 4 years for thinning, 5 for everything else.
  window_n = 5L,
  window_min = 3L,
  event_drop = c(pa = 0L, fire = 0L, thin = 1L, offset = 0L),

  # NO TRUNCATION ON IN-WINDOW DETECTIONS.
  #
  # The earlier code censored the pre-treatment window at the last COLD/CCDC
  # detection inside it and returned NA below window_min years. That rule
  # was measuring the two sides differently, since controls carry no detections
  # by construction and always receive a clean five-year mean while treated
  # pixels receive a truncated or absent value. The asymmetry biased the
  # augmented arm's covariate and discarded three quarters of the thinning
  # treatment group.
  #
  # With c-1 dropped, the remaining in-window detections at gaps +2 through +5
  # are at background rate (0.5 to 1.4 percent of all detections). The 374
  # units with a genuine prior detection showed a value shift of -0.9 SD in
  # AGB, which is real contamination on those units, but 374 out of 60,357 is
  # 0.6 percent, so the effect on the population mean is 0.006 times 0.9 or
  # roughly 0.005 SD, below the resolution of any matching config in the sweep.
  #
  # Genuine prior recorded events are already excluded at stage 6 through the
  # multi-event screen at a 32.26 percent exclusion rate.
  truncate_on_detection = FALSE,

  pa_tiers = data.frame(
    tier = c("post_1990", "1985_1990", "pre_1985"),
    rule = c("true pre-establishment 5-year mean",
             "available pre-establishment years, 3-year minimum",
             "1985 to 1989 baseline-period mean"),
    stringsAsFactors = FALSE
  ),
  pa_sensitivity = "post-1990 establishment subset run separately",
  # CPAD YR_EST uses 0 as an unknown marker, not as a year. The submitted
  # pipeline collapsed it with genuine pre-1986 establishment through
  # ifelse(YR_EST < 1986, 1985, YR_EST), which silently asserts that every
  # unknown-year PA is old. It is kept as its own category here, assigned the
  # same 1985 cohort as the submitted version so the main result is comparable,
  # and dropped in a sensitivity run so the assertion can be tested rather than
  # assumed.
  #
  # Stage 3c writes the unknown as NA, not as 0, because 0 is a null token in
  # the same sense as the CAL FIRE 1899-12-30 placeholder and reading it as a
  # year is exactly the error being corrected. The stage 6 rule therefore keys
  # on NA. The 3c probe confirms the share at 73.5 percent of superunits, which
  # is large enough that the sensitivity run is not optional.
  cpad_year_unknown = NA_integer_,
  cpad_year_unknown_rule = "NA assigned 1985, dropped in sensitivity",
  cpad_year_floor = 1985L,
  # Offset cohorts run 2012 to 2019, not all from 2012, so a fixed 2007 to 2011
  # window would be pre-treatment only for the earliest projects. Offsets use
  # the same per-cohort rule as everything else: the five years before the
  # project's own start year. The fixed window in an earlier draft of this file
  # was carried from the assumption that every project began in 2012.
  offset_window = NULL,
  control_rule = "same window as the matched treatment pixel"
)

# ---------------------------------------------------------------------------
# 0.10  TREATMENT PAIRINGS
# ---------------------------------------------------------------------------
# SIX PAIRINGS. Every treated group is matched to undisturbed forest of its own
# ownership class.
#
# The submitted design matched burned PA against burned non-PA and thinned PA
# against thinned non-PA, 8_2 L57 and L69, 8_3 L57 and L69. But the estimation
# never contrasted those two arms in a regression. 10_2_..._sub_1_1_short.R L73
# subsets the matched data to treat==1 and runs ETWFE with cgroup = "notyet",
# and sub_1_0 L73 does the same for treat==0. Each ownership arm was estimated
# separately against its own not-yet-treated pixels and the ownership
# comparison was made between the two sets of ATTs. The matched counterpart was
# discarded before estimation.
#
# So cross-arm matching imposed the full cost of 1:1 without replacement while
# doing no work in the regression. Both arms collapse to the smaller side. The
# stage 6 feasibility section measured what that costs. Low-intensity thinning
# yields 317 matched pairs from 14,450 treated and 1,171 control units, so the
# stratum cannot be estimated at all.
#
# Matching each treated group to undisturbed forest of the same ownership
# removes the coupling. Control pools are UP at roughly 24 million and UNP at
# roughly 14 million, so no treated group is capped by control scarcity. It
# also supplies genuine never-treated controls, which lets every pairing use
# cgroup = "never" as UP_UNP already did at 10_1 L194 rather than leaning on
# not-yet-treated units with no never-treated anchor.
#
# Every pairing now estimates additionality against a business-as-usual
# baseline of undisturbed forest. The ownership contrast is recovered by
# comparing the P and NP pairings of the same arm, which is what the submitted
# analysis did with its two sets of estimates.
ca_pairings <- data.frame(
  pairing    = c("UP_UNP", "ONP_UNP",
                 "FP_UP", "FNP_UNP",
                 "TP_UP", "TNP_UNP"),
  treat      = c("UP", "ONP",
                 "F2P,F3P,F4P", "F2NP,F3NP,F4NP",
                 "T1P,T2P,T3P", "T1NP,T2NP,T3NP"),
  control    = c("UNP", "UNP",
                 "UP", "UNP",
                 "UP", "UNP"),
  arm        = c("pa", "offset", "fire", "fire", "thin", "thin"),
  ownership  = c(NA, "NP", "P", "NP", "P", "NP"),

  # NO SUBSET ON pa_est_year, IN EITHER ROLE.
  #
  # An earlier draft of this block restricted UP to establishment in 1990 or
  # later on the grounds that a pixel protected before the panel opens has no
  # observed pre-period. That reasoning imports event time into an arm that
  # does not use it, and it contradicts ca_pa_rules$specification.
  #
  # The PA arm is specified in calendar year. Treatment is the time-varying
  # indicator protected_it = 1 where est_year <= t, so a PA established in 1955
  # is treated in every year of the panel and its effect is identified against
  # never-protected controls through the matched level difference and common
  # time effects, not through a within-unit break. Dropping those pixels would
  # discard the ongoing effect of every PA established before the panel opens,
  # which is 98.2 percent of protected forest on Evergreen and the great
  # majority of California's protected carbon.
  #
  # The submitted analysis made the same choice deliberately. The restriction
  # to establishment after 1985 sits in 10_1_DiD_analysis_UP_UNP_all_1.R at
  # L66 to L69 and is commented out. The filter that is active, at L77, keeps
  # subclasses whose treated member carries a valid establishment year in both
  # CPAD sources, which is a data-quality filter and not a date window.
  #
  # Consequence for the fire and thinning arms. UP as a control pool needs no
  # establishment filter either, for the same reason. What those arms do need
  # is that the protection label hold at the cohort year, which is a separate
  # question recorded in ca_pending().

  # The control is undisturbed and carries no severity or intensity, so the
  # stratum cannot enter the exact set as it did in the submitted design. Each
  # pairing runs once and the strata are recovered by subsetting the matched
  # data at stage 10, which is where 10_2 L63 and 10_3 L62 already do it.
  strata     = c(NA, NA,
                 "Low severity,Moderate severity,High severity",
                 "Low severity,Moderate severity,High severity",
                 "Low intensity,Medium intensity,High intensity",
                 "Low intensity,Medium intensity,High intensity"),
  baseline_type = "business as usual",
  cgroup     = "never",
  staggered  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Exact-matching set. Three variables for every pairing now, since the stratum
# has left the exact set along with cross-arm matching.
ca_match_exact <- function(pairing) {
  if (!pairing %in% ca_pairings$pairing) stop("Unknown pairing: ", pairing)
  ca_match$exact
}

# Treated and control group codes for a pairing.
ca_pairing_groups <- function(pairing, side = c("treat", "control")) {
  side <- match.arg(side)
  i <- match(pairing, ca_pairings$pairing)
  if (is.na(i)) stop("Unknown pairing: ", pairing)
  trimws(strsplit(ca_pairings[[side]][i], ",")[[1]])
}

# UP is a treated group in UP_UNP and a control pool in FP_UP and TP_UP. Those
# are different roles for the same pixels in different pairings, which is
# legitimate, but it means the 100,000 cap cannot be applied at stage 7. A
# group cannot be both capped as a treatment and left whole as a control pool
# in one sampled file. The cap therefore moves to stage 8, where it is applied
# per pairing to the treated side only.

# ---------------------------------------------------------------------------
# 0.11  MATCHING SPECIFICATION
# ---------------------------------------------------------------------------
# THREE covariate arms run the full sweep and all three are reported. Baseline
# is the submitted specification, 8 continuous plus 3 exact. Augmented adds the
# three pre-treatment means, 11 continuous plus 3 exact. augmented_noch drops
# canopy height for 10 continuous plus 3 exact. All three run simultaneously
# rather than augmented_noch being held back as a conditional fallback, because
# the collinearity between ch_pre5 and agb_pre5 is the question and running the
# comparison is the only way to answer it.
#
# Aspect in degrees is replaced by northness and eastness, see ca_derive_rules.
# Baseline is therefore 8 continuous rather than the submitted 7, which is a
# correction to a circular covariate rather than a change of covariate set.
#
# COHORT STRATIFICATION
#
# ca_pretrt$control_rule sets the control window as the same window as the
# matched treatment pixel. That is circular unless matching runs cohort by
# cohort, because the augmented arms need the control covariate value before
# the match exists. The alternative, one fixed window for every control,
# balances a 2018 treated unit against a control measured in 1985 to 1989 and
# reports balance that does not exist, since carbon stocks trend over the
# record. So the stratification is the only implementation faithful to the
# locked rule rather than a competing option chosen for convenience.
#
# Applied to all three arms. Applying it to the augmented arms only would
# confound the covariate set with the matching structure, and the arm contrast
# is the reason all three run.
#
# NO CROSS-COHORT BOOKKEEPING. Each matchit() call is self-contained and 1:1
# without replacement holds inside that call. A control pixel drawn by the 1995
# cohort remains available to the 2015 cohort, as it does across pairings. The
# duplication rate is counted by the stage 8 verifier and reported rather than
# prevented, because with control pools in the tens of millions against treated
# sides capped in the hundreds of thousands the collision rate is expected to
# be small and an exhaustion rule would be a cost paid against a problem that
# has not been measured.
#
# THE MATCHABLE POOL
#
# Exact matching on lulc, ecoregion_l3, and huc8 means a control in a cell that
# holds no treated unit for this cohort can never be matched. Those rows are
# dropped before the distance model is fitted. The reduction is exact for the
# match itself. It does change the fitted propensity score, since the glm sees
# the matchable pool rather than the whole pool, and that is stated in the
# Methods rather than left implicit. Without it the distance model for every
# cohort is fitted on 14 to 24 million rows and bart and elasticnet do not run
# at all.

# EVERY LOGGED OUTCOME TAKES log1p, NOT log.
#
# An earlier version of this function reserved log1p for GPP, on the
# grounds that GPP is clipped at zero by ca_clean_rules while the biomass
# products are strictly positive. The stage 9 panel shows they are not.
# Zero aboveground biomass appears in every pairing, at 33,016 pixel-years
# in FP_UP, 26,462 in FNP_UNP, and 766 in UP_UNP, and eMapR carries zeros
# too. Those are not missing values. They are stand-replacing fire, which
# is the observation this manuscript is most interested in.
#
# log() returns -Inf there and R drops the row without a warning, so the
# plain log would silently delete the post-fire floor from the fire arms
# and bias the ATT toward zero. log1p costs nothing at the scale of these
# values, since the outcome is gC m-2 and the median sits above 6,000, so
# log1p and log differ by under 0.02 percent everywhere the two are both
# defined.
#
# This also restores the submitted specification, which applied log1p to
# all four outcomes, 10_1_DiD_analysis_UP_UNP_all_1.R L28 to L36.
ca_log_transform <- function(outcome) "log1p"

ca_match <- list(
  exact = c("lulc", "ecoregion_l3", "huc8"),
  continuous_baseline = c(
    "elevation", "slope", "northness", "eastness", "ppt_normal",
    "tmean_normal", "pop_density", "city_travel_time"
  ),
  continuous_augmented = c(
    "elevation", "slope", "northness", "eastness", "ppt_normal",
    "tmean_normal", "pop_density", "city_travel_time",
    "gpp_pre5", "agb_pre5", "ch_pre5"
  ),
  continuous_augmented_noch = c(
    "elevation", "slope", "northness", "eastness", "ppt_normal",
    "tmean_normal", "pop_density", "city_travel_time",
    "gpp_pre5", "agb_pre5"
  ),
  arms_enabled = c("baseline", "augmented", "augmented_noch"),
  ratio = 1L,
  replace = FALSE,
  caliper_sd = 0.25,
  cluster = ~subclass,
  large_subclass_guard = 50000L,  # present in 10_1, absent in 10_2, now global

  # WHERE THE BALANCE NUMBERS COME FROM
  #
  # summary.matchit, not cobalt::bal.tab. The MS3 template reads bal.tab and
  # keeps Diff and KS, which is SMD and a Kolmogorov-Smirnov statistic. That
  # drops variance ratio and both eCDF statistics, and those are two of the
  # three balance criteria named at Methods L272 to L274 and two of the four
  # columns of Table S3. summary() returns sum.all and sum.matched carrying
  # means by group, standardized mean difference, variance ratio, eCDF mean,
  # and eCDF max, which is the set the submitted scripts used.
  metric_source = "summary.matchit",

  # Cohort stratification, see the block above.
  cohort_stratified = TRUE,
  cohort_replacement_across_cohorts = TRUE,
  matchable_pool_only = TRUE,

  # NO GATE ON MATCHED PAIRS.
  #
  # Every stratum is estimated and reported with its pair count. A hard floor
  # is a researcher degree of freedom and a wide interval is an honest result.
  # Stage 8 writes the counts and makes no decision.
  min_pairs_report = 0L,

  # THE PRIMARY ARM.
  #
  # augmented. The reviewer asked for pre-treatment covariates, so all three go
  # in. augmented_noch existed only to catch collinearity between ch_pre5 and
  # agb_pre5 degrading balance, and the sweep shows it does not. Including
  # canopy height moves the mean absolute SMD of the other ten covariates by at
  # most 0.003 in any pairing, and ch_pre5 balances on its own terms at |SMD|
  # 0.041 or better. The fallback is not needed and is not reported.
  primary_arm = "augmented",

  # SELECTION RULE, applied per pairing at stage 8b.
  #
  # THE RULE SUGGESTS, IT DOES NOT DECIDE. Stage 8b ranks the sweep and marks
  # what this rule would pick, then prints a snippet for ca_match$selected
  # below. Nothing downstream reads the ranking. Stage 9 reads ca_match$selected
  # only, which is written by hand after the ranking has been checked by eye.
  #
  # BALANCE DECIDES, RETENTION GATES. The point of matching is a credible
  # dynamic baseline, and balance is what measures it. Match rate is not a thing
  # to maximise, it is a thing to have enough of. A config retaining 90 percent
  # with every covariate under SMD 0.10 beats one retaining 97 percent with
  # three covariates over it, because the second has not produced the comparison
  # the paper claims to make.
  #
  # RETENTION IS MEASURED AGAINST THE STRUCTURAL CEILING, NOT AGAINST 1.
  #
  # A treated pixel goes unmatched for two unrelated reasons. Its exact cell
  # holds no available control, which is a property of the data and identical
  # for every config in the pairing. Or the config's caliper discarded it, which
  # is a property of the config. Only the second is a choice, so only the second
  # should be judged.
  #
  # nn_exact_1to1 separates them. It carries no caliper and no distance model,
  # so its match rate is exactly the structural ceiling, and every config's
  # shortfall below that ceiling is caliper-induced. The gate is therefore
  #
  #   retention_ratio = match_rate / match_rate(nn_exact_1to1)
  #
  # which reads as the share of the structurally matchable treated units the
  # caliper kept. An absolute floor cannot do this. A pairing whose cells only
  # support 60 percent matching would fail an absolute 0.80 floor on every
  # config including the uncalipered ones, and the rule would return nothing.
  #
  # WHERE 0.90 COMES FROM. Not from the data, and no longer from a gap.
  #
  # An earlier version of this comment put 0.90 inside an empty band, on the
  # grounds that every surviving config landed at 0.907 or above and every
  # failing one at 0.780 or below. That band was an artefact of stage 8b, which
  # was silently excluding every config that had errored on a single cohort.
  # The configs it excluded were the calipered elastic and BART runs, which are
  # exactly the ones that sit in the band. With the full sweep aggregated the
  # observed ratios run 0.788, 0.820, 0.862, 0.880, 0.898, 0.907 with no gap
  # wider than 0.042, and the previous justification does not survive.
  #
  # The threshold is therefore stated a priori and on the estimand rather than
  # derived from the observed distribution. A caliper that discards a tenth or
  # more of the structurally matchable treated units is no longer estimating
  # the ATT on the treated population. It is estimating it on the subset that
  # happened to have a close control, and that subset is selected on covariate
  # proximity, which is correlated with everything the matching is meant to
  # balance. 0.90 is the point at which the paper is willing to call the
  # matched sample representative of the treated group. It is a round number
  # chosen for a reason that does not depend on this dataset, which is the only
  # kind of threshold that can be defended once the data have been seen.
  #
  # The configs the gate excludes are named in the SI with their balance, so a
  # reader can see what was given up. In UP_UNP the gate is doing real work,
  # since bart_cal025 balances better at 0.820 retention than the selected
  # glm_cal025 does at 0.928. Elsewhere it excludes nothing that would have won.
  #
  #   1. gate    retention_ratio >= retention_min
  #   2. pass    every covariate under smd_max and inside the vr window
  #   3. choose  lowest mean absolute SMD among those that pass
  #   4. else    lowest maximum absolute SMD above the gate, flagged
  #
  # A pairing whose structural ceiling is itself low is a data limitation to
  # report in Table S3, not a selection problem. The ceiling is written to the
  # summary so it can be read directly.
  selection = list(
    ceiling_config = "nn_exact_1to1",
    retention_min  = 0.90,
    smd_max        = 0.10,
    vr_window      = c(0.5, 2.0),
    tie_break      = "mean_abs_smd",
    
    # THE MATCH RATE DENOMINATOR IS FIXED PER PAIRING AND ARM.
    #
    # match_rate is matched pairs over matchable treated, where matchable
    # treated is the capped treated draw after the arm's complete-case filter.
    # n_treat is set before any config is fitted and is written on failed rows
    # too, so the level-wise value summed over every level gives a denominator
    # that does not move when a config fails a cohort. An earlier version of
    # stage 8b summed n_treat over surviving rows only, which handed a failing
    # config a smaller denominator and so a higher rate, the opposite of what
    # the failure means.
    #
    # The spaced pool and the cap sit upstream and are not in the spec table.
    # They are constants per pairing and arm, they are in the stage 8 logs, and
    # the Methods state the cap. Nothing imputes them.
    rate_denominator = "matchable_treated",
    
    # A CONFIG THAT FAILS A COHORT IS COMPLETE, NOT PARTIAL.
    #
    # glmnet cannot fit a binomial model to a cohort holding one treated pixel,
    # and a tight caliper can match nothing in a cohort of six. Those are
    # properties of the cohort. The config still ran, still wrote its matched
    # parquet, and the treated units in the failed cohorts are lost to it. They
    # count against its rate and are reported in meta/match_coverage.csv. Only
    # a cohort with NO spec row at all, meaning the task is still running or
    # was level subset through CA_LEVELS, excludes a config from the ranking.
    failed_levels_count_against_rate = TRUE,
    absent_levels_exclude_config     = TRUE,
    
    # THE eCDF MAXIMUM IS AVERAGED ACROSS COHORTS, NOT MAXIMISED.
    #
    # The eCDF maximum between two samples of one unit each is 1 whenever the
    # values differ. Taking a maximum across cohorts of wildly different size
    # therefore returns the smallest cohort rather than the worst balance, and
    # it returned exactly 1.000 for every covariate in every pairing holding a
    # singleton cohort, including covariates at SMD 0.003. Both eCDF statistics
    # are now aggregated as treated-weighted averages, which is the same rule
    # the rest of the table uses. The raw maximum is retained beside them with
    # the cohort it came from and that cohort's pair count.
    ecdf_aggregation = "treated_weighted_mean"
  ),

  # THE CHOSEN SPECIFICATION, ONE PER PAIRING, WRITTEN BY HAND.
  #
  # Filled from the snippet stage 8b prints, after the ranking has been read.
  # One entry per pairing naming both the arm and the config, since the arm
  # decision and the config decision are made together. Table S3 reports the
  # primary arm only.
  #
  # Empty until the sweep is complete. Stage 9 stops rather than guessing.
  selected = list(
    UP_UNP  = list(arm = "augmented", config = "mahal_pscal"),
    ONP_UNP = list(arm = "augmented", config = "elastic_cal025"),
    FP_UP   = list(arm = "augmented", config = "mahal_nocal"),
    FNP_UNP = list(arm = "augmented", config = "elastic_cal025"),
    TP_UP   = list(arm = "augmented", config = "mahal_nocal"),
    TNP_UNP = list(arm = "augmented", config = "mahal_pscal")
  ),

  # Balance is summarised across cohorts as a treated-weighted average, which
  # is how the ATT itself aggregates at stage 10. Means and variances pool
  # exactly from the n and sd columns the balance table carries. eCDF
  # statistics do not pool, so they are reported as the weighted average and
  # the maximum across cohorts.
  summary_weight = "n_treat_matched",

  # ONE ARM IS REPORTED. Three arms run so the covariate set can be chosen on
  # evidence, not so three sets of results are published. Table S3 and Figure
  # S6 carry the primary arm alone, and the arms that lost are named in the
  # Methods without their tables.
  report_arms = "primary_only",

  # Plots are not written at stage 8. match_balance carries everything a love
  # plot draws, and 6 pairings by 3 arms by 11 configs by 4 plot types is 792
  # files produced before anything has been selected.
  write_plots = FALSE
)

# The chosen specification for a pairing. Stops rather than falling back to a
# default, because a silent default here would put an unexamined matched dataset
# into the DiD panel and nothing downstream would notice.
ca_match_selected <- function(pairing) {
  s <- ca_match$selected[[pairing]]
  if (is.null(s)) {
    stop("No selected specification for ", pairing, ". Run 8b_aggregate_match.R, ",
         "read the ranking, and write the choice into ca_match$selected.")
  }
  s
}

# Continuous covariate vector for an arm.
ca_match_covariates <- function(arm) {
  key <- switch(arm,
                baseline = "continuous_baseline",
                augmented = "continuous_augmented",
                augmented_noch = "continuous_augmented_noch",
                stop("Unknown matching arm: ", arm))
  ca_match[[key]]
}

# Which arms need the pre-treatment reads. baseline does not, so its tasks skip
# the cohort carbon reads entirely and still stratify, keeping the matching
# structure identical across arms.
ca_match_arm_needs_pretrt <- function(arm) {
  any(ca_pretrt_vars$covariate %in% ca_match_covariates(arm))
}

# The stage 8 array. Pairing crossed with arm, 6 by 3.
ca_match_tasks <- function() {
  g <- expand.grid(pairing = ca_pairings$pairing,
                   arm = ca_match$arms_enabled,
                   stringsAsFactors = FALSE)
  g <- g[order(match(g$arm, ca_match$arms_enabled),
               match(g$pairing, ca_pairings$pairing)), ]
  rownames(g) <- NULL
  g
}

# The variable that defines a stratification level.
#
# UP_UNP is the awkward case and it is not an inconsistency. ca_pa_rules
# specifies the PA arm in calendar year, so UP has no cohort in the event-time
# sense. What it has is a pre-treatment window assigned by establishment date
# through ca_pretrt$pa_tiers, and the window is what the stratification needs.
# So UP_UNP strata are windows keyed on pa_est_year and every other pairing
# keys on cohort_year.
ca_match_cohort_var <- function(pairing) {
  if (identical(pairing, "UP_UNP")) "pa_est_year" else "cohort_year"
}

# Pre-treatment window for one unit, returned as a year vector.
#
# The rule lives here rather than in the stage 8 script because it is a
# decision, not an implementation detail, and because the stage 8 verifier has
# to test the same rule the script applied.
#
#   fire, thin, offset   the five years before the unit's own cohort year
#   pa, est >= 1990      true pre-establishment five-year mean
#   pa, est 1985 to 1989 available pre-establishment years, three-year minimum,
#                        falling back to the baseline period where fewer than
#                        three years are available
#   pa, est < 1985 or NA the 1985 to 1989 baseline period
#
# NA establishment is CPAD YR_EST 0, kept as its own category and assigned the
# baseline window here so the main result stays comparable with the submitted
# version, then dropped in the sensitivity run named by
# ca_pretrt$cpad_year_unknown_rule.
ca_pretrt_window <- function(pairing, year) {
  base <- ca_years$baseline
  n <- ca_pretrt$window_n
  if (identical(pairing, "UP_UNP")) {
    if (is.na(year) || year <= min(base)) return(base)
    if (year > max(base) + 1L) return(seq(year - n, year - 1L))
    yrs <- seq(min(base), year - 1L)
    if (length(yrs) < ca_pretrt$window_min) return(base)
    return(yrs)
  }
  if (is.na(year)) return(integer(0))
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  drop <- ca_pretrt$event_drop[[arm]]
  if (is.null(drop) || is.na(drop)) drop <- 0L
  end   <- year - 1L - drop          # c-1 for fire/offset/pa, c-2 for thin
  start <- year - n                   # c-5 always, anchor does not move
  yrs <- seq(start, end)
  yrs <- yrs[yrs >= min(ca_years$study)]
  if (length(yrs) < ca_pretrt$window_min) return(integer(0))
  yrs
}

# Full 11-config sweep from the MS3 template. Best config selected per pairing
# per arm on SMD, VR, eCDF, and match retention.
ca_match_configs <- list(
  list(id = "nn_exact_1to1",  distance = "glm",                caliper = FALSE, mahvars = FALSE),
  list(id = "glm_nocal",      distance = "glm",                caliper = FALSE, mahvars = TRUE),
  list(id = "glm_cal025",     distance = "glm",                caliper = TRUE,  mahvars = FALSE),
  list(id = "elastic_nocal",  distance = "elasticnet",         caliper = FALSE, mahvars = TRUE),
  list(id = "elastic_cal025", distance = "elasticnet",         caliper = TRUE,  mahvars = TRUE),
  list(id = "bart_nocal",     distance = "bart",               caliper = FALSE, mahvars = TRUE),
  list(id = "bart_cal025",    distance = "bart",               caliper = TRUE,  mahvars = TRUE),
  list(id = "mahal_nocal",    distance = "mahalanobis",        caliper = FALSE, mahvars = FALSE),
  list(id = "mahal_pscal",    distance = "glm",                caliper = TRUE,  mahvars = TRUE),
  list(id = "rmahal_nocal",   distance = "robust_mahalanobis", caliper = FALSE, mahvars = FALSE),
  list(id = "rmahal_cal025",  distance = "robust_mahalanobis", caliper = TRUE,  mahvars = FALSE)
)

ca_balance_metrics <- c("smd", "vr", "ecdf_mean", "ecdf_max", "n_matched",
                        "retention_rate")

# ---------------------------------------------------------------------------
# 0.11b  DiD PANEL
# ---------------------------------------------------------------------------
# Stage 9 turns the selected matched dataset into the panel stage 10 estimates
# on. It reads ca_match$selected and nothing else from the sweep, so a
# specification that was never chosen cannot reach the results.
#
# THE DiD UNIT IS A MATCHED ROW, NOT A PIXEL
#
# Matching runs cohort by cohort and 1:1 without replacement holds inside a
# call, not across calls, so a control pixel drawn by the 1997 cohort stays
# available to the 2003 cohort. It then appears twice in the matched data under
# two subclasses. Keying the panel on pixel_id would give att_gt() and etwfe()
# duplicate unit-year rows and both would treat the second copy as new
# information.
#
# unit_id is therefore one integer per matched row, assigned at stage 9a and
# stable thereafter. pixel_id is carried beside it so the reuse is measurable.
# The reuse rate per pairing goes to meta/did_units_summary.csv, which is the
# number ca_pending() asks for in the cohort-stratified matching sentence.
#
# CLUSTERING IS ON THE PAIR
#
# vcov = ~ subclass in the submitted version, which is pair_id here. That is
# the correct level for a 1:1 matched design rather than a convention carried
# forward. Matching induces dependence between the two members of a pair, and
# treating them as independent observations understates the variance of the
# ATT. Clustering on the pixel instead would split every pair into two clusters
# and discard exactly that dependence, so it is not an alternative and is not
# reported as a sensitivity.
#
# The one case where the pair is not enough is a control pixel that entered two
# pairs, since its two rows are perfectly dependent and sit in different
# clusters. That is a cross-pair dependence neither pair nor pixel clustering
# captures alone, and it needs two-way clustering on pair and pixel. Whether it
# is worth the sentence is an empirical question and control_reuse in the
# units summary answers it. At a reuse ratio near 1 there is nothing to correct.

ca_did <- list(

  # PANEL YEARS. 1985 to 2025, the full extraction window, not the 1990 to 2025
  # analysis window.
  #
  # ca_years$analysis governs which events may establish a cohort, through
  # ca_window_rules, and that is settled at stages 5 and 6 and must not move.
  # The panel is a different question. Opening it to 1985 gives the 1990 to
  # 1994 cohorts pre-treatment periods they otherwise lack entirely, which is
  # what aggte(type = "dynamic") needs to test parallel trends, and it lets the
  # PA establishments dated 1986 to 1989 enter att_gt() as identified cohorts
  # rather than being dropped as already treated in the first period.
  #
  # eMapR and LEMMA start in 1990 and NBP in 1986 regardless, through
  # ca_layer_years(), so the extension reaches the Almanac and NCSDA series
  # only. The panel is ragged by construction either way.
  #
  # THE CIRCULARITY THIS CREATES, AND WHY IT IS NOT NEW. Units treated before
  # 1990 carry pre-treatment covariates computed over 1985 to 1989, so matching
  # has balanced those very years and an event study over them would show flat
  # pre-trends by construction. That is already true of every staggered cohort,
  # whose window sits inside the panel: a 1995 fire cohort is matched on 1990
  # to 1994. The property is inherent to matching on pre-treatment outcomes and
  # belongs in the Methods rather than in a shorter panel.
  panel_years = ca_years$study,

  # The unit key, and what is clustered on.
  unit_key      = c("pixel_id", "subclass"),
  cluster       = "pair_id",

  # WHERE PAIR CLUSTERING STOPS BEING SUFFICIENT.
  #
  # A control pixel drawn by two cohorts sits in two clusters and its two
  # rows are the same forest, which pair clustering cannot see. Measured
  # at stage 9a, the reuse ratio runs 1.0018 on UP_UNP to 1.0548 on
  # FP_UP, so the dependence reaches at most 5 percent of control rows
  # and vcov = ~ pair_id stands without qualification. The threshold is
  # here so a rerun that moves the number is caught rather than assumed
  # away, and section 9 of the verifier tests against it.
  cluster_reuse_threshold = 1.10,

  # Where each attribute is authoritative.
  cohort_source   = "stage 6 groups file",
  stratum_scope   = "subclass",
  control_gvar    = 0L,
  pa_gvar_source  = "pa_est_year",
  pa_unknown_year = ca_pretrt$cpad_year_floor,

  # THE ESTIMATOR PAIRING, CARRIED FORWARD FROM THE SUBMITTED SCRIPTS.
  #
  # Read from the files rather than remembered. The PA arm is specified in
  # calendar year and reports calendar-time ATTs, 10_1 L232. Every event-time
  # arm reports normalised event time, 10_5 L187, 10_2 L143, 10_3 L140. The
  # parallel trends test is aggte(type = "dynamic") in normalised event time
  # for all four arms, 10_1 L162, 10_2 L122, 10_3 L109, 10_5 L114.
  #
  # Stage 10 reads this table. Stage 9 writes gvar so both estimators can run
  # and takes no position on either.
  emfx_type = c(pa = "calendar", offset = "event",
                fire = "event", thin = "event"),
  attgt_aggregation = c("group", "dynamic"),

  # TWO PER-PANEL QUESTIONS AT STAGE 10. BOTH ARE SCOPED TO ONE OUTCOME.
  #
  # Each panel is its own file and stage 10 runs once per panel, so anything
  # either of these removes is removed from that outcome alone. A pair that
  # LEMMA cannot see is untouched in the Almanac, eMapR, NPP, NEP, and NBP
  # panels. Nothing here propagates across outcomes.
  #
  # complete_pairs_only. FALSE. When one member of a matched pair has no value
  # for this outcome, 9b writes the other member alone. The estimators pool
  # every unit carrying gvar = 0 into one control group rather than reading the
  # pair, so a half pair costs nothing mechanically and deleting its surviving
  # member would discard a real observation to enforce a symmetry the estimator
  # does not use. Row-level complete cases have already run at 9b, which is the
  # filter the estimation actually needs.
  #
  # The share is reported rather than acted on. It is under 0.5 percent for
  # every product except LEMMA, whose static mask breaks 2.7 to 13.6 percent of
  # pairs and reaches 49.7 percent of Mixed pairs in FP_UP. That belongs in the
  # LEMMA data note, not in a deletion rule.
  #
  # viable_cohorts_only. TRUE, and mechanical rather than a choice. A cohort
  # whose gvar falls after the last year of this outcome has no post-treatment
  # observation, so it yields no ATT, and carrying a non-zero gvar it is not in
  # the never-treated control group either. did::att_gt() errors on such a
  # group rather than skipping it. eMapR ends 2017 and LEMMA 2016 while fire
  # cohorts run to 2025, so the retained share falls to 42.6 and 34.0 percent
  # on FP_UP and 42.4 and 29.7 percent on FNP_UNP. The Almanac panels retain
  # 100 percent in every pairing and are unaffected.
  #
  # The retained share is per outcome and goes in Table S3, since it is what
  # the secondary biomass products can speak to rather than a sample loss.
  complete_pairs_only  = FALSE,
  viable_cohorts_only  = TRUE,

  # Rows whose outcome is NA are not written. Stage 10 would drop them at
  # complete.cases() anyway and LEMMA alone carries a 7.4 percent static mask,
  # so writing them costs storage and buys nothing. The share dropped per
  # pairing and outcome is recorded in the manifest.
  drop_missing_outcome = TRUE,

  # NO TRANSFORM AT STAGE 9. The panel carries the converted value once and
  # stage 10 applies ca_log_transform() where ca_outcomes$log is TRUE. Writing
  # both columns would double the panel to encode a decision that already lives
  # in the outcome registry.
  transform_stage = 10L,

  # THE MANIFEST REPLACES csv.list[i].
  #
  # The submitted pipeline labelled outcomes by the position of their file in
  # list.files(), 9_DiD_input_all.R L34 and L53, so adding a product silently
  # renamed every column after it. Stage 9b writes one manifest row per panel
  # naming the pairing, outcome, layer, arm, config, year range, and row count,
  # and stage 10 reads the manifest rather than a directory listing.
  manifest = "did_panel_manifest.csv"
)

# Stage 9b array. Pairing crossed with outcome, filtered by ca_outcomes$scope,
# so NBP appears for the two fire pairings only.
ca_did_tasks <- function() {
  out <- do.call(rbind, lapply(seq_len(nrow(ca_outcomes)), function(i) {
    sc <- ca_outcomes$scope[i]
    pr <- if (identical(sc, "all")) ca_pairings$pairing else
      ca_pairings$pairing[ca_pairings$arm == sc]
    if (!length(pr)) return(NULL)
    data.frame(pairing = pr, outcome = ca_outcomes$outcome[i],
               stringsAsFactors = FALSE)
  }))
  out <- out[order(match(out$pairing, ca_pairings$pairing),
                   match(out$outcome, ca_outcomes$outcome)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Years written for one panel. Keyed on the extraction layer rather than on the
# source, because NBP starts in 1986 where the rest of NCSDA starts in 1985 and
# ca_sources carries only the source-level bound.
ca_did_years <- function(outcome) {
  intersect(ca_layer_years(ca_outcome_layer(outcome)$layer), ca_did$panel_years)
}

# gvar, vectorised over the matched rows of one pairing. Controls are never
# treated and take 0. The rule lives here rather than in 9a because the PA
# asymmetry is a decision the verifier has to test against the same source the
# script applied.
ca_did_gvar <- function(pairing, treat, cohort_year, pa_est_year) {
  if (!pairing %in% ca_pairings$pairing) stop("Unknown pairing: ", pairing)
  g  <- rep(as.integer(ca_did$control_gvar), length(treat))
  tr <- which(treat == 1L)
  if (!length(tr)) return(g)
  if (identical(pairing, "UP_UNP")) {
    y <- as.integer(pa_est_year[tr])
    # ESTABLISHMENT BEFORE THE PANEL COLLAPSES TO ONE ALWAYS-TREATED COHORT.
    #
    # CPAD carries establishment back to 1889, and 98.3 percent of the UP
    # treated sample predates the panel. Left raw, that is roughly sixty
    # distinct cohorts which all mean "protected in every observed year"
    # and are indistinguishable inside the panel. Each one costs att_gt() a
    # dropped group and etwfe() a coefficient absorbed by its own group
    # fixed effect, for no gain, since nothing in data opening in 1985
    # separates 1889 from 1984.
    #
    # The floor also reproduces the submitted specification, which
    # collapsed them through ifelse(YR_EST < 1986, 1985, YR_EST). Unknown
    # establishment, CPAD YR_EST 0 arriving as NA, takes the same value and
    # stays separable through pa_year_known for the sensitivity run.
    y[is.na(y)] <- as.integer(ca_did$pa_unknown_year)
    y <- pmax(y, as.integer(min(ca_did$panel_years)))
    g[tr] <- y
  } else {
    g[tr] <- as.integer(cohort_year[tr])
  }
  g
}

# Admissible gvar range for a pairing, asserted at stage 9a. The PA arm reaches
# below the treatment window because establishment predating the panel is
# treated throughout, so its floor is the CPAD floor rather than 1990.
ca_did_gvar_range <- function(pairing) {
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  i   <- match(arm, ca_treat_years$arm)
  lo  <- if (identical(pairing, "UP_UNP")) ca_pretrt$cpad_year_floor else
    ca_treat_years$first[i]
  c(as.integer(lo), as.integer(ca_treat_years$last[i]))
}

# Which emfx aggregation a pairing reports, from the arm.
ca_did_emfx_type <- function(pairing) {
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  if (is.na(arm)) stop("Unknown pairing: ", pairing)
  unname(ca_did$emfx_type[[arm]])
}

# EVERY emfx AGGREGATION A PAIRING REPORTS. The PA arm reports two.
#
# calendar stays the PA headline. It is the only aggregation that uses the
# always-treated block floored at 1985, which holds the pre-1990 CPAD units
# and therefore most of the protected forest carbon, and that is why the arm
# was given a calendar specification in the first place.
#
# event was added on 2026-08-13 so protection carries the same pre-treatment
# panel as fire, thinning, and offsets, which the submitted etwfe version did
# not return at all. Its negative event times rest on the post-1985 cohorts
# alone, since the always-treated block has no pre-treatment period when the
# panel opens in 1985. That is the same subset did::att_gt() speaks to, which
# drops the always-treated group outright and reports 37 pre-treatment event
# times on 89,153 treated units across 31 cohorts. So the two estimators
# describe the same units in the pre-window, which is what makes the SI
# comparison a comparison. The figure caption has to say which units those
# are, because the post-treatment side pools all 31 cohorts.
#
# ONP_UNP is not given a calendar aggregation. Offset cohorts are all dated
# and the event specification answers the question that arm asks.
ca_did_emfx_types <- function(pairing) {
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  if (is.na(arm)) stop("Unknown pairing: ", pairing)
  if (identical(unname(arm), "pa")) {
    return(c("calendar", "event"))
  }
  ca_did_emfx_type(pairing)
}

# EVERY aggte AGGREGATION A PAIRING REPORTS, the att_gt counterpart of the
# function above.
#
# The SI comparison is only a comparison if the two estimators are asked for
# the same aggregations, so the PA arm gets calendar here for the same reason
# it gets calendar there, and no other arm does. Giving calendar to fire,
# thinning, and offsets would write rows no figure reads and would leave the
# att_gt set carrying an aggregation the reporting estimator does not have.
#
# What the two calendar aggregations describe is NOT the same population.
# did::att_gt() drops the always-treated group outright, so the att_gt
# calendar series rests on the post-1985 dated cohorts alone, while the etwfe
# calendar series uses the always-treated block floored at 1985 that holds the
# pre-1990 CPAD units and most of the protected forest carbon. That is the
# whole reason the arm was given a calendar specification. The SI caption has
# to say so, or the two Figure 2 panels read as one estimator disagreeing with
# another when they are answering questions about different units.
ca_did_aggte_types <- function(pairing) {
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  if (is.na(arm)) stop("Unknown pairing: ", pairing)
  if (identical(unname(arm), "pa")) {
    return(c(ca_did_est$aggte_types, "calendar"))
  }
  ca_did_est$aggte_types
}

# The array task table. For etwfe every run is its own task. For att_gt the
# ecoregion cells of a run ride inside their parent statewide task, so the
# array table is the statewide subset and ca_did_runs()$task points each cell
# at the task that computes it.
# Named ca_did_array() and not ca_did_tasks(), which is already the stage 9b
# pairing by outcome table at section 0.11a.
ca_did_array <- function(estimator = c("etwfe", "attgt")) {
  estimator <- match.arg(estimator)
  r <- ca_did_runs(estimator)
  if (identical(estimator, "attgt")) {
    r <- r[r$ecoregion == "", , drop = FALSE]
  }
  rownames(r) <- NULL
  r
}

# ---------------------------------------------------------------------------
# 0.11c  STAGE 10 ESTIMATION
# ---------------------------------------------------------------------------
# The submitted pipeline held roughly forty near-identical files differing only
# in pairing, log against non-log, ecoregion, and severity or intensity. Adding
# three outcomes and five years meant the same four edits in every one of them,
# and any file missed produced publishable-looking output computed on the old
# specification. Stage 10 is two scripts plus this table.
#
# THE SPLIT BETWEEN THE TWO SCRIPTS IS THE SUBMITTED ONE, AND IT IS BY COST.
#
#   10a_did_attgt.R   att_gt then aggte("simple") and aggte("dynamic"), plus
#                     aggte("calendar") on the PA arm alone, through
#                     ca_did_aggte_types(). 10_1 L112, 10_5 L73, 10_2 L85,
#                     10_3 L72.
#   10b_did_etwfe.R   etwfe then emfx.
#   10b_did_etwfe.R   etwfe then emfx. Statewide and per ecoregion.
#                     10_1 L188, 10_5 L143, and every 10_*_sub_2 file.
#
# BOTH ESTIMATORS NOW RUN AT ECOREGION LEVEL, BY DIFFERENT MEANS.
#
# The submitted scripts never ran att_gt inside an ecoregion file, and this
# pipeline followed them until 2026-08-16. The reason given was cost, and the
# reason was wrong. The heaviest att_gt cell fits in 15 to 20 minutes inside
# 100 GB, so a statewide cell and its six ecoregions together sit well under a
# day on acpu with cpu-normal. Withholding the ecoregion comparison to save
# time that was never at risk left the SI able to compare the two estimators
# statewide and nowhere else.
#
# What is kept from that reasoning is the array shape. One task per ecoregion
# would turn a 15 minute fit into a scheduling problem, so the ecoregions ride
# inside their statewide task, which also means one panel read serves all of
# them. See the task column in ca_did_runs() and ca_did_array().
#
# WHICH CELLS ARE RUN.
#
#   pa and offset arms      statewide, and one run per ecoregion
#   fire and thinning arms  one run per stratum, and one per stratum and
#                           ecoregion
#
# The same cell set for both estimators. They differ in how the cells are
# distributed across array tasks, not in which cells exist.
#
# Fire and thinning are always stratified, as submitted, where the severity and
# intensity loops sit outside every estimator call, 10_2 L61 and 10_3 L59.
# There is no severity-pooled estimate in the submitted results and none here.
#
# The ecoregions are not a fixed list. The submitted files hardcoded
# c(78, 1, 9, 4, 5, 8) for UP_UNP and c(78, 1, 4) for ONP_UNP, which is the
# same empirical fact written twice by hand. ca_did_runs() reads the levels
# actually present in the stage 9 units files, through meta/did_subset_levels.
#
# THE THINNING BASE PERIOD IS c-2, THROUGH A COHORT SHIFT.
#
# Stage 8c put 68.9 percent of canopy loss detections on thinned pixels at c-1,
# one year before the FACTS accomplishment date, which is administrative date
# offset rather than prior disturbance. ca_pretrt$event_drop["thin"] already
# drops c-1 from the matching window for that reason.
#
# If c-1 is also the base period then the reference year already contains the
# cut for most treated pixels, and the short-term stock reduction, which is the
# quantity this manuscript rests on, is estimated as a partly cut stand against
# a fully cut stand.
#
# Neither did::att_gt() nor etwfe::etwfe() exposes a per-cohort offset that
# both estimators read the same way, so the shift is applied to gvar itself, at
# estimation only. gvar_est = gvar - 1 on the thinning pairings. Both
# estimators then take the year before gvar_est as the base, which is c-2, and
# the stage 8 window c-5 to c-2 reads as four pre-treatment years ending at
# t-1 rather than a window with a hole in it.
#
# The shifted value is what stage 8 already wrote, since stage 8 sets the
# cohort to max(window) + 1 and the thinning window ends at c-2. Stage 9a
# restored it to c from the stage 6 groups file so the panel carries the
# administrative year, and this puts it back for estimation alone. The panel is
# untouched and stage 9 does not rerun.
#
# Consequence for stage 11. Thinning event time 0 is the detection year, one
# year before the FACTS accomplishment year. The axis label carries it and the
# Methods states it. Fire, offset, and protection are unshifted.
#
# WHAT THIS STAGE APPLIES THAT STAGE 9 DELIBERATELY DID NOT.
#
# viable_cohorts_only. A cohort whose gvar_est falls after the last year of
# this outcome has no post-treatment observation, so it yields no ATT, and
# carrying a non-zero gvar it is not in the never-treated group either.
# did::att_gt() errors on such a group rather than skipping it. The filter is
# per cell, since the year range is a property of the outcome and the cohorts
# present are a property of the subset.
#
# complete_pairs_only is FALSE and stays FALSE. The half-pair share is recorded
# per run and acted on nowhere.
#
# NO RUN-TIME SUBSAMPLE. 10_2 L69 drew 50,000 subclasses on the fire pairings
# and nothing else, so the fire estimates rested on a sample size rule no other
# arm used. The stage 8 cap of 100,000 treated units per analysis group
# replaces it.

ca_did_est <- list(

  # ESTIMATOR ARGUMENTS, CARRIED FROM THE SUBMITTED SCRIPTS.
  #
  # No covariates in either estimator. Matching balanced them, and the
  # submitted files that carried climate and terrain columns in col.list left
  # the xformla line commented out. 10_1 L116 and L189.
  xformla       = ~ 1,
  control_group = "nevertreated",
  cgroup        = "never",
  # The base set, for every arm. The PA arm adds calendar through
  # ca_did_aggte_types(), which is the aggregation Figure 2 is drawn on.
  #
  # simple replaces group. group returns one ATT per cohort, which no figure
  # reads, and simple returns the overall ATT, which the manuscript has never
  # carried for any arm. 10b does not run emfx type simple because it is a
  # second full pass through marginaleffects and cost 5 hours 48 minutes on
  # task 1217. aggte() is not that. The att_gt fit is already done and simple
  # is a reweighting of the group-time cells it holds, so the overall ATT is
  # nearly free here and expensive there. That asymmetry is deliberate.
  aggte_types   = c("simple", "dynamic"),

  # aggte() returns NA for the whole aggregation whenever any group-time cell
  # it averages is NA, so na.rm is required rather than optional. Submitted at
  # 10_1 L144 and L153. What the submitted pipeline left silent was how many
  # cells it removed, so that count goes to the run manifest.
  aggte_na_rm   = TRUE,

  # compress, not collapse. etwfe 0.6.2 supersedes the collapse argument
  # with compress, identical in effect. FALSE is the submitted setting,
  # 10_1 L239. TRUE pre-aggregates the panel before the marginal effects
  # and trades exact standard errors for speed, which is a change to the
  # reported inference rather than to the runtime alone.
  emfx_compress = FALSE,
  emfx_vcov     = TRUE,

  # att_gt() balances the panel internally when allow_unbalanced_panel is
  # FALSE, which is the submitted default and is kept. Stage 10 does the
  # balancing in the open first and reports the units it costs, so the number
  # is visible rather than inferred from a warning. etwfe needs no balance and
  # keeps every row, which is also what the submitted scripts did.
  attgt_balance          = TRUE,
  attgt_allow_unbalanced = FALSE,
  attgt_base_period      = "varying",
  attgt_bstrap           = TRUE,
  attgt_cband            = TRUE,

  seed = 160617L,

  # THE FLOOR ON TREATED UNITS PER CELL.
  #
  # meta/did_subset_levels.csv holds 189 cells and the ecoregion tail runs
  # down to a single treated unit, which no estimator can speak to. The
  # submitted analysis handled this by hardcoding six ecoregions for UP_UNP
  # and three for ONP_UNP, which is the same fact expressed as a literal.
  # The floor comes from the gap in the observed distribution and is set
  # once, here, after inspecting the levels table. Zero runs every cell.
  # Set to 48 on 2026-08-11 from the 189-cell distribution. Sorted on
  # treated units the ecoregion cells run 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3,
  # 4, 4, 4, 5, 5, 7, 9, 10, 10, 10, 12, 13, 16, 17, 17, 18, 20, 24, 25,
  # 27, 27, 30, then jump to 48, 51, 52, 52, 59, and climb without another
  # break. Thirty-three cells sit below the jump and none between 31 and
  # 47. The floor is the lower edge of the upper cluster, which is how the
  # matching config retention ratio was set against its own gap.
  #
  # The gap also lands where cluster-robust inference stops being
  # credible. Every standard error in this stage clusters on pair_id, and
  # a cell with 30 treated pairs has 30 treated clusters. The two
  # criteria, one from the data and one from the estimator, agree.
  #
  # It touches ecoregion cells only. The smallest statewide cell is
  # TNP_UNP Low intensity at 219 treated units, so 10a is unaffected.
  min_treated_units = 48L,

  # The c-2 rule, by arm.
  base_shift = c(pa = 0L, offset = 0L, fire = 0L, thin = 1L),

  # THE SUBMISSION TIER. Nothing analytical reads it.
  #
  # Tiered on rows times cohorts rather than on rows. emfx is the binding
  # cost in stage 10 and marginaleffects builds a jacobian with one row
  # per observation and one column per parameter, so a cell with twice
  # the rows and twice the cohorts is four times the peak, not twice.
  #
  # Calibrated on the probes. UP_UNP agb_almanac eco01 at 4.6 million
  # cost finished in 25 minutes under 1 GB of R heap. FNP_UNP agb_lemma
  # Moderate severity at 123 million cost was killed inside emfx event
  # after spending 3 hours 48 minutes on emfx simple. FP_UP agb_almanac
  # Low severity, the largest cell in the stage, is 287 million.
  tier_breaks = c(small = 2e7, medium = 8e7),

  levels_file = "did_subset_levels.csv",
  manifest    = c(etwfe = "did_etwfe_manifest.csv",
                  attgt = "did_attgt_manifest.csv")
)

# The base-period shift for a pairing, in years.
ca_did_base_shift <- function(pairing) {
  arm <- ca_pairings$arm[match(pairing, ca_pairings$pairing)]
  if (is.na(arm)) stop("Unknown pairing: ", pairing)
  as.integer(ca_did_est$base_shift[[arm]])
}

# gvar as the estimators see it. Controls stay at 0. Treated cohorts move by
# the arm shift, which is zero everywhere except thinning.
ca_did_gvar_est <- function(pairing, gvar) {
  s <- ca_did_base_shift(pairing)
  g <- as.integer(gvar)
  if (s) g[!is.na(g) & g > 0L] <- g[!is.na(g) & g > 0L] - s
  g
}

# Which transforms an outcome runs. Logged outcomes run both, since the
# submitted pipeline reported both and 11_5_1 and 11_5_2 read them separately.
# Signed outcomes run non-log only, which is a property of the registry.
ca_did_transforms <- function(outcome) {
  i <- match(outcome, ca_outcomes$outcome)
  if (is.na(i)) stop("Unknown outcome: ", outcome)
  c(if (isTRUE(ca_outcomes$log[i])) ca_log_transform(outcome),
    if (isTRUE(ca_outcomes$nonlog[i])) "none")
}

ca_did_apply_transform <- function(value, transform) {
  v <- as.numeric(value)
  if (identical(transform, "none")) return(v)
  if (identical(transform, "log1p")) {
    if (any(v < 0, na.rm = TRUE)) {
      stop("log1p requested on an outcome carrying negative values. Signed ",
           "outcomes are non-log only, see ca_outcomes$log.")
    }
    return(log1p(v))
  }
  stop("Unknown transform: ", transform)
}

# Filename-safe token. Stratum labels carry spaces, ecoregion codes are
# zero-padded strings.
ca_slug <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

# The run key. Every output file, every manifest row, and every stage 11 read
# resolves through this rather than through a filename pattern.
ca_did_run_id <- function(pairing, outcome, transform, stratum, ecoregion) {
  s <- ca_slug(stratum)
  e <- ca_slug(ecoregion)
  s[!nzchar(s)] <- "all"
  e <- ifelse(nzchar(e), paste0("eco", e), "statewide")
  # The field separator is a hyphen and the within-field separator is an
  # underscore, since ca_slug() turns every non-alphanumeric character
  # into an underscore and pairing labels carry underscores of their own.
  # So the five fields split unambiguously on the hyphen.
  paste(pairing, outcome, transform, s, e, sep = "-")
}

# THE STAGE 10 TASK TABLE.
#
# Panel manifest crossed with transform crossed with subset cell. One row is
# one array task, and nothing loops inside a task, so a wall clock hit costs
# one cell rather than a whole pairing. That is the stage 8 lesson, where a
# 72 hour timeout at 32 of 36 levels wrote nothing.
#
# Rows are ordered by tier first, so each tier is a contiguous array range and
# the three submissions carry their own wall clock. Take the ranges from the
# tasks mode rather than from memory, since they move when a level is added.
ca_did_runs <- function(estimator = c("etwfe", "attgt")) {
  estimator <- match.arg(estimator)

  mp <- ca_did_manifest_path()
  if (!file.exists(mp)) {
    stop("Missing ", mp, ". Stage 9b has not completed.")
  }
  lp <- ca_did_levels_path()
  if (!file.exists(lp)) {
    stop("Missing ", lp, ". Build it first:\n",
         "  Rscript 10b_did_etwfe.R levels")
  }

  man <- utils::read.csv(mp, stringsAsFactors = FALSE)
  lev <- utils::read.csv(lp, stringsAsFactors = FALSE,
                         colClasses = c(ecoregion = "character",
                                        stratum = "character"))
  lev$stratum[is.na(lev$stratum)]     <- ""
  lev$ecoregion[is.na(lev$ecoregion)] <- ""

  # The statewide all-stratum cell is the denominator for the row estimate. It
  # is a run for the pa and offset arms and a reference row for the others.
  base <- lev[lev$stratum == "" & lev$ecoregion == "",
              c("pairing", "n_units")]
  names(base)[2] <- "n_units_all"

  arm  <- ca_pairings$arm[match(lev$pairing, ca_pairings$pairing)]
  keep <- ifelse(arm %in% c("fire", "thin"), lev$stratum != "",
                 lev$stratum == "")
  lev <- lev[which(keep), , drop = FALSE]

  out <- do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
    tf <- ca_did_transforms(man$outcome[i])
    cl <- lev[lev$pairing == man$pairing[i], , drop = FALSE]
    if (!length(tf) || !nrow(cl)) return(NULL)
    g <- expand.grid(ti = seq_along(tf), ci = seq_len(nrow(cl)),
                     KEEP.OUT.ATTRS = FALSE)
    data.frame(
      pairing     = man$pairing[i],
      outcome     = man$outcome[i],
      transform   = tf[g$ti],
      stratum     = cl$stratum[g$ci],
      ecoregion   = cl$ecoregion[g$ci],
      n_units     = cl$n_units[g$ci],
      n_treated   = cl$n_treated[g$ci],
      n_pairs     = cl$n_pairs[g$ci],
      n_cohorts   = cl$n_cohorts[g$ci],
      panel_rows  = man$n_rows[i],
      year_first  = man$year_first[i],
      year_last   = man$year_last[i],
      path        = man$path[i],
      stringsAsFactors = FALSE)
  }))
  if (is.null(out) || !nrow(out)) {
    stop("No stage 10 runs. Check ", basename(mp), " and ", basename(lp), ".")
  }

  out$n_units_all <- base$n_units_all[match(out$pairing, base$pairing)]
  # as.numeric first. panel_rows on a statewide Almanac cell runs to tens of
  # millions and n_units to hundreds of thousands, and the product of two
  # integers overflows at 2.1 billion, which returns NA with a warning and
  # leaves every large cell untiered.
  out$est_rows <- round(as.numeric(out$panel_rows) *
                        as.numeric(out$n_units) /
                        as.numeric(out$n_units_all))
  out$est_cost <- out$est_rows * as.numeric(out$n_cohorts)
  out$tier <- ifelse(
    out$est_cost <= ca_did_est$tier_breaks[["small"]], "small",
    ifelse(out$est_cost <= ca_did_est$tier_breaks[["medium"]],
           "medium", "large"))
  out$emfx_type <- vapply(out$pairing, ca_did_emfx_type, character(1))
  out$base_shift <- vapply(out$pairing, ca_did_base_shift, integer(1))
  out$run_id <- ca_did_run_id(out$pairing, out$outcome, out$transform,
                              out$stratum, out$ecoregion)
  out$estimator <- estimator

  # The floor on treated units, applied to every cell. Below it a cell is
  # not a small estimate, it is no estimate, and it costs an array task and
  # an error row in the manifest. Set from the observed distribution in
  # meta/did_subset_levels.csv rather than chosen, following
  # screen_threshold.csv and the config retention ratio. Zero disables it.
  if (ca_did_est$min_treated_units > 0L) {
    out <- out[out$n_treated >= ca_did_est$min_treated_units, ,
               drop = FALSE]
  }

  out <- out[order(match(out$tier, c("small", "medium", "large")),
                   match(out$pairing, ca_pairings$pairing),
                   match(out$outcome, ca_outcomes$outcome),
                   out$transform, out$stratum, out$ecoregion), , drop = FALSE]
  rownames(out) <- NULL

  # WHICH ARRAY TASK COMPUTES A CELL.
  #
  # etwfe submits one task per cell. att_gt does not, because the fit is
  # minutes rather than hours and one task per ecoregion would spend more
  # scheduler time than estimator time. An att_gt task is a statewide cell and
  # it computes its own ecoregions in a loop, reusing the one panel read it
  # has already paid for. 10b re-reads the panel for every ecoregion because
  # every ecoregion is its own task there.
  #
  # The task index is the row number in ca_did_array(), which is the statewide
  # subset of this table in this same order, so the two cannot drift.
  if (identical(estimator, "attgt")) {
    key <- paste(out$pairing, out$outcome, out$transform, out$stratum)
    out$task <- match(key, key[out$ecoregion == ""])
    # An ecoregion cell whose statewide parent fell below the treated-unit
    # floor has no task to ride in. It cannot happen while the floor is a
    # count, since the statewide cell contains the ecoregion, but it is
    # dropped rather than left pointing at nothing.
    out <- out[!is.na(out$task), , drop = FALSE]
    rownames(out) <- NULL
  } else {
    out$task <- seq_len(nrow(out))
  }

  out
}

# ---------------------------------------------------------------------------
# 0.12  SAMPLING
# ---------------------------------------------------------------------------
# Stage 7 does one thing, spacing. It reduces every group that appears on the
# treated side of a pairing to at most one pixel per 150 m block, subject to a
# 150 m minimum separation. The 100,000 cap is a stage 8 operation and the
# order between them is not free, since capping first would leave roughly
# 4,000 units after spacing rather than 100,000.
#
# decimate_control is retired. It was a switch on whether control pools were
# spaced as well, which made sense while a group held one role only. Under the
# six-pairing structure UP is treated in UP_UNP and a control pool in FP_UP and
# TP_UP, so the question is not answered by a flag on the file, it is answered
# by which file stage 8 reads for which role. Spaced treated groups come from
# ca_sample_path(). Whole control pools come from ca_group_path().
#
# ---------------------------------------------------------------------------
ca_sample <- list(

  grid_spacing_m   = 150,
  min_distance_m   = 150,
  block_pixels     = 5L,
  block_offset     = 2L,

  # WHO IS SPACED
  #
  # Treated-role groups only. Spacing exists to stop adjacent 30 m pixels being
  # counted as independent observations of the same stand, and the submitted
  # run applied it to the treated side only. 8_2 filters the treated set by
  # all.smpl.V1 and leaves the control set at full density.
  #
  # The group list is derived from ca_pairings$treat and is never written out
  # here, so UNP and OP fall out by construction. UP is spaced in its treated
  # role for UP_UNP and read whole from ca_group_path() in its control role for
  # FP_UP and TP_UP. The two versions coexist and nothing reconciles them.
  #
  # The asymmetry between the arms remains open and is settled by measurement,
  # not assertion. The stage 8 verifier computes the nearest-neighbour distance
  # distribution among retained control pixels per pairing. pixel_id is the
  # terra cell index, so row and column follow from integer arithmetic and no
  # spatial package is needed.
  spaced_role      = "treatment",

  # HOW THEY ARE SPACED
  #
  # One representative per occupied 150 m block, best centred first and
  # ascending pixel_id on ties, then a second pass enforcing the 150 m minimum
  # between representatives on either side of a shared block edge.
  #
  # This replaces the fixed-offset lattice of the first draft, which kept a
  # pixel only where it was the centre of its own block. That rule is exact for
  # a contiguous surface like UP and severe for a fragmented group like T1NP,
  # where a block holding group pixels but not the centre pixel contributes
  # nothing. Retention now scales with fragmentation, which is correct. A group
  # whose pixels already lie more than 150 m apart is not spatially redundant
  # and is left intact.
  #
  # Deterministic in both passes. No seed, no distance search, no spatial
  # package, bit-identical on rerun. This is what retires
  # spatialEco::subsample.distance in 7_grid_based_sampling.R and the R to
  # ArcMap to R detour in 7_2 and 7_3.
  method           = "block_representative",

  # THE EXCEPTION LIST
  #
  # min_group_n is retired and the per-cell floor that briefly replaced it is
  # retired too. The probe of 3 August measured both and neither survives.
  #
  # A threshold on group size has no bearing on the reason spacing exists.
  # A threshold on units per occupied exact cell measures the wrong side of the
  # match. Matching is treated to control, so what gates feasibility is control
  # abundance in the cell, and control pools are not spaced. Treated density per
  # cell tracks class size instead. Applying a floor of 30 would have taken
  # every Deciduous fire group whole, at a median of 10 to 18 units per cell,
  # while the same groups in Evergreen pass at 97 to 298. The difference is that
  # Deciduous holds 2 percent of the grid, not that Deciduous fire pixels are
  # less spatially redundant.
  #
  # The estimation-relevant quantity is the pooled spaced count across classes,
  # since lulc is an exact matching variable and matching runs within class
  # while estimation pools. On that basis only T1NP is small, at 219 units
  # pooled from 1,174, and T1NP is underpowered whether it is spaced or not.
  #
  # So spacing is uniform and any exemption is named here rather than derived
  # from a threshold. An empty vector means every treated group is spaced,
  # which is one Methods sentence with no exception to defend. Groups too small
  # to estimate are handled at stage 8 by a reporting gate on matched pairs,
  # which is where sample size belongs.
  whole_groups     = character(0),

  # THE CAP
  #
  # Stage 8, treated side, after spacing. Capping first and spacing second
  # would leave roughly 4,000 units rather than 100,000, so the order is not
  # free. The draw is random and carries ca_seed, and it is the only stochastic
  # step in the pipeline before matching.
  #
  # SCOPE. Per analysis group, not per pairing. The cap is a limit on each
  # severity and intensity level, so the treated side of FP_UP reaches up to
  # 300,000 across F2P, F3P, and F4P rather than being held to 100,000 in
  # total. An earlier comment in this file read it as per pairing and was
  # wrong.
  #
  # POOLED ACROSS CLASSES. One draw of 100,000 per analysis group covering
  # Decid, Everg, and Mixed together, matching the submitted run, where 8_1 L52
  # subsets one combined input against All_sampling_sum.csv. Spacing ran per
  # class because it is a geometric operation, the cap does not.
  #
  # ORDER AGAINST COMPLETE CASES. Cap first, complete-case filtering second and
  # per arm. The draw is then one seeded draw shared by all three arms, and
  # each arm loses only what its own covariates cost it. Filtering on the
  # augmented covariates first would give the three arms an identical treated
  # set at the price of restricting baseline by covariates it does not use,
  # which would break the claim that baseline reproduces the submitted
  # specification. The per-arm loss is written to match_specs so the difference
  # is visible rather than hidden.
  cap_per_group        = 100000L,
  cap_stage            = 8L,
  cap_applies_to       = "treatment",
  cap_scope            = "analysis_group",
  cap_pooled_across_lulc = TRUE,
  cap_before_complete_cases = TRUE
)

ca_seed <- 160617L

# ---------------------------------------------------------------------------
# 0.13  IO CONVENTIONS
# ---------------------------------------------------------------------------
# Tabular intermediates are Parquet. Spatial layers are GeoPackage. Shapefile
# is retired for its 2 GB ceiling and 10-character field-name truncation.
# CURC advises that all filesystems perform poorly on many small files, so each
# stage consolidates its part files on completion.

ca_io <- list(
  tabular_format = "parquet",
  spatial_format = "gpkg",
  compression    = "zstd",
  consolidate    = TRUE
)

# ---------------------------------------------------------------------------
# 0.13b  DURABILITY
# ---------------------------------------------------------------------------
# Destination is a property of the artifact, not a decision made at write time.
# Scratch purges on an access-based schedule and the GCB deadline is 18
# September 2026, so anything produced now and read again in September has to
# be on /projects rather than merely still present on scratch.
#
#   scratch   regenerable by rerunning one array from an input that is itself
#             archived. Chunk intermediates, staged copies, temporary reads
#   archive   input to matching or DiD, expensive to regenerate, or quoted in
#             the manuscript. Everything downstream of stage 4
#
# The split inside stage 3a is by whether the upstream is external. The Almanac
# is streamed from source.coop, so re-extraction depends on an outside service
# and those extracts are archived. NCSDA, eMapR, and LEMMA masters already sit
# on /projects, so re-extraction from them is local and cheap and those
# extracts stay regenerable.

ca_durability <- data.frame(
  artifact = c("grid", "extract_almanac", "extract_screens", "extract_mtbs",
               "extract_local", "extract_vector", "extract_static",
               "covariates", "verify", "treat", "screen", "groups",
               "sample", "matched", "did_panel", "model", "figure"),
  tier     = c("archive", "archive", "archive", "archive",
               "scratch", "archive", "archive",
               "archive", "archive", "archive", "archive", "archive",
               "archive", "archive", "archive", "archive", "archive"),
  note     = c("frozen, stage 1",
               "streamed from source.coop, external dependency",
               "small, read repeatedly by stages 5 and 6",
               "small, sparse, read by stages 5 and 6",
               "NCSDA, eMapR, LEMMA. Masters are local, re-extraction is cheap",
               "stage 3c long records",
               "stage 3d one row per pixel",
               "stage 4b matching covariate table",
               "every verification CSV",
               "stage 5a", "stage 5b", "stage 6",
               "stage 7, deterministic, must be byte-identical on rerun",
               "stage 8, cohort-stratified, seed 160617 on the cap draw only",
               "stage 9", "stage 10", "stage 11"),
  stringsAsFactors = FALSE
)

ca_tier <- function(artifact) {
  i <- match(artifact, ca_durability$artifact)
  if (is.na(i)) stop("Unknown artifact class: ", artifact)
  ca_durability$tier[i]
}

# Archive if the artifact class says so, otherwise leave it on scratch. Called
# at the end of a stage, never inside an array task.
ca_persist <- function(path, artifact, subdir = "work") {
  if (ca_tier(artifact) != "archive") {
    ca_log("Tier scratch, not archived: ", basename(path))
    return(invisible(NULL))
  }
  ca_archive(path, subdir = subdir)
}

# ---------------------------------------------------------------------------
# 0.13c  STAGE PATHS
# ---------------------------------------------------------------------------
# Stages 5a, 5b, and 4b each defined their output path locally. Stages 6, 7,
# and 8 all read them, so they move here. The local definitions in 5a and 5b
# can be deleted, they resolve identically.

ca_event_path <- function(arm, lulc) {
  dir <- ca_work("recode", "events")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("events_%s_%s.parquet", arm, lulc))
}

ca_screen_path <- function(lulc) {
  dir <- ca_work("recode", "screen")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("screen_%s.parquet", lulc))
}

ca_covariate_path <- function(lulc) {
  dir <- ca_work("clean", "covariates")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("covariates_%s.parquet", lulc))
}

ca_group_path <- function(lulc) {
  dir <- ca_work("groups")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("groups_%s.parquet", lulc))
}

ca_sample_path <- function(lulc) {
  dir <- ca_work("sample")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, sprintf("sample_%s.parquet", lulc))
}

# Which file stage 8 reads, by role. Structural rather than remembered, because
# the one error this design can produce is reading the spaced file for a control
# pool, and that error is silent. It would not fail, it would quietly shrink UP
# from 25 million to 1.2 million in FP_UP and TP_UP and inflate every standard
# error in the fire and thinning arms on protected land.
#
# Treated groups are spaced, stage 7. Control pools are whole, stage 6. UP holds
# both roles and resolves to a different file in each, which is the whole point
# of separating the two stages.
ca_pool_path <- function(lulc, role) {
  role <- match.arg(role, c("treatment", "control"))
  if (role == "treatment") ca_sample_path(lulc) else ca_group_path(lulc)
}

# Stage 8 outputs. One matched file per pairing, arm, and config. Classes are
# pooled inside a file because lulc is an exact matching variable, which is
# what the submitted run did when it fed one combined input to matchit().
ca_match_dir <- function(...) {
  dir <- ca_work("match", ...)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

ca_match_path <- function(pairing, arm, config) {
  file.path(ca_match_dir(pairing),
            sprintf("matched_%s_%s_%s.parquet", pairing, arm, config))
}

# Two stream-appended tables, written by every task under a file lock. These
# are the sole source for Table S3 and Figure S6, which is why they are single
# files rather than one per task.
ca_match_specs_path   <- function() file.path(ca_match_dir(),
                                              "match_specs.csv.gz")
ca_match_balance_path <- function() file.path(ca_match_dir(),
                                              "match_balance.csv.gz")

# Stage 9 outputs. One units file per pairing at the top level, one panel per
# pairing and outcome in a per-pairing subdirectory, matching the stage 8
# layout so a pairing is one directory at every stage downstream of matching.
ca_did_dir <- function(...) {
  dir <- ca_work("did", ...)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

ca_did_units_path <- function(pairing) {
  file.path(ca_did_dir(), sprintf("did_units_%s.parquet", pairing))
}

ca_did_panel_path <- function(pairing, outcome) {
  file.path(ca_did_dir(pairing),
            sprintf("did_panel_%s_%s.parquet", pairing, outcome))
}

# The outcome-to-file map stage 10 reads. Written by 9b, never by hand.
ca_did_manifest_path <- function() ca_meta(ca_did$manifest)

# Stage 10 outputs. One CSV per run under a per-pairing subdirectory, matching
# the stage 8 and stage 9 layout, plus two meta tables. The levels table is a
# data artifact in the same sense as screen_threshold.csv, built once from the
# stage 9 units files and read by both estimator scripts, so the ecoregions and
# strata a pairing actually holds are recorded rather than hardcoded.
ca_did_est_dir <- function(...) {
  dir <- ca_work("did_est", ...)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

ca_did_est_path <- function(estimator, pairing, run_id) {
  file.path(ca_did_est_dir(pairing),
            sprintf("%s_%s.csv.gz", estimator, run_id))
}

ca_did_levels_path <- function() ca_meta(ca_did_est$levels_file)

ca_did_est_manifest_path <- function(estimator) {
  estimator <- match.arg(estimator, names(ca_did_est$manifest))
  ca_meta(unname(ca_did_est$manifest[[estimator]]))
}

# Stage 10c outputs. Three small tables in meta/, which is on /projects and is
# backed up, so they survive the scratch purge without an archive step and are
# what gets copied to the laptop for stage 11.
#
# results   every estimate from both estimators, one row per term, with the
#           log1p rows back-transformed and se_ratio_prev attached
# coverage  completeness by estimator and treatment arm at the moment results
#           was written, so a figure built on a draining array is identifiable
# missing   the runs that are todo, error, or still queued, with the array
#           index so a resubmission targets exactly those tasks
ca_did_results <- list(
  results  = "did_results.csv.gz",
  coverage = "did_results_coverage.csv",
  missing  = "did_results_missing.csv"
)

ca_did_results_path <- function(what = c("results", "coverage", "missing")) {
  what <- match.arg(what)
  ca_meta(ca_did_results[[what]])
}

# pixel_id is the terra cell index of the NLCD 2001 grid, so grid position is
# integer arithmetic on it. The column count comes from the dim element of the
# frozen grid reference written at stage 1, 1_forest_point_grid.R L235.
ca_grid_ncol <- function() {
  if (!file.exists(ca_grid$meta_rds)) {
    stop("Grid reference missing: ", ca_grid$meta_rds)
  }
  as.integer(readRDS(ca_grid$meta_rds)$dim[["ncol"]])
}

# ---------------------------------------------------------------------------
# 0.14  HELPERS
# ---------------------------------------------------------------------------

ca_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

ca_task_id <- function() {
  id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nzchar(id)) return(as.integer(id))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) as.integer(args[1]) else NA_integer_
}

ca_threads <- function() {
  n <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (!nzchar(n)) n <- Sys.getenv("SLURM_NTASKS", unset = "1")
  max(1L, as.integer(n))
}

ca_log <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

ca_stamp <- function(script) {
  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  out <- ca_meta(paste0("session_", script, "_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  writeLines(capture.output(sessionInfo()), out)
  invisible(out)
}

# Copy a completed intermediate from scratch to the backed-up archive. Scratch
# purges 90 days after file creation and the policy must not be circumvented,
# so anything expensive to rebuild is archived explicitly.
#
# file.copy() copies a directory INTO an existing destination directory. Passing
# a not-yet-existing full destination path silently drops the copy with only a
# warning, so the destination directory is created first and the return value is
# checked.
ca_archive <- function(path, subdir = "work") {
  if (!file.exists(path)) stop("Nothing to archive at ", path)

  dst_dir <- file.path(ca_paths$archive, subdir)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)

  ok <- if (dir.exists(path)) {
    file.copy(path, dst_dir, overwrite = TRUE, recursive = TRUE)
  } else {
    file.copy(path, file.path(dst_dir, basename(path)), overwrite = TRUE)
  }

  dst <- file.path(dst_dir, basename(path))
  if (!all(ok) || !file.exists(dst)) {
    stop("Archive copy failed: ", path, " -> ", dst)
  }

  n_src <- length(list.files(path, recursive = TRUE))
  n_dst <- length(list.files(dst, recursive = TRUE))
  if (dir.exists(path) && n_src != n_dst) {
    stop("Archive incomplete. Source files ", n_src, ", destination ", n_dst)
  }

  ca_log("Archived ", basename(path), " to ", dst_dir,
         " (", n_dst, " files, ",
         round(sum(file.size(list.files(dst, recursive = TRUE,
                                        full.names = TRUE))) / 1e9, 2), " GB)")
  invisible(dst)
}

ca_pending <- function() {
  data.frame(
    item = c(
      "Aspect replaced by northness and eastness, response letter note",
      "Screen threshold derived at 0, replaces the 5 percent Methods text",
      "Joint CPAD units, 81 features, in protected_other rather than treatment",
      "CAFR5170 Early Action predecessor beginning 2004 against a 2015 start",
      "Private thinning record starts 1991, public starts in the 1960s",
      "Carbon_AGB capped at 800 t/ha by the producer, upper tail truncated",
      "Methods text correction, 120 m spacing should read 150 m",
      "Methods text correction, thinning intensity is a crosswalk not basal area",
      "Detection requirement retention per group, Methods number",
      "OP group count, audit of the private-land-only offset claim",
      "Control pool decimation, decided from the stage 8 probe distance run",
      "Cohort-stratified matching, Methods sentence and the control reuse rate",
      "Distance model fitted on the matchable pool, Methods sentence",
      "ONP lost 25,225 tribal-land offset pixels, response letter number",
      "pa_after_cohort count, Methods sentence on the time-varying PA label"
    ),
    blocks = c("manuscript", "manuscript", "stage 6", "stage 8",
               "manuscript", "manuscript", "manuscript", "manuscript",
               "stage 6", "stage 6", "stage 8", "manuscript", "manuscript",
               "manuscript", "manuscript"),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 0.15  SESSION OPTIONS
# ---------------------------------------------------------------------------

options(stringsAsFactors = FALSE, scipen = 999)

Sys.setenv(TMPDIR = ca_tmpdir())

if (requireNamespace("terra", quietly = TRUE)) {
  terra::terraOptions(tempdir = ca_tmpdir(), memfrac = 0.6, progress = 0)
}

if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::set_cpu_count(ca_threads())
}

invisible(TRUE)

# ---------------------------------------------------------------------------
# 0.16  FIGURES
# ---------------------------------------------------------------------------
# Everything the stage 11 scripts would otherwise hardcode. The submitted
# figure scripts carried the ecoregion set as c(78, 1, 9, 4, 5, 8) at
# 11_4_2 L27, the outcome set as c("GPP","eMapR","NCSDA","LEMMA") at
# 11_4_1 L27, the caps as mng.yr at 11_4_1 L25, and the y limits as
# c(-113, 100) at 11_4_1 L144. All four are here instead.
#
# WHERE THE FIGURES RUN
#
# Stage 11 runs on the laptop against did_results.csv.gz. The exceptions are
# the Figure 1 and S1 inputs, which touch the 104 million pixel grid and are
# reduced on Alpine by 11a. Set CA_FIG_DIR on the laptop to point the output
# somewhere local. Everything else resolves through the usual paths and only
# ever reads.

ca_fig <- list(
  dir    = Sys.getenv("CA_FIG_DIR", unset = ""),
  
  ## THE WINDOW FILE STAGE 11 READS
  ##
  ## fig_caps.csv is what 10d section 7 writes from the three support floors.
  ## fig_caps_rev.csv is the hand-checked copy, cut at horizons the automatic
  ## rule allowed and the panels could not carry. The revised file is preferred
  ## where it exists, so both records survive and the difference between the
  ## two files is the record of what was cut. CA_FIG_CAPS forces one or the
  ## other without editing code.
  caps_file = Sys.getenv("CA_FIG_CAPS", unset = ""),
  
  device = "svg",
  dpi    = 400,
  width  = list(single = 3.35, double = 6.9, full = 9.5),
  
  ecoregions = c(
    "01" = "Coast Range",
    "04" = "Cascades",
    "05" = "Sierra Nevada",
    "08" = "S. California Mountains",
    "09" = "Eastern Cascades",
    "78" = "Klamath Mountains"
  ),
  
  ecoregion_colors = c(
    "Statewide"                     = "black",
    "Coast Range"                   = "#e41a1c",
    "Cascades"                      = "#377eb8",
    "Sierra Nevada"                 = "#4daf4a",
    "S. California Mountains"       = "#984ea3",
    "Eastern Cascades"              = "#ff7f00",
    "Klamath Mountains"             = "#e6ab02"
  ),
  
  panel_rows = list(
    flux  = list(outcomes  = c("gpp", "npp_tree", "nep"),
                 transform = "none",
                 axis      = "gC m-2 yr-1"),
    stock = list(outcomes  = c("agb_almanac", "agb_emapr", "agb_lemma"),
                 transform = "log1p",
                 axis      = "Change (%)")
  ),
  
  outcome_short = c(
    gpp         = "GPP",
    npp_tree    = "NPP",
    nep         = "NEP",
    nbp         = "NBP",
    agb_almanac = "AGB (Almanac)",
    agb_emapr   = "AGB (eMapR)",
    agb_lemma   = "AGB (LEMMA)"
  ),
  
  caps = list(
    pa = list(
      type = "calendar",
      lo   = c(gpp = 1990, npp_tree = 1990, nep = 1990,
               agb_almanac = 1990, agb_emapr = 1990, agb_lemma = 1990),
      hi   = c(gpp = 2025, npp_tree = 2021, nep = 2021,
               agb_almanac = 2025, agb_emapr = 2017, agb_lemma = 2016)),
    # att_gt has no calendar aggregation. Its dynamic aggregation is event
    # time, so the att_gt companion to Figure 2 is drawn on event time and
    # needs an event cap the calendar entry above cannot supply. The hi values
    # are built the same way the fire and thin entries are, as the outcome's
    # last data year minus the first treated cohort year, which is 1990 under
    # the calendar-year PA specification because everything earlier is
    # always-treated and carries no event time. 2025 - 1990, 2021 - 1990,
    # 2017 - 1990, 2016 - 1990.
    pa_event = list(
      type = "event",
      lo   = c(gpp = 0, npp_tree = 0, nep = 0,
               agb_almanac = 0, agb_emapr = 0, agb_lemma = 0),
      hi   = c(gpp = 35, npp_tree = 31, nep = 31,
               agb_almanac = 35, agb_emapr = 27, agb_lemma = 26)),
    offset = list(
      type = "event",
      lo   = c(gpp = 0, npp_tree = 0, nep = 0,
               agb_almanac = 0, agb_emapr = 0, agb_lemma = 0),
      hi   = c(gpp = 10, npp_tree = 6, nep = 6,
               agb_almanac = 10, agb_emapr = 3, agb_lemma = 2)),
    fire = list(
      type = "event",
      lo   = c(gpp = 0, npp_tree = 0, nep = 0,
               agb_almanac = 0, agb_emapr = 0, agb_lemma = 0),
      hi   = c(gpp = 34, npp_tree = 30, nep = 30,
               agb_almanac = 34, agb_emapr = 26, agb_lemma = 25)),
    thin = list(
      type = "event",
      lo   = c(gpp = 0, npp_tree = 0, nep = 0,
               agb_almanac = 0, agb_emapr = 0, agb_lemma = 0),
      hi   = c(gpp = 34, npp_tree = 30, nep = 30,
               agb_almanac = 34, agb_emapr = 26, agb_lemma = 25))
  ),
  
  pretrend = c(lo = -9, hi = -2),
  
  ownership_colors = c("Protected" = "#2C6E49", "Non-protected" = "#A15C38"),
  severity_colors  = c("Low severity"      = "#F6C667",
                       "Moderate severity" = "#E07A5F",
                       "High severity"     = "#8C2F39"),
  intensity_colors = c("Low intensity"    = "#A8DADC",
                       "Medium intensity" = "#457B9D",
                       "High intensity"   = "#1D3557"),
  stage_colors     = c("Before matching" = "#4C72B0",
                       "After matching"  = "#C44E52"),
  
  covariate_labels = c(
    elevation        = "Elevation",
    slope            = "Slope",
    northness        = "Northness",
    eastness         = "Eastness",
    ppt_normal       = "Precipitation",
    tmean_normal     = "Air temperature",
    pop_density      = "Population density",
    city_travel_time = "Travel time to cities",
    gpp_pre5         = "Pre-treatment GPP",
    agb_pre5         = "Pre-treatment AGB",
    ch_pre5          = "Pre-treatment canopy height"
  ),
  
  smd_threshold = 0.1,
  
  pairing_labels = c(
    UP_UNP  = "Protected areas",
    ONP_UNP = "Offset programs",
    FP_UP   = "Wildfire, protected",
    FNP_UNP = "Wildfire, non-protected",
    TP_UP   = "Thinning, protected",
    TNP_UNP = "Thinning, non-protected"
  ),
  
  coarsen_factor = 33L,
  disturbance_labels = c(pa      = "Undisturbed",
                         control = "Undisturbed",
                         offset  = "Offset",
                         fire    = "Wildfire",
                         thin    = "Thinning"),
  disturbance_order  = c("Undisturbed", "Wildfire", "Thinning", "Offset"),
  disturbance_colors = c("Undisturbed" = "#B7B7A4",
                         "Wildfire"    = "#BC4749",
                         "Thinning"    = "#6A994E",
                         "Offset"      = "#386641")
)

# Everything below is measured off a submitted file. The two reference maps
# were decoded and sampled rather than eyeballed, so the numbers here are what
# the published figures actually contain.
#
#   11_2_dist_chng_bar.R          L65 labels, L69 strips, L96 and L131 fills,
#                                 L98 and L106 type sizes
#   ca_eco_l3_prj.svg             map type 15.84 pt bold, line spacing 18.5 pt,
#                                 label positions, land fill #F7F7F7, ecoregion
#                                 polygons filled with the strip colours and no
#                                 border, state outline stroke 0.48 pt
#   CA_undist_dist_all.svg        the fifteen group colours, legend swatch
#                                 28.1 by 13.9 pt stepping 20.64 pt, text
#                                 13.92 pt bold
#   CA_undist_dist_all_eco.svg    the framing, the ecoregion boundary colours
#                                 and their 9 pixel width, and the labels
#
# TWO LABEL SETS, NOT ONE
#
# The bar strips carry "Eastern Cascades\nSlopes and Foothills" at 11_2 L65.
# The maps carry "Eastern\nCascades" and nothing more, because four words do
# not fit inside the polygon. Both are correct in their own figure.

ca_fig1 <- list(
  eco_strip_colors = c("#e41a1c", "#377eb8", "#4daf4a",
                       "#984ea3", "#ff7f00", "#e6ab02"),
  eco_facet_labels = c(
    "01" = "Coast Range",
    "04" = "Cascades",
    "05" = "Sierra Nevada",
    "08" = "Southern California\nMountains",
    "09" = "Eastern Cascades\nSlopes and Foothills",
    "78" = "Klamath Mountains"),
  eco_map_labels = c(
    "01" = "Coast\nRange",
    "04" = "Cascades",
    "05" = "Sierra\nNevada",
    "08" = "Southern California\nMountains",
    "09" = "Eastern\nCascades",
    "78" = "Klamath\nMountains"),
  map_frame = c(left = 0.1435, right = 0.0706, top = 0.0160, bottom = 0.0781),
  eco_label_pos = data.frame(
    code = c("01", "04", "05", "08", "09", "78"),
    fx   = c(0.1510, 0.3753, 0.4457, 0.4318, 0.4631, 0.2358),
    fy   = c(0.7010, 0.8208, 0.5883, 0.2670, 0.8896, 0.8796),
    stringsAsFactors = FALSE),
  map_text_pt    = 15.84,
  map_text_face  = "bold",
  map_lineheight = 1.168,
  map_dev_width  = 8.27,
  map_dev_height = 11.69,
  lulc_stack  = c("Decid", "Mixed", "Everg"),
  lulc_colors = c(Decid = "#68ab5f", Mixed = "#b5c58f", Everg = "#1c5f2c"),
  lulc_legend = c(Decid = "Deciduous", Mixed = "Mixed", Everg = "Coniferous"),
  dist_stack  = c("F2", "F3", "F4", "T1", "T2", "T3", "O", "U"),
  dist_colors = c(F2 = "#fee0d2", F3 = "#fc9272", F4 = "#de2d26",
                  T1 = "#e5f5e0", T2 = "#a1d99b", T3 = "#31a354",
                  O  = "#edf8b1", U  = "#a6cee3"),
  dist_legend = c(F2 = "Fire, low",  F3 = "Fire, mid",  F4 = "Fire, high",
                  T1 = "Thin, low",  T2 = "Thin, mid",  T3 = "Thin, high",
                  O  = "Offset",     U  = "Undisturbed"),
  y_limit     = c(0, 20),
  y_breaks    = seq(0, 20, by = 5),
  base_text   = 24,
  strip_text  = 16,
  dev_width   = 14,
  dev_height  = 11,
  rel_heights = c(1.1, 1.3),
  group_colors = c(
    UNP  = "#FFFF99", UP   = "#8DD3C7",
    ONP  = "#1F78B4",
    F2NP = "#FDBE85", F2P  = "#D7B5D8",
    F3NP = "#FD8D3C", F3P  = "#DF65B0",
    F4NP = "#D94701", F4P  = "#CE1256",
    T1NP = "#BAE4B3", T1P  = "#CBC9E2",
    T2NP = "#74C476", T2P  = "#9E9AC8",
    T3NP = "#238B45", T3P  = "#6A51A3"),
  group_order = c("UNP", "UP", "ONP",
                  "F2NP", "F2P", "F3NP", "F3P", "F4NP", "F4P",
                  "T1NP", "T1P", "T2NP", "T2P", "T3NP", "T3P"),
  s1_legend_pos    = c(0.0778, 0.469),
  s1_legend_just   = c(0, 1),
  s1_legend_width  = 27.4,
  s1_legend_height = 20.64,
  s1_legend_text   = 13.92,
  land_fill = "#F7F7F7",
  boundary_color = "#000000",
  boundary_width = 0.225,
  eco_line_width = 1.0,
  simplify_m = 300,
  eco_shapefile = Sys.getenv(
    "CA_ECO_SHP",
    unset = "D:/CA_data/Ecoregion/ca_eco_l3_prj.shp"),
  boundary_shapefile = Sys.getenv(
    "CA_BOUNDARY_SHP",
    unset = "D:/CA_data/Boundary/CA_State_TIGER2016_proj.shp")
)

# The DiD trajectory figures. Everything here is read out of the submitted
# scripts rather than chosen again.
#
#   11_4_0_DiD_eco_log.R    L46 ecoregion colours, L47 shapes, L138 to L167
#                           the panel theme, L163 the in-panel title, L385
#                           the cowplot assembly and its 28 by 7 device
#   11_4_1_DiD_PA_NP_log.R  L92 ownership colours, fills and shapes, L98 to
#                           L112 the theme, L133 to L149 the geometry sizes,
#                           L156 the 26 by 7 device
#
# WHAT CHANGED, AND WHY
#
# Four outcomes became six, so the panel grid is two rows of three rather than
# one row of four. The two rows carry different units, which no single y axis
# title can express, so each row is built as its own plot and the two are
# stacked with cowplot. That is the same assembly 11_4_0 L385 already used to
# put one shared legend under a grid of panels.
#
# The in-panel outcome label moves from the fixed x = 1 and y = 95 at
# 11_4_1 L131 to a fixed fraction of each panel's own axis window, placed with
# hjust = 0. A fixed coordinate lands off the panel under free scales, and
# -Inf with a non-zero hjust indents each label by a fraction of its own string
# width, which is a different indent per panel. See ca_fig_panel_label().
#
# panel_label_grid is separate from panel_label because the ecoregion grids
# draw a two or three line label inside a panel a fraction of a main-text
# panel's width. At PANEL_W = 7 the two sizes coincide. They stop coinciding
# the moment the grid geometry changes, which is why the entry exists now
# rather than when it is first needed.
#
# Fire and thinning now carry three strata and two ownerships in one figure
# rather than one figure per stratum, so six lines share a panel. The ribbon
# alpha of 0.5 at 11_4_1 L136 hides five of them, hence a separate alpha for
# the stratified figures. Everything else is unchanged.

ca_fig2 <- list(
  eco_shapes = c("Statewide"                     = 18,
                 "Coast Range"                   = 16,
                 "Cascades"                      = 17,
                 "Sierra Nevada"                 = 15,
                 "Southern California Mountains" = 3,
                 "Eastern Cascades"              = 4,
                 "Klamath Mountains"             = 8),
  
  own_labels = c("Protected" = "PAs", "Non-protected" = "non-PAs"),
  own_colors = c("Protected" = "#ca0020", "Non-protected" = "#0571b0"),
  own_fills  = c("Protected" = "#f4a582", "Non-protected" = "#92c5de"),
  own_shapes = c("Protected" = 19, "Non-protected" = 17),
  own_lines  = c("Protected" = "solid", "Non-protected" = "22"),
  
  strat_order = c(
    "High_Protected", "Moderate_Protected", "Low_Protected",
    "High_Non-protected", "Moderate_Non-protected", "Low_Non-protected"
  ),
  
  strat_labels = c(
    "High_Protected"         = "High, PAs",
    "Moderate_Protected"     = "Moderate, PAs",
    "Low_Protected"          = "Low, PAs",
    "High_Non-protected"     = "High, non-PAs",
    "Moderate_Non-protected" = "Moderate, non-PAs",
    "Low_Non-protected"      = "Low, non-PAs"
  ),
  
  strat_colors = c(
    "High_Protected"         = "#67000d",
    "Moderate_Protected"     = "#cb181d",
    "Low_Protected"          = "#fb6a4a",
    "High_Non-protected"     = "#08306b",
    "Moderate_Non-protected" = "#2171b5",
    "Low_Non-protected"      = "#6baed6"
  ),
  
  strat_fills = c(
    "High_Protected"         = "#fcbba1",
    "Moderate_Protected"     = "#fcbba1",
    "Low_Protected"          = "#fcbba1",
    "High_Non-protected"     = "#c6dbef",
    "Moderate_Non-protected" = "#c6dbef",
    "Low_Non-protected"      = "#c6dbef"
  ),
  
  strat_shapes = c(
    "High_Protected"         = 1,
    "Moderate_Protected"     = 2,
    "Low_Protected"          = 5,
    "High_Non-protected"     = 1,
    "Moderate_Non-protected" = 2,
    "Low_Non-protected"      = 5
  ),
  
  strat_lines = c(
    "High_Protected"         = "solid",
    "Moderate_Protected"     = "dashed",
    "Low_Protected"          = "dotted",
    "High_Non-protected"     = "solid",
    "Moderate_Non-protected" = "dashed",
    "Low_Non-protected"      = "dotted"
  ),
  
  axis_text    = 30,
  axis_title   = 30,
  legend_text  = 30,
  legend_grid  = 32,
  panel_label  = 10,
  panel_label_grid = 10,
  
  hline_width   = 1.5,
  line_width    = 1.2,
  point_size    = 3.5,
  point_stroke  = 1.2,
  border_width  = 1.2,
  ribbon_alpha  = 0.5,
  ribbon_alpha_strat = 0.18,
  
  dev_width    = 21,
  dev_height   = 15,
  rel_heights  = c(1, 1, 0.25),
  plot_margin  = c(0.2, 0.6, 0.2, 0.2),
  
  calendar_breaks = seq(1990, 2020, by = 10),
  y_axis_angle    = 90,
  
  flux_ylab    = expression("gC m"^-2 ~ "yr"^-1),
  flux_ylab_kg = expression("kgC m"^-2 ~ "yr"^-1),
  stock_ylab   = "%",
  
  x_titles = c(
    pa     = "Year",
    offset = "Years since project commencement",
    fire   = "Years since wildfire",
    thin   = "Years since thinning detection")
)

# Points to the millimetres ggplot sizes text in.
ca_fig_pt <- function(pt) pt / 2.845276

# THE PLOT FRAME, reconstructed from the reference margins.
# Takes the bounding box of the California polygon and returns the extent the
# reference map drew, so the legend has the same space beside the state and the
# label fractions land where they land in the published figure.
ca_fig_frame <- function(bnd_sf) {
  bb <- sf::st_bbox(bnd_sf)
  f  <- ca_fig1$map_frame
  w  <- unname((bb[["xmax"]] - bb[["xmin"]]) / (1 - f[["left"]] - f[["right"]]))
  h  <- unname((bb[["ymax"]] - bb[["ymin"]]) / (1 - f[["top"]] - f[["bottom"]]))
  x0 <- bb[["xmin"]] - f[["left"]] * w
  y0 <- bb[["ymin"]] - f[["bottom"]] * h
  list(xlim = c(x0, x0 + w), ylim = c(y0, y0 + h), w = w, h = h,
       x0 = x0, y0 = y0)
}

# Figure output directory. Defaults under the archive so a run on Alpine has
# somewhere to write, and is overridden by CA_FIG_DIR on the laptop.
ca_fig_dir <- function(...) {
  root <- if (nzchar(ca_fig$dir)) ca_fig$dir else ca_archive_out("figures")
  dir <- file.path(root, ...)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

ca_fig_path <- function(name, ext = ca_fig$device) {
  file.path(ca_fig_dir(), sprintf("%s.%s", name, ext))
}

# WHERE A STAGE 11 INPUT LIVES. On Alpine everything sits where its stage put
# it. On the laptop the files have been copied down into one directory, which
# is meta/ by convention. Searched rather than assumed, so the same script runs
# on both machines without an edit.
ca_fig_find <- function(file) {
  # CA_DATA is the laptop case, where the files have been copied into one
  # directory rather than into the archive layout. Searched first, because on
  # the laptop the archive paths resolve to folders that may be empty. Both the
  # flat form and a meta subfolder are tried, so either copy layout works.
  data_dir <- Sys.getenv("CA_DATA", unset = "")
  cand <- c(if (nzchar(data_dir)) file.path(data_dir, file),
            if (nzchar(data_dir)) file.path(data_dir, "meta", file),
            ca_meta(file),
            file.path(ca_match_dir(), file),
            file.path(ca_did_est_dir(), file),
            file)
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("Cannot find ", file, ". Looked in:\n  ",
         paste(cand, collapse = "\n  "))
  }
  hit[1]
}

# The window file, resolved once for every caller. The revised copy wins where
# it exists, so stage 11 and stage 10d read the same windows without either
# naming a file. Returns a path, never a frame, so the caching in ca_fig.R and
# the one-off reads in 11d both work off it.
ca_fig_caps_path <- function() {
  f <- ca_fig$caps_file
  if (!is.null(f) && nzchar(f)) return(ca_fig_find(f))
  p <- tryCatch(ca_fig_find("fig_caps_rev.csv"), error = function(e) NA_character_)
  if (!is.na(p)) return(p)
  ca_fig_find("fig_caps.csv")
}

# The window file, read and normalised. Excel drops the leading zero from the
# ecoregion codes when a hand-edited copy is saved, so "01" returns as "1" and
# every merge on ecoregion silently loses five of the six regions. Klamath
# survives because 78 has no zero to lose. Padding once here is the only place
# it can be fixed for both stage 10d and stage 11.
ca_fig_caps_read <- function() {
  p <- ca_fig_caps_path()
  k <- utils::read.csv(p, colClasses = c(stratum = "character",
                                         ecoregion = "character"),
                       stringsAsFactors = FALSE)
  k$stratum[is.na(k$stratum)] <- ""
  e <- trimws(ifelse(is.na(k$ecoregion), "", k$ecoregion))
  k$ecoregion <- ifelse(nzchar(e),
                        formatC(as.integer(e), width = 2, flag = "0"), "")
  attr(k, "path") <- p
  k
}

# Ecoregion code to label, with the empty stage 10 code mapping to Statewide.
# Anything outside the six returns NA, which is how the figures drop the
# ecoregions that were estimated but are not shown.
ca_fig_ecoregion <- function(code) {
  code <- as.character(code)
  code[is.na(code)] <- ""
  out <- unname(ca_fig$ecoregions[code])
  out[code == ""] <- "Statewide"
  out
}

# The axis window for one arm and outcome, as c(lo, hi).
ca_fig_cap <- function(arm, outcome) {
  k <- ca_fig$caps[[arm]]
  if (is.null(k)) stop("No cap defined for arm: ", arm)
  if (!outcome %in% names(k$hi)) {
    stop("No cap defined for ", arm, " and ", outcome)
  }
  c(lo = unname(k$lo[[outcome]]), hi = unname(k$hi[[outcome]]))
}