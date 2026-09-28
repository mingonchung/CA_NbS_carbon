### CA carbon revision pipeline
### Stage 3c. Probe the management, disturbance, and ownership layers
###
### Runs before any 3c extraction is written. Reports what is actually in each
### file rather than what ca_vector_layers assumes, because shapefile attribute
### names are capped at ten characters and differ between vintages.
###
### Three outputs in meta/.
###   probe_3c_layers_<date>.csv    one row per layer: features, CRS, geometry,
###                                 date field coverage against the study window
###   probe_3c_fields_<date>.csv    one row per field: name, type, missing,
###                                 distinct count, example values
###   probe_3c_knight_<date>.csv    every distinct activity string in the four
###                                 thinning layers, matched against Knight
###                                 Tables S4 and S5
###
### The Knight report is the one that matters. An unmatched record is worse than
### a misclassified one, because it silently leaves a disturbed pixel eligible
### for the undisturbed control pool.
###
### Each layer is wrapped so that one bad file cannot lose the whole run, and
### the CSVs are written on exit whether or not the script completes. The
### previous version lost all output when a single layer failed.
###
### Usage
###   Rscript 3c_probe_vectors.R

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
ca_require(c("sf", "terra"))

suppressPackageStartupMessages({
  library(sf)
  library(terra)
})

ca_stamp("stage3c_probe")

stamp <- format(Sys.Date(), "%Y%m%d")
ref_crs <- ca_ref_crs()

layer_rows <- list()
field_rows <- list()
knight_rows <- list()

# ---------------------------------------------------------------------------
# WRITE ON EXIT
# ---------------------------------------------------------------------------
# Registered before any work starts, so a failure part way through still leaves
# the rows gathered up to that point on disk.

bind <- function(lst) {
  lst <- lst[vapply(lst, function(d) is.data.frame(d) && nrow(d) > 0, logical(1))]
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
}

write_all <- function() {
  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  # Referenced directly rather than through get(). The lists are rebuilt in the
  # global environment as the loop runs, and lexical scoping means this function
  # picks up their current values at call time.
  objs <- list(layers = layer_rows, fields = field_rows, knight = knight_rows)
  for (nm in names(objs)) {
    obj <- bind(objs[[nm]])
    if (is.null(obj)) next
    write.csv(obj, ca_meta(paste0("probe_3c_", nm, "_", stamp, ".csv")),
              row.names = FALSE)
    ca_log("wrote probe_3c_", nm, "_", stamp, ".csv  (", nrow(obj), " rows)")
  }
}
# on.exit() is a no-op at top level in Rscript, which is why an earlier run
# produced no CSVs despite completing. write_all() is called after every layer
# instead, so output exists from the first layer onward whatever happens next.

# ---------------------------------------------------------------------------
# 3C.1  KNIGHT CROSSWALK
# ---------------------------------------------------------------------------

if (!file.exists(ca_knight_csv)) {
  stop("Knight crosswalk not found at ", ca_knight_csv)
}

knight <- read.csv(ca_knight_csv, stringsAsFactors = FALSE)
knight$activity_key <- ca_norm_activity(knight$activity)

req <- c("table", "activity", "intensity", "method", "record_class")
if (!all(req %in% names(knight))) {
  stop("Crosswalk is missing columns: ",
       paste(setdiff(req, names(knight)), collapse = ", "))
}
if (!"override_reason" %in% names(knight)) knight$override_reason <- ""

ca_log("Knight crosswalk: ", nrow(knight), " activities.  S4 ",
       sum(knight$table == "S4"), ", S5 ", sum(knight$table == "S5"))
ca_log("Record classes: ",
       paste(names(table(knight$record_class)),
             as.integer(table(knight$record_class)),
             sep = "=", collapse = "  "))

# Alias table. Activity strings observed in the shapefiles that are absent from
# the published tables, resolved either as a wording variant of a published
# activity or as an assignment with a stated rationale. Kept separate from the
# verbatim transcription so provenance stays clean.
alias <- if (file.exists(ca_alias_csv)) {
  read.csv(ca_alias_csv, stringsAsFactors = FALSE)
} else {
  ca_log("No alias table at ", ca_alias_csv)
  data.frame(table = character(0), observed_key = character(0),
             intensity = character(0), method = character(0),
             record_class = character(0), resolution = character(0),
             stringsAsFactors = FALSE)
}
if (nrow(alias)) {
  # The missing-value token must not go through ca_norm_activity, which strips
  # the angle brackets and yields MISSING, while the layer side keeps <MISSING>.
  # That mismatch left four THP polygons unresolved.
  # The missing-value token must not go through ca_norm_activity, which strips
  # the angle brackets and yields MISSING while the layer side keeps <MISSING>.
  # An empty or NA observed_activity is also treated as the token, so an older
  # alias file that encoded it as a blank still resolves.
  is_miss <- is.na(alias$observed_activity) |
    !nzchar(trimws(alias$observed_activity)) |
    alias$observed_activity == "<MISSING>"
  alias$observed_key <- ifelse(is_miss, "<MISSING>",
                               ca_norm_activity(alias$observed_activity))
  ca_log("Alias table: ", nrow(alias), " resolutions, ",
         sum(alias$resolution == "alias"), " aliases, ",
         sum(alias$resolution == "assigned"), " assignments")
}

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Nearest published activity by edit distance, for anything that failed to
# match. Guarded on two counts: rows whose key is NA or empty are skipped, and
# a row whose distances are all NA yields NA rather than a zero-length result.
# apply(d, 1, which.min) returns a list in that case, which is what crashed the
# previous run on THP, where SILVI_1 carries missing values.
nearest_activity <- function(keys, kt) {
  out <- rep(NA_character_, length(keys))
  ok <- which(!is.na(keys) & nzchar(keys))
  if (!length(ok) || !nrow(kt)) return(out)
  d <- utils::adist(keys[ok], kt$activity_key, ignore.case = TRUE)
  out[ok] <- vapply(seq_len(nrow(d)), function(r) {
    w <- which.min(d[r, ])
    if (length(w) == 1L) kt$knight_activity[w] else NA_character_
  }, character(1))
  out
}

# Year from a date field of unknown class. Handles Date, POSIXct, character,
# and numeric encodings without assuming any of them.
year_of <- function(x) {
  if (inherits(x, "Date") || inherits(x, "POSIXt")) {
    return(as.integer(format(x, "%Y")))
  }
  y <- suppressWarnings(
    as.integer(substr(gsub("[^0-9]", "", as.character(x)), 1, 4))
  )
  y[!is.na(y) & (y < 1900L | y > 2100L)] <- NA_integer_
  y
}

# ---------------------------------------------------------------------------
# 3C.2  VECTOR LAYERS
# ---------------------------------------------------------------------------

for (i in seq_len(nrow(ca_vector_layers))) {

  lay <- ca_vector_layers$layer[i]

  res <- try({

    path <- ca_find_input(ca_vector_layers$subdir[i], ca_vector_layers$file[i])
    loc <- attr(path, "location")

    if (is.na(path)) {

      layer_rows[[lay]] <- data.frame(
        layer = lay, status = "file missing", location = "missing",
        stringsAsFactors = FALSE)
      ca_log(lay, "  MISSING on both scratch and /projects")

    } else {

    # Winding-order warnings on NTO and THP are autocorrected by GDAL on read.
    # st_make_valid is applied anyway, because invalid rings would otherwise
    # surface later during rasterisation, where the failure is harder to trace.
    v <- suppressWarnings(st_read(as.character(path), quiet = TRUE))
    n_invalid_before <- sum(!st_is_valid(v))
    if (n_invalid_before > 0) v <- st_make_valid(v)

    geom <- as.character(unique(st_geometry_type(v)))
    att <- st_drop_geometry(v)

    lrow <- data.frame(
      layer      = lay,
      status     = "ok",
      location   = loc,
      n_features = nrow(v),
      geom_type  = paste(geom, collapse = " "),
      n_fields   = ncol(att),
      same_crs   = terra::same.crs(ref_crs, st_crs(v)$wkt),
      crs_name   = st_crs(v)$Name,
      invalid_before = n_invalid_before,
      invalid_after  = sum(!st_is_valid(v)),
      empty          = sum(st_is_empty(v)),
      stringsAsFactors = FALSE
    )

    ca_log(lay, "  [", loc, "]  features=", nrow(v),
           "  geom=", paste(geom, collapse = "/"),
           "  same_crs=", lrow$same_crs,
           "  invalid=", n_invalid_before)

    # Field inventory.
    field_rows[[lay]] <- data.frame(
      layer    = lay,
      field    = names(att),
      type     = vapply(att, function(x) class(x)[1], character(1)),
      n_miss   = vapply(att, function(x) sum(is.na(x)), integer(1)),
      n_unique = vapply(att, function(x) length(unique(x[!is.na(x)])), integer(1)),
      example  = vapply(att, function(x) {
        u <- unique(x[!is.na(x)])
        if (!length(u)) "" else
          substr(paste(utils::head(u, 3), collapse = " | "), 1, 120)
      }, character(1)),
      is_key = names(att) %in% c(ca_vector_layers$date_field[i],
                                 ca_vector_layers$key_field[i],
                                 ca_vector_layers$code_field[i]),
      stringsAsFactors = FALSE
    )
    rownames(field_rows[[lay]]) <- NULL

    # Plan status. Unlogged, Withdrawn, and Approved plans are not treatment,
    # and they carry the 1899-12-30 placeholder that inflates the out-of-window
    # count. This distribution is what determines the real treatment sample.
    sfield <- ca_vector_layers$status_field[i]
    if (!is.na(sfield) && sfield %in% names(att)) {
      tb <- table(att[[sfield]], useNA = "ifany")
      lrow$n_status_keep <- sum(att[[sfield]] %in%
                                  ca_vector_layers$status_keep[i], na.rm = TRUE)
      ca_log("  ", sfield, ": ",
             paste(names(tb), as.integer(tb), sep = "=", collapse = "  "))
    }

    # Plan or establishment year, distinct from the completion date. The two
    # CAL FIRE files are split on this, so their spans must not overlap.
    yfield <- ca_vector_layers$year_field[i]
    if (!is.na(yfield) && yfield %in% names(att)) {
      yv <- suppressWarnings(as.integer(att[[yfield]]))
      yv <- yv[!is.na(yv) & yv > 1900]
      if (length(yv)) {
        lrow$year_field <- yfield
        lrow$year_field_min <- min(yv)
        lrow$year_field_max <- max(yv)
        ca_log("  ", yfield, ": ", min(yv), " to ", max(yv),
               "  (", length(unique(yv)), " distinct)")
      }
    }

    # Completion date. Only completed records enter the treatment dataset, and
    # only those inside the study window are usable.
    dfield <- ca_vector_layers$date_field[i]
    if (!is.na(dfield)) {
      if (!dfield %in% names(att)) {
        ca_log("  WARNING date field ", dfield, " absent. Fields: ",
               paste(names(att), collapse = " "))
        lrow$date_status <- "field absent"
      } else {
        yr <- year_of(att[[dfield]])
        in_win <- !is.na(yr) & yr >= min(ca_years$study) & yr <= max(ca_years$study)
        lrow$date_field   <- dfield
        lrow$date_class   <- class(att[[dfield]])[1]
        lrow$date_missing <- sum(is.na(att[[dfield]]))
        lrow$year_min     <- if (all(is.na(yr))) NA_integer_ else min(yr, na.rm = TRUE)
        lrow$year_max     <- if (all(is.na(yr))) NA_integer_ else max(yr, na.rm = TRUE)
        lrow$n_in_window  <- sum(in_win)
        lrow$n_out_window <- sum(!in_win)
        ca_log("  ", dfield, ": class=", lrow$date_class,
               "  missing=", lrow$date_missing,
               "  years ", lrow$year_min, " to ", lrow$year_max,
               "  in window ", format(lrow$n_in_window, big.mark = ","),
               "  outside ", format(lrow$n_out_window, big.mark = ","))
      }
    }

    layer_rows[[lay]] <- lrow

    # Knight match report for the four thinning layers.
    kfield <- ca_vector_layers$key_field[i]
    tab <- ca_vector_layers$crosswalk[i]

    if (!is.na(kfield) && !is.na(tab)) {
      if (!kfield %in% names(att)) {
        ca_log("  WARNING key field ", kfield, " absent. Fields: ",
               paste(names(att), collapse = " "))
      } else {

        # Missing values get an explicit token rather than an empty key, so a
        # blank prescription can be resolved through the alias table like any
        # other value instead of falling out as unresolved.
        raw <- as.character(att[[kfield]])
        raw[is.na(raw) | !nzchar(trimws(raw))] <- "<MISSING>"
        cnt <- as.data.frame(table(raw), stringsAsFactors = FALSE)
        names(cnt) <- c("activity_raw", "n_features")

        # Features that also pass the plan-status filter. Approved, Unlogged,
        # and Withdrawn plans are permits, not operations, so the unfiltered
        # count overstates the treatment sample.
        if (!is.na(sfield) && sfield %in% names(att)) {
          done <- att[[sfield]] %in% ca_vector_layers$status_keep[i]
          cnt_done <- as.data.frame(table(raw[done]), stringsAsFactors = FALSE)
          names(cnt_done) <- c("activity_raw", "n_completed")
          cnt <- merge(cnt, cnt_done, by = "activity_raw", all.x = TRUE)
          cnt$n_completed[is.na(cnt$n_completed)] <- 0L
        } else {
          cnt$n_completed <- NA_integer_
        }
        cnt$activity_key <- ifelse(cnt$activity_raw == "<MISSING>", "<MISSING>",
                                   ca_norm_activity(cnt$activity_raw))

        kcols <- c("activity_key", "activity", "intensity", "method",
                   "record_class", "override_reason")
        kt <- knight[knight$table == tab, kcols]
        names(kt)[2] <- "knight_activity"

        m <- merge(cnt, kt, by = "activity_key", all.x = TRUE, sort = FALSE)

        # Second pass. Anything the published table did not resolve is looked up
        # in the alias table. resolved_by records which pass supplied the class,
        # so the response letter can state how many activities came from the
        # published tables and how many were resolved here.
        m$resolved_by <- ifelse(is.na(m$record_class), NA_character_, "knight")
        al <- alias[alias$table == tab, ]
        if (nrow(al)) {
          hit <- match(m$activity_key, al$observed_key)
          fill <- is.na(m$record_class) & !is.na(hit)
          if (any(fill)) {
            m$intensity[fill]    <- al$intensity[hit[fill]]
            m$method[fill]       <- al$method[hit[fill]]
            m$record_class[fill] <- al$record_class[hit[fill]]
            m$resolved_by[fill]  <- al$resolution[hit[fill]]
          }
        }

        m$layer <- lay
        m$crosswalk_table <- tab
        m$matched <- !is.na(m$record_class)
        m$is_treatment <- m$record_class %in% "treatment_LMH"
        m$disqualifies_control <- m$record_class %in%
          ca_record_classes$record_class[ca_record_classes$disqualifies_control]
        m$nearest <- NA_character_
        un <- !m$matched
        if (any(un)) m$nearest[un] <- nearest_activity(m$activity_key[un], kt)

        knight_rows[[lay]] <- m[order(-m$n_features),
                                c("layer", "crosswalk_table", "activity_raw",
                                  "activity_key", "matched", "intensity",
                                  "method", "record_class", "is_treatment",
                                  "disqualifies_control", "n_features",
                                  "n_completed", "resolved_by",
                                  "knight_activity", "override_reason",
                                  "nearest")]

        ca_log("  ", kfield, ": ", nrow(m), " distinct values, ",
               sum(m$resolved_by %in% "knight"), " from Knight, ",
               sum(m$resolved_by %in% c("alias", "assigned")), " from alias, ",
               sum(!m$matched), " unresolved")
        ca_log("    features  treatment=",
               format(sum(m$n_features[m$is_treatment]), big.mark = ","),
               "  disqualify_control=",
               format(sum(m$n_features[m$disqualifies_control &
                                         !m$is_treatment]), big.mark = ","),
               "  non_disturbing=",
               format(sum(m$n_features[m$record_class %in% "non_disturbing"]),
                      big.mark = ","),
               "  unmatched=",
               format(sum(m$n_features[!m$matched]), big.mark = ","),
               "  total=", format(sum(m$n_features), big.mark = ","))
        if (!all(is.na(m$n_completed))) {
          ca_log("    completed only  treatment=",
                 format(sum(m$n_completed[m$is_treatment]), big.mark = ","),
                 "  disqualify_control=",
                 format(sum(m$n_completed[m$disqualifies_control &
                                            !m$is_treatment]), big.mark = ","),
                 "  total=", format(sum(m$n_completed), big.mark = ","))
        }
      }
    }

    rm(v, att)
    gc(verbose = FALSE)
    }
    TRUE
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    ca_log(lay, "  FAILED: ", conditionMessage(attr(res, "condition")))
    layer_rows[[lay]] <- data.frame(
      layer = lay, status = "error",
      note = conditionMessage(attr(res, "condition")),
      stringsAsFactors = FALSE)
  }

  # After the error row, so a failed layer is recorded in the same pass rather
  # than only appearing once the next layer writes.
  write_all()
}

# ---------------------------------------------------------------------------
# 3C.2b  OFFSET START YEARS
# ---------------------------------------------------------------------------
# The offset shapefile has no date field. Cohort years come from the CARB
# issuance workbook, joined on ARB_id. Reported here so an unmatched project is
# caught before stage 5 rather than silently losing a cohort.

off_path <- ca_find_input(ca_vector_layers$subdir[ca_vector_layers$layer == "offset"],
                          ca_vector_layers$file[ca_vector_layers$layer == "offset"])

if (!is.na(off_path) && file.exists(ca_offset_years_csv)) {

  off <- suppressWarnings(st_drop_geometry(st_read(as.character(off_path),
                                                   quiet = TRUE)))
  oy <- read.csv(ca_offset_years_csv, stringsAsFactors = FALSE)

  if (!ca_offset_join_field %in% names(off)) {
    ca_log("OFFSET: join field ", ca_offset_join_field, " absent. Fields: ",
           paste(names(off), collapse = " "))
  } else {
    # Compliance rows only. Early Action rows stay in the lookup for the
    # predecessor check below but must not join, or an Early Action id sharing
    # a number with a compliance id would attach the wrong start year.
    oy_cop <- oy[oy$ea_cop == ca_offset_program, ]

    keep <- intersect(c(ca_offset_join_field, "P_name", "OPO", "Notes"),
                      names(off))
    j <- merge(off[, keep, drop = FALSE],
               oy_cop, by.x = ca_offset_join_field, by.y = "arb_id",
               all.x = TRUE, sort = FALSE)
    names(j)[1] <- "arb_id"
    ca_log("OFFSET: ", nrow(j), " projects, ",
           sum(!is.na(j$start_year)), " matched to CARB, ",
           sum(is.na(j$start_year)), " unmatched")
    if (any(!is.na(j$start_year))) {
      ca_log("  start years ", min(j$start_year, na.rm = TRUE), " to ",
             max(j$start_year, na.rm = TRUE))
      print(table(j$start_year, useNA = "ifany"))
    }
    if (any(is.na(j$start_year))) {
      miss <- j$arb_id[is.na(j$start_year)]
      ca_log("  unmatched ARB ids: ", paste(miss, collapse = " "))
      in_ea <- miss[miss %in% oy$arb_id[oy$ea_cop != ca_offset_program]]
      if (length(in_ea)) {
        ca_log("  of those, present as Early Action: ",
               paste(in_ea, collapse = " "))
      }
    }

    # Predecessor check. A compliance project whose management already changed
    # under an Early Action project has treated years inside its pre-treatment
    # window. The workbook carries no link between the two, so the only signals
    # available here are the free-text Notes field and the registry URLs.
    if ("Notes" %in% names(off)) {
      nt <- off$Notes[!is.na(off$Notes) & nzchar(off$Notes)]
      idx_nt <- which(!is.na(off$Notes) & nzchar(off$Notes))
      if (length(idx_nt)) {
        ca_log("  Notes present on ", length(idx_nt), " project(s), check for ",
               "an Early Action predecessor:")
        for (k in idx_nt) {
          ca_log("    ", off[[ca_offset_join_field]][k], "  ",
                 substr(off$Notes[k], 1, 140))
        }
      }
    }
    write.csv(j, ca_meta(paste0("probe_3c_offset_", stamp, ".csv")),
              row.names = FALSE)
    ca_log("wrote probe_3c_offset_", stamp, ".csv")
  }
} else if (is.na(off_path)) {
  ca_log("OFFSET: shapefile missing")
} else {
  ca_log("OFFSET: no CARB lookup at ", ca_offset_years_csv)
}

# ---------------------------------------------------------------------------
# 3C.2c  CPAD ESTABLISHMENT YEARS
# ---------------------------------------------------------------------------
# YR_EST uses 0 as an unknown marker. Counted here because the count determines
# how much rests on the assumption that unknown means old.

cp_path <- ca_find_input(ca_vector_layers$subdir[ca_vector_layers$layer == "cpad"],
                         ca_vector_layers$file[ca_vector_layers$layer == "cpad"])

if (!is.na(cp_path)) {
  cp <- suppressWarnings(st_drop_geometry(st_read(as.character(cp_path),
                                                   quiet = TRUE)))
  if (all(c("YR_EST", "MNG_AG_LEV") %in% names(cp))) {
    y <- as.numeric(cp$YR_EST)
    ca_log("CPAD YR_EST: unknown(0)=", sum(y == 0, na.rm = TRUE),
           "  pre-1985=", sum(y > 0 & y < 1985, na.rm = TRUE),
           "  1985-1989=", sum(y >= 1985 & y <= 1989, na.rm = TRUE),
           "  1990+=", sum(y >= 1990, na.rm = TRUE),
           "  NA=", sum(is.na(y)))
    ca_log("CPAD acres by establishment class:")
    cls <- cut(y, c(-Inf, 0, 1984, 1989, Inf),
               labels = c("unknown", "pre_1985", "1985_1989", "post_1990"))
    print(tapply(as.numeric(cp$ACRES), cls, sum))
    ca_log("CPAD MNG_AG_LEV:")
    print(table(cp$MNG_AG_LEV))

    role <- ca_cpad_agency$role[match(cp$MNG_AG_LEV, ca_cpad_agency$label)]
    ca_log("CPAD by role, features. pa_treatment enters UP and the P groups, ",
           "protected_other is removed from the control pool only:")
    print(table(role, useNA = "ifany"))
    ca_log("CPAD by role, acres:")
    print(round(tapply(as.numeric(cp$ACRES), role, sum)))
    if (any(is.na(role))) {
      ca_log("  WARNING unmapped MNG_AG_LEV values: ",
             paste(unique(cp$MNG_AG_LEV[is.na(role)]), collapse = " | "))
    }
  }
}

# ---------------------------------------------------------------------------
# 3C.2d  OWNERSHIP LEVELS
# ---------------------------------------------------------------------------
# The private non-protected control is the complement of this layer and CPAD.
# Any level in here that is not actually public would wrongly remove ground
# from the baseline, so the levels are tabulated rather than assumed.

ow_path <- ca_find_input(
  ca_vector_layers$subdir[ca_vector_layers$layer == "ownership"],
  ca_vector_layers$file[ca_vector_layers$layer == "ownership"])

if (!is.na(ow_path)) {
  ow <- suppressWarnings(st_drop_geometry(st_read(as.character(ow_path),
                                                   quiet = TRUE)))
  for (f in intersect(c("Own_Level", "Own_Group"), names(ow))) {
    ca_log("OWNERSHIP ", f, ":")
    print(sort(table(ow[[f]]), decreasing = TRUE))
  }
  if ("Own_Level" %in% names(ow)) {
    orole <- ca_ownership_level$role[match(ow$Own_Level,
                                           ca_ownership_level$label)]
    ca_log("OWNERSHIP by role:")
    print(table(orole, useNA = "ifany"))
    if (any(is.na(orole))) {
      ca_log("  WARNING unmapped Own_Level values: ",
             paste(unique(ow$Own_Level[is.na(orole)]), collapse = " | "))
    }
  }
}

# ---------------------------------------------------------------------------
# 3C.3  MTBS RASTER STACK
# ---------------------------------------------------------------------------

mtbs_dir <- ca_find_input(ca_mtbs_layers$subdir[1])

if (!is.na(mtbs_dir)) {

  ca_log("MTBS directory [", attr(mtbs_dir, "location"), "] ", mtbs_dir)
  files <- list.files(as.character(mtbs_dir),
                      pattern = ca_mtbs_layers$file_pattern[1],
                      full.names = TRUE, recursive = TRUE)

  if (length(files)) {
    yrs <- sort(as.integer(sub(".*mtbs_CA_([0-9]{4})_prj\\.tif$", "\\1",
                               basename(files))))
    r <- rast(files[1])
    ca_log("MTBS: ", length(files), " files, ", min(yrs), " to ", max(yrs),
           "  gaps=", paste(setdiff(seq(min(yrs), max(yrs)), yrs),
                            collapse = " "),
           "  same_crs=", terra::same.crs(ref_crs, crs(r)),
           "  datatype=", datatype(r)[1],
           "  NAflag=", paste(NAflag(r), collapse = " "))
    ca_log("MTBS severity class counts, ", min(yrs), ":")
    fq <- freq(r)
    fq$role <- ca_fire_severity$role[match(fq$value, ca_fire_severity$code)]
    fq$group <- ca_fire_severity$group[match(fq$value, ca_fire_severity$code)]
    print(fq)
    if (any(is.na(fq$role))) {
      ca_log("  WARNING MTBS values with no role in ca_fire_severity: ",
             paste(fq$value[is.na(fq$role)], collapse = " "))
    }

    layer_rows[["MTBS"]] <- data.frame(
      layer = "MTBS", status = "ok",
      location = attr(mtbs_dir, "location"),
      n_features = length(files),
      year_min = min(yrs), year_max = max(yrs),
      same_crs = terra::same.crs(ref_crs, crs(r)),
      stringsAsFactors = FALSE)
  } else {
    ca_log("MTBS: no files matching ", ca_mtbs_layers$file_pattern[1])
  }
} else {
  ca_log("MTBS: directory missing on both scratch and /projects")
}

# ---------------------------------------------------------------------------
# 3C.4  SUMMARY
# ---------------------------------------------------------------------------

kr <- bind(knight_rows)
if (!is.null(kr)) {
  ca_log("Unmatched activity values across all thinning layers: ",
         sum(!kr$matched), " distinct, ",
         format(sum(kr$n_features[!kr$matched]), big.mark = ","), " features")
  if (any(!kr$matched)) {
    print(kr[!kr$matched, c("layer", "activity_raw", "n_features", "nearest")])
  }
}

write_all()
ca_log("Stage 3c probe complete.")
