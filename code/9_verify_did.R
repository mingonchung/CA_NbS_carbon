### CA carbon revision pipeline
### Stage 9 verifier
###
### Nine sections, one per array task, following the stage 3 to 6 convention.
### Sections run in separate processes rather than sequentially, because Arrow
### allocates outside the R heap and a long single process accumulates it.
###
###   1  units against the matched parquet and the stage 6 groups
###   2  pair structure in the units files
###   3  gvar recomputed independently from stage 6
###   4  manifest completeness and per-outcome year coverage
###   5  pair integrity inside each panel
###   6  cohort viability inside each panel
###   7  outcome distributions and the cross-product AGB comparison
###   8  sentinel residue at the integer type limits
###   9  control reuse and cohort concentration
###
### Every section writes one CSV to meta/ and prints a PASS or FAIL line. A
### FAIL is a statement about the data, not always a bug, and sections 5, 6,
### and 7 are expected to report non-zero counts that stage 10 filters rather
### than errors to fix here.
###
### Usage
###   Rscript 9_verify_did.R          run every section in one process
###   sbatch --array=1-9 9_verify_did.sub

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

SECTION <- ca_task_id()
SECTIONS <- if (is.na(SECTION)) 1:9 else SECTION

PAIRINGS <- ca_pairings$pairing
TASKS    <- ca_did_tasks()

pass <- function(section, ok, msg) {
  ca_log("[", section, "] ", if (ok) "PASS" else "FAIL", "  ", msg)
  invisible(ok)
}

wr <- function(d, name) {
  p <- ca_meta(paste0("verify_did_", name, ".csv"))
  utils::write.csv(d, p, row.names = FALSE)
  ca_log("  wrote meta/", basename(p))
  invisible(p)
}

read_units <- function(pairing, cols = NULL) {
  p <- ca_did_units_path(pairing)
  if (!file.exists(p)) stop("Missing units file: ", p)
  ds <- arrow::open_dataset(p)
  if (is.null(cols)) cols <- names(ds)
  as.data.frame(dplyr::collect(dplyr::select(ds, dplyr::all_of(
    intersect(cols, names(ds))))))
}

read_panel <- function(pairing, outcome, cols) {
  p <- ca_did_panel_path(pairing, outcome)
  if (!file.exists(p)) return(NULL)
  ds <- arrow::open_dataset(p)
  as.data.frame(dplyr::collect(dplyr::select(ds, dplyr::all_of(
    intersect(cols, names(ds))))))
}


# ---------------------------------------------------------------------------
# 1  UNITS AGAINST THEIR SOURCES
# ---------------------------------------------------------------------------
# The units file is a join of two things it does not own. Every pixel must come
# from the matched parquet stage 8 wrote for the selected specification, and
# every attribute must reproduce from the stage 6 groups file.

if (1L %in% SECTIONS) {

  out <- lapply(PAIRINGS, function(pr) {

    sel <- ca_match_selected(pr)
    u <- read_units(pr, c("unit_id", "pair_id", "pixel_id", "lulc", "treat",
                          "gvar", "cohort_year", "pa_est_year", "stratum",
                          "analysis_group", "ecoregion_l3", "arm", "config"))

    pat <- sprintf("^matched_%s_%s_%s(_lev[0-9]+-[0-9]+)?\\.parquet$",
                   pr, sel$arm, sel$config)
    f <- list.files(ca_match_dir(pr), pattern = pat, full.names = TRUE)
    m <- do.call(rbind, lapply(f, function(p) {
      as.data.frame(dplyr::collect(dplyr::select(
        arrow::open_dataset(p), pixel_id, treat)))
    }))

    target_ids <- unique(as.integer(u$pixel_id))
    
    g <- do.call(rbind, lapply(ca_lulc$label, function(cl) {
      as.data.frame(dplyr::collect(dplyr::select(dplyr::filter(
        arrow::open_dataset(ca_group_path(cl)),
        pixel_id %in% target_ids),
        pixel_id, cohort_year, pa_est_year, analysis_group)))
    }))
    
    k <- match(u$pixel_id, g$pixel_id)

    tr <- u$treat == 1L
    data.frame(
      pairing = pr, arm = sel$arm, config = sel$config,
      n_units = nrow(u), n_matched_rows = nrow(m),
      rows_agree = nrow(u) == nrow(m),
      unit_id_unique = !anyDuplicated(u$unit_id),
      arm_agrees = all(as.character(u$arm) == sel$arm),
      config_agrees = all(as.character(u$config) == sel$config),
      pixels_absent_from_stage6 = sum(is.na(k)),
      cohort_year_mismatch =
        sum(tr & !is.na(k) & u$cohort_year != g$cohort_year[k], na.rm = TRUE),
      pa_est_mismatch =
        sum(!is.na(k) & !is.na(u$pa_est_year) &
              u$pa_est_year != g$pa_est_year[k], na.rm = TRUE),
      n_na_ecoregion = sum(is.na(u$ecoregion_l3)),
      n_na_gvar = sum(is.na(u$gvar)),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d)
  wr(d, "01_units_sources")
  pass(1, all(d$rows_agree) && all(d$unit_id_unique) &&
          !sum(d$pixels_absent_from_stage6) && !sum(d$n_na_gvar) &&
          !sum(d$cohort_year_mismatch),
       "units reproduce from the matched parquet and stage 6")
}


# ---------------------------------------------------------------------------
# 2  PAIR STRUCTURE
# ---------------------------------------------------------------------------
# 1:1 without replacement means every pair_id holds exactly one treated and one
# control row. A singleton is two level ranges stacked twice or a subclass that
# lost its partner, and both are silent downstream.

if (2L %in% SECTIONS) {

  out <- lapply(PAIRINGS, function(pr) {
    u <- read_units(pr, c("pair_id", "treat", "pixel_id", "stratum", "gvar"))
    n <- tapply(u$treat, u$pair_id, length)
    t1 <- tapply(u$treat, u$pair_id, function(z) sum(z == 1L))
    # stratum and gvar are pair-level by construction, so a pair carrying two
    # values of either means the propagation at 9a did not take.
    ns <- tapply(as.character(u$stratum), u$pair_id,
                 function(z) length(unique(z)))
    ng <- tapply(u$gvar, u$pair_id, function(z) length(unique(z[z > 0L])))
    data.frame(
      pairing = pr, n_pairs = length(n),
      pairs_not_two = sum(n != 2L),
      pairs_not_one_treated = sum(t1 != 1L),
      pairs_split_stratum = sum(ns > 1L),
      pairs_split_gvar = sum(ng > 1L),
      duplicate_treated_pixels =
        sum(duplicated(u$pixel_id[u$treat == 1L])),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d)
  wr(d, "02_pair_structure")
  pass(2, !sum(d$pairs_not_two) && !sum(d$pairs_not_one_treated) &&
          !sum(d$pairs_split_stratum) && !sum(d$pairs_split_gvar) &&
          !sum(d$duplicate_treated_pixels),
       "every pair holds one treated and one control with a shared stratum")
}


# ---------------------------------------------------------------------------
# 3  gvar RECOMPUTED
# ---------------------------------------------------------------------------
# The rule lives in ca_did_gvar(). This applies it again to the stage 6 columns
# and compares, so a config edit that changes the rule without a 9a rerun
# surfaces here rather than in the estimates.

if (3L %in% SECTIONS) {

  out <- lapply(PAIRINGS, function(pr) {
    u <- read_units(pr, c("treat", "gvar", "cohort_year", "pa_est_year",
                          "pair_id"))
    # pa_est_year on the units file is the pixel's own. The rule takes the
    # treated member's, so it is rebuilt at pair level before comparison.
    tr <- u[u$treat == 1L, , drop = FALSE]
    j <- match(u$pair_id, tr$pair_id)
    g2 <- ca_did_gvar(pr, u$treat, tr$cohort_year[j], tr$pa_est_year[j])
    rng <- ca_did_gvar_range(pr)
    t1 <- u$treat == 1L
    data.frame(
      pairing = pr,
      n_units = nrow(u),
      gvar_mismatch = sum(g2 != u$gvar, na.rm = TRUE),
      control_gvar_nonzero = sum(u$gvar[!t1] != 0L),
      treated_gvar_below_range = sum(u$gvar[t1] < rng[1]),
      treated_gvar_above_range = sum(u$gvar[t1] > rng[2]),
      n_cohorts = length(unique(u$gvar[t1])),
      always_treated_pairs =
        sum(u$gvar[t1] <= min(ca_did$panel_years)),
      pct_always_treated =
        round(100 * mean(u$gvar[t1] <= min(ca_did$panel_years)), 2),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d)
  wr(d, "03_gvar")
  pass(3, !sum(d$gvar_mismatch) && !sum(d$control_gvar_nonzero) &&
          !sum(d$treated_gvar_below_range) && !sum(d$treated_gvar_above_range),
       "gvar reproduces from stage 6 and sits inside the admissible range")
}


# ---------------------------------------------------------------------------
# 4  MANIFEST AND YEAR COVERAGE
# ---------------------------------------------------------------------------
# One row per task, the file it names present, and the years written equal to
# ca_did_years() unless a shard was absent, in which case 9b trimmed and the
# manifest records the shortfall.

if (4L %in% SECTIONS) {

  mp <- ca_did_manifest_path()
  if (!file.exists(mp)) stop("No manifest at ", mp)
  man <- utils::read.csv(mp, stringsAsFactors = FALSE)

  key_task <- paste(TASKS$pairing, TASKS$outcome)
  key_man  <- paste(man$pairing, man$outcome)

  chk <- lapply(seq_len(nrow(TASKS)), function(i) {
    pr <- TASKS$pairing[i]; oc <- TASKS$outcome[i]
    j <- match(paste(pr, oc), key_man)
    p <- ca_did_panel_path(pr, oc)
    exp_y <- ca_did_years(oc)
    yrs <- if (file.exists(p)) {
      sort(unique(as.data.frame(dplyr::collect(dplyr::distinct(dplyr::select(
        arrow::open_dataset(p), year))))$year))
    } else integer(0)
    data.frame(
      pairing = pr, outcome = oc,
      in_manifest = !is.na(j),
      file_present = file.exists(p),
      size_mb = if (file.exists(p)) round(file.size(p) / 1e6, 1) else NA_real_,
      years_expected = length(exp_y),
      years_present = length(yrs),
      years_missing = paste(setdiff(exp_y, yrs), collapse = ";"),
      manifest_rows = if (is.na(j)) NA_integer_ else man$n_rows[j],
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, chk)
  print(d[!d$file_present | d$years_present != d$years_expected, ])
  wr(d, "04_manifest")
  pass(4, all(d$in_manifest) && all(d$file_present) &&
          all(d$years_present == d$years_expected) &&
          nrow(man) == nrow(TASKS),
       paste0(sum(d$file_present), " of ", nrow(TASKS),
              " panels present and complete"))
}


# ---------------------------------------------------------------------------
# 5  PAIR INTEGRITY INSIDE EACH PANEL
# ---------------------------------------------------------------------------
# 9b drops rows with no value, so a pair whose treated or control member has no
# coverage for this outcome survives as a half pair. It is no longer a matched
# comparison and ca_did$complete_pairs_only removes it at stage 10. This counts
# what that filter will take, per pairing, outcome, and class.
#
# Expect zero to negligible everywhere except LEMMA, whose static mask is class
# dependent and reaches double digits on the fire pairings.

if (5L %in% SECTIONS) {

  out <- lapply(seq_len(nrow(TASKS)), function(i) {
    pr <- TASKS$pairing[i]; oc <- TASKS$outcome[i]
    d <- read_panel(pr, oc, c("pair_id", "treat", "lulc"))
    if (is.null(d)) return(NULL)
    u <- unique(d[, c("pair_id", "treat", "lulc")])
    sides <- tapply(u$treat, u$pair_id, function(z) length(unique(z)))
    n_pair <- length(sides)
    n_half <- sum(sides < 2L)
    by_cl <- do.call(rbind, lapply(split(u, u$lulc), function(z) {
      s <- tapply(z$treat, z$pair_id, function(w) length(unique(w)))
      data.frame(lulc = as.character(z$lulc[1]), n_pairs = length(s),
                 n_half = sum(s < 2L), stringsAsFactors = FALSE)
    }))
    data.frame(
      pairing = pr, outcome = oc,
      n_pairs = n_pair, n_half_pairs = n_half,
      pct_half = round(100 * n_half / max(n_pair, 1), 3),
      worst_class = if (is.null(by_cl)) NA_character_ else
        by_cl$lulc[which.max(by_cl$n_half / pmax(by_cl$n_pairs, 1))],
      worst_class_pct = if (is.null(by_cl)) NA_real_ else
        round(100 * max(by_cl$n_half / pmax(by_cl$n_pairs, 1)), 2),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d[order(-d$pct_half), ][1:12, ])
  wr(d, "05_pair_integrity")
  # Reported, not acted on. ca_did$complete_pairs_only is FALSE, and the count
  # is a LEMMA data note rather than a deletion rule. The threshold below only
  # decides whether the note is needed for products other than LEMMA.
  non_lemma <- d[d$outcome != "agb_lemma", ]
  pass(5, max(non_lemma$pct_half) < 1,
       paste0("outside LEMMA the worst half-pair share is ",
              round(max(non_lemma$pct_half), 3), " percent. LEMMA reaches ",
              round(max(d$pct_half[d$outcome == "agb_lemma"]), 2),
              " percent, which is its static mask and is reported, not removed"))
}


# ---------------------------------------------------------------------------
# 6  COHORT VIABILITY INSIDE EACH PANEL
# ---------------------------------------------------------------------------
# A cohort whose gvar exceeds the last year of the outcome has no
# post-treatment observation. It cannot contribute an ATT and, carrying a
# non-zero gvar, it is not in the never-treated control group either.
# ca_did$viable_cohorts_only removes it at stage 10.
#
# This is the number that decides how much of the fire and thinning treated
# sample the eMapR and LEMMA panels can actually speak to, and it belongs in
# Table S3 rather than in a log.

if (6L %in% SECTIONS) {

  out <- lapply(seq_len(nrow(TASKS)), function(i) {
    pr <- TASKS$pairing[i]; oc <- TASKS$outcome[i]
    d <- read_panel(pr, oc, c("pair_id", "treat", "gvar", "year"))
    if (is.null(d)) return(NULL)
    last <- max(d$year)
    tr <- unique(d[d$treat == 1L, c("pair_id", "gvar")])
    ok <- tr$gvar <= last
    # post-treatment years available to the cohorts that do survive
    post <- last - tr$gvar[ok] + 1L
    data.frame(
      pairing = pr, outcome = oc, year_last = last,
      n_treated_pairs = nrow(tr),
      n_viable_pairs = sum(ok),
      pct_viable = round(100 * mean(ok), 2),
      n_cohorts = length(unique(tr$gvar)),
      n_viable_cohorts = length(unique(tr$gvar[ok])),
      median_post_years = if (any(ok)) stats::median(post) else NA_real_,
      min_post_years = if (any(ok)) min(post) else NA_real_,
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d[order(d$pct_viable), ][1:12, ])
  wr(d, "06_cohort_viability")
  pass(6, TRUE,
       paste0("lowest viable share ", round(min(d$pct_viable), 2),
              " percent, at ", d$pairing[which.min(d$pct_viable)], " ",
              d$outcome[which.min(d$pct_viable)],
              ". Removed at stage 10 by viable_cohorts_only"))
}


# ---------------------------------------------------------------------------
# 7  OUTCOME DISTRIBUTIONS
# ---------------------------------------------------------------------------
# Units are gC m-2 for stocks and gC m-2 yr-1 for fluxes after ca_convert().
# The three biomass products should agree to within a product difference rather
# than a factor, and the stage 4 verification put their medians at 8,084,
# 8,811, and 7,573 gC m-2 across the whole grid.
#
# The upper tail is the thing to look at. Carbon_AGB is capped by the producer
# at 800 t/ha, which is 37,600 gC m-2 after conversion, and eMapR lands near
# it. LEMMA carries no cap and no clipping rule, so its maximum is unbounded by
# construction and needs reading against what California forest can hold.

if (7L %in% SECTIONS) {

  out <- lapply(seq_len(nrow(TASKS)), function(i) {
    pr <- TASKS$pairing[i]; oc <- TASKS$outcome[i]
    d <- read_panel(pr, oc, c("value", "treat", "lulc"))
    if (is.null(d)) return(NULL)
    v <- d$value
    q <- stats::quantile(v, c(0, .001, .01, .5, .99, .999, 1), na.rm = TRUE)
    data.frame(
      pairing = pr, outcome = oc, n = length(v),
      min = round(q[[1]], 2), p001 = round(q[[2]], 2), p01 = round(q[[3]], 2),
      median = round(q[[4]], 2), p99 = round(q[[5]], 2),
      p999 = round(q[[6]], 2), max = round(q[[7]], 2),
      n_negative = sum(v < 0), n_zero = sum(v == 0),
      signed_expected = ca_outcomes$signed[match(oc, ca_outcomes$outcome)],
      # A logged outcome with a non-positive value returns NaN or -Inf and R
      # drops the row without comment, so the count that log() would lose is
      # reported here rather than discovered at stage 10.
      log_loss = if (isTRUE(ca_outcomes$log[match(oc, ca_outcomes$outcome)]))
        (if (identical(ca_log_transform(oc), "log1p")) sum(v <= -1) else
          sum(v <= 0)) else NA_integer_,
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d)
  wr(d, "07_distributions")

  agb <- d[d$outcome %in% c("agb_almanac", "agb_emapr", "agb_lemma"), ]
  sp <- tapply(agb$max, agb$outcome, max)
  ca_log("  AGB maxima by product, gC m-2: ",
         paste(names(sp), round(sp), collapse = "  "))
  bad_sign <- d$n_negative > 0 & !d$signed_expected
  pass(7, !any(bad_sign) && !sum(d$log_loss, na.rm = TRUE),
       "no unsigned outcome carries negative values and no logged outcome loses rows")
}


# ---------------------------------------------------------------------------
# 8  SENTINEL RESIDUE
# ---------------------------------------------------------------------------
# ca_sentinel_rules masks the type limits it knows about before conversion.
# NBP is declared at +32767 only, and the observed minima sit within tens of
# units of the INT2S floor, so the negative limit is checked directly against
# the stored extracts rather than inferred from the panel.

if (8L %in% SECTIONS) {

  lims <- c(32767, -32767, 32768, -32768)
  layers <- unique(ca_outcomes$layer[ca_outcomes$signed])

  out <- do.call(rbind, lapply(layers, function(ly) {
    yrs <- ca_layer_years(ly)
    yrs <- yrs[seq(1, length(yrs), length.out = min(6, length(yrs)))]
    do.call(rbind, lapply(ca_lulc$label, function(cl) {
      do.call(rbind, lapply(round(yrs), function(y) {
        p <- ca_extract_path(ca_extract_store(ly), ly, y, cl)
        if (!file.exists(p)) return(NULL)
        v <- arrow::read_parquet(p, col_select = "value")$value
        z <- vapply(lims, function(L) sum(!is.na(v) & abs(v - L) < 1e-6),
                    numeric(1))
        data.frame(layer = ly, lulc = cl, year = y, n = length(v),
                   at_p32767 = z[1], at_m32767 = z[2],
                   at_p32768 = z[3], at_m32768 = z[4],
                   observed_min = round(min(v, na.rm = TRUE), 2),
                   observed_max = round(max(v, na.rm = TRUE), 2),
                   stringsAsFactors = FALSE)
      }))
    }))
  }))
  print(out)
  wr(out, "08_sentinel_residue")
  hit <- sum(out$at_m32767) + sum(out$at_m32768)
  pass(8, hit == 0,
       paste0(hit, " stored values sit at the negative INT2S limit. ",
              "Non-zero means ca_sentinel_rules needs the negative limit ",
              "added for that layer and stage 9b rerun."))
}


# ---------------------------------------------------------------------------
# 9  CONTROL REUSE AND COHORT CONCENTRATION
# ---------------------------------------------------------------------------
# Reuse decides whether pair-level clustering is sufficient. A control pixel
# drawn by two cohorts sits in two clusters and its two rows are the same
# forest, which pair clustering cannot see. At a ratio near 1 there is nothing
# to correct and vcov = ~ pair_id stands unqualified.
#
# Cohort concentration is the companion number. A cohort holding one pair
# contributes a group-time ATT estimated from two units, which is what made
# ecdf_max unusable at stage 8b and will widen intervals at stage 10.

if (9L %in% SECTIONS) {

  out <- lapply(PAIRINGS, function(pr) {
    u <- read_units(pr, c("pixel_id", "treat", "gvar", "pair_id", "stratum"))
    ct <- u[u$treat == 0L, , drop = FALSE]
    tr <- u[u$treat == 1L, , drop = FALSE]
    rep_n <- table(ct$pixel_id)
    coh <- table(tr$gvar)
    data.frame(
      pairing = pr,
      n_control_rows = nrow(ct),
      n_control_pixels = length(rep_n),
      control_reuse = round(nrow(ct) / length(rep_n), 4),
      pct_controls_reused = round(100 * mean(rep_n > 1L), 2),
      max_reuse = max(rep_n),
      n_cohorts = length(coh),
      cohorts_under_10_pairs = sum(coh < 10L),
      pairs_in_those = sum(coh[coh < 10L]),
      smallest_cohort = min(coh),
      largest_cohort_share = round(100 * max(coh) / sum(coh), 2),
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, out)
  print(d)
  wr(d, "09_reuse_cohorts")
  pass(9, max(d$control_reuse) < ca_did$cluster_reuse_threshold,
       paste0("highest control reuse ", max(d$control_reuse), " against a ",
              ca_did$cluster_reuse_threshold, " threshold. Above it, report ",
              "two-way clustering on pair and pixel as an SI line."))
}

ca_stamp(sprintf("9_verify_did_%s", paste(SECTIONS, collapse = "_")))
ca_log("Stage 9 verification section(s) ", paste(SECTIONS, collapse = ", "),
       " complete.")
