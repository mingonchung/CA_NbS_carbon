### CA carbon revision pipeline
### Stage 6 verification
###
### Sections
###
###   integrity    grid coverage, uniqueness, completeness, and agreement
###                between the assigned codes and the config declarations
###   ledger       the exclusion sequence reconciled to the grid total, pooled
###                and per class
###   stage5       every pixel carrying a stage 5a record accounted for in the
###                stage 6 output, in both directions
###   groups       group counts, cohort distribution, detection retention and
###                the never-disturbed against wrong-window split
###   control      UNP composition. The verification-log item. UNP is the
###                complement of every excluded status, so a pixel with no
###                ownership record at all falls into it silently
###   spatial      ecoregion and HUC8 distribution per group, replacing
###                6_4_group_sel_ecoregion.R
###   taxonomy     class code tail, how much of the grid sits in codes too rare
###                to interpret
###   feasibility  exact-cell matching feasibility per pairing. The section
###                that decides whether a pairing can be estimated at all,
###                before stage 7 allocates anything to it
###   sample       stage 7. Referential integrity against stage 6, the 150 m
###                separation guarantee tested rather than assumed, retention
###                reconciled to the summary files, and feasibility recomputed
###                on what stage 8 will actually read, spaced treated against
###                whole control
###   all          every section
###
### Usage
###   Rscript 6_verify_groups.R integrity   [class]
###   Rscript 6_verify_groups.R feasibility [class]
###   Rscript 6_verify_groups.R all         [class]

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

if (!CLASS %in% ca_lulc$label) stop("Unknown class: ", CLASS)

STAMP <- format(Sys.time(), "%Y%m%d")
SECTIONS <- c("integrity", "ledger", "stage5", "groups", "control",
              "spatial", "taxonomy", "feasibility", "sample")
if (!section %in% c(SECTIONS, "all")) stop("Unknown section: ", section)

report <- function(df, name) {
  p <- ca_meta(sprintf("verify_stage6_%s_%s_%s.csv", name, CLASS, STAMP))
  utils::write.csv(df, p, row.names = FALSE)
  ca_log("Wrote ", basename(p))
  print(utils::head(as.data.frame(df), 40))
  invisible(p)
}

FAIL <- character(0)
check <- function(ok, msg) {
  if (isTRUE(ok)) {
    ca_log("  PASS  ", msg)
  } else {
    ca_log("  FAIL  ", msg)
    FAIL <<- c(FAIL, msg)
  }
  invisible(ok)
}

read_groups <- function(cols = NULL) {
  p <- ca_group_path(CLASS)
  if (!file.exists(p)) stop("Stage 6a output missing: ", p)
  if (is.null(cols)) {
    as.data.frame(arrow::read_parquet(p))
  } else {
    as.data.frame(arrow::read_parquet(p, col_select = dplyr::all_of(cols)))
  }
}

read_arm <- function(arm, cols) {
  p <- ca_event_path(arm, CLASS)
  if (!file.exists(p)) stop("Missing stage 5a arm: ", p)
  as.data.frame(arrow::read_parquet(p, col_select = dplyr::all_of(cols)))
}

run <- function(x) section == "all" || section == x

# ---------------------------------------------------------------------------
# 1. INTEGRITY
# ---------------------------------------------------------------------------

if (run("integrity")) {
  ca_log("Section integrity, class ", CLASS)
  g <- read_groups(c("pixel_id", "lulc", "class_code", "event_chain",
                     "status_code", "analysis_group", "arm", "stratum",
                     "eligibility", "cohort_year", "pa_est_year",
                     "n_events_effective"))

  part <- file.path(ca_grid$parquet_dir, paste0("lulc=", CLASS))
  gid <- as.integer(as.data.frame(arrow::Scanner$create(
    arrow::open_dataset(part), projection = "pixel_id")$ToTable())$pixel_id)

  check(nrow(g) == length(gid), sprintf(
    "row count matches the frozen grid, %s against %s",
    format(nrow(g), big.mark = ","), format(length(gid), big.mark = ",")))
  check(!anyDuplicated(g$pixel_id), "pixel_id unique")
  check(setequal(g$pixel_id, gid), "pixel_id set identical to the grid")
  check(all(g$lulc == CLASS), "lulc constant and correct")

  # class_code is the complete taxonomy, so a missing value is a defect rather
  # than a category.
  check(!any(is.na(g$class_code)), "class_code complete, no missing values")
  check(!any(is.na(g$status_code)), "status_code complete")
  check(!any(is.na(g$eligibility)), "eligibility complete")
  check(all(g$status_code %in% ca_status_levels$status),
        "every status declared in ca_status_levels")
  check(all(na.omit(g$analysis_group) %in% ca_analysis_groups$group),
        "every group declared in ca_analysis_groups")
  check(all(g$eligibility %in% c(ca_exclusion_order, "treatment", "control",
                                 "unassigned")),
        "every eligibility value declared in ca_exclusion_order")

  # A pixel is either in a group or excluded, never both and never neither.
  assigned <- !is.na(g$analysis_group)
  check(all(g$eligibility[assigned] %in% c("treatment", "control")),
        "assigned pixels carry a treatment or control eligibility")
  check(!any(g$eligibility[!assigned] %in% c("treatment", "control")),
        "unassigned pixels carry no treatment or control eligibility")
  check(sum(g$eligibility == "unassigned") == 0,
        sprintf("no pixel falls through every rule, %s unassigned",
                format(sum(g$eligibility == "unassigned"), big.mark = ",")))

  # Group, status, arm, and stratum must agree with the declaration table
  # rather than with each other by accident.
  i <- match(g$analysis_group, ca_analysis_groups$group)
  check(all(g$status_code[assigned] == ca_analysis_groups$status[i][assigned]),
        "group and status agree with ca_analysis_groups")
  check(all(g$arm[assigned] == ca_analysis_groups$arm[i][assigned]),
        "group and arm agree with ca_analysis_groups")
  check(identical(g$stratum[assigned], ca_analysis_groups$stratum[i][assigned]),
        "group and stratum agree with ca_analysis_groups")

  # Cohort years exist exactly where a cohort is defined.
  treat_arm <- assigned & g$arm %in% c("fire", "thin", "offset")
  check(!any(is.na(g$cohort_year[treat_arm])),
        "every fire, thinning, and offset group pixel carries a cohort year")
  check(all(is.na(g$cohort_year[assigned & g$arm %in% c("pa", "control")])),
        "no cohort year on the PA and undisturbed control groups")
  ch <- g$cohort_year[treat_arm]
  check(all(ch >= min(ca_years$analysis) & ch <= max(ca_years$analysis)),
        sprintf("cohort years inside the analysis window, observed %s to %s",
                min(ch), max(ch)))

  # The undisturbed groups must carry no effective event, by definition.
  und <- assigned & g$analysis_group %in% c("UP", "UNP", "OP", "ONP")
  check(all(g$n_events_effective[und] == 0L),
        "undisturbed and offset groups carry no effective event")

  # Protection is time-varying, so a P-group fire or thinning pixel must have
  # been protected at its own cohort year rather than merely today.
  if (isTRUE(ca_pa_rules$label_at_cohort_year)) {
    pg <- assigned & g$status_code == "P" & g$arm %in% c("fire", "thin")
    late <- pg & !is.na(g$pa_est_year) & g$pa_est_year > g$cohort_year
    check(sum(late) == 0, sprintf(
      "no P-group event pixel protected after its cohort year, found %s",
      format(sum(late), big.mark = ",")))
  }

  # UP carries every establishment era. The calendar-year specification treats
  # a pre-1990 PA as protected in every panel year, so dropping those pixels
  # would discard the ongoing effect of 98 percent of protected forest.
  up <- assigned & g$analysis_group == "UP"
  check(sum(up & !is.na(g$pa_est_year) & g$pa_est_year < min(ca_years$analysis)) > 0,
        "UP retains pixels established before the analysis window")

  report(as.data.frame(table(eligibility = g$eligibility)), "integrity_elig")
  rm(g, gid); ca_gc()
}

# ---------------------------------------------------------------------------
# 2. LEDGER
# ---------------------------------------------------------------------------

if (run("ledger")) {
  ca_log("Section ledger")
  led <- lapply(ca_lulc$label, function(l) {
    p <- ca_meta(sprintf("groups_ledger_%s.csv", l))
    if (!file.exists(p)) return(NULL)
    d <- utils::read.csv(p, stringsAsFactors = FALSE)
    d$lulc <- l
    d
  })
  led <- do.call(rbind, led[!vapply(led, is.null, logical(1))])

  grid_n <- sum(led$remaining[led$rule == "start"])
  removed <- sum(led$removed[led$rule != "start"])
  check(removed == grid_n, sprintf(
    "removals and assignments reconcile to the grid, %s against %s",
    format(removed, big.mark = ","), format(grid_n, big.mark = ",")))

  ord <- led$rule[led$lulc == ca_lulc$label[1]]
  ord <- ord[!ord %in% c("start", "assigned")]
  check(identical(ord, ca_exclusion_order[ca_exclusion_order %in% ord]),
        "rules fire in the order declared in ca_exclusion_order")

  pooled <- led %>%
    filter(rule != "start") %>%
    group_by(rule) %>%
    summarise(removed = sum(removed), .groups = "drop") %>%
    mutate(pct_of_grid = round(100 * removed / grid_n, 2)) %>%
    arrange(match(rule, c(ca_exclusion_order, "assigned")))
  report(pooled, "ledger_pooled")

  report(led[, c("lulc", "step", "rule", "removed", "remaining")],
         "ledger_by_class")
}

# ---------------------------------------------------------------------------
# 3. STAGE 5 RECONCILIATION
# ---------------------------------------------------------------------------

if (run("stage5")) {
  ca_log("Section stage5, class ", CLASS)
  g <- read_groups(c("pixel_id", "n_events", "n_thin", "n_fire",
                     "undated_thin", "status_code", "offset_start",
                     "n_detect"))

  # Direction one. Every pixel carrying a record must show it in stage 6.
  for (arm in c("thin", "fire")) {
    a <- read_arm(arm, c("pixel_id"))
    px <- unique(a$pixel_id)
    px <- px[px %in% g$pixel_id]
    i <- match(px, g$pixel_id)
    has <- if (arm == "thin") {
      g$n_thin[i] > 0L | g$undated_thin[i]
    } else {
      g$n_fire[i] > 0L
    }
    # Fu records arrive through the thinning layer and are counted as fire, so
    # a thinning pixel whose only record resolves to Fu shows up on the fire
    # side. The test allows either.
    if (arm == "thin") has <- has | g$n_fire[i] > 0L
    check(all(has), sprintf(
      "%s arm, every recorded pixel present in stage 6, %s missing",
      arm, format(sum(!has), big.mark = ",")))
    rm(a, px, i, has); ca_gc()
  }

  # Direction two. Stage 6 must not invent events.
  ev_px <- g$pixel_id[g$n_events > 0L]
  src <- unique(c(read_arm("thin", "pixel_id")$pixel_id,
                  read_arm("fire", "pixel_id")$pixel_id))
  check(all(ev_px %in% src), sprintf(
    "no pixel carries an event without a stage 5a record, %s invented",
    format(sum(!ev_px %in% src), big.mark = ",")))
  rm(ev_px, src); ca_gc()

  # Offsets and tribal, exact set agreement.
  off <- read_arm("offset", c("pixel_id", "event_year"))
  off <- unique(off$pixel_id[!is.na(off$event_year)])
  check(setequal(off[off %in% g$pixel_id],
                 g$pixel_id[!is.na(g$offset_start)]),
        "offset pixel set identical to stage 5a")
  # Tribal status comes from two sources. The dedicated tribal layer and the
  # Tribal role inside the ownership layer, which the polygon layer does not
  # fully cover. TR is the union of the two, so the test is containment in both
  # directions rather than equality against one source.
  tr <- unique(read_arm("tribal", "pixel_id")$pixel_id)
  own <- read_arm("own", c("pixel_id", "class_value"))
  own_tr <- unique(own$pixel_id[
    ca_ownership_level$role[
      match(own$class_value, ca_ownership_level$label)] %in% "tribal"])
  rm(own); ca_gc()

  tr_status <- g$pixel_id[g$status_code == "TR"]
  check(all(tr[tr %in% g$pixel_id] %in% tr_status),
        "every tribal-layer pixel carries TR status")
  check(all(own_tr[own_tr %in% g$pixel_id] %in% tr_status),
        "every ownership-tribal pixel carries TR status")
  check(setequal(tr_status,
                 unique(c(tr, own_tr))[unique(c(tr, own_tr)) %in% g$pixel_id]),
        sprintf("TR is exactly the union of the two tribal sources, %s pixels",
                format(length(tr_status), big.mark = ",")))
  rm(off, tr, own_tr, tr_status); ca_gc()

  # Screen, pixel count agreement.
  sc <- as.data.frame(arrow::read_parquet(
    ca_screen_path(CLASS), col_select = dplyr::all_of("pixel_id")))
  n_sc <- length(unique(sc$pixel_id[sc$pixel_id %in% g$pixel_id]))
  check(n_sc == sum(g$n_detect > 0L), sprintf(
    "detected pixel count matches stage 5b, %s against %s",
    format(n_sc, big.mark = ","),
    format(sum(g$n_detect > 0L), big.mark = ",")))
  rm(g, sc); ca_gc()
}

# ---------------------------------------------------------------------------
# 4. GROUPS
# ---------------------------------------------------------------------------

if (run("groups")) {
  ca_log("Section groups, class ", CLASS)
  g <- read_groups(c("analysis_group", "arm", "stratum", "cohort_year",
                     "pa_est_year", "detect_at_event", "status_code"))

  gs <- g %>%
    filter(!is.na(analysis_group)) %>%
    group_by(analysis_group, arm, stratum) %>%
    summarise(n = n(),
              n_cohorts = length(unique(cohort_year[!is.na(cohort_year)])),
              first_cohort = suppressWarnings(min(cohort_year, na.rm = TRUE)),
              last_cohort = suppressWarnings(max(cohort_year, na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(analysis_group)
  gs$first_cohort[!is.finite(gs$first_cohort)] <- NA
  gs$last_cohort[!is.finite(gs$last_cohort)] <- NA
  report(gs, "group_counts")

  # Cohort concentration. A flat per-group cap lets one dominant cohort absorb
  # the allocation and then drive the event study, so the share held by the
  # largest cohort is the number to watch at stage 7.
  cc <- g %>%
    filter(!is.na(analysis_group), !is.na(cohort_year)) %>%
    group_by(analysis_group, cohort_year) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(share = n / sum(n)) %>%
    summarise(n_cohorts = n(),
              max_cohort_share = round(100 * max(share), 1),
              cohorts_under_1pct = sum(share < 0.01),
              .groups = "drop")
  report(cc, "cohort_concentration")

  # PA establishment, the calendar-year specification rather than event time.
  pa <- g %>%
    filter(analysis_group == "UP") %>%
    mutate(era = ifelse(is.na(pa_est_year), "undated",
                 ifelse(pa_est_year < min(ca_years$analysis), "pre-window",
                        "in-window"))) %>%
    group_by(era) %>%
    summarise(n = n(), .groups = "drop")
  report(pa, "pa_establishment")

  miss <- ca_meta(sprintf("groups_detection_miss_%s.csv", CLASS))
  if (file.exists(miss)) {
    d <- utils::read.csv(miss, stringsAsFactors = FALSE)
    ca_log("  detection failures, never disturbed against wrong window")
    print(d)
    # A group whose failures are mostly disturbed off-window has a date problem
    # rather than an over-assignment problem, and the remedy is the window.
    flag <- d$group[d$pct_disturbed_off_window > 50]
    if (length(flag)) {
      ca_log("  NOTE archive dates unreliable for: ",
             paste(flag, collapse = ", "),
             ". Detection is validating the cohort year for these groups.")
    }
  }
  rm(g); ca_gc()
}

# ---------------------------------------------------------------------------
# 5. CONTROL POOL COMPOSITION
# ---------------------------------------------------------------------------

if (run("control")) {
  ca_log("Section control, class ", CLASS)
  g <- read_groups(c("pixel_id", "analysis_group", "status_code"))
  unp <- g$pixel_id[g$analysis_group %in% "UNP"]

  # NP is the complement of every excluded status, so a pixel with no
  # ownership record at all lands in the control pool without ever being
  # examined. That share is the number this section exists for.
  own <- read_arm("own", c("pixel_id", "class_value"))
  own <- own[own$pixel_id %in% unp, ]
  covered <- length(unique(own$pixel_id))
  ca_log("  UNP pixels ", format(length(unp), big.mark = ","),
         ", carrying an ownership record ", format(covered, big.mark = ","),
         " (", round(100 * covered / length(unp), 2), " percent)")

  comp <- own %>%
    group_by(class_value) %>%
    summarise(n_pixels = length(unique(pixel_id)), .groups = "drop") %>%
    mutate(role = ca_ownership_level$role[
             match(class_value, ca_ownership_level$label)],
           pct_of_unp = round(100 * n_pixels / length(unp), 3)) %>%
    arrange(desc(n_pixels))
  comp <- rbind(comp, data.frame(
    class_value = "NO RECORD", n_pixels = length(unp) - covered,
    role = "defaults to control",
    pct_of_unp = round(100 * (length(unp) - covered) / length(unp), 3)))
  report(comp, "unp_ownership")

  # No excluded role may survive in the control pool.
  bad <- comp$class_value[comp$role %in% c("public", "exclude_both") &
                            comp$n_pixels > 0]
  check(length(bad) == 0, paste0(
    "no excluded ownership role inside UNP",
    if (length(bad)) paste0(", found ", paste(bad, collapse = ", ")) else ""))

  # And no CPAD record of any level.
  pa <- read_arm("pa", c("pixel_id", "class_value"))
  pa_in_unp <- length(unique(pa$pixel_id[pa$pixel_id %in% unp]))
  check(pa_in_unp == 0, sprintf(
    "no CPAD record inside UNP, found %s", format(pa_in_unp, big.mark = ",")))
  rm(g, own, pa); ca_gc()
}

# ---------------------------------------------------------------------------
# 6. SPATIAL DISTRIBUTION
# ---------------------------------------------------------------------------

if (run("spatial")) {
  ca_log("Section spatial, class ", CLASS)
  g <- read_groups(c("pixel_id", "analysis_group"))
  g <- g[!is.na(g$analysis_group), ]

  eco <- ca_read_static("ecoregion_l3", CLASS)
  g$ecoregion <- eco$value[match(g$pixel_id, eco$pixel_id)]
  rm(eco); ca_gc()

  ec <- g %>%
    group_by(analysis_group, ecoregion) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(pct = round(100 * n / sum(n), 2)) %>%
    ungroup() %>%
    arrange(analysis_group, desc(n))
  report(ec, "ecoregion_by_group")

  ca_log("  pixels with no ecoregion ",
         format(sum(is.na(g$ecoregion)), big.mark = ","))

  huc <- ca_read_static("huc8", CLASS)
  g$huc8 <- huc$value[match(g$pixel_id, huc$pixel_id)]
  rm(huc); ca_gc()

  hs <- g %>%
    group_by(analysis_group) %>%
    summarise(n_huc8 = length(unique(huc8[!is.na(huc8)])),
              n_missing_huc8 = sum(is.na(huc8)), .groups = "drop")
  report(hs, "huc8_by_group")
  rm(g); ca_gc()
}

# ---------------------------------------------------------------------------
# 7. TAXONOMY TAIL
# ---------------------------------------------------------------------------

if (run("taxonomy")) {
  ca_log("Section taxonomy, class ", CLASS)
  p <- ca_meta(sprintf("groups_taxonomy_%s.csv", CLASS))
  if (!file.exists(p)) stop("Taxonomy artifact missing: ", p)
  tx <- utils::read.csv(p, stringsAsFactors = FALSE)

  cc <- tx %>%
    group_by(class_code) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n))
  tot <- sum(cc$n)
  ca_log("  distinct class codes ", format(nrow(cc), big.mark = ","))
  for (k in c(10, 100, 1000)) {
    ca_log("    codes with fewer than ", k, " pixels: ",
           format(sum(cc$n < k), big.mark = ","), " codes, ",
           round(100 * sum(cc$n[cc$n < k]) / tot, 3), " percent of the grid")
  }
  report(utils::head(cc, 200), "class_code_top")

  # Codes the analysis never uses but a follow-up study would want. Repeat
  # entry, reburn, and thinning followed by wildfire.
  seq_codes <- cc[nchar(sub("_.*$", "", cc$class_code)) > 2 &
                    !grepl("^U", cc$class_code), ]
  ca_log("  multi-token histories ", format(nrow(seq_codes), big.mark = ","),
         " codes covering ", format(sum(seq_codes$n), big.mark = ","),
         " pixels")
  report(utils::head(seq_codes, 100), "class_code_sequences")
}

# ---------------------------------------------------------------------------
# 8. MATCHING FEASIBILITY
# ---------------------------------------------------------------------------
# The section that decides whether a pairing can be estimated. Matching is 1:1
# without replacement inside exact cells, so the number of pairs a pairing can
# yield is the sum over cells of the smaller side, not the size of either side.
# A pairing with a million treated units and a thousand controls yields a
# thousand pairs at best, and stage 7 should not allocate to it as though it
# were large.

if (run("feasibility")) {
  ca_log("Section feasibility, class ", CLASS)
  g <- read_groups(c("pixel_id", "analysis_group", "stratum"))
  g <- g[!is.na(g$analysis_group), ]

  eco <- ca_read_static("ecoregion_l3", CLASS)
  g$ecoregion <- eco$value[match(g$pixel_id, eco$pixel_id)]
  rm(eco); ca_gc()
  huc <- ca_read_static("huc8", CLASS)
  g$huc8 <- huc$value[match(g$pixel_id, huc$pixel_id)]
  rm(huc); ca_gc()

  g <- g[!is.na(g$ecoregion) & !is.na(g$huc8), ]

  out <- lapply(seq_len(nrow(ca_pairings)), function(i) {
    tr <- ca_pairing_groups(ca_pairings$pairing[i], "treat")
    co <- ca_pairing_groups(ca_pairings$pairing[i], "control")
    d <- g[g$analysis_group %in% c(tr, co), ]
    if (!nrow(d)) return(NULL)
    is_tr <- d$analysis_group %in% tr
    cell <- paste(d$ecoregion, d$huc8, sep = "|")
    lev <- unique(cell)
    k <- match(cell, lev)
    nt <- tabulate(k[is_tr], nbins = length(lev))
    nc <- tabulate(k[!is_tr], nbins = length(lev))
    pairs <- pmin(nt, nc)
    data.frame(
      pairing = ca_pairings$pairing[i],
      n_treat = sum(nt),
      n_control = sum(nc),
      n_cells = length(lev),
      cells_with_both = sum(nt > 0 & nc > 0),
      max_pairs = sum(pairs),
      pct_treat_matchable = round(100 * sum(pairs) / max(sum(nt), 1), 1),
      binding_side = ifelse(sum(nc) < sum(nt), "control", "treatment"),
      stringsAsFactors = FALSE)
  })
  feas <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  report(feas, "pairing_feasibility")

  # Per stratum, because a pairing can look healthy in aggregate while one
  # severity or intensity level inside it cannot match at all.
  strat <- lapply(seq_len(nrow(ca_pairings)), function(i) {
    if (is.na(ca_pairings$strata[i])) return(NULL)
    tr <- ca_pairing_groups(ca_pairings$pairing[i], "treat")
    co <- ca_pairing_groups(ca_pairings$pairing[i], "control")
    d <- g[g$analysis_group %in% c(tr, co), ]
    if (!nrow(d)) return(NULL)
    is_tr <- d$analysis_group %in% tr
    # The control is undisturbed and carries no stratum, so each treated
    # stratum is measured against the whole control pool rather than against a
    # matching stratum.
    do.call(rbind, lapply(trimws(strsplit(ca_pairings$strata[i], ",")[[1]]),
                          function(s) {
      j <- (!is_tr) | (d$stratum %in% s)
      cell <- paste(d$ecoregion[j], d$huc8[j], sep = "|")
      lev <- unique(cell); k <- match(cell, lev)
      nt <- tabulate(k[is_tr[j]], nbins = length(lev))
      nc <- tabulate(k[!is_tr[j]], nbins = length(lev))
      data.frame(pairing = ca_pairings$pairing[i], stratum = s,
                 n_treat = sum(nt), n_control = sum(nc),
                 max_pairs = sum(pmin(nt, nc)),
                 pct_treat_matchable = round(
                   100 * sum(pmin(nt, nc)) / max(sum(nt), 1), 1),
                 stringsAsFactors = FALSE)
    }))
  })
  strat <- do.call(rbind, strat[!vapply(strat, is.null, logical(1))])
  if (!is.null(strat)) {
    report(strat[order(strat$pairing, strat$stratum), ], "stratum_feasibility")
    thin <- strat$max_pairs < 5000
    if (any(thin)) {
      ca_log("  NOTE strata yielding fewer than 5,000 pairs: ",
             paste(paste0(strat$pairing[thin], " ", strat$stratum[thin],
                          " (", strat$max_pairs[thin], ")"), collapse = "; "))
    }
  }
  rm(g); ca_gc()
}

# ---------------------------------------------------------------------------
# 9. SAMPLE
# ---------------------------------------------------------------------------
# Stage 7 verification. Four questions.
#
# Does the sample file describe the same pixels as stage 6. Is the 150 m
# guarantee actually met, measured rather than trusted, since the spacing pass
# is the one piece of stage 7 that could be wrong without failing. Does
# retention match what stage 7 reported. And what does stage 8 face once the
# spaced treated side meets the whole control pool, which is the number that
# decides whether a pairing survives.

if (run("sample")) {

  ca_log("=== 9. SAMPLE")

  sp <- ca_sample_path(CLASS)
  if (!file.exists(sp)) stop("Stage 7 output missing: ", sp)
  smp <- as.data.frame(arrow::read_parquet(sp))
  ca_log("  sample rows ", format(nrow(smp), big.mark = ","))

  # -- referential integrity against stage 6
  gg <- read_groups(c("pixel_id", "analysis_group", "arm", "stratum",
                      "cohort_year"))
  j <- match(smp$pixel_id, gg$pixel_id)

  check(anyDuplicated(smp$pixel_id) == 0L, "sample pixel_id unique")
  check(!anyNA(j), "every sampled pixel present in the stage 6 output")
  jj <- j[!is.na(j)]
  ok <- !is.na(j)
  check(identical(smp$analysis_group[ok], gg$analysis_group[jj]),
        "analysis_group agrees with stage 6")
  check(identical(smp$arm[ok], gg$arm[jj]), "arm agrees with stage 6")
  check(identical(smp$stratum[ok], gg$stratum[jj]),
        "stratum agrees with stage 6")
  check(identical(smp$cohort_year[ok], gg$cohort_year[jj]),
        "cohort_year agrees with stage 6")

  # -- role. A control-only group in this file is the error ca_pool_path exists
  # to prevent, caught here as well in case stage 7 is rerun by hand.
  treated <- trimws(unique(unlist(strsplit(ca_pairings$treat, ",",
                                           fixed = TRUE))))
  control_only <- setdiff(unique(gg$analysis_group[!is.na(gg$analysis_group)]),
                          treated)
  check(all(smp$analysis_group %in% treated),
        "sample holds treated-role groups only")
  check(!any(smp$analysis_group %in% control_only),
        paste0("control-only groups absent from the sample file (",
               paste(control_only, collapse = " "), ")"))
  rm(gg); ca_gc()

  # -- the 150 m guarantee, tested
  # One retained pixel per 150 m block, and no retained pair closer than the
  # block width. Conflicts can only arise between adjacent blocks, so the test
  # is eight offset joins per group rather than a distance matrix. Groups named
  # in ca_sample$whole_groups are exempt by construction and are reported, not
  # failed.
  ncol_grid <- ca_grid_ncol()
  blk <- as.integer(ca_sample$block_pixels)
  bc_span <- as.double(ncol_grid %/% blk + 1L)
  offs <- expand.grid(dbr = -1:1, dbc = -1:1)
  offs <- offs[!(offs$dbr == 0 & offs$dbc == 0), ]

  spacing <- do.call(rbind, lapply(sort(unique(smp$analysis_group)),
                                   function(grp) {
    pid <- smp$pixel_id[smp$analysis_group == grp]
    z <- as.double(pid) - 1
    row <- as.integer(z %/% ncol_grid)
    col <- as.integer(z %% ncol_grid)
    bkey <- as.double(row %/% blk) * bc_span + (col %/% blk)
    n_dup <- sum(duplicated(bkey))
    n_bad <- 0L
    if (length(pid) > 1L) {
      for (k in seq_len(nrow(offs))) {
        m <- match(bkey + offs$dbr[k] * bc_span + offs$dbc[k], bkey)
        a <- which(!is.na(m))
        if (!length(a)) next
        dr <- row[m[a]] - row[a]
        dc <- col[m[a]] - col[a]
        n_bad <- n_bad + sum((dr * dr + dc * dc) < blk * blk)
      }
    }
    data.frame(analysis_group = grp,
               n = length(pid),
               n_blocks = length(unique(bkey)),
               extra_per_block = n_dup,
               close_pairs = as.integer(n_bad %/% 2L),
               exempt = grp %in% ca_sample$whole_groups,
               stringsAsFactors = FALSE)
  }))
  report(spacing, "sample_spacing")

  enforced <- spacing[!spacing$exempt, ]
  check(all(enforced$extra_per_block == 0L),
        "at most one retained pixel per 150 m block")
  check(all(enforced$close_pairs == 0L),
        "no retained pair closer than 150 m")
  if (any(spacing$exempt)) {
    ca_log("  exempt from spacing by ca_sample$whole_groups: ",
           paste(spacing$analysis_group[spacing$exempt], collapse = " "))
  }

  # -- retention reconciled to what stage 7 reported
  mp <- ca_meta(sprintf("sample_summary_%s.csv", CLASS))
  if (file.exists(mp)) {
    ms <- utils::read.csv(mp, stringsAsFactors = FALSE)
    tb <- as.data.frame(table(smp$analysis_group), stringsAsFactors = FALSE)
    names(tb) <- c("analysis_group", "n_file")
    rec <- merge(ms[, c("analysis_group", "n_group", "n_spaced", "n_written",
                        "retention_pct")], tb, by = "analysis_group",
                 all = TRUE)
    rec$delta <- rec$n_file - rec$n_written
    report(rec, "sample_retention")
    check(all(!is.na(rec$delta) & rec$delta == 0L),
          "row counts match the stage 7 summary")
  } else {
    ca_log("  NOTE summary absent, retention not reconciled: ", basename(mp))
  }

  # -- pooled across classes, since lulc is an exact matching variable and
  # matching runs within class while estimation pools. Per-class counts
  # understate every group.
  pooled <- do.call(rbind, lapply(ca_lulc$label, function(cl) {
    f <- ca_meta(sprintf("sample_summary_%s.csv", cl))
    if (!file.exists(f)) return(NULL)
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d[, c("analysis_group", "n_group", "n_spaced")]
  }))
  if (!is.null(pooled)) {
    agg <- stats::aggregate(cbind(n_group, n_spaced) ~ analysis_group,
                            data = pooled, FUN = sum)
    agg$retention_pct <- round(100 * agg$n_spaced / agg$n_group, 2)
    agg$over_cap <- agg$n_spaced > ca_sample$cap_per_group
    report(agg[order(-agg$n_spaced), ], "sample_pooled")
    ca_log("  groups above the stage 8 cap of ",
           format(ca_sample$cap_per_group, big.mark = ","), ": ",
           paste(agg$analysis_group[agg$over_cap], collapse = " "))
  }

  # -- feasibility on what stage 8 actually reads
  # Section 8 measured the whole grid against itself. This measures the spaced
  # treated side against the whole control pool, which is the pairing stage 8
  # faces. ca_pool_path resolves the file per role so the two sides cannot be
  # crossed.
  eco <- ca_read_static("ecoregion_l3", CLASS)
  huc <- ca_read_static("huc8", CLASS)

  # Read every group that serves as a control in any pairing, which is UNP and
  # UP. Filtering to control_only here was wrong and reported FP_UP and TP_UP
  # as having no control pool at all. UP is a control in those two pairings and
  # a treated group in UP_UNP, which is the whole reason ca_pool_path exists,
  # and the file it resolves to for the control role holds UP whole.
  control_groups <- trimws(unique(unlist(strsplit(ca_pairings$control, ",",
                                                  fixed = TRUE))))
  cpool <- as.data.frame(arrow::read_parquet(
    ca_pool_path(CLASS, "control"),
    col_select = dplyr::all_of(c("pixel_id", "analysis_group"))))
  cpool <- cpool[!is.na(cpool$analysis_group) &
                   cpool$analysis_group %in% control_groups, ]
  ca_log("  control pool rows ", format(nrow(cpool), big.mark = ","),
         " across ", paste(control_groups, collapse = " "))

  cell_of <- function(pid) {
    paste(eco$value[match(pid, eco$pixel_id)],
          huc$value[match(pid, huc$pixel_id)], sep = "|")
  }
  smp$cell <- cell_of(smp$pixel_id)
  cpool$cell <- cell_of(cpool$pixel_id)
  rm(eco, huc); ca_gc()

  feas <- do.call(rbind, lapply(seq_len(nrow(ca_pairings)), function(i) {
    tr <- ca_pairing_groups(ca_pairings$pairing[i], "treat")
    co <- ca_pairing_groups(ca_pairings$pairing[i], "control")
    t_cell <- smp$cell[smp$analysis_group %in% tr]
    c_cell <- cpool$cell[cpool$analysis_group %in% co]
    if (!length(t_cell)) return(NULL)
    lev <- unique(c(t_cell, c_cell))
    nt <- tabulate(match(t_cell, lev), nbins = length(lev))
    nc <- tabulate(match(c_cell, lev), nbins = length(lev))
    data.frame(
      pairing = ca_pairings$pairing[i],
      n_treat_spaced = sum(nt),
      n_control_whole = sum(nc),
      n_cells = sum(nt > 0),
      cells_with_both = sum(nt > 0 & nc > 0),
      max_pairs = sum(pmin(nt, nc)),
      pct_treat_matchable = round(100 * sum(pmin(nt, nc)) / max(sum(nt), 1), 1),
      binding_side = ifelse(sum(nc) < sum(nt), "control", "treatment"),
      stringsAsFactors = FALSE)
  }))
  report(feas, "sample_feasibility")
  check(all(feas$max_pairs > 0), "every pairing yields at least one pair")
  thin <- feas$pct_treat_matchable < 90
  if (any(thin)) {
    ca_log("  NOTE pairings losing more than 10 percent of treated units to ",
           "empty control cells: ",
           paste(paste0(feas$pairing[thin], " (",
                        feas$pct_treat_matchable[thin], "%)"),
                 collapse = "; "))
  }

  rm(smp, cpool); ca_gc()
}


# ---------------------------------------------------------------------------

if (length(FAIL)) {
  ca_log("VERIFICATION FAILURES: ", length(FAIL))
  for (f in FAIL) ca_log("  ", f)
  quit(save = "no", status = 1)
}
ca_log("Section ", section, " complete, class ", CLASS, ", no failures")
ca_stamp(paste0("6_verify_groups_", section))
