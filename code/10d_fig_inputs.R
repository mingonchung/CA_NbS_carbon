## CA carbon revision pipeline
## Stage 10d. Figure inputs. Event-time support, and the estimate set the
## figures draw
##
## WHY 10d AND NOT 11g
##
## Both outputs are inputs to stage 11, so neither belongs inside it. This runs
## on Alpine after 10c, reads 10c results and the stage 9 artifacts, and writes
## two files to meta/ that the figure scripts then read. 11c cannot choose its
## x windows before this has run.
##
## THE PROBLEM
##
## ca_fig$caps sets one hi per arm and outcome, built as the outcome's last
## data year minus the first treated cohort year. That is the largest event
## time any unit in the arm can reach, and exactly one cohort reaches it. Every
## event time near the cap rests on the earliest cohorts alone, and in an
## ecoregion that can be a handful of units, which is where the tails come
## from: AGB Almanac in the Coast Range at +22, NEP under low-intensity
## thinning at +19, AGB eMapR under high-intensity thinning at +20. Those are
## not effects that appear late. They are the estimate changing what it is an
## estimate of.
##
## Two things vary that the current caps cannot see. Cohort years differ by
## ecoregion and by severity or intensity, so the same event time is supported
## by 90 percent of units in one cell and 4 percent in another. And the outcome
## windows differ, Almanac 1985 to 2025 against LEMMA 1990 to 2016, so the same
## cohort reaches a different maximum event time in each product.
##
## WHAT THIS WRITES
##
## meta/event_support.csv, one row per cell and event time:
##
##   pairing arm outcome stratum ecoregion year_first year_last n_treated
##   event_time n_units n_cohorts cohort_first cohort_last
##   n_units_0 n_units_ref e_ref n_units_den share_units share_at_0 retention
##
## n_units is the treated units whose event time is inside the outcome's year
## window, share_units is that over the cell's treated total, retention is that
## over the best-supported horizon already passed on the same side of zero, and
## n_cohorts is how many distinct cohorts contribute. A cap is the largest
## event time at which n_units, n_cohorts, and retention all clear a floor.
##
## The floors are not set here. Section 7b prints what each candidate floor
## would cut, for the two count floors and the retention rule that gate and for
## the interval-width rule that does not, so the choice is made against the
## table rather than before it, which is how screen_threshold.csv and
## nondisturbing_decision.csv were settled.
##
## WHY THIS IS CHEAP
##
## No panel is read. The cohort distribution comes from the stage 9a units
## files, which carry unit_id, treat, gvar, stratum, and ecoregion_l3 and are
## small, and the year windows come from meta/did_panel_manifest.csv, which
## stage 9b already wrote per pairing and outcome after trimming. Nothing here
## re-derives anything stage 9 decided.
##
## THREE OUTPUTS
##
##   meta/event_support.csv        section 4, units and cohorts at each event
##                                 time, per cell
##   meta/fig_caps.csv             section 7, the x window per cell, per
##                                 estimator, with the reason it stops there,
##                                 plus cap_ci, the same window under the
##                                 interval-width rule, as a sensitivity column
##   meta/did_results_fig.csv.gz   section 8, the rows the figures draw, from
##                                 the pre-trend floor to each cell's cap_use,
##                                 carrying cap_use and cap_lo_use with them
##
## Usage
##   Rscript 10d_fig_inputs.R
##
## Copy both to the laptop.

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})


## ---------------------------------------------------------------------------
## 2. THE YEAR WINDOWS
## ---------------------------------------------------------------------------
## Trimmed windows, per pairing and outcome, as stage 9b recorded them. Taking
## these from ca_did_years() instead would use the window the outcome was meant
## to have rather than the one its panel ended up with.

mpath <- ca_did_manifest_path()
if (!file.exists(mpath)) {
  stop("meta/", basename(mpath), " not found. Stage 9b writes it.")
}
PM <- utils::read.csv(mpath, stringsAsFactors = FALSE)
PM <- PM[, c("pairing", "outcome", "year_first", "year_last")]
PM <- unique(PM)
ca_log("panel manifest: ", nrow(PM), " pairing by outcome windows")


## ---------------------------------------------------------------------------
## 3. THE COHORT DISTRIBUTION
## ---------------------------------------------------------------------------
## Treated units only, one row per unit. gvar_est applies the thinning base
## shift, so the event time here is the event time the estimator reports.
##
## The always-treated block is dropped rather than floored. Under att_gt it has
## no event time at all, and under etwfe it appears in the calendar aggregation
## and not the event one, so it belongs in neither an event-time window nor the
## denominator that window is judged against.

read_treated <- function(pairing) {
  p <- ca_did_units_path(pairing)
  if (!file.exists(p)) {
    ca_log("  ", pairing, ": units file missing, skipped")
    return(NULL)
  }
  d <- as.data.frame(arrow::read_parquet(
    p, col_select = c("unit_id", "treat", "gvar", "stratum", "ecoregion_l3")))
  d <- d[d$treat == 1L, , drop = FALSE]
  d$stratum   <- as.character(d$stratum)
  d$ecoregion <- as.character(d$ecoregion_l3)
  d$stratum[is.na(d$stratum)]     <- ""
  d$ecoregion[is.na(d$ecoregion)] <- ""
  d$gvar_est  <- ca_did_gvar_est(pairing, d$gvar)
  d$pairing   <- pairing
  d[d$gvar_est > 0L, c("pairing", "unit_id", "gvar_est", "stratum",
                       "ecoregion"), drop = FALSE]
}

TR <- do.call(rbind, lapply(ca_pairings$pairing, function(p) {
  z <- read_treated(p)
  if (!is.null(z)) {
    ca_log("  ", p, ": ", format(nrow(z), big.mark = ","), " treated units, ",
           length(unique(z$gvar_est)), " cohorts, ",
           length(unique(z$ecoregion)), " ecoregions")
  }
  z
}))

## ---------------------------------------------------------------------------
## 4. SUPPORT AT EACH EVENT TIME
## ---------------------------------------------------------------------------
## A unit supports event time e for an outcome when gvar_est + e falls inside
## that outcome's trimmed window. Negative e is included, because the
## pre-trend window in ca_fig$pretrend has the same problem in the other
## direction and is currently set by hand.
##
## Cells are the figure's cells. Statewide is ecoregion "", the pooled stratum
## is stratum "", and both are computed as their own rows rather than summed
## from the parts, so a cell's denominator is its own treated total.
##
## TWO DENOMINATORS, BECAUSE THEY ANSWER DIFFERENT QUESTIONS
##
## share_units is over all treated units in the cell. It falls for two
## unrelated reasons, units aging out of the window as e grows, and outcomes
## whose window never covered most treated units at all. eMapR ends in 2017 and
## LEMMA in 2016, so most units burned after those dates have no post-period at
## any e, which made share_units report no supportable event time for those
## products.
##
## retention isolates attrition with event time. share_at_0 reports the other
## component separately, as the fraction of the cell an outcome can see at all.
##
## THE RETENTION DENOMINATOR IS THE BEST-SUPPORTED HORIZON ALREADY PASSED,
## NOT e = 0.
##
## n_units(e) counts treated units whose gvar_est + e falls inside the
## outcome's window, so as a function of e it is the cohort histogram slid
## under a fixed-length box. It rises, peaks, and falls. Where every cohort
## year sits inside the outcome window the peak is at e = 0, which is every
## Almanac cell, and the old denominator n_units(0) was correct there.
##
## It is not correct for UP_UNP against eMapR and LEMMA. Both products open in
## 1990 while CPAD establishment years run earlier, so at e = 0 only the
## post-1990 establishments are inside the window and n_units(0) is a tail
## count rather than the units at risk. Retention was then computed against
## that tail and returned values above 1, which is not a share of anything.
##
## The denominator is therefore taken as the running maximum of n_units from
## the anchor outward, e = 0 upward on the post side and e = -1 downward on the
## pre side. On the rising limb the ratio is 1, because nothing has attrited
## yet. On the falling limb it is the surviving share of the largest set the
## cell ever supported. It equals n_units(e) / n_units(0) wherever the peak is
## at the anchor, so no Almanac cell changes, and it is bounded by 1 by
## construction everywhere.
##
## n_units_ref and e_ref record the global peak and where it sits, so a cell
## whose support peaks away from zero is visible in the file rather than only
## in its retention column.
##
## Section 7 gates on n_units, n_cohorts, and retention. share_at_0 and
## share_units are reported and gate nothing, because a low share_at_0 is a
## property of the outcome window rather than of drift along the curve.

E_LO <- -15L
E_HI <- 45L

support_one <- function(g, year_first, year_last) {
  e <- E_LO:E_HI
  vapply(e, function(ee) {
    ok <- (g + ee) >= year_first & (g + ee) <= year_last
    c(n_units = sum(ok), n_cohorts = length(unique(g[ok])),
      cohort_first = if (any(ok)) min(g[ok]) else NA_real_,
      cohort_last  = if (any(ok)) max(g[ok]) else NA_real_)
  }, numeric(4))
}

## The running maximum of n_units from the anchor outward, e = 0 upward on the
## post side and e = -1 downward on the pre side. The two sides are disjoint,
## so one column carries both and every row is gated against the anchor its own
## side uses. See the denominator note in the block above.

support_denominator <- function(n, e) {
  den <- numeric(length(e))
  i_post <- which(e >= 0)                 # ascending from 0
  i_pre  <- rev(which(e <  0))            # descending from -1
  if (length(i_post)) den[i_post] <- cummax(n[i_post])
  if (length(i_pre))  den[i_pre]  <- cummax(n[i_pre])
  den
}

cell_rows <- function(d, pairing, outcome, stratum, ecoregion,
                      year_first, year_last) {
  if (!nrow(d)) return(NULL)
  s <- support_one(d$gvar_est, year_first, year_last)

  e   <- E_LO:E_HI
  n   <- s["n_units", ]
  n0  <- n[match(0L, e)]
  den <- support_denominator(n, e)

  ## The global peak, reported so a cell whose support does not peak at the
  ## anchor is identifiable without recomputing anything.
  n_ref <- max(n)
  e_ref <- e[which.max(n)]

  data.frame(
    pairing      = pairing,
    arm          = ca_pairings$arm[match(pairing, ca_pairings$pairing)],
    outcome      = outcome,
    stratum      = stratum,
    ecoregion    = ecoregion,
    year_first   = year_first,
    year_last    = year_last,
    n_treated    = nrow(d),
    event_time   = e,
    n_units      = as.integer(n),
    n_cohorts    = as.integer(s["n_cohorts", ]),
    cohort_first = as.integer(s["cohort_first", ]),
    cohort_last  = as.integer(s["cohort_last", ]),
    n_units_0    = as.integer(n0),
    n_units_ref  = as.integer(n_ref),
    e_ref        = as.integer(e_ref),
    n_units_den  = as.integer(den),
    share_units  = round(n / nrow(d), 4),
    share_at_0   = round(n0 / nrow(d), 4),
    retention    = ifelse(den > 0, round(n / den, 4), NA_real_),
    stringsAsFactors = FALSE)
}

OUT <- list()
k <- 0L

for (i in seq_len(nrow(PM))) {

  pg <- PM$pairing[i]; oc <- PM$outcome[i]
  yf <- PM$year_first[i]; yl <- PM$year_last[i]

  d <- TR[TR$pairing == pg, , drop = FALSE]
  if (!nrow(d)) next

  strata <- unique(d$stratum)
  strata <- if (length(strata) == 1L && strata == "") "" else c("", strata)
  strata <- unique(strata)

  for (st in strata) {
    ds <- if (nzchar(st)) d[d$stratum == st, , drop = FALSE] else d
    if (!nrow(ds)) next
    for (ec in c("", sort(unique(ds$ecoregion)))) {
      de <- if (nzchar(ec)) ds[ds$ecoregion == ec, , drop = FALSE] else ds
      k <- k + 1L
      OUT[[k]] <- cell_rows(de, pg, oc, st, ec, yf, yl)
    }
  }
}

SUP <- do.call(rbind, OUT)
SUP <- SUP[order(match(SUP$pairing, ca_pairings$pairing),
                 match(SUP$outcome, ca_outcomes$outcome),
                 SUP$stratum, SUP$ecoregion, SUP$event_time), ]
rownames(SUP) <- NULL

utils::write.csv(SUP, ca_meta("event_support.csv"), row.names = FALSE)
ca_log("Wrote event_support.csv  ", format(nrow(SUP), big.mark = ","),
       " rows, ",
       length(unique(paste(SUP$pairing, SUP$outcome, SUP$stratum,
                           SUP$ecoregion))), " cells")


## ---------------------------------------------------------------------------
## 5. HOW THE COHORTS ARE DISTRIBUTED
## ---------------------------------------------------------------------------
## Context for section 7, not a decision. The horizon at which half the units
## have aged out is a property of when treatment happened, and it differs by an
## order of magnitude between the two stratified arms, so no single event-time
## cap can serve both.

half_life <- function(z) {
  z <- z[z$event_time >= 0, , drop = FALSE]
  z <- z[order(z$event_time), ]
  i <- which(z$retention < 0.5)
  if (!length(i)) return(max(z$event_time))
  z$event_time[i[1]]
}

sw <- SUP[SUP$ecoregion == "" & SUP$outcome == "agb_almanac", , drop = FALSE]
hl <- do.call(rbind, lapply(
  split(sw, paste(sw$pairing, sw$stratum)), function(z) {
    data.frame(cell = paste(z$pairing[1], z$stratum[1]),
               n_treated = z$n_treated[1],
               share_at_0 = z$share_at_0[1],
               half_gone_at = half_life(z),
               stringsAsFactors = FALSE)
  }))
cat("\nHALF OF TREATED UNITS AGED OUT BY EVENT TIME, Almanac window\n")
print(hl, row.names = FALSE)

## WHERE THE SUPPORT DOES NOT PEAK AT e = 0.
##
## The audit for the retention denominator. A cell listed here has cohort years
## outside its outcome window, so n_units(0) is a tail count and the old
## denominator would have returned retention above 1. Expected for UP_UNP
## against eMapR and LEMMA, whose windows open in 1990, and expected to be
## empty for every Almanac cell. Anything else here is worth reading before the
## caps are taken.
cat("\nCELLS WHOSE SUPPORT PEAKS AWAY FROM e = 0\n")
pk <- SUP[SUP$event_time == 0 & SUP$e_ref != 0L,
          c("pairing", "outcome", "stratum", "ecoregion", "n_treated",
            "n_units_0", "n_units_ref", "e_ref")]
if (nrow(pk)) print(pk, row.names = FALSE) else cat("  none\n")


## ---------------------------------------------------------------------------
## 6. THE FIGURE ROW SET, BEFORE ANY TIME WINDOW
## ---------------------------------------------------------------------------
## did_results.csv.gz holds every aggregation, every transform, and every
## ecoregion, of which the figures draw a subset. This selects that subset on
## everything except time. The time window comes from section 7, which needs
## this table to compute it, so the two cannot be done in one pass.
##
## BOTH TRANSFORMS ARE KEPT, THE FIGURES STILL DRAW ONE
##
## ca_fig$panel_rows fixes one transform per panel row, none for the flux row
## and log1p for the stock row. Selecting on that pair would write a file
## holding GPP in gC m-2 yr-1 and no percent, and the three AGB products in
## percent and no gC m-2. The Results section quotes both scales in the same
## sentence, "1.9 +/- 0.5% in GPP, corresponding to 96.7 +/- 21.3 gC/m2/yr" at
## L449 to L450 of the submitted text, so a file carrying one scale forces the
## other to be recomputed by hand, which is how a quoted number stops matching
## the figure it came from.
##
## Rows are therefore selected on outcome, and every transform stage 10
## produced for that outcome is kept. fig_transform marks the one the figures
## draw. ca_fig_series() merges on outcome and transform against panel_rows, so
## the extra rows are invisible to 11c, 11d, and 11e and no figure script
## changes.
##
## NPP, NEP, and NBP are signed and were never estimated on a log scale, so
## they carry one transform for a data reason rather than a filtering one. The
## count per outcome and transform is printed below so that absence is visible
## rather than assumed.
##
## Ecoregions are the six forested codes plus statewide; stage 10 computes up
## to twelve and no figure draws them.

RES_IN  <- ca_did_results_path("results")
if (!file.exists(RES_IN)) stop("Not found: ", RES_IN, ". Run 10c aggregate.")

R <- utils::read.csv(RES_IN, colClasses = c(stratum = "character",
                                            ecoregion = "character"),
                     stringsAsFactors = FALSE)
R$stratum[is.na(R$stratum)]     <- ""
R$ecoregion[is.na(R$ecoregion)] <- ""
ca_log("read ", format(nrow(R), big.mark = ","), " result rows")

WANT <- do.call(rbind, lapply(names(ca_fig$panel_rows), function(rn) {
  data.frame(outcome       = ca_fig$panel_rows[[rn]]$outcomes,
             block         = rn,
             fig_transform = ca_fig$panel_rows[[rn]]$transform,
             stringsAsFactors = FALSE)
}))

FIG_SPEC <- rbind(
  data.frame(figure = "Figure 2",  estimator = "etwfe", arm = "pa",
             pairing = "UP_UNP",  aggregation = "calendar",
             ecoregions = "all",       strata = "pooled"),
  data.frame(figure = "Figure 3",  estimator = "etwfe", arm = "fire",
             pairing = "FP_UP",   aggregation = "event",
             ecoregions = "statewide", strata = "stratified"),
  data.frame(figure = "Figure 3",  estimator = "etwfe", arm = "fire",
             pairing = "FNP_UNP", aggregation = "event",
             ecoregions = "statewide", strata = "stratified"),
  data.frame(figure = "Figure 4",  estimator = "etwfe", arm = "thin",
             pairing = "TP_UP",   aggregation = "event",
             ecoregions = "statewide", strata = "stratified"),
  data.frame(figure = "Figure 4",  estimator = "etwfe", arm = "thin",
             pairing = "TNP_UNP", aggregation = "event",
             ecoregions = "statewide", strata = "stratified"),
  data.frame(figure = "Figure 5",  estimator = "etwfe", arm = "offset",
             pairing = "ONP_UNP", aggregation = "event",
             ecoregions = "all",       strata = "pooled"),
  data.frame(figure = "Figure S2", estimator = "etwfe", arm = "pa",
             pairing = "UP_UNP",  aggregation = "calendar",
             ecoregions = "eco",       strata = "pooled"),
  data.frame(figure = "Figure S3", estimator = "etwfe", arm = "fire",
             pairing = "FP_UP",   aggregation = "event",
             ecoregions = "eco",       strata = "stratified"),
  data.frame(figure = "Figure S3", estimator = "etwfe", arm = "fire",
             pairing = "FNP_UNP", aggregation = "event",
             ecoregions = "eco",       strata = "stratified"),
  data.frame(figure = "Figure S4", estimator = "etwfe", arm = "thin",
             pairing = "TP_UP",   aggregation = "event",
             ecoregions = "eco",       strata = "stratified"),
  data.frame(figure = "Figure S4", estimator = "etwfe", arm = "thin",
             pairing = "TNP_UNP", aggregation = "event",
             ecoregions = "eco",       strata = "stratified"),
  data.frame(figure = "Figure S5", estimator = "etwfe", arm = "offset",
             pairing = "ONP_UNP", aggregation = "event",
             ecoregions = "eco",       strata = "pooled"),
  stringsAsFactors = FALSE)

ATT <- FIG_SPEC[FIG_SPEC$figure %in% paste("Figure", 2:5), ]
ATT$estimator   <- "attgt"
ATT$figure      <- paste0(ATT$figure, " attgt")
ATT$aggregation <- ifelse(ATT$aggregation == "event", "dynamic",
                          ATT$aggregation)
FIG_SPEC <- rbind(FIG_SPEC, ATT)

PRE <- data.frame(
  figure      = "Figure S6",
  estimator   = rep(c("etwfe", "attgt"), each = nrow(ca_pairings)),
  arm         = rep(ca_pairings$arm, 2),
  pairing     = rep(ca_pairings$pairing, 2),
  aggregation = rep(c("event", "dynamic"), each = nrow(ca_pairings)),
  ecoregions  = "statewide",
  strata      = ifelse(rep(ca_pairings$arm, 2) %in% c("fire", "thin"),
                       "stratified", "pooled"),
  stringsAsFactors = FALSE)
FIG_SPEC <- rbind(FIG_SPEC, PRE)

pick <- function(sp) {
  z <- R[R$estimator == sp$estimator &
           R$pairing == sp$pairing &
           R$aggregation == sp$aggregation, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  z <- merge(z, WANT, by = "outcome")
  if (!nrow(z)) return(NULL)
  z <- switch(sp$ecoregions,
    statewide = z[z$ecoregion == "", , drop = FALSE],
    eco       = z[z$ecoregion %in% names(ca_fig$ecoregions), , drop = FALSE],
    all       = z[z$ecoregion %in% c("", names(ca_fig$ecoregions)), ,
                  drop = FALSE])
  if (!nrow(z)) return(NULL)
  z <- if (identical(sp$strata, "stratified")) {
    z[z$stratum != "", , drop = FALSE]
  } else z[z$stratum == "", , drop = FALSE]
  if (!nrow(z)) return(NULL)
  z$figure <- sp$figure
  z
}

FIG <- do.call(rbind, lapply(seq_len(nrow(FIG_SPEC)), function(i) {
  z <- pick(FIG_SPEC[i, ])
  if (is.null(z)) {
    ca_log("  EMPTY  ", FIG_SPEC$figure[i], "  ", FIG_SPEC$estimator[i], " ",
           FIG_SPEC$pairing[i], " ", FIG_SPEC$aggregation[i])
  }
  z
}))
if (is.null(FIG) || !nrow(FIG)) stop("No rows selected. Check FIG_SPEC.")

## One row per estimate, with the figures listed, so the same number cannot be
## quoted twice with two different provenances.
key_cols <- c("estimator", "pairing", "outcome", "transform", "stratum",
              "ecoregion", "aggregation", "time")
key  <- do.call(paste, c(FIG[, key_cols], sep = "|"))
figs <- tapply(FIG$figure, key, function(x) paste(sort(unique(x)),
                                                  collapse = "; "))
FIG <- FIG[!duplicated(key), , drop = FALSE]
FIG$figures <- unname(figs[do.call(paste, c(FIG[, key_cols], sep = "|"))])
FIG$figure  <- NULL

## THE DISPLAY SCALE, PER ROW RATHER THAN PER PANEL ROW.
##
## The same selection ca_fig_series() makes, applied to whichever transform the
## row carries rather than to the one its panel row draws. A log1p row is a
## percent change through expm1(), a non-log row is in the outcome's native
## units. Both are written so a Results sentence can quote the two scales for
## one cell and one horizon without recomputing either.
##
## estimate_bt is NA on every non-log row by construction at 10c, so a log1p
## row that lost its back transformation reads NA rather than reading as a
## native value in the wrong units.

is_log <- FIG$transform == "log1p"

FIG$value <- ifelse(is_log, FIG$estimate_bt  * 100, FIG$estimate)
FIG$lo    <- ifelse(is_log, FIG$conf.low_bt  * 100, FIG$conf.low)
FIG$hi    <- ifelse(is_log, FIG$conf.high_bt * 100, FIG$conf.high)

FIG$unit <- ifelse(is_log, "percent",
                   ifelse(FIG$block == "stock", "gC m-2", "gC m-2 yr-1"))

## Which rows the figures draw. A row whose interval is missing is not
## drawable, so it is not a figure row, regardless of transform.
FIG$fig_row <- FIG$transform == FIG$fig_transform &
  !is.na(FIG$value) & !is.na(FIG$lo) & !is.na(FIG$hi)

ne <- FIG$transform == FIG$fig_transform & !FIG$fig_row
if (any(ne)) {
  cat("\nNON-ESTIMABLE ROWS, demoted from fig_row\n")
  print(FIG[ne, c("pairing", "outcome", "stratum", "ecoregion", "transform",
                  "time", "estimate", "n_cohorts")], row.names = FALSE)
}

## THE INTERVAL THE FIGURE ACTUALLY DRAWS.
##
## Taken on the display scale of the row, so a width is the width a reader
## would see. Section 7 uses it on the figure rows only.
FIG$ci_width <- FIG$hi - FIG$lo

ca_log("figure row set before windows: ", format(nrow(FIG), big.mark = ","),
       ", of which ", format(sum(FIG$fig_row), big.mark = ","),
       " are drawn")

## Both scales present, per outcome. A blank column is an outcome stage 10
## never estimated on that scale, which is the three signed fluxes and is not
## a filtering loss.
cat("\nROWS BY OUTCOME AND TRANSFORM\n")
print(table(FIG$outcome, FIG$transform))


## ---------------------------------------------------------------------------
## 7. THE X WINDOW, PER CELL AND PER ESTIMATOR
## ---------------------------------------------------------------------------
## An event-study estimate at horizon e averages over the cohorts still
## observable at e. As e grows the cohort set shrinks toward the earliest
## treatments, so the curve mixes recovery with cohort heterogeneity. A visible
## break appears when few cohorts remain, because one cohort leaving moves the
## average by roughly 1/k.
##
## THREE GATING CONDITIONS, ALL ON DESIGN QUANTITIES
##
##   n_units   >= N_MIN     the stage 10 design floor, applied again at e
##   n_cohorts >= COH_MIN   an average over two cohorts is two events
##   retention >= RET_MIN   the surviving share against the best-supported
##                          horizon already passed, section 4
##
## Attrition here is deterministic. A unit supports e when gvar_est + e is
## inside the outcome window, so the retained set at e is exactly the cohorts
## with gvar_est <= year_last - e. Composition therefore drifts along one axis,
## cohort year, and retention is the empirical CDF at that cutoff. All three
## conditions are functions of cohort dates, panel windows, and unit counts,
## fixed before any estimate exists, so the window is pre-specifiable.
##
## The cap is the leading run of satisfied conditions from the anchor, so a
## horizon that rebounds after a cut cannot restore a window. Retention falls
## monotonically once attrition begins and is 1 before it, which is why no
## single noisy horizon can set the cap. The two count floors do the work
## retention cannot: FP_UP High severity and FNP_UNP Low severity hold 133 and
## 2,352 units at e = 33, which retention treats alike and N_MIN separates.
##
## THE INTERVAL RULE IS REPORTED, NOT GATED
##
## cap_ci applies ci_width <= CI_M * ci_ref, where ci_ref is the median width
## over e in 0 to CI_REF_HI. It is carried as a sensitivity column so the
## windows can be compared against a precision-based alternative, and it gates
## nothing, for two reasons. Interval width is already drawn in the ribbon, so
## the reader can judge it, whereas composition drift is invisible and must be
## ruled on. And for the log1p outcomes ci_width is taken on the back
## transformed percent scale, where it scales with exp(beta) as well as with
## the standard error, so it moves with the effect level during recovery.
##
## e = 0 is exempt from the interval condition because ci_ref is anchored on a
## window that contains it, so a wide burn year would otherwise void the cell.
##
## PER ESTIMATOR, BECAUSE THE BANDS ARE NOT THE SAME OBJECT
##
## attgt_cband is TRUE, so att_gt returns uniform bands and etwfe returns
## pointwise intervals. cap is identical across estimators, since it depends on
## no estimate, but cap_ci is not. cap_use is the smaller of the two, so a
## paired panel shares one window and a difference between the two figures is a
## difference in estimates.

N_MIN      <- ca_did_est$min_treated_units
COH_MIN    <- 3L
RET_MIN    <- 0.05

## Sensitivity only. CI_M sets cap_ci, which is written to fig_caps.csv and
## printed against cap in section 7b. It does not enter cap_use.
CI_REF_HI  <- 5L
CI_M       <- 2.0

## Calendar arms are not capped here. Protected areas is drawn on calendar time
## in both estimators, where the always-treated block contributes at every year
## and there is no cohort attrition to measure.
CAL_AGG <- c("calendar")

sup_key  <- c("pairing", "outcome", "stratum", "ecoregion")
cell_key <- c(sup_key, "estimator", "aggregation")

## The leading run from e = 0. A cell cannot regain a window after losing it,
## because the composition never recovers.
run_cap <- function(e, ok) {
  if (!length(ok) || !isTRUE(ok[1])) return(NA_integer_)
  r <- rle(ok)
  as.integer(e[r$lengths[1]])
}

## The same going backwards from e = -1, for Figure S6.
run_cap_lo <- function(e, ok) {
  o <- order(e, decreasing = TRUE)
  e <- e[o]; ok <- ok[o]
  if (!length(ok) || !isTRUE(ok[1])) return(NA_integer_)
  r <- rle(ok)
  as.integer(e[r$lengths[1]])
}

## ret, m, and k default to the globals, so every call outside the sweeps is
## unchanged. They exist so section 7b can sweep each floor without a second
## copy of this function.

cap_one <- function(z, ret = RET_MIN, m = CI_M, k = COH_MIN) {

  s <- SUP[SUP$pairing == z$pairing[1] & SUP$outcome == z$outcome[1] &
             SUP$stratum == z$stratum[1] & SUP$ecoregion == z$ecoregion[1], ,
           drop = FALSE]
  if (!nrow(s)) return(NULL)

  d <- merge(z[, c("time", "ci_width")], s, by.x = "time",
             by.y = "event_time")
  if (!nrow(d)) return(NULL)
  d <- d[order(d$time), ]

  post <- d[d$time >= 0, , drop = FALSE]

  ## The gating rule. Design quantities only.
  ok_n <- post$n_units >= N_MIN
  ok_k <- post$n_cohorts >= k
  ok_r <- !is.na(post$retention) & post$retention >= ret

  cap <- run_cap(post$time, ok_n & ok_k & ok_r)

  ## Which condition stopped it, taken at the first failing event time.
  bind <- NA_character_
  if (!is.na(cap)) {
    j <- which(post$time == cap) + 1L
    if (j <= nrow(post)) {
      bind <- paste(c("n_units", "n_cohorts", "retention")[
        !c(ok_n[j], ok_k[j], ok_r[j])], collapse = "+")
    } else bind <- "end_of_data"
  }

  ## The sensitivity rule. Same two count floors, interval width in place of
  ## retention. e = 0 is exempt, since ci_ref is anchored on a window that
  ## contains it.
  ref <- d$ci_width[d$time >= 0 & d$time <= CI_REF_HI]
  ci_ref <- if (length(ref)) stats::median(ref, na.rm = TRUE) else NA_real_

  ok_c <- if (is.na(ci_ref) || ci_ref <= 0) rep(TRUE, nrow(post)) else
    post$ci_width <= m * ci_ref
  ok_c[post$time == 0] <- TRUE
  ok_c[is.na(ok_c)] <- FALSE

  cap_ci <- run_cap(post$time, ok_n & ok_k & ok_c)

  bind_ci <- NA_character_
  if (!is.na(cap_ci)) {
    j <- which(post$time == cap_ci) + 1L
    if (j <= nrow(post)) {
      bind_ci <- paste(c("n_units", "n_cohorts", "ci_width")[
        !c(ok_n[j], ok_k[j], ok_c[j])], collapse = "+")
    } else bind_ci <- "end_of_data"
  }

  ## Pre-treatment side, for Figure S6. Same gating rule, run backwards.
  pre <- d[d$time < 0, , drop = FALSE]
  cap_lo <- if (nrow(pre)) {
    ok <- pre$n_units >= N_MIN & pre$n_cohorts >= k &
      !is.na(pre$retention) & pre$retention >= ret
    run_cap_lo(pre$time, ok)
  } else NA_integer_

  at <- if (!is.na(cap)) post[post$time == cap, ] else post[0, ]

  data.frame(
    pairing      = z$pairing[1],
    arm          = s$arm[1],
    outcome      = z$outcome[1],
    stratum      = z$stratum[1],
    ecoregion    = z$ecoregion[1],
    estimator    = z$estimator[1],
    aggregation  = z$aggregation[1],
    year_first   = s$year_first[1],
    year_last    = s$year_last[1],
    n_treated    = s$n_treated[1],
    share_at_0   = s$share_at_0[1],
    cap          = cap,
    cap_lo       = cap_lo,
    bind         = bind,
    cap_ci       = cap_ci,
    bind_ci      = bind_ci,
    n_units      = if (nrow(at)) at$n_units[1] else NA_integer_,
    n_cohorts    = if (nrow(at)) at$n_cohorts[1] else NA_integer_,
    retention    = if (nrow(at)) at$retention[1] else NA_real_,
    cohort_first = if (nrow(at)) at$cohort_first[1] else NA_integer_,
    cohort_last  = if (nrow(at)) at$cohort_last[1] else NA_integer_,
    ci_ref       = round(ci_ref, 4),
    ci_at_cap    = if (nrow(at)) round(at$ci_width[1], 4) else NA_real_,
    stringsAsFactors = FALSE)
}

## CAPS ARE COMPUTED ON THE FIGURE TRANSFORM ONLY.
##
## cap depends on no estimate, so it is identical on both scales and computing
## it twice would only duplicate rows. cap_ci does depend on the estimate, and
## the window it describes is the window a reader sees, so it is taken on the
## scale the figure draws. The reduction below then spreads cap_use and
## cap_lo_use over every row sharing pairing, outcome, stratum, and ecoregion,
## which is both transforms, so the two scales are windowed identically and a
## quoted number and a plotted point come from the same horizons.
##
## This also keeps time unique inside each group. Splitting on cell_key with
## both transforms present would give two rows per event time, and the merge
## against SUP in cap_one() would duplicate them.

EV <- FIG[!(FIG$aggregation %in% CAL_AGG) & FIG$fig_row, , drop = FALSE]
sp <- split(EV, do.call(paste, c(EV[, cell_key], sep = "|")))

CAPS <- do.call(rbind, lapply(sp, cap_one))
rownames(CAPS) <- NULL

## cap_use, the smaller of the two estimators for the same cell, so paired
## panels share a window. A cell only one estimator has keeps its own.
ck <- do.call(paste, c(CAPS[, sup_key], sep = "|"))
CAPS$cap_use <- ave(CAPS$cap, ck, FUN = function(x)
  if (all(is.na(x))) NA_integer_ else min(x, na.rm = TRUE))
CAPS$cap_lo_use <- ave(CAPS$cap_lo, ck, FUN = function(x)
  if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE))

## The same reduction on the sensitivity cap, so the two are comparable at the
## level the figures are drawn. Reported only.
CAPS$cap_ci_use <- ave(CAPS$cap_ci, ck, FUN = function(x)
  if (all(is.na(x))) NA_integer_ else min(x, na.rm = TRUE))

CAPS <- CAPS[order(match(CAPS$pairing, ca_pairings$pairing),
                   match(CAPS$outcome, ca_outcomes$outcome),
                   CAPS$stratum, CAPS$ecoregion, CAPS$estimator), ]
rownames(CAPS) <- NULL

utils::write.csv(CAPS, ca_meta("fig_caps.csv"), row.names = FALSE)
ca_log("Wrote fig_caps.csv  ", nrow(CAPS), " cell by estimator rows")


## ---------------------------------------------------------------------------
## 7b. CHOOSING THE FLOOR AGAINST THE DATA
## ---------------------------------------------------------------------------
## RET_MIN is the only gating threshold here without an external justification,
## so it is printed rather than argued. Read the statewide table, set RET_MIN
## at the top, rerun. The CI_M sweep is printed beside it so the cost of the
## rule that gates can be read against the rule that does not.
##
## The RET_MIN = 0 row is the counts-only rule, which is what n_units and
## n_cohorts alone would give. The CI_M = Inf row is the same thing, so the two
## sweeps meet there and any difference between the columns at that setting is
## a coding error.

sw_ev <- EV[EV$ecoregion == "", , drop = FALSE]
sw_sp <- split(sw_ev, do.call(paste, c(sw_ev[, cell_key], sep = "|")))

sweep_print <- function(vals, arg, col) {
  for (v in vals) {
    a <- switch(arg,
      ret = do.call(rbind, lapply(sw_sp, cap_one, ret = v)),
      coh = do.call(rbind, lapply(sw_sp, cap_one, k   = v)),
            do.call(rbind, lapply(sw_sp, cap_one, m   = v)))
    a$val <- a[[col]]
    tb <- stats::aggregate(val ~ arm + estimator, data = a,
                           FUN = function(x) c(min = min(x),
                                               med = stats::median(x),
                                               max = max(x)))
    cat("\n  ",
        switch(arg, ret = "RET_MIN = ", coh = "COH_MIN = ", "CI_M = "), v,
        "  (", sum(!is.na(a[[col]])), " of ", nrow(a), " cells capped)\n",
        sep = "")
    print(do.call(data.frame, tb), row.names = FALSE)
  }
}

cat("\nSTATEWIDE CAPS BY RET_MIN, the gating rule\n")
sweep_print(c(0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50), "ret", "cap")

## COH_MIN is the floor that binds most often, so it is swept on the same
## terms as RET_MIN rather than asserted. Two is the arithmetic minimum for an
## average over cohorts to mean anything. Read the fire and offset rows: those
## are the arms where the floor binds, because their cohorts are concentrated
## late and few survive to long horizons.
cat("\nSTATEWIDE CAPS BY COH_MIN, the second count floor\n")
sweep_print(c(2L, 3L, 4L, 5L, 6L), "coh", "cap")

## How many windows each candidate floor binds, across every cell rather than
## the statewide ones alone. This is the number the SI sentence quotes.
cat("\nWINDOWS BOUND BY COH_MIN, all cells\n")
coh_bind <- do.call(rbind, lapply(c(2L, 3L, 4L, 5L, 6L), function(v) {
  a <- do.call(rbind, lapply(
    split(EV, do.call(paste, c(EV[, cell_key], sep = "|"))), cap_one, k = v))
  data.frame(COH_MIN = v,
             n_cells   = nrow(a),
             n_capped  = sum(!is.na(a$cap)),
             n_bind_coh = sum(grepl("n_cohorts", a$bind)),
             med_cap   = stats::median(a$cap, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
print(coh_bind, row.names = FALSE)

cat("\nSTATEWIDE CAPS BY CI_M, sensitivity only\n")
sweep_print(c(1.5, 2.0, 3.0, 5.0, Inf), "m", "cap_ci")

cat("\nWHICH CONDITION BINDS, at RET_MIN =", RET_MIN, "\n")
print(table(CAPS$arm, CAPS$bind))

cat("\nSTATEWIDE CAPS, etwfe against attgt, at RET_MIN =", RET_MIN, "\n")
w <- reshape(CAPS[CAPS$ecoregion == "",
                  c("pairing", "outcome", "stratum", "estimator", "cap")],
             idvar = c("pairing", "outcome", "stratum"),
             timevar = "estimator", direction = "wide")
print(w, row.names = FALSE)

## cap does not depend on any estimate, so the two estimator columns above
## should agree wherever both are present. A disagreement means one estimator
## is missing event times the other has, which is a 10c coverage problem and
## not a windowing one.

cat("\nGATING CAP AGAINST SENSITIVITY CAP, statewide, RET_MIN =", RET_MIN,
    " CI_M =", CI_M, "\n")
v <- CAPS[CAPS$ecoregion == "",
          c("arm", "pairing", "outcome", "stratum", "estimator",
            "cap", "cap_ci", "bind", "bind_ci", "n_units", "retention")]
v$delta <- v$cap_ci - v$cap
print(v[, setdiff(names(v), "arm")], row.names = FALSE)

cat("\nDELTA SUMMARY, cap_ci minus cap, by arm and estimator\n")
ds <- stats::aggregate(delta ~ arm + estimator, data = v,
                       FUN = function(x) c(min = min(x),
                                           med = stats::median(x),
                                           max = max(x)))
print(do.call(data.frame, ds), row.names = FALSE)

## A positive delta is a horizon the interval rule would have kept and the
## retention rule cuts. A negative delta is the reverse.

cat("\nCELLS WITH NO SUPPORTABLE WINDOW\n")
nw <- CAPS[is.na(CAPS$cap),
           c("pairing", "outcome", "stratum", "ecoregion", "estimator",
             "n_treated", "share_at_0", "cap_lo", "cap_ci")]
if (nrow(nw)) print(nw, row.names = FALSE) else cat("  none\n")


## ---------------------------------------------------------------------------
## 8. APPLY THE WINDOWS AND WRITE
## ---------------------------------------------------------------------------
## THE SUBSET IS THE FIGURE SET, ON TWO SCALES
##
## did_results_fig.csv.gz holds the rows the stage 11 scripts draw plus their
## companions on the other transform, over the same horizons. A number in a
## panel and a number in the file are the same number, no figure script
## re-derives a window, and a Results sentence quoting percent and gC for one
## cell reads two rows of one file rather than two files. fig_row marks the
## drawn rows, and ca_fig_series() selects them on outcome and transform
## without reading that column. One rule per aggregation.
##
##   event      from the pre-trend floor to that cell's cap_use
##   calendar   the ca_fig$caps$pa window, a data coverage bound
##
## The event span runs from ca_fig$pretrend["lo"], currently -9, rather than
## from 0. Figures 2 to 5 and S2 to S5 read the part from 0 up, Figure S6 reads
## the part below 0, and the rows at -1 and 0 sit in both. Earlier this file
## kept negative event time only for rows tagged Figure S6 and only out to -2,
## which meant a second read of the full results to draw anything else.
##
## PRE_GATED holds the pre side to what the cell supports. A cell whose earliest
## supported pre horizon is -5 starts at -5, not at -9, so 11e never draws a
## point the three floors reject. Set it FALSE for a flat -9 on every cell.
##
## UP_UNP appears twice and should. Its calendar rows are Figures 2 and S2, its
## event rows are Figure S6, and the two are different aggregations of the same
## cell rather than duplicates.
##
## WHICH WINDOW FILE
##
## Section 7 always writes fig_caps.csv, so the automatic windows stay on
## record. Where meta/fig_caps_rev.csv exists it is the hand-checked copy, cut
## at horizons the three floors allowed and the panels could not carry, and it
## is what this section applies. The difference between the two files is the
## record of what was cut.
##
## The file is resolved through ca_fig_caps_path(), the same call stage 11
## makes, so the rows written here and the lines drawn later cannot come from
## different windows. The full path is logged rather than the basename, because
## ca_fig_find() searches CA_DATA before meta and a stale copy there would
## otherwise be indistinguishable in the log.
##
## A cell the revised file omits has no cap_use and is dropped by the keep rule
## below rather than falling back to its automatic window. A stale revised file
## therefore shows up as missing rows and never as a silently wider panel,
## which is why the unlisted count is printed and should be zero.

CAPS_USE <- ca_fig_caps_read()
ca_log("windows read from ", attr(CAPS_USE, "path"))

need <- c(sup_key, "cap_use", "cap_lo_use")
miss <- setdiff(need, names(CAPS_USE))
if (length(miss)) {
  stop(basename(attr(CAPS_USE, "path")), " is missing: ",
       paste(miss, collapse = ", "))
}

## The same reduction ca_fig_caps() applies. cap_use is replicated across the
## two estimator rows of a cell by construction, so taking the smaller keeps
## the window independent of row order and a cell edited on one estimator row
## only still cuts both panels of the pair.
rk <- do.call(paste, c(CAPS_USE[, sup_key], sep = "|"))
CAPS_USE$cap_use <- ave(CAPS_USE$cap_use, rk, FUN = function(x)
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
CAPS_USE$cap_lo_use <- ave(CAPS_USE$cap_lo_use, rk, FUN = function(x)
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

## What the applied windows differ from, against the rule section 7 just ran.
ak <- do.call(paste, c(CAPS[, sup_key], sep = "|"))
cmp <- merge(CAPS[!duplicated(ak), c(sup_key, "cap_use", "cap_lo_use", "bind")],
             CAPS_USE[!duplicated(rk), c(sup_key, "cap_use", "cap_lo_use")],
             by = sup_key, all.x = TRUE, suffixes = c("_auto", "_use"))

changed <- !is.na(cmp$cap_use_use) & !is.na(cmp$cap_use_auto) &
  cmp$cap_use_use != cmp$cap_use_auto
unlisted <- is.na(cmp$cap_use_use) & !is.na(cmp$cap_use_auto)

cat("\nAPPLIED WINDOWS AGAINST THE SECTION 7 RULE\n")
cat("  ", nrow(cmp), " cells, ", sum(changed), " with a different cap_use, ",
    sum(unlisted), " not listed in the applied file\n", sep = "")
if (any(changed)) {
  print(cmp[changed, c(sup_key, "cap_use_auto", "cap_use_use", "bind")],
        row.names = FALSE)
}
if (any(unlisted)) {
  cat("  unlisted cells, dropped from did_results_fig.csv.gz\n")
  print(cmp[unlisted, c(sup_key, "cap_use_auto")], row.names = FALSE)
}

lim <- CAPS_USE[!duplicated(rk), c(sup_key, "cap_use", "cap_lo_use")]

FIG <- merge(FIG, lim, by = sup_key, all.x = TRUE)

PRE_GATED <- TRUE

pre_lo <- unname(ca_fig$pretrend[["lo"]])

cal_lo <- vapply(FIG$outcome, function(o)
  unname(ca_fig$caps$pa$lo[[o]]), numeric(1))
cal_hi <- vapply(FIG$outcome, function(o)
  unname(ca_fig$caps$pa$hi[[o]]), numeric(1))

is_cal <- FIG$aggregation %in% CAL_AGG

## The event lower bound. A cell with no supported pre horizon at all takes 0,
## so it contributes to the post-treatment figures and not to Figure S6.
ev_lo <- if (PRE_GATED) {
  ifelse(is.na(FIG$cap_lo_use), 0, pmax(pre_lo, FIG$cap_lo_use))
} else {
  rep(pre_lo, nrow(FIG))
}

keep <- ifelse(
  is_cal,
  FIG$time >= cal_lo & FIG$time <= cal_hi,
  !is.na(FIG$cap_use) & FIG$time >= ev_lo & FIG$time <= FIG$cap_use)
keep[is.na(keep)] <- FALSE

ca_log("windows applied: kept ", format(sum(keep), big.mark = ","), " of ",
       format(nrow(FIG), big.mark = ","), " figure rows")

FIG <- FIG[keep, , drop = FALSE]
FIG <- FIG[order(FIG$estimator,
                 match(FIG$pairing, ca_pairings$pairing),
                 match(FIG$outcome, ca_outcomes$outcome),
                 FIG$stratum, FIG$ecoregion, FIG$aggregation, FIG$time), ]
rownames(FIG) <- NULL

RES_OUT <- ca_meta("did_results_fig.csv.gz")
con <- gzfile(RES_OUT, "w")
utils::write.csv(FIG, con, row.names = FALSE)
close(con)

ca_log("Wrote ", RES_OUT, "  ", format(nrow(FIG), big.mark = ","),
       " of ", format(nrow(R), big.mark = ","), " result rows")

cat("\nROWS BY FIGURE\n")
print(sort(table(unlist(strsplit(FIG$figures, "; "))), decreasing = TRUE))

cat("\nROWS BY ESTIMATOR AND AGGREGATION\n")
print(table(FIG$estimator, FIG$aggregation))

## BOTH SCALES, AFTER WINDOWING.
##
## The check that a Results sentence can be written. paired is the horizons
## carrying percent and native for the same cell, and it should equal drawn for
## GPP and the three AGB products and be zero for the three signed fluxes,
## which have no log scale to pair with.
cat("\nSCALES HELD PER OUTCOME, after windowing\n")
cell_time <- do.call(paste, c(FIG[, c("estimator", "pairing", "outcome",
                                      "stratum", "ecoregion", "aggregation",
                                      "time")], sep = "|"))
sc <- do.call(rbind, lapply(split(seq_len(nrow(FIG)), FIG$outcome),
                            function(i) {
  z <- FIG[i, , drop = FALSE]
  k <- cell_time[i]
  data.frame(outcome = z$outcome[1],
             drawn   = sum(z$fig_row),
             percent = sum(z$transform == "log1p"),
             native  = sum(z$transform != "log1p"),
             paired  = sum(table(k) == 2L),
             stringsAsFactors = FALSE)
}))
print(sc[order(match(sc$outcome, ca_outcomes$outcome)), ], row.names = FALSE)

cat("\nEVENT ROWS EITHER SIDE OF TREATMENT\n")
ev <- FIG[!(FIG$aggregation %in% CAL_AGG), , drop = FALSE]
print(table(ev$estimator, ifelse(ev$time < 0, "pre", "post")))

## Every pairing that carries both aggregations, which UP_UNP must.
cat("\nAGGREGATIONS HELD PER PAIRING\n")
print(table(FIG$pairing, FIG$aggregation))

## The span written per statewide cell, the file's own record of what the
## figures can draw. A first_time short of the pre-trend floor is PRE_GATED
## holding the cell to its supported pre window.
cat("\nSPAN WRITTEN, statewide cells\n")
sw <- FIG[FIG$ecoregion == "", , drop = FALSE]
sp <- stats::aggregate(time ~ estimator + pairing + outcome + stratum +
                         aggregation, data = sw,
                       FUN = function(x) c(first = min(x), last = max(x),
                                           n = length(x)))
print(do.call(data.frame, sp), row.names = FALSE)
