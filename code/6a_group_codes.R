### CA carbon revision pipeline
### Stage 6a. Class codes, analysis groups, exclusion ledger
###
### Stage 5 classified nothing. Every analytical rule lives in ca_thin_rules,
### ca_fire_rules, ca_offset_rules, ca_pa_rules, ca_window_rules, ca_screen,
### ca_cpad_agency, and ca_ownership_level, and every one of them is applied
### here, once, in the open.
###
### TWO OUTPUTS FROM ONE PASS
###
### class_code is complete. Every pixel receives one, there are no missing
### values, and repeat entry, reburn, and thinning-then-fire histories all
### appear as ordered token chains. Nothing in this manuscript is estimated on
### them.
###
### analysis_group is the 16 codes this manuscript estimates on, plus NA.
### eligibility carries the name of the rule that removed a pixel, so the
### exclusion ledger is a crosstab of class_code against analysis_group rather
### than a separate accounting step.
###
### Tokens, status levels, group definitions, and rule order are all in
### ca_config.R sections 0.7b and 0.8. This script applies them and counts.
###
### Usage
###   Rscript 6a_group_codes.R tasks      print the array size and exit
###   Rscript 6a_group_codes.R probe      non-disturbing test, writes meta
###   Rscript 6a_group_codes.R            run one array task
###   Rscript 6a_group_codes.R archive    copy completed output to /projects

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
mode <- if (length(args) && args[1] %in% c("tasks", "probe", "archive")) {
  args[1]
} else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

YR_FIRST <- min(ca_years$analysis)
YR_LAST <- max(ca_years$analysis)

read_events <- function(arm, lulc, cols) {
  p <- ca_event_path(arm, lulc)
  if (!file.exists(p)) stop("Missing stage 5a arm: ", p)
  as.data.frame(arrow::read_parquet(p, col_select = dplyr::all_of(cols)))
}

grid_ids <- function(lulc) {
  part <- file.path(ca_grid$parquet_dir, paste0("lulc=", lulc))
  if (!dir.exists(part)) stop("Grid partition not found: ", part)
  ds <- arrow::open_dataset(part)
  v <- as.data.frame(
    arrow::Scanner$create(ds, projection = "pixel_id")$ToTable())$pixel_id
  v <- as.integer(v)
  if (is.unsorted(v)) v <- sort(v)
  v
}

# Pixel-year key. Numeric rather than a pasted string, because the thinning
# table runs to fifty million rows on Evergreen and string keys would triple it.
ckey <- function(idx, year) as.double(idx) * 1000 + (as.double(year) - 1800)

# Guard against labels in the data that have no match in the config table. A
# label that falls through silently is how 42,425 tribal pixels ended up in the
# control pool on the first verified run.
undeclared <- function(values, roles, layer_name) {
  bad <- unique(values[is.na(roles)])
  if (length(bad)) {
    stop(layer_name, " labels not declared in config: ",
         paste(bad, collapse = ", "),
         ". Add them to the config table before rerunning.")
  }
}


# ---------------------------------------------------------------------------
# MODE tasks / archive
# ---------------------------------------------------------------------------

if (mode == "tasks") {
  cat("Array size: ", nrow(ca_lulc), "\n", sep = "")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_work("groups"), "groups", subdir = "work")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# MODE probe. The non-disturbing test
# ---------------------------------------------------------------------------
# ca_thin_rules keeps non_disturbing in the control pool by default and states
# that the default is tested rather than assumed. Six categories are the
# reason. Three alter standing vegetation by definition, Range Cover Type
# Conversion, Range Cover Manipulation, and Range Forage Improvement. Three are
# indirect markers of operations rather than impacts, Range Piling Slash, Cover
# brush pile for burning, and Leave Trees for Wildlife Reasons, since slash and
# retention designations only exist where cutting occurred.
#
# The test. For pixels whose only record is non_disturbing, measure the
# Disturbance_TreeFrac flag rate in the record year plus or minus the thinning
# detection window against reference pixels carrying no record of any kind,
# weighted to the activity's own ecoregion and year distribution.
#
# Evergreen only. It carries 92 percent of the grid and Deciduous and Mixed
# hold too few non-disturbing records to test separately. The decision applies
# to all three classes.

if (mode == "probe") {

  lulc <- ca_lulc$label[2]
  ca_log("Stage 6a probe, non-disturbing test, class ", lulc)

  gid <- grid_ids(lulc)
  n <- length(gid)

  th <- read_events("thin", lulc,
                    c("pixel_id", "event_year", "record_class", "key_value"))
  th$idx <- match(th$pixel_id, gid)
  th <- th[!is.na(th$idx), ]

  n_rec <- tabulate(th$idx, nbins = n)
  n_nd <- tabulate(th$idx[th$record_class %in% "non_disturbing"], nbins = n)
  nd_only <- n_rec > 0L & n_rec == n_nd

  fi <- read_events("fire", lulc, c("pixel_id"))
  has_fire <- logical(n)
  j <- match(fi$pixel_id, gid)
  has_fire[j[!is.na(j)]] <- TRUE
  rm(fi, j); ca_gc()

  nd_only <- nd_only & !has_fire
  ref <- n_rec == 0L & !has_fire
  ca_log("  non-disturbing-only pixels ", format(sum(nd_only), big.mark = ","),
         ", reference pixels ", format(sum(ref), big.mark = ","))

  eco <- ca_read_static("ecoregion_l3", lulc)
  eco_v <- rep(NA_character_, n)
  j <- match(eco$pixel_id, gid)
  eco_v[j[!is.na(j)]] <- eco$value[!is.na(j)]
  rm(eco, j); ca_gc()

  sc <- as.data.frame(arrow::read_parquet(
    ca_screen_path(lulc), col_select = dplyr::all_of(c("pixel_id", "year"))))
  sc$idx <- match(sc$pixel_id, gid)
  sc <- sc[!is.na(sc$idx), c("idx", "year")]

  # A detection at year y covers the windows centred on y-w to y+w.
  w <- ca_screen$thin_detection_window
  cov_idx <- rep(sc$idx, times = 2L * w + 1L)
  cov_yr <- as.vector(vapply(seq(-w, w), function(d) sc$year + d,
                             numeric(nrow(sc))))
  cov_key <- ckey(cov_idx, cov_yr)
  keep <- !duplicated(cov_key)
  cov_idx <- cov_idx[keep]; cov_yr <- cov_yr[keep]; cov_key <- cov_key[keep]
  rm(sc, keep); ca_gc()

  ref_n <- table(eco_v[ref][!is.na(eco_v[ref])])
  in_ref <- ref[cov_idx]
  re <- eco_v[cov_idx[in_ref]]
  ry <- cov_yr[in_ref]
  ok <- !is.na(re)
  ref_hit <- table(re[ok], ry[ok])
  rm(in_ref, re, ry, ok); ca_gc()

  ref_rate <- function(e, y) {
    i <- match(e, rownames(ref_hit))
    j <- match(as.character(y), colnames(ref_hit))
    d <- as.numeric(ref_n[e])
    out <- rep(NA_real_, length(e))
    good <- !is.na(i) & !is.na(j) & !is.na(d) & d > 0
    out[good] <- ref_hit[cbind(i[good], j[good])] / d[good]
    out
  }

  a <- th[nd_only[th$idx] & th$record_class %in% "non_disturbing" &
            !is.na(th$event_year), c("idx", "event_year", "key_value")]
  a$eco <- eco_v[a$idx]
  a <- a[!is.na(a$eco), ]
  a$hit <- ckey(a$idx, a$event_year) %in% cov_key
  a$exp <- ref_rate(a$eco, a$event_year)
  a <- a[!is.na(a$exp), ]

  res <- a %>%
    group_by(activity = key_value) %>%
    summarise(n_records = n(),
              n_pixels = length(unique(idx)),
              obs_rate = mean(hit),
              exp_rate = mean(exp),
              .groups = "drop") %>%
    mutate(ratio = ifelse(exp_rate > 0, obs_rate / exp_rate, NA_real_),
           p_value = mapply(function(k, m, p) {
             if (is.na(p) || p <= 0 || p >= 1) return(NA_real_)
             stats::binom.test(k, m, p, alternative = "greater")$p.value
           }, round(obs_rate * n_records), n_records, exp_rate),
           decision = ifelse(!is.na(ratio) & ratio >= 2 & n_pixels >= 100 &
                               !is.na(p_value) & p_value < 0.001,
                             "disturbing", "non_disturbing")) %>%
    arrange(desc(ratio))

  dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res, ca_nondist_file, row.names = FALSE)
  ca_log("Wrote ", ca_nondist_file)
  print(as.data.frame(res))
  ca_log("Demoted to disturbing: ",
         paste(res$activity[res$decision == "disturbing"], collapse = "; "))

  ca_stamp("6a_group_codes_probe")
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# MODE run
# ---------------------------------------------------------------------------

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id must be 1 to ", nrow(ca_lulc), ", got ", task)
}
lulc <- ca_lulc$label[task]

out_path <- ca_group_path(lulc)
if (file.exists(out_path) && !OVERWRITE) {
  ca_log("Present, skipping: ", basename(out_path))
  quit(save = "no")
}

DEMOTED <- ca_nondisturbing_classes()
ca_log("Stage 6a groups, class ", lulc, ". Activities demoted to disturbing: ",
       if (length(DEMOTED)) paste(DEMOTED, collapse = "; ") else "none")

gid <- grid_ids(lulc)
n <- length(gid)
ca_log("  grid pixels ", format(n, big.mark = ","))

ledger <- list()
note <- function(rule, removed, remaining) {
  ledger[[length(ledger) + 1]] <<- data.frame(
    step = length(ledger) + 1L, rule = rule,
    removed = as.numeric(removed), remaining = as.numeric(remaining),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 1. THINNING TOKENS
# ---------------------------------------------------------------------------

th <- read_events("thin", lulc, c("pixel_id", "event_year", "class_value",
                                  "record_class", "key_value"))
th$idx <- match(th$pixel_id, gid)
th <- th[!is.na(th$idx), ]
ca_log("  thin records ", format(nrow(th), big.mark = ","))

# A record with no record_class would be dropped silently by every rule below.
# Section 8d of 3_verify_extract.R reports pct_unresolved and it must be zero,
# so this is a guard rather than a case to handle.
if (any(is.na(th$record_class))) {
  stop("Thinning records with no record_class: ",
       sum(is.na(th$record_class)),
       ". Resolve the crosswalk before running stage 6.")
}

tok <- rep(NA_character_, nrow(th))
tok[th$record_class %in% "wildfire"] <- "Fu"
tok[th$record_class %in% "disturbance_only"] <- "Td"
tok[th$record_class %in% "non_disturbing"] <- "Tn"

eli <- th$record_class %in% "treatment_LMH"
m <- match(th$class_value, ca_thin_intensity$label)
tok[eli & !is.na(m)] <- ca_thin_intensity$group[m[eli & !is.na(m)]]
tok[eli & is.na(m)] <- "Tv"

# The probe outcome, applied here and nowhere else. The sensitivity flag is the
# conservative bound on the categories the probe could not isolate, since a
# category that never occurs alone cannot be tested against a clean reference.
tok[tok %in% "Tn" & th$key_value %in% DEMOTED] <- "Td"
if (isTRUE(ca_thin_rules$non_disturbing_as_disturbing)) {
  ca_log("  SENSITIVITY: all non_disturbing records treated as disturbing")
  tok[tok %in% "Tn"] <- "Td"
}
th$token <- tok
rm(tok, eli, m); ca_gc()

print(table(th$token, useNA = "ifany", dnn = "thin tokens"))

# Undated records carry no position in a chain. Flagged and set aside.
undated_thin <- logical(n)
u <- is.na(th$event_year)
undated_thin[unique(th$idx[u])] <- TRUE
ca_log("  undated thinning records ", format(sum(u), big.mark = ","),
       " on ", format(sum(undated_thin), big.mark = ","), " pixels")
th <- th[!u, c("idx", "event_year", "token")]
rm(u); ca_gc()

# Fu arrives through the thinning layer but it is fire evidence, so it is
# carried to the fire arm rather than competing with thinning tokens. A pixel
# with a thinning record and a wildfire record in one year has one entry and
# one fire, which is two events and is meant to be.
th <- th[, c("idx", "event_year", "token")]
th$arm_id <- ifelse(substr(th$token, 1, 1) == "F", 2L, 1L)


# ---------------------------------------------------------------------------
# 2. FIRE TOKENS
# ---------------------------------------------------------------------------

fi <- read_events("fire", lulc, c("pixel_id", "event_year", "class_value"))
fi$idx <- match(fi$pixel_id, gid)
fi <- fi[!is.na(fi$idx), ]

sev <- suppressWarnings(as.integer(fi$class_value))
g <- ca_fire_severity$group[match(sev, ca_fire_severity$code)]
r <- ca_fire_severity$role[match(sev, ca_fire_severity$code)]
ftok <- rep(NA_character_, nrow(fi))
ftok[r %in% "treatment"] <- g[r %in% "treatment"]
ftok[r %in% "disturbance_only"] <- "Fu"
fi$token <- ftok
fi <- fi[!is.na(fi$token) & !is.na(fi$event_year),
         c("idx", "event_year", "token")]
fi$arm_id <- 2L
rm(sev, g, r, ftok); ca_gc()

ca_log("  fire records ", format(nrow(fi), big.mark = ","))
print(table(fi$token, dnn = "fire tokens"))

# ---------------------------------------------------------------------------
# 3. RECORDS TO ENTRIES
# ---------------------------------------------------------------------------
# Both source layers are activity-level rather than event-level. FACTS writes
# one row per activity, so a single stand entry produces a harvest row, a site
# preparation row, and a slash row. MTBS writes one row per perimeter year, and
# a fire recorded through a FACTS wildfire activity line appears again there.
# Counting rows as events reads one entry as several and discards the pixel,
# which is what removed 18,646,561 pixels on the first run. 6,253,412 of those
# carried at most one measured severity and were lost to a duplicate Fu.
#
# Records at the same pixel, in the same year, in the same arm are one event.
# The surviving token is the one carrying the most information about magnitude,
# ca_token_info, so a measured severity beats Fu and a measured intensity beats
# Td. Records in different years stay separate events, because a stand entered
# in 2005 and again in 2012 was entered twice.
#
# This also retires the pixel-year intensity conflict. Two eligible records
# disagreeing in one year resolve to the higher intensity rather than to Tv.

ev <- rbind(th, fi)
rm(th, fi); ca_gc()

n_rec_in <- nrow(ev)
ev$info <- ca_token_info[ev$token]
ev$info[is.na(ev$info)] <- 99L
key <- ckey(ev$idx, ev$event_year) * 10 + ev$arm_id
o <- order(key, ev$info)
ev <- ev[o, ]; key <- key[o]
ev <- ev[!duplicated(key), ]
rm(o, key); ca_gc()
ca_log("  records ", format(n_rec_in, big.mark = ","),
       " collapsed to entries ", format(nrow(ev), big.mark = ","),
       " (", round(100 * (1 - nrow(ev) / n_rec_in), 1), " percent merged)")

# Thinning before fire within a year, arbitrary and stated.
ev <- ev[order(ev$idx, ev$event_year, ev$arm_id, ev$token), ]

# YEAR GAP DIAGNOSTIC
# Same-year merging is the defensible floor. Whether the window should widen
# depends on how far apart the surviving same-arm entries actually sit, which
# is a fact rather than a preference. A FACTS project whose activity rows span
# a fiscal year boundary shows up as a spike at gap 1.
gap <- diff(ev$event_year)
same <- ev$idx[-1] == ev$idx[-nrow(ev)] & ev$arm_id[-1] == ev$arm_id[-nrow(ev)]
gp <- gap[same]
ap <- ev$arm_id[-1][same]
ca_log("  consecutive same-arm entry gaps, thinning then fire")
print(table(ifelse(ap == 1L, "thin", "fire"),
            pmin(gp, 6L), dnn = c("arm", "gap years, 6 is 6 or more")))
rm(gap, same, gp, ap); ca_gc()

# ADJACENT-YEAR MERGE, TOKENS CARRYING NO MAGNITUDE ONLY
# The gap profile above justifies this for fire and not for thinning. Evergreen
# fire gaps 2 to 5 average 121,636 and gap 1 is 678,391, an excess near 557,000
# that is the fiscal-year and perimeter-year artifact. Thinning decays smoothly
# at a ratio near 0.70 across gaps 2 to 5, which is what genuine repeat entry
# looks like.
#
# So the merge is restricted rather than applied by arm. A token carrying no
# magnitude, Fu for fire and Td for thinning, merges into a measured neighbour
# one year away. Two measured tokens never merge, so a 2005 low-severity fire
# followed by a 2006 moderate-severity fire stays two events, and a 2005 medium
# thin followed by a 2006 high thin stays two entries. Fu against Fu and Td
# against Td collapse to one.
#
# Single pass. A chain of three no-magnitude rows leaves one behind, which is
# rare and errs toward keeping events rather than merging them away.
nm <- ev$token %in% ca_tokens_no_magnitude
m1 <- nrow(ev)
adj_prev <- c(FALSE, ev$idx[-1] == ev$idx[-m1] & ev$arm_id[-1] == ev$arm_id[-m1] &
                (ev$event_year[-1] - ev$event_year[-m1]) <= 1L)
adj_next <- c(adj_prev[-1], FALSE)
mag_prev <- c(FALSE, !nm[-m1])
mag_next <- c(!nm[-1], FALSE)
drop_nm <- nm & ((adj_prev & mag_prev) | (adj_next & mag_next) |
                   (adj_prev & !mag_prev))
ev <- ev[!drop_nm, ]
ca_log("  adjacent-year merge removed ", format(sum(drop_nm), big.mark = ","),
       " no-magnitude entries, ", format(nrow(ev), big.mark = ","), " remain")
rm(nm, adj_prev, adj_next, mag_prev, mag_next, drop_nm); ca_gc()

cnt <- tabulate(ev$idx, nbins = n)
chain <- rep("U", n)

# Most pixels carry exactly one event, and those need no string operation.
one <- cnt == 1L
pos <- match(which(one), ev$idx)
chain[one] <- ev$token[pos]
rm(pos)

many <- which(cnt > 1L)
if (length(many)) {
  sub <- ev[ev$idx %in% many, ]
  agg <- sub %>%
    group_by(idx) %>%
    summarise(chain = paste0(token, collapse = ""), .groups = "drop")
  chain[agg$idx] <- agg$chain
  rm(sub, agg); ca_gc()
}
ca_log("  pixels with events, single ", format(sum(one), big.mark = ","),
       ", multiple ", format(length(many), big.mark = ","))
rm(one, many)

first_year <- rep(NA_integer_, n)
last_year <- rep(NA_integer_, n)
fp <- !duplicated(ev$idx)
lp <- !duplicated(ev$idx, fromLast = TRUE)
first_year[ev$idx[fp]] <- ev$event_year[fp]
last_year[ev$idx[lp]] <- ev$event_year[lp]
rm(fp, lp)

n_thin <- tabulate(ev$idx[substr(ev$token, 1, 1) == "T"], nbins = n)
n_fire <- tabulate(ev$idx[substr(ev$token, 1, 1) == "F"], nbins = n)

# Eligibility counts tokens flagged counts_as_event in ca_class_tokens and
# dated at or before the last analysis year. Future records describe intent
# rather than history and do not disqualify a control, per ca_window_rules.
counts <- ca_class_tokens$token[ca_class_tokens$counts_as_event]
eff <- ev$event_year <= YR_LAST & ev$token %in% counts
n_eff <- tabulate(ev$idx[eff], nbins = n)

sole_token <- rep(NA_character_, n)
sole_year <- rep(NA_integer_, n)
sub <- ev[eff, ]
so <- tabulate(sub$idx, nbins = n) == 1L
pos <- match(which(so), sub$idx)
sole_token[so] <- sub$token[pos]
sole_year[so] <- sub$event_year[pos]
rm(sub, so, pos, eff, ev, counts); ca_gc()

# ---------------------------------------------------------------------------
# 4. STATUS
# ---------------------------------------------------------------------------
# Evaluated in the order of ca_status_levels, lowest priority first, so the
# later assignment overwrites the earlier one.

status <- rep("NP", n)

# An undeclared label must not fall through to NP. match() on a value absent
# from ca_ownership_level or ca_cpad_agency returns NA, every role test then
# fails, and the pixel lands in the control pool without ever being examined.
# The stage 6 verifier found 5,807 UNP pixels carrying a CPAD record this way.
# Undeclared records now take status UN, which is excluded from both arms, and
# the offending values are logged so they can be declared in config.
undeclared <- function(vals, roles, layer) {
  u <- vals[is.na(roles)]
  if (!length(u)) return(invisible(NULL))
  tb <- sort(table(u), decreasing = TRUE)
  ca_log("  UNDECLARED ", layer, " labels, assigned status UN and excluded:")
  for (k in seq_along(tb)) {
    ca_log("    ", names(tb)[k], "  ", format(as.integer(tb[k]),
                                              big.mark = ","), " records")
  }
}

own <- read_events("own", lulc, c("pixel_id", "class_value"))
own$idx <- match(own$pixel_id, gid)
own <- own[!is.na(own$idx), ]
r <- ca_ownership_level$role[match(own$class_value, ca_ownership_level$label)]
undeclared(own$class_value, r, "ownership")
status[unique(own$idx[r %in% "exclude_both"])] <- "NGO"
status[unique(own$idx[r %in% "public"])] <- "PUB"

# The ownership layer carries its own tribal role, and the separate tribal
# layer does not cover all of it. 36,618 tribal-owned pixels sat in the control
# pool on the first run because only the tribal layer was consulted.
own_tribal <- unique(own$idx[r %in% "tribal"])
rm(own, r); ca_gc()

pa <- read_events("pa", lulc, c("pixel_id", "event_year", "class_value"))
pa$idx <- match(pa$pixel_id, gid)
pa <- pa[!is.na(pa$idx), ]
role <- ca_cpad_agency$role[match(pa$class_value, ca_cpad_agency$label)]
undeclared(pa$class_value, role, "CPAD agency")

# CPAD tribal is protected at a non-Federal non-State level, so it is
# protected_other by the same logic as County, City, and Non Profit. 5,807
# pixels sat in the control pool on the first run because only pa_treatment and
# protected_other were handled.
status[unique(pa$idx[role %in% "tribal"])] <- "PO"
status[unique(pa$idx[role %in% "protected_other"])] <- "PO"
status[unique(pa$idx[role %in% "pa_treatment"])] <- "P"

# Earliest establishment year across superunits, which is when protection
# began. NA is protection of unknown date, not absence of protection.
pa_est <- rep(NA_integer_, n)
trt <- pa[role %in% "pa_treatment", ]
if (nrow(trt)) {
  d <- trt[!is.na(trt$event_year), ]
  d <- d[order(d$idx, d$event_year), ]
  f <- !duplicated(d$idx)
  pa_est[d$idx[f]] <- d$event_year[f]
  rm(d, f)
}
pa_undated <- logical(n)
pa_undated[unique(trt$idx[is.na(trt$event_year)])] <- TRUE
pa_undated <- pa_undated & is.na(pa_est)
rm(pa, role, trt); ca_gc()

# Tribal last, from both sources. It is an exclusion mask and nothing overrides
# it.
tr <- read_events("tribal", lulc, c("pixel_id"))
j <- match(tr$pixel_id, gid)
status[unique(c(j[!is.na(j)], own_tribal))] <- "TR"
rm(tr, j, own_tribal); ca_gc()

print(table(status, dnn = "status"))

off <- read_events("offset", lulc, c("pixel_id", "event_year"))
off$idx <- match(off$pixel_id, gid)
off <- off[!is.na(off$idx) & !is.na(off$event_year), ]
off <- off[order(off$idx, off$event_year), ]
f <- !duplicated(off$idx)
offset_start <- rep(NA_integer_, n)
offset_start[off$idx[f]] <- off$event_year[f]
rm(off, f); ca_gc()

# ---------------------------------------------------------------------------
# 5. SCREEN AND DETECTION
# ---------------------------------------------------------------------------

sc <- as.data.frame(arrow::read_parquet(
  ca_screen_path(lulc), col_select = dplyr::all_of(c("pixel_id", "year"))))
sc$idx <- match(sc$pixel_id, gid)
sc <- sc[!is.na(sc$idx), ]

n_detect <- tabulate(sc$idx, nbins = n)

# Window is arm-specific. Archival completion dates should sit close to the
# detection, whereas fire attribution can fall a year or more away.
w_vec <- ifelse(substr(sole_token, 1, 1) == "F",
                ca_screen$fire_detection_window,
                ca_screen$thin_detection_window)
detect_at_event <- logical(n)
d <- abs(sc$year - sole_year[sc$idx])
hit <- !is.na(d) & d <= w_vec[sc$idx]
detect_at_event[unique(sc$idx[hit])] <- TRUE
rm(sc, d, hit, w_vec); ca_gc()

ca_log("  pixels with any detection ",
       format(sum(n_detect > 0L), big.mark = ","), " (",
       round(100 * mean(n_detect > 0L), 2), " percent)")

# ---------------------------------------------------------------------------
# 6. CLASS CODE
# ---------------------------------------------------------------------------
# The undated marker is written last because an undated record has no position
# in time and cannot be ordered against dated ones.

class_code <- paste0(chain, ifelse(undated_thin, "Tu", ""), "_", status)

# ---------------------------------------------------------------------------
# 7. ELIGIBILITY
# ---------------------------------------------------------------------------
# Order is ca_exclusion_order. The first rule that fires owns the pixel, so the
# ledger reads as a sequence of disjoint removals rather than overlapping
# counts.

elig <- rep(NA_character_, n)
open <- rep(TRUE, n)

fire_rule <- function(name, hit) {
  if (!name %in% ca_exclusion_order) {
    stop("Rule not declared in ca_exclusion_order: ", name)
  }
  hit <- hit & !is.na(hit)
  h <- open & hit
  elig[h] <<- name
  open[h] <<- FALSE
  note(name, sum(h), sum(open))
}

note("start", 0, n)

fire_rule("tribal", status == "TR")
fire_rule("pa_undated", pa_undated)
fire_rule("status_not_analysed", !status %in% c("P", "NP"))
fire_rule("undated_thin_record", undated_thin)
fire_rule("multi_event", n_eff > 1L)
fire_rule("unclassifiable_event",
          n_eff == 1L & !sole_token %in% ca_tokens_treatment)
fire_rule("pre_window_event", n_eff == 1L & sole_year < YR_FIRST)

# Protection is a time-varying state, per ca_pa_rules$treatment_time_varying,
# so a pixel that burned in 1995 inside a PA established in 2010 was not in a
# PA when it burned. Calling it F2P would describe its status today rather than
# at the event. Fire and thinning arms only, since UP_UNP carries no cohort
# year and is specified in calendar time.
if (isTRUE(ca_pa_rules$label_at_cohort_year)) {
  fire_rule("pa_after_cohort",
            status == "P" & n_eff == 1L & sole_token %in% ca_tokens_treatment &
              sole_year >= YR_FIRST & sole_year <= YR_LAST &
              !is.na(pa_est) & pa_est > sole_year)
}

# Offsets sit on undisturbed forest by definition, so an offset pixel carrying
# an event is a conflict rather than a treatment. Offsets on protected land are
# not excluded, they become OP, which is the audit of the Methods claim that
# projects were restricted to non-protected private land.
fire_rule("offset_event_conflict", !is.na(offset_start) & n_eff == 1L)
fire_rule("offset_out_of_window",
          !is.na(offset_start) &
            (offset_start < YR_FIRST | offset_start > YR_LAST))

# A polygon record asserts work inside a boundary, not work on every pixel
# inside it. The detection requirement is the correction for that
# over-assignment and is the submitted specification.
cand_group <- rep(NA_character_, n)
if (isTRUE(ca_screen$require_treatment_detection)) {
  # Captured before the rule fires, so the summary can report how many recorded
  # events survived the detection requirement per group. That retention is the
  # number the Methods quotes.
  cg <- open & n_eff == 1L & sole_token %in% ca_tokens_treatment &
    sole_year >= YR_FIRST & sole_year <= YR_LAST
  cand_group[cg] <- paste0(sole_token[cg],
                           ifelse(status[cg] == "P", "P", "NP"))
  # WHY A CANDIDATE FAILED THE DETECTION TEST
  # Two causes with opposite remedies. A pixel inside a polygon that shows no
  # canopy loss in any year of the record was probably never treated, and
  # dropping it is the correction the requirement exists for. A pixel that
  # shows loss in some other year was treated on a date the archive records
  # badly, and the remedy there is a wider window, not exclusion.
  miss <- !is.na(cand_group) & !detect_at_event
  dg <- data.frame(group = cand_group[miss],
                   ever = n_detect[miss] > 0L, stringsAsFactors = FALSE)
  dg <- dg %>%
    group_by(group) %>%
    summarise(n_missed = n(),
              pct_never_disturbed = round(100 * mean(!ever), 1),
              pct_disturbed_off_window = round(100 * mean(ever), 1),
              .groups = "drop")
  ca_log("  detection failures, never disturbed against wrong window")
  print(as.data.frame(dg))
  utils::write.csv(dg, ca_meta(sprintf("groups_detection_miss_%s.csv", lulc)),
                   row.names = FALSE)
  rm(cg, miss, dg); ca_gc()

  fire_rule("no_detection_at_event",
            n_eff == 1L & sole_token %in% ca_tokens_treatment &
              !detect_at_event)
}

# The control screen. Applies only to pixels with no event.
fire_rule("screen_detection", n_eff == 0L & n_detect > 0L)

# ---------------------------------------------------------------------------
# 8. ANALYSIS GROUP
# ---------------------------------------------------------------------------

grp <- rep(NA_character_, n)
suffix <- ifelse(status == "P", "P", "NP")

is_treat <- open & n_eff == 1L & sole_token %in% ca_tokens_treatment &
  sole_year >= YR_FIRST & sole_year <= YR_LAST
grp[is_treat] <- paste0(sole_token[is_treat], suffix[is_treat])

is_offset <- open & n_eff == 0L & !is.na(offset_start)
grp[is_offset] <- paste0("O", suffix[is_offset])

is_undist <- open & n_eff == 0L & is.na(offset_start)
grp[is_undist] <- paste0("U", suffix[is_undist])

bad <- !is.na(grp) & !grp %in% ca_analysis_groups$group
if (any(bad)) {
  stop("Group codes not declared in ca_analysis_groups: ",
       paste(unique(grp[bad]), collapse = ", "))
}

i <- match(grp, ca_analysis_groups$group)
role <- ca_analysis_groups$role[i]
arm <- ca_analysis_groups$arm[i]
stratum <- ca_analysis_groups$stratum[i]

elig[!is.na(role)] <- role[!is.na(role)]
elig[is.na(elig)] <- "unassigned"

note("assigned", sum(!is.na(grp)), sum(is.na(grp) & open))

cohort_year <- rep(NA_integer_, n)
cohort_year[is_treat] <- sole_year[is_treat]
cohort_year[is_offset] <- offset_start[is_offset]

# ---------------------------------------------------------------------------
# 9. WRITE
# ---------------------------------------------------------------------------

out <- data.frame(
  pixel_id = gid,
  lulc = lulc,
  class_code = class_code,
  event_chain = chain,
  status_code = status,
  analysis_group = grp,
  arm = arm,
  stratum = stratum,
  eligibility = elig,
  cohort_year = cohort_year,
  pa_est_year = pa_est,
  offset_start = offset_start,
  n_events = as.integer(n_thin + n_fire),
  n_thin = as.integer(n_thin),
  n_fire = as.integer(n_fire),
  n_events_effective = as.integer(n_eff),
  first_event_year = first_year,
  last_event_year = last_year,
  undated_thin = undated_thin,
  n_detect = as.integer(n_detect),
  detect_at_event = detect_at_event,
  stringsAsFactors = FALSE
)

arrow::write_parquet(arrow::as_arrow_table(out), out_path,
                     compression = ca_io$compression)
ca_log("Wrote ", basename(out_path), "  ",
       format(nrow(out), big.mark = ","), " rows  ",
       round(file.size(out_path) / 1e6, 1), " MB")

led <- do.call(rbind, ledger)
lp <- ca_meta(sprintf("groups_ledger_%s.csv", lulc))
utils::write.csv(led, lp, row.names = FALSE)
ca_log("Exclusion ledger: ", lp)
print(led)

# The taxonomy. class_code against analysis_group is the ledger in full detail
# and the artifact the Methods and the response letter quote from.
tax <- out %>%
  group_by(class_code, status_code, eligibility,
           group = ifelse(is.na(analysis_group), "-", analysis_group)) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
tp <- ca_meta(sprintf("groups_taxonomy_%s.csv", lulc))
utils::write.csv(tax, tp, row.names = FALSE)
ca_log("Class code taxonomy: ", tp, "  ",
       format(nrow(tax), big.mark = ","), " distinct rows")

cand_n <- as.data.frame(table(cand_group), stringsAsFactors = FALSE)
names(cand_n) <- c("analysis_group", "n_candidate")

gs <- out %>%
  filter(!is.na(analysis_group)) %>%
  group_by(analysis_group, arm) %>%
  summarise(n = n(),
            n_cohorts = length(unique(cohort_year[!is.na(cohort_year)])),
            .groups = "drop") %>%
  left_join(cand_n, by = "analysis_group") %>%
  mutate(pct_detected = ifelse(is.na(n_candidate), NA_real_,
                               round(100 * n / n_candidate, 1))) %>%
  arrange(analysis_group)
print(as.data.frame(gs))
utils::write.csv(gs, ca_meta(sprintf("groups_summary_%s.csv", lulc)),
                 row.names = FALSE)

ca_stamp("6a_group_codes")
ca_log("Done, class ", lulc)
