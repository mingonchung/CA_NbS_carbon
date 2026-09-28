### CA carbon revision pipeline
### Stage 8. Matching
###
### One script replacing 8_1_matching_UP_UNP_NN.R, 8_2_matching_FP_FNP_NN.R,
### 8_3_matching_TP_TNP_NN.R, and 8_5_matching_ONP_UNP_NN.R.
###
### WHAT THIS STAGE DOES
###
### Applies the 100,000 cap, assembles the matching input, computes the
### pre-treatment covariates per cohort, runs the 11-config sweep for one
### covariate arm, and writes matched data plus balance metrics. It estimates
### nothing. Outcomes for matched pixels are read at stage 9.
###
### THE ARRAY
###
### 18 tasks, pairing crossed with arm, from ca_match_tasks(). Classes pool
### inside a task because lulc is an exact matching variable, which is what the
### submitted run did when it fed one combined input to matchit().
###
### THE TWO SIDES
###
### ca_pool_path() resolves the file per role. Treated groups come spaced from
### stage 7, control pools come whole from stage 6. UP holds both roles and
### resolves to a different file in each. Reading the spaced file for a control
### pool would not fail, it would quietly shrink UP by a factor near 25 and
### inflate every standard error in the protected fire and thinning arms, so no
### path in this script is built any other way.
###
### COHORT STRATIFICATION
###
### ca_pretrt$control_rule sets the control window as the same window as the
### matched treatment pixel, which is circular unless matching runs cohort by
### cohort. See section 0.11 of ca_config.R for why the alternatives are worse
### rather than merely slower. Applied to all three arms so the arm contrast is
### about the covariate set and nothing else.
###
### Each matchit() call is self-contained. 1:1 without replacement holds inside
### the call, and a control pixel drawn by one cohort stays available to the
### next, as it does across pairings. The verifier counts the reuse rate.
###
### THE MATCHABLE POOL
###
### Controls in exact cells that hold no treated unit for this cohort can never
### be matched and are dropped before the distance model is fitted. Exact for
### the match, and it changes the fitted propensity score, which is a Methods
### sentence rather than a silent choice. Without it every cohort fits its
### distance model on 14 to 24 million rows and bart and elasticnet do not run.
###
### RUNTIME
###
### Dominated by the per-cohort carbon reads in the two augmented arms, five
### years times three layers times the matchable pool, repeated per cohort.
### The baseline arm skips them entirely. Probe mode benchmarks a single fit at
### graduated input sizes for glm, elasticnet, and bart before the array is
### submitted.
###
### Usage
###   Rscript 8_match_pairings.R tasks     print the array size and exit
###   Rscript 8_match_pairings.R probe     sizes, cells, distances, benchmark
###   Rscript 8_match_pairings.R           run one array task
###   Rscript 8_match_pairings.R archive   copy completed output to /projects

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_extract.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr", "MatchIt"))

# MatchIt delegates two distances to other packages. A missing backend does not
# fail at load, it fails inside matchit() where try() swallows it, so the sweep
# would run to completion with two configs quietly absent. Checked here.
for (pkg in c("glmnet", "dbarts")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is not installed in this environment. ",
         "MatchIt needs glmnet for distance = 'elasticnet' and dbarts for ",
         "distance = 'bart'. Install it into MC0617 before running stage 8, ",
         "since bart is the distance Table S3 selected for the fire pairing.")
  }
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(MatchIt)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("tasks", "probe", "archive")) {
  args[1]
} else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

TASKS <- ca_match_tasks()

# CONFIG SUBSET.
#
# Each config writes its own matched file and appends to the two shared tables,
# so a task can run any subset with no merge step afterwards. This exists
# because bart is roughly 80 times slower than glm at the same input size, and
# a pairing whose largest level carries millions of matchable controls can put
# one task over the wall clock on the strength of two configs out of eleven.
# Splitting is a scheduling decision and changes no result.
#
# BEWARE sbatch --export. It splits its own argument on commas, so
#
#   sbatch --export=ALL,CA_CONFIGS=glm_nocal,glm_cal025 ...
#
# sets CA_CONFIGS to glm_nocal alone and then tries to export a variable named
# glm_cal025 from the calling environment. The value arrives silently truncated
# to its first element, which is exactly how a full sweep turns into a
# one-config run without a single error. Export the variable first instead.
#
#   export CA_CONFIGS="glm_nocal,glm_cal025"
#   sbatch --export=ALL --array=1-18 8_match_pairings.sub
#
# A colon separator also survives sbatch --export intact, so both work.
#
#   sbatch --export=ALL,CA_CONFIGS=glm_nocal:glm_cal025 --array=1-18 ...
#
CONFIGS <- ca_match_configs
sel <- Sys.getenv("CA_CONFIGS")
if (nzchar(sel)) {
  want <- trimws(strsplit(sel, "[,:;+[:space:]]+")[[1]])
  want <- want[nzchar(want)]
  have <- vapply(ca_match_configs, `[[`, character(1), "id")
  bad <- setdiff(want, have)
  if (length(bad)) {
    stop("Unknown config id: ", paste(bad, collapse = ", "),
         ". Known ids: ", paste(have, collapse = ", "))
  }
  CONFIGS <- ca_match_configs[match(want, have)]

  # A single config is legitimate for the two bart submissions and is also the
  # signature of the sbatch --export truncation above, so it is called out
  # rather than logged quietly.
  if (length(want) == 1L) {
    ca_log("NOTE CA_CONFIGS resolved to one config, '", want,
           "'. If a longer list was intended, sbatch --export split it on ",
           "commas. Export CA_CONFIGS in the shell first, or use colons.")
  }
}


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Read one side of a pairing across all three classes and stack them. Column
# sets differ between the stage 6 and stage 7 files, so the caller names what
# it needs and anything absent is filled with NA rather than silently dropped.
read_side <- function(role, groups, cols) {
  want <- unique(c("pixel_id", cols))
  parts <- lapply(ca_lulc$label, function(cl) {
    p <- ca_pool_path(cl, role)
    if (!file.exists(p)) stop("Missing ", role, " pool: ", p)
    ds <- arrow::open_dataset(p)
    have <- intersect(want, names(ds))
    if (!"pixel_id" %in% have) {
      stop("No pixel_id column in ", basename(p))
    }
    d <- as.data.frame(dplyr::collect(dplyr::select(
      dplyr::filter(ds, analysis_group %in% groups),
      dplyr::all_of(have))))
    # rep() rather than a scalar, because a class partition can legitimately
    # hold no rows for a group and scalar assignment fails on a zero-row frame.
    for (v in setdiff(want, have)) d[[v]] <- rep(NA, nrow(d))
    d$lulc <- rep(cl, nrow(d))
    d[, unique(c("pixel_id", "lulc", cols)), drop = FALSE]
  })
  d <- do.call(rbind, parts)
  # THE ROW ORDER MATTERS. Arrow's filtered scan is multithreaded and does not
  # guarantee a stable row order, so an unsorted frame makes sample() draw a
  # different treated set on every run and in every arm. The cap is meant to be
  # one seeded draw shared by all three arms and byte-identical on rerun, which
  # ca_durability requires of an archived artifact.
  d[order(d$pixel_id, method = "radix"), , drop = FALSE]
}

# Covariate join. Only the columns this arm needs, plus the exact set, and only
# the pixels the two sides actually hold. The full table is one row for each of
# 104 million grid pixels, so the filter is pushed into the read rather than
# applied afterwards.
read_covariates <- function(vars, ids) {
  need <- unique(c("pixel_id", ca_match$exact, vars))
  need <- setdiff(need, c("lulc", ca_pretrt_vars$covariate))
  ids <- as.integer(ids)
  parts <- lapply(ca_lulc$label, function(cl) {
    p <- ca_covariate_path(cl)
    if (!file.exists(p)) stop("Missing covariate table: ", p)
    ds <- arrow::open_dataset(p)
    as.data.frame(dplyr::collect(dplyr::select(
      dplyr::filter(ds, pixel_id %in% ids), dplyr::all_of(need))))
  })
  do.call(rbind, parts)
}

# The cap. Per analysis group, pooled across classes, seeded, treated side only.
apply_cap <- function(tr) {
  cap <- ca_sample$cap_per_group
  set.seed(ca_seed)
  keep <- unlist(lapply(split(seq_len(nrow(tr)), tr$analysis_group),
                        function(idx) {
                          if (length(idx) <= cap) return(idx)
                          sort(sample(idx, cap))
                        }), use.names = FALSE)
  tr[sort(keep), , drop = FALSE]
}

# Pre-treatment mean of one layer over a window, per pixel, no truncation.
#
# The window itself avoids the contaminated year (c-1 for thinning) through
# ca_pretrt_window(), so no in-window detection check is needed here. See the
# ca_pretrt$truncate_on_detection comment and the 8c diagnostic for why.
pretrt_mean <- function(covariate, pixel_id, lulc_vec, window) {
  i <- match(covariate, ca_pretrt_vars$covariate)
  if (is.na(i)) stop("Unknown pre-treatment covariate: ", covariate)
  src <- ca_pretrt_vars$source[i]
  lyr <- ca_pretrt_vars$layer[i]

  out <- rep(NA_real_, length(pixel_id))
  if (!length(window)) return(out)

  for (cl in unique(lulc_vec)) {
    k <- which(lulc_vec == cl)
    if (!length(k)) next
    m <- ca_read_years(src, lyr, window, cl, pixel_id = pixel_id[k])
    n_ok <- rowSums(!is.na(m))
    v <- rowMeans(m, na.rm = TRUE)
    v[n_ok < ca_pretrt$window_min] <- NA_real_
    out[k] <- v
    rm(m, n_ok, v); ca_gc()
  }
  out
}

# One matchit() call. Argument dispatch follows the MS3 template.
#
# ONE COMPLETE-CASE TIER, NOT TWO. The MS3 template gives the exact-only
# ablation its own tier, complete on the exact variables alone, so it is not
# penalised by missingness in covariates it never uses. That made sense there,
# where the building and weather covariates carried real missingness. Here the
# two tiers differ by about 0.03 percent of rows, so the second tier buys
# almost nothing and costs two things. summary.matchit refuses NA in
# addlvariables, which is how the difference surfaced. And the ablation stops
# being a controlled comparison, because it would run on a different sample
# from the ten configs it is meant to be compared against.
#
# So every config reads dat_nn. The exact-only ablation keeps its treat ~ 1
# formula and reports balance on the continuous covariates through
# addlvariables, which is the comparison it exists to provide.
fit_one_config <- function(cfg, dat_nn, xvars, exact_form) {

  is_exact_only <- cfg$id == "nn_exact_1to1"
  dat <- dat_nn
  xform <- as.formula(paste("treat ~", paste(xvars, collapse = " + ")))
  mform <- as.formula(paste("~", paste(xvars, collapse = " + ")))

  a <- list(data = dat, method = "nearest", exact = exact_form,
            ratio = ca_match$ratio, replace = ca_match$replace,
            estimand = "ATT")

  if (is_exact_only) {
    a$formula  <- treat ~ 1
    a$distance <- "glm"
    a$link     <- "linear.logit"
  } else if (cfg$distance %in% c("mahalanobis", "robust_mahalanobis")) {
    a$formula  <- xform
    a$distance <- cfg$distance
    if (cfg$caliper) {
      a$caliper     <- setNames(rep(ca_match$caliper_sd, length(xvars)), xvars)
      a$std.caliper <- TRUE
    }
  } else if (cfg$id == "mahal_pscal") {
    a$formula     <- xform
    a$distance    <- "glm"
    a$link        <- "linear.logit"
    a$mahvars     <- mform
    a$caliper     <- ca_match$caliper_sd
    a$std.caliper <- TRUE
  } else {
    a$formula  <- xform
    a$distance <- cfg$distance
    if (cfg$distance == "glm") a$link <- "linear.logit"
    if (isTRUE(cfg$mahvars)) a$mahvars <- mform
    if (cfg$caliper) {
      a$caliper     <- ca_match$caliper_sd
      a$std.caliper <- TRUE
    }
  }

  try(do.call(matchit, a), silent = TRUE)
}

# BALANCE COMPUTED DIRECTLY, NOT THROUGH summary.matchit.
#
# Two failures forced this and both are properties of MatchIt's reporting path
# rather than of the match. summary.matchit refuses covariates carrying NA in
# addlvariables, which killed the first baseline run. It also builds a model
# frame over the exact variables, which dies with "contrasts can be applied
# only to factors with 2 or more levels" the moment a stratification level
# contains a single lulc, ecoregion, or HUC8. Level 1985_1988 of UP_UNP has two
# exact cells and 18 treated units, so that is not an edge case, it is a third
# of the levels in the thin pairings.
#
# Definitions follow MatchIt so the numbers stay comparable with Table S3 as
# submitted. SMD standardises on the treated standard deviation in the
# unadjusted sample and reuses that same denominator after matching, so the two
# rows are comparable. Var. Ratio is treated over control. eCDF Mean and Max are
# the mean and maximum absolute difference between the two empirical
# distribution functions.
#
# The eCDF pair is evaluated on a grid of every treated value plus up to
# ecdf_grid_max quantiles of the control distribution. The CDFs themselves are
# exact, only the evaluation grid is thinned, and it is thinned only where the
# control side runs to millions of rows.
ecdf_grid_max <- 100000L

balance_one <- function(x_t, x_c, denom) {
  x_t <- x_t[is.finite(x_t)]
  x_c <- x_c[is.finite(x_c)]
  if (!length(x_t) || !length(x_c)) {
    return(c(n_treat = length(x_t), n_ctrl = length(x_c),
             mean_treat = NA_real_, mean_ctrl = NA_real_,
             sd_treat = NA_real_, sd_ctrl = NA_real_,
             smd = NA_real_, vr = NA_real_,
             ecdf_mean = NA_real_, ecdf_max = NA_real_))
  }
  mt <- mean(x_t); mc <- mean(x_c)
  vt <- stats::var(x_t); vc <- stats::var(x_c)

  st <- sort(x_t); sc <- sort(x_c)
  grid <- if (length(sc) > ecdf_grid_max) {
    stats::quantile(sc, probs = seq(0, 1, length.out = ecdf_grid_max),
                    names = FALSE, type = 1L)
  } else sc
  grid <- sort(unique(c(st, grid)))
  dd <- abs(findInterval(grid, st) / length(st) -
            findInterval(grid, sc) / length(sc))

  # n and sd travel with every row so stage 8b can pool means and variances
  # across cohorts exactly, rather than averaging summary statistics. The
  # pooled variance needs the within-level variance and the between-level
  # spread of the means, and neither is recoverable from smd and vr alone.
  c(n_treat    = length(x_t),
    n_ctrl     = length(x_c),
    mean_treat = mt,
    mean_ctrl  = mc,
    sd_treat   = sqrt(vt),
    sd_ctrl    = sqrt(vc),
    smd        = if (is.finite(denom) && denom > 0) (mt - mc) / denom
                 else NA_real_,
    vr         = if (is.finite(vc) && vc > 0) vt / vc else NA_real_,
    ecdf_mean  = mean(dd),
    ecdf_max   = max(dd))
}

# One row per covariate for one sample. denoms is the treated standard
# deviation from the unadjusted sample, passed in so the adjusted row shares it.
balance_table <- function(d, xvars, stage, denoms) {
  tr <- d$treat == 1L
  out <- lapply(xvars, function(v) {
    s <- balance_one(d[[v]][tr], d[[v]][!tr], denoms[[v]])
    data.frame(covariate = v, stage = stage, as.list(s),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# Match counts, read off the matchit object rather than its summary.
match_counts <- function(mobj) {
  tr <- mobj$treat == 1L
  w <- mobj$weights
  c(dropped_treated = sum(tr & w == 0),
    dropped_control = sum(!tr & w == 0))
}

# File-system lock and stream append, carried from the MS3 template. Eighteen
# tasks write two shared tables and the append has to be atomic.
acquire_lock <- function(lock_dir, timeout = 900, stale = 1800, interval = 0.5) {
  t0 <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) return(invisible(TRUE))
    mtime <- file.info(lock_dir)$mtime
    if (!is.na(mtime) &&
        difftime(Sys.time(), mtime, units = "secs") > stale) {
      warning("Removing stale lock: ", lock_dir)
      unlink(lock_dir, recursive = TRUE, force = TRUE)
      next
    }
    if (difftime(Sys.time(), t0, units = "secs") > timeout) {
      stop("Timed out acquiring lock: ", lock_dir)
    }
    Sys.sleep(interval + runif(1, 0, 0.3))
  }
}

stream_append_gz <- function(d, path, lock_dir) {
  if (is.null(d) || !nrow(d)) return(invisible(NULL))
  acquire_lock(lock_dir)
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  hdr <- !file.exists(path)
  con <- gzfile(path, open = "at")
  on.exit(close(con), add = TRUE, after = FALSE)
  utils::write.table(d, con, sep = ",", row.names = FALSE,
                     col.names = hdr, qmethod = "double")
  invisible(NULL)
}

# Nearest-neighbour distance among control pixels, on the terra cell index. No
# spatial package, and the answer decides whether the control pool needs
# spacing. Sampled, because the exact all-pairs answer on 14 million rows is
# not needed to read a distribution.
control_nn_distance <- function(pixel_id, n_sample = 20000L) {
  ncol_grid <- ca_grid_ncol()
  z <- as.double(pixel_id) - 1
  row <- z %/% ncol_grid
  col <- z %% ncol_grid
  set.seed(ca_seed)
  idx <- if (length(pixel_id) > n_sample) sample.int(length(pixel_id), n_sample)
         else seq_along(pixel_id)
  # Blocked scan. A control's nearest neighbour is almost always within a few
  # rows, so the search is bounded rather than all-pairs.
  ord <- order(row, col, method = "radix")
  r <- row[ord]; c0 <- col[ord]
  pos <- match(idx, ord)
  d <- vapply(pos, function(p) {
    lo <- max(1L, p - 200L); hi <- min(length(r), p + 200L)
    k <- setdiff(lo:hi, p)
    if (!length(k)) return(NA_real_)
    sqrt(min((r[k] - r[p])^2 + (c0[k] - c0[p])^2))
  }, numeric(1))
  d * 30   # pixel units to metres
}


# ---------------------------------------------------------------------------
# TASKS MODE
# ---------------------------------------------------------------------------

if (mode == "tasks") {
  ca_log("Stage 8 array size: ", nrow(TASKS))
  print(TASKS)
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_match_dir(), "matched", subdir = "work")
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# INPUT ASSEMBLY, SHARED BY PROBE AND RUN
# ---------------------------------------------------------------------------

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(TASKS)) {
  stop("Task id must be 1 to ", nrow(TASKS), ", got ", task)
}

PAIRING <- TASKS$pairing[task]
ARM     <- TASKS$arm[task]
XVARS   <- ca_match_covariates(ARM)
NEEDS_PRETRT <- ca_match_arm_needs_pretrt(ARM)
COHORT_VAR   <- ca_match_cohort_var(PAIRING)

ca_log("Stage 8, task ", task, ".  pairing ", PAIRING, "  arm ", ARM)
ca_log("  continuous covariates (", length(XVARS), "): ",
       paste(XVARS, collapse = ", "))
ca_log("  exact: ", paste(ca_match$exact, collapse = ", "))
ca_log("  cohort variable: ", COHORT_VAR)
ca_log("  configs this task (", length(CONFIGS), "): ",
       paste(vapply(CONFIGS, `[[`, character(1), "id"), collapse = ", "))

tr_groups <- ca_pairing_groups(PAIRING, "treat")
ct_groups <- ca_pairing_groups(PAIRING, "control")

treated <- read_side("treatment", tr_groups,
                     c("analysis_group", "stratum", "cohort_year",
                       "pa_est_year"))
ca_log("  treated, spaced: ", format(nrow(treated), big.mark = ","))

control <- read_side("control", ct_groups,
                     c("analysis_group", "eligibility", "pa_est_year"))
ca_log("  control pool, whole: ", format(nrow(control), big.mark = ","))

# Stage 6 marks a pixel eligible when it survives every exclusion. A control
# pool row that is not eligible is a coding error upstream rather than a
# routine drop, so it is counted rather than filtered silently.
if ("eligibility" %in% names(control)) {
  bad <- sum(!is.na(control$eligibility) &
               control$eligibility == "unassigned")
  if (bad) ca_log("  NOTE control rows flagged unassigned: ",
                  format(bad, big.mark = ","))
}

treated <- apply_cap(treated)
ca_log("  treated after the cap, per analysis group pooled across classes: ",
       format(nrow(treated), big.mark = ","))
print(table(treated$analysis_group))

static_vars <- setdiff(XVARS, ca_pretrt_vars$covariate)
covs <- read_covariates(XVARS, c(treated$pixel_id, control$pixel_id))
ca_log("  covariate rows joined: ", format(nrow(covs), big.mark = ","))

# pixel_id is the terra cell index of one frozen grid, so it is unique across
# the three class partitions and the join needs no composite key.
join_covs <- function(d) {
  k <- match(d$pixel_id, covs$pixel_id)
  if (anyNA(k)) {
    stop(sum(is.na(k)), " pixels absent from the stage 4b covariate table")
  }
  for (v in c(ca_match$exact, static_vars)) {
    if (identical(v, "lulc")) next
    d[[v]] <- covs[[v]][k]
  }
  d
}
treated <- join_covs(treated)
control <- join_covs(control)
rm(covs); ca_gc("after covariate join")

treated$cell <- paste(treated$lulc, treated$ecoregion_l3, treated$huc8,
                      sep = "|")
control$cell <- paste(control$lulc, control$ecoregion_l3, control$huc8,
                      sep = "|")

# Stratification level and the window that defines it.
lev_year <- treated[[COHORT_VAR]]
win_list <- lapply(lev_year, function(y) ca_pretrt_window(PAIRING, y))
win_key  <- vapply(win_list, function(w) {
  if (!length(w)) NA_character_ else sprintf("%d_%d", min(w), max(w))
}, character(1))

drop_nowin <- is.na(win_key)
if (any(drop_nowin)) {
  ca_log("  treated dropped for an unusable pre-treatment window: ",
         format(sum(drop_nowin), big.mark = ","))
  treated <- treated[!drop_nowin, , drop = FALSE]
  lev_year <- lev_year[!drop_nowin]
  win_key <- win_key[!drop_nowin]
}
treated$win_key <- win_key
treated$lev_year <- lev_year

levels_tbl <- unique(data.frame(win_key = treated$win_key,
                                stringsAsFactors = FALSE))
levels_tbl <- levels_tbl[order(levels_tbl$win_key), , drop = FALSE]
rownames(levels_tbl) <- NULL
ca_log("  stratification levels: ", nrow(levels_tbl))

# LEVEL SUBSET.
#
# Companion to CA_CONFIGS, for the pairings where one config on one task still
# exceeds the wall clock. Levels are independent of each other, so a subset run
# is exact. The level order is sorted rather than encounter order, so the same
# range always means the same levels.
#
#   CA_LEVELS=1-18    the first eighteen levels
#
# When a subset is active the matched filename carries the range, so two halves
# of the same config cannot overwrite each other. Stage 9 reads them by glob.
LEVEL_TAG <- ""
lev_sel <- Sys.getenv("CA_LEVELS")
if (nzchar(lev_sel)) {
  p <- as.integer(strsplit(lev_sel, "-", fixed = TRUE)[[1]])
  if (length(p) != 2L || anyNA(p) || p[1] < 1L || p[2] < p[1]) {
    stop("CA_LEVELS must look like 1-18, got: ", lev_sel)
  }
  p[2] <- min(p[2], nrow(levels_tbl))
  if (p[1] > nrow(levels_tbl)) {
    ca_log("  CA_LEVELS starts past the last level, nothing to do")
    quit(save = "no")
  }
  levels_tbl <- levels_tbl[p[1]:p[2], , drop = FALSE]
  LEVEL_TAG <- sprintf("_lev%d-%d", p[1], p[2])
  ca_log("  running levels ", p[1], " to ", p[2], ", ",
         nrow(levels_tbl), " of them")
}


# ---------------------------------------------------------------------------
# PROBE MODE
# ---------------------------------------------------------------------------

if (mode == "probe") {

  ca_log("Probe, pairing ", PAIRING, ", arm ", ARM)

  cells_tr <- unique(treated$cell)
  matchable <- control[control$cell %in% cells_tr, , drop = FALSE]

  per_lev <- do.call(rbind, lapply(levels_tbl$win_key, function(w) {
    tr <- treated[treated$win_key == w, , drop = FALSE]
    ct <- control[control$cell %in% unique(tr$cell), , drop = FALSE]
    data.frame(win_key = w, n_treat = nrow(tr),
               n_cells = length(unique(tr$cell)),
               n_control_matchable = nrow(ct),
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(per_lev,
                   ca_meta(sprintf("match_probe_levels_%s_%s.csv",
                                   PAIRING, ARM)), row.names = FALSE)
  print(utils::head(per_lev[order(-per_lev$n_treat), ], 20))

  nnd <- control_nn_distance(control$pixel_id)
  qs <- stats::quantile(nnd, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
  ca_log("  control nearest-neighbour distance, metres:")
  print(round(qs, 1))
  ca_log("  share of controls with a neighbour closer than 150 m: ",
         round(100 * mean(nnd < 150, na.rm = TRUE), 2), " percent")

  # Timed single fit at graduated sizes. Skippable, since the benchmark is the
  # slow part of the probe and the numbers do not move between pairings.
  #   CA_SKIP_BENCH=1
  bench <- NULL
  if (length(static_vars) >= 2 && !nzchar(Sys.getenv("CA_SKIP_BENCH"))) {
    for (n in c(5e4, 2e5, 1e6)) {
      if (nrow(matchable) < n) next
      set.seed(ca_seed)
      bcols <- unique(c(ca_match$exact, static_vars))
      a <- treated[sample.int(nrow(treated), min(nrow(treated), 5000L)),
                   bcols, drop = FALSE]
      b <- matchable[sample.int(nrow(matchable), n), bcols, drop = FALSE]
      a$treat <- 1L
      b$treat <- 0L
      d <- rbind(a, b)
      rm(a, b)
      d <- d[stats::complete.cases(d[, c(static_vars, ca_match$exact)]), ]
      for (v in ca_match$exact) d[[v]] <- factor(d[[v]])
      ef <- as.formula(paste("~", paste(ca_match$exact, collapse = " + ")))
      for (dist in c("glm", "elasticnet", "bart")) {
        t0 <- Sys.time()
        fit <- try(matchit(
          as.formula(paste("treat ~", paste(static_vars, collapse = " + "))),
          data = d, method = "nearest", exact = ef, distance = dist,
          ratio = 1L, replace = FALSE, estimand = "ATT"), silent = TRUE)
        el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        bench <- rbind(bench, data.frame(
          pairing = PAIRING, n_control = n, distance = dist,
          seconds = round(el, 1),
          ok = !inherits(fit, "try-error"),
          message = if (inherits(fit, "try-error"))
            gsub("[\r\n]+", " ", trimws(as.character(fit))) else "",
          stringsAsFactors = FALSE))
        rm(fit); ca_gc()
      }
      rm(d); ca_gc()
    }
  }
  if (!is.null(bench)) {
    print(bench)
    utils::write.csv(bench, ca_meta(sprintf("match_probe_bench_%s_%s.csv",
                                            PAIRING, ARM)),
                     row.names = FALSE)
  }

  # The augmented arms read five years of three carbon layers per level, over
  # the matchable control pool. That read, not the matchit() fit, is what sets
  # the wall clock, and nothing above measures it. One timed call on the
  # largest level answers it.
  if (NEEDS_PRETRT) {
    big <- per_lev$win_key[which.max(per_lev$n_control_matchable)]
    tr_b <- treated[treated$win_key == big, , drop = FALSE]
    ct_b <- control[control$cell %in% unique(tr_b$cell), , drop = FALSE]
    window <- {
      p <- as.integer(strsplit(big, "_", fixed = TRUE)[[1]])
      seq(p[1], p[2])
    }
    ca_log("  timing one pre-treatment read, level ", big, ", ",
           format(nrow(ct_b), big.mark = ","), " controls, ",
           length(window), " years")
    v <- intersect(XVARS, ca_pretrt_vars$covariate)[1]
    t0 <- Sys.time()
    z <- pretrt_mean(v, ct_b$pixel_id, ct_b$lulc, window)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ca_log("  ", v, " on the largest level took ", round(el, 1), " seconds, ",
           round(100 * mean(!is.na(z)), 2), " percent non-missing")
    ca_log("  projected for ", length(intersect(XVARS,
           ca_pretrt_vars$covariate)), " covariates across ",
           nrow(per_lev), " levels, order ",
           round(el * length(intersect(XVARS, ca_pretrt_vars$covariate)) *
                 nrow(per_lev) / 3600, 1), " hours if every level were this ",
           "large, which is an upper bound")
    rm(tr_b, ct_b, z); ca_gc()
  }

  ca_stamp(sprintf("8_probe_%s_%s", PAIRING, ARM))
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# RUN
# ---------------------------------------------------------------------------


exact_form <- as.formula(paste("~", paste(ca_match$exact, collapse = " + ")))

matched_acc <- setNames(vector("list", length(CONFIGS)),
                        vapply(CONFIGS, `[[`, character(1), "id"))
specs_acc   <- list()
balance_acc <- list()

for (w in levels_tbl$win_key) {

  tr <- treated[treated$win_key == w, , drop = FALSE]
  window <- {
    p <- as.integer(strsplit(w, "_", fixed = TRUE)[[1]])
    seq(p[1], p[2])
  }
  cohort_year <- max(window) + 1L

  ct <- control[control$cell %in% unique(tr$cell), , drop = FALSE]

  # FP_UP and TP_UP. A UP control must have been protected before the treated
  # unit's cohort year, the symmetric half of ca_pa_rules$label_at_cohort_year.
  # Unknown establishment is CPAD YR_EST 0, assigned the floor as elsewhere.
  if (identical(ct_groups, "UP")) {
    est <- ifelse(is.na(ct$pa_est_year), ca_pretrt$cpad_year_floor,
                  ct$pa_est_year)
    n0 <- nrow(ct)
    ct <- ct[est <= cohort_year, , drop = FALSE]
    ca_log("  [", w, "] UP controls dropped for late protection: ",
           format(n0 - nrow(ct), big.mark = ","))
  }

  ca_log("  [", w, "] treated ", format(nrow(tr), big.mark = ","),
         "  matchable controls ", format(nrow(ct), big.mark = ","),
         "  cells ", length(unique(tr$cell)))

  if (!nrow(tr) || !nrow(ct)) {
    ca_log("    no matchable pair in this level, skipping")
    next
  }

  if (NEEDS_PRETRT) {
    for (v in intersect(XVARS, ca_pretrt_vars$covariate)) {
      tr[[v]] <- pretrt_mean(v, tr$pixel_id, tr$lulc, window)
      ct[[v]] <- pretrt_mean(v, ct$pixel_id, ct$lulc, window)
    }
    ca_gc(paste0("after pre-treatment reads, level ", w))
  }

  keep_cols <- unique(c("pixel_id", "lulc", "analysis_group", "stratum",
                        ca_match$exact, XVARS))
  tr$stratum <- if ("stratum" %in% names(tr)) tr$stratum else NA_character_
  ct$stratum <- NA_character_
  tr2 <- tr[, intersect(keep_cols, names(tr)), drop = FALSE]
  ct2 <- ct[, intersect(keep_cols, names(ct)), drop = FALSE]
  tr2$treat <- 1L
  ct2$treat <- 0L
  input <- rbind(tr2, ct2)
  input$cohort_year <- cohort_year
  input$win_key <- w
  rm(tr2, ct2)

  for (v in ca_match$exact) input[[v]] <- factor(input[[v]])

  dat_ex <- input[stats::complete.cases(
    input[, c("treat", ca_match$exact)]), , drop = FALSE]
  dat_nn <- input[stats::complete.cases(
    input[, c("treat", ca_match$exact, XVARS)]), , drop = FALSE]

  ca_log("    complete cases, exact only ",
         format(nrow(dat_ex), big.mark = ","),
         "  exact plus continuous ", format(nrow(dat_nn), big.mark = ","),
         "  loss to the continuous covariates ",
         format(nrow(dat_ex) - nrow(dat_nn), big.mark = ","),
         " rows, ",
         round(100 * (nrow(dat_ex) - nrow(dat_nn)) / max(nrow(dat_ex), 1), 3),
         " percent")

  # The matching input is the peak, not the covariate join that gets logged
  # earlier. Cores on acpu are bought for their RAM at 3.8 GB each, so this
  # line is what tells us whether 32 cores is right or wasteful.
  ca_gc(paste0("matching input built, level ", w))

  # Unadjusted balance is a property of the level, not of a config, so it is
  # computed once here rather than eleven times. The treated standard deviation
  # from this sample becomes the SMD denominator for every adjusted row, which
  # is what makes the before and after numbers comparable.
  denoms <- lapply(XVARS, function(v) {
    z <- dat_nn[[v]][dat_nn$treat == 1L]
    stats::sd(z[is.finite(z)])
  })
  names(denoms) <- XVARS
  bal_un <- balance_table(dat_nn, XVARS, "unadjusted", denoms)

  for (cfg in CONFIGS) {

    dat <- dat_nn
    n_tr <- sum(dat$treat == 1L)
    n_ct <- sum(dat$treat == 0L)

    spec <- data.frame(
      pairing = PAIRING, arm = ARM, config = cfg$id, win_key = w,
      cohort_year = cohort_year,
      n_input = nrow(dat), n_treat = n_tr, n_control = n_ct,
      n_matched = NA_integer_, match_rate = NA_real_,
      dropped_treated = NA_real_, dropped_control = NA_real_,
      status = "ok", stringsAsFactors = FALSE)

    if (n_tr == 0L || n_ct == 0L) {
      spec$status <- "no units"
      specs_acc[[length(specs_acc) + 1L]] <- spec
      next
    }

    mobj <- fit_one_config(cfg, dat_nn, XVARS, exact_form)
    if (inherits(mobj, "try-error")) {
      spec$status <- paste("failed:", trimws(as.character(mobj)))
      specs_acc[[length(specs_acc) + 1L]] <- spec
      warning("matchit failed [", PAIRING, "/", ARM, "/", cfg$id, "/", w, "]")
      next
    }

    md <- match.data(mobj)
    # match.data() returns a matchdata object, and MatchIt defines an
    # rbind.matchdata method that fails when the parts are stacked with
    # do.call(). The class is dropped here so the accumulated parts stack as
    # ordinary data frames. Its renumbering behaviour is not wanted either,
    # since subclass is namespaced by level immediately below.
    class(md) <- "data.frame"
    md$subclass <- paste(w, md$subclass, sep = ":")
    md$pairing <- PAIRING
    md$arm <- ARM
    md$config <- cfg$id
    matched_acc[[cfg$id]][[w]] <- md

    b <- rbind(bal_un, balance_table(md, XVARS, "adjusted", denoms))
    b$pairing <- PAIRING; b$arm <- ARM; b$config <- cfg$id
    b$win_key <- w; b$cohort_year <- cohort_year
    balance_acc[[length(balance_acc) + 1L]] <- b

    cnt <- match_counts(mobj)

    spec$n_matched <- nrow(md)
    spec$match_rate <- (nrow(md) / 2) / n_tr
    spec$dropped_treated <- cnt[["dropped_treated"]]
    spec$dropped_control <- cnt[["dropped_control"]]
    specs_acc[[length(specs_acc) + 1L]] <- spec

    rm(mobj, md, b, cnt); ca_gc(paste0("after ", cfg$id, ", level ", w))
  }

  rm(tr, ct, input, dat_ex, dat_nn, bal_un, denoms)
  ca_gc(paste0("level ", w, " done"))
}


# ---------------------------------------------------------------------------
# WRITE
# ---------------------------------------------------------------------------

for (id in names(matched_acc)) {
  parts <- matched_acc[[id]]
  if (!length(parts)) next
  d <- do.call(rbind, parts)
  p <- ca_match_path(PAIRING, ARM, paste0(id, LEVEL_TAG))
  if (file.exists(p) && !OVERWRITE) {
    ca_log("Present, not overwritten: ", basename(p))
    next
  }
  arrow::write_parquet(arrow::as_arrow_table(d), p,
                       compression = ca_io$compression)
  ca_log("Wrote ", basename(p), "  ", format(nrow(d), big.mark = ","),
         " rows  ", round(file.size(p) / 1e6, 1), " MB")
  rm(d); ca_gc()
}

specs <- if (length(specs_acc)) do.call(rbind, specs_acc) else NULL
bal   <- if (length(balance_acc)) do.call(rbind, balance_acc) else NULL

stream_append_gz(specs, ca_match_specs_path(),
                 file.path(ca_match_dir(), ".lock_specs"))
stream_append_gz(bal, ca_match_balance_path(),
                 file.path(ca_match_dir(), ".lock_balance"))

if (!is.null(specs)) {
  ok <- specs[specs$status == "ok", , drop = FALSE]
  if (nrow(ok)) {
    agg <- stats::aggregate(cbind(n_matched, n_treat) ~ config, data = ok,
                            FUN = sum)
    agg$match_rate <- round((agg$n_matched / 2) / agg$n_treat, 4)
    ca_log("Pooled across levels, ", PAIRING, " ", ARM, ":")
    print(agg[order(-agg$n_matched), ])
  }
}

ca_stamp(sprintf("8_match_%s_%s", PAIRING, ARM))
ca_log("Stage 8 task ", task, " complete.  pairing ", PAIRING, "  arm ", ARM)
