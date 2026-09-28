### CA carbon revision pipeline
### Stage 10c. Completeness check and results aggregation
###
### WHAT THIS STAGE DOES
###
### Reads every file 10a and 10b wrote, binds them into one tidy table, and
### records what is not there. It is the boundary between the cluster and the
### laptop. Everything upstream needs Alpine. Everything downstream reads
### did_results.csv.gz and needs nothing else.
###
### It runs on partial output by design. 10b is a 1,245 task array and the
### large tier can take days, so the aggregate is regenerated as tasks land
### rather than waited for. did_results_coverage.csv records completeness at
### the moment of writing, so a figure built from an incomplete aggregate is
### identifiable after the fact rather than mistaken for a finished one.
###
### WHAT IT DOES NOT DO
###
### No filtering on ecoregion, no filtering on outcome, no selection of an
### estimator. Every cell that ran is in the output. The six forested
### ecoregions and the etwfe-only figure rule are stage 11 decisions and live
### in ca_fig, not here. If a reviewer asks for a cell the figures do not show,
### it is already in this file.
###
### ONE CONNECTION PER CALL
###
### gz_has_rows() opens, reads, and closes inside a single function call, so
### its on.exit fires once and against the connection it was registered for.
### The first version of this script registered on.exit inside a loop inside a
### function, which defers every close to the function's return and evaluates
### them all against the last connection. That produced a wall of "invalid
### connection" errors.
###
### The same block at top level in 10b_did_etwfe.R does not error, because
### on.exit has no function context there and is ignored. Those connections
### leak and the garbage collector closes them, which is where the "closing
### unused connection" warnings in the 10b missing mode come from. Harmless,
### and not worth editing a script whose array is still running.
###
### THE COLUMN THE CAPS COME FROM
###
### se_ratio_prev is the standard error divided by the standard error at the
### preceding time point within the same run and aggregation. It answers two
### questions with one column. A ratio above two marks where the contributing
### cohort set has thinned enough that the estimate describes a different
### population, which is the empirical basis for the event-time caps. A ratio
### below one half marks the terminal-year variance collapse seen in UP_UNP
### GPP at 2025 and NEP at 2021, which has to be understood before it appears
### in a figure.
###
### THE TWO MEANINGS OF arm
###
### ca_pairings$arm is the treatment arm, pa or offset or fire or thin. The
### manifest arm column is ca_match_selected()$arm, the matching arm, baseline
### or augmented. They are different things with the same name, so the manifest
### one is renamed match_arm on read and arm is the treatment arm throughout.
###
### Usage
###   Rscript 10c_aggregate_did.R            completeness check, writes only
###                                          did_results_missing.csv
###   Rscript 10c_aggregate_did.R aggregate  bind and write the results table
###   Rscript 10c_aggregate_did.R report     print the last coverage table
###
### Delete the three outputs before rerunning rather than relying on overwrite.
###
###   rm /projects/mich9173/CA_carbon/meta/did_results*.csv*

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("check", "aggregate", "report")) {
  args[1]
} else "check"

ESTIMATORS <- c("attgt", "etwfe")

# ONE ROW PER CELL, BUT NOT ONE TASK PER CELL.
#
# ca_did_runs() carries every cell for both estimators and a task column saying
# which array task computes it. For etwfe that is one to one. For att_gt a task
# is a statewide cell that also computes its ecoregions in a loop, so several
# rows share a task and the array strings below are the unique task ids rather
# than the row numbers. array_string() already sorts and deduplicates.
JOB_NAME   <- c(attgt = "ca_10a_attgt", etwfe = "ca_10b_etwfe")

dir.create(ca_meta(), recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------
# STATUS
# ---------------------------------------------------------------------------
# WHAT COUNTS AS DONE. The same definition the array uses, because a file this
# script accepts is a file 10b will never recompute. A present file is not
# enough. A task killed mid-write leaves a truncated csv.gz, and both run paths
# skip any cell whose output already exists. Done means the gzip opens and
# carries a header plus at least one row.
#
# WHAT IS NOT MISSING. A cell with no viable cohort writes a manifest row with
# status skipped and no output, by design. It must not sit in the retry array
# forever.
#
# WHAT IS STILL RUNNING. squeue is read where it exists, so a task in flight is
# not handed back for resubmission on another node. This is the case that
# matters while an array drains.

# One connection, opened and closed inside one call. on.exit is scoped to this
# invocation, so it fires once and against the connection it was registered
# for.
gz_has_rows <- function(path) {
  con <- try(gzfile(path, "r"), silent = TRUE)
  if (inherits(con, "try-error")) return(FALSE)
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  ln <- suppressWarnings(try(readLines(con, n = 2L), silent = TRUE))
  !inherits(ln, "try-error") && length(ln) >= 2L
}

valid_output <- function(paths) {
  present <- file.exists(paths)
  ok <- rep(FALSE, length(paths))
  for (i in which(present)) {
    sz <- file.info(paths[i])$size
    if (is.na(sz) || sz < 50) next
    ok[i] <- gz_has_rows(paths[i])
  }
  list(present = present, valid = ok)
}

queued_tasks <- function(job_name) {
  if (!nzchar(Sys.which("squeue"))) return(integer(0))
  txt <- tryCatch(
    system2("squeue", c("-h", "-u", Sys.getenv("USER"), "-n", job_name,
                        "-o", "%K"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0))
  v <- suppressWarnings(as.integer(trimws(txt)))
  v[!is.na(v)]
}

read_manifest <- function(estimator) {
  p <- ca_did_est_manifest_path(estimator)
  if (!file.exists(p)) return(NULL)
  # Character in, character out. An ecoregion of "08" returns as the number 8
  # under default colClasses and an empty stratum returns as NA, which is the
  # corruption 10a and 10b guard against on every append.
  utils::read.csv(p, colClasses = "character", stringsAsFactors = FALSE)
}

# One row per expected run, for one estimator.
status_table <- function(estimator) {

  runs <- ca_did_runs(estimator)
  ca_log(estimator, ": ", nrow(runs), " expected runs, checking output")

  paths <- vapply(seq_len(nrow(runs)), function(i) {
    ca_did_est_path(estimator, runs$pairing[i], runs$run_id[i])
  }, character(1))

  v   <- valid_output(paths)
  man <- read_manifest(estimator)

  st <- if (is.null(man)) rep(NA_character_, nrow(runs)) else {
    man$status[match(runs$run_id, man$run_id)]
  }
  msg <- if (is.null(man)) rep("", nrow(runs)) else {
    m <- man$message[match(runs$run_id, man$run_id)]
    m[is.na(m)] <- ""
    m
  }

  # A RESOURCE ERROR IS NOT AN ESTIMATOR ERROR. The fit is wrapped in tryCatch,
  # so a task that ran out of address space inside feols writes a manifest row
  # with status error and exits zero. Those rows belong in the retry array on a
  # bigger node, not in a list of things no node can fix.
  res <- grepl("cannot allocate|out of memory|memory exhausted|std::bad_alloc",
               msg, ignore.case = TRUE)

  running <- queued_tasks(JOB_NAME[[estimator]])

  cls <- rep("todo", nrow(runs))
  cls[!v$valid & !is.na(st) & st == "skipped"]      <- "skipped"
  cls[!v$valid & !is.na(st) & st == "error" & !res] <- "error"
  cls[runs$task %in% running & !v$valid]            <- "running"
  cls[v$valid]                                      <- "done"

  data.frame(
    estimator = estimator,
    task      = runs$task,
    run_id    = runs$run_id,
    pairing   = runs$pairing,
    arm       = ca_pairings$arm[match(runs$pairing, ca_pairings$pairing)],
    outcome   = runs$outcome,
    transform = runs$transform,
    stratum   = runs$stratum,
    ecoregion = runs$ecoregion,
    tier      = runs$tier,
    n_treated = runs$n_treated,
    n_cohorts = runs$n_cohorts,
    est_cost  = runs$est_cost,
    class     = cls,
    mem_error = res,
    truncated = v$present & !v$valid,
    path      = paths,
    message   = msg,
    stringsAsFactors = FALSE)
}

STATUS <- do.call(rbind, lapply(ESTIMATORS, status_table))
rownames(STATUS) <- NULL


# ---------------------------------------------------------------------------
# REPORT THE STATUS
# ---------------------------------------------------------------------------

# Contiguous runs collapsed to an sbatch array string.
array_string <- function(v) {
  if (!length(v)) return("")
  v <- sort(unique(v))
  g <- cumsum(c(1L, as.integer(diff(v) != 1L)))
  paste(vapply(split(v, g), function(x) {
    if (length(x) == 1L) as.character(x) else paste0(min(x), "-", max(x))
  }, character(1)), collapse = ",")
}

print_status <- function() {
  for (es in ESTIMATORS) {
    s <- STATUS[STATUS$estimator == es, , drop = FALSE]
    cat("\n", es, ", ", nrow(s), " expected runs\n", sep = "")
    print(table(s$class, s$tier))

    tr <- which(s$truncated)
    if (length(tr)) {
      cat("\nTRUNCATED OUTPUT, delete these before resubmitting:\n")
      cat(paste0("  ", s$path[tr], collapse = "\n"), "\n")
    }

    e <- which(s$class == "error")
    if (length(e)) {
      cat("\nESTIMATOR ERRORS, not a resource problem:\n")
      for (i in e) cat("  ", s$run_id[i], "  ", s$message[i], "\n", sep = "")
    }

    todo <- which(s$class == "todo")
    if (length(todo)) {
      if (any(s$mem_error[todo])) {
        cat("\n", sum(s$mem_error[todo]),
            " runs failed on memory inside the fit and are in the retry ",
            "array.\n", sep = "")
      }
      cat("\nTO SUBMIT, one array string per tier:\n")
      for (ti in c("small", "medium", "large")) {
        v <- s$task[todo][s$tier[todo] == ti]
        if (!length(v)) next
        cat(sprintf("  %-7s %d runs   --array=%s\n",
                    ti, length(v), array_string(v)))
      }
    } else {
      cat("\nNothing to submit.\n")
    }
  }
  cat("\n")
}

write_status <- function() {
  keep <- STATUS[STATUS$class %in% c("todo", "error", "running"), ,
                 drop = FALSE]
  utils::write.csv(keep, ca_did_results_path("missing"), row.names = FALSE)
  ca_log("Wrote ", ca_did_results_path("missing"), "  ", nrow(keep), " rows")
  invisible(keep)
}

if (mode == "check") {
  print_status()
  write_status()
  quit(save = "no")
}

if (mode == "report") {
  p <- ca_did_results_path("coverage")
  if (!file.exists(p)) {
    stop("No coverage table. Run: Rscript 10c_aggregate_did.R aggregate")
  }
  cov <- utils::read.csv(p, colClasses = "character", stringsAsFactors = FALSE)
  print(cov, right = FALSE)
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# READ ONE ESTIMATE FILE
# ---------------------------------------------------------------------------
# EVERY COLUMN AS CHARACTER, THEN COERCED BY NAME. Three reasons, all seen in
# the written files. Ecoregion "08" parses to the integer 8. An empty stratum
# parses to NA rather than "". And p.value is written in fixed notation with
# several hundred decimal places on the largest fire cells, which is valid but
# is not something to leave to type inference.
#
# The two estimators write different columns and both are handled here rather
# than by two readers, because the difference is three columns and a time
# variable, not two schemas.
#
#   attgt   aggregation term egt estimate std.error crit_val conf.low
#           conf.high statistic p.value + identity + wald_pval. aggregation
#           is group, dynamic, or calendar, the last on the PA arm alone.
#   etwfe   term contrast .Dtreat estimate std.error statistic p.value
#           s.value conf.low conf.high aggregation [year] [event] + identity
#           + emfx_type
#
# year and event are present only for the aggregations that produce them, so a
# PA run carries both and an offset run carries event alone. They are coalesced
# into one time column with time_type recording which it was.

CANON <- c("run_id", "estimator", "pairing", "outcome", "transform",
           "stratum", "ecoregion", "base_shift", "aggregation", "term",
           "time_type", "time", "estimate", "std.error", "statistic",
           "p.value", "conf.low", "conf.high", "crit_val")

read_est <- function(path, estimator) {

  d <- tryCatch(
    utils::read.csv(path, colClasses = "character", check.names = FALSE,
                    stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)

  col <- function(nm) {
    if (nm %in% names(d)) d[[nm]] else rep(NA_character_, nrow(d))
  }
  num <- function(x) suppressWarnings(as.numeric(x))

  agg  <- col("aggregation")
  term <- col("term")

  time_type <- rep(NA_character_, nrow(d))
  time      <- rep(NA_real_, nrow(d))

  if (identical(estimator, "attgt")) {
    # egt is the cohort year under the group aggregation, event time under the
    # dynamic one, and the calendar year under the calendar one, which the PA
    # arm alone reports. Same column, three meanings, so time_type carries the
    # distinction rather than the reader assuming one of them.
    egt <- num(col("egt"))
    time_type[!is.na(agg) & agg == "group"]    <- "cohort"
    time_type[!is.na(agg) & agg == "dynamic"]  <- "event"
    time_type[!is.na(agg) & agg == "calendar"] <- "calendar"
    time <- egt
    ov <- !is.na(term) & term == "overall"
    time_type[ov] <- "overall"
    time[ov]      <- NA_real_
  } else {
    yr <- num(col("year"))
    ev <- num(col("event"))
    time_type[!is.na(agg) & agg == "calendar"] <- "calendar"
    time_type[!is.na(agg) & agg == "event"]    <- "event"
    time_type[!is.na(agg) & agg == "simple"]   <- "overall"
    time <- ifelse(!is.na(agg) & agg == "calendar", yr,
                   ifelse(!is.na(agg) & agg == "event", ev, NA_real_))
  }

  out <- data.frame(
    run_id      = col("run_id"),
    estimator   = estimator,
    pairing     = col("pairing"),
    outcome     = col("outcome"),
    transform   = col("transform"),
    stratum     = col("stratum"),
    ecoregion   = col("ecoregion"),
    base_shift  = as.integer(num(col("base_shift"))),
    aggregation = agg,
    term        = term,
    time_type   = time_type,
    time        = time,
    estimate    = num(col("estimate")),
    std.error   = num(col("std.error")),
    statistic   = num(col("statistic")),
    p.value     = num(col("p.value")),
    conf.low    = num(col("conf.low")),
    conf.high   = num(col("conf.high")),
    crit_val    = num(col("crit_val")),
    stringsAsFactors = FALSE)

  # An empty stratum written as "" comes back as "" and an unwritten one as NA.
  # Both mean statewide across strata, so they are made the same thing.
  out$stratum[is.na(out$stratum)]     <- ""
  out$ecoregion[is.na(out$ecoregion)] <- ""

  out[, CANON, drop = FALSE]
}


# ---------------------------------------------------------------------------
# BIND
# ---------------------------------------------------------------------------

done <- STATUS[STATUS$class == "done", , drop = FALSE]
if (!nrow(done)) {
  stop("No completed runs found. Run the check mode first.")
}

ca_log("Reading ", nrow(done), " files")

parts <- vector("list", nrow(done))
for (i in seq_len(nrow(done))) {
  parts[[i]] <- read_est(done$path[i], done$estimator[i])
  if (i %% 200L == 0L) ca_log("  ", i, " of ", nrow(done))
}

bad <- vapply(parts, is.null, logical(1))
if (any(bad)) {
  ca_log("WARNING, ", sum(bad), " files read as empty and were skipped")
  for (i in which(bad)) ca_log("  ", done$path[i])
}

RES <- do.call(rbind, parts[!bad])
rownames(RES) <- NULL
rm(parts); ca_gc("after bind")

ca_log("Bound ", format(nrow(RES), big.mark = ","), " estimate rows")


# ---------------------------------------------------------------------------
# TREATMENT ARM AND MANIFEST CONTEXT
# ---------------------------------------------------------------------------
# The manifest carries what the estimate file cannot, which is what the cell
# was built from. A caption that states how many pairs a line rests on reads it
# from here rather than from a separate lookup.

RES$arm <- ca_pairings$arm[match(RES$pairing, ca_pairings$pairing)]

MAN_KEEP <- c("run_id", "status", "arm", "config", "tier",
              "year_first", "year_last", "n_rows", "n_units", "n_pairs",
              "n_treated_units", "n_cohorts", "gvar_min", "gvar_max")

man_all <- do.call(rbind, lapply(ESTIMATORS, function(es) {
  m <- read_manifest(es)
  if (is.null(m)) return(NULL)
  for (nm in MAN_KEEP) if (!nm %in% names(m)) m[[nm]] <- NA_character_
  m <- m[, MAN_KEEP, drop = FALSE]
  m$estimator <- es
  m
}))

MAN_NUM <- c("year_first", "year_last", "n_rows", "n_units", "n_pairs",
             "n_treated_units", "n_cohorts", "gvar_min", "gvar_max")

if (!is.null(man_all)) {
  # The manifest arm is the matching arm, baseline or augmented. RES$arm is the
  # treatment arm. Renamed rather than resolved at read time so the collision is
  # visible in one place.
  names(man_all)[names(man_all) == "arm"] <- "match_arm"
  idx <- match(paste(RES$estimator, RES$run_id),
               paste(man_all$estimator, man_all$run_id))
  for (nm in c("status", "match_arm", "config", "tier", MAN_NUM)) {
    RES[[nm]] <- man_all[[nm]][idx]
  }
  for (nm in MAN_NUM) {
    RES[[nm]] <- suppressWarnings(as.numeric(RES[[nm]]))
  }
}


# ---------------------------------------------------------------------------
# BACK-TRANSFORM
# ---------------------------------------------------------------------------
# exp(beta) - 1 on the log1p rows, through expm1 for accuracy near zero. The
# transform is monotone, so the confidence bounds map directly and are not
# recomputed. The result is a proportion. Multiply by 100 at the figure, not
# here, so nothing in this table sits on a plotting scale.
#
# The submitted figures plotted estimate * 100, which is the log-point
# approximation. It is close at small effects and is not close at the
# magnitudes the fire arm reports, where a coefficient near -1 is -63 percent
# and not -100. This is the correction that moves the headline numbers.
#
# Non-log rows keep NA in the back-transformed columns rather than a copy of
# the estimate. A figure that reaches for estimate_bt on a native-scale outcome
# should fail visibly rather than plot the wrong units.

is_log <- !is.na(RES$transform) & RES$transform == "log1p"

RES$scale        <- ifelse(is_log, "proportion", "native")
RES$estimate_bt  <- ifelse(is_log, expm1(RES$estimate),  NA_real_)
RES$conf.low_bt  <- ifelse(is_log, expm1(RES$conf.low),  NA_real_)
RES$conf.high_bt <- ifelse(is_log, expm1(RES$conf.high), NA_real_)


# ---------------------------------------------------------------------------
# SE RATIO
# ---------------------------------------------------------------------------
# Standard error against the standard error at the preceding time point, within
# one run and one aggregation. The cap on an event-time axis is the last point
# before this exceeds two, which is where the contributing cohort set has
# thinned enough that the estimate describes a different population. Measured
# on the offset arm in the Coast Range, GPP breaks between event 10 and 11 at
# 2.3 and the NCSDA fluxes break between 6 and 7 at 2.6, which is where the
# submitted caps already sat.
#
# The same column catches the opposite case. A terminal ratio below one half
# marks the variance collapse seen in UP_UNP GPP at 2025 and NEP at 2021, which
# is inspected before it reaches a figure rather than after a reviewer finds it.
#
# NA on the first point of a series, on the overall rows, and wherever the
# preceding standard error is zero or absent, which includes the emfx event
# reference period where the estimate is zero by construction.

RES <- RES[order(RES$estimator, RES$run_id, RES$aggregation, RES$time,
                 na.last = TRUE), , drop = FALSE]
rownames(RES) <- NULL

key      <- paste(RES$estimator, RES$run_id, RES$aggregation)
prev_se  <- c(NA_real_, utils::head(RES$std.error, -1L))
same_run <- c(FALSE, key[-1L] == utils::head(key, -1L))

RES$se_ratio_prev <- ifelse(
  same_run & !is.na(prev_se) & prev_se > 0 & !is.na(RES$std.error) &
    !is.na(RES$time_type) & RES$time_type != "overall",
  RES$std.error / prev_se, NA_real_)


# ---------------------------------------------------------------------------
# WRITE
# ---------------------------------------------------------------------------

OUT_COLS <- c("run_id", "estimator", "arm", "pairing", "outcome", "transform",
              "scale", "stratum", "ecoregion", "aggregation", "time_type",
              "time", "estimate", "std.error", "conf.low", "conf.high",
              "estimate_bt", "conf.low_bt", "conf.high_bt", "statistic",
              "p.value", "crit_val", "se_ratio_prev", "base_shift",
              "match_arm", "config", "tier", "year_first", "year_last",
              "n_rows", "n_units", "n_pairs", "n_treated_units", "n_cohorts",
              "gvar_min", "gvar_max")
OUT_COLS <- OUT_COLS[OUT_COLS %in% names(RES)]

res_path <- ca_did_results_path("results")
con <- gzfile(res_path, "w")
utils::write.csv(RES[, OUT_COLS, drop = FALSE], con, row.names = FALSE)
close(con)
ca_log("Wrote ", res_path, "  ",
       format(nrow(RES), big.mark = ","), " rows, ",
       round(file.info(res_path)$size / 1e6, 1), " MB")


# ---------------------------------------------------------------------------
# COVERAGE
# ---------------------------------------------------------------------------
# Completeness at the moment the results table was written, one row per
# estimator and treatment arm. A figure produced from an incomplete aggregate
# is identifiable from this file rather than from memory, which matters while
# the large tier is still draining.

cov <- do.call(rbind, lapply(
  split(STATUS, list(STATUS$estimator, STATUS$arm), drop = TRUE),
  function(s) {
    data.frame(
      estimator = s$estimator[1],
      arm       = s$arm[1],
      expected  = nrow(s),
      done      = sum(s$class == "done"),
      running   = sum(s$class == "running"),
      todo      = sum(s$class == "todo"),
      skipped   = sum(s$class == "skipped"),
      error     = sum(s$class == "error"),
      pct_done  = round(100 * sum(s$class == "done") / nrow(s), 1),
      stringsAsFactors = FALSE)
  }))
cov <- cov[order(match(cov$estimator, ESTIMATORS), cov$arm), , drop = FALSE]
cov$written_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
rownames(cov) <- NULL

utils::write.csv(cov, ca_did_results_path("coverage"), row.names = FALSE)
ca_log("Wrote ", ca_did_results_path("coverage"))

write_status()

cat("\n")
print(cov, right = FALSE)
cat("\n")

flag <- RES[!is.na(RES$se_ratio_prev) & RES$se_ratio_prev < 0.5, , drop = FALSE]
if (nrow(flag)) {
  cat("SE COLLAPSE, ratio below 0.5 against the preceding time point.\n")
  cat("Inspect before plotting.\n")
  print(utils::head(flag[, c("run_id", "aggregation", "time", "std.error",
                             "se_ratio_prev")], 40), right = FALSE)
  cat("\n", nrow(flag), " rows flagged in total.\n\n", sep = "")
}

ca_stamp("10c_aggregate_did")
ca_log("Stage 10c complete")
