### CA carbon revision pipeline
### Stage 10b. ETWFE and emfx, one file per run
###
### WHAT THIS STAGE DOES
###
### Reads one stage 9 panel, subsets it to one cell, applies the transform,
### fits etwfe(), and writes the emfx aggregation stage 11 plots. One array
### task is one run and nothing loops inside it.
###
### It replaces every 10_*_sub_* file. Those differed only in pairing, log
### against non-log, ecoregion, and severity or intensity, so they are the
### columns of ca_did_runs() rather than forty files carrying the same four
### edits.
###
### THE ARRAY
###
### ca_did_runs("etwfe"). Panel manifest crossed with transform crossed with
### subset cell, ordered by tier so each tier is a contiguous range. Cells come
### from meta/did_subset_levels.csv, which this script builds in levels mode
### from the stage 9 units files, so the ecoregions and strata are the ones a
### pairing actually holds rather than the hardcoded c(78, 1, 9, 4, 5, 8) of
### 10_1_..._sub_2.R L43.
###
### THE AGGREGATION IS THE ARM'S, NOT A CHOICE MADE HERE
###
### ca_did_emfx_types() returns calendar and event for the PA arm and event
### alone for offset, fire, and thinning. The single-type version is read from
### the submitted files at 10_1 L239, 10_5 L187, 10_2 L146, and 10_3 L140. The
### PA arm gained event on 2026-08-13 so protection carries the same
### pre-treatment panel as every other arm. See ca_config.R section 0.11b.
###
### THE THINNING BASE PERIOD
###
### gvar_est = gvar - 1 on TP_UP and TNP_UNP, through ca_did_gvar_est(), so the
### base period is c-2 rather than the detection year. See ca_config.R section
### 0.11c. Estimation only, the panel is untouched.
###
### NO COVARIATES, NO SUBSAMPLE
###
### fml = y ~ 0. Matching balanced the covariates and the submitted files that
### carried them in col.list left the xformla line commented out. The 50,000
### subclass draw at 10_2 L69 is gone, replaced by the stage 8 cap.
###
### Usage
###   Rscript 10b_did_etwfe.R levels    build meta/did_subset_levels.csv
###   Rscript 10b_did_etwfe.R tasks     print the array and tier ranges
###   Rscript 10b_did_etwfe.R           run one array task
###   Rscript 10b_did_etwfe.R missing   print an sbatch array of unfinished runs
###   Rscript 10b_did_etwfe.R report    wall time against cost, per tier
###   Rscript 10b_did_etwfe.R archive   copy completed output to /projects
###
###   CA_OVERWRITE=1                    rerun cells that already have output

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr", "etwfe", "broom"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# THREADS.
#
# fixest parallelises through OpenMP and, left alone, counts the physical cores
# of the machine rather than the cores Slurm allocated. On a 128 core node with
# eight cores allocated that is sixteen times oversubscription inside a cgroup
# that will not give the extra cores, so the threads fight each other for the
# eight that exist. Every thread count in this stage comes from
# SLURM_CPUS_PER_TASK through ca_threads().
if (requireNamespace("fixest", quietly = TRUE)) {
  fixest::setFixest_nthreads(ca_threads(), save = FALSE)
}
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(ca_threads())
}
ca_log("threads ", ca_threads())

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("levels", "tasks", "archive", "missing", "report")) {
  args[1]
} else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

ESTIMATOR <- "etwfe"


# ---------------------------------------------------------------------------
# LEVELS
# ---------------------------------------------------------------------------
# One row per pairing and subset cell, built from the stage 9 units files. The
# units files are small, carry the stratum already propagated from the treated
# member to its control, and carry ecoregion_l3, which is an exact matching
# variable and therefore identical within a pair. So a cell is a clean subset
# of whole pairs in both dimensions.
#
# The statewide all-stratum row is written for every pairing. It is a run for
# the pa and offset arms and the denominator for the row estimate on the
# others, which is why it exists even where it is not estimated.

build_levels <- function() {

  UNIT_COLS <- c("unit_id", "pair_id", "treat", "gvar", "stratum",
                 "ecoregion_l3")

  cell_row <- function(u, pairing, stratum, ecoregion) {
    tr <- u$treat == 1L
    g  <- u$gvar[tr]
    g  <- g[!is.na(g) & g > 0L]
    data.frame(
      pairing   = pairing,
      stratum   = stratum,
      ecoregion = ecoregion,
      n_units   = nrow(u),
      n_treated = sum(tr),
      n_pairs   = length(unique(u$pair_id)),
      n_cohorts = length(unique(g)),
      gvar_min  = if (length(g)) min(g) else NA_integer_,
      gvar_max  = if (length(g)) max(g) else NA_integer_,
      stringsAsFactors = FALSE)
  }

  out <- list()

  for (pr in ca_pairings$pairing) {

    p <- ca_did_units_path(pr)
    if (!file.exists(p)) {
      stop("Missing units file: ", p, ". Run 9a_did_units.R first.")
    }

    u <- as.data.frame(dplyr::collect(dplyr::select(
      arrow::open_dataset(p), dplyr::all_of(UNIT_COLS))))
    u$stratum      <- as.character(u$stratum)
    u$ecoregion_l3 <- as.character(u$ecoregion_l3)
    u$gvar         <- as.integer(u$gvar)

    tr <- u$treat == 1L

    st <- unique(u$stratum[tr])
    st <- st[!is.na(st) & nzchar(st)]
    # Config order where the pairing declares strata, so Low, Moderate, High
    # rather than alphabetical.
    dec <- ca_pairings$strata[match(pr, ca_pairings$pairing)]
    if (!is.na(dec)) {
      dec <- trimws(strsplit(dec, ",", fixed = TRUE)[[1]])
      st <- st[order(match(st, dec))]
    } else {
      st <- sort(st)
    }

    ec <- unique(u$ecoregion_l3[tr])
    ec <- sort(ec[!is.na(ec) & nzchar(ec)])

    rows <- list(cell_row(u, pr, "", ""))

    for (s in st) {
      us <- u[!is.na(u$stratum) & u$stratum == s, , drop = FALSE]
      if (!sum(us$treat == 1L)) next
      rows[[length(rows) + 1L]] <- cell_row(us, pr, s, "")
    }

    for (e in ec) {
      ue <- u[!is.na(u$ecoregion_l3) & u$ecoregion_l3 == e, , drop = FALSE]
      if (!sum(ue$treat == 1L)) next
      rows[[length(rows) + 1L]] <- cell_row(ue, pr, "", e)
      for (s in st) {
        us <- ue[!is.na(ue$stratum) & ue$stratum == s, , drop = FALSE]
        if (!sum(us$treat == 1L)) next
        rows[[length(rows) + 1L]] <- cell_row(us, pr, s, e)
      }
    }

    d <- do.call(rbind, rows)
    ca_log(pr, ": ", nrow(d), " cells, ", length(st), " strata, ",
           length(ec), " ecoregions")
    out[[length(out) + 1L]] <- d
    rm(u); ca_gc(paste0("after ", pr))
  }

  d <- do.call(rbind, out)
  d$built_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  utils::write.csv(d, ca_did_levels_path(), row.names = FALSE)
  ca_log("Wrote ", ca_did_levels_path(), "  ", nrow(d), " cells")
  invisible(d)
}

if (mode == "levels") {
  build_levels()
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_did_est_dir(), "model", subdir = "work")
  quit(save = "no")
}

RUNS <- ca_did_runs(ESTIMATOR)

if (mode == "tasks") {
  show <- RUNS[, c("tier", "pairing", "outcome", "transform", "stratum",
                   "ecoregion", "n_treated", "n_cohorts", "est_rows",
                   "est_cost", "emfx_type")]
  print(show, right = FALSE)
  cat("\narray size:", nrow(RUNS), "\n")
  for (tr in c("small", "medium", "large")) {
    i <- which(RUNS$tier == tr)
    if (!length(i)) next
    cat(sprintf("%-7s %d-%d  (%d runs, cost %s to %s)\n", tr,
                min(i), max(i), length(i),
                format(min(RUNS$est_cost[i]), big.mark = ","),
                format(max(RUNS$est_cost[i]), big.mark = ",")))
  }
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# MISSING AND REPORT
# ---------------------------------------------------------------------------
# WHAT COUNTS AS DONE. A present file is not enough. A task killed mid-write
# leaves a truncated csv.gz, and the run path skips any cell whose output
# already exists, so a truncated file would never be recomputed and would go
# on to be read by stage 11. Done means the gzip opens and carries a header
# plus at least one row.
#
# WHAT IS NOT MISSING. A cell with no viable cohort writes a manifest row with
# status skipped and no output, by design, and must not sit in the retry array
# forever. A cell whose estimator errored is also not a memory problem and is
# listed separately rather than resubmitted blindly.
#
# WHAT IS STILL RUNNING. squeue is read where it exists, so a task in flight
# is not handed back for resubmission on another node. This is the case that
# matters while an array is draining.

JOB_NAME <- "ca_10b_etwfe"

if (mode %in% c("missing", "report")) {

  paths <- vapply(seq_len(nrow(RUNS)), function(i) {
    ca_did_est_path(ESTIMATOR, RUNS$pairing[i], RUNS$run_id[i])
  }, character(1))

  present <- file.exists(paths)
  valid   <- rep(FALSE, length(paths))
  for (i in which(present)) {
    sz <- file.info(paths[i])$size
    if (is.na(sz) || sz < 50) next
    valid[i] <- tryCatch({
      con <- gzfile(paths[i], "r")
      on.exit(close(con), add = TRUE)
      length(readLines(con, n = 2L)) >= 2L
    }, error = function(e) FALSE)
  }

  mpath <- ca_did_est_manifest_path(ESTIMATOR)
  man <- if (file.exists(mpath)) {
    utils::read.csv(mpath, colClasses = "character")
  } else NULL
  st <- if (is.null(man)) rep(NA_character_, nrow(RUNS)) else {
    man$status[match(RUNS$run_id, man$run_id)]
  }

  running <- integer(0)
  if (nzchar(Sys.which("squeue"))) {
    txt <- tryCatch(
      system2("squeue", c("-h", "-u", Sys.getenv("USER"), "-n", JOB_NAME,
                          "-o", "%K"), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0))
    running <- suppressWarnings(as.integer(trimws(txt)))
    running <- running[!is.na(running)]
  }

  msg <- if (is.null(man)) rep("", nrow(RUNS)) else {
    m <- man$message[match(RUNS$run_id, man$run_id)]
    m[is.na(m)] <- ""
    m
  }

  # A RESOURCE ERROR IS NOT AN ESTIMATOR ERROR.
  #
  # The fit is wrapped in tryCatch, so a task that ran out of address space
  # inside feols writes a manifest row with status error and exits zero, which
  # the scheduler records as COMPLETED in two minutes. Those rows say "cannot
  # allocate vector of size 88 Gb" and they belong in the retry array, not in
  # a list of things a bigger node cannot fix. Only an error that is not about
  # memory is a real estimator failure.
  res <- grepl("cannot allocate|out of memory|memory exhausted|std::bad_alloc",
               msg, ignore.case = TRUE)

  cls <- rep("todo", nrow(RUNS))
  cls[!valid & !is.na(st) & st == "skipped"] <- "skipped"
  cls[!valid & !is.na(st) & st == "error" & !res] <- "error"
  cls[seq_len(nrow(RUNS)) %in% running & !valid] <- "running"
  cls[valid] <- "done"

  truncated <- which(present & !valid)

  cat(nrow(RUNS), "runs\n")
  print(table(cls, RUNS$tier))

  if (length(truncated)) {
    cat("\nTRUNCATED OUTPUT, delete before resubmitting these tasks:\n")
    cat(paste0("  ", paths[truncated], collapse = "\n"), "\n")
  }

  if (any(res & cls == "todo")) {
    cat("\n", sum(res & cls == "todo"),
        " runs failed on memory inside the fit and are in the retry array.\n",
        "  largest request seen: ",
        max(as.numeric(sub(".*size ([0-9.]+) Gb.*", "\\1",
            grep("size .* Gb", msg[res & cls == "todo"], value = TRUE))),
            na.rm = TRUE), " Gb\n", sep = "")
  }

  if (any(cls == "error")) {
    cat("\nESTIMATOR ERRORS, not a resource problem:\n")
    e <- which(cls == "error")
    for (i in e) {
      cat("  ", RUNS$run_id[i], "  ",
          if (is.null(man)) "" else man$message[match(RUNS$run_id[i],
                                                     man$run_id)], "\n")
    }
  }

  if (mode == "report") {
    if (!is.null(man)) {
      d <- man[man$status == "ok", , drop = FALSE]
      d$wall_min <- as.numeric(d$wall_min)
      d$cost <- RUNS$est_cost[match(d$run_id, RUNS$run_id)]
      d$tier <- RUNS$tier[match(d$run_id, RUNS$run_id)]
      d <- d[!is.na(d$cost), , drop = FALSE]
      cat("\nCOMPLETED RUNS, wall minutes against cost\n")
      for (tr in c("small", "medium", "large")) {
        x <- d[d$tier == tr, , drop = FALSE]
        if (!nrow(x)) next
        cat(sprintf("%-7s n %4d  cost %s to %s  wall %.1f / %.1f / %.1f min",
                    tr, nrow(x),
                    format(min(x$cost), big.mark = ","),
                    format(max(x$cost), big.mark = ","),
                    stats::median(x$wall_min), stats::quantile(x$wall_min,
                    0.9, names = FALSE), max(x$wall_min)))
        cat(sprintf("  minutes per million cost %.3f\n",
                    stats::median(x$wall_min / (x$cost / 1e6))))
      }
    }
    quit(save = "no")
  }

  idx <- which(cls == "todo")
  cat("\n", length(idx), " runs to submit",
      if (length(running)) paste0(", ", sum(cls == "running"),
                                  " still in the queue") else "", "\n",
      sep = "")
  # One array string per tier. A retry sizes its memory on the tier it is
  # retrying, and mixing a 20 million cost cell into a 950 GB request wastes a
  # node that something else needs.
  arr <- function(v) {
    g <- cumsum(c(1L, as.integer(diff(v) != 1L)))
    paste(vapply(split(v, g), function(x) {
      if (length(x) == 1L) as.character(x) else paste0(min(x), "-", max(x))
    }, character(1)), collapse = ",")
  }
  for (tr in c("small", "medium", "large")) {
    v <- idx[RUNS$tier[idx] == tr]
    if (!length(v)) next
    cat("\n", tr, ", ", length(v), " runs, cost ",
        format(min(RUNS$est_cost[v]), big.mark = ","), " to ",
        format(max(RUNS$est_cost[v]), big.mark = ","), "\n", sep = "")
    cat("--array=", arr(v), "\n", sep = "")
  }
  quit(save = "no")
}

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(RUNS)) {
  stop("Task id must be 1 to ", nrow(RUNS), ", got ", task)
}

R <- RUNS[task, ]
t_start <- Sys.time()

ca_log("Stage 10b, task ", task, ".  ", R$run_id)
ca_log("  tier ", R$tier, "  est rows ", format(R$est_rows, big.mark = ","),
       "  cost ", format(R$est_cost, big.mark = ","),
       "  emfx type ", R$emfx_type, "  base shift ", R$base_shift)

out_path <- ca_did_est_path(ESTIMATOR, R$pairing, R$run_id)
if (file.exists(out_path) && !OVERWRITE) {
  ca_log("Present, not overwritten: ", basename(out_path))
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# READ THE CELL
# ---------------------------------------------------------------------------
# The subset is pushed into the Arrow scan so a one-ecoregion run does not
# decompress the whole panel. Dictionary-encoded columns compare against a
# string, but a version that refuses falls back to an in-memory filter rather
# than failing the task.

PANEL_COLS <- c("unit_id", "pair_id", "treat", "gvar", "year", "value",
                "stratum", "ecoregion_l3")

read_cell <- function() {
  ds <- arrow::open_dataset(R$path)
  q  <- dplyr::select(ds, dplyr::all_of(PANEL_COLS))
  st <- R$stratum
  ec <- R$ecoregion
  pushed <- tryCatch({
    if (nzchar(st)) q <- dplyr::filter(q, stratum == st)
    if (nzchar(ec)) q <- dplyr::filter(q, ecoregion_l3 == ec)
    as.data.frame(dplyr::collect(q))
  }, error = function(e) {
    ca_log("  scan pushdown refused (", conditionMessage(e),
           "), filtering in memory")
    NULL
  })
  if (!is.null(pushed)) return(pushed)
  d <- as.data.frame(dplyr::collect(dplyr::select(ds,
        dplyr::all_of(PANEL_COLS))))
  if (nzchar(st)) d <- d[as.character(d$stratum) == st, , drop = FALSE]
  if (nzchar(ec)) d <- d[as.character(d$ecoregion_l3) == ec, , drop = FALSE]
  d
}

d <- read_cell()
d$stratum <- NULL
d$ecoregion_l3 <- NULL
ca_gc("cell read")

status  <- "ok"
message <- ""

fail <- function(msg) {
  status  <<- "skipped"
  message <<- msg
  ca_log("  SKIP ", msg)
}

if (!nrow(d)) fail("no rows in this cell")


# ---------------------------------------------------------------------------
# COHORTS
# ---------------------------------------------------------------------------

n_cohorts_dropped <- 0L
n_always_treated  <- 0L
n_half_pairs      <- 0L
yr_first <- NA_integer_
yr_last  <- NA_integer_

if (status == "ok") {

  d$year     <- as.integer(d$year)
  d$gvar_est <- ca_did_gvar_est(R$pairing, d$gvar)

  yr_first <- min(d$year)
  yr_last  <- max(d$year)

  # viable_cohorts_only. A cohort starting after the last year of this outcome
  # has no post-treatment observation. It contributes no ATT and it is not in
  # the never-treated group either, so it is removed rather than carried.
  coh <- sort(unique(d$gvar_est[d$treat == 1L]))
  bad <- coh[coh > yr_last]
  n_cohorts_dropped <- length(bad)
  if (n_cohorts_dropped) {
    drop_id <- unique(d$unit_id[d$treat == 1L & d$gvar_est %in% bad])
    ca_log("  cohorts after ", yr_last, " dropped: ", n_cohorts_dropped,
           " (", format(length(drop_id), big.mark = ","), " treated units)")
    d <- d[!(d$unit_id %in% drop_id), , drop = FALSE]
  }

  # Reported, not acted on. A cohort at the first panel year is treated in
  # every observed period, which is the PA floor at 1985 and is the submitted
  # specification. Its effect is identified against never-treated controls
  # through the matched level difference.
  n_always_treated <- sum(unique(d$gvar_est[d$treat == 1L]) == yr_first)

  # complete_pairs_only is FALSE. Counted so the LEMMA data note has its
  # number, and nothing is deleted.
  p1 <- unique(d$pair_id[d$treat == 1L])
  p0 <- unique(d$pair_id[d$treat == 0L])
  n_half_pairs <- length(setdiff(p1, p0)) + length(setdiff(p0, p1))
  rm(p1, p0)

  if (!sum(d$treat == 1L)) fail("no treated units after the cohort filter")
  if (status == "ok" && !any(d$gvar_est == 0L)) {
    fail("no never-treated units in this cell")
  }
  if (status == "ok" && !length(unique(d$gvar_est[d$treat == 1L]))) {
    fail("no viable cohort in this cell")
  }
}


# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------

est        <- NULL
n_terms    <- NA_integer_
att_simple <- NA_real_
se_simple  <- NA_real_
n_rows   <- if (status == "ok") nrow(d) else 0L
n_units  <- if (status == "ok") length(unique(d$unit_id)) else 0L
n_pairs  <- if (status == "ok") length(unique(d$pair_id)) else 0L
n_treat  <- if (status == "ok") length(unique(d$unit_id[d$treat == 1L])) else 0L
n_coh    <- if (status == "ok") length(unique(d$gvar_est[d$treat == 1L])) else 0L
gv_min   <- NA_integer_
gv_max   <- NA_integer_

if (status == "ok") {

  d$y <- ca_did_apply_transform(d$value, R$transform)
  d$value <- NULL
  d$gvar <- NULL
  d$treat <- NULL
  # unit_id is not passed to etwfe, since ivar stays NULL as in the submitted
  # calls, and every count that needs it has already been taken. It is dropped
  # here because emfx holds the model frame alongside its own jacobian and the
  # peak is what gets a task killed.
  d$unit_id <- NULL
  ca_gc("transform applied")

  g <- d$gvar_est[d$gvar_est > 0L]
  gv_min <- min(g); gv_max <- max(g)

  ca_log("  rows ", format(n_rows, big.mark = ","),
         "  units ", format(n_units, big.mark = ","),
         "  pairs ", format(n_pairs, big.mark = ","),
         "  treated ", format(n_treat, big.mark = ","),
         "  cohorts ", n_coh, " (", gv_min, " to ", gv_max, ")")

  set.seed(ca_did_est$seed)

  fit <- tryCatch({
    mod <- etwfe::etwfe(
      fml    = y ~ 0,
      tvar   = year,
      gvar   = gvar_est,
      data   = d,
      vcov   = stats::as.formula(paste0("~", ca_did$cluster)),
      cgroup = ca_did_est$cgroup)
    ca_log("  etwfe fitted")

    # THE ARM'S OWN AGGREGATIONS, AND NOTHING ELSE.
    #
    # ca_did_emfx_types() gives calendar and event for protection and event
    # alone for offset, fire, and thinning. The task table column emfx_type
    # carries the arm's primary aggregation and drives tiering and logging.
    # The aggregation column in the output is the row-level truth and is what
    # stage 11 filters on.
    #
    # type = "simple" is not run. It was added because 10_2 L137 and 10_5 L181
    # called emfx untyped, but 10_1 L233 and 10_3 L126 had that call commented
    # out, so the submitted results never carried an overall ATT for protection
    # or thinning and the manuscript does not rest on one. It is a full second
    # pass through marginaleffects, and on task 1217 it cost 5 hours 48 minutes
    # of an 18 hour run. The claims in this paper are about trajectories, the
    # short-term trough and the recovery path, which the event and calendar
    # aggregations carry directly.
    #
    # The att_simple and se_simple manifest columns stay, holding NA, because
    # the manifest already written for 1,162 completed runs has them and the
    # append path refuses a changed column set.
    parts <- lapply(ca_did_emfx_types(R$pairing), function(ty) {
      m <- etwfe::emfx(mod, type = ty,
                       compress = ca_did_est$emfx_compress,
                       vcov = ca_did_est$emfx_vcov)
      ca_log("  emfx ", ty, " done")
      out <- as.data.frame(broom::tidy(m))
      out$aggregation <- ty
      # The marginal effects object carries a per-observation jacobian and on a
      # three million row cell it is gigabytes. Released before the next
      # aggregation starts rather than at the end of the lapply, since the
      # second call is where the peak lands.
      rm(m)
      ca_gc(paste0("after emfx ", ty))
      out
    })
    dplyr::bind_rows(parts)
  }, error = function(e) {
    status  <<- "error"
    message <<- conditionMessage(e)
    ca_log("  ERROR ", conditionMessage(e))
    NULL
  })

  if (!is.null(fit)) {
    est <- as.data.frame(fit)
    n_terms <- nrow(est)
    sm <- est[est$aggregation == "simple", , drop = FALSE]
    if (nrow(sm) == 1L) {
      att_simple <- as.numeric(sm$estimate)
      se_simple  <- as.numeric(sm$std.error)
      ca_log("  overall ATT ", signif(att_simple, 6), " (SE ",
             signif(se_simple, 4), ")")
    }
  }
}


# ---------------------------------------------------------------------------
# WRITE
# ---------------------------------------------------------------------------
# Identity columns ride on every row, so a stage 11 read is a bind of the
# whole directory followed by a filter rather than a filename parse.

if (!is.null(est) && nrow(est)) {
  est$run_id     <- R$run_id
  est$estimator  <- ESTIMATOR
  est$pairing    <- R$pairing
  est$outcome    <- R$outcome
  est$transform  <- R$transform
  est$stratum    <- R$stratum
  est$ecoregion  <- R$ecoregion
  est$emfx_type  <- R$emfx_type
  est$base_shift <- R$base_shift

  # Written to a temporary name in the same directory and renamed on success.
  # rename is atomic within a filesystem, so a reader never sees a partial
  # file and two tasks racing on the same cell cannot interleave into one
  # corrupt gzip. That matters because a retry can be queued on a second
  # partition while the first is still running.
  tmp_path <- paste0(out_path, ".tmp", Sys.getpid())
  con <- gzfile(tmp_path, "w")
  utils::write.csv(est, con, row.names = FALSE)
  close(con)
  if (!file.rename(tmp_path, out_path)) {
    unlink(tmp_path)
    stop("Could not move ", tmp_path, " into place")
  }
  ca_log("Wrote ", basename(out_path), "  ", nrow(est), " terms")
} else {
  out_path <- NA_character_
}


# ---------------------------------------------------------------------------
# MANIFEST
# ---------------------------------------------------------------------------

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

sel <- ca_match_selected(R$pairing)

man <- data.frame(
  run_id            = R$run_id,
  estimator         = ESTIMATOR,
  pairing           = R$pairing,
  outcome           = R$outcome,
  transform         = R$transform,
  stratum           = R$stratum,
  ecoregion         = R$ecoregion,
  arm               = sel$arm,
  config            = sel$config,
  emfx_type         = R$emfx_type,
  base_shift        = R$base_shift,
  status            = status,
  tier              = R$tier,
  year_first        = yr_first,
  year_last         = yr_last,
  n_rows            = n_rows,
  n_units           = n_units,
  n_pairs           = n_pairs,
  n_treated_units   = n_treat,
  n_cohorts         = n_coh,
  n_cohorts_dropped = n_cohorts_dropped,
  n_always_treated  = n_always_treated,
  n_half_pairs      = n_half_pairs,
  gvar_min          = gv_min,
  gvar_max          = gv_max,
  n_terms           = n_terms,
  att_simple        = att_simple,
  se_simple         = se_simple,
  wall_min          = round(as.numeric(difftime(Sys.time(), t_start,
                                                units = "mins")), 2),
  path              = out_path,
  written_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  message           = message,
  stringsAsFactors  = FALSE
)

lock <- file.path(ca_did_est_dir(), ".lock_manifest_etwfe")
acquire_lock(lock)
on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

# READ THE MANIFEST BACK AS CHARACTER.
#
# read.csv() with default colClasses reparses every column of every row already
# written. An ecoregion of "08" returns as the number 8 and an empty stratum
# returns as NA, so the rows written by earlier tasks silently lose their
# zero padding and their empty strings. With 1,245 tasks appending in turn the
# damage compounds. Character in, character out, and the numbers are numbers
# again wherever they are read for analysis.
mpath <- ca_did_est_manifest_path(ESTIMATOR)
man <- data.frame(lapply(man, as.character), stringsAsFactors = FALSE)
if (file.exists(mpath)) {
  old <- utils::read.csv(mpath, colClasses = "character")
  old <- old[old$run_id != R$run_id, , drop = FALSE]
  if (nrow(old) && !identical(sort(names(old)), sort(names(man)))) {
    stop("meta/", basename(mpath), " has a different column set from this ",
         "script. Delete it and rerun the array rather than appending.")
  }
  man <- rbind(old[, names(man), drop = FALSE], man)
}
man <- man[order(match(man$pairing, ca_pairings$pairing),
                 match(man$outcome, ca_outcomes$outcome),
                 man$transform, man$stratum, man$ecoregion), , drop = FALSE]
utils::write.csv(man, mpath, row.names = FALSE)
unlink(lock, recursive = TRUE, force = TRUE)

ca_log("manifest ", nrow(man), " of ", nrow(RUNS), " runs recorded")

ca_stamp(sprintf("10b_etwfe_%s", R$run_id))
ca_log("Stage 10b task ", task, " ", status, ".  ", R$run_id)
