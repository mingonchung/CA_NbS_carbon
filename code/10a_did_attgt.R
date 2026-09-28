### CA carbon revision pipeline
### Stage 10a. att_gt and aggte, one file per cell, several cells per task
###
### WHAT THIS STAGE DOES
###
### Callaway and Sant'Anna group-time ATTs, aggregated to an overall ATT, to
### dynamic event time, and on the PA arm to calendar time. This is the SI
### comparison and the parallel trends evidence. The headline estimate is
### etwfe, stage 10b.
###
### ONE TASK, SEVERAL CELLS
###
### An array task is a statewide cell. It reads its panel once and then fits
### that cell and every ecoregion inside it, writing one output file and one
### manifest row per cell. ca_did_runs("attgt") carries every cell and a task
### column saying which task computes it, ca_did_array("attgt") is the task
### table, and the array is sized on the latter.
###
### The submitted pipeline ran att_gt statewide only, 10_1 L112, 10_5 L73,
### 10_2 L85, 10_3 L72, and this script followed it until 2026-08-16. The
### reason given was the cost of the multiplier bootstrap. The heaviest cell
### fits in 15 to 20 minutes inside 100 GB, so the cost was never the
### constraint and the ecoregion comparison was withheld for nothing. What is
### kept from that reasoning is the array shape. One task per ecoregion would
### spend more scheduler time than estimator time, and would re-read the panel
### once per ecoregion the way 10b has to.
###
### A CELL WITH OUTPUT IS SKIPPED, NOT RECOMPUTED
###
### The loop checks each cell's own file, so a task that died on its fourth
### ecoregion resumes at the fourth rather than refitting the first three.
### CA_OVERWRITE=1 turns that off. This is why the output is one file per cell
### rather than one file per task.
###
### WHAT LEAVES THE SUBMITTED SPECIFICATION, AND WHY
###
### 1. clustervars = "pair_id". The submitted call left it unset, which
###    clusters on idname, the unit. ca_did$cluster is the matched pair, and
###    etwfe already uses vcov = ~ pair_id, so the two estimators now report
###    standard errors on the same clustering.
###
### 2. The panel is balanced here rather than inside att_gt(). The submitted
###    default allow_unbalanced_panel = FALSE is kept, so did would balance it
###    anyway and report the cost in a warning. Doing it first puts the number
###    in the manifest.
###
### 3. gvar_est = gvar - 1 on the thinning pairings, so the base period is c-2.
###    ca_config.R section 0.11c.
###
### na.rm = TRUE, as submitted at 10_1 L144 and L153, and required rather than
### optional, since aggte() returns NA for the whole aggregation when any
### group-time cell it averages is NA. The count of NA cells goes to the
### manifest, which is what the submitted pipeline left silent.
###
### Usage
###   Rscript 10a_did_attgt.R tasks     print the array and tier ranges
###   Rscript 10a_did_attgt.R           run one array task
###   Rscript 10a_did_attgt.R missing   print an sbatch array of unfinished runs
###   Rscript 10a_did_attgt.R report    wall time against cost, per tier
###   Rscript 10a_did_attgt.R archive   copy completed output to /projects
###
###   meta/did_subset_levels.csv must exist. Build it once with
###   Rscript 10b_did_etwfe.R levels
###
###   CA_OVERWRITE=1                    rerun cells that already have output

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
source(file.path(Sys.getenv("HOME"), "ca_clean.R"))

ca_require(c("arrow", "dplyr", "did"))

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) && args[1] %in% c("tasks", "archive", "missing", "report")) args[1] else "run"
OVERWRITE <- "overwrite" %in% args || nzchar(Sys.getenv("CA_OVERWRITE"))

ESTIMATOR <- "attgt"

if (mode == "archive") {
  ca_persist(ca_did_est_dir(), "model", subdir = "work")
  quit(save = "no")
}

RUNS  <- ca_did_runs(ESTIMATOR)
ARRAY <- ca_did_array(ESTIMATOR)

# Cells per task, for the log line and the tasks report.
CELLS_PER_TASK <- as.integer(table(factor(RUNS$task,
                                          levels = seq_len(nrow(ARRAY)))))

if (mode == "tasks") {
  show <- ARRAY[, c("tier", "pairing", "outcome", "transform", "stratum",
                    "n_treated", "n_cohorts", "est_rows", "est_cost")]
  show$cells <- CELLS_PER_TASK
  print(show, right = FALSE)
  cat("\narray size:", nrow(ARRAY), " tasks,", nrow(RUNS), "cells\n")
  for (tr in c("small", "medium", "large")) {
    i <- which(ARRAY$tier == tr)
    if (!length(i)) next
    cat(sprintf("%-7s %d-%d  (%d tasks, %d cells, cost %s to %s)\n", tr,
                min(i), max(i), length(i), sum(CELLS_PER_TASK[i]),
                format(min(ARRAY$est_cost[i]), big.mark = ","),
                format(max(ARRAY$est_cost[i]), big.mark = ",")))
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

JOB_NAME <- "ca_10a_attgt"

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

  # Classification is per cell. Submission is per task, because a task is what
  # the scheduler takes. A task with one unfinished ecoregion is resubmitted
  # whole and skips the cells that already have output.
  cls <- rep("todo", nrow(RUNS))
  cls[!valid & !is.na(st) & st == "skipped"] <- "skipped"
  cls[!valid & !is.na(st) & st == "error" & !res] <- "error"
  cls[RUNS$task %in% running & !valid] <- "running"
  cls[valid] <- "done"

  truncated <- which(present & !valid)

  cat(nrow(RUNS), "cells in", nrow(ARRAY), "tasks\n")
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
  cat("\n", length(idx), " cells to compute in ",
      length(unique(RUNS$task[idx])), " tasks",
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
  # Tier comes from the task, not the cell. A task is sized by its statewide
  # cell, which is the largest cell it holds.
  todo_task <- sort(unique(RUNS$task[idx]))
  for (tr in c("small", "medium", "large")) {
    v <- todo_task[ARRAY$tier[todo_task] == tr]
    if (!length(v)) next
    cat("\n", tr, ", ", length(v), " tasks, ",
        sum(cls[RUNS$task %in% v] == "todo"), " cells, cost ",
        format(min(ARRAY$est_cost[v]), big.mark = ","), " to ",
        format(max(ARRAY$est_cost[v]), big.mark = ","), "\n", sep = "")
    cat("--array=", arr(v), "\n", sep = "")
  }
  quit(save = "no")
}

task <- ca_task_id()
if (is.na(task) || task < 1L || task > nrow(ARRAY)) {
  stop("Task id must be 1 to ", nrow(ARRAY), ", got ", task)
}

RT    <- ARRAY[task, ]
CELLS <- RUNS[RUNS$task == task, , drop = FALSE]

# Statewide first, then the ecoregions in code order. The statewide cell is the
# one every figure needs, so it is computed before anything that could run the
# clock out.
CELLS <- CELLS[order(CELLS$ecoregion != "", CELLS$ecoregion), , drop = FALSE]
rownames(CELLS) <- NULL

t_start <- Sys.time()

ca_log("Stage 10a, task ", task, ".  ", RT$run_id)
ca_log("  tier ", RT$tier, "  est rows ", format(RT$est_rows, big.mark = ","),
       "  cost ", format(RT$est_cost, big.mark = ","),
       "  base shift ", RT$base_shift)
ca_log("  ", nrow(CELLS), " cells: statewide",
       if (nrow(CELLS) > 1L) paste0(" plus ", nrow(CELLS) - 1L,
                                    " ecoregions (",
                                    paste(CELLS$ecoregion[-1], collapse = ", "),
                                    ")") else "")

CELLS$out_path <- vapply(seq_len(nrow(CELLS)), function(i) {
  ca_did_est_path(ESTIMATOR, CELLS$pairing[i], CELLS$run_id[i])
}, character(1))

# NOTHING LEFT TO DO IS NOT A FAILURE.
#
# A resubmitted task whose cells all have output exits before the panel read,
# which is the whole point of resuming at cell level rather than task level.
todo <- OVERWRITE | !file.exists(CELLS$out_path)
if (!any(todo)) {
  ca_log("All ", nrow(CELLS), " cells present, not overwritten.")
  quit(save = "no")
}
if (any(!todo)) {
  ca_log("  ", sum(!todo), " cells already present, skipping: ",
         paste(ifelse(nzchar(CELLS$ecoregion[!todo]),
                      CELLS$ecoregion[!todo], "statewide"), collapse = ", "))
}


# ---------------------------------------------------------------------------
# READ THE TASK PANEL
# ---------------------------------------------------------------------------
# Once, for the whole task. ecoregion_l3 is kept so the loop can subset in
# memory rather than rescanning the dataset per ecoregion, which is the saving
# that makes the loop worth having.

PANEL_COLS <- c("unit_id", "pair_id", "treat", "gvar", "year", "value",
                "stratum", "ecoregion_l3")

read_task_panel <- function() {
  ds <- arrow::open_dataset(RT$path)
  q  <- dplyr::select(ds, dplyr::all_of(PANEL_COLS))
  st <- RT$stratum
  pushed <- tryCatch({
    if (nzchar(st)) q <- dplyr::filter(q, stratum == st)
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
  d
}

D_TASK <- read_task_panel()
D_TASK$stratum <- NULL
D_TASK$ecoregion_l3 <- as.character(D_TASK$ecoregion_l3)
ca_gc("task panel read")

ca_log("  panel ", format(nrow(D_TASK), big.mark = ","), " rows, ",
       length(unique(D_TASK$ecoregion_l3)), " ecoregions present")


# ---------------------------------------------------------------------------
# ONE CELL
# ---------------------------------------------------------------------------
# Everything from the cohort filter to the manifest row, for one cell. The
# cohort filter, the balance step, and the never-treated check are all
# properties of the subset, so none of them can be done once on the statewide
# panel and reused by an ecoregion.

tidy_aggte <- function(a, aggregation) {
  egt <- as.numeric(a$egt)
  cv  <- a$crit.val.egt
  if (is.null(cv) || !is.finite(cv)) cv <- stats::qnorm(0.975)
  out <- data.frame(
    aggregation = aggregation,
    term        = c("overall", as.character(a$egt)),
    egt         = c(NA_real_, egt),
    estimate    = c(as.numeric(a$overall.att), as.numeric(a$att.egt)),
    std.error   = c(as.numeric(a$overall.se), as.numeric(a$se.egt)),
    crit_val    = c(stats::qnorm(0.975), rep(cv, length(egt))),
    stringsAsFactors = FALSE)
  out$conf.low  <- out$estimate - out$crit_val * out$std.error
  out$conf.high <- out$estimate + out$crit_val * out$std.error
  out$statistic <- out$estimate / out$std.error
  out$p.value   <- 2 * stats::pnorm(-abs(out$statistic))
  out
}

run_cell <- function(R) {

  t_cell <- Sys.time()
  lab <- if (nzchar(R$ecoregion)) R$ecoregion else "statewide"
  ca_log("Cell ", lab, "  ", R$run_id)

  status  <- "ok"
  message <- ""
  fail <- function(msg) {
    status  <<- "skipped"
    message <<- msg
    ca_log("  SKIP ", msg)
  }

  d <- if (nzchar(R$ecoregion)) {
    D_TASK[D_TASK$ecoregion_l3 == R$ecoregion, , drop = FALSE]
  } else D_TASK
  d$ecoregion_l3 <- NULL

  if (!nrow(d)) fail("no rows in this cell")

  n_cohorts_dropped <- 0L
  n_always_treated  <- 0L
  n_unbalanced_drop <- 0L
  yr_first <- NA_integer_
  yr_last  <- NA_integer_

  if (status == "ok") {

    d$year     <- as.integer(d$year)
    d$unit_id  <- as.integer(d$unit_id)
    d$gvar_est <- ca_did_gvar_est(R$pairing, d$gvar)

    yr_first <- min(d$year)
    yr_last  <- max(d$year)

    coh <- sort(unique(d$gvar_est[d$treat == 1L]))
    bad <- coh[coh > yr_last]
    n_cohorts_dropped <- length(bad)
    if (n_cohorts_dropped) {
      drop_id <- unique(d$unit_id[d$treat == 1L & d$gvar_est %in% bad])
      ca_log("  cohorts after ", yr_last, " dropped: ", n_cohorts_dropped,
             " (", format(length(drop_id), big.mark = ","), " treated units)")
      d <- d[!(d$unit_id %in% drop_id), , drop = FALSE]
    }

    # BALANCE.
    #
    # att_gt() with allow_unbalanced_panel = FALSE converts to a balanced panel
    # internally. Stage 9 dropped rows whose outcome was NA, so a unit missing
    # any year would be removed there without a number reaching the manifest.
    # Done here it is one count. LEMMA's static mask is the case that matters,
    # since it removes whole pixels rather than scattered years.
    if (isTRUE(ca_did_est$attgt_balance)) {
      ny <- length(unique(d$year))
      cnt <- table(d$unit_id)
      keep <- as.integer(names(cnt)[cnt == ny])
      n_unbalanced_drop <- length(cnt) - length(keep)
      if (n_unbalanced_drop) {
        ca_log("  unbalanced units dropped: ",
               format(n_unbalanced_drop, big.mark = ","), " of ",
               format(length(cnt), big.mark = ","))
        d <- d[d$unit_id %in% keep, , drop = FALSE]
      }
    }

    n_always_treated <- sum(unique(d$gvar_est[d$treat == 1L]) == yr_first)

    if (!nrow(d)) fail("no rows after balancing")
    if (status == "ok" && !sum(d$treat == 1L)) {
      fail("no treated units after the cohort filter")
    }
    if (status == "ok" && !any(d$gvar_est == 0L)) {
      fail("no never-treated units in this cell")
    }
  }

  # THE FLOOR, APPLIED AGAIN HERE.
  #
  # ca_did_runs() applies it to the cell counts stage 9 recorded. Those counts
  # are before the cohort filter and the balance step, both of which can take a
  # thin ecoregion below the floor. A cell that arrives at the fit with too few
  # treated units is a skip with a reason, not an estimate with a wide band.
  if (status == "ok" && ca_did_est$min_treated_units > 0L) {
    n_tr <- length(unique(d$unit_id[d$treat == 1L]))
    if (n_tr < ca_did_est$min_treated_units) {
      fail(paste0("treated units after balancing ", n_tr, " below floor ",
                  ca_did_est$min_treated_units))
    }
  }

  est        <- NULL
  wald_p     <- NA_real_
  n_gt       <- NA_integer_
  n_gt_na    <- NA_integer_
  n_pre      <- NA_integer_
  n_pre_sig  <- NA_integer_
  max_pre    <- NA_real_
  att_simple <- NA_real_
  se_simple  <- NA_real_
  n_rows  <- if (status == "ok") nrow(d) else 0L
  n_units <- if (status == "ok") length(unique(d$unit_id)) else 0L
  n_pairs <- if (status == "ok") length(unique(d$pair_id)) else 0L
  n_treat <- if (status == "ok") length(unique(d$unit_id[d$treat == 1L])) else 0L
  n_coh   <- if (status == "ok") length(unique(d$gvar_est[d$treat == 1L])) else 0L
  gv_min  <- NA_integer_
  gv_max  <- NA_integer_

  if (status == "ok") {

    d$y <- ca_did_apply_transform(d$value, R$transform)
    d$value <- NULL
    d$gvar <- NULL

    g <- d$gvar_est[d$gvar_est > 0L]
    gv_min <- min(g); gv_max <- max(g)

    ca_log("  rows ", format(n_rows, big.mark = ","),
           "  units ", format(n_units, big.mark = ","),
           "  pairs ", format(n_pairs, big.mark = ","),
           "  treated ", format(n_treat, big.mark = ","),
           "  cohorts ", n_coh, " (", gv_min, " to ", gv_max, ")")

    set.seed(ca_did_est$seed)

    fit <- tryCatch({
      mp <- did::att_gt(
        yname                  = "y",
        tname                  = "year",
        idname                 = "unit_id",
        gname                  = "gvar_est",
        xformla                = ca_did_est$xformla,
        data                   = d,
        control_group          = ca_did_est$control_group,
        clustervars            = ca_did$cluster,
        allow_unbalanced_panel = ca_did_est$attgt_allow_unbalanced,
        base_period            = ca_did_est$attgt_base_period,
        anticipation           = 0,
        bstrap                 = ca_did_est$attgt_bstrap,
        cband                  = ca_did_est$attgt_cband)
      ca_log("  att_gt fitted, ", length(mp$att), " group-time cells")

      wald_p  <<- if (is.null(mp$Wpval)) NA_real_ else as.numeric(mp$Wpval)[1]
      n_gt    <<- length(mp$att)
      n_gt_na <<- sum(is.na(mp$att))
      ca_log("  parallel trends Wald p = ", signif(wald_p, 4),
             ", NA group-time cells removed by na.rm = ", n_gt_na,
             " of ", n_gt)

      # AN AGGREGATION THAT FAILS DOES NOT TAKE THE CELL WITH IT.
      #
      # aggte() can error on a thin cell, most often calendar or dynamic where
      # a year holds no usable group-time comparison. Losing the overall ATT
      # because one event time could not be formed would throw away a fit that
      # costs minutes to rebuild, so each type is wrapped and a failure is
      # logged and dropped.
      parts <- lapply(ca_did_aggte_types(R$pairing), function(ty) {
        a <- tryCatch(did::aggte(mp, type = ty,
                                 na.rm = ca_did_est$aggte_na_rm),
                      error = function(e) {
                        ca_log("  aggte ", ty, " FAILED: ",
                               conditionMessage(e))
                        NULL
                      })
        if (is.null(a)) return(NULL)
        if (identical(ty, "simple")) {
          att_simple <<- as.numeric(a$overall.att)
          se_simple  <<- as.numeric(a$overall.se)
        }
        ca_log("  aggte ", ty, " done")
        tidy_aggte(a, ty)
      })
      parts <- parts[!vapply(parts, is.null, logical(1))]
      if (!length(parts)) stop("every aggregation failed")
      do.call(rbind, parts)
    }, error = function(e) {
      status  <<- "error"
      message <<- conditionMessage(e)
      ca_log("  ERROR ", conditionMessage(e))
      NULL
    })

    if (!is.null(fit)) {
      est <- as.data.frame(fit)

      # THE PRE-TREND EVIDENCE, SINCE THE WALD PRE-TEST IS NOT AVAILABLE.
      #
      # did declines to report Wpval when clustervars goes beyond the unit,
      # because its analytical variance matrix does not carry between-cluster
      # correlation. The choice is a Wald statistic computed on an understated
      # variance or the bootstrap simultaneous bands computed on the right one.
      # The bands are also what the submitted figures plotted, 10_1 L163, so
      # the pre-trend claim is summarised from the dynamic aggregation at
      # negative event time and wald_pval stays in the manifest as NA to record
      # why.
      pre <- est[est$aggregation == "dynamic" & !is.na(est$egt) &
                   est$egt < 0, , drop = FALSE]
      n_pre     <- nrow(pre)
      n_pre_sig <- sum(pre$conf.low > 0 | pre$conf.high < 0, na.rm = TRUE)
      max_pre   <- if (n_pre) max(abs(pre$estimate), na.rm = TRUE) else NA_real_
      ca_log("  pre-treatment event times ", n_pre, ", simultaneous band ",
             "excludes zero in ", n_pre_sig, ", largest |estimate| ",
             signif(max_pre, 4))
    }
  }

  rm(d)
  ca_gc(paste0("cell ", lab, " done"))

  out_path <- R$out_path

  if (!is.null(est) && nrow(est)) {
    est$run_id     <- R$run_id
    est$estimator  <- ESTIMATOR
    est$pairing    <- R$pairing
    est$outcome    <- R$outcome
    est$transform  <- R$transform
    est$stratum    <- R$stratum
    est$ecoregion  <- R$ecoregion
    est$base_shift <- R$base_shift
    est$wald_pval  <- wald_p

    con <- gzfile(out_path, "w")
    utils::write.csv(est, con, row.names = FALSE)
    close(con)
    ca_log("  wrote ", basename(out_path), "  ", nrow(est), " rows")
  } else {
    out_path <- NA_character_
  }

  sel <- ca_match_selected(R$pairing)

  data.frame(
    run_id            = R$run_id,
    estimator         = ESTIMATOR,
    task              = task,
    pairing           = R$pairing,
    outcome           = R$outcome,
    transform         = R$transform,
    stratum           = R$stratum,
    ecoregion         = R$ecoregion,
    arm               = sel$arm,
    config            = sel$config,
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
    n_unbalanced_drop = n_unbalanced_drop,
    n_gt_cells        = n_gt,
    n_gt_na           = n_gt_na,
    wald_pval         = wald_p,
    n_pre             = n_pre,
    n_pre_sig         = n_pre_sig,
    max_abs_pre       = max_pre,
    att_simple        = att_simple,
    se_simple         = se_simple,
    gvar_min          = gv_min,
    gvar_max          = gv_max,
    wall_min          = round(as.numeric(difftime(Sys.time(), t_cell,
                                                  units = "mins")), 2),
    path              = out_path,
    written_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    message           = message,
    stringsAsFactors  = FALSE)
}


# ---------------------------------------------------------------------------
# THE LOOP
# ---------------------------------------------------------------------------

MAN <- vector("list", nrow(CELLS))
for (i in which(todo)) {
  MAN[[i]] <- tryCatch(run_cell(CELLS[i, ]), error = function(e) {
    # A cell that fails outside the fit still gets a manifest row, or the task
    # would look complete to 10c while one of its cells has neither output nor
    # a record of why.
    ca_log("  CELL ERROR ", conditionMessage(e))
    data.frame(run_id = CELLS$run_id[i], estimator = ESTIMATOR, task = task,
               pairing = CELLS$pairing[i], status = "error",
               ecoregion = CELLS$ecoregion[i],
               message = conditionMessage(e), stringsAsFactors = FALSE)
  })
}
MAN <- MAN[!vapply(MAN, is.null, logical(1))]


# ---------------------------------------------------------------------------
# MANIFEST
# ---------------------------------------------------------------------------
# One lock per task, not one per cell. Seven cells appending in turn would take
# the lock seven times and rewrite the whole manifest each time.

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

if (length(MAN)) {

  man <- do.call(rbind, lapply(MAN, function(x) {
    data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
  }))

  lock <- file.path(ca_did_est_dir(), ".lock_manifest_attgt")
  acquire_lock(lock)
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  # READ THE MANIFEST BACK AS CHARACTER.
  #
  # read.csv() with default colClasses reparses every column of every row
  # already written. An ecoregion of "08" returns as the number 8 and an empty
  # stratum returns as NA, so the rows written by earlier tasks silently lose
  # their zero padding and their empty strings. With hundreds of tasks
  # appending in turn the damage compounds. Character in, character out, and
  # the numbers are numbers again wherever they are read for analysis.
  mpath <- ca_did_est_manifest_path(ESTIMATOR)
  if (file.exists(mpath)) {
    old <- utils::read.csv(mpath, colClasses = "character")
    old <- old[!(old$run_id %in% man$run_id), , drop = FALSE]
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

  ca_log("manifest ", nrow(man), " of ", nrow(RUNS), " cells recorded")
}

ca_stamp(sprintf("10a_attgt_task_%s", task))
ca_log("Stage 10a task ", task, " done, ", sum(todo), " cells computed in ",
       round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
       " min.  ", RT$run_id)
