### CA carbon revision pipeline
### Stage 7. Spatial thinning of treated-role groups to 150 m
###
### WHAT THIS STAGE DOES, AND WHAT IT DOES NOT
###
### One operation. Reduce each group that appears on the treated side of a
### pairing to at most one pixel per 150 m block, subject to a 150 m minimum
### separation. No cap, no random draw, no covariate join, no carbon read.
###
### The 100,000 cap is a stage 8 operation, applied per pairing to the treated
### side after the pairing is known. Spacing must run first, because capping
### first and spacing second would leave roughly 4,000 units rather than
### 100,000.
###
### TREATED ROLE ONLY
###
### Spacing exists to stop adjacent 30 m pixels being counted as independent
### observations of the same stand. The submitted run applied it to treated
### pixels only, 8_2 filters the treated set by all.smpl.V1 and leaves the
### control set at full density, and this stage restores that specification.
###
### Consequence for UP, which is the only group holding both roles. UP is
### treated in UP_UNP and a control pool in FP_UP and TP_UP, so it exists in
### two forms. Spaced, in sample_<lulc>.parquet, read by stage 8 as the treated
### side of UP_UNP. Whole, in groups_<lulc>.parquet, read by stage 8 as the
### control pool for FP_UP and TP_UP. Nothing reconciles them because nothing
### needs to. UNP and OP never appear on a treated side and are not processed
### here at all.
###
### THE RETENTION RULE
###
### Block representative, not fixed lattice. A fixed-offset lattice keeps a
### pixel only where it is the centre of its own 150 m block, so a block that
### holds group pixels but not the centre pixel contributes nothing. That is
### harmless for a contiguous surface like UP and severe for a fragmented group
### like T1NP, where retention collapses for a reason that is an artifact of
### the rule rather than a property of the data.
###
### Instead, every occupied 150 m block contributes its best-centred pixel, and
### a second pass enforces the 150 m minimum between representatives sitting on
### either side of a shared block edge. Deterministic, no seed, and retention
### now scales with fragmentation, which is the correct behaviour. A group whose
### pixels are naturally more than 150 m apart is not spatially redundant and is
### left intact.
###
### EXCEPTIONS
###
### None by default. Spacing is uniform across all fourteen treated-role groups.
### ca_sample$whole_groups names any exemption explicitly and is empty, because
### the probe showed that neither a floor on group size nor a floor on treated
### units per exact cell measures anything that gates matching. Groups too small
### to estimate are handled at stage 8 by a gate on matched pairs.
###
### Probe mode remains, as the diagnostic behind that decision and as the source
### of the retention numbers the Methods reports.
###
### Usage
###   Rscript 7_sample_grid.R tasks      print the array size and exit
###   Rscript 7_sample_grid.R probe      per-group diagnostics, writes meta
###   Rscript 7_sample_grid.R            run one array task
###   Rscript 7_sample_grid.R archive    copy completed output to /projects

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))
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

BLK <- as.integer(ca_sample$block_pixels)      # 5, from 150 m over a 30 m cell
OFF <- as.integer(ca_sample$block_offset)      # 2, the centre of a 0 to 4 cycle
MIND2 <- as.integer(BLK * BLK)                 # squared 150 m in pixel units


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# pixel_id is the terra cell index of the frozen stage 1 grid, so grid position
# is integer arithmetic on it and no spatial package is involved. Arithmetic
# runs in double because the full California grid exceeds the integer range in
# the intermediate product, then returns to integer.
grid_rc <- function(pixel_id, ncol_grid) {
  z <- as.double(pixel_id) - 1
  list(row = as.integer(z %/% ncol_grid),
       col = as.integer(z %% ncol_grid))
}

# One representative per occupied 150 m block. Best centred wins, ties broken on
# ascending pixel_id, so the result is bit-identical on rerun.
block_representatives <- function(pixel_id, ncol_grid, bc_span) {

  rc <- grid_rc(pixel_id, ncol_grid)
  br <- rc$row %/% BLK
  bc <- rc$col %/% BLK
  dr <- (rc$row %% BLK) - OFF
  dc <- (rc$col %% BLK) - OFF
  d2 <- as.integer(dr * dr + dc * dc)
  bkey <- as.double(br) * bc_span + bc

  ord <- order(bkey, d2, pixel_id, method = "radix")
  first <- !duplicated(bkey[ord])
  keep <- ord[first]

  data.frame(pixel_id = pixel_id[keep],
             row = rc$row[keep],
             col = rc$col[keep],
             bkey = bkey[keep],
             d2 = d2[keep],
             stringsAsFactors = FALSE)
}

# Second pass. Representatives in diagonally or orthogonally adjacent blocks can
# sit closer than 150 m when both lie against the shared edge. Blocks two apart
# cannot conflict, since the minimum separation is then 2 * BLK - (BLK - 1)
# pixels, which exceeds BLK. So the conflict graph is confined to the eight
# surrounding blocks and the whole pass is local.
#
# Resolution is greedy in priority order, best centred first and ascending
# pixel_id on ties, executed as repeated rounds of a deterministic independent
# set rather than a sequential loop. Each round accepts every representative
# that is a local priority minimum among its undecided conflicting neighbours,
# which is exactly what the sequential greedy would have accepted, and then
# rejects their conflicting neighbours. The result is identical to the
# sequential version and runs vectorised.
enforce_spacing <- function(rp, bc_span) {

  n <- nrow(rp)
  if (n <= 1L) return(rep(TRUE, n))

  offs <- expand.grid(dbr = -1:1, dbc = -1:1)
  offs <- offs[!(offs$dbr == 0 & offs$dbc == 0), ]
  k_n <- nrow(offs)

  nb <- matrix(NA_integer_, n, k_n)
  cf <- matrix(FALSE, n, k_n)

  for (k in seq_len(k_n)) {
    j <- match(rp$bkey + offs$dbr[k] * bc_span + offs$dbc[k], rp$bkey)
    ok <- which(!is.na(j))
    if (!length(ok)) next
    jj <- j[ok]
    dr <- rp$row[jj] - rp$row[ok]
    dc <- rp$col[jj] - rp$col[ok]
    hit <- (dr * dr + dc * dc) < MIND2
    nb[ok, k] <- jj
    cf[ok[hit], k] <- TRUE
  }

  # Priority. Lower is better.
  prio <- integer(n)
  prio[order(rp$d2, rp$pixel_id, method = "radix")] <- seq_len(n)

  status <- integer(n)                       # 0 undecided, 1 keep, -1 drop
  round <- 0L

  repeat {
    if (!any(status == 0L)) break
    round <- round + 1L
    if (round > 100L) stop("Spacing pass failed to converge")

    blocked <- logical(n)
    for (k in seq_len(k_n)) {
      idx <- which(cf[, k])
      if (!length(idx)) next
      jj <- nb[idx, k]
      sel <- status[idx] == 0L & status[jj] == 0L & prio[jj] < prio[idx]
      if (any(sel)) blocked[idx[sel]] <- TRUE
    }

    acc <- status == 0L & !blocked
    if (!any(acc)) stop("Spacing pass stalled at round ", round)
    status[acc] <- 1L

    for (k in seq_len(k_n)) {
      idx <- which(cf[, k])
      if (!length(idx)) next
      jj <- nb[idx, k]
      sel <- status[idx] == 0L & status[jj] == 1L
      if (any(sel)) status[idx[sel]] <- -1L
    }
  }

  ca_log("    spacing pass converged in ", round, " rounds, dropped ",
         format(sum(status == -1L), big.mark = ","), " of ",
         format(n, big.mark = ","), " representatives")

  status == 1L
}

# Units per occupied exact cell. Matching is 1:1 without replacement inside
# cells defined by lulc, ecoregion, and HUC8, so this is the quantity that
# decides whether a group can be matched at all.
per_cell <- function(cell) {
  cell <- cell[!is.na(cell)]
  if (!length(cell)) return(c(n_cells = 0, median = NA_real_, q25 = NA_real_))
  tb <- tabulate(match(cell, unique(cell)))
  c(n_cells = length(tb),
    median = as.numeric(stats::median(tb)),
    q25 = as.numeric(stats::quantile(tb, 0.25, names = FALSE)))
}


# ---------------------------------------------------------------------------
# MODE tasks / archive
# ---------------------------------------------------------------------------

if (mode == "tasks") {
  cat("Array size: ", nrow(ca_lulc), "\n", sep = "")
  quit(save = "no")
}

if (mode == "archive") {
  ca_persist(ca_work("sample"), "sample", subdir = "work")
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# 1. CLASS, GROUPS, GEOMETRY
# ---------------------------------------------------------------------------

task <- ca_task_id()
if (mode == "probe" && is.na(task)) task <- 2L          # Evergreen by default
if (is.na(task) || task < 1L || task > nrow(ca_lulc)) {
  stop("Task id out of range: ", task)
}
lulc <- ca_lulc$label[task]

out_path <- ca_sample_path(lulc)
if (mode == "run" && file.exists(out_path) && !OVERWRITE) {
  ca_log("Output exists, nothing to do: ", basename(out_path))
  quit(save = "no")
}

ca_log("Stage 7, class ", lulc, ", mode ", mode)

# The group list is derived, never hand maintained. Anything on the treated
# side of a pairing is spaced, anything that only ever serves as a control pool
# is not. UNP and OP fall out of ca_pairings by construction.
treated_groups <- unique(unlist(strsplit(ca_pairings$treat, ",", fixed = TRUE)))
treated_groups <- trimws(treated_groups)
ca_log("Treated-role groups, ", length(treated_groups), ": ",
       paste(treated_groups, collapse = " "))

ncol_grid <- ca_grid_ncol()
bc_span <- as.double(ncol_grid %/% BLK + 1L)
ca_log("Grid ncol ", format(ncol_grid, big.mark = ","),
       ", block ", BLK, " pixels, ", ca_sample$grid_spacing_m, " m")


# ---------------------------------------------------------------------------
# 2. READ STAGE 6, KEEP TREATED-ROLE PIXELS
# ---------------------------------------------------------------------------

g <- as.data.frame(arrow::read_parquet(
  ca_group_path(lulc),
  col_select = dplyr::all_of(c("pixel_id", "analysis_group", "arm", "stratum",
                               "status_code", "cohort_year", "pa_est_year"))))
g <- g[!is.na(g$analysis_group) & g$analysis_group %in% treated_groups, ]
g <- g[order(g$pixel_id, method = "radix"), ]
ca_log("Treated-role pixels ", format(nrow(g), big.mark = ","))

missing_grp <- setdiff(treated_groups, unique(g$analysis_group))
if (length(missing_grp)) {
  ca_log("Groups absent in this class: ", paste(missing_grp, collapse = " "))
}

# Exact matching cell, for the per-cell diagnostic and the floor.
eco <- ca_read_static("ecoregion_l3", lulc)
g$ecoregion <- eco$value[match(g$pixel_id, eco$pixel_id)]
rm(eco); ca_gc()
huc <- ca_read_static("huc8", lulc)
g$huc8 <- huc$value[match(g$pixel_id, huc$pixel_id)]
rm(huc); ca_gc()
g$cell <- paste(g$ecoregion, g$huc8, sep = "|")


# ---------------------------------------------------------------------------
# 3. SPACE EACH GROUP
# ---------------------------------------------------------------------------

groups <- sort(unique(g$analysis_group))
stats <- vector("list", length(groups))
keep_ids <- vector("list", length(groups))

for (i in seq_along(groups)) {

  grp <- groups[i]
  sub <- g[g$analysis_group == grp, ]
  ca_log("  ", grp, "  n = ", format(nrow(sub), big.mark = ","))

  rep_i <- block_representatives(sub$pixel_id, ncol_grid, bc_span)
  ok <- enforce_spacing(rep_i, bc_span)
  kept <- rep_i$pixel_id[ok]
  keep_ids[[i]] <- kept

  cell_before <- per_cell(sub$cell)
  cell_after <- per_cell(sub$cell[match(kept, sub$pixel_id)])

  stats[[i]] <- data.frame(
    lulc = lulc,
    analysis_group = grp,
    arm = sub$arm[1],
    n_group = nrow(sub),
    n_blocks = nrow(rep_i),
    n_spaced = length(kept),
    retention_pct = round(100 * length(kept) / nrow(sub), 3),
    n_cells = cell_before[["n_cells"]],
    per_cell_before = cell_before[["median"]],
    per_cell_after = cell_after[["median"]],
    q25_cell_after = cell_after[["q25"]],
    cells_emptied = cell_before[["n_cells"]] - cell_after[["n_cells"]],
    stringsAsFactors = FALSE
  )

  rm(sub, rep_i); ca_gc()
}

stats <- do.call(rbind, stats)
names(keep_ids) <- groups


# ---------------------------------------------------------------------------
# 4. MODE probe. Write the diagnostic and stop
# ---------------------------------------------------------------------------
# Diagnostic only, run mode does not read it. It carries the retention and
# per-cell numbers the Methods and SI report, and it is what showed that a
# threshold-based floor could not be defended. Keep the three class files with
# the run.

if (mode == "probe") {
  p <- ca_meta(sprintf("sample_floor_%s.csv", lulc))
  utils::write.csv(stats, p, row.names = FALSE)
  ca_log("Floor diagnostic: ", p)
  print(stats)
  ca_log("Diagnostic written. Exemptions, if any, go in ca_sample$whole_groups")
  ca_stamp("7_sample_grid_probe")
  quit(save = "no")
}


# ---------------------------------------------------------------------------
# 5. EXEMPTIONS
# ---------------------------------------------------------------------------
# Spacing is uniform. ca_sample$whole_groups is empty by default and names any
# group exempted from spacing, which is the only mechanism for an exception.
# There is no threshold, because the probe showed that a floor on group size
# measures nothing relevant and a floor on treated units per exact cell
# measures the wrong side of the match. rule_applied rides on every row rather
# than living only in the summary, so stage 8 and the Methods both read the
# rule from the data.

stats$rule_applied <- ifelse(
  stats$analysis_group %in% ca_sample$whole_groups, "whole", "spaced")

whole <- stats$analysis_group[stats$rule_applied == "whole"]
if (length(whole)) {
  ca_log("Exempted from spacing, taken whole: ", paste(whole, collapse = " "))
} else {
  ca_log("Spacing applied to all ", nrow(stats), " treated-role groups")
}


# ---------------------------------------------------------------------------
# 6. WRITE
# ---------------------------------------------------------------------------

sel <- logical(nrow(g))
for (i in seq_along(groups)) {
  in_grp <- g$analysis_group == groups[i]
  if (stats$rule_applied[i] == "whole") {
    sel[in_grp] <- TRUE
  } else {
    sel[in_grp & g$pixel_id %in% keep_ids[[i]]] <- TRUE
  }
}

out <- data.frame(
  pixel_id = g$pixel_id[sel],
  lulc = lulc,
  analysis_group = g$analysis_group[sel],
  arm = g$arm[sel],
  stratum = g$stratum[sel],
  status_code = g$status_code[sel],
  cohort_year = g$cohort_year[sel],
  pa_est_year = g$pa_est_year[sel],
  rule_applied = stats$rule_applied[match(g$analysis_group[sel],
                                          stats$analysis_group)],
  stringsAsFactors = FALSE
)
out <- out[order(out$pixel_id, method = "radix"), ]

if (anyDuplicated(out$pixel_id)) stop("Duplicate pixel_id in stage 7 output")

arrow::write_parquet(arrow::as_arrow_table(out), out_path,
                     compression = ca_io$compression)
ca_log("Wrote ", basename(out_path), "  ",
       format(nrow(out), big.mark = ","), " rows  ",
       round(file.size(out_path) / 1e6, 1), " MB")

stats$n_written <- as.integer(table(factor(out$analysis_group,
                                           levels = stats$analysis_group)))
sp <- ca_meta(sprintf("sample_summary_%s.csv", lulc))
utils::write.csv(stats, sp, row.names = FALSE)
ca_log("Sample summary: ", sp)
print(stats[, c("analysis_group", "n_group", "n_blocks", "n_spaced",
                "retention_pct", "per_cell_before", "per_cell_after",
                "rule_applied", "n_written")])

ca_stamp("7_sample_grid")
ca_log("Done, class ", lulc)
