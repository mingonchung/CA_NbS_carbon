### CA carbon revision pipeline
### Stage 8b. Aggregate the matching sweep and recommend a specification
###
### WHY THIS EXISTS
###
### Stage 8 matches cohort by cohort, so it writes one row per stratification
### level. FNP_UNP baseline glm_nocal produces 36 rows in match_specs.csv.gz and
### 36 by 11 by 2 rows in match_balance.csv.gz. Table S3 needs one row per
### pairing and Figure S6 needs one point per covariate. This script closes that
### gap and nothing else.
###
### IT RECOMMENDS, IT DOES NOT DECIDE
###
### The ranking here is a reading aid. It prints a per-pairing table, marks the
### config the rule in ca_match$selection would pick, and emits an R snippet to
### paste into ca_config.R once the choice has been checked by eye. Nothing
### downstream reads the recommendation, only ca_match$selected, which is
### written by hand.
###
### ---------------------------------------------------------------------------
### FAILED IS NOT THE SAME AS ABSENT
### ---------------------------------------------------------------------------
###
### The previous version of this script dropped every spec row with a status
### other than ok, then excluded any config that did not carry the full set of
### levels. Those two rules together deleted whole configs from the ranking
### because glmnet cannot fit a propensity model to a cohort holding one treated
### pixel. UP_UNP augmented lost elastic_nocal, elastic_cal025, and
### rmahal_cal025 that way, all three of which ran to completion and wrote their
### matched parquet.
###
### Two kinds of missing level, and they call for opposite treatment.
###
###   failed   A row exists and its status is not ok. The config was attempted
###            on that level and could not be fitted or matched nothing. The run
###            is complete. The treated units in that level are lost to this
###            config and count against its match rate.
###
###   absent   No row exists. The task is still running, or the run was level
###            subset through CA_LEVELS and the other half has not landed. The
###            config is excluded from the ranking, because a rate built from
###            whichever cohorts finished first is not a rate.
###
### Only absence excludes. Failure is counted, reported, and carried into the
### denominator.
###
### ---------------------------------------------------------------------------
### THE MATCH RATE DENOMINATOR
### ---------------------------------------------------------------------------
###
### match_rate is matched pairs over MATCHABLE TREATED, and matchable treated is
### FIXED PER PAIRING AND ARM.
###
### n_treat is set from the complete-case frame before any config is fitted, is
### identical across configs, and is written on failed rows too. Summing the
### level-wise value over every level therefore gives a denominator that does not
### move when a config fails a cohort. The previous version summed n_treat over
### surviving rows only, which handed a failing config a smaller denominator and
### so a higher rate, the opposite of what the failure means.
###
### Upstream of that sit the spaced treated pool from stage 7 and the
### cap_per_group draw, neither of which stage 8 writes to the spec table. They
### are constants per pairing and arm, they are in the stage 8 logs, and stage 8
### is not rerun to capture them. The Methods state the cap. Nothing here
### imputes them.
###
### ---------------------------------------------------------------------------
### HOW BALANCE POOLS
### ---------------------------------------------------------------------------
###
### Means and variances pool exactly, because stage 8 writes n and sd on every
### balance row. The pooled variance is the within-level variance plus the
### between-level spread of the means, which is why sd alone would not do.
###
###   N   = sum(n_l)
###   mu  = sum(n_l * m_l) / N
###   var = [ sum((n_l - 1) * s_l^2) + sum(n_l * (m_l - mu)^2) ] / (N - 1)
###
### SMD and VR then follow from pooled moments and are exact. The denominator is
### the pooled unadjusted treated standard deviation, held fixed across stages so
### the before and after numbers sit on one scale, which is what MatchIt does
### within a level and what cobalt does with s.d.denom = "treated".
###
### The pooled SMD is the "as if one sample" statistic, which is what Table S3
### reports and what the submitted analysis reported. Per-level balance is not
### discarded, it stays in match_balance.csv.gz where stage 8 wrote it.
###
### ---------------------------------------------------------------------------
### THE eCDF STATISTICS
### ---------------------------------------------------------------------------
###
### An empirical distribution function cannot be reconstructed from a mean and a
### standard deviation, so the eCDF statistics do not pool exactly. They are
### aggregated across levels the same way everything else in the table is, as a
### treated-weighted average, and that is the only change.
###
### The previous version reported the eCDF maximum as a maximum across levels.
### For two samples of one unit each the Kolmogorov-Smirnov statistic is 1
### whenever the values differ, so any pairing containing a singleton cohort
### returned exactly 1.000 for every covariate, including covariates at SMD
### 0.003. UP_UNP, FNP_UNP, and FP_UP all hold cohorts of one treated pixel and
### all read 1.000 throughout. TNP_UNP, TP_UP, and ONP_UNP hold none and read
### normal values. The column was measuring the size of the smallest cohort.
###
###   ecdf_mean_wavg   treated-weighted mean of the within-level mean eCDF
###                    difference. The eCDF Mean column of Table S3.
###   ecdf_max_wavg    treated-weighted mean of the within-level eCDF maximum.
###                    The eCDF Max column of Table S3.
###   ecdf_max_worst   the raw maximum across levels, kept for reference with
###                    the level it came from and that level's pair count
###                    beside it, so a one-pair cohort reads as a one-pair
###                    cohort rather than as imbalance.
###
### ---------------------------------------------------------------------------
### THE UNADJUSTED SIDE
### ---------------------------------------------------------------------------
###
### Balance before matching is a property of the pairing, arm, and level, not of
### a config, so the eleven identical copies collapse to one and the pooled rows
### are written once with config = "(before)". Figure S6 joins its before points
### on pairing and arm. The previous version tried to expand them to every config
### and, because the deduplicated frame still carried a config column, wrote them
### with config NA instead, which broke that join silently.
###
### A control pixel matchable in several levels is counted once per level in the
### pooled unadjusted figures. That is the cohort-weighted before picture, which
### is the right comparator for a cohort-weighted after picture, and it is not a
### deduplicated census of the control pool. One SI sentence.
###
### Usage
###   Rscript 8b_aggregate_match.R            aggregate, rank, write, print
###   Rscript 8b_aggregate_match.R quiet      write without the console tables

rm(list = ls())

source(file.path(Sys.getenv("HOME"), "ca_config.R"))

ca_require(c("dplyr"))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
QUIET <- "quiet" %in% args

SEL <- ca_match$selection

# max and min over an empty or all-NA vector return -Inf and Inf with a warning.
# Both would propagate into the ranking as a finite-looking extreme.
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}
wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}


# ---------------------------------------------------------------------------
# READ
# ---------------------------------------------------------------------------

specs_path <- ca_match_specs_path()
bal_path   <- ca_match_balance_path()

for (p in c(specs_path, bal_path)) {
  if (!file.exists(p)) stop("Missing stage 8 output: ", p)
}

specs <- utils::read.csv(gzfile(specs_path), stringsAsFactors = FALSE)
bal   <- utils::read.csv(gzfile(bal_path),   stringsAsFactors = FALSE)

ca_log("specs rows ", format(nrow(specs), big.mark = ","),
       "  balance rows ", format(nrow(bal), big.mark = ","))

# A rerun that was not preceded by clearing the work directory appends rather
# than replaces, so identical rows would be double counted. Dropped here rather
# than trusted, since the append is by design and the duplication is not.
key_s <- c("pairing", "arm", "config", "win_key")
key_b <- c("pairing", "arm", "config", "win_key", "covariate", "stage")
dup_s <- sum(duplicated(specs[, key_s]))
dup_b <- sum(duplicated(bal[, key_b]))
if (dup_s || dup_b) {
  # LAST WRITE WINS. stream_append_gz opens the shared tables in append mode and
  # never replaces, so rerunning one config leaves the old rows in place above
  # the new ones on the same keys. Keeping the first occurrence would silently
  # return the run that was being replaced, which is the opposite of what a
  # rerun is for. fromLast is a backstop, not a substitute for stripping the
  # config out of the tables before resubmitting.
  ca_log("NOTE duplicate keys found, keeping the last occurrence. specs ",
         dup_s, ", balance ", dup_b, ". If these are not from a deliberate ",
         "rerun, clear the work directory and start again.")
  specs <- specs[!duplicated(specs[, key_s], fromLast = TRUE), , drop = FALSE]
  bal   <- bal[!duplicated(bal[, key_b], fromLast = TRUE), , drop = FALSE]
}

specs <- specs |>
  mutate(
    ok = status == "ok",
    fail_kind = case_when(
      status == "ok"       ~ NA_character_,
      status == "no units" ~ "no units",
      grepl("^failed:", status) ~ sub("^failed:\\s*", "", status),
      TRUE ~ status))


# ---------------------------------------------------------------------------
# THE LEVEL LEDGER. One row per pairing, arm, and level, config independent
# ---------------------------------------------------------------------------
#
# n_treat is set from the complete-case frame before any config is fitted, so it
# is a property of the level. Taking it from any row that carries the level, ok
# or failed, gives the fixed denominator that every config is scored against.

levels_ledger <- specs |>
  group_by(pairing, arm, win_key) |>
  summarise(
    cohort_year   = dplyr::first(cohort_year),
    n_treat_level = safe_max(n_treat),
    n_ctrl_level  = safe_max(n_control),
    n_configs_run = n_distinct(config),
    .groups = "drop")

# A level whose n_treat differs between configs would mean stage 8 filtered
# inside the config loop, which it does not. Checked rather than assumed.
inconsistent <- specs |>
  group_by(pairing, arm, win_key) |>
  summarise(n_distinct_treat = n_distinct(n_treat[is.finite(n_treat)]),
            .groups = "drop") |>
  filter(n_distinct_treat > 1L)
if (nrow(inconsistent)) {
  ca_log("WARNING n_treat is not constant across configs in ",
         nrow(inconsistent), " levels. The fixed denominator assumes it is.")
  print(as.data.frame(inconsistent), row.names = FALSE)
}

denom_tbl <- levels_ledger |>
  group_by(pairing, arm) |>
  summarise(
    n_levels_total     = n(),
    n_treat_matchable  = sum(n_treat_level, na.rm = TRUE),
    n_levels_singleton = sum(n_treat_level <= 2, na.rm = TRUE),
    .groups = "drop")

# ---------------------------------------------------------------------------
# COVERAGE. Failed levels are counted, absent levels exclude
# ---------------------------------------------------------------------------

coverage <- specs |>
  group_by(pairing, arm, config) |>
  summarise(
    n_levels_ok     = sum(ok),
    n_levels_failed = sum(!ok),
    n_treat_ok      = sum(n_treat[ok], na.rm = TRUE),
    n_treat_failed  = sum(n_treat[!ok], na.rm = TRUE),
    fail_reasons    = paste(sort(unique(fail_kind[!ok])), collapse = " | "),
    .groups = "drop") |>
  left_join(denom_tbl, by = c("pairing", "arm")) |>
  mutate(
    n_levels_present = n_levels_ok + n_levels_failed,
    n_levels_absent  = n_levels_total - n_levels_present,
    run_complete     = n_levels_absent == 0L,
    pct_treat_failed = ifelse(n_treat_matchable > 0,
                              100 * n_treat_failed / n_treat_matchable, NA_real_),
    coverage_note = case_when(
      !run_complete ~ sprintf("%d of %d levels absent, still running or level subset",
                              n_levels_absent, n_levels_total),
      n_levels_failed > 0 ~ sprintf("complete, %d levels failed carrying %.2f%% of treated units",
                                    n_levels_failed, pct_treat_failed),
      TRUE ~ "complete"))

n_incomplete <- sum(!coverage$run_complete)
n_withfail   <- sum(coverage$run_complete & coverage$n_levels_failed > 0)
ca_log("configs ranked ", sum(coverage$run_complete), " of ", nrow(coverage),
       ", of which ", n_withfail, " carry failed levels. ",
       n_incomplete, " excluded for absent levels.")

keep <- coverage |> filter(run_complete) |> select(pairing, arm, config)
if (!nrow(keep)) {
  utils::write.csv(coverage, ca_meta("match_coverage.csv"), row.names = FALSE)
  stop("No pairing, arm, and config combination has every level yet. ",
       "Coverage written to meta/match_coverage.csv.")
}


# ---------------------------------------------------------------------------
# POOL THE BALANCE TABLE
# ---------------------------------------------------------------------------

# Exact pooling of a mean and a variance across levels.
pool_moments <- function(n, m, s) {
  n <- as.numeric(n); m <- as.numeric(m); s <- as.numeric(s)
  ok <- is.finite(n) & n > 0 & is.finite(m)
  n <- n[ok]; m <- m[ok]; s <- s[ok]
  if (!length(n)) return(c(n = 0, mean = NA_real_, sd = NA_real_))
  s[!is.finite(s)] <- 0
  N  <- sum(n)
  mu <- sum(n * m) / N
  v  <- if (N > 1) {
    (sum((n - 1) * s^2) + sum(n * (m - mu)^2)) / (N - 1)
  } else NA_real_
  c(n = N, mean = mu, sd = sqrt(v))
}

bal <- bal |> semi_join(keep, by = c("pairing", "arm", "config"))

bal_ad <- bal |> filter(stage == "adjusted")

# The unadjusted comparison belongs to the level. The config column is dropped
# before deduplication so it cannot leak into the grouping or collide in a join.
bal_un <- bal |>
  filter(stage == "unadjusted") |>
  select(-config) |>
  distinct(pairing, arm, win_key, covariate, .keep_all = TRUE)

pool_stage <- function(d, by_config) {
  grp <- if (by_config) c("pairing", "arm", "config", "covariate") else
    c("pairing", "arm", "covariate")
  d |>
    group_by(across(all_of(grp))) |>
    summarise(
      n_levels       = n(),
      t              = list(pool_moments(n_treat, mean_treat, sd_treat)),
      c              = list(pool_moments(n_ctrl,  mean_ctrl,  sd_ctrl)),
      ecdf_mean_wavg = wmean(ecdf_mean, n_treat),
      ecdf_max_wavg  = wmean(ecdf_max,  n_treat),
      ecdf_max_worst = safe_max(ecdf_max),
      worst_level    = {
        i <- which.max(ifelse(is.finite(ecdf_max), ecdf_max, -Inf))
        if (length(i)) win_key[i] else NA_character_
      },
      worst_level_pairs = {
        i <- which.max(ifelse(is.finite(ecdf_max), ecdf_max, -Inf))
        if (length(i)) n_treat[i] else NA_real_
      },
      .groups = "drop") |>
    mutate(
      n_treat    = vapply(t, `[[`, numeric(1), "n"),
      n_ctrl     = vapply(c, `[[`, numeric(1), "n"),
      mean_treat = vapply(t, `[[`, numeric(1), "mean"),
      mean_ctrl  = vapply(c, `[[`, numeric(1), "mean"),
      sd_treat   = vapply(t, `[[`, numeric(1), "sd"),
      sd_ctrl    = vapply(c, `[[`, numeric(1), "sd")) |>
    select(-t, -c)
}

un <- pool_stage(bal_un, by_config = FALSE)
ad <- pool_stage(bal_ad, by_config = TRUE)

# The SMD denominator is the pooled unadjusted treated standard deviation, held
# fixed across stages so before and after sit on one scale. A property of the
# pairing, arm, and covariate.
denom <- un |> select(pairing, arm, covariate, denom = sd_treat)

finish <- function(d, stage_label) {
  d |>
    left_join(denom, by = c("pairing", "arm", "covariate")) |>
    mutate(
      stage = stage_label,
      smd = ifelse(is.finite(denom) & denom > 0,
                   (mean_treat - mean_ctrl) / denom, NA_real_),
      vr  = ifelse(is.finite(sd_ctrl) & sd_ctrl > 0 & is.finite(sd_treat),
                   (sd_treat / sd_ctrl)^2, NA_real_)) |>
    select(-denom)
}

cov_summary <- bind_rows(
  finish(un, "unadjusted") |> mutate(config = "(before)"),
  finish(ad, "adjusted")) |>
  arrange(pairing, arm, config, covariate, desc(stage)) |>
  select(pairing, arm, config, covariate, stage, n_levels,
         n_treat, n_ctrl, mean_treat, mean_ctrl, sd_treat, sd_ctrl,
         smd, vr, ecdf_mean_wavg, ecdf_max_wavg, ecdf_max_worst,
         worst_level, worst_level_pairs)


# ---------------------------------------------------------------------------
# CONFIG LEVEL SUMMARY
# ---------------------------------------------------------------------------
#
# The denominator is n_treat_matchable, fixed per pairing and arm. It does not
# move when a config fails a level, so rates across configs are comparable and a
# failure is a cost rather than a discount.

rates <- specs |>
  semi_join(keep, by = c("pairing", "arm", "config")) |>
  group_by(pairing, arm, config) |>
  summarise(
    n_matched_pairs = sum(n_matched[ok], na.rm = TRUE) / 2,
    .groups = "drop") |>
  left_join(coverage, by = c("pairing", "arm", "config")) |>
  mutate(
    match_rate       = ifelse(n_treat_matchable > 0,
                              n_matched_pairs / n_treat_matchable, NA_real_),
    n_treat_unmatched = n_treat_matchable - n_matched_pairs)

bal_cfg <- cov_summary |>
  filter(stage == "adjusted") |>
  group_by(pairing, arm, config) |>
  summarise(
    n_cov          = n(),
    max_abs_smd    = safe_max(abs(smd)),
    mean_abs_smd   = mean(abs(smd), na.rm = TRUE),
    n_smd_over     = sum(abs(smd) > SEL$smd_max, na.rm = TRUE),
    vr_min         = safe_min(vr),
    vr_max         = safe_max(vr),
    n_vr_outside   = sum(vr < SEL$vr_window[1] | vr > SEL$vr_window[2],
                         na.rm = TRUE),
    ecdf_mean_wavg = safe_max(ecdf_mean_wavg),
    ecdf_max_wavg  = safe_max(ecdf_max_wavg),
    ecdf_max_worst = safe_max(ecdf_max_worst),
    .groups = "drop")

cfg_summary <- rates |>
  left_join(bal_cfg, by = c("pairing", "arm", "config")) |>
  mutate(passes_balance = n_smd_over == 0 & n_vr_outside == 0)

# The structural ceiling. nn_exact_1to1 carries no caliper and no distance
# model, so its match rate is what the exact cells and the control pool allow
# and nothing else. Every other config's shortfall below it is caliper-induced,
# which is the only part that is a choice. Falls back to the best rate in the
# pairing if the ablation is missing.
ceilings <- cfg_summary |>
  group_by(pairing, arm) |>
  summarise(
    ceiling_rate = {
      z <- match_rate[config == SEL$ceiling_config]
      if (length(z) == 1L && is.finite(z)) z else safe_max(match_rate)
    },
    ceiling_from = if (any(config == SEL$ceiling_config)) SEL$ceiling_config
                   else "max_observed",
    .groups = "drop")

cfg_summary <- cfg_summary |>
  left_join(ceilings, by = c("pairing", "arm")) |>
  mutate(
    retention_ratio = ifelse(is.finite(ceiling_rate) & ceiling_rate > 0,
                             match_rate / ceiling_rate, NA_real_),
    meets_floor = is.finite(retention_ratio) &
                  retention_ratio >= SEL$retention_min,
    eligible    = meets_floor & passes_balance)

# Balance decides, retention gates.
cfg_summary <- cfg_summary |>
  group_by(pairing, arm) |>
  arrange(desc(eligible), desc(meets_floor),
          ifelse(eligible, mean_abs_smd, max_abs_smd),
          .by_group = TRUE) |>
  mutate(
    rank = row_number(),
    suggested = rank == 1L,
    note = case_when(
      !meets_floor ~ sprintf(
        "caliper kept %.3f of the structurally matchable units, below %.2f",
        retention_ratio, SEL$retention_min),
      meets_floor & !passes_balance ~ sprintf(
        "above the floor, %d covariates over SMD %.2f, %d outside the VR window",
        n_smd_over, SEL$smd_max, n_vr_outside),
      TRUE ~ "above the floor and passing")) |>
  ungroup() |>
  arrange(pairing, arm, rank)

# WHAT THE PICK COST IN RETENTION, AND WHAT IT BOUGHT.
best_ret <- cfg_summary |>
  filter(eligible) |>
  group_by(pairing, arm) |>
  slice_max(retention_ratio, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(pairing, arm,
         alt_config = config,
         alt_retention = retention_ratio,
         alt_mean_smd = mean_abs_smd,
         alt_max_smd = max_abs_smd)

cfg_summary <- cfg_summary |>
  left_join(best_ret, by = c("pairing", "arm")) |>
  mutate(
    retention_cost = ifelse(suggested & eligible,
                            alt_retention - retention_ratio, NA_real_),
    balance_gain   = ifelse(suggested & eligible,
                            alt_mean_smd - mean_abs_smd, NA_real_),
    trade_flag = ifelse(
      suggested & eligible & is.finite(retention_cost) &
        retention_cost > 0.02 & balance_gain < 0.005,
      sprintf("gives up %.3f retention against %s for %.4f mean SMD, check this",
              retention_cost, alt_config, balance_gain),
      ""))

any_pass <- cfg_summary |>
  group_by(pairing, arm) |>
  summarise(any_eligible = any(eligible), .groups = "drop")

cfg_summary <- cfg_summary |>
  left_join(any_pass, by = c("pairing", "arm")) |>
  mutate(note = ifelse(suggested & !any_eligible,
                       paste0(note, ". Nothing above the gate passes, ",
                              "this is the least bad and needs a human look"),
                       note))

# ---------------------------------------------------------------------------
# WRITE
# ---------------------------------------------------------------------------

utils::write.csv(coverage,  ca_meta("match_coverage.csv"),  row.names = FALSE)
utils::write.csv(cov_summary, file.path(ca_match_dir(),
                 "match_summary_covariate.csv"), row.names = FALSE)
utils::write.csv(cfg_summary, file.path(ca_match_dir(),
                 "match_summary_config.csv"), row.names = FALSE)

ca_log("wrote match_summary_covariate.csv, ",
       format(nrow(cov_summary), big.mark = ","), " rows")
ca_log("wrote match_summary_config.csv, ",
       format(nrow(cfg_summary), big.mark = ","), " rows")
ca_log("wrote meta/match_coverage.csv")


# ---------------------------------------------------------------------------
# PRINT
# ---------------------------------------------------------------------------

if (!QUIET) {

  cat("\n==== COVERAGE ====\n")
  cat("Failed levels are counted against the config and reported. Only absent\n")
  cat("levels exclude a config from the ranking.\n")
  print(as.data.frame(coverage |>
    filter(!run_complete | n_levels_failed > 0) |>
    transmute(pairing, arm, config,
              ok = n_levels_ok, failed = n_levels_failed,
              absent = n_levels_absent,
              pct_treat_lost = round(pct_treat_failed, 2),
              ranked = run_complete,
              why = substr(fail_reasons, 1, 58)) |>
    arrange(pairing, arm, config)), row.names = FALSE)

  cat("\n==== SWEEP, RANKED WITHIN PAIRING AND ARM ====\n")
  cat("Balance decides. Retention only has to clear ")
  cat(sprintf("%.2f of the structural ceiling,\n", SEL$retention_min))
  cat("which is the match rate of ")
  cat(sprintf("%s, the config with no caliper and no distance model.\n",
              SEL$ceiling_config))
  for (pr in unique(cfg_summary$pairing)) {
    for (ar in unique(cfg_summary$arm[cfg_summary$pairing == pr])) {
      d <- cfg_summary |> filter(pairing == pr, arm == ar)
      cat(sprintf("\n-- %s / %s   %d levels, %s matchable treated, ceiling %.3f\n",
                  pr, ar, safe_max(d$n_levels_total),
                  format(safe_max(d$n_treat_matchable), big.mark = ","),
                  safe_max(d$ceiling_rate)))
      print(as.data.frame(d |>
        transmute(config,
                  rate = round(match_rate, 4),
                  ret = round(retention_ratio, 3),
                  gate = meets_floor,
                  fail_lev = n_levels_failed,
                  max_smd = round(max_abs_smd, 3),
                  mean_smd = round(mean_abs_smd, 3),
                  smd_over = n_smd_over,
                  vr_out = n_vr_outside,
                  ecdf = round(ecdf_mean_wavg, 3),
                  ecdf_max = round(ecdf_max_wavg, 3),
                  pick = ifelse(suggested, "  <==", ""))),
        row.names = FALSE)
    }
  }

  cat("\n==== eCDF MAXIMUM, WHERE IT COMES FROM ====\n")
  cat("A cohort of one matched pair returns an eCDF maximum of 1 by\n")
  cat("construction. ecdf_max_wavg is the treated-weighted average across\n")
  cat("cohorts and is the Table S3 column. The raw worst is printed with its\n")
  cat("cohort and pair count so that distinction is visible.\n")
  print(as.data.frame(cov_summary |>
    filter(stage == "adjusted",
           config %in% cfg_summary$config[cfg_summary$suggested]) |>
    group_by(pairing, arm) |>
    slice_max(ecdf_max_worst, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(pairing, arm, covariate,
              ecdf_mean = round(ecdf_mean_wavg, 3),
              ecdf_max = round(ecdf_max_wavg, 3),
              raw_worst = round(ecdf_max_worst, 3),
              at_level = worst_level,
              pairs = worst_level_pairs) |>
    arrange(pairing, arm)), row.names = FALSE)

  cat("\n==== ARM COMPARISON, LIKE FOR LIKE ====\n")
  cat("Arms cannot be compared on their own max SMD, because each is taken over\n")
  cat("a different covariate count. These columns restrict every arm to the\n")
  cat("covariates all arms share.\n")

  shared_cov <- cov_summary |>
    filter(stage == "adjusted") |>
    distinct(arm, covariate) |>
    group_by(covariate) |>
    summarise(n_arms = n(), .groups = "drop") |>
    filter(n_arms == n_distinct(cov_summary$arm)) |>
    pull(covariate)

  if (length(shared_cov)) {
    cat(sprintf("\nShared covariates (%d): %s\n", length(shared_cov),
                paste(shared_cov, collapse = ", ")))
    common <- cov_summary |>
      filter(stage == "adjusted", covariate %in% shared_cov) |>
      group_by(pairing, arm, config) |>
      summarise(max_smd_shared  = safe_max(abs(smd)),
                mean_smd_shared = mean(abs(smd), na.rm = TRUE),
                .groups = "drop")

    added <- cov_summary |>
      filter(stage == "adjusted", !covariate %in% shared_cov) |>
      group_by(pairing, arm, config) |>
      summarise(added_cov = paste(covariate, collapse = ", "),
                max_smd_added = safe_max(abs(smd)),
                .groups = "drop")

    print(as.data.frame(common |>
      left_join(added, by = c("pairing", "arm", "config")) |>
      semi_join(cfg_summary |> filter(suggested) |>
                  select(pairing, arm, config),
                by = c("pairing", "arm", "config")) |>
      transmute(pairing, arm, config,
                max_shared  = round(max_smd_shared, 3),
                mean_shared = round(mean_smd_shared, 3),
                added_cov   = ifelse(is.na(added_cov), "", added_cov),
                max_added   = round(max_smd_added, 3)) |>
      arrange(pairing, arm)), row.names = FALSE)
  }

  flagged <- cfg_summary |> filter(nzchar(trade_flag))
  if (nrow(flagged)) {
    cat("\n==== TRADES WORTH A SECOND LOOK ====\n")
    print(as.data.frame(flagged |>
      transmute(pairing, arm, suggested = config,
                ret = round(retention_ratio, 3),
                mean_smd = round(mean_abs_smd, 4),
                alternative = alt_config,
                alt_ret = round(alt_retention, 3),
                alt_mean_smd = round(alt_mean_smd, 4),
                alt_max_smd = round(alt_max_smd, 3))), row.names = FALSE)
  }

  cat("\n==== SUGGESTED CONFIG IN EACH ARM ====\n")
  print(as.data.frame(cfg_summary |>
    filter(suggested) |>
    transmute(pairing, arm, config,
              rate = round(match_rate, 4),
              max_smd = round(max_abs_smd, 3),
              mean_smd = round(mean_abs_smd, 3),
              ecdf = round(ecdf_mean_wavg, 3)) |>
    arrange(pairing, arm)), row.names = FALSE)

  cat("\n==== SNIPPET FOR ca_config.R, AFTER YOU HAVE CHECKED IT ====\n")
  cat("One block per arm. Copy the block for the arm you are making primary.\n")
  sel_rows <- cfg_summary |> filter(suggested) |> arrange(arm, pairing)
  for (ar in unique(sel_rows$arm)) {
    d <- sel_rows |> filter(arm == ar)
    cat(sprintf("\n  # arm: %s\n  selected = list(\n", ar))
    for (i in seq_len(nrow(d))) {
      cat(sprintf('    %-9s = list(arm = "%s", config = "%s")%s   # rate %.3f, max SMD %.3f\n',
                  d$pairing[i], d$arm[i], d$config[i],
                  if (i < nrow(d)) "," else "",
                  d$match_rate[i], d$max_abs_smd[i]))
    }
    cat("  ),\n")
  }
  cat("\nNothing downstream reads the ranking, only what is written by hand.\n")
}

ca_stamp("8b_aggregate_match")
ca_log("Stage 8b complete.")
