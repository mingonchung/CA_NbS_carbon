# CA_NbS_carbon

Code and derived data for Chung et al., "Promoting long-term carbon durability in fire-prone forests requires reducing stocks", *Global Change Biology* (in revision).

Statistical matching and difference-in-differences estimates of carbon additionality and durability for protected areas (PAs), forest offsets, forest thinning, and wildfire in California forests, 1990–2025, at 30 m resolution.

Version 2.0.0 reproduces the revised manuscript. Version 1.0.0 (tag `v.1.0.0`, DOI 10.5281/zenodo.17665766) holds the code for the original submission.

## Changes from version 1.0.0

- AGB and GPP from the Wildland Almanac (1985–2025), replacing NCSDA total live biomass and GPP.
- Tree NPP and NEP from NCSDA (1985–2021) added as outcomes.
- Pre-treatment AGB, GPP, and canopy height added as matching covariates.
- FACTS and CAL FIRE thinning records and CARB offsets updated through 2025.
- Wildfire and thinning effects estimated against matched undisturbed forests of the same protection status.
- Ecoregion-level estimates added for thinning.
- Full pipeline from data extraction to estimation, replacing the matching and DiD scripts of version 1.0.0.

## Repository layout

```
code/   stages 1–10, shared library, SLURM submission files
meta/   input tables, decision tables, and results
```

`code/` is flat because every script sources `file.path(Sys.getenv("HOME"), "ca_config.R")`. `meta/` is the directory `ca_meta()` resolves to when `CA_ARCHIVE` is set to the repository root.

Figure scripts (stage 11) are not included. The estimates shown in the figures are in `meta/did_results_fig.csv.gz`.

## Pairings

| Code | Treatment | Business-as-usual control |
|---|---|---|
| UP_UNP | Undisturbed PAs | Undisturbed non-PAs |
| ONP_UNP | Forest offsets, non-PAs | Undisturbed non-PAs |
| FP_UP | Wildfire, PAs | Undisturbed PAs |
| FNP_UNP | Wildfire, non-PAs | Undisturbed non-PAs |
| TP_UP | Thinning, PAs | Undisturbed PAs |
| TNP_UNP | Thinning, non-PAs | Undisturbed non-PAs |

## Stages

| Stage | Scripts | Output |
|---|---|---|
| Library | `ca_config.R`, `ca_extract.R`, `ca_clean.R` | Analytical rules, paths, and shared functions |
| 0 | `download_almanac.sh` | Wildland Almanac layers |
| 1 | `1_forest_point_grid.R` | Forest point grid on NLCD 2001 |
| 3 | `3a_*`, `3b_*`, `3c_*`, `3d_*` | Carbon, disturbance, treatment, and covariate extraction |
| 4 | `4b_clean_covariates.R` | Unit conversion, cleaning, matching covariates |
| 5 | `5a_recode_treatment.R`, `5b_recode_screen.R` | Event tables, disturbance screen |
| 6 | `6a_group_codes.R` | Treatment and control groups |
| 7 | `7_sample_grid.R` | Grid-based sampling of treated pixels |
| 8 | `8_match_pairings.R`, `8b_aggregate_match.R` | Matched samples, covariate balance, specification selection |
| 9 | `9a_did_units.R`, `9b_did_panel.R` | Difference-in-differences panels |
| 10 | `10a_did_attgt.R`, `10b_did_etwfe.R`, `10c_aggregate_did.R`, `10d_fig_inputs.R` | Pre-treatment trends, ATT estimates, figure windows |

Treatment effects are estimated with extended two-way fixed effects (`etwfe`). Pre-treatment trends are assessed with `att_gt` and `aggte` from the `did` package. Scripts named `*_verify_*` write a verification CSV to `meta/`. Each script opens with a header stating its inputs and outputs.

## Running

R 4.3. Packages: `terra`, `sf`, `arrow`, `dplyr`, `tidyverse`, `MatchIt`, `glmnet`, `dbarts`, `did`, `etwfe` (0.6.2), `marginaleffects` (0.32.0).

| Variable | Set to |
|---|---|
| `HOME` | `code/`, or copy `code/` into `$HOME` |
| `CA_ARCHIVE` | The repository root, so `meta/` resolves |
| `CA_SCRATCH` | A working directory for extracted data and panels |

Stages 1–10c run as SLURM array jobs. Most scripts accept `tasks` (print the array size), `probe` (checks before a full run), or no argument (run one array task). Stage 10d runs interactively. The `.sub` files carry cluster-specific partitions and log paths.

## Input data

Input datasets are publicly available and are not redistributed. Sources, periods, and resolutions are listed in Table S1 of the Supporting Information.

## meta/

| File | Contents |
|---|---|
| `knight_intensity_crosswalk.csv` | Thinning activity to intensity, from Knight et al. (2022) Tables S4–S5 |
| `activity_alias.csv` | Thinning activities absent from Knight et al. (2022), with assigned intensity and rationale |
| `carb_offset_start_years.csv` | CARB offset project start years |
| `screen_threshold.csv`, `nondisturbing_decision.csv`, `did_subset_levels.csv` | Decision tables written at stages 4, 6, and 10 and read by later stages |
| `match_summary_config.csv`, `match_summary_covariate.csv` | Match rates and covariate balance (SMD, VR, eCDF) |
| `did_results.csv.gz` | ATT estimates written at stage 10c |
| `did_results_coverage.csv` | Completeness of stage 10 estimation |
| `event_support.csv`, `fig_caps.csv`, `fig_caps_rev.csv` | Event-time windows for the figures |
| `did_results_fig.csv.gz` | Estimates shown in the figures |

AGB was analyzed with and without log-transformation to estimate proportional and physical changes. GPP, tree NPP, and NEP were analyzed in physical units. AGB is in gC/m² and fluxes in gC/m²/yr.

## Contact

Min Gon Chung, Cooperative Institute for Research in Environmental Sciences, University of Colorado Boulder.
