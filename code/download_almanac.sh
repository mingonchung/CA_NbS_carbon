#!/bin/bash
# CA carbon revision pipeline
# Download the Wildland Almanac California v2026.1 layers this pipeline needs.
#
# The stage 3 scripts do NOT download anything. They read from disk only.
# Run this once, before the first probe.
#
# Downloads 245 files into the folder layout the registries expect. Only the
# six layers used by 3a and 3b are pulled, not all 18 properties. From Fire_LCP
# only the CH band and the three static topographic bands are taken, not all
# eight bands.
#
# Filenames follow WildlandAlmanac_CA_{Property}_{WaterYear}.tif, so the file
# list is deterministic and no listing call is needed.
#
# Where to run. This is network and disk bound, not CPU bound. Use a compile
# node (`acompile`) or a short job, not a login node, and not an array task.
# Existing complete files are skipped, so the script is safe to rerun after an
# interruption.
#
# Usage
#   bash download_almanac.sh

set -u

ROOT="/scratch/alpine/mich9173/CA_carbon/input/almanac"
BASE="https://data.source.coop/wildland-almanac/california/v2026.1"

mkdir -p "$ROOT"/{Carbon_AGB,Carbon_GPP,Fire_LCP,Disturbance_TreeFrac,Disturbance_AGB,Veg_TreeFrac}

# Scaffold folders for a future release that carries these. Empty for now.
mkdir -p "$ROOT"/{Carbon_NPP,Carbon_NEP,Carbon_NBP}

fetch () {
  local url="$1" out="$2"
  if [ -s "$out" ]; then
    echo "skip  $(basename "$out")"
    return 0
  fi
  # -C - resumes a partial file, --fail turns an HTTP error into a nonzero exit
  curl -sS --fail --retry 5 --retry-delay 5 -C - -o "$out" "$url" \
    && echo "ok    $(basename "$out")" \
    || echo "FAIL  $(basename "$out")"
}

echo "=== Carbon_AGB and Carbon_GPP and Veg_TreeFrac, 1985-2025 ==="
for YEAR in $(seq 1985 2025); do
  for LAYER in Carbon_AGB Carbon_GPP Veg_TreeFrac; do
    FILE="WildlandAlmanac_CA_${LAYER}_${YEAR}.tif"
    fetch "${BASE}/${LAYER}/${FILE}" "${ROOT}/${LAYER}/${FILE}"
  done
done

echo "=== Disturbance layers, 1986-2024 ==="
# Boundary years 1985 and 2025 are not produced. Each annual estimate needs a
# bracketing pre and post year, which is unavailable at the ends of the series.
for YEAR in $(seq 1986 2024); do
  for LAYER in Disturbance_TreeFrac Disturbance_AGB; do
    FILE="WildlandAlmanac_CA_${LAYER}_${YEAR}.tif"
    fetch "${BASE}/${LAYER}/${FILE}" "${ROOT}/${LAYER}/${FILE}"
  done
done

echo "=== Fire_LCP canopy height, 1985-2025 ==="
# Canopy height only. Stage 3a reads this band alone, matched by the file
# pattern in ca_almanac_layers. Units are decimetres, scaled to metres at
# extraction. Nodata is 0 here, not -9999.
for YEAR in $(seq 1985 2025); do
  FILE="WildlandAlmanac_CA_Fire_LCP_CH_${YEAR}.tif"
  fetch "${BASE}/Fire_LCP/${FILE}" "${ROOT}/Fire_LCP/${FILE}"
done

echo "=== Fire_LCP static topography ==="
# No year in the filename. Candidates to replace the separately derived
# elevation, slope, and aspect covariates at stage 3d, which would put every
# matching covariate on one co-registered stack.
for BAND in Elevation Slope Aspect; do
  FILE="WildlandAlmanac_CA_Fire_LCP_${BAND}.tif"
  fetch "${BASE}/Fire_LCP/${FILE}" "${ROOT}/Fire_LCP/${FILE}"
done

echo
echo "=== Summary ==="
for LAYER in Carbon_AGB Carbon_GPP Veg_TreeFrac Disturbance_TreeFrac Disturbance_AGB Fire_LCP; do
  N=$(ls -1 "${ROOT}/${LAYER}" 2>/dev/null | wc -l)
  SZ=$(du -sh "${ROOT}/${LAYER}" 2>/dev/null | cut -f1)
  printf "%-22s %4s files  %8s\n" "$LAYER" "$N" "$SZ"
done
echo
echo "Expected: Carbon_AGB 41, Carbon_GPP 41, Veg_TreeFrac 41,"
echo "          Disturbance_TreeFrac 39, Disturbance_AGB 39, Fire_LCP 44"
echo
echo "Any FAIL lines above mean a retry is needed. Rerun this script."
