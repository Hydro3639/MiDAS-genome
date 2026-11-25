#!/usr/bin/env bash
# Supply metadata for the two catalogs, or use the main script's input options.
CONFIG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$CONFIG_DIR}"
OUT_DIR="${OUT_DIR:-results}"
CATALOG_A_METADATA="${CATALOG_A_METADATA:-}"
CATALOG_B_METADATA="${CATALOG_B_METADATA:-}"
LABEL_A="${LABEL_A:-Catalog_A}"
LABEL_B="${LABEL_B:-Catalog_B}"

SKANI="${SKANI:-skani}"
# PATH uses the active shell. Otherwise supply a Conda environment name or prefix.
SKANI_ENV="${SKANI_ENV:-PATH}"
PYTHON="${PYTHON:-python3}"
THREADS="${THREADS:-8}"
TAXONOMY_RELEASE="${TAXONOMY_RELEASE:-unspecified}"

# ANI and aligned fractions use percentages.
ANI_THRESHOLD="${ANI_THRESHOLD:-95}"
AF_THRESHOLD="${AF_THRESHOLD:-65}"
AF_MODE="${AF_MODE:-max}"  # max: either direction; both: both directions

# Broader search settings retain candidates for the final and sensitivity rules.
SCREEN_ANI="${SCREEN_ANI:-80}"
SEARCH_MIN_AF="${SEARCH_MIN_AF:-15}"
SKETCH_C="${SKETCH_C:-125}"
MARKER_C="${MARKER_C:-1000}"
