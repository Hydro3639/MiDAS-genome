#!/usr/bin/env bash
# Set your genome folder and the three database paths before running the workflow.
# Environment variables and command-line options can override these defaults.

CONFIG_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GENOMES_DIR="${GENOMES_DIR:-}"
OUTDIR="${OUTDIR:-$CONFIG_DIR/results}"

# CheckM2 requires its DIAMOND database file. The other two values are folders.
CHECKM2_DB="${CHECKM2_DB:-}"
GTDBTK_DB="${GTDBTK_DB:-${GTDBTK_DATA_PATH:-}}"
BAKTA_DB="${BAKTA_DB:-}"

# Use PATH for tools available in the active shell, or provide Conda names/prefixes.
CHECKM2_ENV="${CHECKM2_ENV:-PATH}"
GTDBTK_ENV="${GTDBTK_ENV:-PATH}"
BAKTA_ENV="${BAKTA_ENV:-PATH}"
BARRNAP_ENV="${BARRNAP_ENV:-PATH}"
TRNASCAN_ENV="${TRNASCAN_ENV:-PATH}"

CHECKM2_THREADS="${CHECKM2_THREADS:-8}"
GTDBTK_THREADS="${GTDBTK_THREADS:-8}"
GTDBTK_PPLACER_THREADS="${GTDBTK_PPLACER_THREADS:-1}"
CHECKM2_CHUNK_SIZE="${CHECKM2_CHUNK_SIZE:-400}"
CHECKM2_LOWMEM="${CHECKM2_LOWMEM:-0}"
PARALLEL_GENOMES="${PARALLEL_GENOMES:-4}"
BAKTA_THREADS_PER_GENOME="${BAKTA_THREADS_PER_GENOME:-1}"
BARRNAP_THREADS_PER_GENOME="${BARRNAP_THREADS_PER_GENOME:-1}"
TRNASCAN_THREADS_PER_GENOME="${TRNASCAN_THREADS_PER_GENOME:-1}"
GTDBTK_SCRATCH_DIR="${GTDBTK_SCRATCH_DIR:-}"
GTDBTK_EXTRA_ARGS="${GTDBTK_EXTRA_ARGS:-}"

# GTDB-Tk domain assignment selects the bacterial or archaeal RNA models.
# Unassigned genomes are evaluated with both model sets.
BAKTA_RNA_ONLY="${BAKTA_RNA_ONLY:-1}"            # 1 for RNA modules, 0 for full Bakta
RRNA_EDGE_MARGIN_BP="${RRNA_EDGE_MARGIN_BP:-10}"

# Inclusive CheckM2 thresholds and minimum number of amino-acid tRNA types.
HQ_COMPLETENESS_GE="${HQ_COMPLETENESS_GE:-90}"
HQ_CONTAMINATION_LE="${HQ_CONTAMINATION_LE:-5}"
MQ_COMPLETENESS_GE="${MQ_COMPLETENESS_GE:-50}"
MQ_CONTAMINATION_LE="${MQ_CONTAMINATION_LE:-10}"
MIN_UNIQUE_TRNA="${MIN_UNIQUE_TRNA:-18}"

# Resolve user-supplied relative paths against this configuration file.
for path_variable in GENOMES_DIR OUTDIR CHECKM2_DB GTDBTK_DB BAKTA_DB GTDBTK_SCRATCH_DIR; do
    path_value="${!path_variable}"
    if [[ -n "$path_value" && "$path_value" != /* ]]; then
        printf -v "$path_variable" '%s' "$CONFIG_DIR/$path_value"
    fi
done
for environment_variable in CHECKM2_ENV GTDBTK_ENV BAKTA_ENV BARRNAP_ENV TRNASCAN_ENV; do
    environment_value="${!environment_variable}"
    if [[ "$environment_value" == */* && "$environment_value" != /* ]]; then
        printf -v "$environment_variable" '%s' "$CONFIG_DIR/$environment_value"
    fi
done
