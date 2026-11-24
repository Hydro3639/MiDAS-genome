#!/usr/bin/env bash
# Bash entry point for the GTDB-domain-guided genome QC workflow.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${QC_CONFIG:-$ROOT/config.sh}"
ACTION=help

usage() {
    cat <<'HELP'
Run a QC stage using config.sh or explicit input arguments.

  bash genome_qc_pipeline.sh full
  bash genome_qc_pipeline.sh full --config config.local.sh
  bash genome_qc_pipeline.sh prepare
  bash genome_qc_pipeline.sh checkm2
  bash genome_qc_pipeline.sh gtdbtk
  bash genome_qc_pipeline.sh domains
  bash genome_qc_pipeline.sh bakta
  bash genome_qc_pipeline.sh barrnap
  bash genome_qc_pipeline.sh trnascan
  bash genome_qc_pipeline.sh summarize
  bash genome_qc_pipeline.sh aggregate
  bash genome_qc_pipeline.sh rna
  bash genome_qc_pipeline.sh status

Options
  --config FILE        Configuration file, default config.sh beside this script
  --genomes DIR        Input genome folder
  --output DIR         Output directory
  --checkm2-db FILE     CheckM2 DIAMOND database
  --gtdbtk-db DIR       GTDB-Tk database
  --bakta-db DIR        Bakta database
  --threads INT        CheckM2 and GTDB-Tk threads
  --jobs INT           Concurrent genome annotation jobs
  -h, --help           Show help

Full stage order
  prepare, checkm2, gtdbtk, domains, bakta, barrnap, trnascan, summarize, aggregate

The rna action runs bakta, barrnap, trnascan and summarize.
GTDB-Tk domain assignments select the Barrnap and tRNAscan-SE models.
Run gtdbtk and domains before running the RNA stages separately.
Additional resources, software environments and thresholds are set in config.sh.
Command-line paths use the launch directory. Configuration paths use the
directory containing the selected configuration file.
HELP
}

cli_path() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$PWD" "$1" ;; esac; }
need_value() { (( $# >= 2 )) && [[ -n "$2" ]] || { printf 'Missing value after %s\n' "$1" >&2; exit 2; }; }
while (( $# > 0 )); do
    case "$1" in
        help|-h|--help) usage; exit 0 ;;
        prepare|checkm2|gtdbtk|domains|domain|bakta|barrnap|trnascan|trnascan-se|summarize|summary|aggregate|finalise|finalize|rna|full|status|check-config)
            ACTION="$1"; shift ;;
        --config) need_value "$@"; CONFIG="$(cli_path "$2")"; shift 2 ;;
        --genomes) need_value "$@"; export QC_OVERRIDE_GENOMES_DIR="$(cli_path "$2")"; shift 2 ;;
        --output) need_value "$@"; export QC_OVERRIDE_OUTDIR="$(cli_path "$2")"; shift 2 ;;
        --checkm2-db) need_value "$@"; export QC_OVERRIDE_CHECKM2_DB="$(cli_path "$2")"; shift 2 ;;
        --gtdbtk-db) need_value "$@"; export QC_OVERRIDE_GTDBTK_DB="$(cli_path "$2")"; shift 2 ;;
        --bakta-db) need_value "$@"; export QC_OVERRIDE_BAKTA_DB="$(cli_path "$2")"; shift 2 ;;
        --threads) need_value "$@"; export QC_OVERRIDE_CHECKM2_THREADS="$2" QC_OVERRIDE_GTDBTK_THREADS="$2"; shift 2 ;;
        --jobs) need_value "$@"; export QC_OVERRIDE_PARALLEL_GENOMES="$2"; shift 2 ;;
        *) printf 'Unknown argument %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done
if [[ "$ACTION" == help ]]; then usage; exit 0; fi
[[ -r "$CONFIG" ]] || { printf 'Cannot read configuration %s\n' "$CONFIG" >&2; exit 1; }
export QC_CONFIG="$(cli_path "$CONFIG")"
source "$ROOT/lib/common.sh"

run_module() { bash "$ROOT/modules/$1.sh"; }

status() {
    printf 'OUTDIR\t%s\n' "$OUTDIR"
    for stage in checkm2 gtdbtk domains bakta barrnap trnascan summarize aggregate; do
        state=PENDING
        [[ ! -f "$MARKERS/$stage.done" ]] || state=DONE
        printf '%s\t%s\n' "$stage" "$state"
    done
    if [[ -s "$WORK/gtdb_domain_counts.tsv" ]]; then cat "$WORK/gtdb_domain_counts.tsv"; fi
    if [[ -s "$RESULTS/final/status_counts.tsv" ]]; then cat "$RESULTS/final/status_counts.tsv"; fi
}

full() {
    run_module prepare
    run_module checkm2
    run_module gtdbtk
    run_module domains
    run_module bakta
    run_module barrnap
    run_module trnascan
    run_module summarize
    run_module aggregate
}

rna() {
    if [[ ! -s "$DOMAIN_MAP" ]]; then
        die 'Run gtdbtk and domains before the RNA stages'
    fi
    run_module bakta
    run_module barrnap
    run_module trnascan
    run_module summarize
}

if [[ "$ACTION" == status ]]; then status; exit 0; fi
case "$ACTION" in full|prepare|check-config) validate_config ;; esac
mkdir -p "$WORK"
exec 9>"$WORK/.run.lock"
flock -n 9 || die "Another QC workflow is running in $OUTDIR"

case "$ACTION" in
    prepare) run_module prepare ;;
    checkm2) run_module checkm2 ;;
    gtdbtk) run_module gtdbtk ;;
    domains|domain) run_module domains ;;
    bakta) run_module bakta ;;
    barrnap) run_module barrnap ;;
    trnascan|trnascan-se) run_module trnascan ;;
    summarize|summary) run_module summarize ;;
    aggregate|finalise|finalize) run_module aggregate ;;
    rna) rna ;;
    full) full ;;
    check-config) run_module prepare; msg 'Configuration OK' ;;
esac
