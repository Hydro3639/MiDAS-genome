#!/usr/bin/env bash
# Compare two genome catalogs using skani and configurable ANI/AF thresholds.
set -Eeuo pipefail
export LC_ALL=C
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${OVERLAP_CONFIG:-$SCRIPT_DIR/config.sh}"
HELPER="$SCRIPT_DIR/scripts/overlap_tables.py"
STAGE=run
declare -A CLI_OVERRIDES=()

usage() {
    cat <<'HELP'
Usage
  bash run_overlap.sh [run|prepare|sketch|search|summarize] [options]

Stages
  run        Validate inputs, sketch catalog B, search catalog A and summarize
  prepare    Validate the two metadata tables and prepare FASTA lists
  sketch     Prepare inputs and build or reuse catalog B sketches
  search     Prepare inputs, sketch catalog B and search catalog A
  summarize  Rebuild reports from compatible completed raw search results

Options
  --config FILE       Configuration file, default config.sh beside this script
  --catalog-a FILE    Catalog A metadata TSV, used as the query catalog
  --catalog-b FILE    Catalog B metadata TSV, used as the reference catalog
  --label-a NAME      Display label for catalog A, default Catalog_A
  --label-b NAME      Display label for catalog B, default Catalog_B
  --root DIR          Base for relative FASTA paths inside metadata tables
  --output DIR        Output directory
  --threads INT       Number of skani threads
  --ani NUMBER        Minimum ANI percentage, default 95
  --af NUMBER         Minimum aligned fraction percentage, default 65
  --af-mode MODE      max or both, default max
  -h, --help          Show help

Edit config.sh or supply input options. Command-line paths use the launch
directory. Configuration paths use PROJECT_ROOT, which defaults to the selected
configuration file's directory. FASTA paths inside metadata use PROJECT_ROOT.
Set SKANI_ENV in config.sh to select a Conda environment for skani, or use PATH
to run the executable available in the launching shell.
HELP
}

cli_path() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$PWD" "$1" ;; esac; }
need_value() { (( $# >= 2 )) && [[ -n "$2" ]] || { printf 'Missing value after %s\n' "$1" >&2; exit 2; }; }
while (( $# > 0 )); do
    case "$1" in
        -h|--help|help) usage; exit 0 ;;
        run|prepare|sketch|search|summarize) STAGE="$1"; shift ;;
        --config) need_value "$@"; CONFIG="$(cli_path "$2")"; shift 2 ;;
        --catalog-a) need_value "$@"; CLI_OVERRIDES[CATALOG_A_METADATA]="$(cli_path "$2")"; shift 2 ;;
        --catalog-b) need_value "$@"; CLI_OVERRIDES[CATALOG_B_METADATA]="$(cli_path "$2")"; shift 2 ;;
        --label-a) need_value "$@"; CLI_OVERRIDES[LABEL_A]="$2"; shift 2 ;;
        --label-b) need_value "$@"; CLI_OVERRIDES[LABEL_B]="$2"; shift 2 ;;
        --root) need_value "$@"; CLI_OVERRIDES[PROJECT_ROOT]="$(cli_path "$2")"; shift 2 ;;
        --output) need_value "$@"; CLI_OVERRIDES[OUT_DIR]="$(cli_path "$2")"; shift 2 ;;
        --threads) need_value "$@"; CLI_OVERRIDES[THREADS]="$2"; shift 2 ;;
        --ani) need_value "$@"; CLI_OVERRIDES[ANI_THRESHOLD]="$2"; shift 2 ;;
        --af) need_value "$@"; CLI_OVERRIDES[AF_THRESHOLD]="$2"; shift 2 ;;
        --af-mode) need_value "$@"; CLI_OVERRIDES[AF_MODE]="$2"; shift 2 ;;
        *) printf 'Unknown argument %s\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
done
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
say() { printf '[%s] %s\n' "$(date -u +%FT%TZ)" "$*"; }
[[ -r "$CONFIG" ]] || die "Cannot read configuration $CONFIG"
CONFIG_DIR="$(cd -- "$(dirname -- "$CONFIG")" && pwd)"
source "$CONFIG"
for variable in "${!CLI_OVERRIDES[@]}"; do
    printf -v "$variable" '%s' "${CLI_OVERRIDES[$variable]}"
done
[[ "$PROJECT_ROOT" == /* ]] || PROJECT_ROOT="$CONFIG_DIR/$PROJECT_ROOT"
[[ "$OUT_DIR" == /* ]] || OUT_DIR="$PROJECT_ROOT/$OUT_DIR"
SKANI_ENV="${SKANI_ENV:-PATH}"
if [[ "$SKANI_ENV" != PATH && "$SKANI_ENV" == */* && "$SKANI_ENV" != /* ]]; then
    SKANI_ENV="$PROJECT_ROOT/$SKANI_ENV"
fi
for variable in CATALOG_A_METADATA CATALOG_B_METADATA; do
    [[ -n "${!variable}" ]] || die "Set $variable in config.sh or supply its input option"
done
for dependency in "$PYTHON" cmp awk flock; do
    command -v "$dependency" >/dev/null || die "Missing dependency: $dependency"
done
for variable in THREADS SKETCH_C MARKER_C; do
    [[ "${!variable}" =~ ^[1-9][0-9]*$ ]] || die "$variable must be a positive integer"
done
[[ "$AF_MODE" == max || "$AF_MODE" == both ]] || die 'AF_MODE must be max or both'
"$PYTHON" "$HELPER" check-parameters --ani "$ANI_THRESHOLD" --af "$AF_THRESHOLD" \
    --screen "$SCREEN_ANI" --search-af "$SEARCH_MIN_AF" --label-a "$LABEL_A" --label-b "$LABEL_B"
mkdir -p "$OUT_DIR"/{01_input,02_sketch,03_search,04_summary,logs,record}
OUT_DIR="$(cd -- "$OUT_DIR" && pwd)"
exec 9>"$OUT_DIR/.run.lock"
flock -n 9 || die "Another workflow is running in $OUT_DIR"
trap 'printf "ERROR at line %s. Check %s/logs/.\n" "$LINENO" "$OUT_DIR" >&2' ERR

prepare() {
    say '[1/4] Validate the two catalog inputs'
    if ! "$PYTHON" "$HELPER" prepare --root "$PROJECT_ROOT" --out "$OUT_DIR/01_input" \
        --catalog-a "$CATALOG_A_METADATA" --catalog-b "$CATALOG_B_METADATA" \
        --label-a "$LABEL_A" --label-b "$LABEL_B" >"$OUT_DIR/logs/01_prepare.log" 2>&1; then
        cat "$OUT_DIR/logs/01_prepare.log" >&2
        return 1
    fi
    cat "$OUT_DIR/01_input/catalog_counts.tsv"
    REFERENCE_N="$(awk -F '\t' '$1=="b" {print $3}' "$OUT_DIR/01_input/catalog_counts.tsv")"
    [[ "$REFERENCE_N" =~ ^[1-9][0-9]*$ ]] || die 'Catalog B has no reference genomes'
}

run_skani() {
    if [[ "$SKANI_ENV" == PATH ]]; then
        "$SKANI" "$@"
    elif [[ "$SKANI_ENV" == */* ]]; then
        conda run --no-capture-output -p "$SKANI_ENV" -- "$SKANI" "$@"
    else
        conda run --no-capture-output -n "$SKANI_ENV" -- "$SKANI" "$@"
    fi
}

setup_skani() {
    if [[ "$SKANI_ENV" == PATH ]]; then
        command -v "$SKANI" >/dev/null || die "Cannot find skani: $SKANI"
    else
        command -v conda >/dev/null || die 'conda must be available to use SKANI_ENV'
        if [[ "$SKANI_ENV" == */* ]]; then
            [[ -d "$SKANI_ENV" ]] || die "Conda environment prefix does not exist: $SKANI_ENV"
        fi
    fi
    VERSION="$(run_skani --version)"
    "$PYTHON" "$HELPER" check-version "$VERSION"
    printf '%s\n' "$VERSION" >"$OUT_DIR/record/skani.version.txt"
    run_skani sketch --help >"$OUT_DIR/record/skani.sketch.help.txt"
    run_skani search --help >"$OUT_DIR/record/skani.search.help.txt"
    for flag in --ql --min-af; do
        "$PYTHON" -c 'import sys; sys.exit(sys.argv[2] not in open(sys.argv[1]).read())' \
            "$OUT_DIR/record/skani.search.help.txt" "$flag" || die "Installed skani lacks $flag"
    done
    DB_DIR="$OUT_DIR/02_sketch/catalog_b"
    DB_INPUTS="$OUT_DIR/record/sketch_inputs.tsv"
    {
        cat "$OUT_DIR/01_input/catalog_b.files.tsv"
        printf 'skani_version\t%s\nsketch_c\t%s\nmarker_c\t%s\n' "$VERSION" "$SKETCH_C" "$MARKER_C"
    } >"$DB_INPUTS"
    {
        printf 'parameter\tvalue\n'
        printf 'catalog_a_label\t%s\ncatalog_b_label\t%s\n' "$LABEL_A" "$LABEL_B"
        printf 'catalog_a_metadata\t%s\ncatalog_b_metadata\t%s\n' "$CATALOG_A_METADATA" "$CATALOG_B_METADATA"
        printf 'project_root\t%s\ntaxonomy_release\t%s\n' "$PROJECT_ROOT" "$TAXONOMY_RELEASE"
        printf 'skani_version\t%s\nthreads\t%s\n' "$VERSION" "$THREADS"
        printf 'skani_executable\t%s\nskani_environment\t%s\n' "$SKANI" "$SKANI_ENV"
        printf 'ANI_threshold_pct\t%s\nAF_threshold_pct\t%s\nAF_mode\t%s\n' "$ANI_THRESHOLD" "$AF_THRESHOLD" "$AF_MODE"
        printf 'screen_ANI_pct\t%s\nsearch_min_AF_pct\t%s\n' "$SCREEN_ANI" "$SEARCH_MIN_AF"
        printf 'sketch_c\t%s\nmarker_c\t%s\n' "$SKETCH_C" "$MARKER_C"
        printf 'active_reference_database\t%s\nrun_time_UTC\t%s\n' "$DB_DIR" "$(date -u +%FT%TZ)"
    } >"$OUT_DIR/record/parameters.tsv"
}

sketch() {
    say '[2/4] Build or reuse catalog B sketches'
    if [[ -d "$DB_DIR" && -s "$DB_DIR.complete" ]] && cmp -s "$DB_INPUTS" "$DB_DIR.complete"; then
        say "Reuse reference sketches: $DB_DIR"
        return
    fi
    [[ ! -e "$DB_DIR" ]] || die "Reference inputs/sketch settings changed, or database is incomplete; use a new OUT_DIR"
    local temp_root temp_db
    temp_root="$(mktemp -d "$OUT_DIR/02_sketch/.building.XXXXXX")"
    temp_db="$temp_root/database"
    # skani creates its own output directory; only reserve the private parent.
    run_skani sketch -l "$OUT_DIR/01_input/catalog_b.paths.txt" -o "$temp_db" \
        -t "$THREADS" -c "$SKETCH_C" -m "$MARKER_C" \
        >"$OUT_DIR/logs/02_sketch.log" 2>&1 || die 'skani sketch failed; see 02_sketch.log'
    "$PYTHON" "$HELPER" audit-log "$OUT_DIR/logs/02_sketch.log"
    [[ -d "$temp_db" && -n "$(ls -A -- "$temp_db")" ]] || die 'skani did not create a nonempty sketch database'
    mv -- "$temp_db" "$DB_DIR"
    rmdir -- "$temp_root"
    cp -- "$DB_INPUTS" "$DB_DIR.complete"
}

search_inputs() {
    cat "$DB_INPUTS" "$OUT_DIR/01_input/catalog_a.files.tsv"
    printf 'screen_ANI\t%s\nsearch_min_AF\t%s\nreference_count\t%s\n' "$SCREEN_ANI" "$SEARCH_MIN_AF" "$REFERENCE_N"
}

search_catalogs() {
    say '[3/4] Search catalog A against catalog B'
    RAW="$OUT_DIR/03_search/catalog_a_vs_catalog_b.skani.tsv"
    SEARCH_INPUTS="$OUT_DIR/record/search_inputs.tsv"
    search_inputs >"$SEARCH_INPUTS"
    if [[ -s "$RAW" && -s "$RAW.inputs.tsv" ]] && cmp -s "$SEARCH_INPUTS" "$RAW.inputs.tsv"; then
        say "Reuse raw search: $RAW"
        return
    fi
    local temp_raw
    temp_raw="$(mktemp "$OUT_DIR/03_search/.search.XXXXXX.tsv")"
    # Retain all reported reference candidates instead of truncating to a top hit.
    run_skani search -d "$DB_DIR" --ql "$OUT_DIR/01_input/catalog_a.paths.txt" \
        -o "$temp_raw" -t "$THREADS" -n "$REFERENCE_N" -s "$SCREEN_ANI" --min-af "$SEARCH_MIN_AF" \
        >"$OUT_DIR/logs/03_search.log" 2>&1 || die 'skani search failed; see 03_search.log'
    "$PYTHON" "$HELPER" audit-log "$OUT_DIR/logs/03_search.log"
    "$PYTHON" "$HELPER" check-raw "$temp_raw"
    mv -- "$temp_raw" "$RAW"
    cp -- "$SEARCH_INPUTS" "$RAW.inputs.tsv"
}

summarize() {
    say '[4/4] Write overlap reports for both catalogs'
    local raw="$OUT_DIR/03_search/catalog_a_vs_catalog_b.skani.tsv"
    local current_inputs="$OUT_DIR/record/search_inputs.tsv"
    search_inputs >"$current_inputs"
    [[ -s "$raw" && -s "$raw.inputs.tsv" ]] || die 'No completed raw output; run the search stage first'
    cmp -s "$current_inputs" "$raw.inputs.tsv" || die 'Raw search inputs/settings changed; rerun the search stage'
    "$PYTHON" "$HELPER" summarize --catalog-a "$OUT_DIR/01_input/catalog_a.metadata.tsv" \
        --catalog-b "$OUT_DIR/01_input/catalog_b.metadata.tsv" --raw "$raw" --out "$OUT_DIR/04_summary" \
        --label-a "$LABEL_A" --label-b "$LABEL_B" --release "$TAXONOMY_RELEASE" \
        --ani "$ANI_THRESHOLD" --af "$AF_THRESHOLD" --af-mode "$AF_MODE"
    cat "$OUT_DIR/04_summary/summary.tsv"
}

prepare
[[ "$STAGE" == prepare ]] && exit 0
setup_skani
case "$STAGE" in
    sketch) sketch ;;
    search) sketch; search_catalogs ;;
    summarize) summarize ;;
    run) sketch; search_catalogs; summarize ;;
esac
say "Finished: $OUT_DIR"
