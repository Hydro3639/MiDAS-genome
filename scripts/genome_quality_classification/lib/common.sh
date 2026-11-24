#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${QC_CONFIG:-$PIPELINE_ROOT/config.sh}"
for option_variable in GENOMES_DIR OUTDIR CHECKM2_DB GTDBTK_DB BAKTA_DB CHECKM2_THREADS GTDBTK_THREADS PARALLEL_GENOMES; do
    override_variable="QC_OVERRIDE_$option_variable"
    if [[ -n "${!override_variable-}" ]]; then
        printf -v "$option_variable" '%s' "${!override_variable}"
    fi
done
export LC_ALL=C

WORK="$OUTDIR/work"
RESULTS="$OUTDIR/results"
LOGS="$OUTDIR/logs"
MANIFEST="$WORK/genomes.tsv"
ENV_FILE="$WORK/env_prefixes.sh"
DOMAIN_MAP="$WORK/gtdb_domain.tsv"
MARKERS="$WORK/markers"

msg() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { msg "ERROR: $*"; exit 1; }
require_file() { [[ -f "$1" ]] || die "File not found: $1"; }
require_dir() { [[ -d "$1" ]] || die "Directory not found: $1"; }

resolve_env_prefix() {
    local spec="$1" prefix=""
    if [[ "$spec" == PATH ]]; then printf '%s\n' PATH; return; fi
    [[ -n "$spec" ]] || die "Empty conda environment specification"
    if [[ -d "$spec/bin" ]]; then
        (cd "$spec" && pwd)
        return
    fi
    command -v conda >/dev/null 2>&1 || die "conda is not on PATH; use a full env prefix for '$spec'"
    prefix="$(conda env list | awk -v n="$spec" '$1==n {print $NF; exit}')"
    [[ -n "$prefix" && -d "$prefix/bin" ]] || die "Cannot resolve conda env '$spec'"
    printf '%s\n' "$prefix"
}

with_env() {
    local prefix="$1"; shift
    if [[ "$prefix" == PATH ]]; then "$@"; else env PATH="$prefix/bin:$PATH" "$@"; fi
}

load_envs() {
    require_file "$ENV_FILE"
    # shellcheck source=/dev/null
    source "$ENV_FILE"
}

check_binary() {
    local prefix="$1" exe="$2"
    if [[ "$prefix" == PATH ]]; then
        command -v "$exe" >/dev/null 2>&1 || die "$exe is unavailable on PATH"
    else
        [[ -x "$prefix/bin/$exe" ]] || die "$exe not found in $prefix/bin"
    fi
}

make_manifest() {
    require_dir "$GENOMES_DIR"
    mkdir -p "$WORK" "$RESULTS" "$LOGS" "$MARKERS"
    local tmp="$MANIFEST.tmp"
    printf 'genome\tpath\n' > "$tmp"

    declare -A seen=()
    local n=0 path base gid
    while IFS= read -r -d '' path; do
        base="$(basename "$path")"
        case "$base" in
            *.fna.gz) gid="${base%.fna.gz}" ;;
            *.fa.gz) gid="${base%.fa.gz}" ;;
            *.fasta.gz) gid="${base%.fasta.gz}" ;;
            *.fna) gid="${base%.fna}" ;;
            *.fa) gid="${base%.fa}" ;;
            *.fasta) gid="${base%.fasta}" ;;
            *) continue ;;
        esac
        gid="$(printf '%s' "$gid" | sed -E 's/[^A-Za-z0-9_.-]+/_/g; s/^_+//; s/_+$//')"
        [[ -n "$gid" ]] || die "Genome ID became empty: $path"
        [[ -z "${seen[$gid]+x}" ]] || die "Duplicate genome ID '$gid': ${seen[$gid]} and $path"
        seen[$gid]="$path"
        printf '%s\t%s\n' "$gid" "$path" >> "$tmp"
        ((n+=1))
    done < <(find -L "$GENOMES_DIR" -maxdepth 1 -type f \
        \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' -o \
           -name '*.fna.gz' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) \
        -print0 | sort -z)

    (( n > 0 )) || die "No genomes found in $GENOMES_DIR"
    # Associate reusable annotation outputs with their input files and settings.
    local snapshot="$WORK/input_snapshot.tmp"
    cat "$tmp" > "$snapshot"
    local setting
    for setting in CHECKM2_DB GTDBTK_DB BAKTA_DB CHECKM2_ENV GTDBTK_ENV BAKTA_ENV BARRNAP_ENV TRNASCAN_ENV BAKTA_RNA_ONLY GTDBTK_EXTRA_ARGS CHECKM2_CHUNK_SIZE; do
        printf '%s=%q\n' "$setting" "${!setting}" >> "$snapshot"
    done
    printf 'RNA_MODEL_SELECTION=GTDB_DOMAIN\n' >> "$snapshot"
    while IFS=$'\t' read -r gid path; do
        printf '%s\t%s\n' "$path" "$(stat -Lc '%s:%y' -- "$path")" >> "$snapshot"
    done < <(tail -n +2 "$tmp")
    if [[ -s "$WORK/input_snapshot.tsv" ]] && ! cmp -s "$snapshot" "$WORK/input_snapshot.tsv" && \
        [[ -n "$(find "$RESULTS" -type f -print -quit)" ]]; then
        rm -f "$tmp" "$snapshot"
        die 'Genome inputs or annotation settings changed. Use a separate OUTDIR for those annotations'
    fi
    mv "$tmp" "$MANIFEST"
    mv "$snapshot" "$WORK/input_snapshot.tsv"
    msg "Manifest: $n genomes -> $MANIFEST"
}

write_env_prefixes() {
    mkdir -p "$WORK"
    local c g b r t
    c="$(resolve_env_prefix "$CHECKM2_ENV")"
    g="$(resolve_env_prefix "$GTDBTK_ENV")"
    b="$(resolve_env_prefix "$BAKTA_ENV")"
    r="$(resolve_env_prefix "$BARRNAP_ENV")"
    t="$(resolve_env_prefix "$TRNASCAN_ENV")"
    {
        printf 'CHECKM2_PREFIX=%q\n' "$c"
        printf 'GTDBTK_PREFIX=%q\n' "$g"
        printf 'BAKTA_PREFIX=%q\n' "$b"
        printf 'BARRNAP_PREFIX=%q\n' "$r"
        printf 'TRNASCAN_PREFIX=%q\n' "$t"
    } > "$ENV_FILE"
    # shellcheck source=/dev/null
    source "$ENV_FILE"
    check_binary "$CHECKM2_PREFIX" checkm2
    check_binary "$GTDBTK_PREFIX" gtdbtk
    check_binary "$BAKTA_PREFIX" bakta
    check_binary "$BARRNAP_PREFIX" barrnap
    check_binary "$TRNASCAN_PREFIX" tRNAscan-SE
    msg "Resolved software environments -> $ENV_FILE"
}

validate_config() {
    for variable in GENOMES_DIR OUTDIR CHECKM2_DB GTDBTK_DB BAKTA_DB; do
        [[ -n "${!variable}" ]] || die "Set $variable in config.sh or supply its command-line option"
    done
    require_dir "$GENOMES_DIR"
    require_file "$CHECKM2_DB"
    require_dir "$GTDBTK_DB"
    require_dir "$BAKTA_DB"
    command -v python3 >/dev/null 2>&1 || die "python3 not found"
    command -v flock >/dev/null 2>&1 || die 'flock not found'
    command -v cmp >/dev/null 2>&1 || die 'cmp not found'
    for variable in CHECKM2_THREADS GTDBTK_THREADS GTDBTK_PPLACER_THREADS CHECKM2_CHUNK_SIZE PARALLEL_GENOMES BAKTA_THREADS_PER_GENOME BARRNAP_THREADS_PER_GENOME TRNASCAN_THREADS_PER_GENOME MIN_UNIQUE_TRNA; do
        [[ "${!variable}" =~ ^[1-9][0-9]*$ ]] || die "$variable must be a positive integer"
    done
    (( MIN_UNIQUE_TRNA <= 20 )) || die 'MIN_UNIQUE_TRNA must be at most 20 amino-acid types'
}

ensure_prepared() {
    if [[ ! -s "$MANIFEST" || ! -s "$ENV_FILE" ]]; then
        msg "Preparation files missing; running prepare first."
        bash "$PIPELINE_ROOT/modules/prepare.sh"
    fi
}

materialize_fasta() {
    local fasta="$1" tmpdir="$2"
    if [[ "$fasta" == *.gz ]]; then
        local out="$tmpdir/$(basename "${fasta%.gz}")"
        gzip -cd "$fasta" > "$out"
        printf '%s\n' "$out"
    else
        printf '%s\n' "$fasta"
    fi
}

build_domain_tasks() {
    local outfile="$1"
    require_file "$MANIFEST"
    require_file "$DOMAIN_MAP"
    awk -F '\t' 'BEGIN{OFS="\t"} NR==FNR {if(FNR>1)d[$1]=$2; next} FNR>1 {print $1,$2,(($1 in d)?d[$1]:"Unknown")}' \
        "$DOMAIN_MAP" "$MANIFEST" > "$outfile"
}

# Run worker function for a 2- or 3-column TSV: gid, path [, domain].
# Successful jobs continue; failures are collected and the module returns nonzero.
parallel_tsv() {
    local tsv="$1" worker="$2" max_jobs="$3" failure_file="$4"
    : > "$failure_file"
    local failed=0 gid path third pid
    local -a worker_pids=()
    while IFS=$'\t' read -r gid path third; do
        [[ -n "$gid" ]] || continue
        (
            set -euo pipefail
            "$worker" "$gid" "$path" "${third:-}"
        ) || { printf '%s\n' "$gid" >> "$failure_file"; exit 1; } &
        worker_pids+=("$!")
        if (( ${#worker_pids[@]} >= max_jobs )); then
            for pid in "${worker_pids[@]}"; do wait "$pid" || failed=1; done
            worker_pids=()
        fi
    done < "$tsv"
    for pid in "${worker_pids[@]}"; do wait "$pid" || failed=1; done
    if (( failed != 0 )); then
        msg "Failures recorded in $failure_file"
        return 1
    fi
    rm -f "$failure_file"
}
