#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
load_envs

marker="$MARKERS/gtdbtk.done"
out="$RESULTS/gtdbtk"
batch="$WORK/gtdbtk_batch.tsv"
if [[ -f "$marker" ]]; then
    msg "gtdbtk: already complete; remove $marker to rerun"
    exit 0
fi

# GTDB-Tk batchfile format: FASTA_PATH<TAB>GENOME_ID.
awk -F '\t' 'BEGIN{OFS="\t"} FNR>1 {print $2,$1}' "$MANIFEST" > "$batch"
rm -rf "$out"
mkdir -p "$out"

extra=()
if [[ -n "$GTDBTK_SCRATCH_DIR" ]]; then
    mkdir -p "$GTDBTK_SCRATCH_DIR"
    extra+=(--scratch_dir "$GTDBTK_SCRATCH_DIR")
fi
if [[ -n "$GTDBTK_EXTRA_ARGS" ]]; then
    # Intended for simple whitespace-separated flags/options configured by the user.
    read -r -a user_extra <<< "$GTDBTK_EXTRA_ARGS"
    extra+=("${user_extra[@]}")
fi

msg "gtdbtk: classify_wf for $(($(wc -l < "$batch"))) genomes"
GTDBTK_DATA_PATH="$GTDBTK_DB" with_env "$GTDBTK_PREFIX" gtdbtk classify_wf \
    --batchfile "$batch" \
    --out_dir "$out" \
    --prefix gtdbtk \
    --cpus "$GTDBTK_THREADS" \
    --pplacer_cpus "$GTDBTK_PPLACER_THREADS" \
    --force \
    "${extra[@]}" \
    >"$out/stdout.log" 2>"$out/stderr.log"

touch "$marker"
msg "gtdbtk: complete -> $out"
