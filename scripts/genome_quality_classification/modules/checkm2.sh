#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
load_envs

marker="$MARKERS/checkm2.done"
summary="$RESULTS/checkm2/checkm2_summary.tsv"
if [[ -f "$marker" && -s "$summary" ]]; then
    msg "checkm2: already complete; remove $marker to rerun"
    exit 0
fi

chunkdir="$WORK/checkm2_chunks"
outroot="$RESULTS/checkm2/chunks"
logroot="$LOGS/checkm2"
rm -rf "$chunkdir"
mkdir -p "$chunkdir" "$outroot" "$logroot"
tail -n +2 "$MANIFEST" | split -l "$CHECKM2_CHUNK_SIZE" -d -a 5 --additional-suffix=.tsv - "$chunkdir/chunk_"

for chunk in "$chunkdir"/chunk_*.tsv; do
    idx="$(basename "$chunk" .tsv | sed 's/^chunk_//')"
    out="$outroot/$idx"
    if [[ -s "$out/quality_report.tsv" && -f "$out/.complete" ]]; then
        msg "checkm2: chunk $idx already complete; skip"
        continue
    fi
    # CheckM2 requires --output-directory to not exist (or to be empty).
    # Do not pre-create $out; keep logs in the separate pipeline log directory.
    rm -rf "$out"
    logdir="$logroot/$idx"
    mkdir -p "$logdir"
    mapfile -t files < <(cut -f2 "$chunk")
    opts=()
    [[ "$CHECKM2_LOWMEM" == 1 ]] && opts+=(--lowmem)
    msg "checkm2: chunk $idx (${#files[@]} genomes)"
    with_env "$CHECKM2_PREFIX" checkm2 predict \
        --threads "$CHECKM2_THREADS" \
        --input "${files[@]}" \
        --output-directory "$out" \
        --database_path "$CHECKM2_DB" \
        "${opts[@]}" \
        >"$logdir/stdout.log" 2>"$logdir/stderr.log"
    [[ -s "$out/quality_report.tsv" ]] || die "CheckM2 chunk $idx produced no quality_report.tsv"
    touch "$out/.complete"
done

python3 "$ROOT/scripts/merge_checkm2.py" --root "$outroot" --manifest "$MANIFEST" --output "$summary"
touch "$marker"
msg "checkm2: complete -> $summary"
