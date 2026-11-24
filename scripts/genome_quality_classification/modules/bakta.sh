#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
load_envs

outroot="$RESULTS/rna_raw/bakta"
logroot="$LOGS/bakta"
mkdir -p "$outroot" "$logroot" "$WORK/tmp"
rm -f "$MARKERS/bakta.done"

tasks="$WORK/bakta_tasks.tsv"
tail -n +2 "$MANIFEST" > "$tasks"

bakta_one() {
    local gid="$1" fasta="$2" _unused="${3:-}"
    local out="$outroot/$gid" log="$logroot/$gid"
    mkdir -p "$log"
    if [[ -s "$out/$gid.tsv" && -s "$out/$gid.json" ]]; then
        return 0
    fi
    rm -rf "$out"; mkdir -p "$out"
    opts=(
        --db "$BAKTA_DB" --prefix "$gid" --output "$out"
        --keep-contig-headers --threads "$BAKTA_THREADS_PER_GENOME"
        --meta --force
    )
    if [[ "$BAKTA_RNA_ONLY" == 1 ]]; then
        opts+=(--skip-tmrna --skip-ncrna --skip-ncrna-region --skip-crispr
              --skip-cds --skip-pseudo --skip-sorf --skip-gap --skip-ori --skip-plot)
    fi
    if ! with_env "$BAKTA_PREFIX" bakta "${opts[@]}" "$fasta" \
        >"$log/stdout.log" 2>"$log/stderr.log"; then
        rm -rf -- "$out"
        return 1
    fi
    [[ -s "$out/$gid.tsv" && -s "$out/$gid.json" ]]
}

msg "bakta: all genomes; PARALLEL_GENOMES=$PARALLEL_GENOMES, threads/genome=$BAKTA_THREADS_PER_GENOME"
if parallel_tsv "$tasks" bakta_one "$PARALLEL_GENOMES" "$WORK/bakta_failed.txt"; then
    touch "$MARKERS/bakta.done"
    msg "bakta: complete"
else
    die "bakta: one or more genomes failed"
fi
