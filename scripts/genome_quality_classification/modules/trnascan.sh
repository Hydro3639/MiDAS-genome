#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
load_envs

[[ -s "$DOMAIN_MAP" ]] || die "Run: genome_qc_pipeline.sh gtdbtk && genome_qc_pipeline.sh domains first"

outroot="$RESULTS/rna_raw/trnascan"
logroot="$LOGS/trnascan"
doneroot="$outroot/.done"
mkdir -p "$outroot" "$logroot" "$doneroot" "$WORK/tmp"
rm -f "$MARKERS/trnascan.done"

tasks="$WORK/trnascan_tasks.tsv"
build_domain_tasks "$tasks"

trnascan_one() {
    local gid="$1" fasta="$2" domain="${3:-Unknown}"
    local td seq mode label outfile stats donef logdir
    td="$(mktemp -d "$WORK/tmp/${gid}.trnascan.XXXXXX")"
    trap 'rm -rf "$td"' RETURN
    seq="$(materialize_fasta "$fasta" "$td")"
    logdir="$logroot/$gid"; mkdir -p "$logdir"

    modes=()
    case "$domain" in
        Bacteria)
            modes=(B)
            # Remove stale archaeal fallback files from an earlier Unknown run.
            rm -f "$outroot/$gid.arc.tsv" "$outroot/$gid.arc.stats.txt" "$doneroot/$gid.arc.done"
            ;;
        Archaea)
            modes=(A)
            # Remove stale bacterial fallback files from an earlier Unknown run.
            rm -f "$outroot/$gid.bac.tsv" "$outroot/$gid.bac.stats.txt" "$doneroot/$gid.bac.done"
            ;;
        *) modes=(B A) ;;
    esac

    for mode in "${modes[@]}"; do
        case "$mode" in B) label=bac ;; A) label=arc ;; esac
        outfile="$outroot/$gid.$label.tsv"
        stats="$outroot/$gid.$label.stats.txt"
        donef="$doneroot/$gid.$label.done"
        [[ -f "$donef" ]] && continue
        rm -f "$outfile" "$stats"
        with_env "$TRNASCAN_PREFIX" tRNAscan-SE "-$mode" \
            -o "$outfile" -m "$stats" -d \
            --thread "$TRNASCAN_THREADS_PER_GENOME" "$seq" \
            >"$logdir/$label.stdout.log" 2>"$logdir/$label.stderr.log" || return 1
        touch "$donef"
    done
}

msg "trnascan: models selected by GTDB-Tk domain"
if parallel_tsv "$tasks" trnascan_one "$PARALLEL_GENOMES" "$WORK/trnascan_failed.txt"; then
    touch "$MARKERS/trnascan.done"
    msg "trnascan: complete"
else
    die "trnascan: one or more genomes failed"
fi
