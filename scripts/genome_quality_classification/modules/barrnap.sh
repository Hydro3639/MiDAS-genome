#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
load_envs

[[ -s "$DOMAIN_MAP" ]] || die "Run: genome_qc_pipeline.sh gtdbtk && genome_qc_pipeline.sh domains first"

outroot="$RESULTS/rna_raw/barrnap"
logroot="$LOGS/barrnap"
doneroot="$outroot/.done"
mkdir -p "$outroot" "$logroot" "$doneroot" "$WORK/tmp"
rm -f "$MARKERS/barrnap.done"

tasks="$WORK/barrnap_tasks.tsv"
build_domain_tasks "$tasks"

barrnap_one() {
    local gid="$1" fasta="$2" domain="${3:-Unknown}"
    local td seq model outfile donef logdir
    td="$(mktemp -d "$WORK/tmp/${gid}.barrnap.XXXXXX")"
    trap 'rm -rf "$td"' RETURN
    seq="$(materialize_fasta "$fasta" "$td")"
    logdir="$logroot/$gid"; mkdir -p "$logdir"

    models=()
    case "$domain" in
        Bacteria)
            models=(bac)
            # Remove stale archaeal fallback files from an earlier Unknown run.
            rm -f "$outroot/$gid.arc.gff3" "$doneroot/$gid.arc.done"
            ;;
        Archaea)
            models=(arc)
            # Remove stale bacterial fallback files from an earlier Unknown run.
            rm -f "$outroot/$gid.bac.gff3" "$doneroot/$gid.bac.done"
            ;;
        *) models=(bac arc) ;;
    esac

    for model in "${models[@]}"; do
        outfile="$outroot/$gid.$model.gff3"
        donef="$doneroot/$gid.$model.done"
        [[ -f "$donef" ]] && continue
        rm -f "$outfile"
        with_env "$BARRNAP_PREFIX" barrnap \
            --threads "$BARRNAP_THREADS_PER_GENOME" \
            --kingdom "$model" --quiet "$seq" \
            >"$outfile" 2>"$logdir/$model.stderr.log" || return 1
        touch "$donef"
    done
}

msg "barrnap: models selected by GTDB-Tk domain"
if parallel_tsv "$tasks" barrnap_one "$PARALLEL_GENOMES" "$WORK/barrnap_failed.txt"; then
    touch "$MARKERS/barrnap.done"
    msg "barrnap: complete"
else
    die "barrnap: one or more genomes failed"
fi
