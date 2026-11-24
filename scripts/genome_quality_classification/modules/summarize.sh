#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
rm -f "$MARKERS/summarize.done"

[[ -f "$MARKERS/bakta.done" ]] || die "Bakta step is not complete"
[[ -f "$MARKERS/barrnap.done" ]] || die "Barrnap step is not complete"
[[ -f "$MARKERS/trnascan.done" ]] || die "tRNAscan-SE step is not complete"

[[ -s "$DOMAIN_MAP" ]] || die "Run: genome_qc_pipeline.sh gtdbtk && genome_qc_pipeline.sh domains first"

mkdir -p "$RESULTS/rna"
python3 "$ROOT/scripts/summarize_rna.py" \
    --manifest "$MANIFEST" \
    --domain-map "$DOMAIN_MAP" \
    --outroot "$RESULTS/rna_raw" \
    --output "$RESULTS/rna/rna_summary.tsv" \
    --hits-output "$RESULTS/rna/rrna_hits.tsv" \
    --min-trna "$MIN_UNIQUE_TRNA" \
    --edge-margin "$RRNA_EDGE_MARGIN_BP"
touch "$MARKERS/summarize.done"
msg "summarize: $RESULTS/rna/rna_summary.tsv"
msg "summarize: hit-level rRNA QC -> $RESULTS/rna/rrna_hits.tsv"
