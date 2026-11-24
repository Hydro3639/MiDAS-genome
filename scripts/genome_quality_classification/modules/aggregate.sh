#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
rm -f "$MARKERS/aggregate.done"

check="$RESULTS/checkm2/checkm2_summary.tsv"
rna="$RESULTS/rna/rna_summary.tsv"
require_file "$check"
require_file "$rna"
mkdir -p "$RESULTS/final"

python3 "$ROOT/scripts/aggregate.py" \
    --manifest "$MANIFEST" \
    --checkm2 "$check" \
    --rna "$rna" \
    --output "$RESULTS/final/all_genomes_quality.tsv" \
    --counts-output "$RESULTS/final/status_counts.tsv" \
    --hq-comp-ge "$HQ_COMPLETENESS_GE" \
    --hq-cont-le "$HQ_CONTAMINATION_LE" \
    --mq-comp-ge "$MQ_COMPLETENESS_GE" \
    --mq-cont-le "$MQ_CONTAMINATION_LE" \
    --min-trna "$MIN_UNIQUE_TRNA"
touch "$MARKERS/aggregate.done"
msg "aggregate: $RESULTS/final/all_genomes_quality.tsv"
cat "$RESULTS/final/status_counts.tsv"
