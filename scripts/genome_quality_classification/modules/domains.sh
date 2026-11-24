#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/lib/common.sh"
ensure_prepared
[[ -f "$MARKERS/gtdbtk.done" ]] || die "Run: genome_qc_pipeline.sh gtdbtk first"
python3 "$ROOT/scripts/parse_gtdb_domains.py" \
    --manifest "$MANIFEST" \
    --gtdb-root "$RESULTS/gtdbtk" \
    --output "$DOMAIN_MAP"

# Print domain counts. The final quality table retains this domain assignment.
printf 'domain\tn_genomes\n' > "$WORK/gtdb_domain_counts.tsv"
tail -n +2 "$DOMAIN_MAP" | cut -f2 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2,$1}' >> "$WORK/gtdb_domain_counts.tsv"

touch "$MARKERS/domains.done"
msg "domains: internal Bacteria/Archaea map -> $DOMAIN_MAP"
cat "$WORK/gtdb_domain_counts.tsv" >&2
