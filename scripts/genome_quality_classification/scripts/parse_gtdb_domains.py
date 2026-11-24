#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_manifest(path):
    with open(path, newline='') as fh:
        return [r['genome'] for r in csv.DictReader(fh, delimiter='\t')]


def domain_from_row(row, file_domain):
    """Return (domain, status) without exporting taxonomy to downstream tables."""
    classification = (row.get('classification') or '').strip()
    low = classification.lower()

    if classification.startswith('d__Bacteria'):
        return 'Bacteria', 'classified'
    if classification.startswith('d__Archaea'):
        return 'Archaea', 'classified'

    # GTDB-Tk >=2.1 may retain genomes filtered during alignment as
    # "Unclassified Bacteria" / "Unclassified Archaea".  Domain is still known.
    if 'unclassified bacteria' in low:
        return 'Bacteria', 'unclassified_bacteria'
    if 'unclassified archaea' in low:
        return 'Archaea', 'unclassified_archaea'

    # A plain "Unclassified" can be emitted when no domain can be assigned.
    if low == 'unclassified' or low.startswith('unclassified;'):
        return 'Unknown', 'unclassified_domain'

    # Older GTDB-Tk versions may not provide a useful classification string in
    # all summary rows.  In that case the summary-file marker set is the best
    # domain signal available.
    if not classification:
        return file_domain, 'summary_file_domain'

    # Taxonomy strings from older versions may still make the domain explicit.
    if 'd__bacteria' in low:
        return 'Bacteria', 'classified'
    if 'd__archaea' in low:
        return 'Archaea', 'classified'

    return file_domain, 'summary_file_domain'


def ingest(path, file_domain, mapping):
    if not path.exists() or path.stat().st_size == 0:
        return
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames:
            return
        col = 'user_genome' if 'user_genome' in reader.fieldnames else reader.fieldnames[0]
        for row in reader:
            gid = (row.get(col) or '').strip()
            if not gid:
                continue
            domain, status = domain_from_row(row, file_domain)
            prev = mapping.get(gid)
            if prev is None:
                mapping[gid] = (domain, status, path.name)
            elif prev[1] == 'conflicting_domain_summaries':
                # Later summary files must not erase a previously detected conflict.
                mapping[gid] = ('Unknown', prev[1], f'{prev[2]};{path.name}')
            elif prev[0] == 'Unknown' and domain != 'Unknown':
                mapping[gid] = (domain, status, path.name)
            elif domain != 'Unknown' and prev[0] != domain:
                # Conflicting non-Unknown domain evidence should never silently
                # choose one side.
                mapping[gid] = ('Unknown', 'conflicting_domain_summaries', f'{prev[2]};{path.name}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--gtdb-root', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()

    ids = read_manifest(a.manifest)
    root = Path(a.gtdb_root)
    mapping = {}

    # Prefix may be changed by the user; do not assume literal "gtdbtk".
    for p in root.rglob('*.bac120.summary.tsv'):
        ingest(p, 'Bacteria', mapping)
    for p in root.rglob('*.ar53.summary.tsv'):
        ingest(p, 'Archaea', mapping)

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['genome', 'domain', 'gtdb_domain_status', 'gtdb_summary_file'])
        for gid in ids:
            domain, status, source = mapping.get(gid, ('Unknown', 'missing_from_summary', ''))
            w.writerow([gid, domain, status, source])


if __name__ == '__main__':
    main()
