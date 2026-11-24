#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


def read_table(path, required):
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fields = reader.fieldnames or []
        missing = set(required) - set(fields)
        if missing:
            raise ValueError(f'{path}: missing columns {", ".join(sorted(missing))}; regenerate the input table with the current workflow')
        rows = {}
        for row in reader:
            genome = row['genome']
            if not genome or genome in rows:
                raise ValueError(f'{path}: empty or duplicate genome ID {genome!r}')
            rows[genome] = row
        return fields, rows


def as_float(x):
    try:
        value = float(x)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) and 0 <= value <= 100 else None


def as_int(x):
    try: return int(float(x))
    except (TypeError, ValueError, OverflowError): return 0


def b01(x): return str(x).strip().lower() in {'1','true','yes','y'}


def classify(comp, cont, rrna_ok, trna_ok, args):
    if comp is None or cont is None:
        return 'UNASSESSED'
    if comp >= args.hq_comp_ge and cont <= args.hq_cont_le and rrna_ok and trna_ok:
        return 'HQ'
    if comp >= args.mq_comp_ge and cont <= args.mq_cont_le:
        return 'MQ'
    if cont <= args.mq_cont_le:
        return 'LQ'
    return 'Contaminated'


def hq_reasons(comp, cont, r, trna_ok, args):
    if comp is None or cont is None:
        return 'missing_or_invalid_CheckM2'
    reasons = []
    if comp < args.hq_comp_ge: reasons.append(f'completeness_lt_{args.hq_comp_ge:g}')
    if cont > args.hq_cont_le: reasons.append(f'contamination_gt_{args.hq_cont_le:g}')
    if not b01(r.get('merged_5S')): reasons.append('missing_5S')
    if not b01(r.get('merged_16S')): reasons.append('missing_16S')
    if not b01(r.get('merged_23S')): reasons.append('missing_23S')
    if not trna_ok: reasons.append('tRNA_requirement_failed')
    return ';'.join(reasons)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--checkm2', required=True)
    ap.add_argument('--rna', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--counts-output', required=True)
    ap.add_argument('--hq-comp-ge', type=float, default=90.0)
    ap.add_argument('--hq-cont-le', type=float, default=5.0)
    ap.add_argument('--mq-comp-ge', type=float, default=50.0)
    ap.add_argument('--mq-cont-le', type=float, default=10.0)
    ap.add_argument('--min-trna', type=int, default=18)
    args = ap.parse_args()

    thresholds = (args.hq_comp_ge, args.hq_cont_le, args.mq_comp_ge, args.mq_cont_le)
    if not all(math.isfinite(value) and 0 <= value <= 100 for value in thresholds):
        ap.error('Completeness and contamination thresholds must be finite percentages from 0 to 100')
    if args.hq_comp_ge < args.mq_comp_ge or args.hq_cont_le > args.mq_cont_le:
        ap.error('HQ completeness must be at least the MQ threshold and HQ contamination at most the MQ threshold')
    if not 1 <= args.min_trna <= 20:
        ap.error('--min-trna must be between 1 and 20 amino-acid types')

    try:
        _, check = read_table(args.checkm2, ('genome', 'completeness', 'contamination'))
        rna_columns, rna = read_table(args.rna, (
            'genome', 'gtdb_domain', 'bakta_trna_unique_standard20',
            'trnascan_trna_unique_standard20', 'trna_complete',
            'merged_5S', 'merged_16S', 'merged_23S', 'rrna_complete'))
        _, genomes = read_table(args.manifest, ('genome', 'path'))
    except ValueError as error:
        ap.error(str(error))

    rna_fields = [name for name in rna_columns if name != 'genome']
    fields = ['genome','completeness','contamination'] + rna_fields + [
        'bin_status', 'hq_fail_reasons'
    ]
    statuses = ('HQ','MQ','LQ','Contaminated','UNASSESSED')
    counts = {status: 0 for status in statuses}

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter='\t'); w.writeheader()
        for gid in genomes:
            c, r = check.get(gid, {}), rna.get(gid, {})
            comp, cont = as_float(c.get('completeness')), as_float(c.get('contamination'))
            rrna_ok = all(b01(r.get(name)) for name in ('merged_5S', 'merged_16S', 'merged_23S'))
            bakta_types = as_int(r.get('bakta_trna_unique_standard20'))
            trnascan_types = as_int(r.get('trnascan_trna_unique_standard20'))
            trna_ok = bakta_types >= args.min_trna or trnascan_types >= args.min_trna
            status = classify(comp, cont, rrna_ok, trna_ok, args)
            counts[status] += 1

            w.writerow({
                'genome':gid, 'completeness':c.get('completeness',''), 'contamination':c.get('contamination',''),
                **{k:r.get(k,'') for k in rna_fields},
                'rrna_complete':int(rrna_ok), 'trna_complete':int(trna_ok),
                'bin_status':status,
                'hq_fail_reasons':hq_reasons(comp,cont,r,trna_ok,args),
            })

    Path(args.counts_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.counts_output, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['bin_status','n_genomes'])
        for status in statuses:
            w.writerow([status,counts[status]])


if __name__ == '__main__':
    main()
