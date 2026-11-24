#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
import re
from pathlib import Path
from urllib.parse import unquote

STANDARD_AA = {
    'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile',
    'Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val'
}
AA_CODES = {
    'A':'Ala', 'R':'Arg', 'N':'Asn', 'D':'Asp', 'C':'Cys',
    'Q':'Gln', 'E':'Glu', 'G':'Gly', 'H':'His', 'I':'Ile',
    'L':'Leu', 'K':'Lys', 'M':'Met', 'F':'Phe', 'P':'Pro',
    'S':'Ser', 'T':'Thr', 'W':'Trp', 'Y':'Tyr', 'V':'Val'
}
AA_LABELS = {aa.lower(): aa for aa in STANDARD_AA}
AA_LABELS.update({code.lower(): aa for code, aa in AA_CODES.items()})
AA_LABELS.update({'imet':'Met', 'fmet':'Met', 'ini':'Met', 'im':'Met', 'fm':'Met', 'ile2':'Ile'})


def standard_amino_acid(label, gene=False):
    """Normalize a tRNA type, Bakta gene or product to one of 20 amino acids."""
    prefix = r'^trn' if gene else r'^tRNA[-_ ]+'
    label = re.sub(prefix, '', label.strip(), flags=re.I)
    label = re.split(r'[(-]', label, maxsplit=1)[0].strip()
    return AA_LABELS.get(label.lower())


def read_manifest(path):
    d = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[row['genome']] = row['path']
    return d


def read_domain(path):
    d = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            d[row['genome']] = row['domain']
    return d


def attrs(s):
    out = {}
    for item in s.strip().split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            out[k] = unquote(v)
    return out


def rrna_class(text):
    t = text.lower().replace('_', ' ')
    if re.search(r'(^|\W)5s(\W|$)', t): return '5S'
    if re.search(r'(^|\W)16s(\W|$)', t): return '16S'
    if re.search(r'(^|\W)23s(\W|$)', t): return '23S'
    return None


def open_text(path):
    return gzip.open(path, 'rt', errors='replace') if str(path).endswith('.gz') else open(path, 'rt', errors='replace')


def fasta_lengths(path):
    lengths, name, n = {}, None, 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = n
                name = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
        if name is not None:
            lengths[name] = n
    return lengths


def edge_fields(start, end, contig_len, margin):
    if contig_len in (None, 0):
        return '', '', '', ''
    s, e = min(start, end), max(start, end)
    left = max(0, s - 1)
    right = max(0, contig_len - e)
    at_edge = int(left == 0 or right == 0)
    near = int(min(left, right) <= margin)
    return left, right, at_edge, near


def parse_bakta_tsv(path):
    counts = {'5S': 0, '16S': 0, '23S': 0}
    trna_total, raw_types, standard_types = 0, set(), set()
    if not path.exists():
        return None
    with open(path, errors='replace') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 7:
                continue
            ftype, gene = p[1], p[6]
            if ftype == 'tRNA':
                trna_total += 1
                product = p[7] if len(p) > 7 else ''
                if gene and gene != '-':
                    raw_types.add(gene)
                aa = standard_amino_acid(gene, gene=True) or standard_amino_acid(product)
                if aa:
                    standard_types.add(aa)
            elif ftype == 'rRNA':
                if gene == 'rrs': counts['16S'] += 1
                elif gene == 'rrl': counts['23S'] += 1
                elif gene == 'rrf': counts['5S'] += 1
    return counts, trna_total, raw_types, standard_types


def bakta_hits(json_path, genome, domain, fallback_lengths, margin):
    hits = []
    if not json_path.exists():
        return hits
    with open(json_path, errors='replace') as fh:
        data = json.load(fh)
    lens = dict(fallback_lengths)
    for seq in data.get('sequences', []):
        sid = seq.get('id') or seq.get('name')
        if sid and seq.get('length') is not None:
            lens[sid] = int(seq['length'])
    for feat in data.get('features', []):
        if str(feat.get('type', '')).lower() != 'rrna':
            continue
        cls = {'rrf':'5S','rrs':'16S','rrl':'23S'}.get(feat.get('gene')) or rrna_class(feat.get('product',''))
        if not cls:
            continue
        contig = feat.get('sequence') or feat.get('contig') or ''
        start, end = int(feat['start']), int(feat.get('stop', feat.get('end')))
        clen = lens.get(contig)
        left, right, edge, near = edge_fields(start, end, clen, margin)
        trunc_detail = feat.get('truncated', '')
        trunc_flag = str(trunc_detail).strip().lower() not in {'', 'none', 'false', '0'}
        length_bp = len(feat.get('nt','')) or abs(end - start) + 1
        hits.append({
            'genome': genome, 'gtdb_domain': domain, 'predictor': 'Bakta', 'model': 'Infernal/Rfam',
            'rrna_type': cls, 'contig': contig, 'start': start, 'end': end,
            'strand': feat.get('strand',''), 'length_bp': length_bp, 'contig_length_bp': clen or '',
            'distance_left_edge_bp': left, 'distance_right_edge_bp': right,
            'contig_edge': edge, 'near_contig_edge': near,
            'truncated': int(trunc_flag), 'truncated_detail': trunc_detail,
            'partial': int(trunc_flag), 'coverage': feat.get('coverage',''),
            'score': feat.get('score',''), 'evalue': feat.get('evalue','')
        })
    return hits


def parse_barrnap_file(path, genome, domain, contig_lengths, margin):
    model = Path(path).name.split('.')[-2]
    counts = {'5S': 0, '16S': 0, '23S': 0}
    hits = []
    with open(path, errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2].lower() != 'rrna':
                continue
            a = attrs(p[8])
            text = ' '.join([a.get('Name',''), a.get('product',''), p[8]])
            cls = rrna_class(text)
            if not cls:
                continue
            counts[cls] += 1
            start, end = int(p[3]), int(p[4])
            contig = p[0]
            clen = contig_lengths.get(contig)
            left, right, edge, near = edge_fields(start, end, clen, margin)
            partial = int('partial' in text.lower())
            hits.append({
                'genome': genome, 'gtdb_domain': domain, 'predictor': 'Barrnap', 'model': model,
                'rrna_type': cls, 'contig': contig, 'start': start, 'end': end, 'strand': p[6],
                'length_bp': abs(end-start)+1, 'contig_length_bp': clen or '',
                'distance_left_edge_bp': left, 'distance_right_edge_bp': right,
                'contig_edge': edge, 'near_contig_edge': near,
                'truncated': partial, 'truncated_detail': 'partial' if partial else '',
                'partial': partial, 'coverage': '', 'score': '', 'evalue': '' if p[5] == '.' else p[5]
            })
    return counts, hits, model


def parse_barrnap(paths, genome, domain, contig_lengths, margin):
    results, all_hits = [], []
    for p in paths:
        if not p.exists():
            continue
        c, h, model = parse_barrnap_file(p, genome, domain, contig_lengths, margin)
        results.append((model, c))
        all_hits.extend(h)
    if not results:
        return {'5S':0,'16S':0,'23S':0}, '', all_hits
    # Unknown-domain fallback may have two models. Keep raw hit sets for QC but
    # take the per-rRNA maximum count so the same biological locus is not added twice.
    merged = {k:max(r[1][k] for r in results) for k in ('5S','16S','23S')}
    selected = results[0][0] if len(results) == 1 else 'max(bac,arc)'
    return merged, selected, all_hits


def parse_trnascan_file(path):
    total, uniq_raw, uniq_standard = 0, set(), set()
    with open(path, errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('Sequence') or set(s) <= {'-',' '}:
                continue
            p = s.split()
            if len(p) < 6 or not p[1].isdigit():
                continue
            total += 1
            typ = p[4]
            if typ not in {'Undet','Sup'}:
                uniq_raw.add(typ)
            aa = standard_amino_acid(typ)
            if aa:
                uniq_standard.add(aa)
    return total, uniq_raw, uniq_standard


def parse_trnascan(paths):
    rows = []
    for p in paths:
        if not p.exists():
            continue
        label = Path(p).name.split('.')[-2]
        total, raw, standard = parse_trnascan_file(p)
        rows.append((label,total,raw,standard))
    if not rows:
        return '', 0, set(), set()
    if len(rows) == 1:
        return rows[0]
    # Unknown-domain fallback only. Prefer the model covering the largest number
    # of the standard 20 amino acids; raw labels (e.g. Met/iMet) are secondary.
    best = max(rows, key=lambda x:(len(x[3]), len(x[2]), x[1]))
    label, total, raw, standard = best
    return f'max(bac,arc)->{label}', total, raw, standard


def selected_barrnap_paths(root, gid, domain):
    d = root / 'barrnap'
    if domain == 'Bacteria':
        labels = ['bac']
    elif domain == 'Archaea':
        labels = ['arc']
    else:
        labels = ['bac', 'arc']
    return [d / f'{gid}.{label}.gff3' for label in labels]


def selected_trnascan_paths(root, gid, domain):
    d = root / 'trnascan'
    if domain == 'Bacteria':
        labels = ['bac']
    elif domain == 'Archaea':
        labels = ['arc']
    else:
        labels = ['bac', 'arc']
    return [d / f'{gid}.{label}.tsv' for label in labels]


def yn(x): return 1 if x else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--domain-map', required=True)
    ap.add_argument('--outroot', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--hits-output', required=True)
    ap.add_argument('--min-trna', type=int, default=18)
    ap.add_argument('--edge-margin', type=int, default=10)
    args = ap.parse_args()
    if not 1 <= args.min_trna <= 20:
        ap.error('--min-trna must be between 1 and 20 amino-acid types')
    if args.edge_margin < 0:
        ap.error('--edge-margin must be nonnegative')

    genomes = read_manifest(args.manifest)
    domains = read_domain(args.domain_map)
    root = Path(args.outroot)

    fields = [
        'genome','gtdb_domain',
        'bakta_trna_total','bakta_trna_unique_raw','bakta_trna_unique_standard20',
        'bakta_16s','bakta_23s','bakta_5s', 'barrnap_16s','barrnap_23s','barrnap_5s',
        'trnascan_trna_total','trnascan_trna_unique_raw','trnascan_trna_unique_standard20',
        'barrnap_model_used','trnascan_model_used',
        'merged_5S','merged_16S','merged_23S','rrna_complete',
        'trna_complete','rrna_hit_count','rrna_truncated_hit_count','rrna_contig_edge_hit_count'
    ]
    hit_fields = [
        'genome','gtdb_domain','predictor','model','rrna_type','contig','start','end','strand',
        'length_bp','contig_length_bp','distance_left_edge_bp','distance_right_edge_bp',
        'contig_edge','near_contig_edge','truncated','truncated_detail','partial','coverage','score','evalue'
    ]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.hits_output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w', newline='') as out, open(args.hits_output, 'w', newline='') as hout:
        w = csv.DictWriter(out, fieldnames=fields, delimiter='\t'); w.writeheader()
        hw = csv.DictWriter(hout, fieldnames=hit_fields, delimiter='\t', extrasaction='ignore'); hw.writeheader()

        for gid, fasta in genomes.items():
            domain = domains.get(gid, 'Unknown')
            lengths = fasta_lengths(fasta)

            btsv = root / 'bakta' / gid / f'{gid}.tsv'
            bj = root / 'bakta' / gid / f'{gid}.json'
            bres = parse_bakta_tsv(btsv)
            if bres is None:
                bc = {'5S':0,'16S':0,'23S':0}; btotal = 0; braw = set(); bstandard = set()
            else:
                bc, btotal, braw, bstandard = bres
            bhits = bakta_hits(bj, gid, domain, lengths, args.edge_margin)

            # IMPORTANT: do not glob every historical model output. Select only
            # the model(s) appropriate for the current GTDB domain.
            barr_paths = selected_barrnap_paths(root, gid, domain)
            br, br_model, br_hits = parse_barrnap(
                barr_paths, gid, domain, lengths, args.edge_margin
            )

            trna_paths = selected_trnascan_paths(root, gid, domain)
            trna_model, ttotal, traw, tstandard = parse_trnascan(trna_paths)

            all_hits = bhits + br_hits
            for hit in all_hits:
                hw.writerow(hit)

            m5 = bc['5S'] >= 1 or br['5S'] >= 1
            m16 = bc['16S'] >= 1 or br['16S'] >= 1
            m23 = bc['23S'] >= 1 or br['23S'] >= 1
            rrna_ok = m5 and m16 and m23

            # MIMAG means tRNAs for >=18 of the 20 standard amino acids.
            # Do not let raw labels such as Met + iMet count as two amino acids.
            trna_ok = len(bstandard) >= args.min_trna or len(tstandard) >= args.min_trna

            w.writerow({
                'genome':gid, 'gtdb_domain':domain,
                'bakta_trna_total':btotal, 'bakta_trna_unique_raw':len(braw),
                'bakta_trna_unique_standard20':len(bstandard),
                'bakta_16s':bc['16S'], 'bakta_23s':bc['23S'], 'bakta_5s':bc['5S'],
                'barrnap_16s':br['16S'], 'barrnap_23s':br['23S'], 'barrnap_5s':br['5S'],
                'trnascan_trna_total':ttotal,
                'trnascan_trna_unique_raw':len(traw),
                'trnascan_trna_unique_standard20':len(tstandard),
                'barrnap_model_used':br_model, 'trnascan_model_used':trna_model,
                'merged_5S':yn(m5), 'merged_16S':yn(m16), 'merged_23S':yn(m23),
                'rrna_complete':yn(rrna_ok), 'trna_complete':yn(trna_ok),
                'rrna_hit_count':len(all_hits),
                'rrna_truncated_hit_count':sum(int(h['truncated']) for h in all_hits),
                'rrna_contig_edge_hit_count':sum(int(h['contig_edge']) for h in all_hits if h['contig_edge'] != ''),
            })


if __name__ == '__main__':
    main()
