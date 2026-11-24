#!/usr/bin/env python3
import argparse, csv, re
from pathlib import Path

SUFFIXES = ('.fna.gz','.fa.gz','.fasta.gz','.fna','.fa','.fasta','.gz')

def norm(name):
    s = Path(name).name
    for x in SUFFIXES:
        if s.endswith(x):
            s = s[:-len(x)]; break
    return re.sub(r'[^A-Za-z0-9_.-]+','_',s).strip('_')


def resolve_genome(name, genome_ids):
    base = Path(name).name
    safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', base).strip('_')
    # A report can already contain the manifest ID, including a terminal .fa.
    for candidate in (name, base, safe, norm(name)):
        if candidate in genome_ids:
            return candidate
    raise ValueError(f'CheckM2 report genome {name!r} does not match the input manifest')


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',required=True); ap.add_argument('--output',required=True); ap.add_argument('--manifest',required=True); a=ap.parse_args()
    with open(a.manifest, newline='') as fh:
        genome_ids = {row['genome'] for row in csv.DictReader(fh, delimiter='\t')}
    rows=[]
    seen=set()
    for p in sorted(Path(a.root).rglob('quality_report.tsv')):
        with open(p,newline='') as fh:
            r=csv.DictReader(fh,delimiter='\t')
            if not r.fieldnames: continue
            low={x.lower():x for x in r.fieldnames}
            name=low.get('name') or low.get('genome') or r.fieldnames[0]
            comp=next((x for x in r.fieldnames if x.lower()=='completeness'),None)
            cont=next((x for x in r.fieldnames if x.lower()=='contamination'),None)
            if not comp or not cont: raise SystemExit(f'Cannot find completeness/contamination in {p}')
            for q in r:
                try:
                    genome = resolve_genome(q[name].strip(), genome_ids)
                except ValueError as error:
                    ap.error(str(error))
                if genome in seen:
                    ap.error(f'Duplicate CheckM2 report for {genome!r}')
                seen.add(genome)
                rows.append({'genome':genome,'completeness':q[comp],'contamination':q[cont]})
    rows.sort(key=lambda x:x['genome'])
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with open(out,'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['genome','completeness','contamination'],delimiter='\t'); w.writeheader(); w.writerows(rows)
if __name__=='__main__': main()
