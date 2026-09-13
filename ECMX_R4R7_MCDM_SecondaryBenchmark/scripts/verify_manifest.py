#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, sys
ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'sha256_manifest.csv'
rows=list(csv.DictReader(manifest.open(encoding='utf-8')))
bad=[]
for r in rows:
    p=ROOT/r['relative_path']
    if not p.exists():
        bad.append((r['relative_path'],'missing')); continue
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h!=r['sha256'] or p.stat().st_size!=int(r['bytes']):
        bad.append((r['relative_path'],'mismatch'))
print(f'checked={len(rows)} bad={len(bad)}')
for x in bad: print(x)
sys.exit(1 if bad else 0)
