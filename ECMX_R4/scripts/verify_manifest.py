#!/usr/bin/env python3
"""Verify outputs/sha256_manifest.csv against the current ECMX-R4 tree."""
from pathlib import Path
import csv, hashlib, sys
ROOT=Path(__file__).resolve().parents[1]
MAN=ROOT/'outputs/sha256_manifest.csv'

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

rows=list(csv.DictReader(MAN.open(encoding='utf-8')))
fail=[]
for r in rows:
    p=ROOT/r['relative_path']
    if not p.is_file():
        fail.append((r['relative_path'],'MISSING',''))
        continue
    got=sha256(p)
    if got!=r['sha256']:
        fail.append((r['relative_path'],r['sha256'],got))
print(f'Checked {len(rows)} manifested files.')
if fail:
    for x in fail: print('FAIL',*x)
    raise SystemExit(1)
print('PASS: manifest matches the ECMX-R4 tree.')
