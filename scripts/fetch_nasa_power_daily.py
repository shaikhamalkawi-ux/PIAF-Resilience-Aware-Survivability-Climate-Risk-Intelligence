#!/usr/bin/env python3
"""Fetch the NASA POWER daily point files listed in data/nasa_power_query_urls.csv.

This script is intentionally separate from run_analysis.py. It requires internet
access and reproduces the upstream daily source acquisition. The manuscript's
version-controlled analysis uses the committed climate-summary table so that
later updates to the live NASA service cannot silently change the revision.
"""
from pathlib import Path
import hashlib
import re
import urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
q = pd.read_csv(ROOT / "data" / "nasa_power_query_urls.csv")
out = ROOT / "data" / "raw_nasa_power"
out.mkdir(parents=True, exist_ok=True)
rows=[]
for _, r in q.iterrows():
    slug=re.sub(r"[^A-Za-z0-9]+","_",r.location).strip("_").lower()
    path=out/f"{slug}_2014_2023.csv"
    print("Fetching", r.location)
    with urllib.request.urlopen(r.url, timeout=120) as resp:
        data=resp.read()
    path.write_bytes(data)
    rows.append({"location":r.location,"file":str(path.relative_to(ROOT)),"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),"url":r.url})
pd.DataFrame(rows).to_csv(ROOT/"data"/"raw_nasa_power_manifest.csv",index=False)
print("Done")
