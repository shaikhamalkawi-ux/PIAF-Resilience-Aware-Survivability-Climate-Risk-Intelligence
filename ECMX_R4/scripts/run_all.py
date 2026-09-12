#!/usr/bin/env python3
"""End-to-end ECMX-R4 workflow.

Modes:
  1) With --daily-dir: rebuild annual metrics from an existing NASA POWER daily
     directory and run the analysis without network access.
  2) With --acquire: call download_nasa_power.py first, then rebuild annual
     metrics from the acquired daily files and run the analysis.
"""
from __future__ import annotations
import argparse, shutil, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def run(cmd):
    print('+',' '.join(map(str,cmd)),flush=True)
    subprocess.run([str(x) for x in cmd],check=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--daily-dir',type=Path); ap.add_argument('--acquire',action='store_true'); ap.add_argument('--acquisition-dir',type=Path,default=ROOT/'acquired_nasa_power')
    a=ap.parse_args()
    daily=a.daily_dir
    if a.acquire:
        run([sys.executable,ROOT/'scripts/download_nasa_power.py','--output',a.acquisition_dir,'--overwrite'])
        daily=a.acquisition_dir/'daily_csv'
        # carry acquisition provenance into the R4 data folder
        for src,dst in [
            (a.acquisition_dir/'provenance/nasa_power_query_urls.csv',ROOT/'data/nasa_power_query_urls.csv'),
            (a.acquisition_dir/'provenance/acquisition_metadata.json',ROOT/'data/acquisition_metadata.json'),
            (a.acquisition_dir/'provenance/qa_checks.csv',ROOT/'data/qa_checks.csv'),
            (a.acquisition_dir/'processed/ghi_seasonality_diagnostics.csv',ROOT/'data/ghi_seasonality_diagnostics.csv')]:
            if src.exists(): shutil.copy2(src,dst)
    if daily is None:
        raise SystemExit('Provide --daily-dir PATH or use --acquire')
    run([sys.executable,ROOT/'scripts/build_annual_site_metrics.py','--daily-dir',daily,'--output',ROOT/'data/annual_site_metrics.csv','--seasonality-output',ROOT/'data/ghi_seasonality_diagnostics.csv'])
    run([sys.executable,ROOT/'scripts/run_analysis.py'])

if __name__=='__main__': main()
