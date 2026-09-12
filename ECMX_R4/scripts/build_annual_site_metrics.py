#!/usr/bin/env python3
"""Build annual ECMX-R4 site metrics from NASA POWER daily CSV files.

This closes the reproducibility step between raw/daily NASA POWER evidence and
`data/annual_site_metrics.csv`, the direct input to `run_analysis.py`.

Definitions (per location and calendar year):
- mean_ghi: arithmetic mean of daily ALLSKY_SFC_SW_DWN.
- heat_p95_tmax: 0.95 sample quantile of daily T2M_MAX using pandas' default
  linear interpolation.
- mean_rh: arithmetic mean of daily RH2M within the year.
- ghi_within_month_cv: arithmetic mean of the 12 within-month daily GHI
  coefficients of variation, each using sample SD (ddof=1) / monthly mean.
- wind_p95: 0.95 sample quantile of daily WS10M using linear interpolation.

Rows with missing values in a metric are excluded only from that metric's
calculation. The admitted ECMX-R4 snapshot has no missing daily values in the
five required variables, as recorded by QA.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

REQUIRED = ["date","location","country","latitude","longitude",
            "ALLSKY_SFC_SW_DWN","T2M_MAX","RH2M","WS10M"]


def cv_sample(s: pd.Series) -> float:
    x = pd.to_numeric(s, errors="coerce").dropna()
    if len(x) < 2:
        return float("nan")
    m = float(x.mean())
    if abs(m) < 1e-12:
        return float("nan")
    return float(x.std(ddof=1) / m)


def build_one(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["date"])
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing columns {missing}")
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month
    out = []
    for year, g in df.groupby("year", sort=True):
        monthly_cv = g.groupby("month")["ALLSKY_SFC_SW_DWN"].apply(cv_sample)
        out.append({
            "location": str(g["location"].iloc[0]),
            "country": str(g["country"].iloc[0]),
            "latitude": float(g["latitude"].iloc[0]),
            "longitude": float(g["longitude"].iloc[0]),
            "year": int(year),
            "mean_ghi": float(g["ALLSKY_SFC_SW_DWN"].mean()),
            "heat_p95_tmax": float(g["T2M_MAX"].quantile(0.95, interpolation="linear")),
            "mean_rh": float(g["RH2M"].mean()),
            "ghi_within_month_cv": float(monthly_cv.mean()),
            "wind_p95": float(g["WS10M"].quantile(0.95, interpolation="linear")),
        })
    return pd.DataFrame(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--daily-dir", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--seasonality-output", type=Path, default=None)
    args = ap.parse_args()
    files = sorted(args.daily_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {args.daily_dir}")
    result = pd.concat([build_one(p) for p in files], ignore_index=True)
    result = result.sort_values(["location","year"]).reset_index(drop=True)
    if result.location.nunique() != 10 or result.year.nunique() != 10 or len(result) != 100:
        raise ValueError(f"Expected 100 rows (10 sites x 10 years); got {len(result)}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    # Rebuild the site-level GHI diagnostics from the same daily files.
    diag_rows=[]
    for path in files:
        df=pd.read_csv(path,parse_dates=["date"])
        df["year"]=df["date"].dt.year; df["month"]=df["date"].dt.month
        ghi=df[["date","year","month","ALLSKY_SFC_SW_DWN"]].dropna().copy()
        annual=ghi.groupby("year")["ALLSKY_SFC_SW_DWN"].mean()
        ym=ghi.groupby(["year","month"])["ALLSKY_SFC_SW_DWN"].mean().reset_index()
        month_inter=ym.groupby("month")["ALLSKY_SFC_SW_DWN"].apply(cv_sample)
        within=ghi.groupby(["year","month"])["ALLSKY_SFC_SW_DWN"].apply(cv_sample)
        diag_rows.append({
            "raw_daily_ghi_cv":cv_sample(ghi["ALLSKY_SFC_SW_DWN"]),
            "interannual_annual_mean_ghi_cv":cv_sample(annual),
            "mean_month_conditioned_interannual_cv":float(month_inter.mean()),
            "max_month_conditioned_interannual_cv":float(month_inter.max()),
            "mean_within_month_daily_cv":float(within.mean()),
            "n_valid_ghi_days":int(len(ghi)),
            "n_years":int(len(annual)),
            "location":str(df["location"].iloc[0]),
            "country":str(df["country"].iloc[0]),
            "latitude":float(df["latitude"].iloc[0]),
            "longitude":float(df["longitude"].iloc[0]),
        })
    seasonality_path=args.seasonality_output or args.output.with_name("ghi_seasonality_diagnostics.csv")
    pd.DataFrame(diag_rows).sort_values("location").to_csv(seasonality_path,index=False)
    print(f"Wrote {args.output} ({len(result)} rows)")
    print(f"Wrote {seasonality_path} ({len(diag_rows)} rows)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
