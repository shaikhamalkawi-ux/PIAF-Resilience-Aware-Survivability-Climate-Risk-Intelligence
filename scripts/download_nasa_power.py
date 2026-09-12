#!/usr/bin/env python3
"""
PIAF / ECMX-D-26-01378
NASA POWER raw-data acquisition and reproducible climate-summary builder.

Downloads the exact 2014-2023 daily NASA POWER inputs used by the revised
photovoltaic screening study, preserves raw JSON responses, writes clean daily
CSV files, exact query URLs, SHA-256 manifests, site summaries, GHI seasonality
diagnostics, and a thermal-redundancy audit.

Interpretation boundary: raw daily GHI CV contains seasonal plus day-to-day
variation. Interannual and month-conditioned diagnostics are generated
separately so this is not misrepresented as pure stochastic volatility.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import requests

API_URL = "https://power.larc.nasa.gov/api/temporal/daily/point"
PARAMETERS = ["ALLSKY_SFC_SW_DWN", "T2M", "T2M_MAX", "RH2M", "WS10M"]
SITES = [
    {"location": "Abu Dhabi", "country": "United Arab Emirates", "latitude": 24.4539, "longitude": 54.3773},
    {"location": "NEOM Bay/Sharma", "country": "Saudi Arabia", "latitude": 28.1120, "longitude": 35.1940},
    {"location": "Jodhpur", "country": "India", "latitude": 26.2389, "longitude": 73.0243},
    {"location": "Las Vegas", "country": "United States", "latitude": 36.1699, "longitude": -115.1398},
    {"location": "Phoenix", "country": "United States", "latitude": 33.4484, "longitude": -112.0740},
    {"location": "Antofagasta", "country": "Chile", "latitude": -23.6509, "longitude": -70.3975},
    {"location": "Alice Springs", "country": "Australia", "latitude": -23.6980, "longitude": 133.8807},
    {"location": "Dunhuang", "country": "China", "latitude": 40.1421, "longitude": 94.6618},
    {"location": "Aswan", "country": "Egypt", "latitude": 24.0889, "longitude": 32.8998},
    {"location": "Ouarzazate", "country": "Morocco", "latitude": 30.9335, "longitude": -6.9370},
]
FILL_VALUES = {-999, -999.0, -99, -99.0}


def slugify(text: str) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", text.lower().strip())
    return text.strip("_")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def cv(series: pd.Series) -> float:
    x = pd.to_numeric(series, errors="coerce").dropna()
    if len(x) < 2:
        return float("nan")
    m = float(x.mean())
    if abs(m) < 1e-12:
        return float("nan")
    return float(x.std(ddof=1) / m)


def build_query(site: Dict, start: str, end: str) -> Dict[str, str]:
    return {
        "parameters": ",".join(PARAMETERS),
        "community": "RE",
        "longitude": str(site["longitude"]),
        "latitude": str(site["latitude"]),
        "start": start,
        "end": end,
        "format": "JSON",
    }


def request_with_retry(session: requests.Session, params: Dict[str, str], retries: int = 5, timeout: int = 120) -> requests.Response:
    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            r = session.get(API_URL, params=params, timeout=timeout)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                wait = min(60, 2 ** attempt)
                print(f"  transient HTTP {r.status_code}; retrying in {wait}s")
                time.sleep(wait)
                continue
            r.raise_for_status()
        except requests.RequestException as exc:
            last_exc = exc
            if attempt == retries:
                break
            wait = min(60, 2 ** attempt)
            print(f"  request error: {exc}; retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"NASA POWER request failed after {retries} attempts: {last_exc}")


def parse_power_json(payload: Dict) -> pd.DataFrame:
    try:
        block = payload["properties"]["parameter"]
    except KeyError as exc:
        raise ValueError("Unexpected NASA POWER JSON: properties.parameter not found") from exc
    all_dates = set()
    for p in PARAMETERS:
        if p not in block:
            raise ValueError(f"NASA POWER response missing parameter: {p}")
        all_dates.update(block[p].keys())
    dates = sorted(all_dates)
    data = {"date": pd.to_datetime(dates, format="%Y%m%d", errors="coerce")}
    for p in PARAMETERS:
        vals = []
        for d in dates:
            v = block[p].get(d, np.nan)
            if v in FILL_VALUES:
                v = np.nan
            vals.append(v)
        data[p] = pd.to_numeric(pd.Series(vals), errors="coerce")
    return pd.DataFrame(data).dropna(subset=["date"]).sort_values("date").reset_index(drop=True)


def ghi_seasonality_diagnostics(df: pd.DataFrame) -> Dict[str, float]:
    ghi = df[["date", "ALLSKY_SFC_SW_DWN"]].dropna().copy()
    ghi["year"] = ghi["date"].dt.year
    ghi["month"] = ghi["date"].dt.month
    annual_mean_daily = ghi.groupby("year")["ALLSKY_SFC_SW_DWN"].mean()
    monthly_year = ghi.groupby(["year", "month"])["ALLSKY_SFC_SW_DWN"].mean().rename("monthly_mean_daily_ghi").reset_index()
    month_cvs = monthly_year.groupby("month")["monthly_mean_daily_ghi"].apply(cv)
    return {
        "raw_daily_ghi_cv": cv(ghi["ALLSKY_SFC_SW_DWN"]),
        "interannual_annual_mean_ghi_cv": cv(annual_mean_daily),
        "mean_month_conditioned_interannual_cv": float(month_cvs.mean()),
        "max_month_conditioned_interannual_cv": float(month_cvs.max()),
        "min_month_conditioned_interannual_cv": float(month_cvs.min()),
        "n_years": int(annual_mean_daily.shape[0]),
        "n_valid_days_ghi": int(ghi.shape[0]),
    }


def summarize_site(site: Dict, df: pd.DataFrame) -> Dict[str, float]:
    return {
        "location": site["location"], "country": site["country"],
        "latitude": site["latitude"], "longitude": site["longitude"],
        "mean_ghi_kwh_m2_day": float(df["ALLSKY_SFC_SW_DWN"].mean()),
        "cv_ghi_raw_daily": cv(df["ALLSKY_SFC_SW_DWN"]),
        "mean_t2m_c": float(df["T2M"].mean()),
        "mean_t2m_max_c": float(df["T2M_MAX"].mean()),
        "mean_rh_pct": float(df["RH2M"].mean()),
        "mean_ws10m_ms": float(df["WS10M"].mean()),
        "cv_ws10m_raw_daily": cv(df["WS10M"]),
        "corr_t2m_t2m_max": float(df[["T2M", "T2M_MAX"]].corr().iloc[0, 1]),
        "n_rows": int(df.shape[0]),
    }


def build_url(params: Dict[str, str]) -> str:
    return requests.Request("GET", API_URL, params=params).prepare().url


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="20140101")
    ap.add_argument("--end", default="20231231")
    ap.add_argument("--output", default="NASA_POWER_ECMX_R1")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--pause", type=float, default=0.5)
    args = ap.parse_args()

    out = Path(args.output)
    raw_dir, daily_dir = out / "raw_json", out / "daily_csv"
    prov_dir, proc_dir = out / "provenance", out / "processed"
    for d in (raw_dir, daily_dir, prov_dir, proc_dir):
        d.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "ECMX-D-26-01378-PIAF-reproducibility/1.0"})
    summaries: List[Dict] = []
    seasonal_rows: List[Dict] = []
    thermal_rows: List[Dict] = []
    query_rows: List[Dict] = []
    manifest_rows: List[Dict] = []
    retrieved_at = datetime.now(timezone.utc).isoformat()

    for idx, site in enumerate(SITES, 1):
        slug = slugify(site["location"])
        raw_path = raw_dir / f"{slug}_{args.start}_{args.end}.json"
        daily_path = daily_dir / f"{slug}_{args.start}_{args.end}.csv"
        params = build_query(site, args.start, args.end)
        print(f"[{idx}/{len(SITES)}] {site['location']}")
        if raw_path.exists() and not args.overwrite:
            payload = json.loads(raw_path.read_text(encoding="utf-8"))
        else:
            r = request_with_retry(session, params)
            raw_path.write_bytes(r.content)
            payload = r.json()
            time.sleep(args.pause)
        df = parse_power_json(payload)
        df.insert(1, "location", site["location"])
        df.insert(2, "country", site["country"])
        df.insert(3, "latitude", site["latitude"])
        df.insert(4, "longitude", site["longitude"])
        df.to_csv(daily_path, index=False)
        summary = summarize_site(site, df)
        summaries.append(summary)
        seas = ghi_seasonality_diagnostics(df)
        seas.update({"location": site["location"], "country": site["country"], "latitude": site["latitude"], "longitude": site["longitude"]})
        seasonal_rows.append(seas)
        thermal_rows.append({"location": site["location"], "country": site["country"], "mean_t2m_c": summary["mean_t2m_c"], "mean_t2m_max_c": summary["mean_t2m_max_c"], "daily_corr_t2m_t2m_max": summary["corr_t2m_t2m_max"]})
        query_rows.append({"location": site["location"], "country": site["country"], "latitude": site["latitude"], "longitude": site["longitude"], "start": args.start, "end": args.end, "parameters": ",".join(PARAMETERS), "community": "RE", "format": "JSON", "url": build_url(params)})
        for kind, path in (("raw_json", raw_path), ("daily_csv", daily_path)):
            manifest_rows.append({"location": site["location"], "artifact_type": kind, "relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "bytes": path.stat().st_size, "retrieved_or_generated_utc": retrieved_at})

    summary_path = proc_dir / "sites_climate_summary.csv"
    seasonal_path = proc_dir / "ghi_seasonality_diagnostics.csv"
    thermal_path = proc_dir / "thermal_redundancy_audit.csv"
    query_path = prov_dir / "nasa_power_query_urls.csv"
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    pd.DataFrame(seasonal_rows).to_csv(seasonal_path, index=False)
    pd.DataFrame(thermal_rows).to_csv(thermal_path, index=False)
    pd.DataFrame(query_rows).to_csv(query_path, index=False)
    for kind, path in (("processed_summary", summary_path), ("seasonality_diagnostic", seasonal_path), ("thermal_redundancy_audit", thermal_path), ("query_registry", query_path)):
        manifest_rows.append({"location": "ALL", "artifact_type": kind, "relative_path": str(path.relative_to(out)), "sha256": sha256_file(path), "bytes": path.stat().st_size, "retrieved_or_generated_utc": retrieved_at})
    pd.DataFrame(manifest_rows).to_csv(prov_dir / "acquisition_manifest.csv", index=False)
    (prov_dir / "acquisition_metadata.json").write_text(json.dumps({"manuscript": "ECMX-D-26-01378", "project": "PIAF Resilience-Aware Survivability / Climate-Risk Intelligence", "retrieved_utc": retrieved_at, "n_sites": len(SITES), "period": {"start": args.start, "end": args.end}, "parameters": PARAMETERS, "api_endpoint": API_URL, "community": "RE", "interpretation_note": "Raw daily GHI CV includes seasonal and day-to-day variability; use the seasonality diagnostic for interannual and month-conditioned sensitivity."}, indent=2), encoding="utf-8")
    print(f"Done: {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
