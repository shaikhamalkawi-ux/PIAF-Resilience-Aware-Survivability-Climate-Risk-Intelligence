# ECMX-R4 data provenance

## Public source
NASA POWER Daily Point API, 2014-01-01 through 2023-12-31, for ten fixed geographic points.

Variables:
- `ALLSKY_SFC_SW_DWN`
- `T2M`
- `T2M_MAX`
- `RH2M`
- `WS10M`

The exact coordinates and query URLs are stored in `nasa_power_query_urls.csv`. The separately archived raw-evidence ZIP is `ECMX_D_26_01378_NASA_POWER_2014_2023_COMPLETE.zip`, SHA-256 `8dd9b7acfbe92ca83b2fb1c8c35271858444e6cf10f015f15fd04f837bfb57ec`.

## Reproducible transformation chain
1. `scripts/download_nasa_power.py` reacquires raw NASA POWER JSON and daily CSV files and records provenance.
2. `scripts/build_annual_site_metrics.py` transforms daily records into the 100 location-year rows used by the paper.
3. `scripts/run_analysis.py` regenerates all R4 numerical tables and publication figures from the committed annual metrics and seasonality diagnostics.
4. `scripts/run_all.py` wires the steps together. It can operate from an existing daily-data directory or reacquire NASA POWER first.

## Annual metric definitions
For each location and calendar year:
- `mean_ghi`: mean daily `ALLSKY_SFC_SW_DWN`.
- `heat_p95_tmax`: 0.95 quantile of daily `T2M_MAX` with linear interpolation.
- `mean_rh`: mean daily `RH2M`.
- `ghi_within_month_cv`: mean of the 12 within-month daily GHI coefficients of variation, each using sample SD (`ddof=1`) divided by the monthly mean.
- `wind_p95`: 0.95 quantile of daily `WS10M` with linear interpolation.

The admitted raw package contained 3,652 calendar days for every site and no missing/fill values in the five required variables. The transformation was replayed independently during the R4 editorial audit with numerical agreement to machine precision.

## Interpretation boundary
NASA POWER is used as a public, analysis-ready satellite/model/assimilation product rather than as a substitute for site instrumentation. The climate-context summary is benchmark-relative and descriptive. It is not a calibrated failure probability, damage function, energy-yield model, or climate projection.
