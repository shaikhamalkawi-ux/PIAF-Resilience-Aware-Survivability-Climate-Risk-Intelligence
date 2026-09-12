# Data dictionary

## `data/sites_climate_summary.csv`

Each row is one fixed benchmark point. The locations are geographic points and are **not** treated as nationally representative samples.

- `location`, `country`: benchmark-point label and country.
- `latitude`, `longitude`: query coordinates in decimal degrees.
- `mean_ghi_kwh_m2_day`: 2014–2023 daily mean of NASA POWER `ALLSKY_SFC_SW_DWN`.
- `cv_ghi`: sample coefficient of variation of the raw daily GHI series.
- `mean_t2m_c`: 2014–2023 daily mean of `T2M`.
- `mean_t2m_max_c`: 2014–2023 mean of daily `T2M_MAX`.
- `mean_rh_pct`: 2014–2023 daily mean of `RH2M`.
- `mean_ws10m_ms`: 2014–2023 daily mean of `WS10M`.
- `cv_ws10m`: sample coefficient of variation of the raw daily `WS10M` series.

## `data/pv_scenario_intervals.csv`

Declared early-stage screening envelopes for four PV technology families.

- `s0_lower`, `s0_central`, `s0_upper`: baseline-suitability screening values.
- `rho_lower`, `rho_central`, `rho_upper`: resilience screening values.
- `evidence_status`: qualitative evidence class/claim status.
- `interpretation`: explicit claim boundary.

These values are **not** field-calibrated probability distributions and the lower/upper bounds are **not** confidence intervals.

## `data/technology_evidence_register.csv`

Literature-to-claim map recording the role of each source and what it does not establish.

## `data/nasa_power_query_urls.csv`

Exact point-query URLs, coordinates, period, variables, community, and output format used for independent NASA POWER retrieval.
