# Data dictionary

## `data/sites_climate_summary.csv`

Each row is one fixed benchmark point. The locations are geographic points and are **not** treated as nationally representative samples.

- `location`, `country`: benchmark-point label and country.
- `latitude`, `longitude`: NASA POWER point-query coordinates in decimal degrees.
- `mean_ghi_kwh_m2_day`: 2014-2023 daily mean of NASA POWER `ALLSKY_SFC_SW_DWN`.
- `cv_ghi`: sample coefficient of variation of the raw daily GHI series. The raw daily series was not deseasonalized, so this descriptor combines seasonal and day-to-day variation.
- `mean_t2m_c`: 2014-2023 daily mean of `T2M`.
- `mean_t2m_max_c`: 2014-2023 mean of daily `T2M_MAX`.
- `mean_rh_pct`: 2014-2023 daily mean of `RH2M`.
- `mean_ws10m_ms`: 2014-2023 daily mean of `WS10M`.
- `cv_ws10m`: sample coefficient of variation of the raw daily `WS10M` series; seasonality is retained.

## `data/pv_screening_components.csv`

Declared early-stage scenario inputs for the four comparable PV families. The values are retained as transparent screening assumptions and are **not** field-calibrated family constants.

Baseline-suitability descriptors:
- `E`: performance potential
- `M`: market/technology maturity
- `C`: cost attractiveness
- `OM`: operations-and-maintenance simplicity

Resilience descriptors:
- `RT`: thermal tolerance
- `RD`: degradation resistance
- `RH`: humidity/soiling tolerance
- `RO`: operational robustness

`S0` and `rho` are calculated by `scripts/run_analysis.py` from the declared descriptor values and the weights in `model_settings.csv`.

## `data/model_settings.csv`

Exact declared weights used to reconstruct baseline suitability, resilience, and historical climate context, plus the attenuation setting `lambda=1`. The settings define a transparent screening scenario and are not empirical damage coefficients.

## `data/technology_evidence_register.csv`

Literature-to-claim map recording why each source is relevant and what it does **not** establish. The literature constrains interpretation and plausibility; it does not numerically calibrate the scenario matrix.

## `data/nasa_power_query_urls.csv`

Exact NASA POWER point-query URLs, coordinates, period, variables, community, and output format used to document independent retrieval of the daily upstream data.

## `outputs/tables/table_9_descriptor_reversal_margins.csv`

For each benchmark point, the central winner and runner-up are identified. One corresponding 0-1 technology descriptor is decreased for the winner and increased for the runner-up by the same amount. The table reports the smallest feasible paired change that makes their central scores tie. This is a **local deterministic reversal margin**, not a statistical confidence interval or empirical probability distribution.
