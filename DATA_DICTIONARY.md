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

## `data/pv_screening_components.csv`

Declared early-stage scenario inputs for the four comparable PV families. The component values are retained from the submitted screening design so that the revised calculations are reproducible. They are **not** field-calibrated family constants.

Baseline-suitability components: performance (`E`), market maturity (`M`), cost attractiveness (`C`), dispatchability/flexibility (`D`), and O&M simplicity (`OM`).

Resilience components: thermal tolerance (`RT`), degradation resistance (`RD`), humidity/soiling tolerance (`RH`), and operational robustness (`RO`).

## `data/model_settings.csv`

Exact declared weights used to reconstruct baseline suitability, resilience, and climate context, plus the baseline attenuation setting `lambda=1`. These values define a transparent screening scenario and are not empirical damage coefficients.

## `data/technology_evidence_register.csv`

Literature-to-claim map recording the role of each source and what it does not establish.

## `data/nasa_power_query_urls.csv`

Exact point-query URLs, coordinates, period, variables, community, and output format used for independent NASA POWER retrieval.
