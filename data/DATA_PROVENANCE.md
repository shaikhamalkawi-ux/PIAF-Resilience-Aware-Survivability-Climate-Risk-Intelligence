# Climate-data provenance and reproducibility boundary

The analysis uses a committed 2014-2023 summary table (`sites_climate_summary.csv`) for ten fixed geographic points. That committed table is the numerical input used by `scripts/run_analysis.py` and is covered by the SHA-256 manifest.

For each point, `nasa_power_query_urls.csv` records the exact NASA POWER Daily API request for:

- `ALLSKY_SFC_SW_DWN`
- `T2M`
- `T2M_MAX`
- `RH2M`
- `WS10M`

The query period is 2014-01-01 through 2023-12-31.

The API URLs document the acquisition route to the daily source series. The committed summary is retained as the revision-stage analysis snapshot because upstream gridded products can be reprocessed or revised by data providers over time. The manuscript therefore distinguishes the version-controlled analysis input from the external live API.

NASA POWER solar products are satellite/model-derived products and should not be described as co-located site observations. NASA's own Solar Flux Assessment notes that high-quality surface measurements are generally more accurate than satellite-derived values, while also documenting validation of the satellite/model products.
