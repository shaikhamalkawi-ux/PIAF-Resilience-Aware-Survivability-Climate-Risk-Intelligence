# Papers and data: acquisition/status list

This list separates sources that are openly accessible from sources for which a publisher/institution link may require library access. It is intended to support the ECMX-D-26-01378 revision only.

## Open / directly accessible

- Al-Ali, Olabi & Mahmoud (2024), UAE PV-park MCDM case study: https://doi.org/10.3390/en17174235
- Ranjgar, Niccolai & Leva (2026), GIS-MCDM collinearity/sensitivity: https://doi.org/10.1016/j.renene.2026.125873
- Abdallah, Ali & Kivambe (2023), five-year desert PV reliability: https://doi.org/10.1016/j.solener.2022.11.042
- Sen et al. (2024), TOPCon/PERC damp-heat failure modes: https://doi.org/10.1016/j.solmat.2024.112877
- Voyant et al. (2026), stochastic coefficient of variation for solar irradiance: https://doi.org/10.1016/j.renene.2025.123913
- Mazurek & Strzalka (2022), Monte Carlo weights in MCDM: https://doi.org/10.1371/journal.pone.0268950
- IEA PVPS Task 13 (2025), climate-specific PV optimisation: https://doi.org/10.69766/QSYC8858
- NASA POWER Solar Flux Assessment / Validation: https://power.larc.nasa.gov/docs/methodology/energy-fluxes/solar-validation/
- NASA POWER Daily API: https://power.larc.nasa.gov/docs/services/api/temporal/daily/

## Publisher / library access may be needed

- Almasad et al. (2023), Saudi PV site suitability: https://doi.org/10.1016/j.solener.2022.11.046
- Noorollahi et al. (2022), PV site selection and technical potential: https://doi.org/10.1016/j.renene.2021.12.124
- Habte et al. (2020), long-term solar resource variability: https://doi.org/10.1016/j.rser.2020.110285
- Nan et al. (2024), PERC/SHJ/TOPCon temperature dependence: https://doi.org/10.1016/j.solmat.2023.112649
- Adothu et al. (2024), desert PV reliability review: https://doi.org/10.1002/pip.3827
- Mehdi et al. (2023), hot-arid crystalline-Si vs CdTe experiment: https://doi.org/10.1016/j.renene.2023.01.082
- Kivambe et al. (2025), Qatar PERC/TOPCon/SHJ field comparison: https://doi.org/10.1016/j.solener.2025.113555

## Data

The project records all ten exact NASA POWER point-query URLs in `data/nasa_power_query_urls.csv`. `scripts/fetch_nasa_power_daily.py` retrieves the upstream daily CSVs and writes a SHA-256 acquisition manifest when run with internet access. The revision calculations use the committed summary snapshot in `data/sites_climate_summary.csv` to prevent live-provider updates from silently changing the submitted results.
