# ECMX-D-26-01378 — final major-revision reproducibility layer

Manuscript: **A Resource-Explicit Historical Climate-Stress Framework for Pre-Feasibility Photovoltaic Infrastructure Screening**

This folder contains the R3 analysis used for the final major revision of manuscript **ECMX-D-26-01378** submitted to *Energy Conversion and Management: X*.

## Scientific boundary

The empirical result uses NASA POWER-derived resource and climate evidence only. No numerical PV-family score, market coefficient, or arbitrary perturbation distribution enters the main calculation.

The site screen keeps resource magnitude and historical climate context separate:

- Resource: 2014-2023 mean daily GHI.
- Heat proxy: mean across years of the annual 95th percentile of daily T2M_MAX.
- Humidity proxy: mean daily RH2M.
- Solar-variability proxy: mean within-month daily GHI coefficient of variation.
- Wind-exposure proxy: mean across years of the annual 95th percentile of daily WS10M.

The four climate proxies are min-max normalized inside the ten-site benchmark and averaged with equal weights. Equal weighting is a neutral summary convention, not a calibrated physical importance model. Component values remain visible and component-omission sensitivity is reported.

The final screen is bi-objective: maximize mean GHI and minimize the historical climate-context summary. Pareto nondominance is used instead of a single weighted technology score.

## Main result

Aswan and Antofagasta form the baseline resource-climate Pareto set. Both remain Pareto-nondominated in every leave-one-year-out analysis and in all 10,000 empirical year-block bootstrap resamples. Alternative GHI-variability definitions preserve the same Pareto set.

Four comparable PV families — PERC, TOPCon, SHJ/HJT, and CdTe — remain in an evidence-only table and discussion. No universal family coefficient or quantitative family ranking is claimed because a common calibration boundary is not available across all four families.

## Reproducibility

The exact NASA POWER acquisition queries remain in the repository root at `data/nasa_power_query_urls.csv`. The raw evidence package used for this revision has SHA-256:

`8dd9b7acfbe92ca83b2fb1c8c35271858444e6cf10f015f15fd04f837bfb57ec`

Run `ECMX_R3/scripts/run_analysis.py` with the R3 input files to regenerate the final tables and figures.

This is historical pre-feasibility screening only. It is not a climate projection, field-failure probability model, degradation prediction, energy-yield model, LCOE model, or bankability assessment.
