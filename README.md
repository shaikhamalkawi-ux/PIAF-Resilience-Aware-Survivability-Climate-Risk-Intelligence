# ECMX-D-26-01378 revision reproducibility package

This public repository supports the major revision of:

**A Resource-Explicit Historical Climate-Stress Framework for Pre-Feasibility Photovoltaic Infrastructure Screening**

submitted to *Energy Conversion and Management: X* (ECMX-D-26-01378).

## Final revision layer

The journal-facing final revision is in [`ECMX_R3/`](ECMX_R3/).

R3 was rebuilt to answer the reviewers without retaining unsupported numerical technology-family coefficients. The main empirical analysis uses only NASA POWER-derived historical resource and climate evidence.

### R3 design

- Ten fixed geographic benchmark points; NASA POWER daily data for 2014-2023.
- Resource magnitude is explicit: long-term mean daily GHI.
- Climate context uses four visible proxies: annual p95 daily maximum temperature, mean RH2M, within-month daily GHI coefficient of variation, and annual p95 daily WS10M.
- The four proxies are normalized within the ten-site benchmark and averaged equally. Equal weighting is a neutral summary convention, not physical calibration.
- The final screen is bi-objective: maximize mean GHI and minimize historical climate context.
- Pareto nondominance replaces the old single weighted technology ranking.
- Robustness uses 10,000 empirical year-block bootstrap resamples, leave-one-year-out tests, leave-one-location-out normalization, component omission, and alternative GHI-variability definitions.
- PERC, TOPCon, SHJ/HJT and CdTe are retained only in an evidence-boundary table and discussion. No universal family coefficient or family ranking enters the main calculation.

### Main result

Aswan and Antofagasta are the baseline resource-climate Pareto set. Both remain Pareto-nondominated in every leave-one-year-out analysis and in all 10,000 empirical year-block bootstrap resamples. Alternative GHI-variability definitions preserve the same baseline Pareto set.

## Reproduce R3

```bash
pip install numpy pandas matplotlib
python ECMX_R3/scripts/run_analysis.py
```

The script regenerates the numerical tables, figures, analysis summary and SHA-256 manifest from the committed R3 data.

## Data provenance

Exact NASA POWER acquisition queries are preserved in the repository and in the revision transfer package. The complete reacquired raw evidence package used for the final revision has SHA-256:

`8dd9b7acfbe92ca83b2fb1c8c35271858444e6cf10f015f15fd04f837bfb57ec`

NASA POWER is treated as satellite/model-derived analysis-ready evidence, not direct site observation.

## Interpretation boundary

This is historical comparative pre-feasibility screening. It is not a climate projection, field-failure probability model, degradation prediction, energy-yield model, LCOE model, or bankability assessment. A benchmark Pareto result is not a claim of universal site or PV-family superiority.

Earlier root-level ECMX-R2 materials are retained only as revision history. **Use `ECMX_R3/` for the submitted major revision.**
