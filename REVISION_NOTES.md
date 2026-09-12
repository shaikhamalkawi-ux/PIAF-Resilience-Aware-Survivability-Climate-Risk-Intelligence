# ECMX-D-26-01378 - final reviewer implementation map (R3)

## Reviewer 3 major comments

1. **Technology comparability** - system-level alternatives were removed. PERC, TOPCon, SHJ/HJT and CdTe remain only as comparable PV-family evidence context; the old numerical family ranking is removed because a common four-family calibration boundary is unavailable.
2. **Resource magnitude missing** - long-term mean GHI is now explicit and remains analytically separate from climate context. The final screen is bi-objective rather than another weighted single score.
3. **Insufficient documentation / calibrated wording** - unsupported calibrated-coefficient language and technology constants were removed from the empirical calculation. All NASA-derived metrics and equations are defined explicitly.
4. **Benchmark-dependent normalization** - the 0.33/0.66 regime thresholds were removed. `C=0` is defined only as the lowest value within the ten-site benchmark. Leave-one-location-out normalization is reported; worst Spearman correlation is 0.900 and maximum rank shift is 2.
5. **MCDM reproducibility / Abu Dhabi focus** - WSM/TOPSIS/VIKOR are removed from the central analysis. All ten sites are treated under the same resource-climate Pareto screen.
6. **Narrow Monte Carlo** - the old +/-10% uniform perturbation is removed. Robustness now uses 10,000 empirical year-block bootstrap resamples, leave-one-year-out analysis, component omission and GHI-variability-definition sensitivity. No invented probability distribution is assigned to unsupported technology scores.
7. **Table/figure inconsistency** - final numerical tables and figures are regenerated from one R3 script and one committed R3 input set.
8. **Climate-risk/survivability claims** - journal-facing wording is narrowed to historical climate stress/context and pre-feasibility screening. No future climate, failure-probability, degradation-prediction, energy-yield, LCOE or bankability claim is made.
9. **Cross-sector transferability** - removed as an empirical claim and limited to future work.

## Minor comments

- Countries and coordinates are provided for every benchmark location.
- NASA POWER is described as satellite/model-derived analysis-ready evidence, not direct observation.
- Mean temperature is removed from the main climate-context calculation; heat uses the mean across years of the annual 95th percentile of daily maximum temperature.
- The baseline GHI variability term is mean within-month daily CV, reducing the deterministic annual seasonal cycle. Raw-daily and interannual alternatives are reported as sensitivity tests, and all preserve the same baseline Pareto set.
- Every revised table and figure is introduced and discussed.
- The old journal-name/title-page inconsistency is removed.
- Data/code availability points to the public repository and the exact immutable commit used for the resubmission.

## Reviewer 4

The Discussion now compares the framework with recent solar-site and climate-zone MCDM studies, NASA POWER benchmarking/solar-variability literature, and PV reliability evidence. The comparison is used to state both what the present screening adds and what it does not replace.

## R3 result

The baseline Pareto set is **Antofagasta and Aswan**. Both remain Pareto-nondominated in every leave-one-year-out analysis and all 10,000 year-block bootstrap resamples. The same two-site baseline Pareto set is retained under all four tested GHI-variability definitions.
