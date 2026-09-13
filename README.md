# PIAF Resilience-Aware Survivability / Climate-Risk Intelligence

This repository contains reproducibility materials for manuscript **ECMX-D-26-01378**, *A Resource-Explicit Historical Climate-Stress Framework for Pre-Feasibility Photovoltaic Infrastructure Screening*.

## Locked primary analysis

The primary scientific source of record remains **ECMX_R4/**. It closes the daily NASA POWER -> annual metrics -> final results path, reports temporal and structural sensitivity, and retains PV technology families as evidence context only.

Reproduce from an existing NASA daily-data directory:

```bash
pip install -r ECMX_R4/requirements.txt
python ECMX_R4/scripts/run_all.py --daily-dir /path/to/NASA_POWER_2014_2023/daily_csv
```

The R4 source snapshot contains the exact query registry, annual metrics, evidence registers, analysis code, tables, verification records, and SHA-256 manifest. The exact primary-analysis code snapshot is immutable Git commit:

`89d266e3ea18777ff131f3108690e4fdcd940601`

## Secondary reviewer benchmark

The R4R7 reviewer-closure adds a separately versioned all-site WSM/TOPSIS/VIKOR consistency benchmark under:

`ECMX_R4R7_MCDM_SecondaryBenchmark/`

It uses the same ten fixed sites and the same two decision dimensions as the primary Pareto screen: mean GHI as a benefit and historical climate-context summary C as a cost. Pareto nondominance remains the primary decision rule. The equal 0.5/0.5 resource-climate split is a declared diagnostic reference convention, not a calibrated physical weighting or elicited stakeholder preference.

At that reference, WSM, TOPSIS and VIKOR all rank Antofagasta first and Aswan second; VIKOR retains both in the compromise set. Deterministic weight and normalization audits show that the identity of a single compensatory winner is more specification-sensitive than the two-site shortlist.

The secondary folder contains the executable script, locked ten-site input table, reviewer-facing Supplementary Note S1, rank and normalization outputs, and a manifest verifier. Full grid outputs and the formula-bearing audit workbook are included in the archival R4R7 deposit package.

## Permanent archive

The locked ECMX-R4 reproducibility materials and reacquired NASA POWER raw-evidence archive are permanently archived in Zenodo:

**DOI: 10.5281/zenodo.22730737**

https://doi.org/10.5281/zenodo.22730737

The existing DOI identifies the published primary R4 archive. The R4R7 secondary MCDM benchmark is versioned separately in this repository so the locked primary scientific snapshot is not rewritten. A Zenodo new-version package is prepared separately for archival extension.

Earlier folders/files are retained for provenance but do not supersede the locked ECMX-R4 primary analysis.
