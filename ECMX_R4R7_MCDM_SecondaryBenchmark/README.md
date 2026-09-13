# ECMX R4R7 — Secondary all-site MCDM benchmark

This folder provides the reviewer-facing secondary WSM/TOPSIS/VIKOR benchmark for manuscript ECMX-D-26-01378. The primary scientific result remains the two-objective Pareto screen in the locked ECMX-R4 analysis.

## Scope
- Ten fixed benchmark locations only.
- Same two dimensions as the primary Pareto screen: mean GHI (benefit) and historical climate-context summary C (cost).
- Equal 0.5/0.5 resource-climate weighting is a declared diagnostic reference convention, not a calibrated physical weighting or elicited stakeholder preference.
- Full deterministic weight sensitivity over lambda in [0,1].
- Shared-minmax common-matrix WSM/TOPSIS/VIKOR benchmark plus conventional vector-normalized TOPSIS audit.
- Five-criterion result is a structural diagnostic only.

## Locked primary snapshot
Primary R4 analysis commit: `89d266e3ea18777ff131f3108690e4fdcd940601`.

## Main secondary result
At 0.5/0.5, WSM, TOPSIS, and VIKOR all rank Antofagasta first and Aswan second. VIKOR does not satisfy acceptable advantage because Q2-Q1 = 0.023615 < DQ = 0.111111, so both remain in the compromise set.

## Reproduce
```bash
python scripts/run_mcdm_benchmark.py
```

The script contains branch assertions for the standard VIKOR C1/C2 compromise rule. `run_log.txt` is intentionally excluded from the publishable manifest so shell redirection cannot create a self-inconsistent checksum record.

## Files
- `data/`: locked ten-site input table.
- `scripts/run_mcdm_benchmark.py`: executable benchmark and sensitivity analysis.
- `outputs/`: decision matrix, rankings, VIKOR test, weight sensitivity, normalization audit, and structural diagnostic.
- `ECMX_R4R7_MCDM_SecondaryBenchmark.xlsx`: formula-bearing audit workbook (included in the archival package).
- `Supplementary_Note_S1_MCDM.tex`: reviewer-facing supplementary note.
- `MCDM_Benchmark_Report.md`: concise technical report.
- `sha256_manifest.csv`: checksum register for the archival package.
