# ECMX R4R7 Secondary MCDM Benchmark Report

## Decision boundary
This is a secondary methodological benchmark only. Pareto nondominance remains the primary decision rule. The benchmark does not change any R4 site, NASA POWER input, climate-context calculation, Pareto set, bootstrap result, or conclusion.

## Equal-weight result
At the declared diagnostic reference weighting (resource 0.5, climate context 0.5):
- WSM: Antofagasta #1, Aswan #2.
- TOPSIS: Antofagasta #1, Aswan #2.
- VIKOR: Antofagasta #1, Aswan #2.
- Spearman: WSM-TOPSIS 0.987879; WSM-VIKOR 0.975758; TOPSIS-VIKOR 0.963636.
- VIKOR DQ = 0.111111; Q2-Q1 = 0.023615; C1=False, C2=True; compromise set = Antofagasta + Aswan.

## Weight sensitivity
Shared-minmax common-matrix rank-one switches occur on the 0.001 grid at:
- WSM: Antofagasta through 0.537; Aswan from 0.538.
- TOPSIS: Antofagasta through 0.534; Aswan from 0.535.
- VIKOR: Antofagasta through 0.537; Aswan from 0.538.
The exact Pareto pair remains the top two through lambda=0.977, 0.904, and 0.954, respectively.

## Normalization audit
Conventional vector-normalized TOPSIS also ranks Antofagasta first and Aswan second at 0.5/0.5. Its rank-one switch occurs between lambda=0.732 and 0.733, confirming that precise scalar-switch thresholds are normalization-dependent.

## Five-criterion structural diagnostic
When the climate summary is represented by its four existing normalized components (resource weight 0.5; each climate component 0.125), all three methods rank Aswan first and Antofagasta second. The change in scalar winner arises from representation and normalization geometry, not from new evidence.

## QA repairs in R4R7
- Standard three-branch VIKOR compromise rule implemented, including the C1=True/C2=False branch.
- VIKOR regret notation changed in the Supplement to `\mathcal{R}^{\mathrm{VIKOR}}_l` to avoid collision with the manuscript resource symbol `R_l`.
- Full TOPSIS distances and VIKOR S/regret/Q equations included.
- 0.5/0.5 weights explicitly labeled diagnostic, not calibrated.
- Conventional normalization audit made explicit.
- `run_log.txt` removed from the publishable manifest to eliminate redirect-time checksum inconsistency.
