# Major-revision notes — ECMX-D-26-01378

These notes summarize the methodological changes made after peer review. They are provided for transparency; the formal point-by-point response is submitted separately through Elsevier.

## Scope

The initial submission mixed PV module families with system-level configurations. The revised quantitative analysis is restricted to four comparable PV families: PERC, TOPCon, SHJ/HJT, and CdTe.

## Resource term

Mean solar resource is now explicit in the final screening score through `R = mean GHI / benchmark maximum`. Resource and climate context remain separate terms.

## Climate-context interpretation

The former normalized operational-stress quantity is relabeled as a **benchmark-relative climate-context index**. It is not a physical damage scale. Zero means the lowest composite value in the selected ten-point benchmark, not zero climate stress. The former 0.33/0.66 category thresholds have been removed.

## Coefficients

The climate weights and PV component scores are described as declared screening-scenario values. They are not described as field calibration. The PV baseline-suitability and resilience scores are now reconstructed explicitly from the component matrix and disclosed weights used in the submitted screening design.

## Robustness

The earlier ±10% aggregate Monte Carlo test is removed. The revised analysis uses deterministic leave-one-location-out normalization sensitivity, climate-component sensitivity, mean-temperature-only heat sensitivity, equal-weight technology sensitivity, and leave-one-technology-criterion-out tests. This avoids presenting an arbitrary perturbation width as a field probability model.

## Removed secondary analyses

The WSM/TOPSIS/VIKOR comparison and market-price uplift layer were removed rather than expanded because they were not central to the revised contribution and were not documented strongly enough in the initial submission.

## Reproducibility

The revision uses one committed input set and one end-to-end script for all revised numerical results. The exact NASA POWER point-query URLs are provided for independent upstream retrieval. The live archive can evolve, so the committed climate-summary table is the revision's numerical input snapshot.

## Claim boundary

The paper is an early-stage, historical climate-exposure screening study. It is not a climate projection, field-degradation validation, bankability model, or proof that one PV family is universally superior.
