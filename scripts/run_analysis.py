#!/usr/bin/env python3
"""Deterministic major-revision analysis for ECMX-D-26-01378.

The climate layer uses committed 2014-2023 NASA POWER summary statistics.
Technology descriptors and weights are declared screening scenarios, not field calibration.
All numerical tables and publication figures are regenerated from the same inputs.
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "outputs"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)
for p in TABLES.glob("*.csv"): p.unlink()
for p in FIGS.glob("*.png"): p.unlink()

SITES = pd.read_csv(DATA / "sites_climate_summary.csv")
TECH = pd.read_csv(DATA / "pv_screening_components.csv")
SETTINGS = pd.read_csv(DATA / "model_settings.csv")

SUIT_COLS = ["E_performance", "M_market_maturity", "C_cost_attractiveness", "OM_om_simplicity"]
RES_COLS = ["RT_thermal_tolerance", "RD_degradation_resistance", "RH_humidity_soiling_tolerance", "RO_operational_robustness"]
SUIT_SYMBOLS = ["E", "M", "C", "OM"]
RES_SYMBOLS = ["RT", "RD", "RH", "RO"]
CLIMATE_SYMBOLS = ["H", "CVGHI", "CVwind", "RH2M"]
HEAT_SYMBOLS = ["Tmean", "Tmax"]

# Fail early if repository files drift out of sync.
for col in SUIT_COLS + RES_COLS:
    if col not in TECH.columns:
        raise ValueError(f"Missing technology component column: {col}")

def setting_weights(layer: str, symbols: list[str]) -> np.ndarray:
    d = SETTINGS[SETTINGS.layer.eq(layer)].set_index("symbol")["value"]
    missing = [s for s in symbols if s not in d.index]
    if missing:
        raise ValueError(f"Missing {layer} settings for: {missing}")
    w = np.array([float(d[s]) for s in symbols], dtype=float)
    if np.any(w < 0) or np.isclose(w.sum(), 0):
        raise ValueError(f"Invalid weights for {layer}: {w}")
    return w / w.sum()

WS = setting_weights("baseline_suitability", SUIT_SYMBOLS)
WR = setting_weights("resilience", RES_SYMBOLS)
WC = setting_weights("climate_context", CLIMATE_SYMBOLS)
WH = setting_weights("heat_component", HEAT_SYMBOLS)
LAMBDA = float(SETTINGS[(SETTINGS.layer == "attenuation") & (SETTINGS.symbol == "lambda")].value.iloc[0])


def normalize(s: pd.Series) -> pd.Series:
    lo, hi = float(s.min()), float(s.max())
    if np.isclose(hi, lo):
        return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    return (s - lo) / (hi - lo)


def climate_context(df: pd.DataFrame, weights=None, heat=None) -> pd.DataFrame:
    wt = WC if weights is None else np.asarray(weights, dtype=float)
    ht = WH if heat is None else np.asarray(heat, dtype=float)
    wt = wt / wt.sum(); ht = ht / ht.sum()
    z_tmean = normalize(df["mean_t2m_c"])
    z_tmax = normalize(df["mean_t2m_max_c"])
    z_heat = ht[0] * z_tmean + ht[1] * z_tmax
    z_ghi = normalize(df["cv_ghi"])
    z_wind = normalize(df["cv_ws10m"])
    z_rh = normalize(df["mean_rh_pct"])
    raw = wt[0] * z_heat + wt[1] * z_ghi + wt[2] * z_wind + wt[3] * z_rh
    context = normalize(raw)
    return pd.DataFrame({
        "z_tmean": z_tmean, "z_tmax": z_tmax, "z_heat": z_heat,
        "z_ghi_cv": z_ghi, "z_wind_cv": z_wind, "z_rh": z_rh,
        "context_raw": raw, "climate_context": context
    })


def technology_scores(ws=WS, wr=WR) -> pd.DataFrame:
    t = TECH.copy()
    t["S0"] = t[SUIT_COLS].to_numpy(float) @ np.asarray(ws, dtype=float)
    t["rho"] = t[RES_COLS].to_numpy(float) @ np.asarray(wr, dtype=float)
    return t


def screen_score(resource_ratio, context, s0, rho, lam=LAMBDA):
    return resource_ratio * s0 * np.exp(-lam * context * (1.0 - rho))


def spearman_rank(a, b):
    ra = pd.Series(a).rank(method="average")
    rb = pd.Series(b).rank(method="average")
    return float(ra.corr(rb, method="pearson"))


def normalized_without(w, idx):
    x = np.array(w, dtype=float); x[idx] = 0.0
    return x / x.sum()

# Central benchmark
ctx = climate_context(SITES)
base = pd.concat([SITES.copy(), ctx], axis=1)
base["resource_ratio"] = base["mean_ghi_kwh_m2_day"] / base["mean_ghi_kwh_m2_day"].max()
base["context_rank"] = base["climate_context"].rank(method="min", ascending=False).astype(int)
base["resource_rank"] = base["resource_ratio"].rank(method="min", ascending=False).astype(int)
base.to_csv(TABLES / "table_1_site_climate_resource_context.csv", index=False)

tech0 = technology_scores()
tech0[["technology"] + SUIT_COLS + RES_COLS + ["S0", "rho", "status"]].to_csv(
    TABLES / "table_2_pv_component_matrix_and_scores.csv", index=False)

rows = []
for _, srow in base.iterrows():
    for _, t in tech0.iterrows():
        q = t.S0 * np.exp(-LAMBDA * srow.climate_context * (1.0 - t.rho))
        rows.append({
            "location": srow.location, "country": srow.country, "technology": t.technology,
            "resource_ratio": srow.resource_ratio, "climate_context": srow.climate_context,
            "technical_score_Q": q,
            "resource_resilience_score_A": srow.resource_ratio * q,
        })
scores = pd.DataFrame(rows)
scores["rank_within_location"] = scores.groupby("location")["resource_resilience_score_A"].rank(method="min", ascending=False).astype(int)
scores.to_csv(TABLES / "table_3_resource_resilience_scores.csv", index=False)
mat = scores.pivot(index="location", columns="technology", values="resource_resilience_score_A").loc[base.location]
mat.to_csv(TABLES / "table_4_central_score_matrix.csv")

# Resource-vs-context inversion examples
# Critical resilience value at which the higher-resource location ceases to offset its larger context penalty.
def critical_rho(high_name, low_name, lam=LAMBDA):
    a = base.set_index("location").loc[high_name]
    b = base.set_index("location").loc[low_name]
    gap = float(a.climate_context - b.climate_context)
    rratio = float(a.resource_ratio / b.resource_ratio)
    crit = np.nan if gap <= 0 or rratio <= 1 else 1.0 - np.log(rratio) / (lam * gap)
    return {
        "higher_resource_location": high_name,
        "comparison_location": low_name,
        "resource_advantage_pct": 100.0 * (rratio - 1.0),
        "climate_context_gap": gap,
        "critical_resilience_rho": crit,
    }
pd.DataFrame([
    critical_rho("NEOM Bay/Sharma", "Phoenix"),
    critical_rho("NEOM Bay/Sharma", "Ouarzazate"),
    critical_rho("Jodhpur", "Dunhuang"),
]).to_csv(TABLES / "table_5_resource_context_inversion.csv", index=False)

# Benchmark-dependence audit: leave one location out and recompute normalization.
base_order = base.set_index("location")["climate_context"]
loo_rows = []
for omitted in SITES.location:
    sub = SITES[SITES.location != omitted].reset_index(drop=True)
    values = pd.Series(climate_context(sub).climate_context.values, index=sub.location)
    common = [x for x in SITES.location if x != omitted]
    base_rank = base_order.loc[common].rank(method="min", ascending=False)
    new_rank = values.loc[common].rank(method="min", ascending=False)
    loo_rows.append({
        "omitted_location": omitted,
        "n_remaining": len(common),
        "spearman_rank_correlation": spearman_rank(base_order.loc[common], values.loc[common]),
        "max_absolute_rank_shift": int((base_rank - new_rank).abs().max()),
        "max_absolute_context_change": float((base_order.loc[common] - values.loc[common]).abs().max()),
    })
loo = pd.DataFrame(loo_rows).sort_values(["spearman_rank_correlation", "max_absolute_context_change"])
loo.to_csv(TABLES / "table_6_leave_one_location_out.csv", index=False)

# Structural climate sensitivity.
climate_scenarios = {
    "baseline": (WC, WH),
    "equal_components": (np.ones(4)/4, WH),
    "omit_solar_variability": (normalized_without(WC, 1), WH),
    "omit_wind_variability": (normalized_without(WC, 2), WH),
    "omit_humidity": (normalized_without(WC, 3), WH),
    "omit_heat": (normalized_without(WC, 0), WH),
    "mean_temperature_only_inside_heat": (WC, np.array([1.0, 0.0])),
}
cs_rows = []
for name, (w, h) in climate_scenarios.items():
    c = climate_context(SITES, weights=w, heat=h).climate_context
    cs_rows.append({
        "scenario": name,
        "spearman_vs_baseline": spearman_rank(base.climate_context, c),
        "max_absolute_context_change": float(np.max(np.abs(base.climate_context-c))),
        "max_absolute_rank_shift": int(np.max(np.abs(base.climate_context.rank(method='min',ascending=False) - c.rank(method='min',ascending=False)))),
        "top_context_location": SITES.loc[c.idxmax(), "location"],
    })
pd.DataFrame(cs_rows).to_csv(TABLES / "table_7_climate_context_structural_sensitivity.csv", index=False)

# Technology weighting sensitivity: equal weights and leave-one-criterion-out.
tech_scenarios = [("baseline", WS, WR), ("equal_suitability_and_resilience", np.ones(len(WS))/len(WS), np.ones(len(WR))/len(WR))]
for i, s in enumerate(SUIT_SYMBOLS): tech_scenarios.append((f"omit_suitability_{s}", normalized_without(WS,i), WR))
for i, s in enumerate(RES_SYMBOLS): tech_scenarios.append((f"omit_resilience_{s}", WS, normalized_without(WR,i)))
tw_rows = []
for name, ws, wr in tech_scenarios:
    ts = technology_scores(ws, wr); leaders = []
    for _, srow in base.iterrows():
        vals = {t.technology: screen_score(srow.resource_ratio, srow.climate_context, t.S0, t.rho) for _, t in ts.iterrows()}
        leaders.append(max(vals, key=vals.get))
    counts = pd.Series(leaders).value_counts()
    tw_rows.append({
        "scenario": name,
        "topcon_leader_sites": int(counts.get("TOPCon PV",0)),
        "perc_leader_sites": int(counts.get("PERC PV",0)),
        "shj_hjt_leader_sites": int(counts.get("SHJ/HJT PV",0)),
        "cdte_leader_sites": int(counts.get("CdTe PV",0)),
    })
pd.DataFrame(tw_rows).to_csv(TABLES / "table_8_technology_weight_structural_sensitivity.csv", index=False)

# Deterministic descriptor-uncertainty audit.
# For each site, identify the central winner and runner-up and compute the
# smallest feasible *paired* perturbation to one corresponding 0-1 descriptor
# that makes them tie: the winner descriptor is decreased by delta while the
# runner-up descriptor is increased by the same delta. This is a local
# sensitivity margin, not a probability statement and not an uncertainty
# distribution. Smaller margins indicate greater dependence on the declared
# technology scenario values.
def paired_descriptor_reversal_margin(srow, tech_scores):
    vals = {
        t.technology: screen_score(srow.resource_ratio, srow.climate_context, t.S0, t.rho)
        for _, t in tech_scores.iterrows()
    }
    order = sorted(vals, key=vals.get, reverse=True)
    winner, runner = order[0], order[1]
    tw = tech_scores.set_index("technology").loc[winner]
    tr = tech_scores.set_index("technology").loc[runner]
    c = float(srow.climate_context)
    log_gap = float(np.log(vals[winner] / vals[runner]))
    candidates = []

    # Suitability descriptors: solve the tie condition analytically when the
    # same descriptor is moved by -delta for the winner and +delta for runner.
    k = float(np.exp(LAMBDA * c * (tr.rho - tw.rho)))
    for col, w, sym in zip(SUIT_COLS, WS, SUIT_SYMBOLS):
        delta = float((tw.S0 - k * tr.S0) / (w * (1.0 + k)))
        feasible = delta >= -1e-12 and delta <= float(tw[col]) + 1e-12 and delta <= (1.0 - float(tr[col])) + 1e-12
        if feasible:
            candidates.append((max(0.0, delta), f"{sym} suitability descriptor"))

    # Resilience descriptors: the paired move changes rho_w-rho_r by 2*w*delta.
    if c > 1e-12:
        for col, w, sym in zip(RES_COLS, WR, RES_SYMBOLS):
            delta = float(log_gap / (2.0 * LAMBDA * c * w))
            feasible = delta >= -1e-12 and delta <= float(tw[col]) + 1e-12 and delta <= (1.0 - float(tr[col])) + 1e-12
            if feasible:
                candidates.append((max(0.0, delta), f"{sym} resilience descriptor"))

    if not candidates:
        return {
            "location": srow.location, "winner": winner, "runner_up": runner,
            "score_ratio_winner_to_runner": float(vals[winner] / vals[runner]),
            "minimum_paired_descriptor_delta": np.nan,
            "most_sensitive_descriptor": "No feasible single-descriptor reversal within [0,1]",
        }
    delta, desc = min(candidates, key=lambda x: x[0])
    return {
        "location": srow.location, "winner": winner, "runner_up": runner,
        "score_ratio_winner_to_runner": float(vals[winner] / vals[runner]),
        "minimum_paired_descriptor_delta": float(delta),
        "most_sensitive_descriptor": desc,
    }

margin_rows = [paired_descriptor_reversal_margin(r, tech0) for _, r in base.iterrows()]
margins = pd.DataFrame(margin_rows).sort_values("minimum_paired_descriptor_delta")
margins.to_csv(TABLES / "table_9_descriptor_reversal_margins.csv", index=False)

# Thermal redundancy audit.
temp_corr = float(SITES.mean_t2m_c.corr(SITES.mean_t2m_max_c))
pd.DataFrame({
    "metric": ["Pearson correlation: mean T2M vs mean T2M_MAX", "Interpretation"],
    "value": [f"{temp_corr:.6f}", "The two summaries are combined inside one heat component; mean-temperature-only sensitivity is reported."],
}).to_csv(TABLES / "table_10_temperature_redundancy_audit.csv", index=False)

# Publication figures (all generated from the same inputs/tables).
fig, ax = plt.subplots(figsize=(8.1,5.4))
ax.scatter(base.resource_ratio, base.climate_context, s=52)
off={"Alice Springs":(5,8),"Antofagasta":(5,-10),"Aswan":(5,5),"NEOM Bay/Sharma":(5,5)}
for _,r in base.iterrows(): ax.annotate(r.location,(r.resource_ratio,r.climate_context),xytext=off.get(r.location,(4,4)),textcoords="offset points",fontsize=8)
ax.set_xlabel("Relative mean solar resource, R = mean GHI / benchmark maximum")
ax.set_ylabel("Benchmark-relative climate-context index, C")
ax.set_xlim(0.72,1.02); ax.set_ylim(-0.04,1.04); ax.grid(alpha=0.25); fig.tight_layout()
fig.savefig(FIGS / "fig_1_resource_vs_climate_context.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(9.2,4.8)); img=ax.imshow(mat.T.values,aspect="auto")
ax.set_xticks(range(len(mat.index)),labels=mat.index,rotation=45,ha="right",fontsize=8)
ax.set_yticks(range(len(mat.columns)),labels=mat.columns,fontsize=9)
ax.set_xlabel("Benchmark location"); ax.set_ylabel("PV family")
for i in range(mat.T.shape[0]):
    for j in range(mat.T.shape[1]): ax.text(j,i,f"{mat.T.iloc[i,j]:.3f}",ha="center",va="center",fontsize=7)
fig.colorbar(img,ax=ax,label="Resource-resilience screening score, A"); fig.tight_layout()
fig.savefig(FIGS / "fig_2_central_score_heatmap.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(8.1,4.8)); abu=scores[scores.location=="Abu Dhabi"].sort_values("resource_resilience_score_A",ascending=True)
ax.barh(abu.technology,abu.resource_resilience_score_A); ax.set_xlabel("Resource-resilience screening score, A")
ax.set_title("Abu Dhabi: conditional four-family screening scenario"); ax.grid(axis='x',alpha=0.25); fig.tight_layout()
fig.savefig(FIGS / "fig_3_abu_dhabi_central_scores.png", dpi=300); plt.close(fig)

loo_plot=loo.sort_values("spearman_rank_correlation"); fig, ax = plt.subplots(figsize=(8.0,5.0)); y=np.arange(len(loo_plot))
ax.barh(y,loo_plot.spearman_rank_correlation); ax.set_yticks(y,labels=loo_plot.omitted_location,fontsize=8)
ax.set_xlabel("Spearman rank correlation with full 10-location climate-context ranking"); ax.set_xlim(0,1.02); ax.grid(axis='x',alpha=0.25); fig.tight_layout()
fig.savefig(FIGS / "fig_4_leave_one_out_normalization_sensitivity.png", dpi=300); plt.close(fig)

mean_only=climate_context(SITES, heat=np.array([1.0,0.0])).climate_context; fig, ax = plt.subplots(figsize=(8.4,5.0)); x=np.arange(len(SITES)); width=.38
ax.bar(x-width/2,base.climate_context,width,label="Baseline heat: 0.60 mean + 0.40 max")
ax.bar(x+width/2,mean_only,width,label="Mean-temperature-only heat sensitivity")
ax.set_xticks(x,labels=SITES.location,rotation=45,ha='right',fontsize=8); ax.set_ylabel("Benchmark-relative climate-context index, C")
ax.legend(fontsize=8); ax.grid(axis='y',alpha=0.25); fig.tight_layout()
fig.savefig(FIGS / "fig_5_temperature_component_sensitivity.png", dpi=300); plt.close(fig)

leaders = scores.loc[scores.groupby("location")["resource_resilience_score_A"].idxmax()]
summary = {
    "site_count": int(len(SITES)),
    "technology_count": int(len(TECH)),
    "baseline_context_min_location": base.loc[base.climate_context.idxmin(),"location"],
    "baseline_context_max_location": base.loc[base.climate_context.idxmax(),"location"],
    "temperature_pearson": temp_corr,
    "baseline_leader_counts": leaders.technology.value_counts().to_dict(),
    "worst_leave_one_out_spearman": float(loo.spearman_rank_correlation.min()),
    "worst_leave_one_out_omission": str(loo.iloc[0].omitted_location),
    "minimum_paired_descriptor_reversal_margin": float(margins.minimum_paired_descriptor_delta.min()),
    "minimum_margin_location": str(margins.iloc[0].location),
    "notes": [
        "Climate-context and technology weights are declared screening-scenario weights, not calibrated physical coefficients.",
        "GHI and wind coefficients of variation use the raw daily series and therefore combine seasonal and day-to-day variability.",
        "C=0 denotes the lowest value in this benchmark, not absence of physical climate stress.",
        "PV-family scores are conditional scenario outputs; no universal family superiority is claimed.",
    ],
}
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

manifest=[]
for p in sorted(list(DATA.glob("*.csv")) + list(TABLES.glob("*.csv")) + list(FIGS.glob("*.png")) + [OUT/"analysis_summary.json"]):
    b=p.read_bytes(); manifest.append({"path":str(p.relative_to(ROOT)).replace('\\','/'),"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()})
pd.DataFrame(manifest).to_csv(OUT / "sha256_manifest.csv", index=False)
print(json.dumps(summary, indent=2))
print(f"Wrote {len(list(TABLES.glob('*.csv')))} tables and {len(list(FIGS.glob('*.png')))} figures.")
