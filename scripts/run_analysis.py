#!/usr/bin/env python3
"""Reproduce the deterministic analyses used in the ECMX major revision.

Inputs are committed CSV files. The climate-summary CSV is the locked 2014-2023
summary used for the revision; exact NASA POWER query URLs are supplied separately
for independent retrieval of the daily source series.
"""
from __future__ import annotations
import hashlib
import json
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

SITES = pd.read_csv(DATA / "sites_climate_summary.csv")
TECH = pd.read_csv(DATA / "pv_scenario_intervals.csv")


def normalize(s: pd.Series) -> pd.Series:
    lo, hi = float(s.min()), float(s.max())
    if np.isclose(hi, lo):
        return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    return (s - lo) / (hi - lo)


def climate_context(df: pd.DataFrame, weights=(0.40, 0.20, 0.20, 0.20), heat=(0.60, 0.40)):
    """Benchmark-relative climate-context index.

    Components: heat context, raw-daily GHI CV, raw-daily wind CV, mean RH.
    The coefficients are declared screening-scenario weights, not calibrated
    physical damage coefficients.
    """
    wt = np.asarray(weights, dtype=float)
    wt = wt / wt.sum()
    ht = np.asarray(heat, dtype=float)
    ht = ht / ht.sum()
    z_tmean = normalize(df["mean_t2m_c"])
    z_tmax = normalize(df["mean_t2m_max_c"])
    z_heat = ht[0] * z_tmean + ht[1] * z_tmax
    z_ghi = normalize(df["cv_ghi"])
    z_wind = normalize(df["cv_ws10m"])
    z_rh = normalize(df["mean_rh_pct"])
    raw = wt[0] * z_heat + wt[1] * z_ghi + wt[2] * z_wind + wt[3] * z_rh
    context = normalize(raw)
    return pd.DataFrame({
        "z_tmean": z_tmean,
        "z_tmax": z_tmax,
        "z_heat": z_heat,
        "z_ghi_cv": z_ghi,
        "z_wind_cv": z_wind,
        "z_rh": z_rh,
        "context_raw": raw,
        "climate_context": context,
    })


def score(resource_ratio, context, s0, rho, lam=1.0):
    return resource_ratio * s0 * np.exp(-lam * context * (1.0 - rho))


def spearman_rank(a, b):
    ra = pd.Series(a).rank(method="average")
    rb = pd.Series(b).rank(method="average")
    return float(ra.corr(rb, method="pearson"))

ctx = climate_context(SITES)
base = pd.concat([SITES.copy(), ctx], axis=1)
base["resource_ratio"] = base["mean_ghi_kwh_m2_day"] / base["mean_ghi_kwh_m2_day"].max()
base["context_rank"] = base["climate_context"].rank(method="min", ascending=False).astype(int)
base["resource_rank"] = base["resource_ratio"].rank(method="min", ascending=False).astype(int)
base.to_csv(TABLES / "table_1_site_climate_resource_context.csv", index=False)

rows = []
for _, srow in base.iterrows():
    for _, t in TECH.iterrows():
        lo = score(srow.resource_ratio, srow.climate_context, t.s0_lower, t.rho_lower)
        ce = score(srow.resource_ratio, srow.climate_context, t.s0_central, t.rho_central)
        hi = score(srow.resource_ratio, srow.climate_context, t.s0_upper, t.rho_upper)
        rows.append({
            "location": srow.location,
            "country": srow.country,
            "technology": t.technology,
            "resource_ratio": srow.resource_ratio,
            "climate_context": srow.climate_context,
            "score_lower": lo,
            "score_central": ce,
            "score_upper": hi,
            "central_rank_within_location": np.nan,
        })
scores = pd.DataFrame(rows)
scores["central_rank_within_location"] = scores.groupby("location")["score_central"].rank(method="min", ascending=False).astype(int)
scores.to_csv(TABLES / "table_2_resource_resilience_scores.csv", index=False)
mat = scores.pivot(index="location", columns="technology", values="score_central").loc[base.location]
mat.to_csv(TABLES / "table_3_central_score_matrix.csv")

cert = []
for loc, g in scores.groupby("location", sort=False):
    leader = g.loc[g.score_central.idxmax()]
    comp_max = g.loc[g.technology != leader.technology, "score_upper"].max()
    cert.append({
        "location": loc,
        "central_leader": leader.technology,
        "leader_central": leader.score_central,
        "leader_lower": leader.score_lower,
        "largest_competitor_upper": comp_max,
        "guaranteed_dominance_over_declared_envelopes": bool(leader.score_lower > comp_max),
    })
cert_df = pd.DataFrame(cert)
cert_df.to_csv(TABLES / "table_4_interval_dominance_audit.csv", index=False)

def critical_rho(high_name, low_name, lam=1.0):
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

components = {
    "baseline_40_20_20_20": ((0.40,0.20,0.20,0.20),(0.60,0.40)),
    "equal_components": ((0.25,0.25,0.25,0.25),(0.60,0.40)),
    "heat_only": ((1.0,0.0,0.0,0.0),(0.60,0.40)),
    "omit_solar_variability": ((0.50,0.0,0.25,0.25),(0.60,0.40)),
    "omit_wind_variability": ((0.50,0.25,0.0,0.25),(0.60,0.40)),
    "omit_humidity": ((0.50,0.25,0.25,0.0),(0.60,0.40)),
    "omit_heat": ((0.0,1/3,1/3,1/3),(0.60,0.40)),
    "mean_temperature_only_inside_heat": ((0.40,0.20,0.20,0.20),(1.0,0.0)),
}
sens_rows = []
baseC = base.climate_context
for name,(w,h) in components.items():
    c = climate_context(SITES, weights=w, heat=h).climate_context
    sens_rows.append({
        "scenario": name,
        "spearman_vs_baseline": spearman_rank(baseC, c),
        "max_absolute_context_change": float(np.max(np.abs(baseC-c))),
        "max_absolute_rank_shift": int(np.max(np.abs(baseC.rank(method='min',ascending=False) - c.rank(method='min',ascending=False)))),
        "top_context_location": SITES.loc[c.idxmax(),"location"],
    })
pd.DataFrame(sens_rows).to_csv(TABLES / "table_7_climate_context_structural_sensitivity.csv", index=False)

pd.DataFrame({
    "metric": ["Pearson correlation: mean T2M vs mean T2M_MAX", "Interpretation"],
    "value": [f"{SITES.mean_t2m_c.corr(SITES.mean_t2m_max_c):.6f}", "The two temperature summaries are combined inside one heat component; they are not separate top-level criteria."],
}).to_csv(TABLES / "table_8_temperature_redundancy_audit.csv", index=False)

leader_rows=[]
for name,(w,h) in components.items():
    c = climate_context(SITES, weights=w, heat=h).climate_context
    for j, srow in SITES.iterrows():
        rr = srow.mean_ghi_kwh_m2_day / SITES.mean_ghi_kwh_m2_day.max()
        vals = {t.technology: score(rr,c.iloc[j],t.s0_central,t.rho_central) for _,t in TECH.iterrows()}
        leader=max(vals,key=vals.get)
        leader_rows.append({"scenario":name,"location":srow.location,"central_leader":leader,"leader_score":vals[leader]})
pd.DataFrame(leader_rows).to_csv(TABLES / "table_9_central_leader_structural_sensitivity.csv",index=False)

fig, ax = plt.subplots(figsize=(8.1,5.4))
ax.scatter(base.resource_ratio, base.climate_context, s=52)
label_offsets = {"Alice Springs": (5, 8), "Antofagasta": (5, -10), "Aswan": (5, 5), "NEOM Bay/Sharma": (5, 5)}
for _, r in base.iterrows():
    ax.annotate(r.location, (r.resource_ratio, r.climate_context), xytext=label_offsets.get(r.location,(4,4)), textcoords="offset points", fontsize=8)
ax.set_xlabel("Relative mean solar resource, R = mean GHI / benchmark maximum")
ax.set_ylabel("Benchmark-relative climate-context index, C")
ax.set_xlim(0.72,1.02); ax.set_ylim(-0.04,1.04); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "fig_1_resource_vs_climate_context.png", dpi=300); plt.close(fig)

fig, ax = plt.subplots(figsize=(9.2,4.8))
img = ax.imshow(mat.T.values, aspect="auto")
ax.set_xticks(range(len(mat.index)), labels=mat.index, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(mat.columns)), labels=mat.columns, fontsize=9)
ax.set_xlabel("Benchmark location"); ax.set_ylabel("PV family")
for i in range(mat.T.shape[0]):
    for j in range(mat.T.shape[1]):
        ax.text(j,i,f"{mat.T.iloc[i,j]:.3f}",ha="center",va="center",fontsize=7)
fig.colorbar(img, ax=ax, label="Central resource-resilience screening score")
fig.tight_layout(); fig.savefig(FIGS / "fig_2_central_score_heatmap.png", dpi=300); plt.close(fig)

abu=scores[scores.location=="Abu Dhabi"].set_index("technology").loc[TECH.technology]
fig, ax = plt.subplots(figsize=(7.6,4.6)); y=np.arange(len(abu))
ax.errorbar(abu.score_central, y, xerr=np.vstack([abu.score_central-abu.score_lower,abu.score_upper-abu.score_central]), fmt='o', capsize=4)
ax.set_yticks(y, labels=abu.index); ax.set_xlabel("Resource-resilience screening score")
ax.set_title("Abu Dhabi: central values and declared scenario envelopes"); ax.grid(axis='x',alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "fig_3_abu_dhabi_score_envelopes.png", dpi=300); plt.close(fig)

loo_plot=loo.sort_values("spearman_rank_correlation")
fig, ax = plt.subplots(figsize=(8.0,5.0)); y=np.arange(len(loo_plot))
ax.barh(y, loo_plot.spearman_rank_correlation); ax.set_yticks(y, labels=loo_plot.omitted_location, fontsize=8)
ax.set_xlabel("Spearman rank correlation with full 10-location climate-context ranking"); ax.set_xlim(0,1.02); ax.grid(axis='x',alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "fig_4_leave_one_out_normalization_sensitivity.png", dpi=300); plt.close(fig)

mean_only=climate_context(SITES, heat=(1,0)).climate_context
fig, ax = plt.subplots(figsize=(8.4,5.0)); x=np.arange(len(SITES)); width=.38
ax.bar(x-width/2, baseC, width, label="Baseline heat: 0.60 mean + 0.40 max")
ax.bar(x+width/2, mean_only, width, label="Mean-temperature-only heat sensitivity")
ax.set_xticks(x, labels=SITES.location, rotation=45, ha='right', fontsize=8)
ax.set_ylabel("Benchmark-relative climate-context index"); ax.legend(fontsize=8); ax.grid(axis='y',alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "fig_5_temperature_component_sensitivity.png", dpi=300); plt.close(fig)

summary = {
    "site_count": int(len(SITES)),
    "technology_count": int(len(TECH)),
    "baseline_context_min_location": base.loc[base.climate_context.idxmin(),"location"],
    "baseline_context_max_location": base.loc[base.climate_context.idxmax(),"location"],
    "temperature_pearson": float(SITES.mean_t2m_c.corr(SITES.mean_t2m_max_c)),
    "central_leader_all_locations": bool(cert_df.central_leader.eq("TOPCon PV").all()),
    "guaranteed_dominance_count": int(cert_df.guaranteed_dominance_over_declared_envelopes.sum()),
    "worst_leave_one_out_spearman": float(loo.spearman_rank_correlation.min()),
    "worst_leave_one_out_omission": str(loo.iloc[0].omitted_location),
    "notes": [
        "Climate-context weights are declared screening-scenario weights, not calibrated physical damage coefficients.",
        "The GHI and wind coefficients of variation use the raw daily series and therefore combine seasonal and day-to-day variability.",
        "C=0 denotes the lowest value in this benchmark, not absence of physical climate stress.",
        "Technology lower/upper values are declared screening envelopes, not empirical confidence intervals.",
    ],
}
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
manifest=[]
for p in sorted(list(DATA.glob("*.csv"))+list(TABLES.glob("*.csv"))+list(FIGS.glob("*.png"))+[OUT/"analysis_summary.json"]):
    b=p.read_bytes(); manifest.append({"path":str(p.relative_to(ROOT)).replace('\\','/'),"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()})
pd.DataFrame(manifest).to_csv(OUT / "sha256_manifest.csv", index=False)
print(json.dumps(summary, indent=2))
print(f"Wrote {len(list(TABLES.glob('*.csv')))} tables and {len(list(FIGS.glob('*.png')))} figures.")
