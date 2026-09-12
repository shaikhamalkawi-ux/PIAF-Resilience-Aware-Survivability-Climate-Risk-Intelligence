#!/usr/bin/env python3
"""ECMX-D-26-01378 R3 deterministic analysis.

Main empirical calculation uses NASA POWER-derived annual/site metrics only.
The screen is bi-objective: maximize long-term mean GHI and minimize an
historical climate-context summary. Equal weighting of the four visible climate
components is a neutral summary convention, not a calibrated physical model.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'data'; OUT=ROOT/'outputs'; TABLES=OUT/'tables'; FIGS=OUT/'figures'
TABLES.mkdir(parents=True,exist_ok=True); FIGS.mkdir(parents=True,exist_ok=True)
ANNUAL=pd.read_csv(DATA/'annual_site_metrics.csv'); SEASON=pd.read_csv(DATA/'ghi_seasonality_diagnostics.csv')
LOC_ORDER=ANNUAL.location.drop_duplicates().tolist(); YEARS=sorted(ANNUAL.year.unique().tolist()); SEED=20260912; BOOTSTRAP_N=10000

def norm(x):
    x=np.asarray(x,float); lo,hi=np.nanmin(x),np.nanmax(x)
    return np.zeros_like(x) if np.isclose(hi,lo) else (x-lo)/(hi-lo)

def spear(a,b):
    ra=pd.Series(a).rank(method='average'); rb=pd.Series(b).rank(method='average'); return float(ra.corr(rb,method='pearson'))

def context(m,idx=None):
    idx=[1,2,3,4] if idx is None else idx; z=np.column_stack([norm(m[:,i]) for i in idx]); return z.mean(axis=1),z

def pareto(r,c):
    r=np.asarray(r,float); c=np.asarray(c,float)
    d=((r[:,None]>=r[None,:])&(c[:,None]<=c[None,:])&((r[:,None]>r[None,:])|(c[:,None]<c[None,:])))
    np.fill_diagonal(d,False); by=d.sum(axis=0); does=d.sum(axis=1); return by==0,does,by

cube=np.empty((len(LOC_ORDER),len(YEARS),5)); coords={}
for li,loc in enumerate(LOC_ORDER):
    s=ANNUAL[ANNUAL.location.eq(loc)].set_index('year').loc[YEARS]
    for j,col in enumerate(['mean_ghi','heat_p95_tmax','mean_rh','ghi_within_month_cv','wind_p95']): cube[li,:,j]=s[col].to_numpy(float)
    coords[loc]=(s.country.iloc[0],float(s.latitude.iloc[0]),float(s.longitude.iloc[0]))
M=cube.mean(axis=1); C,Z=context(M); P,does,by=pareto(M[:,0],C)
site=pd.DataFrame({'location':LOC_ORDER,'country':[coords[l][0] for l in LOC_ORDER],'latitude':[coords[l][1] for l in LOC_ORDER],'longitude':[coords[l][2] for l in LOC_ORDER],'mean_ghi_kwh_m2_day':M[:,0],'heat_mean_annual_p95_tmax_c':M[:,1],'mean_rh_pct':M[:,2],'ghi_mean_within_month_daily_cv':M[:,3],'wind_mean_annual_p95_ws10m_ms':M[:,4],'climate_context_equal_weight':C,'dominates_count':does,'dominated_by_count':by,'pareto_nondominated':P})
for j,n in enumerate(['heat','humidity','ghi_variability','wind_exposure']): site[f'normalized_{n}']=Z[:,j]

rng=np.random.default_rng(SEED); pc=np.zeros(len(LOC_ORDER),int); ds=np.zeros(len(LOC_ORDER)); cs=np.zeros(len(LOC_ORDER)); css=np.zeros(len(LOC_ORDER))
for _ in range(BOOTSTRAP_N):
    idx=rng.integers(0,len(YEARS),size=len(YEARS)); mb=cube[:,idx,:].mean(axis=1); cb,_=context(mb); pb,db,_=pareto(mb[:,0],cb)
    pc+=pb; ds+=db; cs+=cb; css+=cb*cb
site['bootstrap_pareto_frequency']=pc/BOOTSTRAP_N; site['bootstrap_mean_dominates_count']=ds/BOOTSTRAP_N; site['bootstrap_context_mean']=cs/BOOTSTRAP_N; site['bootstrap_context_sd']=np.sqrt(np.maximum(css/BOOTSTRAP_N-(cs/BOOTSTRAP_N)**2,0))
ly=np.zeros(len(LOC_ORDER),int); lctx=[]
for yi,_ in enumerate(YEARS):
    ml=np.delete(cube,yi,axis=1).mean(axis=1); cl,_=context(ml); pl,_,_=pareto(ml[:,0],cl); ly+=pl; lctx.append(cl)
lctx=np.stack(lctx); site['leave_one_year_out_pareto_frequency']=ly/len(YEARS); site['leave_one_year_out_context_range']=lctx.max(axis=0)-lctx.min(axis=0)
site.to_csv(TABLES/'table_1_site_resource_climate_screen.csv',index=False)

rows=[]
for oi,om in enumerate(LOC_ORDER):
    keep=[i for i in range(len(LOC_ORDER)) if i!=oi]; sub,_=context(M[keep,:]); base=C[keep]
    rows.append({'omitted_location':om,'spearman_context_rank':spear(base,sub),'max_absolute_rank_shift':int((pd.Series(base).rank(method='average')-pd.Series(sub).rank(method='average')).abs().max()),'max_absolute_context_change':float(np.max(np.abs(base-sub)))})
loo=pd.DataFrame(rows).sort_values(['spearman_context_rank','max_absolute_context_change']); loo.to_csv(TABLES/'table_2_leave_one_location_out_normalization.csv',index=False)

names=['heat','humidity','ghi_variability','wind_exposure']; scenarios=[('baseline',[1,2,3,4])]+[(f'omit_{n}',[x for x in [1,2,3,4] if x!=k+1]) for k,n in enumerate(names)]
rows=[]
for name,idxs in scenarios:
    cc,_=context(M,idxs); pp,_,_=pareto(M[:,0],cc)
    rows.append({'scenario':name,'spearman_vs_baseline_context':spear(C,cc),'max_absolute_rank_shift':int((pd.Series(C).rank()-pd.Series(cc).rank()).abs().max()),'pareto_set':'; '.join([LOC_ORDER[i] for i in range(len(LOC_ORDER)) if pp[i]]),'pareto_count':int(pp.sum())})
pd.DataFrame(rows).to_csv(TABLES/'table_3_component_omission_sensitivity.csv',index=False)

season=SEASON.set_index('location').loc[LOC_ORDER]
alts=[('within_month_daily_CV_baseline',season['mean_within_month_daily_cv'].to_numpy(float)),('raw_daily_CV',season['raw_daily_ghi_cv'].to_numpy(float)),('month_conditioned_interannual_CV',season['mean_month_conditioned_interannual_cv'].to_numpy(float)),('annual_mean_interannual_CV',season['interannual_annual_mean_ghi_cv'].to_numpy(float))]
rows=[]
for label,g in alts:
    mm=M.copy(); mm[:,3]=g; cg,_=context(mm); pg,_,_=pareto(mm[:,0],cg)
    rows.append({'ghi_variability_definition':label,'spearman_vs_baseline_context':spear(C,cg),'max_absolute_rank_shift':int((pd.Series(C).rank()-pd.Series(cg).rank()).abs().max()),'pareto_set':'; '.join([LOC_ORDER[i] for i in range(len(LOC_ORDER)) if pg[i]])})
gs=pd.DataFrame(rows); gs.to_csv(TABLES/'table_4_ghi_variability_definition_sensitivity.csv',index=False)

tech=pd.DataFrame([
{'PV family':'PERC','evidence_basis':'Qatar desert field comparison; damp-heat comparison','supported_interpretation':'Mature reference family; tested samples can be comparatively stable, but product and manufacturing differences remain important.','calculation_role':'Evidence context only; no family coefficient or rank.'},
{'PV family':'TOPCon','evidence_basis':'Qatar desert field comparison; 20-module-type accelerated-aging study; damp-heat handling study','supported_interpretation':'Strong efficiency/temperature-coefficient potential with substantial product-to-product moisture/UV reliability heterogeneity.','calculation_role':'Evidence context only; no family coefficient or rank.'},
{'PV family':'SHJ/HJT','evidence_basis':'Qatar desert field comparison; damp-heat handling study','supported_interpretation':'Favorable temperature coefficients can coexist with product-specific field or humidity-related degradation.','calculation_role':'Evidence context only; no family coefficient or rank.'},
{'PV family':'CdTe','evidence_basis':'DOE/IEA technology and hot-climate evidence','supported_interpretation':'Comparable PV technology with hot-climate relevance, but not covered by the same common four-family field experiment used for the crystalline families.','calculation_role':'Evidence context only; no family coefficient or rank.'}]); tech.to_csv(TABLES/'table_5_four_family_evidence_boundary.csv',index=False)

plt.rcParams.update({'font.size':9})
fig,ax=plt.subplots(figsize=(7.2,5.2)); ax.scatter(site.climate_context_equal_weight,site.mean_ghi_kwh_m2_day,s=55)
for _,r in site.iterrows(): ax.annotate(r.location,(r.climate_context_equal_weight,r.mean_ghi_kwh_m2_day),xytext=(4,4),textcoords='offset points',fontsize=8)
pdf=site[site.pareto_nondominated].sort_values('climate_context_equal_weight'); ax.plot(pdf.climate_context_equal_weight,pdf.mean_ghi_kwh_m2_day,marker='o'); ax.set_xlabel('Historical climate-context summary C (lower is less demanding within benchmark)'); ax.set_ylabel('Mean GHI (kWh m$^{-2}$ day$^{-1}$; higher is better)'); ax.set_title('Resource-climate pre-feasibility screen'); ax.grid(True,alpha=.25); fig.tight_layout(); fig.savefig(FIGS/'fig_1_resource_climate_pareto.png',dpi=300); plt.close(fig)
order=site.sort_values('climate_context_equal_weight').location.tolist(); q=site.set_index('location').loc[order,['normalized_heat','normalized_humidity','normalized_ghi_variability','normalized_wind_exposure']]; fig,ax=plt.subplots(figsize=(7.2,5.4)); im=ax.imshow(q.values,aspect='auto'); ax.set_yticks(range(len(order)),labels=order); ax.set_xticks(range(4),labels=['Heat','Humidity','GHI variability','Wind exposure'],rotation=25,ha='right')
for i in range(q.shape[0]):
    for j in range(q.shape[1]): ax.text(j,i,f'{q.values[i,j]:.2f}',ha='center',va='center',fontsize=7)
fig.colorbar(im,ax=ax,label='Benchmark-normalized component'); ax.set_title('Climate-context components remain visible'); fig.tight_layout(); fig.savefig(FIGS/'fig_2_climate_component_heatmap.png',dpi=300); plt.close(fig)
q=loo.sort_values('spearman_context_rank'); fig,ax=plt.subplots(figsize=(7.2,4.6)); ax.bar(q.omitted_location,q.spearman_context_rank); ax.set_ylim(0,1.05); ax.set_ylabel('Spearman rank correlation'); ax.set_xlabel('Location omitted from normalization benchmark'); ax.set_title('Leave-one-location-out normalization sensitivity'); ax.tick_params(axis='x',rotation=35); ax.grid(axis='y',alpha=.25); fig.tight_layout(); fig.savefig(FIGS/'fig_3_leave_one_location_normalization.png',dpi=300); plt.close(fig)
q=site.sort_values('bootstrap_pareto_frequency',ascending=False); fig,ax=plt.subplots(figsize=(7.2,4.8)); ax.bar(q.location,q.bootstrap_pareto_frequency); ax.set_ylim(0,1.05); ax.set_ylabel('Pareto inclusion frequency'); ax.set_xlabel('Location'); ax.set_title(f'Year-block bootstrap robustness ({BOOTSTRAP_N:,} resamples)'); ax.tick_params(axis='x',rotation=35); ax.grid(axis='y',alpha=.25); fig.tight_layout(); fig.savefig(FIGS/'fig_4_bootstrap_pareto_frequency.png',dpi=300); plt.close(fig)

summary={'site_count':len(LOC_ORDER),'years':YEARS,'bootstrap_iterations':BOOTSTRAP_N,'bootstrap_seed':SEED,'baseline_pareto_set':site.loc[site.pareto_nondominated,'location'].tolist(),'bootstrap_pareto_frequency':dict(zip(site.location,site.bootstrap_pareto_frequency.round(4))),'leave_one_year_out_pareto_frequency':dict(zip(site.location,site.leave_one_year_out_pareto_frequency.round(4))),'worst_leave_one_location_spearman':float(loo.spearman_context_rank.min()),'worst_leave_one_location_max_rank_shift':int(loo.max_absolute_rank_shift.max()),'all_ghi_definitions_same_pareto_set':bool(gs.pareto_set.nunique()==1),'scientific_boundary':['Main results use NASA POWER-derived climate/resource evidence only.','Equal component weighting is a neutral summary convention and is sensitivity-tested.','No numerical PV-family score or market coefficient enters the empirical result.','Historical pre-feasibility screening only; not climate projection, failure prediction, yield modeling, LCOE, or bankability analysis.']}
with open(OUT/'analysis_summary.json','w',encoding='utf-8') as f: json.dump(summary,f,indent=2)
rows=[]
for p in sorted(ROOT.rglob('*')):
    if p.is_file() and 'outputs/sha256_manifest.csv' not in str(p): rows.append({'relative_path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
pd.DataFrame(rows).to_csv(OUT/'sha256_manifest.csv',index=False); print(json.dumps(summary,indent=2))
