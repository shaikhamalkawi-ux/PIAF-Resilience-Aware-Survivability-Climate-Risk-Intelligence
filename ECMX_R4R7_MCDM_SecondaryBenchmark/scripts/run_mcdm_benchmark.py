#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'R4R6_table_1_site_resource_climate_screen.csv'
OUT = ROOT / 'outputs'
OUT.mkdir(exist_ok=True)
BASE_PARETO = {'Antofagasta','Aswan'}
V = 0.5


def ranks(values, higher_better=True):
    return pd.Series(values).rank(ascending=not higher_better, method='min').astype(int).to_numpy()


def mcdm_common(X, weights, v=0.5):
    """WSM, TOPSIS and VIKOR on one shared [0,1] benefit-oriented matrix."""
    X = np.asarray(X, float)
    w = np.asarray(weights, float)
    w = w / w.sum()
    wsm = X @ w

    # TOPSIS: use the same already-normalized matrix; weights are applied once.
    Y = X * w
    ideal_pos = Y.max(axis=0)
    ideal_neg = Y.min(axis=0)
    d_pos = np.linalg.norm(Y - ideal_pos, axis=1)
    d_neg = np.linalg.norm(Y - ideal_neg, axis=1)
    den = d_pos + d_neg
    topsis = np.divide(d_neg, den, out=np.full_like(d_neg, 0.5), where=~np.isclose(den, 0))

    # VIKOR: all X columns are benefits on [0,1].
    best = X.max(axis=0)
    worst = X.min(axis=0)
    span = best - worst
    loss = np.divide(best - X, span, out=np.zeros_like(X), where=~np.isclose(span, 0))
    weighted_loss = loss * w
    S = weighted_loss.sum(axis=1)
    R = weighted_loss.max(axis=1)
    Sstar, Sminus = S.min(), S.max()
    Rstar, Rminus = R.min(), R.max()
    Sterm = np.zeros_like(S) if np.isclose(Sminus, Sstar) else (S-Sstar)/(Sminus-Sstar)
    Rterm = np.zeros_like(R) if np.isclose(Rminus, Rstar) else (R-Rstar)/(Rminus-Rstar)
    Q = v*Sterm + (1-v)*Rterm
    return {'WSM':wsm, 'TOPSIS':topsis, 'VIKOR_S':S, 'VIKOR_R':R, 'VIKOR_Q':Q}


def top_order(values, higher_better=True):
    return np.argsort(-np.asarray(values) if higher_better else np.asarray(values), kind='stable')


def compress_intervals(grid, labels):
    rows=[]; start=0; last=labels[0]
    for i in range(1,len(labels)):
        if labels[i] != last:
            rows.append((float(grid[start]), float(grid[i-1]), last))
            start=i; last=labels[i]
    rows.append((float(grid[start]), float(grid[-1]), last))
    return rows


def vikor_compromise_set(q_values, s_values, r_values, names, dq, tol=1e-12):
    """Return standard VIKOR C1/C2 diagnostics and compromise set."""
    q_values=np.asarray(q_values,float); s_values=np.asarray(s_values,float); r_values=np.asarray(r_values,float)
    order=top_order(q_values,False)
    q1=q_values[order[0]]; q2=q_values[order[1]]
    c1=bool((q2-q1) >= dq - tol)
    qwinner=names[order[0]]
    swinner=names[top_order(s_values,False)[0]]
    rwinner=names[top_order(r_values,False)[0]]
    c2=bool(qwinner in {swinner,rwinner})
    if c1 and c2:
        comp=[qwinner]
    elif c1 and not c2:
        comp=[names[order[0]], names[order[1]]]
    else:
        comp=[names[i] for i in order if (q_values[i]-q1) < dq + tol]
    return {
        'order':order,'q1':q1,'q2':q2,'C1':c1,'C2':c2,
        'Q_best_site':qwinner,'S_best_site':swinner,'R_best_site':rwinner,
        'compromise_set':comp,
    }

# Synthetic branch checks for the standard VIKOR compromise rule.
_test_names=np.array(['A1','A2','A3'])
_t=vikor_compromise_set([0.0,0.6,1.0],[0.0,0.7,1.0],[0.0,0.8,1.0],_test_names,0.5)
assert _t['C1'] and _t['C2'] and _t['compromise_set']==['A1']
_t=vikor_compromise_set([0.0,0.6,1.0],[0.2,0.0,1.0],[0.3,0.1,1.0],_test_names,0.5)
assert _t['C1'] and (not _t['C2']) and _t['compromise_set']==['A1','A2']
_t=vikor_compromise_set([0.0,0.2,0.8],[0.0,0.7,1.0],[0.0,0.8,1.0],_test_names,0.5)
assert (not _t['C1']) and _t['compromise_set']==['A1','A2']


df = pd.read_csv(DATA)
assert len(df)==10 and not df.isna().any().any()
assert set(df.loc[df.pareto_nondominated.astype(bool),'location']) == BASE_PARETO
names = df.location.to_numpy()
Rraw = df.mean_ghi_kwh_m2_day.to_numpy(float)
Craw = df.climate_context_equal_weight.to_numpy(float)

# Primary two-objective benchmark: same two dimensions as the current Pareto screen.
# Convert both to a common [0,1] benefit-oriented scale across the fixed ten-site benchmark.
R = (Rraw-Rraw.min())/(Rraw.max()-Rraw.min())
Cbenefit = (Craw.max()-Craw)/(Craw.max()-Craw.min())
X2 = np.column_stack([R, Cbenefit])

decision = df[['location','country','mean_ghi_kwh_m2_day','climate_context_equal_weight','pareto_nondominated']].copy()
decision['ghi_benefit_minmax'] = R
decision['climate_benefit_minmax'] = Cbenefit
decision.to_csv(OUT/'table_1_mcdm_decision_matrix_2d.csv', index=False)

base = mcdm_common(X2, [0.5,0.5], V)
rank = decision.copy()
rank['WSM_score'] = base['WSM']; rank['WSM_rank'] = ranks(base['WSM'], True)
rank['TOPSIS_closeness'] = base['TOPSIS']; rank['TOPSIS_rank'] = ranks(base['TOPSIS'], True)
rank['VIKOR_S'] = base['VIKOR_S']; rank['VIKOR_R'] = base['VIKOR_R']; rank['VIKOR_Q'] = base['VIKOR_Q']; rank['VIKOR_rank'] = ranks(base['VIKOR_Q'], False)
rank = rank.sort_values(['WSM_rank','TOPSIS_rank','VIKOR_rank','location'])
rank.to_csv(OUT/'table_2_equal_weight_mcdm_rankings.csv', index=False)

# Rank correlations.
rk = pd.DataFrame({
    'WSM': ranks(base['WSM'], True),
    'TOPSIS': ranks(base['TOPSIS'], True),
    'VIKOR': ranks(base['VIKOR_Q'], False),
}, index=names)
corr = rk.corr(method='spearman')
corr.to_csv(OUT/'table_3_spearman_rank_correlations.csv')

# VIKOR compromise conditions at equal weights.
DQ = 1/(len(df)-1)
vd = vikor_compromise_set(base['VIKOR_Q'], base['VIKOR_S'], base['VIKOR_R'], names, DQ)
vik_order=vd['order']; q1=vd['q1']; q2=vd['q2']; advantage=vd['C1']; stability=vd['C2']
qwinner=vd['Q_best_site']; S_winner=vd['S_best_site']; R_winner=vd['R_best_site']; compromise=vd['compromise_set']
vikor_df = pd.DataFrame([{
    'v':V,'DQ':DQ,'Q_best_site':qwinner,'Q_best':q1,'Q_second_site':names[vik_order[1]],'Q_second':q2,
    'acceptable_advantage_C1':advantage,'S_best_site':S_winner,'R_best_site':R_winner,
    'acceptable_stability_C2':stability,'compromise_set':'; '.join(compromise)
}])
vikor_df.to_csv(OUT/'table_4_vikor_compromise_test.csv', index=False)

# Resource-vs-climate weight sensitivity, lambda = resource weight.
grid = np.round(np.linspace(0,1,1001), 3)
rows=[]
winners={m:[] for m in ['WSM','TOPSIS','VIKOR']}
top2sets={m:[] for m in ['WSM','TOPSIS','VIKOR']}
for lam in grid:
    res=mcdm_common(X2,[lam,1-lam],V)
    for method,key,hb in [('WSM','WSM',True),('TOPSIS','TOPSIS',True),('VIKOR','VIKOR_Q',False)]:
        order=top_order(res[key],hb)
        top1=names[order[0]]; top2=names[order[1]]
        winners[method].append(top1); top2sets[method].append(frozenset([top1,top2]))
        rows.append({'resource_weight_lambda':lam,'climate_weight':1-lam,'method':method,
                     'rank_1':top1,'rank_2':top2,'top2_matches_pareto_set':frozenset([top1,top2])==frozenset(BASE_PARETO)})
pd.DataFrame(rows).to_csv(OUT/'table_5_weight_sensitivity_grid.csv', index=False)

summ=[]
for method in winners:
    for lo,hi,label in compress_intervals(grid,winners[method]):
        summ.append({'method':method,'interval_type':'rank_1_site','resource_weight_from':lo,'resource_weight_to':hi,'value':label})
    flags=[x==frozenset(BASE_PARETO) for x in top2sets[method]]
    for lo,hi,label in compress_intervals(grid,flags):
        if label:
            summ.append({'method':method,'interval_type':'top2_exactly_pareto_set','resource_weight_from':lo,'resource_weight_to':hi,'value':'Antofagasta; Aswan'})
pd.DataFrame(summ).to_csv(OUT/'table_6_weight_sensitivity_intervals.csv', index=False)

# Normalization-method audit: standard raw-vector TOPSIS and raw ideal-gap VIKOR.
Xraw=np.column_stack([Rraw,Craw])
w=np.array([0.5,0.5])
N=Xraw/np.sqrt((Xraw**2).sum(axis=0)); Y=N*w
ip=np.array([Y[:,0].max(),Y[:,1].min()]); im=np.array([Y[:,0].min(),Y[:,1].max()])
dp=np.linalg.norm(Y-ip,axis=1); dm=np.linalg.norm(Y-im,axis=1); topsis_raw=dm/(dp+dm)
reg=np.column_stack([(Rraw.max()-Rraw)/(Rraw.max()-Rraw.min()), (Craw-Craw.min())/(Craw.max()-Craw.min())])*w
S=reg.sum(axis=1); RR=reg.max(axis=1)
Q=V*(S-S.min())/(S.max()-S.min())+(1-V)*(RR-RR.min())/(RR.max()-RR.min())
norm_audit=pd.DataFrame({'location':names,
                         'shared_minmax_TOPSIS_rank':ranks(base['TOPSIS'],True),
                         'raw_vector_TOPSIS_rank':ranks(topsis_raw,True),
                         'shared_minmax_VIKOR_rank':ranks(base['VIKOR_Q'],False),
                         'raw_ideal_gap_VIKOR_rank':ranks(Q,False)})
norm_audit.to_csv(OUT/'table_7_normalization_audit.csv', index=False)

# Conventional vector-normalized TOPSIS weight sensitivity.
vector_rows=[]
for lam in grid:
    wv=np.array([lam,1-lam], float)
    Nv=Xraw/np.sqrt((Xraw**2).sum(axis=0))
    Yv=Nv*wv
    ipv=np.array([Yv[:,0].max(),Yv[:,1].min()])
    imv=np.array([Yv[:,0].min(),Yv[:,1].max()])
    dpv=np.linalg.norm(Yv-ipv,axis=1); dmv=np.linalg.norm(Yv-imv,axis=1)
    ccv=np.divide(dmv,dpv+dmv,out=np.full_like(dmv,0.5),where=~np.isclose(dpv+dmv,0))
    order=top_order(ccv,True)
    vector_rows.append({'resource_weight_lambda':lam,'climate_weight':1-lam,
                        'rank_1':names[order[0]],'rank_2':names[order[1]],
                        'top2_matches_pareto_set':frozenset(names[order[:2]])==frozenset(BASE_PARETO)})
vector_df=pd.DataFrame(vector_rows)
vector_df.to_csv(OUT/'table_8_vector_TOPSIS_weight_sensitivity.csv', index=False)

# Five-criterion diagnostic: resource plus four climate components, all benefit-oriented.
Z=df[['normalized_heat','normalized_humidity','normalized_ghi_variability','normalized_wind_exposure']].to_numpy(float)
X5=np.column_stack([R,1-Z])
w5=np.array([0.5,0.125,0.125,0.125,0.125])
res5=mcdm_common(X5,w5,V)
five=pd.DataFrame({'location':names,
                   'resource_weight':0.5,'each_climate_component_weight':0.125,
                   'WSM_score':res5['WSM'],'WSM_rank':ranks(res5['WSM'],True),
                   'TOPSIS_closeness':res5['TOPSIS'],'TOPSIS_rank':ranks(res5['TOPSIS'],True),
                   'VIKOR_Q':res5['VIKOR_Q'],'VIKOR_rank':ranks(res5['VIKOR_Q'],False)})
five.to_csv(OUT/'table_9_five_criterion_diagnostic.csv', index=False)

# Five-criterion lambda sensitivity.
frows=[]
for lam in grid:
    w5=np.array([lam]+[(1-lam)/4]*4)
    rr=mcdm_common(X5,w5,V)
    for method,key,hb in [('WSM','WSM',True),('TOPSIS','TOPSIS',True),('VIKOR','VIKOR_Q',False)]:
        order=top_order(rr[key],hb)
        frows.append({'resource_weight_lambda':lam,'method':method,'rank_1':names[order[0]],'rank_2':names[order[1]],
                      'top2_matches_pareto_set':frozenset(names[order[:2]])==frozenset(BASE_PARETO)})
pd.DataFrame(frows).to_csv(OUT/'table_10_five_criterion_weight_sensitivity_grid.csv', index=False)

# Verification assertions.
for method,col in [('WSM','WSM_rank'),('TOPSIS','TOPSIS_rank'),('VIKOR','VIKOR_rank')]:
    top2=set(rank.nsmallest(2,col).location)
    assert top2 == BASE_PARETO, (method,top2)
assert compromise == ['Antofagasta','Aswan'], compromise
assert corr.loc['WSM','TOPSIS'] > 0.95 and corr.loc['WSM','VIKOR'] > 0.95 and corr.loc['TOPSIS','VIKOR'] > 0.95

vector_winners=vector_df.rank_1.tolist()
vector_intervals=compress_intervals(grid, vector_winners)
vector_topsis_rank1_intervals=[{'from':lo,'to':hi,'site':site} for lo,hi,site in vector_intervals]

summary={
 'scope':'Secondary all-site methodological benchmark supporting R4R7; Pareto remains the primary decision rule.',
 'sites':len(df),
 'primary_criteria':['mean GHI (benefit)','equal-weight climate context C (cost)'],
 'normalization':'Common benchmark min-max normalization to [0,1] benefit orientation for both criteria.',
 'baseline_weights':{'resource':0.5,'climate':0.5},
 'vikor_v':V,
 'pareto_set':sorted(BASE_PARETO),
 'equal_weight_rank_1':{
   'WSM':rank.sort_values('WSM_rank').iloc[0].location,
   'TOPSIS':rank.sort_values('TOPSIS_rank').iloc[0].location,
   'VIKOR':rank.sort_values('VIKOR_rank').iloc[0].location,
 },
 'equal_weight_rank_2':{
   'WSM':rank.sort_values('WSM_rank').iloc[1].location,
   'TOPSIS':rank.sort_values('TOPSIS_rank').iloc[1].location,
   'VIKOR':rank.sort_values('VIKOR_rank').iloc[1].location,
 },
 'vikor_compromise_set':compromise,
 'vector_TOPSIS_equal_weight_rank_1':names[top_order(topsis_raw,True)[0]],
 'vector_TOPSIS_equal_weight_rank_2':names[top_order(topsis_raw,True)[1]],
 'vector_TOPSIS_rank1_intervals':vector_topsis_rank1_intervals,
 'spearman':{a:{b:float(corr.loc[a,b]) for b in corr.columns} for a in corr.index},
 'interpretation':'All three equal-weight methods place the two Pareto sites first and second. The scalar winner is normalization/weight dependent, supporting retention of Pareto as the central non-compensatory screen.'
}
(OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')

rows=[]
for p in sorted(ROOT.rglob('*')):
    if p.is_file() and p.name not in {'sha256_manifest.csv','run_log.txt'}:
        rows.append({'relative_path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
pd.DataFrame(rows).to_csv(ROOT/'sha256_manifest.csv',index=False)
print(json.dumps(summary,indent=2))
