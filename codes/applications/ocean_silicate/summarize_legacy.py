"""Audit and export saved latest fits; no fitting or prediction calls.

Tables are statistical CSV intermediates; plotting requires system matplotlib.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

if len(sys.argv) != 2:
    raise SystemExit("Usage: python3 summarize_legacy.py OCEAN_OUTPUT_DIRECTORY")
ROOT = Path(sys.argv[1]).resolve()
CONT = ROOT / 'continuation_plus20k'
OUT = CONT / 'summary'
OUT.mkdir(parents=True, exist_ok=True)
CELLS = ['O0', 'O20', 'O50']
BASES = ['UC-GP', 'LVGP', 'EzGP']
ORDER = BASES + [f'{o} EIV-GP {m}' for o in CELLS for m in ['xc', 'mixed']]
METRICS = ['RMSE', 'MAE', 'CRPS', 'NLPD', 'Coverage95', 'Width95', 'IntervalScore95']
cohort = pd.read_csv(ROOT / 'data/cohort_manifest.csv')
ys = cohort.Y.std(ddof=1)
us = cohort.U.std(ddof=1)
ym = cohort.Y.mean()
assert len(cohort) == cohort.cohort_id.nunique() == cohort.profile.nunique() == 200

def source(f, o):
    extended = CONT / f'fold{f:02d}' / o
    return extended if (extended / 'predictive_metrics.csv').exists() else ROOT / 'full5_class4' / f'fold{f:02d}' / o

def table(df):
    columns = list(df.columns)
    lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    for row in df.itertuples(index=False, name=None):
        lines.append('| ' + ' | '.join(f'{x:.4f}' if isinstance(x, (float, np.floating)) else str(x) for x in row) + ' |')
    return '\n'.join(lines)

points, u_points, provenance = [], [], []
for f in range(1, 6):
    base = pd.read_csv(ROOT / 'full5_class4' / f'fold{f:02d}' / 'baseline_predictions.csv')
    base['fold'] = f
    base['test_U_observed'] = False
    points.append(base)
    for o in CELLS:
        src = source(f, o)
        provenance.append(dict(fold=f, cell=o, source=str(src.relative_to(ROOT)), continued=src.is_relative_to(CONT)))
        xc = pd.read_csv(src / 'predictions.csv').set_index('cohort_id')
        mix = pd.read_csv(src / 'mixed_predictions.csv').set_index('cohort_id')
        assert set(xc.index) == set(mix.index)
        mix = mix.loc[xc.index]
        known = mix.test_U_observed.astype(str).str.lower().eq('true')
        assert known.sum() == {'O0': 0, 'O20': 8, 'O50': 20}[o]
        assert np.allclose(xc.observed_Y, mix.observed_Y)
        assert np.allclose(xc.loc[~known, 'predicted_Y'], mix.loc[~known, 'predicted_Y'])
        for mode, frame in [('xc', xc), ('mixed', mix)]:
            p = frame.reset_index()
            p['method'] = f'{o} EIV-GP {mode}'
            p['fold'] = f
            if mode == 'xc':
                p['test_U_observed'] = False
            points.append(p)
        if o != 'O0':
            imp = pd.read_csv(src / 'imputations.csv')
            ints = pd.read_csv(src / 'U_posterior_intervals.csv')
            assert ints.scale.eq('raw').all()
            imp = imp.merge(ints[['cohort_id', 'posterior_mean', 'q025', 'q975']], on='cohort_id', validate='one_to_one')
            assert np.allclose(imp.EIV_U_mean, imp.posterior_mean)
            assert len(imp) == {'O20': 128, 'O50': 80}[o]
            imp['fold'], imp['cell'] = f, o
            u_points.append(imp)

points = pd.concat(points, ignore_index=True)
truth = cohort.set_index('cohort_id').Y
assert np.allclose(points.observed_Y, truth.loc[points.cohort_id])
assert points.groupby('method').size().eq(200).all()
assert not points.duplicated(['method', 'cohort_id']).any()
points['absolute_error'] = abs(points.predicted_Y - points.observed_Y)
points['observed_Y_z'] = (points.observed_Y - ym) / ys
points['predicted_Y_z'] = (points.predicted_Y - ym) / ys
points['absolute_error_over_SD_Y'] = points.absolute_error / ys
points.to_csv(OUT / 'point_predictions.csv', index=False)
pd.DataFrame(provenance).to_csv(OUT / 'fit_sources.csv', index=False)

folds = pd.read_csv(CONT / 'latest_fold_predictive_metrics.csv')
folds['label'] = np.where(folds.cell.eq('baseline'), folds.method, folds.cell + ' ' + folds.method)
summary = pd.read_csv(CONT / 'latest_summary_metrics.csv').set_index('method').loc[ORDER].reset_index()
for row in summary.itertuples():
    p = points[points.method.eq(row.method)]
    assert np.isclose(np.sqrt(np.mean(p.absolute_error ** 2)), row.pooled_RMSE)
    assert np.isclose(p.absolute_error.mean(), row.pooled_MAE)
    f = folds[folds.label.eq(row.method)]
    assert len(f) == 5
    for metric in METRICS:
        assert np.isclose(f[metric].mean(), getattr(row, f'fold_{metric}_mean'))
        assert np.isclose(f[metric].std(ddof=1), getattr(row, f'fold_{metric}_sd'))
summary.to_csv(OUT / 'summary_original_scale.csv', index=False)
normalized = summary.copy()
for metric in METRICS:
    for prefix in ['pooled_', 'fold_']:
        cols = [prefix + metric] if prefix == 'pooled_' else [prefix + metric + '_mean', prefix + metric + '_sd']
        for col in cols:
            if metric == 'NLPD':
                if not col.endswith('_sd'):
                    normalized[col] -= np.log(ys)
            elif metric != 'Coverage95':
                normalized[col] /= ys
normalized.to_csv(OUT / 'summary_SD_Y_scale.csv', index=False)
folds.to_csv(OUT / 'fold_metrics_original_scale.csv', index=False)
fn = folds.copy()
for metric in METRICS:
    if metric == 'NLPD':
        fn[metric] -= np.log(ys)
    elif metric != 'Coverage95':
        fn[metric] /= ys
fn['common_Y_SD'] = ys
fn.to_csv(OUT / 'fold_metrics_SD_Y_scale.csv', index=False)
conv = pd.read_csv(CONT / 'latest_convergence_summary.csv')
assert len(conv) == 15 and conv.limited_screen_pass.astype(str).str.lower().eq('true').all()
conv.to_csv(OUT / 'convergence.csv', index=False)

up = pd.concat(u_points, ignore_index=True)
up.to_csv(OUT / 'u_imputation_points.csv', index=False)
uf = []
for (f, o), group in up.groupby(['fold', 'cell']):
    for name, col in [('EIV-GP', 'EIV_U_mean'), ('Calibration class mean', 'class_calibration_mean')]:
        err = group[col] - group.true_U
        uf.append(dict(fold=f, cell=o, method=name, n=len(group), RMSE=np.sqrt(np.mean(err**2)), MAE=np.mean(abs(err)), common_U_SD=us))
uf = pd.DataFrame(uf)
uf.to_csv(OUT / 'u_imputation_fold_metrics.csv', index=False)
u_summary = uf.groupby(['cell', 'method']).agg(RMSE_mean=('RMSE', 'mean'), RMSE_sd=('RMSE', 'std'), MAE_mean=('MAE', 'mean'), MAE_sd=('MAE', 'std')).reset_index()
u_summary['RMSE_mean_over_SD_U'] = u_summary.RMSE_mean / us
u_summary.to_csv(OUT / 'u_imputation_summary.csv', index=False)

improvements = []
for label in ORDER[3:]:
    for b in BASES:
        sr = summary.set_index('method').loc[label, 'pooled_RMSE']
        br = summary.set_index('method').loc[b, 'pooled_RMSE']
        e = folds[folds.label.eq(label)].set_index('fold').RMSE
        bb = folds[folds.label.eq(b)].set_index('fold').RMSE
        improvements.append(dict(method=label, baseline=b, RMSE_reduction_percent=100*(1-sr/br), RMSE_reduction_over_SD_Y=(br-sr)/ys, better_fold_count=int((e < bb).sum())))
pd.DataFrame(improvements).to_csv(OUT / 'paired_RMSE_comparisons.csv', index=False)

sections = ['# U<180 joint-stratified five-fold results after selected +20k continuations',
    f'Frozen cohort: 200 distinct profiles, 160 training / 40 test per fold. Oxygen <180 µmol/kg; pressure 300–1500 dbar; four oxygen classes at 50/100/150. Shared SD(Y)={ys:.8f} µmol/kg; mean(Y)={ym:.8f}; SD(U)={us:.8f} µmol/kg.',
    'Only fold 1 O0, fold 3 O20 and fold 4 O0 were continued (+20,000 transitions per chain). All other fits and the three baselines remain unchanged. The hard-threshold package engine and u_block_size=8 were retained. All 15 EIV fits pass the limited finite-Rhat≤1.01 / finite-bulk-ESS≥400 screen; this is not a proof of convergence for every posterior quantity.',
    'Only three of the five frozen folds passed all historical R1–R3 gates. All five folds are reported without outcome-based rerolling or exclusion. This selected-cohort experiment is exploratory and does not establish that the gates generalize.',
    '## Interpretation',
    'O50 mixed reduces pooled RMSE by 16.80% relative to EzGP, the best baseline in pooled RMSE; O20 mixed reduces it by 12.98%. O50 xc improves by only 2.08%. Mixed prediction observes test U for 0/20/50% of points; the baselines always use (X,C). Its gain is an operational benefit with extra information, not a strictly equal-input comparison. O0 xc and mixed are identical.',
    'Normalization does not change rankings or relative percentage gains, and does not establish statistical significance. The fold SD is sample SD across five held-out-fold scores, not posterior SD or a significance test. The same SD(Y) computed from the full frozen cohort is used descriptively for every method, not as a modeling or selection input.',
    '## Complete pooled metrics — original scale',
    table(summary[['method'] + ['pooled_' + m for m in METRICS]]),
    '## Complete pooled metrics — SD(Y) scale',
    'RMSE, MAE, CRPS, interval width and interval score are divided by shared SD(Y). Coverage is unchanged. Standardized NLPD = original NLPD − log(SD(Y)); it is not divided by SD(Y).',
    table(normalized[['method'] + ['pooled_' + m for m in METRICS]])]
for title, frame in [('Original scale', summary), ('SD(Y) scale', normalized)]:
    formatted = pd.DataFrame({'Method': frame.method})
    for m in METRICS:
        formatted[m + ' mean ± SD'] = [f'{a:.4f} ± {b:.4f}' for a,b in zip(frame[f'fold_{m}_mean'], frame[f'fold_{m}_sd'])]
    sections += [f'## Five-fold mean ± SD — {title}', table(formatted)]
sections += ['## Every fold — original scale', table(folds[['fold','label']+METRICS]),
    '## Every fold — SD(Y) scale', table(fn[['fold','label']+METRICS]),
    '## U-imputation',
    'Evaluation uses missing-U training cases only: 128 per fold for O20 and 80 for O50. Training sets overlap; their concatenation is only for faceted illustration, not 5× independent observations. Calibration class means use observed training U only. EIV raw-U imputation remains worse than the class-mean comparator. O0 has no raw-U calibration anchor; raw-scale error or a raw-scale identity plot cannot be interpreted, and no true-U-based alignment was performed.',
    table(u_summary), '## Convergence', table(conv),
    '## Figures',
    'Color identifies method, never class. Each Y panel overlays the same 200 held-out observations for EIV-GP and all three baselines. Top row uses (X,C) EIV prediction; bottom row uses mixed EIV prediction. Shared axes show all points without selective clipping.',
    '![Y prediction](y_prediction_overlay_original.png)',
    '![Y prediction standardized](y_prediction_overlay_SD_Y.png)',
    '![Absolute point-wise errors](pointwise_absolute_error_overlay_original.png)',
    '![Absolute errors standardized](pointwise_absolute_error_overlay_SD_Y.png)',
    '![U imputation](u_imputation_all5.png)',
    '![U imputation including O0 identifiability warning](u_imputation_fold1.png)',
    'U error bars are posterior 95% intervals, not fold SD. Posterior support in the extreme classes is open-ended; intervals outside the selected U range are retained rather than clipped.',
    '## Reproducibility',
    'fit_sources.csv resolves original versus continued outputs per fit. point_predictions.csv stores keyed predictions and errors; summary_original_scale.csv and summary_SD_Y_scale.csv contain all pooled and five-fold mean/SD metrics. Per-fold tables, U-imputation tables and convergence.csv are also included. No models were refitted for this report.']
(OUT / 'RESULTS.md').write_text('\n\n'.join(sections) + '\n', encoding='utf-8')
(OUT / 'normalization.json').write_text(json.dumps({'Y_SD_ddof1':ys,'Y_mean':ym,'U_SD_ddof1':us,'n_unique_profiles':200},indent=2))
print(normalized[['method','pooled_RMSE','fold_RMSE_mean','fold_RMSE_sd','pooled_MAE','pooled_CRPS','pooled_NLPD','pooled_Coverage95','pooled_Width95','pooled_IntervalScore95']].to_string(index=False))
print(u_summary.to_string(index=False))

if '--tables-only' in sys.argv:
    sys.exit(0)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
COLORS = {'EIV-GP':'#0072B2','UC-GP':'#009E73','LVGP':'#E69F00','EzGP':'#CC79A7','Calibration class mean':'#D55E00'}
MARKERS = {'EIV-GP':'o','UC-GP':'s','LVGP':'^','EzGP':'D','Calibration class mean':'x'}
def save(fig, name):
    for ext in ['png','pdf']:
        fig.savefig(OUT / (name+'.'+ext), dpi=180, bbox_inches='tight', pad_inches=.15)
    plt.close(fig)
legend = [Line2D([],[],ls='',marker=MARKERS[n],color=COLORS[n],label=n,markersize=7) for n in ['EIV-GP']+BASES]
for standardized in [False, True]:
    unit = 'SD(Y) units' if standardized else 'µmol/kg'
    xcol = 'observed_Y_z' if standardized else 'observed_Y'
    ycol = 'predicted_Y_z' if standardized else 'predicted_Y'
    ecol = 'absolute_error_over_SD_Y' if standardized else 'absolute_error'
    suffix = 'SD_Y' if standardized else 'original'
    extent = np.r_[points[xcol], points[ycol]]
    pad = np.ptp(extent)*.035
    limits = (extent.min()-pad, extent.max()+pad)
    for residual in [False, True]:
        fig, axes = plt.subplots(2,3,figsize=(16,10),sharex=True,sharey=True)
        for row, mode in enumerate(['xc','mixed']):
            for col, o in enumerate(CELLS):
                ax = axes[row,col]
                for name in BASES+['EIV-GP']:
                    label = name if name in BASES else f'{o} EIV-GP {mode}'
                    p = points[points.method.eq(label)]
                    ax.scatter(p[xcol],p[ecol if residual else ycol],color=COLORS[name],marker=MARKERS[name],s=25 if name=='EIV-GP' else 23,alpha=.78 if name=='EIV-GP' else .42,linewidths=.4,label=name)
                if not residual:
                    ax.plot(limits,limits,color='.4',ls='--',lw=1,zorder=0)
                    ax.set_ylim(limits)
                else:
                    ax.set_ylim(-.02*points[ecol].max(),1.06*points[ecol].max())
                ax.set_xlim(limits)
                ax.grid(alpha=.17)
                score = summary.set_index('method').loc[f'{o} EIV-GP {mode}', 'pooled_RMSE'] / (ys if standardized else 1)
                ax.set_title(f'{o.replace("O","O=")}% · {"(X,C)" if mode=="xc" else "mixed"}\nEIV RMSE={score:.3f}')
                ax.set_xlabel(f'Observed silicate ({unit})')
                ax.set_ylabel(f'{"Absolute prediction error" if residual else "Predicted silicate"} ({unit})')
        fig.suptitle('Point-wise absolute prediction errors' if residual else 'Held-out Y prediction: EIV-GP and three baselines',fontsize=16,y=1.01)
        fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,-.02),ncol=4,frameon=False)
        fig.tight_layout(rect=(0,.03,1,.99),h_pad=2)
        save(fig, ('pointwise_absolute_error_overlay_' if residual else 'y_prediction_overlay_')+suffix)

def u_panel(ax, p, title, annotate=True):
    ax.errorbar(p.true_U,p.EIV_U_mean,yerr=np.vstack([p.EIV_U_mean-p.q025,p.q975-p.EIV_U_mean]),fmt='none',ecolor=COLORS['EIV-GP'],alpha=.16,lw=.65)
    for name,col in [('Calibration class mean','class_calibration_mean'),('EIV-GP','EIV_U_mean')]:
        ax.scatter(p.true_U,p[col],c=COLORS[name],marker=MARKERS[name],s=21,alpha=.7,label=name)
    ax.plot(ulimits,ulimits,ls='--',color='.4',lw=1)
    ax.set_xlim(15,185)
    ax.set_ylim(ulimits)
    ax.grid(alpha=.16)
    ax.set_title(title)
    ax.set_xlabel('True oxygen (µmol/kg)')
    ax.set_ylabel('Imputed oxygen (µmol/kg)')
    if annotate:
        a=np.sqrt(np.mean((p.EIV_U_mean-p.true_U)**2)); b=np.sqrt(np.mean((p.class_calibration_mean-p.true_U)**2))
        ax.text(.03,.97,f'RMSE: EIV {a:.2f}; class mean {b:.2f}',transform=ax.transAxes,va='top',fontsize=10)
uext = np.r_[up.q025,up.q975,up.true_U]
ulimits = (uext.min()-5,uext.max()+5)
ulegend = [Line2D([],[],ls='',marker=MARKERS[n],color=COLORS[n],label=n,markersize=7) for n in ['EIV-GP','Calibration class mean']]
fig,axes=plt.subplots(2,5,figsize=(23,9),sharex=True,sharey=True)
for row,o in enumerate(['O20','O50']):
    for col,f in enumerate(range(1,6)):
        u_panel(axes[row,col],up[up.fold.eq(f)&up.cell.eq(o)],f'Fold {f} · {o.replace("O","O=")}%')
fig.suptitle('Missing training-U imputation: posterior means and 95% intervals',fontsize=16)
fig.legend(handles=ulegend,loc='lower center',bbox_to_anchor=(.5,-.025),ncol=2,frameon=False)
fig.tight_layout(rect=(0,.035,1,.95))
save(fig,'u_imputation_all5')
fig,axes=plt.subplots(1,3,figsize=(16,5))
axes[0].axis('off')
axes[0].set_title('O=0% · no raw-U anchor')
axes[0].text(.5,.5,'Raw-scale U-imputation is not identified.\n\nPosterior U is on the model scale.\nNo alignment using true U is applied.',ha='center',va='center',transform=axes[0].transAxes,fontsize=12)
for ax,o in zip(axes[1:],['O20','O50']):
    u_panel(ax,up[up.fold.eq(1)&up.cell.eq(o)],f'Fold 1 · {o.replace("O","O=")}%')
fig.suptitle('U-imputation: fixed fold 1 illustration',fontsize=16)
fig.legend(handles=ulegend,loc='lower center',bbox_to_anchor=(.5,-.03),ncol=2,frameon=False)
fig.tight_layout(rect=(0,.04,1,.94))
save(fig,'u_imputation_fold1')
print(f'Wrote tables, report and figures: {OUT}')
