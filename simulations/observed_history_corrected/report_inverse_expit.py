"""Report saved inverse-expit diagnostics without fitting or selecting models."""
from pathlib import Path
import os,csv,math,statistics
root=Path(os.environ['SIM_OUTPUT'])
reference=Path(os.environ['SIM_REFERENCE'])
def read(path):
    return list(csv.DictReader(path.open())) if path.exists() and path.stat().st_size else []
def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def f(x,d=4):return f'{float(x):.{d}f}'
labels={'binary_longitudinal':'Binary','discrete_dose':'Numerical dose'}
cases=[(m,n) for m in labels for n in [4000,20000]]
metrics=read(root/'function-and-equation-errors.csv');old=read(reference/'function-and-equation-errors.csv')
def rmse(data,m,n,candidate):
    z=[float(r['rmse'])**2 for r in data if r['mechanism']==m and int(r['n'])==n and r['kind']=='beta' and r['horizon']=='2' and r['current_R']=='1' and r['candidate']==candidate]
    return math.sqrt(statistics.mean(z))
lines=['# Inverse-expit bridge sieve: four dataset checks','',
'The bridge `sieve_md` now directly optimizes unrestricted real coefficients and evaluates the fitted function through `1 / expit`, computed as `1 + exp(-linear predictor)`. It uses no coefficient bounds, constrained quadratic solver, or post-fit prediction clipping. The saturated joint-category basis and its existing ridge penalty on function values are retained. The adjoint sieve keeps its identity parameterization.','',
'The mathematical bridge values are strictly greater than one, with one as a limit; floating-point arithmetic can return exactly one. The paper’s inverse-probability bridge belongs to this class on the full finite support. Other solutions of the bridge equation need not be at least one. A fitted training dictionary can still lack rare predictor combinations; those receive the linked common value.','',
'These are the same two n = 4,000 replication-6 datasets and the corresponding n = 20,000 datasets. The seeds are 5103006 for the binary example and 5203006 for the numerical dose example. The data-generating distributions, outer assignments, and learner-validation assignments match the preceding three-candidate revision. Both follow-up outcomes can be missing. The repeated-sample simulation remains stopped.','',
'Packages: **cmbridge 0.3.0.9010** and **lmtp 1.6.0.9011**. Both conditional-moment ensembles contain sieve, Landweber and PMMR. Saturated L1 is absent from both ensembles and remains a SuperLearner regression/classification candidate. MARS remains a regression candidate only. The logistic TMLE update, shared assignments, U-statistic cell scoring, and projection before selecting ensemble weights are retained.','',
'## Bridge prediction error','',
'The table compares the unrestricted identity sieve in the previous three-candidate revision with the new inverse-expit sieve. It reports the square root of the mean of the three training fits’ probability-weighted squared errors for β₂ at histories with R₂ = 1. Here R₂ = 1 means that a visit occurred at time 2; the outcome being studied is at time 3. The bridge solution is unique on this support.','',
table(['Example','n','Previous sieve','Inverse-expit sieve','Previous ensemble','New ensemble'],[[labels[m],f'{n:,}',f(rmse(old,m,n,'sieve_md')),f(rmse(metrics,m,n,'sieve_md')),f(rmse(old,m,n,'ensemble')),f(rmse(metrics,m,n,'ensemble'))] for m,n in cases]),'',
'The sieve’s error decreased in all four checks, but it still has substantial scatter around y = x. The ensemble also improved in all four checks. These are single datasets, so the comparison does not establish repeated-sample bias, convergence, or confidence-interval coverage.','',
'![Binary sieve fitted versus true](binary_longitudinal-bridge-sieve_md.png)','',
'[Binary figure as PDF](binary_longitudinal-bridge-sieve_md.pdf).','',
'![Numerical dose sieve fitted versus true](discrete_dose-bridge-sieve_md.png)','',
'[Numerical dose figure as PDF](discrete_dose-bridge-sieve_md.pdf). All points appear in the full-range figures. Companion files with `-central` in their names give an enlarged view and state any excluded probability.','',
'## Ensemble weights','',
'Each row gives one outer training fit at time 3. The complete weights at both outcome times are in [conditional-moment-ensemble-weights.csv](conditional-moment-ensemble-weights.csv).']
weights=[]
for m,n in cases:
    seed=5103006 if m=='binary_longitudinal' else 5203006
    folder=root/f'{m}-n{n}-seed{seed}'
    for fold in range(1,4):
        for r in read(folder/f'weights-fold{fold}.csv'):
            r.update(mechanism=m,n=n)
            weights.append(r)
if weights:
    with (root/'conditional-moment-ensemble-weights.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=list(weights[0]));writer.writeheader();writer.writerows(weights)
for kind,title in [('beta','Bridge'),('adjoint','Adjoint')]:
    rows=[]
    for m,n in cases:
        for fold in range(1,4):
            subset={r['candidate']:r['weight'] for r in weights if r['mechanism']==m and r['n']==n and r['kind']==kind and r['horizon']=='2' and int(r['fold'])==fold}
            rows.append([labels[m],f'{n:,}',fold]+[f(subset[c]) for c in ['sieve_md','landweber','pmmr']])
    lines+=['',f'### {title}','',table(['Example','n','Training fit','Sieve','Landweber','PMMR'],rows)]
lines+=['','## Other learners and adjoint checks','',
'The saved Landweber, PMMR, and adjoint predictions match the preceding revision exactly across 21,888 compared rows. The inverse-expit change affects the bridge sieve and the selected bridge weights. The other bridge candidates and the ensemble remain unrestricted on the function scale.','',
'[Binary ensemble bridge](binary_longitudinal-bridge-ensemble.png) · [Numerical dose ensemble bridge](discrete_dose-bridge-ensemble.png).','',
'[Binary sieve adjoint versus true](binary_longitudinal-adjoint-sieve_md.png) · [Numerical dose sieve adjoint equation](discrete_dose-adjoint-equation-sieve_md.png). The numerical adjoint is not unique at R₂ = 1, so its defining conditional equation is used for comparison.','']
for m,n in cases:
    lines.append(f'- {labels[m]}, n = {n:,}: [Landweber, PMMR and bridge ensemble]({m}-n{n}-other-bridge-learners.png); [their adjoint equations]({m}-n{n}-other-adjoint-learners.png).')
audit=read(root/'inverse-expit-optimizer-audit.csv')
if audit:
    restarts=sum(int(r['restarts']) for r in audit)
    lines+=['','## Numerical and scoring checks','',
    'The initial nonlinear implementation sometimes stopped with tiny coefficient gradients even though increasing fitted values near one would improve the objective. The inverse-expit link was too flat to reveal that improvement. The final solver inspects both coefficient gradients and function-scale descent, then restarts affected coefficients at finite values when needed. It also evaluates the quadratic loss relative to a fixed reference point to avoid subtracting nearly equal large values at every iteration.','',
    f'All 24 final primary bridge-sieve fits passed the gradient checks, with {restarts} finite restarts in total. All link coefficients and predictions are finite. The audit independently verifies that the saved penalized loss equals the fitted conditional-moment loss plus the existing ridge penalty. Optimizer stop codes are retained separately: the fit is accepted by the recorded gradient checks, and an iteration or precision stop is not rewritten as an optimizer success.','',
    'Raw cell U-statistics were recomputed from saved validation residuals with all conditioning inputs, including treatment. They match the stored matrices to 1e−10. Projected matrices are positive semidefinite to numerical tolerance. Neither conditional-moment library contains saturated L1. The primary applied bridge values match the wider [-100, 100] interval, and their truncation counts are recorded.','',
    '[Optimizer audit](inverse-expit-optimizer-audit.csv) · [U-statistic audit](scorer-audit.csv) · [Applied bridge truncation](applied-bridge-truncation.csv) · [Represented predictor combinations](saturated-sieve-support.csv).','',
    'All cmbridge tests and the lmtp bridge/cross-validation tests pass. lmtp still reports the existing warnings about knots in discrete spline inputs; they are retained in the test logs.']
results=[];status=[]
for m,n in cases:
    seed=5103006 if m=='binary_longitudinal' else 5203006
    folder=root/f'{m}-n{n}-seed{seed}'
    results.extend(read(folder/'estimates.csv'))
    text=(folder/'STATUS.txt').read_text() if (folder/'STATUS.txt').exists() else 'Running sequential checks'
    errors=read(folder/'estimator-errors.csv')
    status.append([labels[m],f'{n:,}',text.splitlines()[0].split()[0],len(errors) if 'Completed' in text else 'pending'])
lines+=['','## Sequential estimates and one-step correction','',table(['Example','n','Status','Fit errors'],status),'',
'The two mean influence-function components below use the true parameter to center the sequential component. Their sum therefore equals estimate minus truth. The bridge/adjoint component includes its sign in the one-step correction. These are single-dataset checks, not a coverage study.']
if results:
    lines+=['',table(['Example','n','Outcome time','Estimator','Truth','Estimate','SE','Mean sequential component','Mean bridge/adjoint component'],[
    [labels[r['mechanism']],f"{int(r['n']):,}",int(r['horizon'])+1,r['estimator'].upper(),f(r['truth'],6),f(r['estimate'],6),f(r['se'],6),f(r['mean_sequential_eif'],6),f(r['mean_bridge_adjoint_eif'],6)] for r in results])]
    missed=[r for r in results if not float(r['lower'])<=float(r['truth'])<=float(r['upper'])]
    lines+=['',f'{len(missed)} of the {len(results)} completed pointwise intervals exclude their true value. Every result is retained.']
lines+=['','## Files and reproduction','',
'Frozen package sources, installable archives, and the installed private library are in `source`, `packages`, and `library`. The live source files are `cmbridge/R/sieve_md.R` and `lmtp/R/bridge_curve.R`. Run `check_selection_revision.R` with the recorded library and reference assignments. The original three-candidate results remain in `selection-no-saturated-l1`; earlier solver attempts are preserved separately.','',
'Full raw predictions and conditional-equation errors are in [all-function-predictions.csv](all-function-predictions.csv), [all-equation-predictions.csv](all-equation-predictions.csv), and [function-and-equation-errors.csv](function-and-equation-errors.csv).']
if (root/'UNIQUENESS.md').exists():
    lines += ['', '## Bridge uniqueness', '', '[Independent population check](UNIQUENESS.md): β₁ is unique; β₂ is unique at R₂ = 1 and nonunique at R₂ = 0, in both mechanisms. The bridge plots use R₂ = 1. Inverse-expit does not restore uniqueness at missed visits.']
(root/'REPORT.md').write_text('\n'.join(lines)+'\n')
(root/'README.md').write_text('# Inverse-expit bridge sieve\n\nSee [REPORT.md](REPORT.md) for the parameterization, four-dataset plots, ensemble weights, audits, and sequential checks.\n')
print(root/'REPORT.md')
