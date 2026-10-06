"""Report saved complete-estimator checks, including unfavorable comparisons."""
import csv, math, os
from pathlib import Path
root = Path(os.environ['SIM_OUTPUT'])
reference = Path(os.environ['SIM_REFERENCE'])
def read(path):
    if not path.stat().st_size: return []
    return list(csv.DictReader(path.open()))
def table(headers, rows):
    return '| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---']*len(headers)) + ' |\n' + ''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows) + '\n'
names = {'binary_longitudinal':'Binary treatment','discrete_dose':'Numerical dose'}
estimates, weights, failures, trials = [], [], [], []
for m in names:
    seed = 5103006 if m == 'binary_longitudinal' else 5203006
    folder = root / f'{m}-n4000-seed{seed}'
    assert 'Estimator error rows 0' in (folder/'STATUS.txt').read_text()
    r = read(folder/'estimates.csv'); assert len(r)==4
    assert all(abs(float(x['mean_sequential_eif']) + float(x['mean_bridge_adjoint_eif']) - float(x['estimate']) + float(x['truth'])) < 1e-10 for x in r)
    estimates += r
    weights += [dict(x,mechanism=m) for x in read(folder/'all-weights.csv')]
    failures += read(folder/'nested-candidate-failures.csv')
    trials += read(folder/'nested-penalty-trial-failures.csv')
metrics = read(root/'function-and-equation-errors.csv')
old = read(reference/'function-and-equation-errors.csv')
def error(rows,m,c):
    vals=[float(x['rmse'])**2 for x in rows if x['mechanism']==m and x['candidate']==c and x['kind']=='beta' and x['horizon']=='2' and x['current_R']=='1']
    assert len(vals)==3
    return math.sqrt(sum(vals)/3)
penalties=read(root/'sieve-selected-penalties.csv'); assert len(penalties)==24
scoring=read(root/'scorer-audit.csv'); assert len(scoring)==72
clips=read(root/'applied-bridge-truncation.csv'); assert len(clips)==12 and all(int(x['validation_clipped_n'])==0 for x in clips)
text='''# Full conditional equations: saved n = 4,000 checks

**This revision is not adopted.** The numerical-dose bridge is less accurate: its ensemble RMSE increases from 0.246 to 0.384. The binary bridge also worsens modestly. The original frozen cloud study is unchanged; these diagnostic checks are never pooled with it.

Both datasets finished with four estimates and zero final estimator errors. Binary treatment took 943.598 seconds (15.73 minutes); numerical dose took 900.135 seconds (15.00 minutes). Seeds are 5103006 and 5203006. The data and all outer and learner assignments match the earlier saved checks exactly.

The bridge sieve and Landweber now condition on all observed joint categories. Their target function classes are unchanged. Penalty selection and ensemble weighting use the cell U-statistic, retaining every conditioning predictor, including treatment, and excluding self-products. Only the ensemble Gram is projected to positive semidefinite. Positive sieve penalties are selected inside training samples; boundary minima extend the grid or become recorded failures. Landweber conditioning weights in this revision still have the original fixed ridge, 1e-8.

Packages remain cmbridge 0.3.0.9020 and modified lmtp 1.6.0.9021. All bridge candidates use inverse-expit and unrestricted real coefficients; adjoints retain identity links. SuperLearner regression retains saturated L1 with CV, MARS and mean; treatment classification retains saturated L1 and mean. Saturated L1 is absent from cmbridge. Both SDR and logistic TMLE use the one-step bridge correction and shared training-only splits, including strict rebuilding of sequential responses. No previous nuisance fits were copied into these checks.

## True versus estimated bridge

All combinations and all three training fits appear on equal axes, without clipping. The displayed bridge is unique where both follow-up visits are observed. Elsewhere its full conditional equation is the relevant requirement, rather than agreement with one chosen solution.

![Binary bridge](binary_longitudinal-bridge-ensemble.png)

![Numerical-dose bridge](discrete_dose-bridge-ensemble.png)

'''
text += table(['Example','Candidate','Previous RMSE','Full-equation RMSE'],[[names[m],c,f'{error(old,m,c):.6f}',f'{error(metrics,m,c):.6f}'] for m in names for c in ['ensemble','sieve_md','landweber','pmmr']])
text += 'RMSE is population-probability weighted within the observed intermediate/final-visit group, then averaged in squared units across the three training fits. It measures function error, not confidence-interval coverage.\n\n![Binary candidates](binary_longitudinal-all-bridge-learners.png)\n\n![Dose candidates](discrete_dose-all-bridge-learners.png)\n\n'
for kind,label in [('beta','Bridge'),('adjoint','Adjoint')]:
    rows=[]
    for m in names:
        for fold in range(1,4):
            w={x['candidate']:float(x['weight']) for x in weights if x['mechanism']==m and x['kind']==kind and x['horizon_or_depth']=='2' and x['fold']==str(fold)}
            assert set(w)=={'sieve_md','landweber','pmmr'}
            rows.append([names[m],fold]+[f'{w[c]:.6f}' for c in ['sieve_md','landweber','pmmr']])
    text += f'## {label} weights at time 3\n\n' + table(['Example','Training fit','Sieve','Landweber','PMMR'],rows)
text += '## Estimates and sample mean EIF components\n\n' + table(['Example','Outcome time','Estimator','Truth','Estimate','SE','Sequential EIF mean','Bridge/adjoint EIF mean'],[[names[r['mechanism']],int(r['horizon'])+1,r['estimator'].upper()]+[f'{float(r[c]):.6f}' for c in ['truth','estimate','se','mean_sequential_eif','mean_bridge_adjoint_eif']] for r in estimates])
misses=sum(not float(r['lower'])<=float(r['truth'])<=float(r['upper']) for r in estimates)
assert misses==2
text+='Two of eight pointwise intervals exclude the truth, both for the numerical-dose outcome at time 2 in the same dataset. All time-3 pointwise intervals include it. One dataset per mechanism cannot estimate coverage or establish a convergence rate. EIF component means above are actual sample means; their sum equals estimate minus truth.\n\n'
text+='## Defining equations and audit\n\n![Binary adjoint equations](binary_longitudinal-adjoint-equations.png)\n\n![Dose adjoint equations](discrete_dose-adjoint-equations.png)\n\n'
text+='Adjoint diagnostics check the conditional equation, permitting multiple solutions. Full-support bridge class containment was verified algebraically, but this does not guarantee accurate fitted coefficients. The finite-sample adjoint sieve dictionary still omits combinations absent from its training data, so blanket finite-sample class-containment claims are inappropriate.\n\n'
text+=f'The audits verify unchanged data and assignments, shared training-only penalty splits, positive interior selected penalties, all 24 final and 72 learner-training sieve choices, existing final solver checks, all 24 raw ensemble Grams independently reconstructed with ordered-pair cell formulas, and PSD projected matrices. There is no bridge clipping in cmbridge and no values changed by the approved [-100,100] application bounds. The final estimators record {len(failures)} nested candidate-failure event and {len(trials)} failed penalty-trial records; all are retained, including a tiny-sample boundary failure. These counts are records, not independent statistical events.\n\n'
text+='The source snapshot retains historical audit helpers for reference. Separate current scripts implement the actual cell-scored audit; frozen sources were not silently rewritten to make historical Gaussian assumptions pass.\n\n[All function/equation errors](function-and-equation-errors.csv) · [Score audit](scorer-audit.csv) · [Positive penalties](sieve-selected-penalties.csv) · [Applied bounds audit](applied-bridge-truncation.csv).\n\nThe next diagnostic examines finite-sample conditioning weights and selects any proposed tuning parameter using training-only validation, never true population functions. The running cloud study and all failed/unfavorable attempts remain preserved.\n'
(root/'REPORT.md').write_text(text)
print('Reported eight complete estimates, all failures, and the unfavorable revision.')
