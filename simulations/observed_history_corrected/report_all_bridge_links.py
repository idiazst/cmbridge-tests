from pathlib import Path
import csv,math,collections,json,hashlib,os
root=Path(os.environ.get('SIM_OUTPUT','results/observed-history-corrected/selection-inverse-expit-all-bridges-v4'))
reference=root.parent/'selection-inverse-expit-final'
read=lambda p:list(csv.DictReader(p.open()))
label=lambda m:'Binary' if m=='binary_longitudinal' else 'Numerical dose'
fmt=lambda v:f'{float(v):.6f}'
def table(headers,rows):
 return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in rows)+'\n'
metrics=read(root/'function-and-equation-errors.csv')
selected=[r for r in metrics if r['kind']=='beta' and r['horizon']=='2' and r['current_R']=='1']
groups=collections.defaultdict(list)
for r in selected:groups[(r['mechanism'],int(r['n']),r['candidate'])].append(float(r['rmse'])**2)
errorrows=[]
for m in ['binary_longitudinal','discrete_dose']:
 for n in [4000,20000]:
  errorrows.append([label(m),f'{n:,}']+[f'{math.sqrt(sum(groups[(m,n,c)])/len(groups[(m,n,c)])):.4f}' for c in ['sieve_md','landweber','pmmr','ensemble']])
weights=[];estimates=[];completed=[]
for m in ['binary_longitudinal','discrete_dose']:
 for n in [4000,20000]:
  seed=5103006 if m=='binary_longitudinal' else 5203006
  folder=root/f'{m}-n{n}-seed{seed}'
  status=(folder/'STATUS.txt').read_text();assert 'Estimator error rows 0' in status
  r=read(folder/'estimates.csv');assert len(r)==4;estimates.extend(r)
  weights.extend(dict(r, mechanism=m, n=str(n)) for r in read(folder/'all-weights.csv'))
  completed.append([label(m),f'{n:,}','Completed',0])
with (root/'conditional-moment-ensemble-weights.csv').open('w') as f:
 fields=list(weights[0]);w=csv.DictWriter(f,fields);w.writeheader();w.writerows([r for r in weights if r['kind'] in ['beta','adjoint']])
weightrows={kind:[] for kind in ['beta','adjoint']}
for m in ['binary_longitudinal','discrete_dose']:
 for n in [4000,20000]:
  for kind in weightrows:
   for fold in [1,2,3]:
    z=[r for r in weights if r['mechanism']==m and int(r['n'])==n and r['kind']==kind and int(r['horizon_or_depth'])==2 and int(r['fold'])==fold]
    value={r['candidate']:float(r['weight']) for r in z}
    assert set(value)=={'sieve_md','landweber','pmmr'}
    weightrows[kind].append([label(m),f'{n:,}',fold]+[f'{value[c]:.4f}' for c in ['sieve_md','landweber','pmmr']])
eifrows=[]
for r in estimates:
 assert abs(float(r['mean_sequential_eif'])+float(r['mean_bridge_adjoint_eif'])-(float(r['estimate'])-float(r['truth'])))<1e-10
 eifrows.append([label(r['mechanism']),f'{int(r["n"]):,}',int(r['horizon'])+1,r['estimator'].upper()]+[fmt(r[k]) for k in ['truth','estimate','se','mean_sequential_eif','mean_bridge_adjoint_eif']])
miss=sum(not(float(r['lower'])<=float(r['truth'])<=float(r['upper'])) for r in estimates)
gh='https://github.com/idiazst/cmbridge-tests'
run=gh+'/actions/runs/37398481005';results=gh+'/tree/main/results/run-37398481005'
text='''# Inverse-expit for all bridge candidates\n\nAll three bridge candidates—saturated joint-category sieve, Landweber and PMMR—fit inverse-expit inside their conditional-moment objectives, using unrestricted real coefficients. Their predictions and convex ensemble are at least one. All three adjoint candidates use identity and remain unrestricted, including negative values. Saturated L1 remains absent from both conditional-moment libraries; it remains in SuperLearner regression and classification. MARS remains a regression candidate only. The logistic TMLE update is unchanged.\n\nPackages: **cmbridge 0.3.0.9014** and **lmtp 1.6.0.9015**. Generic standalone cmbridge calls retain their identity defaults for compatibility with the existing continuous recovery specifications; the lmtp integration and this study explicitly select inverse-expit for all three bridge candidates.\n\nThese are the same binary and numerical-dose datasets, with n = 4,000 and n = 20,000, seeds 5103006 and 5203006, and three shared outer and learner-validation groups. Both follow-up outcomes can be missing. Data and split assignments are verified against the preceding validated revision. The data-generating distributions and truths are unchanged. The full repeated-sample simulation remains stopped.\n\n## Fitted versus true bridges\n\nEach figure shows all three candidates and their ensemble, both sample sizes and all three outer training fits. The plots use the support with R₂ = 1, where β₂ is unique. Every reachable combination is retained; no fitted bridge value is clipped for plotting. The axes have equal scales within each panel.\n\n![Binary treatment: all bridge learners](binary_longitudinal-all-bridge-learners.png)\n\n![Numerical dose: all bridge learners](discrete_dose-all-bridge-learners.png)\n\n[Binary PDF](binary_longitudinal-all-bridge-learners.pdf) · [Numerical-dose PDF](discrete_dose-all-bridge-learners.pdf).\n\nThe table uses the same population-probability weighting as the preceding function checks and reports the square root of the mean squared error across the three training fits at R₂ = 1.\n\n'''
text+=table(['Example','n','Sieve','Landweber','PMMR','Ensemble'],errorrows)
text+='''## Does the sieve basis contain the true bridge?\n\nThe plotted sieve uses joint-category indicators in both the target and conditioning variables. These indicators include all interactions among the supplied variables. Inverse-expit assigns a separate unrestricted real coefficient to each represented target combination. An algebraic check using the actual prediction code reproduces the true bridge on every represented combination to an absolute tolerance of 1e-12. No oracle estimator was fitted for this check.\n\nHowever, the fitted dictionary contains only combinations observed in its own training sample. Unrepresented combinations receive a common prediction. Consequently, the realized finite-sample class does not contain the true bridge over its entire population support. Calling this class correctly specified without this qualification was too strong.\n\nThe table summarizes the three outer training samples for the plotted R₂ = 1 stratum. The error fraction uses population-probability-weighted squared error.\n\n'''
support=read(root/'saturated-sieve-support.csv')
membership=read(root/'sieve-representation-audit.csv')
membershiprows=[]
for m in ['binary_longitudinal','discrete_dose']:
 for n in [4000,20000]:
  z=[r for r in support if r['mechanism']==m and int(r['n'])==n and r['kind']=='beta' and r['horizon']=='2' and r['current_R']=='1']
  assert len(z)==3
  represent=[int(r['represented_combinations']) for r in z]
  omitted=[100*float(r['missing_probability']) for r in z]
  error=[100*float(r['error_fraction_in_represented']) for r in z]
  membershiprows.append([label(m),f'{n:,}',z[0]['population_combinations'],f'{min(represent)}–{max(represent)}',f'{min(omitted):.2f}–{max(omitted):.2f}%',f'{min(error):.2f}–{max(error):.2f}%'])
assert max(float(r['maximum_representation_error_seen']) for r in membership)<1e-12
text+=table(['Example','n','Population combinations','Represented combinations','Omitted probability','Error on represented combinations'],membershiprows)
text+='''Most of the sieve error occurs where the class does represent the truth. Missing combinations therefore do not explain most of the scatter. An additional algebraic check of the saved empirical moment matrices found that groups with insufficient rank account for at most about 9.1% of the displayed sieve squared error. This does not establish the cause of the remaining discrepancy: poorly estimated moments, regularization effects, and other implementation issues need to be distinguished before drawing a conclusion. The sieve discrepancy remains unresolved, despite successful recovery tests on the simpler GitHub examples. The saved sieve fits use the default function-value ridge penalty of 1e-8.\n\n[Exact representability audit](sieve-representation-audit.csv) · [Saved empirical moment matrices](sieve-empirical-operator-audit.csv) · [Support and error audit](saturated-sieve-support.csv).\n\n'''
text+='''Scatter remains, especially for the saturated sieve. PMMR remains flattened in these examples; its fixed penalty now acts on the linear predictor inside the link. Passing known-truth recovery tests does not imply that every candidate recovers every function in these longitudinal examples. These are four individual datasets and do not establish confidence-interval coverage.\n\n## Checks and repairs\n\nA nested numerical-dose bridge fit used only 28 training observations. A dose value absent from that subset triggered extreme polynomial spline extrapolation; exponentiating the predictor produced values around 1e189 and overflowed the validation matrix. Fixed iteration steps alone did not resolve that extrapolation. Linked Landweber splines now continue constantly beyond their training boundaries. This is part of their basis definition and does not clip fitted bridge values. Their initial iteration step can decrease through backtracking and never increases. Identity adjoints retain their previous spline continuation. The exact failing split is preserved as a regression fixture; its corrected validation values range from 1.20 to 2.86.\n\nGitHub's Linux runner also exposed insufficient function-scale accuracy near the inverse-expit boundary. The sieve now checks descent on the function scale until the limiting value is approached to the requested accuracy. Finite restarts can approach or leave that limit; coefficients remain unrestricted. The existing test's truth and recovery tolerance are unchanged. Earlier attempts are preserved in the adjacent all-bridge result folders.\n\nThe audits verify finite link coefficients and raw predictions, values at least one for every bridge candidate and ensemble, objective consistency, decreasing Landweber loss, raw cell U-statistics with all conditioning variables including treatment, and positive-semidefinite matrices before weight selection. The applied longitudinal bridge interval remains [-100, 100]; its use and changes are recorded. Both the transformed outcomes and the one-step correction use those same applied values. Adjoint predictions match the preceding validated version.\n\n[All-learner link audit](all-bridge-link-audit.csv) · [Data, split and adjoint invariance](all-bridge-link-audit.txt) · [U-statistic audit](scorer-audit.csv) · [Applied prediction interval](applied-bridge-truncation.csv).\n\n'''
text+=f'''## GitHub validation\n\n[Run 37398481005]({run}) tests package commit **6afde0e6c8757ea8ac611a2a2b181ddf06b73805**. [Tests, numerical results and plots]({results}).\n\n'''
text+=table(['Check','Specification'],[
 ['cmbridge unit tests','All installed package tests, including the saved failing nested split and unchanged boundary-accuracy test'],
 ['Original bridge recovery','Three original candidate specifications; n = 300,000; five seeds'],
 ['Original adjoint recovery','Three original candidate specifications; n = 300,000; five seeds'],
 ['Original ensemble selection','Two original truth scenarios; bridge and adjoint; n = 2,000, 20,000 and 200,000; five seeds'],
 ['Inverse-expit bridge / signed adjoint','All three learners; four-category truth; n = 300,000; five seeds']])
text+='''The existing continuous recovery specifications and thresholds are retained. The additional check explicitly uses inverse-expit for bridges and identity for adjoints, including a negative adjoint truth. Local package and lmtp integration tests also pass. Existing spline-knot warnings are retained in the logs.\n\n![New inverse-expit and signed-adjoint recovery](github-run/linked-bridge/truth_vs_estimate.png)\n\n![Original ensemble recovery](github-run/ensemble/truth_vs_estimate.png)\n\n[Original bridge recovery plots](github-run/bridge/truth_vs_estimate.png) · [Original adjoint recovery plots](github-run/adjoint/truth_vs_estimate.png).\n\n## Ensemble weights at outcome time 3\n\n### Bridge\n\n'''
text+=table(['Example','n','Training fit','Sieve','Landweber','PMMR'],weightrows['beta'])
text+='### Adjoint\n\n'+table(['Example','n','Training fit','Sieve','Landweber','PMMR'],weightrows['adjoint'])
text+='[Complete weights at both outcome times](conditional-moment-ensemble-weights.csv).\n\n## Sequential estimates and one-step correction\n\n'
text+=table(['Example','n','Status','Fit errors'],completed)
text+='The sequential mean below is centered at the true parameter. Adding the signed bridge/adjoint mean gives estimate minus truth.\n\n'
text+=table(['Example','n','Outcome time','Estimator','Truth','Estimate','SE','Mean sequential component','Mean bridge/adjoint component'],eifrows)
text+=f'{miss} of the 16 pointwise intervals exclude their true value. Every result is retained; these individual checks are not a coverage study.\n\n'
text+='''## Reproduction\n\nFrozen package sources, archives, installed private library, and simulation scripts are in `source`, `packages` and `library`. `check_selection_revision.R` records full fitted objects and reuses the reference data/splits; `plot_all_bridge_learners.R` draws the displayed figures. Individual predictions and defining conditional-equation values are in [all-function-predictions.csv](all-function-predictions.csv) and [all-equation-predictions.csv](all-equation-predictions.csv). The adjoint equation, rather than a particular nonunique numerical-dose adjoint solution, is used for its equation plots.\n\nThe earlier [population uniqueness calculation](../selection-inverse-expit-final/UNIQUENESS.md) is unchanged: β₁ is unique in both mechanisms; β₂ is unique at R₂ = 1 and nonunique at R₂ = 0. Both mechanisms admit bridge and adjoint solutions.\n'''
(root/'REPORT.md').write_text(text)
(root/'README.md').write_text('# Inverse-expit for all bridge candidates\n\nSee [REPORT.md](REPORT.md) for the four datasets, plots, all learner weights, EIF means, audits, and successful GitHub validation. The full repeated-sample study remains stopped.\n')
(root/'validation/report-integrity.txt').write_text(f'Four completed datasets; {len(estimates)} estimates; zero fit errors.\nAll EIF component sums agree with estimate minus truth to 1e-10.\n')
print('Wrote report with',len(estimates),'estimates and',len(weightrows['beta']),'bridge-weight rows.')
